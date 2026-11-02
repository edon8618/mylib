#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

// Matroid intersection (unweighted / weighted), with solution reconstruction.
// C++17+, STL only. Elements are indexed 0..n-1.
//
// A user-supplied matroid M implements:
//   void build(const std::vector<char>& selected);
//   bool can_add(int x);              // selected + x is independent
//   bool can_exchange(int x, int y);  // selected - y + x is independent
// Here x is NOT selected and y IS selected. The arguments to build() always
// represent an independent set in BOTH matroids. The two matroids may cache
// arbitrary preprocessing in build(), but must not modify selected.
//
// Public API (all return the chosen element indices in ascending order):
//   matroid_intersection_max_cardinality(n, m1, m2)
//   matroid_intersection_max_weight(weights, m1, m2)
//   matroid_intersection_max_cardinality_weight(weights, m1, m2)
//   matroid_intersection(n, m1, m2)   // alias for max_cardinality
// The max_weight mode allows the empty set, and the max_cardinality_weight
// mode maximizes the cardinality first and the weight second.
// Negative and zero weights are supported. The weight type T is deduced from
// std::vector<T>; no infinity sentinel is needed. Integrals up to 64 bits use
// int64 shortest paths when provably safe and exact signed int128 otherwise.
// Floating-point weights use long double internally; near-ties may differ
// due to floating-point rounding.
//
// Complexity: let n be ground size, r maximal intersection size, and B/A/X
// be total build / add / exchange oracle costs for both matroids.
//   Unweighted: O((r+1) * (B + n*A + n*r*X + n*r)).
//   Weighted:   O((r+1) * (B + n*A + n*r*(X+log(n+1)) + n*n)).
//                The n*n term covers tiny-n dense Dijkstra; for heap phases
//                no n*n term is required.
// Auxiliary memory for both algorithms: O(n), excluding the matroid oracles.
// Weighted exchange edges are generated lazily during Dijkstra, not stored.
//
// Note: all arithmetic using integral weights is exact unless the sum of all
// input magnitudes itself exceeds signed int128 (not possible at practical
// n with <=64-bit weights). The floating variant is approximate.

namespace matroid_intersection_detail {
inline std::vector<int> chosen_indices(const std::vector<char>& S) {
    std::vector<int> result;
    result.reserve(S.size());
    for (int i = 0; i < (int)S.size(); ++i) if (S[i]) result.push_back(i);
    return result;
}

template<class M1, class M2>
std::vector<int> cardinality(int n, M1& m1, M2& m2) {
    assert(n >= 0);
    std::vector<char> S(n, 0);
    std::vector<int> prev(n), que(n), in, out;
    std::vector<char> is_target(n);
    in.reserve(n); out.reserve(n);
    while (true) {
        m1.build(S); m2.build(S);
        in.clear(); out.clear();
        for (int x = 0; x < n; ++x) (S[x] ? in : out).push_back(x);
        std::fill(prev.begin(), prev.end(), -2);
        int tail = 0;
        int direct = -1;
        for (int x : out) {
            const bool source = m1.can_add(x);
            const bool target = m2.can_add(x);
            is_target[x] = target;
            if (source) {
                if (target) { direct = x; break; }
                prev[x] = -1;
                que[tail++] = x;
            }
        }
        // A one-element augmentation is always valid, with no exchange graph.
        if (direct >= 0) { S[direct] = 1; continue; }
        int head = 0, last = -1;
        while (head < tail) {
            const int v = que[head++];
            if (!S[v] && is_target[v]) { last = v; break; }
            if (S[v]) {
                for (int x : out) {
                    if (prev[x] == -2 && m1.can_exchange(x, v)) {
                        prev[x] = v;
                        que[tail++] = x;
                    }
                }
            } else {
                for (int y : in) {
                    if (prev[y] == -2 && m2.can_exchange(v, y)) {
                        prev[y] = v;
                        que[tail++] = y;
                    }
                }
            }
        }
        if (last < 0) return chosen_indices(S);
        for (int x = last; x != -1; x = prev[x]) S[x] ^= 1;
    }
}

template<class Cost, class T, class M1, class M2>
std::optional<std::vector<int>> weighted_impl(const std::vector<T>& weight, M1& m1, M2& m2,
                          bool cardinality_first) {
    static_assert(std::is_arithmetic_v<T> &&
                  (!std::is_integral_v<T> || sizeof(T) <= 8),
                  "Weight type must be integral up to 64 bits or floating point");
    const int n = (int)weight.size(), s = n, t = n + 1, nn = n + 2;
    std::vector<char> S(n, 0);
    std::vector<Cost> potential(nn, 0), dist(nn), enter(n);
    for (int x = 0; x < n; ++x) enter[x] = -Cost(weight[x]);
    std::vector<int> ins, outs, prev(nn), hop(nn);
    std::vector<char> add1(n), add2(n);
    std::vector<char> reached(nn), fixed(nn);
    ins.reserve(n); outs.reserve(n);
    Cost current_weight = 0, best_weight = 0;
    std::vector<int> best;

    __int128_t max_abs_weight = 0;
    if constexpr (std::is_same_v<Cost, long long>) {
        for (T w : weight) {
            __int128_t val = (__int128_t)w;
            if (val < 0) val = -val;
            max_abs_weight = std::max(max_abs_weight, val);
        }
    }
    while (true) {
        // Check a worst-case upper bound for every reduced cost, temporary
        // distance, and potential update of the next augmentation.
        // If even this conservative bound could overflow int64, restart with
        // signed int128. The wrapper checks the input before the first negation.
        if constexpr (std::is_same_v<Cost, long long>) {
            __int128_t pmax = 0;
            for (Cost p : potential) {
                __int128_t val = (__int128_t)p;
                if (val < 0) val = -val;
                pmax = std::max(pmax, val);
            }
            const __int128_t pp = pmax + max_abs_weight;
            const __int128_t edge_bound = max_abs_weight + 2 * pp;
            const __int128_t safety = pp + (__int128_t(nn) + 3) * edge_bound
                                      + (__int128_t(nn) + 2) * max_abs_weight;
            if (safety > (__int128_t)std::numeric_limits<long long>::max() / 4)
                return std::nullopt;
        }
        m1.build(S); m2.build(S);
        ins.clear(); outs.clear();
        for (int i = 0; i < n; ++i) (S[i] ? ins : outs).push_back(i);
        // Safe weighted greedy shortcut: if a maximum-weight unused element
        // is addable to both matroids, simply add it. This yields an optimal
        // size-(k+1) intersection from the optimal size-k intersection.
#ifndef MATROID_INTERSECTION_DISABLE_DIRECT_SHORTCUT
        if (!outs.empty()) {
            int max_x = outs[0];
            for (int x : outs) if (weight[x] > weight[max_x]) max_x = x;
            // If all strictly positive elements are already selected, no
            // feasible set can have strictly greater total weight.
            if (!cardinality_first && weight[max_x] <= T(0)) break;
            int direct = -1;
            for (int x : outs) {
                if (weight[x] == weight[max_x] && m1.can_add(x) && m2.can_add(x)) {
                    direct = x; break;
                }
            }
            if (direct >= 0) {
                const Cost old_l = enter[direct];
                potential[direct] -= old_l;
                S[direct] = 1;
                enter[direct] = -old_l;
                current_weight -= old_l;
                if (cardinality_first || current_weight >= best_weight) {
                    best = chosen_indices(S); best_weight = current_weight;
                }
                if (!cardinality_first && std::is_integral_v<T> && current_weight < best_weight)
                    break;
                continue;
            }
        }
#endif
        // Check whether either side is already saturated before asking the
        // potentially expensive exchange oracle about any pairs.
        bool has_source = false, has_sink = false;
        Cost ps = 0, pt = 0;
        for (int x : outs) {
            add1[x] = m1.can_add(x);
            add2[x] = m2.can_add(x);
            if (add1[x]) {
                if (!has_sink || potential[x] < pt) pt = potential[x];
                has_sink = true;
            }
            if (add2[x]) {
                const Cost p = potential[x] - enter[x];
                if (!has_source || p > ps) ps = p;
                has_source = true;
            }
        }
        if (!has_source || !has_sink) break;
        potential[s] = ps; potential[t] = pt;
        std::fill(reached.begin(), reached.end(), false);
        std::fill(fixed.begin(), fixed.end(), false);
        std::fill(prev.begin(), prev.end(), -1);
        reached[s] = true; dist[s] = 0; hop[s] = 0;
        auto visit_edges = [&](int u, auto&& fn) {
            if (u == s) {
                for (int x : outs) if (add2[x]) fn(x);
            } else if (S[u]) {
                for (int x : outs)
                    if (!add2[x] && m2.can_exchange(x, u)) fn(x);
            } else if (add1[u]) {
                fn(t);
            } else {
                for (int y : ins)
                    if (m1.can_exchange(u, y)) fn(y);
            }
        };
        bool use_dense = nn <= 30;
#ifdef MATROID_INTERSECTION_FORCE_DENSE
        use_dense = true;
#endif
#ifdef MATROID_INTERSECTION_FORCE_HEAP
        use_dense = false;
#endif
        auto relax = [&](int u, int v) {
            if (fixed[v]) return false;
            const Cost l = v == t ? Cost(0) : enter[v];
            Cost rc = l + potential[u] - potential[v];
            if constexpr (std::is_floating_point_v<T>) {
                const Cost mag = std::max({Cost(1), std::fabs(l),
                                          std::fabs(potential[u]), std::fabs(potential[v])});
                assert(rc >= -Cost(64)*std::numeric_limits<Cost>::epsilon()*mag);
                if (rc < 0) rc = 0;
            } else assert(rc >= 0);
            const Cost nd = dist[u] + rc;
            if (!reached[v] || nd < dist[v] || (nd == dist[v] && hop[u]+1 < hop[v])) {
                reached[v] = true; dist[v] = nd; hop[v] = hop[u]+1; prev[v] = u;
                return true;
            }
            return false;
        };
        if (use_dense) {
            for (int z = 0; z < nn; ++z) {
                int u = -1;
                for (int v = 0; v < nn; ++v) {
                    if (fixed[v] || !reached[v]) continue;
                    if (u == -1 || dist[v] < dist[u] ||
                        (dist[v] == dist[u] && hop[v] < hop[u])) u = v;
                }
                if (u < 0 || u == t) break;
                fixed[u] = true;
                visit_edges(u,[&](int v){relax(u,v);});
            }
        } else {
            // Indexed binary heap: decrease-key without duplicate entries.
            std::vector<int> heap, position(nn, -1);
            heap.reserve(nn);
            auto less = [&](int a, int b) {
                return dist[a] < dist[b] ||
                       (dist[a] == dist[b] && hop[a] < hop[b]);
            };
            auto sift_up = [&](int i) {
                while (i > 0) {
                    const int p = (i-1)/2;
                    if (!less(heap[i], heap[p])) break;
                    std::swap(heap[i], heap[p]);
                    position[heap[i]] = i; position[heap[p]] = p;
                    i = p;
                }
            };
            auto update = [&](int v) {
                int p = position[v];
                if (p < 0) {
                    p = (int)heap.size();
                    heap.push_back(v); position[v] = p;
                }
                sift_up(p);
            };
            auto pop = [&]() {
                const int v = heap[0], last = heap.back();
                position[v] = -1;
                heap.pop_back();
                if (!heap.empty()) {
                    heap[0] = last; position[last] = 0;
                    int i = 0;
                    for (;;) {
                        int j = 2*i+1;
                        if (j >= (int)heap.size()) break;
                        if (j+1 < (int)heap.size() && less(heap[j+1],heap[j])) ++j;
                        if (!less(heap[j],heap[i])) break;
                        std::swap(heap[j],heap[i]);
                        position[heap[j]] = j; position[heap[i]] = i;
                        i = j;
                    }
                }
                return v;
            };
            update(s);
            while (!heap.empty()) {
                const int u = pop();
                if (u == t) break;
                fixed[u] = true;
                visit_edges(u,[&](int v){
                    if (relax(u,v)) update(v);
                });
            }
        }
        if (!reached[t]) break;
        const Cost dt = dist[t];
        for (int v = 0; v < nn; ++v) {
            const Cost d = reached[v] && dist[v] < dt ? dist[v] : dt;
            potential[v] += d;
        }
        for (int v = prev[t]; v != s; v = prev[v]) {
            assert(v >= 0);
            const Cost old_l = enter[v];
            potential[v] -= old_l;
            S[v] ^= 1;
            enter[v] = -old_l;
            current_weight -= old_l;
        }
        if (cardinality_first || current_weight >= best_weight) {
            best = chosen_indices(S); best_weight = current_weight;
        }
        // The optimal value for each cardinality is concave in cardinality.
        // After the first strictly negative marginal gain, no later size can
        // improve maximum total weight. Stop early (integer arithmetic only).
        if (!cardinality_first && std::is_integral_v<T> && current_weight < best_weight)
            break;
    }
    return best;
}
// Use int64 for the common case, with a *checked* fallback to exact int128.
// This avoids a static, overly conservative input-weight restriction.
template<class T, class M1, class M2>
std::vector<int> weighted(const std::vector<T>& w, M1& m1, M2& m2,
                          bool cardinality_first) {
    if constexpr (std::is_integral_v<T> && sizeof(T) <= 8) {
        __int128_t max_w = 0;
        for (T x : w) {
            __int128_t y = (__int128_t)x;
            if (y < 0) y = -y;
            max_w = std::max(max_w, y);
        }
        // Initial guard is required before narrowing and negating the weight.
        const __int128_t initial_bound = max_w * (4 * (__int128_t(w.size()) + 2) + 12);
        if (initial_bound <= (__int128_t)std::numeric_limits<long long>::max() / 4) {
            if (auto ans = weighted_impl<long long>(w, m1, m2, cardinality_first))
                return std::move(*ans);
        }
        return *weighted_impl<__int128_t>(w,m1,m2,cardinality_first);
    } else {
        return *weighted_impl<long double>(w,m1,m2,cardinality_first);
    }
}
}  // namespace matroid_intersection_detail

template<class M1, class M2>
std::vector<int> matroid_intersection_max_cardinality(int n, M1& m1, M2& m2) {
    return matroid_intersection_detail::cardinality(n,m1,m2);
}
template<class T, class M1, class M2>
std::vector<int> matroid_intersection_max_weight(const std::vector<T>& w, M1& m1, M2& m2) {
    // Constant weights reduce to unweighted intersection (or the empty set).
    if (!w.empty() && std::all_of(w.begin()+1, w.end(),
                                   [&](const T& x) { return x == w[0]; })) {
        if (w[0] < T(0)) return {};
        return matroid_intersection_detail::cardinality((int)w.size(),m1,m2);
    }
    return matroid_intersection_detail::weighted(w,m1,m2,false);
}
template<class T, class M1, class M2>
std::vector<int> matroid_intersection_max_cardinality_weight(const std::vector<T>& w, M1& m1, M2& m2) {
    if (!w.empty() && std::all_of(w.begin()+1, w.end(),
                                   [&](const T& x) { return x == w[0]; }))
        return matroid_intersection_detail::cardinality((int)w.size(),m1,m2);
    return matroid_intersection_detail::weighted(w,m1,m2,true);
}
template<class M1, class M2>
std::vector<int> matroid_intersection(int n, M1& m1, M2& m2) {
    return matroid_intersection_max_cardinality(n,m1,m2);
}
