#pragma once

// Partial Interval Hitting Set (maximum weighted interval stabbing).
// Dependency: Mincostflow_with_Negative.hpp from mylib (which uses ACL).
// The dependency uses unqualified std names, hence the using-directive below.
#include <bits/stdc++.h>
#include <atcoder/mincostflow>
using namespace std;
#include "Mincostflow_with_Negative.hpp"

// Choose at most k points on a line to maximize the total weight of
// closed intervals [l,r] containing at least one chosen point.
// Nonnegative integral weights only. Duplicate intervals are supported.
// solve(k, true) also returns an optimal set of points.
// Coordinates need only support total ordering (normally long long).
//
// Graph: candidate points = distinct right endpoints p[0..m-1].
// s->t: (W+2,0), s->0: (W+1,k), m->t: (W+1,0),
// j->j+1: (W+1,0), interval [l,r]: (L,R,w,-1),
// where L=lower_bound(p,l), R=upper_bound(p,r).
// Send W+1 units; answer = W+minimum_cost.
//
// Restoration: compute a feasible potential d of the final residual graph.
// The optimal selected point set is {p[j] : d[j]>d[j+1]}.
// The W+1 flow and W+2 dummy capacity guarantee an unsaturated dummy edge;
// W+1 capacities on all other noninterval edges keep their forward residuals.
//
// Time: O(n log n) preprocessing, plus MinCostFlowWithNegative::flow
// (worst-case O(W*(n+1)*log(n+1)) in its integer-flow bound, with an
// initial SPFA). Restoration: one residual-graph SPFA (worst O(n^2)).
// Memory: O(n). W = sum of positive weights, and the graph has O(n) edges.

template<class Coord = long long>
struct PartialIntervalHittingSet {
    struct Result {
        long long value = 0;
        vector<Coord> points; // empty when restore=false
    };

private:
    struct Interval {
        Coord l, r;
        long long w;
    };
    vector<Interval> a_;

public:
    void add(Coord l, Coord r, long long w = 1) {
        assert(!(r < l));
        assert(w >= 0);
        a_.push_back({l, r, w});
    }

    int size() const { return (int)a_.size(); }

    // Minimum number of points hitting ALL added intervals (including w=0).
    // Right-endpoint greedy; NO min-cost flow. O(n log n) time, O(n) space.
    // result.value is the minimum number of points, unlike solve(), where
    // it denotes the maximum covered weight.
    Result min_points_to_cover_all(bool restore = false) const {
        Result result;
        if (a_.empty()) return result;

        auto sorted = a_;
        sort(sorted.begin(), sorted.end(), [](const Interval& x, const Interval& y) {
            if (x.r != y.r) return x.r < y.r;
            return x.l < y.l;
        });

        optional<Coord> last;
        for (const auto& it : sorted) {
            // The previously chosen right endpoint is <= it.r.
            if (!last || *last < it.l) {
                last = it.r;
                ++result.value;
                if (restore) result.points.push_back(it.r);
            }
        }
        return result;
    }

    Result solve(int k, bool restore = false) const {
        assert(k >= 0);
        Result result;
        if (k == 0 || a_.empty()) return result;

        vector<Interval> a;
        a.reserve(a_.size());
        long long W = 0;
        for (auto it : a_) {
            if (it.w == 0) continue;
            // Keep all capacities, accumulated costs, and potentials well
            // within signed 64-bit range for the current MCF dependency.
            assert(it.w <= (LLONG_MAX / 4) - W);
            W += it.w;
            a.push_back(it);
        }
        if (W == 0) return result;

        sort(a.begin(), a.end(), [](const Interval &x, const Interval &y) {
            if (x.r != y.r) return x.r < y.r;
            return x.l < y.l;
        });

        // Fast path: the usual right-endpoint greedy gives the minimum
        // number of points needed to hit ALL positive-weight intervals.
        vector<Coord> all_covered_points;
        for (const auto &it : a) {
            if (all_covered_points.empty() || all_covered_points.back() < it.l)
                all_covered_points.push_back(it.r);
        }
        if (all_covered_points.size() <= (size_t)k) {
            result.value = W;
            if (restore) result.points = std::move(all_covered_points);
            return result;
        }

        vector<Coord> xs;
        xs.reserve(a.size());
        for (auto &it : a) xs.push_back(it.r);
        xs.erase(unique(xs.begin(), xs.end()), xs.end());
        int m = (int)xs.size();
        int s = m + 1, t = m + 2, V = m + 3;

        // The graph encodes the DUAL problem. Flow does not directly
        // represent chosen points; the points come from residual potentials.
        MinCostFlowWithNegative<long long, long long> g(V);
        long long cap = W + 1;
        g.add_edge(s, t, W + 2, 0);  // dummy path, never saturated
        g.add_edge(s, 0, cap, k);
        g.add_edge(m, t, cap, 0);
        for (int j = 0; j < m; ++j) g.add_edge(j, j + 1, cap, 0);
        for (auto &it : a) {
            int L = lower_bound(xs.begin(), xs.end(), it.l) - xs.begin();
            int R = upper_bound(xs.begin(), xs.end(), it.r) - xs.begin();
            assert(L < R);
            g.add_edge(L, R, it.w, -1);
        }

        auto [flow, cost] = g.flow(s, t, cap);
        assert(flow == cap);
        result.value = W + cost;
        assert(0 <= result.value && result.value <= W);
        if (!restore) return result;

        // Recover any feasible residual-graph potential d satisfying
        // d[v] <= d[u] + cost(u,v) on all residual arcs.
        // Super-source with zero-cost arcs to every vertex is represented
        // by initializing every distance to zero and enqueuing all vertices.
        // No negative cycles exist after an optimal min-cost flow.
        vector<vector<pair<int, long long>>> adj(V);
        for (const auto &e : g.edges()) {
            if (e.flow < e.cap) adj[e.from].push_back({e.to, e.cost});
            if (e.flow > 0) adj[e.to].push_back({e.from, -e.cost});
        }

        vector<long long> d(V, 0);
        vector<unsigned char> inq(V, true);
        deque<int> q;
        for (int v = 0; v < V; ++v) q.push_back(v);
        while (!q.empty()) {
            int u = q.front();
            q.pop_front();
            inq[u] = false;
            for (auto [v, c] : adj[u]) {
                if (d[v] <= d[u] + c) continue;
                d[v] = d[u] + c;
                if (!inq[v]) {
                    inq[v] = true;
                    if (!q.empty() && d[v] < d[q.front()]) q.push_front(v);
                    else q.push_back(v);
                }
            }
        }

        for (int j = 0; j < m; ++j) {
            assert(d[j] >= d[j + 1]);
            if (d[j] > d[j + 1]) result.points.push_back(xs[j]);
        }
        assert(result.points.size() <= (size_t)k);
        return result;
    }
};
