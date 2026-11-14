#pragma once
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

// Dynamic forest with vertex values and directional edge values.
// op(a,b): visit a, then b. See docs/LinkCutTree.txt.
template<class S, S (*op)(S, S), S (*e)()>
struct LinkCutTree {
private:
    struct Node {
        int l = -1, r = -1, p = -1;
        int vertices = 1, size = 1, id = -1;
        bool rev = false, backward = false;
        S val = e(), back = e(), prod = e(), rprod = e();
    };
    struct Edge { int u, v, node; };
    int n;
    int plain_edges = 0;
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::vector<int> free_nodes, path;
    std::unordered_map<std::uint64_t, int> edge_ids;

    static std::uint64_t key(int u, int v) {
        if (u > v) std::swap(u, v);
        return (std::uint64_t(std::uint32_t(u)) << 32) | std::uint32_t(v);
    }
    bool is_root(int x) const {
        int p = nodes[x].p;
        return p == -1 || (nodes[p].l != x && nodes[p].r != x);
    }
    void pull(int x) {
        auto& t = nodes[x];
        t.vertices = (x < n); t.size = 1;
        t.prod = t.val; t.rprod = t.back;
        if (t.l != -1) {
            const auto& l = nodes[t.l];
            t.vertices += l.vertices; t.size += l.size;
            t.prod = op(l.prod, t.prod);
            t.rprod = op(t.rprod, l.rprod);
        }
        if (t.r != -1) {
            const auto& r = nodes[t.r];
            t.vertices += r.vertices; t.size += r.size;
            t.prod = op(t.prod, r.prod); t.rprod = op(r.rprod, t.rprod);
        }
    }
    void toggle(int x) {
        auto& t = nodes[x];
        std::swap(t.l, t.r); std::swap(t.prod, t.rprod);
        if (x >= n) { std::swap(t.val, t.back); t.backward = !t.backward; }
        t.rev = !t.rev;
    }
    void push(int x) {
        auto& t = nodes[x];
        if (!t.rev) return;
        if (t.l != -1) toggle(t.l);
        if (t.r != -1) toggle(t.r);
        t.rev = false;
    }
    void rotate(int x) {
        int p = nodes[x].p, g = nodes[p].p;
        if (!is_root(p)) {
            if (nodes[g].l == p) nodes[g].l = x;
            else nodes[g].r = x;
        }
        nodes[x].p = g;
        if (nodes[p].l == x) {
            int b = nodes[x].r;
            nodes[p].l = b; if (b != -1) nodes[b].p = p;
            nodes[x].r = p;
        } else {
            int b = nodes[x].l;
            nodes[p].r = b; if (b != -1) nodes[b].p = p;
            nodes[x].l = p;
        }
        nodes[p].p = x; pull(p); pull(x);
    }
    void splay(int x) {
        path.clear();
        for (int y = x;; y = nodes[y].p) {
            path.push_back(y); if (is_root(y)) break;
        }
        for (auto it = path.rbegin(); it != path.rend(); ++it) push(*it);
        while (!is_root(x)) {
            int p = nodes[x].p;
            if (!is_root(p)) {
                int g = nodes[p].p;
                rotate((nodes[g].l == p) == (nodes[p].l == x) ? p : x);
            }
            rotate(x);
        }
    }
    int access(int x) {
        int last = -1;
        for (int y = x; y != -1; y = nodes[y].p) {
            splay(y); nodes[y].r = last; pull(y); last = y;
        }
        splay(x); return last;
    }
    void make_root(int x) { access(x); toggle(x); }
    int root(int x) {
        access(x); push(x);
        while (nodes[x].l != -1) { x = nodes[x].l; push(x); }
        splay(x); return x;
    }
    void cut_nodes(int u, int v) {
        make_root(u); access(v);
        assert(nodes[v].l == u && nodes[u].r == -1);
        nodes[v].l = -1; nodes[u].p = -1; pull(v);
    }
    int add_edge(int u, int v, S f, S g, bool identity, bool check) {
        assert(0 <= u && u < n && 0 <= v && v < n);
        make_root(u);
        // Keep the old link's no-op on an already connected pair.
        if (check && root(v) == u) return -1;
        int id = int(edges.size());
        edge_ids.emplace(key(u, v), id);
        if (identity) {
            // Keep the vertex-only topology compact until an edge gets a value.
            edges.push_back({u, v, -2}); ++plain_edges;
            nodes[u].p = v;
            return id;
        }
        int x = allocate_node();
        edges.push_back({u, v, x});
        auto& t = nodes[x];
        t.id = id; t.vertices = 0;
        // Represented root remains on v's side: traversal is v -> edge -> u.
        t.val = g; t.back = f; t.prod = g; t.rprod = f; t.backward = true;
        nodes[u].p = x; t.p = v;
        return id;
    }
    int allocate_node() {
        if (free_nodes.empty()) {
            int x = int(nodes.size()); nodes.emplace_back(); return x;
        }
        int x = free_nodes.back(); free_nodes.pop_back(); nodes[x] = Node{}; return x;
    }
    int materialize(int id) {
        int x = edge_node(id);
        if (x != -2) return x;
        int u = edges[id].u, v = edges[id].v;
        int old_root = root(u);
        cut_nodes(u, v);
        // u is now a represented root. Attach it through the new edge to v.
        x = allocate_node();
        auto& t = nodes[x];
        t.id = id; t.vertices = 0; t.backward = true;
        nodes[u].p = x; t.p = v;
        edges[id].node = x; --plain_edges;
        make_root(old_root);
        return x;
    }
    int edge_node(int id) const {
        assert(0 <= id && id < int(edges.size()) && edges[id].node != -1);
        return edges[id].node;
    }
    int kth(int x, int k, bool vertex) {
        assert(0 <= k && k < (vertex ? nodes[x].vertices : nodes[x].size - nodes[x].vertices));
        for (;;) {
            push(x);
            int l = nodes[x].l;
            int cnt = l == -1 ? 0 : (vertex ? nodes[l].vertices : nodes[l].size - nodes[l].vertices);
            if (k < cnt) { x = l; continue; }
            k -= cnt;
            bool self = ((x < n) == vertex);
            if (self && k == 0) { splay(x); return vertex ? x : nodes[x].id; }
            k -= self; x = nodes[x].r;
        }
    }
public:
    LinkCutTree(int n) : n(n) {
        assert(n >= 0);
        nodes.reserve(n ? std::size_t(n) * 2 - 1 : 0);
        nodes.resize(n);
        edges.reserve(n); edge_ids.reserve(n); path.reserve(64);
    }
    void link(int u, int v) { add_edge(u, v, e(), e(), true, true); }
    int link_edge(int u, int v, S f, S g) { return add_edge(u, v, f, g, false, false); }
    void cut(int u, int v) { int id = get_edge_id(u, v); if (id != -1) cut_edge(id); }
    void cut_edge(int id) {
        int x = edge_node(id), u = edges[id].u, v = edges[id].v;
        if (x == -2) { cut_nodes(u, v); --plain_edges; }
        else { cut_nodes(u, x); cut_nodes(v, x); free_nodes.push_back(x); }
        // The surviving v component was rerooted at v, never at an edge node.
        edge_ids.erase(key(u, v)); edges[id].node = -1;
    }
    int get_edge_id(int u, int v) {
        auto it = edge_ids.find(key(u, v)); return it == edge_ids.end() ? -1 : it->second;
    }
    bool is_connected(int u, int v) { return u == v || root(u) == root(v); }
    int find_root(int v) { return root(v); }
    void set(int v, S x) { splay(v); nodes[v].val = nodes[v].back = x; pull(v); }
    S get(int v) { splay(v); return nodes[v].val; }
    void set_edge(int id, S x) { set_edge(id, x, x); }
    void set_edge(int id, S f, S g) {
        int x = materialize(id); splay(x);
        if (nodes[x].backward) std::swap(f, g);
        nodes[x].val = f; nodes[x].back = g; pull(x);
    }
    void set_edge(int u, int v, S f, S g) {
        int id = get_edge_id(u, v); edge_node(id);
        if (edges[id].u != u) std::swap(f, g);
        set_edge(id, f, g);
    }
    std::pair<S, S> get_edge(int id) {
        int x = edge_node(id);
        if (x == -2) return {e(), e()};
        splay(x);
        if (nodes[x].backward) return {nodes[x].back, nodes[x].val};
        return {nodes[x].val, nodes[x].back};
    }
    std::pair<S, S> get_edge(int u, int v) {
        int id = get_edge_id(u, v); auto ans = get_edge(id);
        if (edges[id].u != u) std::swap(ans.first, ans.second);
        return ans;
    }
    void split(int u, int v) { make_root(u); access(v); }
    S prod(int u, int v) { split(u, v); return nodes[v].prod; }
    std::optional<S> try_prod(int u, int v) {
        make_root(u);
        // u starts as an auxiliary root. If connected and distinct,
        // exposing v puts u below v; otherwise u remains isolated in its splay.
        access(v);
        if (u != v && nodes[u].p == -1) return std::nullopt;
        return nodes[v].prod;
    }
    int vertex_count(int u, int v) { split(u, v); return nodes[v].vertices; }
    int edge_count(int u, int v) { return vertex_count(u, v) - 1; }
    int kth_vertex(int u, int v, int k) { split(u, v); return kth(v, k, true); }
    int kth_edge(int u, int v, int k) {
        if (plain_edges == 0) { split(u, v); return kth(v, k, false); }
        int a = kth_vertex(u, v, k), b = kth_vertex(u, v, k + 1);
        return get_edge_id(a, b);
    }
    S get_prod(int v) { return nodes[v].prod; }
    int get_sz(int v) { return nodes[v].vertices; }
    int get_kth(int v, int k) { return kth(v, k, true); }
};
