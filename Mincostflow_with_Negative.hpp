#pragma once
/** https://atcoder.jp/contests/abc247/submissions/73798754
 ** https://judge.yosupo.jp/submission/357342
 *
 * MinCostFlowWithNegative
 * 負の辺があるグラフで最小費用流を解く O(FE logE) Fは流量
 * * 使い方:
 * g.add_edge(from, to, cap, cost)
 * g.flow(s, t, limit) -> {流量, 最小費用}
 * g.slope(s, t, limit) -> {{流量, 最小費用}, ...} の折れ線
 * g.get_cost(slope_data, f) -> 流量 f のときの最小費用 (O(log V))
 * g.edges() -> 全辺の状態 {from, to, cap, flow, cost} を取得
 * g.get_edge(i) -> i番目に追加された辺の状態を取得
 */
template <typename Cap, typename Cost>
struct MinCostFlowWithNegative {
    struct Edge {
        int from, to;
        Cap cap, flow;
        Cost cost;
    };

    int n;
    vector<Edge> _edges;
    const Cost INF_COST = numeric_limits<Cost>::max() / 2;

    MinCostFlowWithNegative(int n) : n(n) {}

    void add_edge(int from, int to, Cap cap, Cost cost) {
        assert(cap >= 0);
        if (cap == 0) return;
        _edges.push_back({from, to, cap, 0, cost});
    }

    pair<Cap, Cost> flow(int s, int t, Cap flow_limit = numeric_limits<Cap>::max()) {
        // SPFAで初期ポテンシャルを計算
        vector<Cost> h(n, INF_COST);
        vector<int> cnt(n, 0);
        vector<bool> in_queue(n, false);
        queue<int> q;

        vector<vector<pair<int, Cost>>> g(n);
        for (const auto& e : _edges) {
            if (e.cap > 0) g[e.from].push_back({e.to, e.cost});
        }

        h[s] = 0;
        q.push(s);
        in_queue[s] = true;

        while (!q.empty()) {
            int u = q.front(); q.pop();
            in_queue[u] = false;

            for (auto& edge : g[u]) {
                int v = edge.first;
                Cost weight = edge.second;
                if (h[v] > h[u] + weight) {
                    h[v] = h[u] + weight;
                    if (!in_queue[v]) {
                        q.push(v);
                        in_queue[v] = true;
                        if (++cnt[v] > n) return {0, -INF_COST}; // 負閉路検出
                    }
                }
            }
        }

        if (h[t] == INF_COST) return {0, 0};

        // ポテンシャルを用いて非負コスト化し、mcf_graphに辺を張る
        atcoder::mcf_graph<Cap, Cost> mcf(n);
        vector<int> mcf_edge_idx(_edges.size(), -1);
        for (size_t i = 0; i < _edges.size(); ++i) {
            if (h[_edges[i].from] == INF_COST) continue;
            Cost new_cost = _edges[i].cost + h[_edges[i].from] - h[_edges[i].to];
            mcf_edge_idx[i] = mcf.add_edge(_edges[i].from, _edges[i].to, _edges[i].cap, new_cost);
        }

        auto result = mcf.flow(s, t, flow_limit);

        // 各辺に流れた量(flow)を元の辺情報に反映
        for (size_t i = 0; i < _edges.size(); ++i) {
            if (mcf_edge_idx[i] != -1) {
                _edges[i].flow = mcf.get_edge(mcf_edge_idx[i]).flow;
            } else {
                _edges[i].flow = 0;
            }
        }

        // 補正したコストから真のコストを復元
        Cost true_cost = result.second + (h[t] - h[s]) * result.first;
        return {result.first, true_cost};
    }

    vector<pair<Cap, Cost>> slope(int s, int t, Cap flow_limit = numeric_limits<Cap>::max()) {
        // SPFAで初期ポテンシャルを計算
        vector<Cost> h(n, INF_COST);
        vector<int> cnt(n, 0);
        vector<bool> in_queue(n, false);
        queue<int> q;

        vector<vector<pair<int, Cost>>> g(n);
        for (const auto& e : _edges) {
            if (e.cap > 0) g[e.from].push_back({e.to, e.cost});
        }

        h[s] = 0;
        q.push(s);
        in_queue[s] = true;

        while (!q.empty()) {
            int u = q.front(); q.pop();
            in_queue[u] = false;

            for (auto& edge : g[u]) {
                int v = edge.first;
                Cost weight = edge.second;
                if (h[v] > h[u] + weight) {
                    h[v] = h[u] + weight;
                    if (!in_queue[v]) {
                        q.push(v);
                        in_queue[v] = true;
                        if (++cnt[v] > n) return {{0, -INF_COST}}; // 負閉路検出
                    }
                }
            }
        }

        if (h[t] == INF_COST) return {{0, 0}};

        // ポテンシャルを用いて非負コスト化し、mcf_graphに辺を張る
        atcoder::mcf_graph<Cap, Cost> mcf(n);
        vector<int> mcf_edge_idx(_edges.size(), -1);
        for (size_t i = 0; i < _edges.size(); ++i) {
            if (h[_edges[i].from] == INF_COST) continue;
            Cost new_cost = _edges[i].cost + h[_edges[i].from] - h[_edges[i].to];
            mcf_edge_idx[i] = mcf.add_edge(_edges[i].from, _edges[i].to, _edges[i].cap, new_cost);
        }

        auto result = mcf.slope(s, t, flow_limit);

        // 各辺に流れた量(flow)を元の辺情報に反映
        for (size_t i = 0; i < _edges.size(); ++i) {
            if (mcf_edge_idx[i] != -1) {
                _edges[i].flow = mcf.get_edge(mcf_edge_idx[i]).flow;
            } else {
                _edges[i].flow = 0;
            }
        }

        // 折れ線の各頂点のコストを真のコストに復元
        vector<pair<Cap, Cost>> true_slope;
        for (auto& p : result) {
            Cap f = p.first;
            Cost mod_cost = p.second;
            Cost true_cost = mod_cost + (h[t] - h[s]) * f;
            true_slope.push_back({f, true_cost});
        }

        return true_slope;
    }

    Cost get_cost(const vector<pair<Cap, Cost>>& slope_data, Cap f) {
        if (slope_data.empty()) return INF_COST; 
        if (f <= 0) return 0;                    
        if (f > slope_data.back().first) return INF_COST;

        auto it = lower_bound(slope_data.begin(), slope_data.end(), f,
            [](const pair<Cap, Cost>& p, Cap val) {
                return p.first < val;
            });

        if (it->first == f) return it->second;

        auto it_prev = prev(it);
        Cap f0 = it_prev->first;
        Cost c0 = it_prev->second;
        Cap f1 = it->first;
        Cost c1 = it->second;

        Cost unit_cost = (c1 - c0) / (f1 - f0);
        return c0 + unit_cost * (f - f0);
    }

    vector<Edge> edges() {
        return _edges;
    }

    Edge get_edge(int i) {
        return _edges[i];
    }
};