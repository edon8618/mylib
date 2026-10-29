#pragma once

#include <algorithm>
#include <utility>
#include <vector>

/** diameter  https://judge.yosupo.jp/submission/365865
 *  lca  https://atcoder.jp/contests/abc014/submissions/74903786
 *  decompose  https://codeforces.com/problemset/problem/321/C
 */
template <typename T>
struct Tree {
    // 内部に隠蔽された Edge 構造体（Tの型をそのまま使える）
    struct Edge {
        int to;
        T weight;
    };

    int n;
    std::vector<std::vector<Edge>> graph;

    Tree(int n) : n(n), graph(n) {}

    void add_edge(int u, int v, T w) {
        graph[u].push_back({v, w});
        graph[v].push_back({u, w});
    }

    void add_edge(int u, int v) {
        graph[u].push_back({v, 1});
        graph[v].push_back({u, 1});
    }


    // diameter
    void dfs(int current, int parent, T depth, std::vector<T>& dist, std::vector<int>& parents) {
        dist[current] = depth;
        parents[current] = parent;
        for (const auto& edge : graph[current]) {
            if (edge.to != parent) {
                dfs(edge.to, current, depth + edge.weight, dist, parents);
            }
        }
    }

    std::pair<T, std::vector<int>> get_diameter_and_path() {
        if (n == 0) return {0, {}};

        std::vector<T> dist(n);
        std::vector<int> parents(n);

        dfs(0, -1, 0, dist, parents);
        int u = 0;
        for (int i = 0; i < n; ++i) if (dist[i] > dist[u]) u = i;

        dfs(u, -1, 0, dist, parents);
        int v = 0;
        for (int i = 0; i < n; ++i) if (dist[i] > dist[v]) v = i;

        std::vector<int> path;
        int current = v;
        while (current != -1) {
            path.push_back(current);
            current = parents[current];
        }
        std::reverse(path.begin(), path.end());

        return {dist[v], path};
    }

    // LCA
    std::vector<int> tour;           // 訪れた頂点の順序
    std::vector<int> first_visit;    // 各頂点がツアーで最初に現れるインデックス
    std::vector<int> depth;          // 各頂点の深さ
    std::vector<T> dist_from_root;   // 根からの距離（和）
    std::vector<T> xor_from_root;    // 根からの距離（xor）
    std::vector<std::vector<int>> st;     // スパーステーブル (tourのインデックスを保持)

    // DFSでオイラーツアーを構築
    void dfs_euler(int v, int p, int d, T sum_dist, T xor_dist) {
        first_visit[v] = tour.size();
        tour.push_back(v);
        depth[v] = d;
        dist_from_root[v] = sum_dist;
        xor_from_root[v] = xor_dist;

        for (const auto& edge : graph[v]) {
            if (edge.to != p) {
                dfs_euler(edge.to, v, d + 1, sum_dist + edge.weight, xor_dist ^ edge.weight);
                tour.push_back(v); // 戻ってきたときも記録
            }
        }
    }

    void build_lca(int root = 0) {
        first_visit.assign(n, -1);
        depth.assign(n, 0);
        dist_from_root.assign(n, 0);
        xor_from_root.assign(n, 0);
        tour.reserve(n * 2);
        
        dfs_euler(root, -1, 0, 0, 0);

        // スパーステーブルの構築
        int m = tour.size();
        int log_m = 0;
        while ((1 << (log_m + 1)) <= m) log_m++;
        st.assign(log_m + 1, std::vector<int>(m));

        for (int i = 0; i < m; i++) st[0][i] = tour[i];

        for (int k = 0; k < log_m; k++) {
            for (int i = 0; i + (1 << (k + 1)) <= m; i++) {
                int left = st[k][i];
                int right = st[k][i + (1 << k)];
                // 深さが小さい方の頂点を選択
                st[k + 1][i] = (depth[left] < depth[right] ? left : right);
            }
        }
    }

    int get_lca(int u, int v) {
        int l = first_visit[u];
        int r = first_visit[v];
        if (l > r) std::swap(l, r);
        int k = 31 - __builtin_clz(r - l + 1); // log2を高速計算
        int left = st[k][l];
        int right = st[k][r - (1 << k) + 1];
        return (depth[left] < depth[right] ? left : right);
    }

    // 通常の距離 (和)
    T get_dist(int u, int v) {
        return dist_from_root[u] + dist_from_root[v] - 2 * dist_from_root[get_lca(u, v)];
    }

    // XOR距離 (LCAすら不要！)
    T get_xor_dist(int u, int v) {
        return xor_from_root[u] ^ xor_from_root[v];
    }

    // --- 重心分解用のメンバ変数 ---
    std::vector<bool> cd_removed;         // 重心として採用済み（削除扱い）フラグ
    std::vector<int> cd_sz;               // 部分木のサイズ
    std::vector<int> centroid_parent;     // 重心木における親（根の場合は -1）

    // ① 部分木のサイズを計算する内部DFS
    void cd_calc_size(int v, int p) {
        cd_sz[v] = 1;
        for (const auto& edge : graph[v]) {
            if (edge.to != p && !cd_removed[edge.to]) {
                cd_calc_size(edge.to, v);
                cd_sz[v] += cd_sz[edge.to];
            }
        }
    }

    // ② 現在の連結成分における重心を探す内部DFS
    int cd_find_centroid(int v, int p, int total_size) {
        for (const auto& edge : graph[v]) {
            if (edge.to != p && !cd_removed[edge.to]) {
                if (cd_sz[edge.to] > total_size / 2) {
                    return cd_find_centroid(edge.to, v, total_size);
                }
            }
        }
        return v; // どこに行っても半分以下なら、現在地が重心
    }

    // ③ 重心分解を再帰的に構築する内部関数
    int build_centroid_tree(int v, int p = -1) {
        cd_calc_size(v, -1);
        int centroid = cd_find_centroid(v, -1, cd_sz[v]);
        
        centroid_parent[centroid] = p;
        cd_removed[centroid] = true;
        
        for (const auto& edge : graph[centroid]) {
            if (!cd_removed[edge.to]) {
                build_centroid_tree(edge.to, centroid);
            }
        }
        return centroid;
    }

    // ④ ユーザーが呼ぶ重心分解の実行関数
    void decompose() {
        // 既存のコンストラクタを汚さないよう、ここで配列を初期化・確保する
        cd_removed.assign(n, false);
        cd_sz.assign(n, 0);
        centroid_parent.assign(n, -1);
        
        // 頂点0を起点に分解スタート
        build_centroid_tree(0, -1);
    }
};