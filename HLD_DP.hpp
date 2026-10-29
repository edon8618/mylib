#pragma once
/** https://atcoder.jp/contests/abc269/submissions/74228386
 * 
 * Heavy-Light Tree DP (Static) / HLD + 分割統治による木DPライブラリ
 * 
 * * 【概要】
 * 一直線に長い木（Linear Graph等）で計算量が O(N^2) に壊れる木DPを、
 * HLD とパス上の分割統治を用いて O(N log^2 N) に落とす。
 * 各頂点で FPS を掛け合わせる問題や、巨大なデータ構造をマージする問題に最適。
 * * 【型の定義】
 * - struct S (Data): 部分木に対するDPの答えの型（例: FPS）
 * - struct F (Func): 頂点やパス上の遷移（一次関数や行列）を表す型（例: pair<FPS, FPS>）
 * * 【実装する3つの関数】
 * 1. F mergelight(ll v, vector<S> &ss)
 * - 頂点 v に集まった「Lightな子たちの答え(ss)」をマージし、v 単体の遷移(F)を構築する。
 * - ※ ss のマージが重い場合（FPS等）は、priority_queue でサイズの小さい順に掛けること！
 * - ※ ss が空（葉）の場合の単位元処理を忘れないこと。
 * * 2. F mergepath(F f, F g)
 * - 同じ Heavy Path 上の遷移関数を合成する。
 * - 引数は f が「上側（根側）」、g が「下側（葉側）」。
 * - 合成順序は f(g(X)) になるように実装すること！（非可換な場合に注意）
 * * 3. S eval(F f)
 * - Heavy Path 全体を合成してできた関数 f の一番下に、木構造の終端の値
 * （X = 1 や X = 0 など）を代入し、最終的な部分木の答え(S)を出力する。
 * * 【使い方】
 * HLD hld(n);
 * rep(i, n-1) hld.add_edge(u, v);
 * hld.build(root); // root を指定して構築
 * * HeavyLightDP<S, F> hlddp(hld, mergelight, mergepath, eval);
 * S ans = hlddp.solve(root);
 * =========================================================================
 */
struct HLD {
    int n;
    vector<vector<int>> g;
    // グラフ構造の基本情報
    vector<int> sz, parent, head, heavy_child;
    // 配列マッピング用（セグ木などに載せるときに使う）
    vector<int> in, out, rev; 
    int timer;

    HLD(int n) : n(n), g(n), sz(n, 1), parent(n, -1), head(n), 
                 heavy_child(n, -1), in(n), out(n), rev(n), timer(0) {}

    // 無向辺を追加
    void add_edge(int u, int v) {
        g[u].push_back(v);
        g[v].push_back(u);
    }

    // 構築（根を指定）
    void build(int root = 0) {
        dfs_sz(root, -1);
        head[root] = root;
        dfs_hld(root, -1);
    }

private:
    // 1回目のDFS：部分木のサイズ(sz)と、最も重い子(heavy_child)を決定する
    void dfs_sz(int v, int p) {
        parent[v] = p;
        int max_sub = 0;
        for (int u : g[v]) {
            if (u == p) continue;
            dfs_sz(u, v);
            sz[v] += sz[u];
            if (max_sub < sz[u]) {
                max_sub = sz[u];
                heavy_child[v] = u;
            }
        }
    }

    // 2回目のDFS：headを決定し、行きがけ順(in)・帰りがけ順(out)を振る
    void dfs_hld(int v, int p) {
        in[v] = timer++;       // 頂点vが配列のどのインデックスに対応するか
        rev[in[v]] = v;        // 配列のインデックスから頂点番号を逆引きする用

        // 【重要】必ず Heavy Child を最初に訪問する！
        // これにより、同じ Heavy Path 上の頂点が `in` 配列上で「連続」するようになる
        if (heavy_child[v] != -1) {
            head[heavy_child[v]] = head[v]; // Heavy Childは親と同じパス
            dfs_hld(heavy_child[v], v);
        }

        // 次に Light Child たちを訪問する
        for (int u : g[v]) {
            if (u == p || u == heavy_child[v]) continue;
            head[u] = u; // Light Child は新しいパスの先頭になる
            dfs_hld(u, v);
        }
        
        out[v] = timer; // 部分木クエリ用の抜け出し時間
    }
};


template<class Data, class Func>
struct HeavyLightDP {
    HLD& hld; // HLD自体は重いので、参照で持つのがおすすめ

    // ユーザーが定義して注入する演算たち
    // 第1引数の int は「現在の頂点番号」です
    using MergeLight = function<Func(int, vector<Data>&)>;
    using MergePath  = function<Func(Func, Func)>;
    using Eval       = function<Data(Func)>;
    
    MergeLight merge_light;
    MergePath merge_path;
    Eval eval;

    HeavyLightDP(HLD& hld, MergeLight ml, MergePath mp, Eval ev)
        : hld(hld), merge_light(ml), merge_path(mp), eval(ev) {}

    // 分割統治（D&C）の本体 [l, r)
    Func dc(vector<Func>& v, int l, int r) {
        if (r - l == 1) return std::move(v[l]);
        int m = (l + r) / 2;
        // 上側(top) と 下側(bot) を合成
        return merge_path(dc(v, l, m), dc(v, m, r));
    }

    // DP本体
    Data solve(int root = 0) {
        int n = hld.n;
        vector<Data> dp_data(n);
        vector<Func> dp_func(n);

        // HLDの「行きがけ順(in)」の逆順 = 完璧なトポロジカルソート(葉から根へ)
        for (int i = n - 1; i >= 0; i--) {
            int u = hld.rev[i];

            // 1. Light子の Data を集める
            vector<Data> light_children;
            for (int v : hld.g[u]) {
                if (v == hld.parent[u] || v == hld.heavy_child[u]) continue;
                
                // 子の Data はここでしか使わないので move してメモリと時間を節約！
                light_children.push_back(std::move(dp_data[v]));
            }

            // 2. 頂点 u の Func（一次関数や行列など） を構築
            dp_func[u] = merge_light(u, light_children);

            // 3. u が Heavy Path の先頭(head)なら、パス全体をマージして Data を確定させる
            if (u == hld.head[u]) {
                vector<Func> path;
                int curr = u;
                // パスを下に向かって収集
                while (curr != -1) {
                    path.push_back(std::move(dp_func[curr])); // これも move でヨシ！
                    curr = hld.heavy_child[curr];
                }
                
                // パス全体を分割統治でドカンと合成
                Func path_func = dc(path, 0, path.size());
                
                // 評価関数を通して最終的な Data に変換し、dp_data に保存
                dp_data[u] = eval(std::move(path_func));
            }
        }
        
        // 根(root)は必ずある Heavy Path の head なので、ここに木全体の答えが入っている
        return dp_data[root];
    }
};
