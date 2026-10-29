#pragma once
/**https://judge.yosupo.jp/submission/379073
 * --- Static 2D Segment Tree (事前構築型 2次元セグメント木) ---
 * * [特徴]
 * - 空間計算量 O(N log N) / クエリ O(log^2 N) の高速・省メモリ仕様
 * - 逆元が不要な演算 (Min, Max, GCDなど) の矩形取得が可能
 * * [使い方]
 * // 1. 乗せるモノイドの定義 (ACLと同じインターフェース)
 * using S = long long;
 * S op(S a, S b) { return a + b; }  // 演算 (Max, Min, + など)
 * S e() { return 0; }               // 単位元
 * * // 2. 座標の事前登録と構築 (オフライン)
 * // ※ 構築後に点が追加されることはないため、クエリで登場する座標もすべて入れる,setとaddで使う座標のみ登録すればよい
 * vector<pair<long long, long long>> pts;
 * pts.push_back({x1, y1});
 * pts.push_back({x2, y2});
 * Static2DSegTree<S, op, e> seg(pts);
 * * // 3. クエリ操作
 * seg.add(x, y, val);  // 座標(x, y) の現在の値に val を作用させる ( op(現在値, val) )
 * seg.set(x, y, val);  // 座標(x, y) の値を val で上書きする
 * S v = seg.get(x, y); // 座標(x, y) の現在の値を取得する (未登録座標なら e() が返る)
 * * // 矩形領域 [min_x, max_x) × [min_y, max_y) の演算結果を取得する (半開区間)
 * S ans = seg.prod(min_x, max_x, min_y, max_y);
 */
template <class S, S (*op)(S, S), S (*e)(), class Coord = long long>
struct Static2DSegTree {
  private:
    int n;
    std::vector<Coord> xs;
    std::vector<std::vector<Coord>> ys;
    std::vector<std::vector<S>> d;

    int get_y_idx(int k, Coord y) const {
        return std::lower_bound(ys[k].begin(), ys[k].end(), y) - ys[k].begin();
    }

    S inner_prod(int k, Coord min_y, Coord max_y) const {
        int l = std::lower_bound(ys[k].begin(), ys[k].end(), min_y) - ys[k].begin();
        int r = std::lower_bound(ys[k].begin(), ys[k].end(), max_y) - ys[k].begin();
        l += ys[k].size();
        r += ys[k].size();
        S sml = e(), smr = e();
        while (l < r) {
            if (l & 1) sml = op(sml, d[k][l++]);
            if (r & 1) smr = op(d[k][--r], smr);
            l >>= 1;
            r >>= 1;
        }
        return op(sml, smr);
    }

  public:
    // コンストラクタ：使う可能性のある座標のペアを全部渡す
    Static2DSegTree(std::vector<std::pair<Coord, Coord>> pts) {
        std::sort(pts.begin(), pts.end());
        pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
        
        for (auto p : pts) xs.push_back(p.first);
        xs.erase(std::unique(xs.begin(), xs.end()), xs.end());

        n = 1;
        while (n < (int)xs.size()) n *= 2;

        ys.resize(2 * n);
        d.resize(2 * n);

        for (auto p : pts) {
            int x_idx = std::lower_bound(xs.begin(), xs.end(), p.first) - xs.begin();
            ys[x_idx + n].push_back(p.second);
        }

        for (int i = 0; i < n; i++) {
            std::sort(ys[i + n].begin(), ys[i + n].end());
            ys[i + n].erase(std::unique(ys[i + n].begin(), ys[i + n].end()), ys[i + n].end());
        }

        for (int i = n - 1; i > 0; i--) {
            std::set_union(ys[2 * i].begin(), ys[2 * i].end(),
                           ys[2 * i + 1].begin(), ys[2 * i + 1].end(),
                           std::back_inserter(ys[i]));
            ys[i].erase(std::unique(ys[i].begin(), ys[i].end()), ys[i].end());
        }

        for (int i = 1; i < 2 * n; i++) {
            d[i].assign(2 * ys[i].size(), e());
        }
    }

    void set(Coord x, Coord y, S val) {
        auto it = std::lower_bound(xs.begin(), xs.end(), x);
        assert(it != xs.end() && *it == x);
        int x_idx = (it - xs.begin()) + n;

        int y_idx = get_y_idx(x_idx, y);
        assert(y_idx < (int)ys[x_idx].size() && ys[x_idx][y_idx] == y);
        int p = y_idx + ys[x_idx].size();
        d[x_idx][p] = val;
        while (p >>= 1) d[x_idx][p] = op(d[x_idx][2 * p], d[x_idx][2 * p + 1]);

        x_idx >>= 1;
        while (x_idx > 0) {
            y_idx = get_y_idx(x_idx, y);
            p = y_idx + ys[x_idx].size();

            S l_val = e(), r_val = e();
            int l_child = 2 * x_idx, r_child = 2 * x_idx + 1;

            int l_y_idx = get_y_idx(l_child, y);
            if (l_y_idx < (int)ys[l_child].size() && ys[l_child][l_y_idx] == y) {
                l_val = d[l_child][l_y_idx + ys[l_child].size()];
            }
            int r_y_idx = get_y_idx(r_child, y);
            if (r_y_idx < (int)ys[r_child].size() && ys[r_child][r_y_idx] == y) {
                r_val = d[r_child][r_y_idx + ys[r_child].size()];
            }

            d[x_idx][p] = op(l_val, r_val);
            while (p >>= 1) d[x_idx][p] = op(d[x_idx][2 * p], d[x_idx][2 * p + 1]);
            x_idx >>= 1;
        }
    }

    // 現在の値に val を作用（加算）させるメソッド
    void add(Coord x, Coord y, S val) {
        auto it = std::lower_bound(xs.begin(), xs.end(), x);
        assert(it != xs.end() && *it == x);
        int x_idx = (it - xs.begin()) + n;

        int y_idx = get_y_idx(x_idx, y);
        assert(y_idx < (int)ys[x_idx].size() && ys[x_idx][y_idx] == y);
        int p = y_idx + ys[x_idx].size();
        
        // --- 上書きではなく op で累積 ---
        d[x_idx][p] = op(d[x_idx][p], val); 
        // --------------------------------
        
        while (p >>= 1) d[x_idx][p] = op(d[x_idx][2 * p], d[x_idx][2 * p + 1]);

        x_idx >>= 1;
        while (x_idx > 0) {
            y_idx = get_y_idx(x_idx, y);
            p = y_idx + ys[x_idx].size();

            S l_val = e(), r_val = e();
            int l_child = 2 * x_idx, r_child = 2 * x_idx + 1;

            int l_y_idx = get_y_idx(l_child, y);
            if (l_y_idx < (int)ys[l_child].size() && ys[l_child][l_y_idx] == y) {
                l_val = d[l_child][l_y_idx + ys[l_child].size()];
            }
            int r_y_idx = get_y_idx(r_child, y);
            if (r_y_idx < (int)ys[r_child].size() && ys[r_child][r_y_idx] == y) {
                r_val = d[r_child][r_y_idx + ys[r_child].size()];
            }

            d[x_idx][p] = op(l_val, r_val);
            while (p >>= 1) d[x_idx][p] = op(d[x_idx][2 * p], d[x_idx][2 * p + 1]);
            x_idx >>= 1;
        }
    }

    S prod(Coord min_x, Coord max_x, Coord min_y, Coord max_y) const {
        if (min_x >= max_x || min_y >= max_y) return e();
        int l = std::lower_bound(xs.begin(), xs.end(), min_x) - xs.begin();
        int r = std::lower_bound(xs.begin(), xs.end(), max_x) - xs.begin();
        l += n;
        r += n;
        S sml = e(), smr = e();
        while (l < r) {
            if (l & 1) sml = op(sml, inner_prod(l++, min_y, max_y));
            if (r & 1) smr = op(inner_prod(--r, min_y, max_y), smr);
            l >>= 1;
            r >>= 1;
        }
        return op(sml, smr);
    }

    S get(Coord x, Coord y) const {
        auto it = std::lower_bound(xs.begin(), xs.end(), x);
        // 登録されていないX座標なら単位元を返す
        if (it == xs.end() || *it != x) return e();
        
        int x_idx = (it - xs.begin()) + n;
        int y_idx = get_y_idx(x_idx, y);
        
        // 登録されていないY座標なら単位元を返す
        if (y_idx >= (int)ys[x_idx].size() || ys[x_idx][y_idx] != y) return e();
        
        // 葉ノードの値を直接返す
        return d[x_idx][y_idx + ys[x_idx].size()];
    }
};
