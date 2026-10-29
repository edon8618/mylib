#pragma once
/**
 * @brief Point Add Rectangle Sum (2次元点追加・矩形和クエリ)
 * @tparam C 座標の型 (int, long long など)
 * @tparam T 重み・和の型 (long long, mint など)
 * * 内部実装: CDQ分割統治 + 動的抽出std::sort + 1次元 Fenwick Tree
 * 時間計算量: O((N + Q) log^2 (N + Q))
 * 空間計算量: O(N + Q)
 */
template <typename C = long long, typename T = long long>
struct PointAddRectangleSum {
private:
    struct RawOp {
        int type; // 0: point, 1: query
        C x1, y1, x2, y2;
        T w;
        int q_id;
    };

    struct Event {
        int type; // 0: query, 1: point (ソート時にクエリを先行させるため 0 < 1)
        C x;
        int ylo, yhi;
        int sign; // +1 or -1 (クエリ用)
        T w;      // 点の重み
        int q_id;

        bool operator<(const Event& e) const {
            if (x != e.x) return x < e.x;
            return type < e.type; // 同じX座標ではクエリ(0)を点(1)より先に処理
        }
    };

    struct FenwickTree {
        int n;
        std::vector<T> data;
        FenwickTree(int n = 0) : n(n), data(n + 1, 0) {}

        void add(int i, T v) {
            for (int idx = i + 1; idx <= n; idx += idx & -idx) data[idx] += v;
        }

        T sum(int i) const {
            T s = 0;
            for (int idx = i + 1; idx > 0; idx -= idx & -idx) s += data[idx];
            return s;
        }

        T query(int l, int r) const {
            if (l > r) return 0;
            return sum(r) - sum(l - 1);
        }
    };

    std::vector<RawOp> ops;
    int q_count = 0;

public:
    PointAddRectangleSum() = default;

    /**
     * @brief 点 (x, y) に重み w を加算するクエリを登録 O(1)
     */
    void add_point(C x, C y, T w) {
        ops.push_back({0, x, y, 0, 0, w, -1});
    }

    /**
     * @brief 矩形領域 [l, r) x [d, u) の重み和を求めるクエリを登録 O(1)
     * l <= x < r 且つ d <= y < u
     */
    void add_query(C l, C r, C d, C u) {
        ops.push_back({1, l, d, r, u, 0, q_count++});
    }

    /**
     * @brief クエリをまとめて処理し、各 add_query の答えを配列で返す
     */
    std::vector<T> solve() {
        std::vector<T> ans(q_count, 0);
        if (ops.empty()) return ans;

        // 1. Y座標の座標圧縮
        std::vector<C> ys;
        ys.reserve(ops.size());
        for (const auto& op : ops) {
            if (op.type == 0) ys.push_back(op.y1);
        }
        std::sort(ys.begin(), ys.end());
        ys.erase(std::unique(ys.begin(), ys.end()), ys.end());

        int y_max = ys.size();
        if (y_max == 0) return ans; // 点が1つもない場合

        // 2. 内部用構造体への変換
        struct PointOp { C x; int y; T w; };
        struct QueryOp { C l, r; int ylo, yhi; bool valid; int q_id; };
        struct InternalOp {
            int type;
            PointOp pt;
            QueryOp q;
        };

        std::vector<InternalOp> items(ops.size());
        for (size_t i = 0; i < ops.size(); ++i) {
            if (ops[i].type == 0) {
                int y_idx = std::lower_bound(ys.begin(), ys.end(), ops[i].y1) - ys.begin();
                items[i] = {0, {ops[i].x1, y_idx, ops[i].w}, {}};
            } else {
                int ylo = std::lower_bound(ys.begin(), ys.end(), ops[i].y1) - ys.begin();
                int yhi = int(std::lower_bound(ys.begin(), ys.end(), ops[i].y2) - ys.begin()) - 1;
                bool valid = (ylo <= yhi);
                items[i] = {1, {}, {ops[i].x1, ops[i].x2, ylo, yhi, valid, ops[i].q_id}};
            }
        }

        FenwickTree ft(y_max);
        std::vector<Event> events;

        // 3. CDQ分割統治 (時間軸 l ~ r)
        auto cdq = [&](auto& self, int l, int r) -> void {
            if (l >= r) return;
            int m = l + (r - l) / 2;
            self(self, l, m);
            self(self, m + 1, r);

            events.clear();

            // 左半分 (過去) の点追加イベントを抽出
            for (int i = l; i <= m; ++i) {
                if (items[i].type == 0) {
                    events.push_back({1, items[i].pt.x, items[i].pt.y, items[i].pt.y, 0, items[i].pt.w, -1});
                }
            }
            if (events.empty()) return;

            int left_pt_cnt = events.size();

            // 右半分 (未来) の矩形和クエリイベントを抽出 (lで-1, rで+1)
            for (int i = m + 1; i <= r; ++i) {
                if (items[i].type == 1 && items[i].q.valid) {
                    events.push_back({0, items[i].q.l, items[i].q.ylo, items[i].q.yhi, -1, 0, items[i].q.q_id});
                    events.push_back({0, items[i].q.r, items[i].q.ylo, items[i].q.yhi,  1, 0, items[i].q.q_id});
                }
            }

            if (int(events.size()) == left_pt_cnt) return;

            // X座標順にソート (Xが同じ場合はクエリ(type=0)を点(type=1)より先に処理することで [l, r) を実現)
            std::sort(events.begin(), events.end());

            // イベントの適用と計算
            for (const auto& e : events) {
                if (e.type == 1) {
                    ft.add(e.ylo, e.w);
                } else {
                    if (e.sign == 1) {
                        ans[e.q_id] += ft.query(e.ylo, e.yhi);
                    } else {
                        ans[e.q_id] -= ft.query(e.ylo, e.yhi);
                    }
                }
            }

            // Fenwick Treeのロールバック
            for (const auto& e : events) {
                if (e.type == 1) {
                    ft.add(e.ylo, -e.w);
                }
            }
        };

        cdq(cdq, 0, items.size() - 1);
        return ans;
    }
};
