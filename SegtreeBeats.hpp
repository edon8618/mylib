#pragma once
//https://judge.yosupo.jp/submission/358536
struct SegTreeBeats {
    static constexpr ll INF = 2e18;

    struct Node {
        ll sum, len;
        ll max_v, smax_v, max_c; // 最大値、2番目の最大値、最大値の個数
        ll min_v, smin_v, min_c; // 最小値、2番目の最小値、最小値の個数
        ll lazy_add;             // 区間加算の遅延タグ
    };

    int n, sz;
    vector<Node> dat;

    SegTreeBeats(const vector<ll>& a) {
        n = a.size();
        sz = 1;
        while (sz < n) sz *= 2;
        dat.resize(2 * sz);
        
        for (int i = 0; i < 2 * sz; ++i) {
            dat[i].lazy_add = 0;
            dat[i].len = 0;
            dat[i].max_v = dat[i].smax_v = -INF;
            dat[i].min_v = dat[i].smin_v = INF;
            dat[i].max_c = dat[i].min_c = 0;
        }
        for (int i = 0; i < n; ++i) {
            dat[sz + i].sum = dat[sz + i].max_v = dat[sz + i].min_v = a[i];
            dat[sz + i].len = dat[sz + i].max_c = dat[sz + i].min_c = 1;
        }
        for (int i = sz - 1; i >= 1; --i) update(i);
    }

    void add_val(int k, ll x) {
        if (x == 0) return;
        dat[k].sum += x * dat[k].len;
        dat[k].max_v += x;
        if (dat[k].smax_v != -INF) dat[k].smax_v += x;
        dat[k].min_v += x;
        if (dat[k].smin_v != INF) dat[k].smin_v += x;
        dat[k].lazy_add += x;
    }

    void update_chmax(int k, ll x) {
        if (dat[k].min_v >= x) return;
        dat[k].sum += (x - dat[k].min_v) * dat[k].min_c;
        if (dat[k].max_v == dat[k].min_v) {
            dat[k].max_v = x;
        } else if (dat[k].smax_v == dat[k].min_v) {
            dat[k].smax_v = x;
        }
        dat[k].min_v = x;
    }

    void update_chmin(int k, ll x) {
        if (dat[k].max_v <= x) return;
        dat[k].sum += (x - dat[k].max_v) * dat[k].max_c;
        if (dat[k].min_v == dat[k].max_v) {
            dat[k].min_v = x;
        } else if (dat[k].smin_v == dat[k].max_v) {
            dat[k].smin_v = x;
        }
        dat[k].max_v = x;
    }

    void push(int k) {
        if (dat[k].len <= 1) return;
        if (dat[k].lazy_add != 0) {
            add_val(2 * k, dat[k].lazy_add);
            add_val(2 * k + 1, dat[k].lazy_add);
            dat[k].lazy_add = 0;
        }
        if (dat[k].max_v < dat[2 * k].max_v) update_chmin(2 * k, dat[k].max_v);
        if (dat[k].min_v > dat[2 * k].min_v) update_chmax(2 * k, dat[k].min_v);
        if (dat[k].max_v < dat[2 * k + 1].max_v) update_chmin(2 * k + 1, dat[k].max_v);
        if (dat[k].min_v > dat[2 * k + 1].min_v) update_chmax(2 * k + 1, dat[k].min_v);
    }

    void update(int k) {
        dat[k].len = dat[2 * k].len + dat[2 * k + 1].len;
        dat[k].sum = dat[2 * k].sum + dat[2 * k + 1].sum;

        // 最大値のマージ
        if (dat[2 * k].max_v > dat[2 * k + 1].max_v) {
            dat[k].max_v = dat[2 * k].max_v;
            dat[k].max_c = dat[2 * k].max_c;
            dat[k].smax_v = max(dat[2 * k].smax_v, dat[2 * k + 1].max_v);
        } else if (dat[2 * k].max_v < dat[2 * k + 1].max_v) {
            dat[k].max_v = dat[2 * k + 1].max_v;
            dat[k].max_c = dat[2 * k + 1].max_c;
            dat[k].smax_v = max(dat[2 * k].max_v, dat[2 * k + 1].smax_v);
        } else {
            dat[k].max_v = dat[2 * k].max_v;
            dat[k].max_c = dat[2 * k].max_c + dat[2 * k + 1].max_c;
            dat[k].smax_v = max(dat[2 * k].smax_v, dat[2 * k + 1].smax_v);
        }

        // 最小値のマージ
        if (dat[2 * k].min_v < dat[2 * k + 1].min_v) {
            dat[k].min_v = dat[2 * k].min_v;
            dat[k].min_c = dat[2 * k].min_c;
            dat[k].smin_v = min(dat[2 * k].smin_v, dat[2 * k + 1].min_v);
        } else if (dat[2 * k].min_v > dat[2 * k + 1].min_v) {
            dat[k].min_v = dat[2 * k + 1].min_v;
            dat[k].min_c = dat[2 * k + 1].min_c;
            dat[k].smin_v = min(dat[2 * k].min_v, dat[2 * k + 1].smin_v);
        } else {
            dat[k].min_v = dat[2 * k].min_v;
            dat[k].min_c = dat[2 * k].min_c + dat[2 * k + 1].min_c;
            dat[k].smin_v = min(dat[2 * k].smin_v, dat[2 * k + 1].smin_v);
        }
    }

    void _apply_chmin(int l, int r, ll x, int k, int sl, int sr) {
        if (r <= sl || sr <= l || dat[k].max_v <= x) return;
        if (l <= sl && sr <= r && dat[k].smax_v < x) {
            update_chmin(k, x);
            return;
        }
        push(k);
        int sm = (sl + sr) / 2;
        _apply_chmin(l, r, x, 2 * k, sl, sm);
        _apply_chmin(l, r, x, 2 * k + 1, sm, sr);
        update(k);
    }

    void _apply_chmax(int l, int r, ll x, int k, int sl, int sr) {
        if (r <= sl || sr <= l || dat[k].min_v >= x) return;
        if (l <= sl && sr <= r && dat[k].smin_v > x) {
            update_chmax(k, x);
            return;
        }
        push(k);
        int sm = (sl + sr) / 2;
        _apply_chmax(l, r, x, 2 * k, sl, sm);
        _apply_chmax(l, r, x, 2 * k + 1, sm, sr);
        update(k);
    }

    void _apply_add(int l, int r, ll x, int k, int sl, int sr) {
        if (r <= sl || sr <= l) return;
        if (l <= sl && sr <= r) {
            add_val(k, x);
            return;
        }
        push(k);
        int sm = (sl + sr) / 2;
        _apply_add(l, r, x, 2 * k, sl, sm);
        _apply_add(l, r, x, 2 * k + 1, sm, sr);
        update(k);
    }

    ll _query_sum(int l, int r, int k, int sl, int sr) {
        if (r <= sl || sr <= l) return 0;
        if (l <= sl && sr <= r) return dat[k].sum;
        push(k);
        int sm = (sl + sr) / 2;
        return _query_sum(l, r, 2 * k, sl, sm) + _query_sum(l, r, 2 * k + 1, sm, sr);
    }

public:
    // [l, r) の各要素 v_i を min(v_i, x) に更新
    void apply_chmin(int l, int r, ll x) { _apply_chmin(l, r, x, 1, 0, sz); }
    
    // [l, r) の各要素 v_i を max(v_i, x) に更新
    void apply_chmax(int l, int r, ll x) { _apply_chmax(l, r, x, 1, 0, sz); }
    
    // [l, r) の各要素に x を加算
    void apply_add(int l, int r, ll x) { _apply_add(l, r, x, 1, 0, sz); }
    
    // [l, r) の和を取得
    ll query_sum(int l, int r) { return _query_sum(l, r, 1, 0, sz); }
};
