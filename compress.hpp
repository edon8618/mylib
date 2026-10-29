#pragma once
template <typename T>
struct Compress {
    vector<T> xs;
    void add(T x) { xs.push_back(x); }
    void build() {
        sort(all(xs));
        xs.erase(unique(all(xs)), xs.end());
    }
    // 値 x が何番目か返す
    int to_idx(T x) {
        return lower_bound(all(xs), x) - xs.begin();
    }
    // k番目の値は何か返す
    T to_val(int k) { return xs[k]; }
    int size() { return xs.size(); }
};