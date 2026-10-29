#pragma once
// https://atcoder.jp/contests/abc244/submissions/73762024
// 有理数
struct yuri {
    long long bunsi; // 分子
    long long bunbo; // 分母

    // コンストラクタ: 引数を __int128_t で受けるのが最大のポイント。
    // 四則演算で発生した巨大な数を受け止め、約分してから long long に安全に落とし込む。
    yuri(__int128_t si = 0, __int128_t bo = 1) {
        assert(bo != 0); // ゼロ除算を検知してREを出させる（デバッグ用）
        
        // 分母は常に正に保つ（比較演算をバグらせないため）
        if (bo < 0) {
            si = -si;
            bo = -bo;
        }

        // __int128_t 用のユークリッド互除法
        __int128_t a = si < 0 ? -si : si;
        __int128_t b = bo;
        while (b) {
            __int128_t temp = a % b;
            a = b;
            b = temp;
        }

        bunsi = (long long)(si / a);
        bunbo = (long long)(bo / a);
    }

    // --- 比較演算子 ---
    // ここも比較の瞬間のクロス掛け算で long long を突破するのでキャスト必須
    bool operator<(const yuri& rhs) const {
        return (__int128_t)bunsi * rhs.bunbo < (__int128_t)rhs.bunsi * bunbo;
    }
    bool operator>(const yuri& rhs) const { return rhs < *this; }
    bool operator<=(const yuri& rhs) const { return !(*this > rhs); }
    bool operator>=(const yuri& rhs) const { return !(*this < rhs); }
    bool operator==(const yuri& rhs) const {
        // コンストラクタで必ず約分＆分母正化されるため、単純比較でOK
        return bunsi == rhs.bunsi && bunbo == rhs.bunbo;
    }
    bool operator!=(const yuri& rhs) const { return !(*this == rhs); }

    // --- 四則演算 ---
    // 計算前に __int128_t にキャストして渡し、コンストラクタ内で約分させる
    yuri operator+(const yuri& rhs) const {
        return yuri((__int128_t)bunsi * rhs.bunbo + (__int128_t)rhs.bunsi * bunbo, (__int128_t)bunbo * rhs.bunbo);
    }
    yuri operator-(const yuri& rhs) const {
        return yuri((__int128_t)bunsi * rhs.bunbo - (__int128_t)rhs.bunsi * bunbo, (__int128_t)bunbo * rhs.bunbo);
    }
    yuri operator*(const yuri& rhs) const {
        return yuri((__int128_t)bunsi * rhs.bunsi, (__int128_t)bunbo * rhs.bunbo);
    }
    yuri operator/(const yuri& rhs) const {
        return yuri((__int128_t)bunsi * rhs.bunbo, (__int128_t)bunbo * rhs.bunsi);
    }

    // --- ユーティリティ ---
    friend yuri abs(const yuri& a) {
        return yuri(a.bunsi < 0 ? -a.bunsi : a.bunsi, a.bunbo);
    }

    // cout で直接出力できるようにするオーバーロード
    friend std::ostream& operator<<(std::ostream& os, const yuri& a) {
        if (a.bunbo == 1) os << a.bunsi;
        else os << a.bunsi << "/" << a.bunbo;
        return os;
    }
};