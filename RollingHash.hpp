#pragma once
// https://mofecoder.com/contests/kyoprocamp2026east/submissions/28274
struct RollingHash {
    static const uint64_t MOD = (1ULL << 61) - 1;
    static uint64_t base;
    static vector<uint64_t> power;

    vector<uint64_t> hash;

    // 実行時に一度だけ呼んで base を乱数初期化する
    static void init_base() {
        if (base > 0) return;
        // 時間をシードにして乱数を生成 (2からMOD-2の範囲)
        uint64_t seed = chrono::steady_clock::now().time_since_epoch().count();
        // 簡易的なXorshift等でバラけさせる
        seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
        base = (seed % (MOD - 3)) + 2;
    }

    RollingHash(const string& s) {
        if (base == 0) init_base();
        int n = s.size();
        hash.assign(n + 1, 0);
        
        while (power.size() <= n) {
            power.push_back(Mul(power.back(), base));
        }

        for (int i = 0; i < n; ++i) {
            hash[i + 1] = Mul(hash[i], base) + s[i];
            if (hash[i + 1] >= MOD) hash[i + 1] -= MOD;
        }
    }

    // [l, r) のハッシュ値を O(1) で取得
    uint64_t get(int l, int r) const {
        uint64_t ret = hash[r] + MOD - Mul(hash[l], power[r - l]);
        if (ret >= MOD) ret -= MOD;
        return ret;
    }

    // O(1) ハッシュ結合
    static uint64_t concat(uint64_t h1, uint64_t h2, int len2) {
        uint64_t ret = Mul(h1, power[len2]) + h2;
        if (ret >= MOD) ret -= MOD;
        return ret;
    }

private:
    // a * b (mod 2^61 - 1) を __uint128_t で高速計算
    static inline uint64_t Mul(uint64_t a, uint64_t b) {
        __uint128_t t = (__uint128_t)a * b;
        // 上位ビットと下位ビットを足し合わせる (mod 2^61 - 1 の性質)
        t = (t >> 61) + (t & MOD);
        if (t >= MOD) t -= MOD;
        return (uint64_t)t;
    }
};

// static メンバ変数の実体定義
uint64_t RollingHash::base = 0;
vector<uint64_t> RollingHash::power = {1};