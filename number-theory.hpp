#pragma once

/*
============================================================
Number Theory Utility
============================================================

[素数判定]
bool nt::is_prime(u64 n);

[mod乗算]
u64 nt::mul_mod(u64 a, u64 b, u64 mod);

[mod累乗]
u64 nt::pow_mod(u64 a, u64 e, u64 mod);
    a^e % mod

[素因数分解]
vector<u64> nt::factor(u64 n);
        nt::factor(12)
        -> {2, 2, 3}

[素因数分解 + 指数]
vector<pair<u64,int>> nt::factor_count(u64 n);
        nt::factor_count(12)
        -> {{2,2}, {3,1}}

[相異なる素因数]
vector<u64> nt::distinct_prime_factors(u64 n);
        nt::distinct_prime_factors(12)
        -> {2,3}

[原始根判定]
bool nt::is_primitive_root(u64 g, u64 p);
    p は素数であること。
    g が mod p の原始根なら true。

[原始根判定 / p-1の素因数既知版]
bool nt::is_primitive_root(
    u64 g,
    u64 p,
    const vector<u64>& pf
);
    pf は p-1 の相異なる素因数。(nt::distinct_prime_factors)
    同じ p で何度も判定するときはこちらを使う。

[原始根探索]
u64 nt::primitive_root(u64 p);

    p は素数であること。
    mod p の原始根を1つ返す。
------------------------------------------------------------
注意
------------------------------------------------------------
・ライブラリ内部は uint64_t 前提。
・long long の非負値ならそのまま渡してよい。
・負の long long を渡すと uint64_t に変換されて壊れる。
・is_primitive_root / primitive_root の p は素数であること。
・factor は Pollard Rho + Miller-Rabin。
・is_prime は 64bit 全域で deterministic。
============================================================
*/

namespace nt {

using u64 = uint64_t;
using u128 = __uint128_t;

// ============================================================
// modular arithmetic
// ============================================================

inline u64 mul_mod(u64 a, u64 b, u64 mod) {
    return (u128)a * b % mod;
}

u64 pow_mod(u64 a, u64 e, u64 mod) {
    assert(mod > 0);
    u64 r = 1 % mod;

    while (e) {
        if (e & 1) {
            r = mul_mod(r, a, mod);
        }

        a = mul_mod(a, a, mod);
        e >>= 1;
    }

    return r;
}


// ============================================================
// deterministic Miller-Rabin for uint64_t
// ============================================================

bool is_prime(u64 n) {
    if (n < 2) return false;

    static constexpr u64 small_primes[] = {
        2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL,
        17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL
    };

    for (u64 p : small_primes) {
        if (n % p == 0) {
            return n == p;
        }
    }

    u64 d = n - 1;
    int s = 0;

    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    auto witness = [&](u64 a) -> bool {
        if (a % n == 0) return false;

        u64 x = pow_mod(a, d, n);

        if (x == 1 || x == n - 1) {
            return false;
        }

        for (int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);

            if (x == n - 1) {
                return false;
            }
        }

        return true;
    };

    // Deterministic for all uint64_t values.
    static constexpr u64 bases[] = {
        2ULL,
        325ULL,
        9375ULL,
        28178ULL,
        450775ULL,
        9780504ULL,
        1795265022ULL
    };

    for (u64 a : bases) {
        if (witness(a)) {
            return false;
        }
    }

    return true;
}


// ============================================================
// splitmix64 RNG
// ============================================================

struct SplitMix64 {
    u64 x;

    explicit SplitMix64(
        u64 seed = chrono::steady_clock::now()
                       .time_since_epoch()
                       .count()
    ) : x(seed) {}

    u64 operator()() {
        u64 z = (
            x += 0x9e3779b97f4a7c15ULL
        );

        z = (
            z ^ (z >> 30)
        ) * 0xbf58476d1ce4e5b9ULL;

        z = (
            z ^ (z >> 27)
        ) * 0x94d049bb133111ebULL;

        return z ^ (z >> 31);
    }

    u64 range(u64 l, u64 r) {
        return l + (*this)() % (r - l + 1);
    }
};

static SplitMix64 rng;


// ============================================================
// Pollard Rho (Brent)
// ============================================================

u64 pollard_rho(u64 n) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    if (n % 5 == 0) return 5;

    constexpr u64 BATCH = 128;

    while (true) {
        u64 y = rng.range(1, n - 1);
        u64 c = rng.range(1, n - 1);

        auto f = [&](u64 x) -> u64 {
            return (
                mul_mod(x, x, n) + c
            ) % n;
        };

        u64 g = 1;
        u64 r = 1;
        u64 q = 1;

        u64 x = 0;
        u64 ys = 0;

        while (g == 1) {
            x = y;

            for (u64 i = 0; i < r; ++i) {
                y = f(y);
            }

            for (
                u64 k = 0;
                k < r && g == 1;
                k += BATCH
            ) {
                ys = y;
                q = 1;

                u64 lim = min(
                    BATCH,
                    r - k
                );

                for (u64 i = 0; i < lim; ++i) {
                    y = f(y);

                    u64 diff = (
                        x > y
                        ? x - y
                        : y - x
                    );

                    q = mul_mod(
                        q,
                        diff,
                        n
                    );
                }

                g = gcd(q, n);
            }

            r <<= 1;
        }

        // Batch gcd returned n itself.
        if (g == n) {
            do {
                ys = f(ys);

                u64 diff = (
                    x > ys
                    ? x - ys
                    : ys - x
                );

                g = gcd(
                    diff,
                    n
                );

            } while (g == 1);
        }

        if (g != n) {
            return g;
        }
    }
}


// ============================================================
// factorization
// ============================================================

void factor_rec(
    u64 n,
    vector<u64>& res
) {
    if (n == 1) {
        return;
    }

    if (is_prime(n)) {
        res.push_back(n);
        return;
    }

    u64 d = pollard_rho(n);

    factor_rec(d, res);
    factor_rec(n / d, res);
}


// 重複込み・昇順
vector<u64> factor(u64 n) {
    vector<u64> res;

    if (n <= 1) {
        return res;
    }

    factor_rec(n, res);

    sort(
        res.begin(),
        res.end()
    );

    return res;
}


// (prime, exponent)
vector<pair<u64, int>> factor_count(
    u64 n
) {
    vector<u64> fs = factor(n);

    vector<pair<u64, int>> res;

    for (u64 p : fs) {
        if (
            res.empty() ||
            res.back().first != p
        ) {
            res.push_back({
                p,
                1
            });
        }
        else {
            ++res.back().second;
        }
    }

    return res;
}


// 相異なる素因数のみ
vector<u64> distinct_prime_factors(
    u64 n
) {
    vector<u64> fs = factor(n);

    fs.erase(
        unique(
            fs.begin(),
            fs.end()
        ),
        fs.end()
    );

    return fs;
}


// ============================================================
// primitive root
// ============================================================

// p must be prime.
// pf = distinct prime factors of p-1.
bool is_primitive_root(
    u64 g,
    u64 p,
    const vector<u64>& pf
) {
    if (p == 2) {
        return g % 2 == 1;
    }

    if (g == 0 || g >= p) {
        return false;
    }

    for (u64 q : pf) {
        if (
            pow_mod(
                g,
                (p - 1) / q,
                p
            ) == 1
        ) {
            return false;
        }
    }

    return true;
}


// p must be prime.
// Factorizes p-1 internally.
bool is_primitive_root(
    u64 g,
    u64 p
) {
    auto pf =
        distinct_prime_factors(
            p - 1
        );

    return is_primitive_root(
        g,
        p,
        pf
    );
}


// p must be prime.
u64 primitive_root(u64 p) {
    if (p == 2) {
        return 1;
    }

    auto pf =
        distinct_prime_factors(
            p - 1
        );

    for (u64 g = 2;; ++g) {
        if (
            is_primitive_root(
                g,
                p,
                pf
            )
        ) {
            return g;
        }
    }
}

} // namespace nt
