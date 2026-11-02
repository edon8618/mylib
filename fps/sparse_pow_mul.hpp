#pragma once
#include <vector>
#include <map>
#include <stdexcept>
#include <utility>
#include <algorithm>

// Prerequisite: qandori fps/fps.hpp (fps, mint).
namespace sparse_pow_mul_detail {
inline mint signed_power(mint a, long long k) {
    unsigned long long e = k < 0 ? 0ULL - (unsigned long long)k : (unsigned long long)k;
    if (k < 0) a = a.inv();
    mint r = 1;
    while (e) { if (e & 1) r *= a; a *= a; e >>= 1; }
    return r;
}
}

// Negative exponents require a nonzero constant coefficient (no Laurent series).
// Zero to the zeroth power is 1. Output has exactly deg coefficients.
inline fps sparse_pow_mul(const fps& f, long long n, const fps& g, long long m, int deg) {
    if (deg < 0 || (unsigned)deg >= mint::mod())
        throw std::invalid_argument("sparse_pow_mul: require 0 <= deg < mod");
    if ((n < 0 && (f.empty() || f[0] == mint(0))) ||
        (m < 0 && (g.empty() || g[0] == mint(0))))
        throw std::invalid_argument("sparse_pow_mul: negative power requires nonzero constant");
    if (!deg) return {};
    auto valuation = [](const fps& a) {
        int s = 0;
        while (s < (int)a.size() && a[s] == mint(0)) ++s;
        return s;
    };
    int sf = n == 0 ? 0 : valuation(f), sg = m == 0 ? 0 : valuation(g);
    if ((n > 0 && sf == (int)f.size()) || (m > 0 && sg == (int)g.size())) return fps(deg);
    __int128 shift = (__int128)sf * n + (__int128)sg * m;
    if (shift >= deg) return fps(deg);
    int off = (int)shift, len = deg - off;
    using Term = std::pair<int, mint>;
    auto collect = [&](const fps& a, long long k, int s) {
        std::vector<Term> v{{0, 1}};
        if (k == 0) return v;
        mint inv = a[s].inv();
        for (int j = 1; j < len && s + j < (int)a.size(); ++j)
            if (a[s+j] != mint(0)) v.emplace_back(j, a[s+j] * inv);
        return v;
    };
    auto a = collect(f, n, sf), b = collect(g, m, sg);
    // A H' = B H, A = f_hat*g_hat, B = n*f_hat'*g_hat + m*g_hat'*f_hat.
    std::map<int, mint> am, bm;
    for (auto [i, x] : a) for (auto [j, y] : b) {
        long long d = (long long)i + j;
        if (d >= len) continue;
        mint z = x*y;
        if (d) {
            am[(int)d] += z;
            bm[(int)d-1] += (mint(n)*i + mint(m)*j)*z;
        }
    }
    std::vector<Term> A, B;
    for (auto [i,x] : am) if (x != mint(0)) A.emplace_back(i,x);
    for (auto [i,x] : bm) if (x != mint(0)) B.emplace_back(i,x);
    std::vector<mint> inv(len);
    if (len > 1) inv[1] = 1;
    for (int i = 2; i < len; ++i) inv[i] = -inv[mint::mod()%i] * (mint::mod()/i);
    fps h(len); h[0] = 1;
    for (int k = 1; k < len; ++k) {
        mint sum = 0;
        for (auto [j,x] : B) { if (j >= k) break; sum += x*h[k-1-j]; }
        for (auto [j,x] : A) { if (j >= k) break; sum -= x*(k-j)*h[k-j]; }
        h[k] = sum*inv[k];
    }
    mint c = 1;
    if (n) c *= sparse_pow_mul_detail::signed_power(f[sf],n);
    if (m) c *= sparse_pow_mul_detail::signed_power(g[sg],m);
    fps res(deg);
    for (int i = 0; i < len; ++i) res[off+i] = h[i]*c;
    return res;
}
