#pragma once
#include <vector>
#include <stdexcept>
// Prerequisite: qandori fps/fps.hpp and fps/taylor_shift.hpp.
inline fps stirling2_n_fix(int n) {
    if (n < 0 || (unsigned)n >= mint::mod()) throw std::invalid_argument("stirling: require 0 <= n < mod");
    if (n == 0) return fps{1};
    fps A(n + 1), B(n + 1);
    std::vector<mint> invf(n + 1);
    mint fact = 1;
    for (int i = 1; i <= n; ++i) fact *= i;
    invf[n] = fact.inv();
    for (int i = n; i > 0; --i) invf[i-1] = invf[i] * i;
    A[0] = 1;
    B[0] = 0;
    for (int i = 1; i <= n; i++) {
        A[i] = (i % 2 == 1 ? -invf[i] : invf[i]);
        B[i] = mint(i).pow(n) * invf[i];
    }
    fps res = A * B;
    res.resize(n + 1);
    return res;
}

inline fps stirling1_n_fix(int n, bool signed_res = false) {
    if (n < 0 || (unsigned)n >= mint::mod()) throw std::invalid_argument("stirling: require 0 <= n < mod");
    if (n == 0) return fps{1};
    fps res = {0, 1};
    int m = 1;
    int h = 31 - __builtin_clz(n);

    for (int i = h - 1; i >= 0; i--) {
        res = res * taylor_shift(res, mint(m));
        m *= 2;

        if ((n >> i) & 1) {
            fps next_res(res.size() + 1, 0);
            for (int j = 0; j < (int)res.size(); j++) {
                next_res[j + 1] += res[j];
                next_res[j] += res[j] * mint(m);
            }
            res = next_res;
            m++;
        }
    }

    if (signed_res) {
        for (int k = 0; k <= n; k++) {
            if ((n - k) % 2 != 0) res[k] = -res[k];
        }
    }
    return res;
}

inline fps stirling2_k_fix(int k, int max_n) {
    if (k < 0 || max_n < 0 || (unsigned)max_n >= mint::mod()) throw std::invalid_argument("stirling: invalid k or max_n");
    if (k > max_n) return fps(max_n + 1, 0);
    if (k == 0) {
        fps res(max_n + 1, 0);
        res[0] = 1;
        return res;
    }
    std::vector<mint> fact(max_n + 1), invf(max_n + 1);
    fact[0] = 1; invf[0] = 1;
    for (int i = 1; i <= max_n; i++) fact[i] = fact[i - 1] * i;
    invf[max_n] = fact[max_n].inv();
    for (int i = max_n; i > 0; i--) invf[i - 1] = invf[i] * i;

    fps f(max_n + 1, 0);
    for (int i = 1; i <= max_n; i++) f[i] = invf[i];

    fps res = f.pow(k);
    mint inv_k_fact = invf[k];

    for (int n = k; n <= max_n; n++) {
        res[n] *= fact[n] * inv_k_fact;
    }
    for (int n = 0; n < k; n++) res[n] = 0;
    return res;
}

inline fps stirling1_k_fix(int k, int max_n, bool signed_res = false) {
    if (k < 0 || max_n < 0 || (unsigned)max_n >= mint::mod()) throw std::invalid_argument("stirling: invalid k or max_n");
    if (k > max_n) return fps(max_n + 1, 0);
    if (k == 0) {
        fps res(max_n + 1, 0);
        res[0] = 1;
        return res;
    }
    std::vector<mint> fact(max_n + 1), invf(max_n + 1);
    fact[0] = 1; invf[0] = 1;
    for (int i = 1; i <= max_n; i++) fact[i] = fact[i - 1] * i;
    invf[max_n] = fact[max_n].inv();
    for (int i = max_n; i > 0; i--) invf[i - 1] = invf[i] * i;

    std::vector<mint> inv(max_n + 1);
    if (max_n > 0) inv[1] = 1;
    for (int i = 2; i <= max_n; i++) inv[i] = -inv[mint::mod() % i] * (mint::mod() / i);

    fps f(max_n + 1, 0);
    for (int i = 1; i <= max_n; i++) f[i] = inv[i];

    fps res = f.pow(k);
    mint inv_k_fact = invf[k];

    for (int n = k; n <= max_n; n++) {
        res[n] *= fact[n] * inv_k_fact;
        if (signed_res && (n - k) % 2 != 0) res[n] = -res[n];
    }
    for (int n = 0; n < k; n++) res[n] = 0;
    return res;
}
