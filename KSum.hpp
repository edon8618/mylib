#pragma once

// Requires qandori's fps/fps.hpp to be included BEFORE this header.
// Uses its mint (mod 998244353), fps, and ntt(vector<mint>&, bool).
//
// Select exactly k DISTINCT INDICES from a nonnegative integer array.
// A may contain duplicates, and zero is allowed.
// Fingerprint[t] is a randomized weighted count of size-k subsets with sum t.
// - fingerprint[t] != 0  => a solution certainly exists.
// - fingerprint[t] == 0  => false negative probability <= k / 998244353
//   (assuming independent uniform random weights).
// A single random set of weights is shared by all queried targets.
//
// Time: O(n*k + k*L*log(L) + k*k*L), L = pow2ceil(2*limit+1).
// Space: O(k*L + n). Optional witness() takes O(k*L + k*limit)
// additional temporary memory and O(k*L*log(L)+n*k*k+k*k*limit) time.
// qandori ntt supports L <= 2^23; k must be < 998244353.
//
// Example:
//   #include "fps/fps.hpp"
//   #include "WeightedKSum.hpp"
//   ksum::WeightedKSum ks(a, 3, X);
//   if (ks.exists(X)) {
//       auto ids = ks.witness(X); // 0-based, all distinct
//   }

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <random>
#include <vector>

namespace ksum {

inline std::uint64_t random_seed() {
    std::random_device rd;
    return (std::uint64_t(rd()) << 32) ^ std::uint64_t(rd()) ^
        std::uint64_t(std::chrono::steady_clock::now().time_since_epoch().count());
}

struct WeightedKSum {
private:
    int k_ = 0, limit_ = 0, length_ = 1;
    std::vector<long long> a_;
    std::vector<mint> weight_;
    fps answer_;
    std::vector<fps> e_fft_;  // NTT(e_m), m=0,...,k-1; for witness reconstruction

    void build(std::uint64_t seed) {
        assert(k_ >= 0 && limit_ >= 0);
        for (auto x : a_) assert(x >= 0);
        assert(static_cast<std::uint64_t>(k_) < mint::mod());
        answer_ = fps(limit_ + 1);
        if (k_ == 0) {
            answer_[0] = 1;
            return;
        }
        if (k_ > static_cast<int>(a_.size())) return;

        const long long need = 2LL * limit_ + 1;
        length_ = 1;
        while (length_ < need && length_ <= (1 << 23)) length_ *= 2;
        assert(length_ <= (1 << 23));

        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<std::uint32_t> uniform(0, mint::mod() - 1);
        weight_.resize(a_.size());
        for (auto &w : weight_) w = mint::raw(uniform(rng));

        // p_j(z) = (-1)^(j-1) sum_i (w_i)^j z^(j*a_i).
        // Store the NTT of each p_j, padded to length_.
        std::vector<fps> p_fft(k_ + 1);
        for (int j = 1; j <= k_; ++j) p_fft[j].resize(length_);
        for (std::size_t i = 0; i < a_.size(); ++i) {
            const long long x = a_[i];
            if (x > limit_) continue;
            mint pw = 1;
            for (int j = 1; j <= k_ && x * j <= limit_; ++j) {
                pw *= weight_[i];
                if (j & 1) p_fft[j][x * j] += pw;
                else p_fft[j][x * j] -= pw;
            }
        }
        for (int j = 1; j <= k_; ++j) ntt(p_fft[j]);

        // e_0(z)=1, so NTT(e_0) is the all-ones vector.
        e_fft_.reserve(k_);
        e_fft_.emplace_back(length_, mint(1));

        // m*e_m = sum_{j=1}^m p_j * e_{m-j}.
        // Combine all products in frequency domain, then inverse NTT only once.
        for (int m = 1; m <= k_; ++m) {
            fps cur(length_);
            for (int j = 1; j <= m; ++j) {
                const fps &p = p_fft[j], &e = e_fft_[m - j];
                for (int t = 0; t < length_; ++t) cur[t] += p[t] * e[t];
            }
            ntt(cur, true);
            const mint inv_m = mint(m).inv();
            for (int t = 0; t <= limit_; ++t) cur[t] *= inv_m;
            if (m == k_) {
                answer_ = cur.pre(limit_ + 1);
                break;
            }
            std::fill(cur.begin() + limit_ + 1, cur.end(), mint(0));
            ntt(cur);
            e_fft_.push_back(std::move(cur));
        }
    }

public:
    template <class T>
    explicit WeightedKSum(const std::vector<T> &a, int k, int limit,
                          std::uint64_t seed = random_seed())
        : k_(k), limit_(limit), a_(a.begin(), a.end()) {
        build(seed);
    }

    int k() const { return k_; }
    int limit() const { return limit_; }
    const fps &fingerprints() const { return answer_; }

    mint fingerprint(int sum) const {
        if (sum < 0 || sum > limit_) return 0;
        return answer_[sum];
    }

    bool exists(int sum) const { return fingerprint(sum) != mint(0); }

    // Returns one actual witness (0-based indices) whenever exists(sum) is true.
    // Returns {} when the fingerprint is zero. For k=0, {} is a valid witness
    // exactly when sum=0, so check exists() to distinguish these cases.
    std::vector<int> witness(int sum) const {
        if (!exists(sum)) return {};
        if (k_ == 0) return {};

        // Recover e_0,...,e_{k-1} in coefficient form.
        std::vector<fps> e(k_);
        e[0] = fps(sum + 1);
        e[0][0] = 1;
        for (int j = 1; j < k_; ++j) {
            e[j] = e_fft_[j];
            ntt(e[j], true);
            e[j].resize(sum + 1);
        }

        std::vector<unsigned char> used(a_.size(), 0);
        std::vector<int> picked;
        picked.reserve(k_);
        int target = sum;
        for (int m = k_; m >= 1; --m) {
            int chosen = -1;
            for (int i = 0; i < static_cast<int>(a_.size()); ++i) {
                if (used[i] || a_[i] > target || weight_[i] == mint(0)) continue;
                const int x = static_cast<int>(a_[i]);
                mint pw = weight_[i], contribution = 0;
                for (int j = 0; j < m && 1LL * (j + 1) * x <= target; ++j) {
                    const mint term = pw * e[m - 1 - j][target - (j + 1) * x];
                    contribution += (j & 1) ? -term : term;
                    pw *= weight_[i];
                }
                if (contribution != mint(0)) {
                    chosen = i;
                    break;
                }
            }
            assert(chosen != -1); // follows from m*e_m[target] != 0
            const int x = static_cast<int>(a_[chosen]);
            const mint w = weight_[chosen];
            used[chosen] = 1;
            picked.push_back(chosen);

            // Divide E(t,z) by (1 + t*w*z^x):
            // e'_j = e_j - w*z^x*e'_{j-1} (update in increasing j).
            for (int j = 1; j < m; ++j) {
                for (int t = x; t <= target; ++t) {
                    e[j][t] -= w * e[j - 1][t - x];
                }
            }
            target -= x;
        }
        assert(target == 0);
        return picked;
    }
};

}  // namespace ksum
