/** https://judge.yosupo.jp/submission/364786
 *  https://atcoder.jp/contests/arc117/submissions/74751959
 * 任意 mod における二項係数
 * mが中くらいのとき、もしくは、nが中くらいでmが合成数のとき　に使う
 * コンストラクタ　ArbitraryModBinomial comb(mod); O(min(M,N_max))
 * クエリ comb.nCr(ll n,ll r) O(logN)
 */
class ArbitraryModBinomial {
private:
    static long long modinv(long long a, long long m) {
        long long b = m, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b; std::swap(a, b);
            u -= t * v; std::swap(u, v);
        }
        u %= m;
        if (u < 0) u += m;
        return u;
    }

    struct Barrett {
        unsigned long long m;
        unsigned long long im;
        explicit Barrett(unsigned long long m) : m(m), im((unsigned long long)(-1) / m + 1) {}
        unsigned long long mul(unsigned long long a, unsigned long long b) const {
            return (unsigned __int128)a * b % m;
            unsigned long long z = a * b;
            unsigned long long x = (unsigned long long)(((unsigned __int128)(z) * im) >> 64);
            unsigned long long v = z - x * m;
            if (m <= v) v += m;
            return v;
        }
    };

    struct BinomialModPe {
        long long p, pe;
        int e_max;
        Barrett br;
        std::vector<long long> fact;
        std::vector<long long> p_pow;

        BinomialModPe(long long p, long long pe, long long n_max) : p(p), pe(pe), br(pe) {
            e_max = 0;
            long long tmp = pe;
            while (tmp > 1) { tmp /= p; e_max++; }

            p_pow.assign(e_max + 1, 1);
            for (int i = 1; i <= e_max; ++i) {
                p_pow[i] = p_pow[i - 1] * p;
            }

            // ここがキモ：n_max と p^e の小さい方だけ前計算する
            long long lim = std::min(n_max, pe);
            
            // 限界チェック (約2000万を超えたら現実的な時間とメモリで解けない)
            if (lim > 20000000) {
                throw std::runtime_error("エラー: n_max と m の素因数(p^e) が両方とも巨大すぎます。別の数学的考察が必要です。");
            }

            fact.assign(lim + 1, 1);
            for (long long i = 1; i <= lim; ++i) {
                if (i % p == 0) fact[i] = fact[i - 1];
                else fact[i] = br.mul(fact[i - 1], i);
            }
        }

        long long fact_no_p(unsigned long long n) const {
            unsigned long long res = 1;
            bool negate = false;
            long long lim = fact.size() - 1; 
            
            while (n > 0) {
                // p^e が lim 以下の場合のみ Wilson の定理の条件を確認可能
                if (pe <= lim && fact[pe - 1] == pe - 1 && ((n / pe) & 1)) {
                    negate = !negate;
                }
                res = br.mul(res, fact[n % pe]);
                n /= p;
            }
            return negate ? pe - res : res;
        }

        long long nCr(unsigned long long n, unsigned long long r) const {
            if (n < r) return 0;
            
            unsigned long long n_copy = n, r_copy = r, nr_copy = n - r;
            long long e = 0;
            while (n_copy) e += (n_copy /= p);
            while (r_copy) e -= (r_copy /= p);
            while (nr_copy) e -= (nr_copy /= p);
            
            if (e >= e_max) return 0; 

            long long num = fact_no_p(n);
            long long den = br.mul(fact_no_p(r), fact_no_p(n - r));
            
            long long res = br.mul(num, modinv(den, pe));
            res = br.mul(res, p_pow[e]);
            return res;
        }
    };

    long long m;
    long long n_max_limit;
    std::vector<long long> M;
    std::vector<long long> crt_coeffs;
    std::vector<BinomialModPe> bims;

public:
    // コンストラクタで n_max を受け取る
    ArbitraryModBinomial(long long m, long long n_max) : m(m), n_max_limit(n_max) {
        if (m <= 0 || n_max < 0) throw std::invalid_argument("mod must be positive; n_max must be nonnegative");
        if (m == 1) return;
        long long tmp = m;
        for (long long i = 2; i <= tmp / i; ++i) {
            if (tmp % i == 0) {
                long long pe = 1;
                while (tmp % i == 0) {
                    pe *= i;
                    tmp /= i;
                }
                M.push_back(pe);
                bims.emplace_back(i, pe, n_max); // n_max を伝搬
            }
        }
        if (tmp > 1) {
            M.push_back(tmp);
            bims.emplace_back(tmp, tmp, n_max); // n_max を伝搬
        }
        
        for (size_t i = 0; i < M.size(); ++i) {
            long long m_i = M[i];
            long long m_other = m / m_i;
            long long inv = modinv(m_other, m_i);
            crt_coeffs.push_back((unsigned __int128)m_other * inv % m);
        }
    }

    long long nCr(long long n, long long r) const {
        if (n < 0 || r < 0 || n < r) return 0;
        if (m == 1) return 0;
        
        // n が申告された n_max を超えた場合は安全のためエラー
        if (n > n_max_limit) {
            throw std::invalid_argument("エラー: クエリの n がコンストラクタで指定した n_max を超えています。");
        }
        
        unsigned __int128 res = 0; 
        for (size_t i = 0; i < M.size(); ++i) {
            unsigned long long b_i = bims[i].nCr(n, r);
            res += (unsigned __int128)b_i * crt_coeffs[i];
        }
        return (long long)(res % m);
    }
};
