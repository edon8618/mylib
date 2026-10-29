#pragma once
//http://atcoder.jp/contests/abc236/submissions/73613608
//https://judge.yosupo.jp/submission/355510
template <class S>
struct Matrix {
    using T = typename S::T;
    int h, w;
    std::vector<std::vector<T>> val;

    // zero() で初期化
    Matrix(int h, int w) : h(h), w(w), val(h, std::vector<T>(w, S::zero())) {}

    // 単位行列の生成（対角成分を one() にする）
    static Matrix eye(int n) {
        Matrix res(n, n);
        for (int i = 0; i < n; i++) res.val[i][i] = S::one();
        return res;
    }

    // 行列同士の積
    Matrix operator*(const Matrix& a) const {
        assert(w == a.h);
        Matrix res(h, a.w);
        for (int i = 0; i < h; i++) {
            for (int k = 0; k < w; k++) {
                for (int j = 0; j < a.w; j++) {
                    // res[i][j] += val[i][k] * a[k][j] の抽象化
                    res.val[i][j] = S::add(res.val[i][j], S::mul(val[i][k], a.val[k][j]));
                }
            }
        }
        return res;
    }

    // ベクトルへの線形変換 (Matrix * Vector)
    std::vector<T> operator*(const std::vector<T>& vec) const {
        assert(w == vec.size());
        std::vector<T> res(h, S::zero());
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                res[i] = S::add(res[i], S::mul(val[i][j], vec[j]));
            }
        }
        return res;
    }

    // 行列の累乗（ダブリング）
    Matrix pow(long long k) const {
        assert(h == w);
        Matrix res = eye(h);
        Matrix a = *this;
        while (k > 0) {
            if (k & 1) res = res * a;
            a = a * a;
            k >>= 1;
        }
        return res;
    }
};
