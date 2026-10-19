// #pragma GCC optimize("O3,unroll-loops")
// #pragma GCC target("avx2")

#include <bits/stdc++.h>
using namespace std;
#include <atcoder/all>
using namespace atcoder;



// #include <boost/multiprecision/cpp_int.hpp>
// using namespace boost::multiprecision;

#define ll long long
#define ld long double
#define rep(i, n) for (ll i = 0; i < (ll)(n); ++i)
#define vi vector<int>
#define vl vector<ll>
#define vd vector<double>
// #define vb vector<bool>
#define vs vector<string>
#define vc vector<char>
#define ull unsigned long long
#define all(a) (a).begin(), (a).end()
#define rall(a) (a).rbegin(), (a).rend()

template<class T, class U>
inline bool chmax(T &a, const U &b) {
    if (a < b) {
        a = b;
        return true;
    }
    return false;
}

template<class T, class U>
inline bool chmin(T &a, const U &b) {
    if (a > b) {
        a = b;
        return true;
    }
    return false;
}


// #define ll int
// #define ll int128_t
// #define ll int256_t
// #define ll cpp_int


constexpr ll inf = (1ll << 60);
// constexpr ll inf = (1 << 30);
// const double PI=3.1415926535897932384626433832795028841971;

uint32_t xor_x = 123456789, xor_y = 362436069, xor_z = 521288629, xor_w = 88675123;
inline uint32_t xor_next() {
    uint32_t t = xor_x ^ (xor_x << 11);
    xor_x = xor_y; xor_y = xor_z; xor_z = xor_w;
    return xor_w = (xor_w ^ (xor_w >> 19)) ^ (t ^ (t >> 8));
}
inline int rnd(int max_val) { return xor_next() % max_val; }


struct Timer {
    std::chrono::steady_clock::time_point start_time;

    Timer() {
        reset();
    }

    // 測定の起点リセット用
    void reset() {
        start_time = std::chrono::steady_clock::now();
    }

    // スタートからの経過時間をミリ秒(msec)で返す
    long long get_ms() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count();
    }
};

struct I128Hash {
    size_t operator()(__int128_t a) const {
        uint64_t x = (uint64_t)a ^ (uint64_t)((__uint128_t)a>>64);
        x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27; x *= 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
};



// vl fact;
// ll kai(ll n){
//     fact.resize(n,1);
//     rep(i,n-1)fact[i+1]=fact[i]*(i+1);
// }

// using mint = ld;
using mint = modint998244353;//static_modint<998244353>
// using mint = modint1000000007;//static_modint<1000000007>
// using mint = static_modint<922267487>;   // 多分落とされにくい NOT ntt-friendly
// using mint = static_modint<469762049>;   // ntt-friendly
// using mint = static_modint<167772161>;   // ntt-friendly
// using mint = modint;//mint::set_mod(mod);

// ll const mod=1000000007ll;
// ll const mod=998244353ll;
// ll modrui(ll a,ll b,ll mod){
//     a%=mod;
//     if(b==0)return 1;
//     if(b%2==1) return a*modrui(a*a%mod,b/2,mod)%mod;
//     return modrui(a*a%mod,b/2,mod)%mod;
// }

// void incr(vl &v,ll n){//　n進法
//     ll k=v.size();
//     v[k-1]++;
//     ll now=k-1;
//     while (v[now]>=n)
//     {
//         v[now]=0;
//         if(now==0)break;
//         v[now-1]++;
//         now--;
//     }
//     return;
// }

// vector<mint> fact,invf;
// void init_modfact(ll sz){
//     fact.resize(sz);
//     invf.resize(sz);
//     fact[0]=1;
//     rep(i,sz-1){
//         fact[i+1]=fact[i]*(i+1);
//     }
//     invf[sz-1]=1/fact[sz-1];
//     for(ll i=sz-2; i>=0; i--){
//         invf[i]=invf[i+1]*(i+1);
//     }
// }
// mint choose(ll n,ll r){
//     if(n<r || r<0 || n<0)return 0;
//     return fact[n]*invf[r]*invf[n-r];
// }

// vector<mint> modpow,invpow;
// void init_modpow(ll x,ll sz){
//     mint inv=1/mint(x);
//     modpow.assign(sz,1);
//     invpow.assign(sz,1);
//     rep(i,sz-1){
//         modpow[i+1]=modpow[i]*x;
//         invpow[i+1]=invpow[i]*inv;
//     }
// }
// long long phi(long long n) {// O(sqrt(n))
//     long long res = n;
//     for (long long i = 2; i * i <= n; i++) {
//         if (n % i == 0) {
//             res -= res / i;
//             while (n % i == 0) n /= i;
//         }
//     }
//     if (n > 1) res -= res / n;
//     return res;
// }



void solve(){

}




int main(){
    ios::sync_with_stdio(false);
    std::cin.tie(nullptr);


    // ll mx=12345678;
    // vc fl(mx+1,0);
    // for(ll d=2;d<=mx;d++){
    //     if(fl[d])continue;
    //     ll x=d;
    //     ps.push_back(x);
    //     while(x<=mx){
    //         fl[x]=1;
    //         x+=d;
    //     }
    // }


    ll t = 1;
    // cin >> t;
    while (t--){
        solve();
        cout << flush;
    }
}