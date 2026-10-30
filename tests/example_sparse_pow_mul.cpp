#include <bits/stdc++.h>
#include "fps/fps.hpp"
#include "sparse_pow_mul.hpp"
using namespace std;
using ll = long long;
#define rep(i,n) for(ll i=0;i<(ll)(n);++i)

void solve(){
    auto f = sparse_pow_mul(fps{1,1},2,fps{1,-1},-1,5);
    rep(i,f.size()) cout << f[i] << (i+1==(ll)f.size()?'\n':' ');
    // 1 3 4 4 4
}
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
}
