#include <bits/stdc++.h>
#include "fps/fps.hpp"
#include "product_of_polynomials.hpp"
using namespace std;
using ll = long long;
#define rep(i,n) for(ll i=0;i<(ll)(n);++i)

void solve(){
    vector<fps> a{{1,1},{1,2},{3}};
    auto f = product_of_polynomials(a);
    rep(i,f.size()) cout << f[i] << (i+1==(ll)f.size()?'\n':' ');
    // 3 9 6
}
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
}
