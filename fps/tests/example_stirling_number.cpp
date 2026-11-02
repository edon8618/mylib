#include <bits/stdc++.h>
#include "fps/taylor_shift.hpp"
#include "stirling_number.hpp"
using namespace std;
using ll = long long;
#define rep(i,n) for(ll i=0;i<(ll)(n);++i)

void solve(){
    auto a = stirling2_n_fix(4);
    rep(i,a.size()) cout << a[i] << (i+1==(ll)a.size()?'\n':' ');
    // 0 1 7 6 1
    auto b = stirling1_k_fix(2,4);
    rep(i,b.size()) cout << b[i] << (i+1==(ll)b.size()?'\n':' ');
    // 0 0 1 3 11
}
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
}
