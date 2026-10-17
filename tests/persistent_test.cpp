#include <bits/stdc++.h>
using namespace std;
using ll=long long;
ll sop(ll x,ll y){return x+y;}
ll se(){return 0;}
string cop(string a,string b){return a+b;}
string ce(){return "";}
#include "../persistant_segtree.txt"
bool threshold_nonzero(ll x){return x<=100000;}
int main(){
 mt19937 rng(128759);
 for(int n=0;n<=80;n++) {
  vector<ll> a(n);for(ll &x:a)x=rng()%8;
  persistent_segtree<ll,sop,se> t(a);
  vector<vector<ll>> vers{a};
  for(int iter=0;iter<50;iter++) {
   int ver=rng()%vers.size();
   if(n && iter%2==0){int idx=rng()%n;ll val=rng()%10;vector<ll> b=vers[ver]; b[idx]=val;int newid=t.set(idx,val,ver);assert(newid==(int)vers.size());vers.push_back(b);ver=newid;}
   const auto&v=vers[ver];
   ll cap=rng()%150;
   auto predicate=[&](ll sm){return sm<=cap;};
   for(int l=0;l<=n;l++) {
    int actual=t.max_right(l,predicate,ver),expected=l;ll sm=0;
    while(expected<n && sm+v[expected]<=cap)sm+=v[expected++];
    if(actual!=expected){cerr<<"MR n="<<n<<" l="<<l<<" ver="<<ver<<" cap="<<cap<<" actual="<<actual<<" expected="<<expected<<endl;return 1;}
   }
   for(int r=0;r<=n;r++) {
    int actual=t.min_left(r,predicate,ver),expected=r;ll sm=0;
    while(expected>0 && sm+v[expected-1]<=cap)sm+=v[--expected];
    if(actual!=expected){cerr<<"ML n="<<n<<" r="<<r<<" ver="<<ver<<" cap="<<cap<<" actual="<<actual<<" expected="<<expected<<endl;return 1;}
   }
   for(int q=0;q<20;q++){int l=n?(rng()%(n+1)):0;int r=n?(rng()%(n+1)):0;if(l>r)swap(l,r);ll expect=accumulate(v.begin()+l,v.begin()+r,0LL);assert(t.prod(l,r,ver)==expect);}
  }
 }
 // Verify preservation of non-commutative monoid order via prefix/suffix predicates.
 for(int n=0;n<=14;n++){
  vector<string> a(n);for(auto &s:a)s=string(1,char('a'+rng()%3));
  persistent_segtree<string,cop,ce> t(a);
  for(int l=0;l<=n;l++)for(int r=l;r<=n;r++){
   string target="";for(int j=l;j<r;j++)target+=a[j];
   auto pf=[&](string s){return s.size()<=target.size() && target.substr(0,s.size())==s;};
   // if next char equals target at boundary, no more target chars: exact r.
   int mr=t.max_right(l,pf);if(mr!=r){cerr<<"noncom MR "<<l<<" "<<r<<" "<<mr<<endl;return 1;}
   string st="";for(int j=l;j<r;j++)st+=a[j];
   auto sf=[&](string s){return s.size()<=st.size() && st.substr(st.size()-s.size())==s;};
   int ml=t.min_left(r,sf);if(ml!=l){cerr<<"noncom ML "<<l<<" "<<r<<" "<<ml<<endl;return 1;}
  }
 }
 // Exercise function-pointer overload, not just lambda overload.
 persistent_segtree<ll,sop,se> t(vector<ll>{1,1,1,1});
 assert(t.max_right<threshold_nonzero>(0)==4);
 assert(t.min_left<threshold_nonzero>(4)==0);
 cerr<<"persistent PASS\n";
}
