#include <bits/stdc++.h>
using namespace std;using u64=unsigned long long;using u128=__uint128_t;
struct Barrett{u64 m,im;explicit Barrett(u64 m):m(m),im(-1ULL/m+1){}u64 mul(u64 a,u64 b)const{u64 z=a*b,x=u128(z)*im>>64,v=z-x*m;if(m<=v)v+=m;return v;}};
int main(int argc,char**argv){u64 m=argc>1?stoull(argv[1]):998244353;assert(m>0 && m<=numeric_limits<unsigned int>::max());Barrett b(m);vector<pair<u64,u64>>v(2000000);mt19937_64 rng(567);for(auto&[x,y]:v){x=rng()%m;y=rng()%m;}vector<double>fast,wide;u64 last=0;
for(int k=0;k<3;k++){for(int mode=0;mode<2;mode++){u64 sum=0;auto t=chrono::steady_clock::now();for(auto[x,y]:v)sum+=mode?u64(u128(x)*y%m):b.mul(x,y);double ms=chrono::duration<double,milli>(chrono::steady_clock::now()-t).count();(mode?wide:fast).push_back(ms);if(k||mode)assert(sum==last);last=sum;}}sort(fast.begin(),fast.end());sort(wide.begin(),wide.end());cout<<"modmul mod="<<m<<" operations=2000000 barrett_ms="<<fast[1]<<" wide_ms="<<wide[1]<<" checksum="<<last<<"\n";}
