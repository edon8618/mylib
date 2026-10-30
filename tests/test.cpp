#include <bits/stdc++.h>
#include "fps/taylor_shift.hpp"
#include "product_of_polynomials.hpp"
#include "stirling_number.hpp"
#include "sparse_pow_mul.hpp"
using namespace std;
using ll = long long;
#define rep(i,n) for(int i=0;i<(int)(n);++i)
int checks=0;
void eq(fps a,fps b){++checks; if(a!=b){cerr<<"FAIL "<<checks<<endl;abort();}}
fps mul(const fps&a,const fps&b,int d){fps c(d);rep(i,a.size())rep(j,b.size())if(i+j<d)c[i+j]+=a[i]*b[j];return c;}
fps inverse(const fps&a,int d){fps b(d);b[0]=a[0].inv();for(int i=1;i<d;++i){for(int j=1;j<=i&&j<(int)a.size();++j)b[i]-=a[j]*b[i-j];b[i]*=b[0];}return b;}
fps power(fps a,ll k,int d){if(k<0)a=inverse(a,d);unsigned long long e=k<0?0ULL-(unsigned long long)k:k;fps b(d);b[0]=1;while(e){if(e&1)b=mul(b,a,d);a=mul(a,a,d);e>>=1;}return b;}
int main(){
 mt19937 rng(123456);
 rep(t,2500){int d=1+rng()%25;fps a(rng()%13),b(rng()%13);for(auto&x:a)x=rng()%4;for(auto&x:b)x=rng()%4;ll n=(int)(rng()%11)-5,m=(int)(rng()%11)-5;if(n<0){if(a.empty())a.resize(1);a[0]=1+rng()%4;}if(m<0){if(b.empty())b.resize(1);b[0]=1+rng()%4;}eq(sparse_pow_mul(a,n,b,m,d),mul(power(a,n,d),power(b,m,d),d));}
 for(ll k:{LLONG_MIN,LLONG_MAX,-1LL,0LL,1LL})eq(sparse_pow_mul(fps{2,1},k,fps{3,0,1},-1,17),mul(power(fps{2,1},k,17),power(fps{3,0,1},-1,17),17));
 eq(sparse_pow_mul({},0,{1,1},1,3),fps{1,1,0});
 eq(sparse_pow_mul({0,1},LLONG_MAX,{1},1,3),fps(3));
 eq(sparse_pow_mul({1},1,{1},1,0),fps{});
 int bad=0;try{sparse_pow_mul({},-1,{0},1,3);}catch(const invalid_argument&){++bad;}try{sparse_pow_mul({1},0,{1},0,-1);}catch(const invalid_argument&){++bad;}assert(bad==2);
 rep(t,300){vector<fps>a(rng()%15);fps want{1};for(auto&f:a){f.resize(rng()%8);for(auto&x:f)x=rng()%6;if(want.empty()||f.empty())want={};else want=mul(want,f,want.size()+f.size()-1);}auto old=a;eq(product_of_polynomials(a),want);assert(a==old);}
 const int N=160;vector<fps>s1(N+1,fps(N+1)),s2=s1;s1[0][0]=s2[0][0]=1;
 for(int n=1;n<=N;++n)for(int k=1;k<=n;++k){s1[n][k]=s1[n-1][k-1]+(n-1)*s1[n-1][k];s2[n][k]=s2[n-1][k-1]+k*s2[n-1][k];}
 for(int n=0;n<=N;++n){fps a=s1[n],b=s2[n];a.resize(n+1);b.resize(n+1);eq(stirling1_n_fix(n),a);eq(stirling2_n_fix(n),b);rep(k,n+1)if((n-k)%2)a[k]=-a[k];eq(stirling1_n_fix(n,true),a);}
 for(int k=0;k<=N+2;++k){fps a(N+1),b(N+1);if(k<=N)rep(n,N+1){a[n]=s1[n][k];b[n]=s2[n][k];}eq(stirling1_k_fix(k,N),a);eq(stirling2_k_fix(k,N),b);rep(n,N+1)if((n-k)%2)a[n]=-a[n];eq(stirling1_k_fix(k,N,true),a);}
 // Larger inputs exercise qandori's NTT and FPS paths.
 vector<fps> factors(512,fps{1,1});fps a(513);a[0]=1;rep(i,512)a[i+1]=a[i]*(512-i)/mint(i+1);eq(product_of_polynomials(factors),a);
 eq(sparse_pow_mul({1,1},200,{1,1},300,1024),fps({1,1}).inv(1024).inv(1024).pow(500));
 cout<<checks<<" comparisons passed\n";
}
