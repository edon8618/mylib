#include <bits/stdc++.h>
using namespace std;
#include "../ArbitraryModBinomial.txt"
#ifdef BINOMIAL_BENCH_COMPARE
namespace original {
#include "bench_reference/binomial_original.hpp"
}
namespace forced_wide {
#include "bench_reference/binomial_forced_wide.hpp"
}
#endif
template<class B> void run(const char*name,long long mod,int trials){
 vector<double> times;uint64_t last=0;
 for(int trial=0;trial<trials;trial++){
  B b(mod,200000);mt19937 rng(678);uint64_t checksum=0;auto start=chrono::steady_clock::now();
  for(int i=0;i<200000;i++){int n=rng()%200001,r=rng()%(n+1);checksum+=b.nCr(n,r);}
  double ms=chrono::duration<double,milli>(chrono::steady_clock::now()-start).count();times.push_back(ms);
  if(trial)assert(checksum==last);last=checksum;
 }
 sort(times.begin(),times.end());cout<<name<<" mod="<<mod<<" median_ms="<<fixed<<setprecision(3)<<times[times.size()/2]<<" checksum="<<last<<"\n";
}
int main(){
#ifdef BINOMIAL_BENCH_COMPARE
 run<original::ArbitraryModBinomial>("original32",998244353,3);
 run<forced_wide::ArbitraryModBinomial>("forced128",998244353,3);
#endif
 run<ArbitraryModBinomial>("hybrid32",998244353,3);
 run<ArbitraryModBinomial>("hybrid64",10460353203LL,3);
}
