#include <bits/stdc++.h>
#include <atcoder/modint>
#include <atcoder/convolution>
using namespace std;using mint=atcoder::modint998244353;
#include "../FPS/poly.txt"
#include "../FPS/sparse_pow_mul.txt"
#include "../FPS/multipoint_evaluation_geometric_.txt"
#include "../FPS/Lagrange_Inversion_Formula.txt"
#include "../FPS/multipoint_evaluation (middle_product).txt"
#include "../FPS/polynomial_interpolation.txt"
#include "../FPS/stirling_number.txt"
FPS naive_mul(const FPS&a,const FPS&b,int size){FPS c(size);for(int i=0;i<(int)a.size();i++)for(int j=0;j<(int)b.size()&&i+j<size;j++)c[i+j]+=a[i]*b[j];return c;}
FPS naive_pow(FPS a,int n,int size){FPS r(size);r[0]=1;while(n--)r=naive_mul(r,a,size);return r;}
int main(){
 assert(sparse_pow_mul(FPS{0},0,FPS{1,1},1,2)==FPS({1,1,0}));
 assert(sparse_pow_mul(FPS{0},0,FPS{0},0,0)==FPS({1}));
 for(int m=0;m<20;m++)for(mint r: {mint(0),mint(1),mint(2)}) assert(multipoint_evaluation_geometric({},1,r,m)==vector<mint>(m,0));
 for(int d=0;d<20;d++){assert(FPS{1}.inv(d).size()==size_t(d+1));assert(FPS{1}.dense_inv(d).size()==size_t(d+1));assert(FPS{1}.sparse_inv(d).size()==size_t(d+1));
  FPS f(2*d+3);f.back()=1;assert(f.sqrt(d)==FPS(d+1,0));assert(f.dense_sqrt(d)==FPS(d+1,0));
 }
 FPS dense(75,1);bool caught=false;try{dense.pow(-1,3);}catch(invalid_argument&){caught=true;}assert(caught);
 caught=false;try{FPS{1,1}.pow(-1,3);}catch(invalid_argument&){caught=true;}assert(caught);
 assert(FPS({1,1}).inv(2).pow(2,2)==FPS({1,-2,3}));
 mt19937 rng(567);int checks=0;
 for(int d=0;d<=32;d++)for(int trial=0;trial<10;trial++){
  FPS f(1+rng()%10),g(1+rng()%10);for(auto&x:f)x=rng()%5;for(auto&x:g)x=rng()%5;int n=rng()%5,m=rng()%5;
  auto ex=naive_mul(naive_pow(f,n,d+1),naive_pow(g,m,d+1),d+1);assert(sparse_pow_mul(f,n,g,m,d)==ex);checks++;
  int index=rng()%(d+1);auto lif=LIF(f,g,index,d);assert(lif.size()==size_t(d+1));FPS power(d+1);power[0]=1;
  for(int k=0;k<=d;k++){assert(lif[k]==naive_mul(f,power,d+1)[index]);power=naive_mul(power,g,d+1);}checks++;
  vector<mint>xs(d+1);for(auto&x:xs)x=rng()%50;auto ev=multipoint_evaluation(f,xs);for(int i=0;i<=d;i++){mint value=0;for(int j=(int)f.size()-1;j>=0;j--)value=value*xs[i]+f[j];assert(ev[i]==value);}checks++;
 }
 // Dense paths and the real ACL NTT path (min operand length > 60).
 for(int size: {75,128,257,1024}){
  FPS f(size);f[0]=1;for(int i=1;i<size;i++)f[i]=rng()%20+1;
  int d=size-1;auto inverse=f.inv(d);assert(naive_mul(f,inverse,size)==FPS({1}).pre(size));assert(f.pow(2,d)==naive_mul(f,f,size));
  auto log=f.log(d);assert(log.exp(d)==f);auto root=naive_mul(f,f,size).sqrt(d);assert(naive_mul(root,root,size)==naive_mul(f,f,size));checks+=4;
 }
 // n/k fixed Stirling callers must retain their old sequence length.
 for(int n=0;n<20;n++)for(int k=0;k<=n;k++){auto s1=stirling1_n_fix(n),s2=stirling2_n_fix(n);assert(s1[k]==stirling1_k_fix(k,n)[n]);assert(s2[k]==stirling2_k_fix(k,n)[n]);}
 cerr<<"FPS extended PASS "<<checks<<" checks; actual ACL\n";
}
