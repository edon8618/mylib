#include <bits/stdc++.h>
using namespace std;
#include <atcoder/modint>
#include <atcoder/convolution>
using mint=atcoder::modint998244353;
#include "../FPS/poly.txt"
#include "../FPS/stirling_number.txt"
#include "../FPS/get_nth_coef_of_inverse.txt"
#include "../FPS/multipoint_evaluation (middle_product).txt"
int main(int argc,char**argv){
 if(argc>1&&string(argv[1])=="negative"){
  FPS f{1,2}; try { f.pow(-1,4); } catch (const invalid_argument&) { cerr<<"negative pow rejected PASS\n"; return 0; } return 1;
 }
 assert((FPS{0,1}.compositional_inverse(0)==FPS{0}));
 mt19937 rng(1234);
 for(int deg=0;deg<16;deg++){
  FPS f{1,1};
  FPS iv=f.inv(deg);assert((int)iv.size()==deg+1);
  assert((f*iv).pre(deg+1)==FPS({1}).pre(deg+1));
  FPS sq=f.pow(2,deg);assert((int)sq.size()==deg+1);
  assert(sq==(f*f).pre(deg+1));
  assert(f.pow(0,deg)==FPS({1}).pre(deg+1));
  FPS lg=f.log(deg);assert((int)lg.size()==deg+1);
  assert(lg.exp(deg)==f.pre(deg+1));
  FPS root=(f*f).sqrt(deg);assert((int)root.size()==deg+1);assert((root*root).pre(deg+1)==(f*f).pre(deg+1));
  assert(f.dense_inv(deg)==iv);
  assert(f.sparse_inv(deg)==iv);
  assert(f.dense_pow(2,deg)==sq);
  assert(f.sparse_pow(2,deg)==sq);
 }
 assert(FPS{1}.inv(0)==FPS({1}));
 assert((FPS{0,0,1}.sqrt(0)==FPS({0})));
 assert((FPS{0,0,1}.sqrt(1)==FPS({0,1})));
 FPS zeros;assert(zeros.pow(0,4)==FPS({1}).pre(5));
 // n > 60 triggers dense dispatch.
 FPS f(75);f[0]=1;for(int i=1;i<75;i++)f[i]=rng()%17+1;
 for(int deg: {0,1,3,16,63,70}){
  FPS iv=f.inv(deg);assert(iv==f.dense_inv(deg));
  assert((f*iv).pre(deg+1)==FPS({1}).pre(deg+1));
  FPS sq=f.pow(2,deg);assert(sq==(f*f).pre(deg+1));
  assert(f.pow(0,deg)==FPS({1}).pre(deg+1));
 }
 for(int n=1;n<=9;n++){
  FPS f(n+1); f[0]=0;f[1]=1;for(int i=2;i<=n;i++)f[i]=rng()%11;
  auto g=f.compositional_inverse(n);assert((int)g.size()==n+1);
  auto fg=f.composition(g,n);assert(fg==FPS({0,1}).pre(n+1));
  assert(get_nth_coef_of_inverse(f,n)==g[n]);
 }
 for(int n=0;n<=15;n++){
  auto a=stirling1_k_fix(3,n);auto b=stirling2_k_fix(3,n);
  assert((int)a.size()==n+1);assert((int)b.size()==n+1);
 }
 cerr<<"fps PASS\n";
}
