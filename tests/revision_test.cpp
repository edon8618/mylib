#include <bits/stdc++.h>
#if __has_include(<boost/multiprecision/cpp_int.hpp>)
#include <boost/multiprecision/cpp_int.hpp>
#define HAVE_BOOST 1
#endif
using namespace std;using ll=long long;
#ifdef HAVE_BOOST
using boost::multiprecision::cpp_int;
#endif
#include "../WaveletMatrix.hpp"
#include "../binarytrie.hpp"
#include "../AVLtree.hpp"
#include "../ArbitraryModBinomial.hpp"
#include "../CRT.hpp"
#include "../floor ceil.hpp"
#include "../Diophantine.hpp"
#include "../DynamicLiChaoTree.hpp"
#include "../Tree/Tree.txt"
#include "../number-theory.hpp"
#include "../geometoric/geo_real.hpp"
namespace lca_euler {
#include "../LCA eulor tour + sparse table.hpp"
}
namespace lca_doubling {
#include "../LCA ダブリング.hpp"
}
namespace treap_lib {
#include "../implicit_treap.hpp"
}
struct TS {ll sum,len;bool operator==(const TS&b)const{return sum==b.sum&&len==b.len;}};
TS top(TS a,TS b){return {a.sum+b.sum,a.len+b.len};}TS te(){return {0,0};}TS tmap(ll f,TS a){return {a.sum+f*a.len,a.len};}ll tc(ll a,ll b){return a+b;}ll ti(){return 0;}TS tr(TS a){return a;}
using Treap=treap_lib::implicit_treap<TS,ll,top,te,tmap,tc,ti,tr>;
void wavelet(){
 mt19937 rng(198);
 for(int n=0;n<=130;n++){
  vector<int>a(n);for(int&i:a)i=rng()%4;WaveletMatrix<int,2>w(a);CompressedWaveletMatrix<int,2>c(a);RangeDistinctQuery<int>rd(a);
  for(int i=0;i<n;i++)assert(w.access(i)==a[i]);
  for(int q=0;q<200;q++){int l=rng()%(n+1),r=rng()%(n+1);if(l>r)swap(l,r);int bound=int(rng()%10)-2;int cnt=0,exact=0;set<int> ss;
   for(int i=l;i<r;i++){cnt+=a[i]<bound;exact+=a[i]==bound;ss.insert(a[i]);}
   assert(w.rank_less_than(l,r,bound)==cnt);assert(c.rank_less_than(l,r,bound)==cnt);assert(w.rank(l,r,bound)==exact);assert(rd.query(l,r)==(int)ss.size());
   if(l<r){auto b=vector<int>(a.begin()+l,a.begin()+r);sort(b.begin(),b.end());int k=rng()%b.size();assert(w.kth_smallest(l,r,k)==b[k]);}
  }
  BitVector bv(n);for(int i=0;i<n;i++)if(a[i]&1)bv.set(i);bv.build();int cnt=0;for(int i=0;i<=n;i++){assert(bv.rank1(i)==cnt);if(i<n)cnt+=a[i]&1;}
 }
 BinaryTrie<int,2> t;vector<int>a;for(int i=0;i<2000;i++){int x=rng()%4;if(a.empty()||rng()%2){t.insert(x);a.push_back(x);}else{int j=rng()%a.size();t.erase(a[j]);a.erase(a.begin()+j);}
  for(int z=0;z<4;z++)for(int limit=-2;limit<=6;limit++){int cnt=0;for(auto v:a)cnt+=(v^z)<limit;assert(t.count_less_than(z,limit)==cnt);}assert(t.count_exact(4)==0);
 }
 cerr<<"B02 B03 PASS\n";
}
void move_ownership(){
 static_assert(!is_copy_constructible_v<AVLtree<int>>);static_assert(!is_copy_assignable_v<AVLtree<int>>);
 static_assert(!is_copy_constructible_v<Treap>);static_assert(!is_copy_assignable_v<Treap>);
 AVLtree<int> a;vector<int>v;mt19937 rng(18);
 for(int i=0;i<5000;i++){int x=rng()%50;if(v.empty()||rng()%2){a.insert(x);v.insert(upper_bound(v.begin(),v.end(),x),x);}else{auto p=lower_bound(v.begin(),v.end(),x);int expected=(p==v.end()||*p!=x)?-1:p-v.begin();assert(a.erase(x)==expected);if(expected!=-1)v.erase(p);}
  assert(a.count(x)==count(v.begin(),v.end(),x));assert(a.to_vector()==v);
 }
 AVLtree<int>b(std::move(a));assert(a.empty());assert(b.to_vector()==v);a.insert(99);AVLtree<int>c;c.insert(-1);c=std::move(b);assert(b.empty());assert(c.to_vector()==v);c=std::move(c);assert(c.to_vector()==v);
 vector<AVLtree<int>> containers;containers.push_back(std::move(c));assert(containers[0].to_vector()==v);
 Treap s;s.build({{1,1},{2,1},{3,1}});s.apply(0,3,2);s.range_reverse(0,3);
 Treap u(std::move(s));assert(s.size()==0);assert(u.get(0).sum==5);s.insert(0,{9,1});Treap z;z.insert(0,{-1,1});z=std::move(u);assert(u.size()==0);assert(z.prod(0,3).sum==12);z=std::move(z);assert(z.size()==3);
 vector<Treap> ts;ts.push_back(std::move(z));assert(ts[0].prod(0,3).sum==12);
 cerr<<"B04 B16 PASS\n";
}
void number_cases(){
 vector<ll> bounds={LLONG_MIN,LLONG_MIN+1,-100,-3,-1,0,1,3,100,LLONG_MAX-1,LLONG_MAX};
 for(ll a:bounds)for(ll b:bounds){if(b==0)continue;if(a==LLONG_MIN&&b==-1){bool caught=false;try{::floor(a,b);}catch(overflow_error&){caught=true;}assert(caught);continue;}
  __int128 q=(__int128)a/b,r=(__int128)a%b;ll f=q-(r!=0&&((r>0)!=(b>0))),c=q+(r!=0&&((r>0)==(b>0)));assert(::floor(a,b)==f);assert(::ceil(a,b)==c);
 }
 for(ll a=-5;a<=5;a++)for(ll b=-5;b<=5;b++)for(ll c=-20;c<=20;c++){
  ll expected=0;for(ll x=-4;x<=4;x++)for(ll y=-4;y<=4;y++)expected+=a*x+b*y==c;assert(Diophantine::count(a,b,c,-4,4,-4,4)==expected);
 }
 assert(Diophantine::count(LLONG_MIN,1,LLONG_MIN,1,1,0,0)==1);
 bool caught=false;try{Diophantine::count(0,0,0,LLONG_MIN,LLONG_MAX,LLONG_MIN,LLONG_MAX);}catch(overflow_error&){caught=true;}assert(caught);
 for(ll m1=1;m1<35;m1++)for(ll m2=1;m2<35;m2++)for(ll r1=-3;r1<4;r1++)for(ll r2=-3;r2<4;r2++){
  auto [x,m]=CRT(r1,m1,r2,m2);ll expected=-1;for(ll v=0;v<lcm(m1,m2);v++)if(v%m1==safe_mod(r1,m1)&&v%m2==safe_mod(r2,m2)){expected=v;break;}
  assert(expected<0?(m==0):(x==expected&&m==lcm(m1,m2)));
 }
 caught=false;try{CRT(0LL,2LL,1000000000000000000LL,4000000000000000001LL);}catch(overflow_error&){caught=true;}assert(caught);
 auto wide=CRT(__int128(0),__int128(2),__int128(1000000000000000000LL),__int128(4000000000000000001LL));assert(wide.first==__int128(1000000000000000000LL));
 #ifdef HAVE_BOOST
 cpp_int p=(cpp_int(1)<<160)+1,q=(cpp_int(1)<<170)+3;auto arbitrary=CRT(cpp_int(7),p,cpp_int(7),q);assert(arbitrary.first==7);assert(arbitrary.second%p==0&&arbitrary.second%q==0);
 #endif
 caught=false;try{CRT((__int128)0,(__int128(1)<<100),(__int128)0,(__int128(1)<<100)+1);}catch(overflow_error&){caught=true;}assert(caught);
 assert(nt::pow_mod(2,0,1)==0);assert(nt::pow_mod(7,0,13)==1);assert(nt::pow_mod(7,3,13)==5);
 cerr<<"B12 B13 B17 B20 PASS\n";
}
void lca_tree(){
 mt19937 rng(87);for(int n=1;n<=70;n++){
  vector<vector<ll>>g(n);vector<int>par(n,-1),depth(n);Tree<ll>t(n);Tree<double>real(n);vector<ll>xs(n);vector<double>sum(n);
  for(int i=1;i<n;i++){par[i]=rng()%i;depth[i]=depth[par[i]]+1;g[i].push_back(par[i]);g[par[i]].push_back(i);ll w=rng()%128;t.add_edge(i,par[i],w);real.add_edge(i,par[i],w/2.0);xs[i]=xs[par[i]]^w;sum[i]=sum[par[i]]+w/2.0;}
  lca_euler::LCA a(g);lca_doubling::LCA b(g);t.build_lca();t.build_xor();real.build_lca();real.build_lca();
  for(int v=0;v<n;v++)for(int k=0;k<=n+10;k++){int expected=v;for(int j=0;j<k&&expected>=0;j++)expected=par[expected];assert(a.kth_ancestor(v,k)==expected);assert(b.kth_ancestor(v,k)==expected);}
  for(int u=0;u<n;u++)for(int v=0;v<n;v++){int x=u,y=v;while(depth[x]>depth[y])x=par[x];while(depth[y]>depth[x])y=par[y];while(x!=y){x=par[x];y=par[y];}assert(real.get_lca(u,v)==x);assert(real.get_dist(u,v)==sum[u]+sum[v]-2*sum[x]);assert(t.get_xor_dist(u,v)==(xs[u]^xs[v]));}
 }
 cerr<<"B14 B19 PASS\n";
}
void binomial(){
 #ifdef HAVE_BOOST
 for(ll mod: {1LL,2LL,4LL,12LL,99991LL,998244353LL,10460353203LL,10000000019LL}){
  ArbitraryModBinomial b(mod,150);vector<cpp_int> row(152);row[0]=1;
  for(int n=0;n<=150;n++){for(int r=0;r<=n;r++)assert(b.nCr(n,r)==(row[r]%mod).convert_to<ll>());for(int r=n+1;r>0;r--)row[r]+=row[r-1];}
 }
 #else
 for(ll mod: {1LL,2LL,4LL,12LL,99991LL,998244353LL,10460353203LL,10000000019LL}){
  ArbitraryModBinomial b(mod,150);vector<ll>row(152);row[0]=1%mod;
  for(int n=0;n<=150;n++){for(int r=0;r<=n;r++)assert(b.nCr(n,r)==row[r]);for(int r=n+1;r>0;r--)row[r]=(row[r]+row[r-1])%mod;}
 }
 #endif
 cerr<<"B11 PASS\n";
}
void tangent(){
 using namespace geo_real;vector<pair<Circle,Circle>> cs={{{{0,0},1},{{2,0},1}},{{{0,0},2},{{1,0},1}},{{{0,0},1},{{4,0},1}},{{{0,0},2},{{3,4},1}}};
 for(auto [a,b]:cs){auto lines=commonTangents(a,b);assert(!lines.empty());for(auto l:lines){assert(norm2(l.b-l.a)>0);assert(fabsl(distance(l,a.c)-a.r)<1e-9L);assert(fabsl(distance(l,b.c)-b.r)<1e-9L);}}
 cerr<<"B10 PASS\n";
}
void li_chao(){
 mt19937 rng(11);for(int n=1;n<55;n++){
  ll lo=-n,hi=n;DynamicLiChaoTree<ll>t(lo,hi,LLONG_MIN);struct Line{ll m,c,l,r;};vector<Line>lines;
  for(int k=0;k<300;k++){ll m=int(rng()%41)-20,c=int(rng()%41)-20,a=int(rng()%(2*n+7))-n-3,b=int(rng()%(2*n+7))-n-3;if(a>b)swap(a,b);if(rng()%2){t.add_line(m,c);lines.push_back({m,c,lo,hi+1});}else{t.add_segment(m,c,a,b);lines.push_back({m,c,a,b});}
   for(ll x=lo;x<=hi;x++){ll best=LLONG_MIN;for(auto f:lines)if(f.l<=x&&x<f.r)best=max(best,f.m*x+f.c);assert(t.query(x)==best);}
  }
 }
 DynamicLiChaoTree<ll>extreme(LLONG_MIN,LLONG_MAX,LLONG_MIN);extreme.add_line(0,7);assert(extreme.query(LLONG_MAX)==7);extreme.add_segment(0,9,LLONG_MIN,LLONG_MIN+1);assert(extreme.query(LLONG_MIN)==9);assert(extreme.query(LLONG_MIN+1)==7);extreme.add_line(LLONG_MAX,0);bool caught=false;try{extreme.query(LLONG_MAX);}catch(overflow_error&){caught=true;}assert(caught);
 for(int iter=0;iter<30;iter++){
  DynamicLiChaoTree<double>t(-10,10,-1e300);struct L{double m,c,l,r;};vector<L>ls;
  for(int k=0;k<100;k++){double m=int(rng()%41)-20,c=int(rng()%41)-20,l=(int(rng()%81)-40)/4.0,r=(int(rng()%81)-40)/4.0;if(l>r)swap(l,r);if(rng()%2){t.add_line(m,c);ls.push_back({m,c,-10,11});}else{t.add_segment(m,c,l,r);ls.push_back({m,c,l,r});}
   for(int j=0;j<=80;j++){double x=-10+j/4.0,best=-1e300;for(auto f:ls)if(f.l<=x&&x<f.r)best=max(best,f.m*x+f.c);assert(t.query(x)==best);}
  }
 }
 double a=1e18,b=nextafter(a,numeric_limits<double>::infinity());DynamicLiChaoTree<double>adj(a,b,-1e300);adj.add_segment(0,5,a,b);assert(adj.query(a)==5);assert(adj.query(b)==-1e300);adj.add_line(0,3);assert(adj.query(b)==3);
 DynamicLiChaoTree<double>sub(-1,1,-1e300);sub.add_segment(0,2,0,numeric_limits<double>::denorm_min());assert(sub.query(0)==2);assert(sub.query(numeric_limits<double>::denorm_min())==-1e300);
 cerr<<"B15 PASS\n";
}
int main(){wavelet();move_ownership();number_cases();lca_tree();binomial();tangent();li_chao();cerr<<"revision PASS\n";}
