#include <bits/stdc++.h>
#include <sys/resource.h>
using namespace std;
#ifdef OLD
#include "old_LinkCutTree.hpp"
#else
#include "../../LinkCutTree.hpp"
#endif
constexpr int mod=998244353;
#ifdef AFFINE
struct S{int a,b;};
S op(S l,S r){return {int(1LL*r.a*l.a%mod),int((1LL*r.a*l.b+r.b)%mod)};}
S e(){return {1,0};} S value(int x){return {x%7,x%101};} long long hashval(S x){return x.a+ x.b;}
#else
using S=long long; S op(S a,S b){return a+b;} S e(){return 0;} S value(int x){return x%101;} long long hashval(S x){return x;}
#endif
int main(int argc,char**argv){
    int n=argc>1?atoi(argv[1]):200000,q=n,mode=argc>2?atoi(argv[2]):0;
    // 0 vertex/random tree; 1 symmetric edges; 2 directed edges;
    // 3 mixed vertices/directed edges; 4 vertex/long path; 5 vertex/link-cut.
    bool weighted=mode>=1&&mode<=3,verts=!weighted||mode==3;
    mt19937 rng(1234567);vector<int> parent(n),ids(n);
    for(int i=1;i<n;i++)parent[i]=mode==4?i-1:rng()%i;
    auto start=chrono::steady_clock::now(); auto cpu_start=clock();
#ifdef OLD
    if(mode==2||mode==3)return 2;
    LinkCutTree<S,op,e> t(weighted?2*n-1:n);
#else
    LinkCutTree<S,op,e> t(n);
#endif
    if(verts)for(int i=0;i<n;i++)t.set(i,value(i));
    for(int i=1;i<n;i++){
#ifdef OLD
        if(weighted){int x=n+i-1;t.set(x,value(i));t.link(parent[i],x);t.link(x,i);}
        else t.link(parent[i],i);
#else
        if(weighted)ids[i]=t.link_edge(parent[i],i,value(i),value(mode==1?i:i+13));
        else t.link(parent[i],i);
#endif
    }
    auto built=chrono::steady_clock::now(); auto cpu_built=clock();long long checksum=0;
    for(int j=0;j<q;j++){
        int u=rng()%n,v=rng()%n;
        if(mode==5&&j%4==0){int x=1+rng()%(n-1);t.cut(x,parent[x]);t.link(parent[x],x);}
        else if(j%4==0){
            if(weighted){int x=1+rng()%(n-1);
#ifdef OLD
                t.set(n+x-1,value(j));
#else
                t.set_edge(ids[x],value(j),value(mode==1?j:j+13));
#endif
            }else t.set(u,value(j));
        }else{
#ifdef OLD
            t.split(u,v);checksum^=hashval(t.get_prod(v));
#else
            checksum^=hashval(t.prod(u,v));
#endif
        }
    }
    auto done=chrono::steady_clock::now(); auto cpu_done=clock();rusage usage{};getrusage(RUSAGE_SELF,&usage);
    cout<<n<<','<<mode<<','<<chrono::duration<double>(built-start).count()<<','<<chrono::duration<double>(done-built).count()<<','<<usage.ru_maxrss<<','<<checksum<<','<<double(cpu_built-cpu_start)/CLOCKS_PER_SEC<<','<<double(cpu_done-cpu_built)/CLOCKS_PER_SEC<<'\n';
}
