#include <bits/stdc++.h>
#include "../../LinkCutTree.hpp"
using namespace std;
constexpr int MOD=998244353;
struct Aff { int a,b; bool operator==(const Aff&) const = default; };
Aff compose(Aff l,Aff r){return {int(1LL*r.a*l.a%MOD),int((1LL*r.a*l.b+r.b)%MOD)};}
Aff unit(){return {1,0};}
int add(int a,int b){return a+b;} int zero(){return 0;}
template<class S,S(*op)(S,S),S(*e)(),class Gen>
void random_test(unsigned seed,Gen gen){
    mt19937 rng(seed); const int n=25;
    LinkCutTree<S,op,e> t(n);
    vector<S> val(n,e());
    struct E{int u,v; S f,g; bool alive;}; vector<E> es;
    vector<vector<int>> adj(n,vector<int>(n,-1));
    auto route=[&](int u,int v){
        vector<int> p(n,-1),q{u}; p[u]=u;
        for(size_t i=0;i<q.size();i++)for(int j=0;j<n;j++)if(adj[q[i]][j]>=0&&p[j]<0){p[j]=q[i];q.push_back(j);}
        vector<int> a; if(p[v]<0)return a;
        for(int x=v;;x=p[x]){a.push_back(x);if(x==u)break;} reverse(a.begin(),a.end());return a;
    };
    for(int step=0;step<10000;step++){
        int u=rng()%n,v=rng()%n,action=rng()%12;
        auto p=route(u,v);
        if(action<2&&p.empty()){
            S f=gen(rng),g=gen(rng); int id;
            if(action==0){t.link(u,v);id=t.get_edge_id(v,u);f=g=e();}
            else id=t.link_edge(u,v,f,g);
            assert(id==int(es.size())); es.push_back({u,v,f,g,true});adj[u][v]=adj[v][u]=id;
        }else if(action==2){
            int id=adj[u][v];t.cut(v,u);
            if(id>=0){es[id].alive=false;adj[u][v]=adj[v][u]=-1;assert(t.get_edge_id(u,v)==-1);}
        }else if(action==3){val[u]=gen(rng);t.set(u,val[u]);assert(t.get(u)==val[u]);}
        else if(action<=6&&!es.empty()){
            int id=rng()%es.size();auto& x=es[id];if(!x.alive)continue;
            if(action==4){t.cut_edge(id);x.alive=false;adj[x.u][x.v]=adj[x.v][x.u]=-1;continue;}
            S f=gen(rng),g=gen(rng);
            if(action==6 && rng()%2){t.set_edge(id,f);x.f=x.g=f;}
            else if(action==5){t.set_edge(x.v,x.u,f,g);x.f=g;x.g=f;}
            else{t.set_edge(id,f,g);x.f=f;x.g=g;}
            assert((t.get_edge(id)==pair<S,S>{x.f,x.g}));
            assert((t.get_edge(x.v,x.u)==pair<S,S>{x.g,x.f}));
        }
        p=route(u,v);
        assert(t.get_edge_id(u,v)==adj[u][v]);
        assert(t.is_connected(u,v)==!p.empty());
        auto result=t.try_prod(u,v);
        assert(bool(result)==!p.empty());
        int r=t.find_root(u);assert(0<=r&&r<n&&!route(u,r).empty());
        if(p.empty())continue;
        S want=val[u];
        for(int i=1;i<int(p.size());i++){
            auto& x=es[adj[p[i-1]][p[i]]];
            want=op(op(want,x.u==p[i-1]?x.f:x.g),val[p[i]]);
        }
        assert(*result==want); assert(t.prod(u,v)==want);
        assert(t.vertex_count(u,v)==int(p.size()));assert(t.edge_count(u,v)==int(p.size())-1);
        int k=rng()%p.size();assert(t.kth_vertex(u,v,k)==p[k]);
        if(p.size()>1){k=rng()%(p.size()-1);assert(t.kth_edge(u,v,k)==adj[p[k]][p[k+1]]);}
        t.split(u,v);assert(t.get_prod(v)==want);assert(t.get_sz(v)==int(p.size()));
        k=rng()%p.size();assert(t.get_kth(v,k)==p[k]);
        if(step%37==0)for(int id=0;id<int(es.size());id++)if(es[id].alive){
            auto x=es[id];assert((t.get_edge(id)==pair<S,S>{x.f,x.g}));
        }
    }
}
void directed(){
    LinkCutTree<Aff,compose,unit> t(3);
    int id=t.link_edge(0,1,{2,3},{0,7});t.link(1,2);
    for(int i=0;i<100;i++){assert((t.prod(0,1)==Aff{2,3}));assert((t.prod(1,0)==Aff{0,7}));}
    t.set_edge(1,0,{4,5},{6,7});assert((t.get_edge(id)==pair<Aff,Aff>{{6,7},{4,5}}));
    t.set_edge(id,{0,9},{8,1});assert((t.prod(0,1)==Aff{0,9}));
    t.set_edge(t.get_edge_id(1,2),{3,2});assert((t.prod(0,2)==Aff{0,29}));
    t.cut_edge(id);int next=t.link_edge(1,0,{5,6},{7,8});assert(next>id);
    assert((t.get_edge(0,1)==pair<Aff,Aff>{{7,8},{5,6}}));
    t.link(0,2);assert(t.get_edge_id(0,2)==-1); // old link no-op
}
void long_path(){
    const int n=200000;LinkCutTree<int,add,zero> t(n);
    for(int i=0;i<n;i++){t.set(i,1);if(i)t.link_edge(i-1,i,2,3);}
    for(int j=0;j<100;j++){
        assert(t.prod(0,n-1)==3*n-2);assert(t.prod(n-1,0)==4*n-3);
        assert(t.kth_vertex(0,n-1,n/2)==n/2);assert(t.kth_edge(n-1,0,n/2)==n-2-n/2);
    }
    for(int i=1;i<n;i+=2)t.cut(i-1,i);
    for(int i=1;i<n;i+=2)t.link(i,i-1);
    assert(t.vertex_count(0,n-1)==n);
}
void corner_cases(){
    LinkCutTree<int,add,zero> empty(0), one(1);
    one.set(0,7); assert(one.prod(0,0)==7 && *one.try_prod(0,0)==7);
    one.link(0,0); one.cut(0,0); assert(one.get_edge_id(0,0)==-1);
    assert(one.vertex_count(0,0)==1 && one.edge_count(0,0)==0);
    assert(one.kth_vertex(0,0,0)==0 && one.find_root(0)==0);
    LinkCutTree<Aff,compose,unit> t(10000);
    Aff forward=unit(),back=unit();
    for(int i=1;i<10000;i++){
        Aff f{1+i%6,i%101},g{1+(i+3)%6,i%103};
        t.link_edge(i-1,i,f,g); forward=compose(forward,f); back=compose(g,back);
    }
    for(int i=0;i<100;i++){
        assert(t.prod(0,9999)==forward); assert(t.prod(9999,0)==back);
    }
    LinkCutTree<int,add,zero> p(3);
    p.link(0,1);p.link(1,2);p.prod(1,0);
    int r=p.find_root(0), id=p.get_edge_id(0,1);
    p.set_edge(id,2,3);assert(p.find_root(0)==r);assert(p.get_edge_id(0,1)==id);
    p.cut(1,0);int next=p.link_edge(1,0,4,5);assert(next>id);
    assert(p.prod(0,2)==5 && p.prod(2,0)==4);
}
int main(int argc,char**argv){
    int seeds=argc>1?atoi(argv[1]):100;
    directed(); corner_cases();
    for(int s=0;s<seeds;s++){
        random_test<int,add,zero>(s,[](auto& r){return int(r()%21)-10;});
        random_test<Aff,compose,unit>(s,[](auto& r){return Aff{int(r()%5),int(r()%MOD)};});
    }
    long_path();cout<<"PASS "<<seeds*20000<<" randomized steps + direction/200000-vertex path\n";
}
