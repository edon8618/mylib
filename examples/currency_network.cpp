// New-API example; the supplied problem summary did not specify input syntax.
// This file uses: N Q, then Q queries (0-based vertices, initially empty forest):
// 0 u v a b c d : add u->v: a*t+b, v->u: c*t+d (different components)
// 1 u v         : delete the edge
// 2 u v a b c d : update in the ORIGINAL opening direction, regardless of u,v order
// 3 x y s       : convert s from x to y, or print -1 if disconnected
// Adjust only this parser to the actual judge's format.
#include <bits/stdc++.h>
#include <atcoder/modint>
#include "../LinkCutTree.hpp"
using namespace std;
using mint = atcoder::modint998244353;
struct S { mint a,b; };
S op(S l,S r){return {r.a*l.a,r.a*l.b+r.b};}
S e(){return {1,0};}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int n,q;if(!(cin>>n>>q))return 0;
    LinkCutTree<S,op,e> lct(n);
    while(q--){
        int type,u,v;cin>>type>>u>>v;
        if(type==0||type==2){
            long long a,b,c,d;cin>>a>>b>>c>>d;
            if(type==0)lct.link_edge(u,v,{a,b},{c,d});
            else lct.set_edge(lct.get_edge_id(u,v),{a,b},{c,d});
        }else if(type==1)lct.cut(u,v);
        else{
            long long s;cin>>s;auto f=lct.try_prod(u,v);
            if(!f)cout<<-1<<'\n';else cout<<(f->a*s+f->b).val()<<'\n';
        }
    }
}
