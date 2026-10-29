#pragma once
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3,unroll-loops")
#endif
// Dynamic Graph Vertex Add Component Sum -- online HDT + Euler Tour Treap
// Library Checker: https://judge.yosupo.jp/problem/dynamic_graph_vertex_add_component_sum
// Component aggregates via a COMMUTATIVE monoid.
// Requirements: S copyable/default-constructible, op associative + commutative, e identity.
// API: add/erase/same/size/find, set(v,x), get(v), apply(v,x), prod(v).
#include <bits/stdc++.h>
#ifndef DC_SCAN_BUDGET
#define DC_SCAN_BUDGET 32
#endif
// Online dynamic connectivity for a simple undirected graph, 0 <= v < n.
// Fully online add/erase/same/size/find. Duplicate add and absent erase return false.
// HDT-style level hierarchy + Euler-tour trees (randomized implicit treaps).
// find returns an opaque component token: compare only within the current state.
// Optional explicit seed is intended for REPRODUCIBLE TESTING ONLY.
// Use the default secret seed against oblivious adversarial inputs.
// NOTE: 600 ms worst-case target is NOT met on all tested cases.
// Optional: DC_SKIP_POOL_FREE=1 intentionally leaks node/edge pools at destruction.
// Use only for a one-instance judge process, to let the OS reclaim memory at exit.
#ifndef DC_SKIP_POOL_FREE
#define DC_SKIP_POOL_FREE 0
#endif
template<class S, S (*op)(S,S), S (*identity)()>
struct DynamicConnectivityHDT {
    struct Node {
        Node *l=nullptr,*r=nullptr,*p=nullptr;
        uint32_t pri;
        int vertex, sz=1, cnt, adj=0, sum=0, mark=-1, marks=0;
        uint64_t xval=0,xsum=0;
        S weight=identity(), weight_sum=identity();
        Node()=default; Node(uint32_t z,int v):pri(z),vertex(v),cnt(v>=0){}
    };
    struct Pool {
        static constexpr int B=1<<15;
        std::vector<std::unique_ptr<Node[]>> blocks;
        std::vector<Node*> reclaimed;
        int at=B;
        Node* make(uint32_t p,int v){
            Node* x;
            if(!reclaimed.empty()){x=reclaimed.back();reclaimed.pop_back();}
            else {if(at==B){blocks.emplace_back(std::make_unique<Node[]>(B));at=0;}x=&blocks.back()[at++];}
            *x=Node(p,v);return x;
        }
        void recycle(Node*x){reclaimed.push_back(x);}
    } pool;
    struct Forest {
        DynamicConnectivityHDT *dc;
        std::vector<Node*> vertex_nodes;
        explicit Forest(DynamicConnectivityHDT*d):dc(d),vertex_nodes(d->n){for(int i=0;i<d->n;i++)vertex_nodes[i]=dc->pool.make(d->random32(),i);}
        ~Forest()=default;
        static int sz(Node*x){return x?x->sz:0;}
        static int cnt(Node*x){return x?x->cnt:0;}
        static int sum(Node*x){return x?x->sum:0;}
        static int marks(Node*x){return x?x->marks:0;}
        static uint64_t xsum(Node*x){return x?x->xsum:0;}
        static S weight_sum(Node*x){return x?x->weight_sum:identity();}
        static void pull(Node*x){if(!x)return;x->sz=1+sz(x->l)+sz(x->r);x->cnt=(x->vertex>=0)+cnt(x->l)+cnt(x->r);x->sum=x->adj+sum(x->l)+sum(x->r);x->marks=(x->mark!=-1)+marks(x->l)+marks(x->r);x->xsum=x->xval^xsum(x->l)^xsum(x->r);x->weight_sum=op(op(weight_sum(x->l),x->weight),weight_sum(x->r));if(x->l)x->l->p=x;if(x->r)x->r->p=x;}
        static Node*merge(Node*a,Node*b){if(!a){if(b)b->p=nullptr;return b;}if(!b){a->p=nullptr;return a;}if(a->pri>b->pri){a->r=merge(a->r,b);pull(a);a->p=nullptr;return a;}b->l=merge(a,b->l);pull(b);b->p=nullptr;return b;}
        static std::pair<Node*,Node*> split(Node*t,int k){if(!t)return {nullptr,nullptr};if(sz(t->l)>=k){auto[a,b]=split(t->l,k);t->l=b;pull(t);t->p=nullptr;if(a)a->p=nullptr;return {a,t};}auto[a,b]=split(t->r,k-sz(t->l)-1);t->r=a;pull(t);t->p=nullptr;if(b)b->p=nullptr;return {t,b};}
        static Node*root(Node*x){while(x->p)x=x->p;return x;}
        static int index(Node*x){int i=sz(x->l);while(x->p){if(x==x->p->r)i+=sz(x->p->l)+1;x=x->p;}return i;}
        Node*reroot(int v){Node*t=root(vertex_nodes[v]);auto[a,b]=split(t,index(vertex_nodes[v]));return merge(b,a);}
        bool same(int u,int v)const{return root(vertex_nodes[u])==root(vertex_nodes[v]);}
        int size(int u)const{return cnt(root(vertex_nodes[u]));}
        uintptr_t find(int u)const{return uintptr_t(root(vertex_nodes[u]));}
        void touch(Node*x){while(x){pull(x);x=x->p;}}
        void bump(int v,int d){vertex_nodes[v]->adj+=d;touch(vertex_nodes[v]);}
        void toggle_xor(int v,uint64_t k){vertex_nodes[v]->xval^=k;touch(vertex_nodes[v]);}
        uint64_t component_xor(int v){return root(vertex_nodes[v])->xsum;}
        void link(int u,int v,Node*&a,Node*&b,int mark=-1){assert(!same(u,v));Node*x=reroot(u),*y=reroot(v);a=dc->pool.make(dc->random32(),-1);b=dc->pool.make(dc->random32(),-1);a->mark=mark;pull(a);merge(merge(merge(x,a),y),b);}
        void cut(Node*a,Node*b){assert(a&&b);Node*t=root(a);int i=index(a),j=index(b);if(i>j)std::swap(i,j);auto[x,p]=split(t,i);auto[del1,q]=split(p,1);auto[y,r]=split(q,j-i-1);auto[del2,z]=split(r,1);merge(z,x);dc->pool.recycle(del1);dc->pool.recycle(del2);}
        Node* next_mark_root(Node*t){if(t->marks==0)return nullptr;while(t){if(marks(t->l)){t=t->l;continue;}if(t->mark!=-1)return t;t=t->r;}return nullptr;}
        int kth_adj(Node*t,int k){
            while(t){
                int l=sum(t->l);
                if(k<l){t=t->l;continue;}
                k-=l;
                if(t->adj){if(k==0)return t->vertex;k--;}
                t=t->r;
            }
            return -1;
        }
        int next_adj_root(Node*t){if(t->sum==0)return -1;while(t){if(sum(t->l)){t=t->l;continue;}if(t->adj)return t->vertex;t=t->r;}return -1;}
    };
    struct Edge {
        int u,v,level=0,id=-1;bool tree=false;
        Edge *prev[2]{nullptr,nullptr}, *next[2]{nullptr,nullptr};
        struct Arcs {
            std::pair<Node*,Node*> in[4]{};
            std::vector<std::pair<Node*,Node*>> tail;
            int cnt=0;
            int size()const{return cnt;}
            auto& operator[](int i){return i<4?in[i]:tail[i-4];}
            const auto& operator[](int i)const{return i<4?in[i]:tail[i-4];}
            void push_back(std::pair<Node*,Node*> a){
                if(cnt<4)in[cnt]=a;
                else tail.push_back(a);
                ++cnt;
            }
            void clear(){cnt=0;tail.clear();}
        } arcs;
        uint64_t label=0;
    };
    static uint64_t entropy(){
        std::random_device rd;
        uint64_t x=(uint64_t(rd())<<32)^rd();
        x ^= uint64_t(std::chrono::steady_clock::now().time_since_epoch().count());
        return x?x:0x9e3779b97f4a7c15ull;
    }
    struct EdgeMap {
        struct Slot {uint64_t key=0;Edge* val=nullptr;};
        std::vector<Slot> a;
        size_t mask=0,live=0,used=0;
        uint64_t salt;
        explicit EdgeMap(uint64_t seed):salt(seed){}
        size_t bucket(uint64_t k)const{
            uint64_t x=k^salt;
            x^=x>>32;
            x*=0xd6e8feb86659fd93ull;
            x^=x>>32;
            return size_t(x)&mask;
        }
        void rehash(size_t cap){
            std::vector<Slot> b(cap);
            mask=cap-1;
            for(auto &z:a)if(z.key && z.key!=UINT64_MAX){
                size_t i=bucket(z.key);
                while(b[i].key)i=(i+1)&mask;
                b[i]=z;
            }
            a.swap(b);used=live;
        }
        void reserve(size_t n){size_t cap=16;while(cap<n*2)cap*=2;rehash(cap);}
        Edge* find(uint64_t k)const{
            if(a.empty())return nullptr;
            size_t i=bucket(k);
            while(a[i].key){if(a[i].key==k)return a[i].val;i=(i+1)&mask;}
            return nullptr;
        }
        void insert(uint64_t k,Edge* e){
            if((used+1)*10 >= a.size()*7)rehash(a.size()?a.size()*2:16);
            size_t i=bucket(k),hole=SIZE_MAX;
            while(a[i].key){if(a[i].key==UINT64_MAX && hole==SIZE_MAX)hole=i;i=(i+1)&mask;}
            if(hole!=SIZE_MAX)i=hole;else used++;
            a[i].key=k;a[i].val=e;live++;
        }
        void erase(uint64_t k){
            size_t i=bucket(k);
            while(a[i].key!=k)i=(i+1)&mask;
            a[i].key=UINT64_MAX;a[i].val=nullptr;live--;
        }
    } es;
    struct EdgePool{
        static constexpr int B=1<<13;
        std::vector<std::unique_ptr<Edge[]>> blocks;
        int at=B;
        Edge* make(int u,int v,int id){
            if(at==B){blocks.emplace_back(std::make_unique<Edge[]>(B));at=0;}
            Edge* e=&blocks.back()[at++];
            e->u=u;e->v=v;e->id=id;e->level=0;e->tree=false;
            e->prev[0]=e->prev[1]=e->next[0]=e->next[1]=nullptr;
            e->arcs.clear();e->label=0;
            return e;
        }
    } epool;
    int n,L; size_t extra_count=0;
    uint64_t pending_key=0;
    struct Undo {uint64_t key=0;Edge* rep=nullptr;int lev=-1;};
    Undo recent[32]{}; unsigned recent_at=0;
    void remember(Edge* r,int lev){recent[recent_at++%32]={pending_key,r,lev};}

    uint64_t rng;
    std::vector<std::unique_ptr<Forest>> f;
    std::vector<std::vector<Edge*>> heads;

    uint32_t random32(){rng^=rng<<7;rng^=rng>>9;return uint32_t(rng);}
    static uint64_t key(int a,int b){if(a>b)std::swap(a,b);return (uint64_t(uint32_t(a))<<32)|uint32_t(b);}
    explicit DynamicConnectivityHDT(int n_, uint64_t seed=entropy()):es(seed^0xc6a4a7935bd1e995ull),n(n_),L(1),rng(seed?seed:0x9e3779b97f4a7c15ull){while((1LL<<L)<=std::max(1,n))++L; f.emplace_back(std::make_unique<Forest>(this));heads.emplace_back(n,nullptr);es.reserve(n+16);ids.reserve(n+16);}
    ~DynamicConnectivityHDT(){
#if DC_SKIP_POOL_FREE
        for(auto &x:pool.blocks)x.release();
        for(auto &x:epool.blocks)x.release();
#endif
    }
    void ensure(int lev){while(int(f.size())<=lev){f.emplace_back(std::make_unique<Forest>(this));heads.emplace_back(n,nullptr);}}
    bool same(int u,int v)const{return f[0]->same(u,v);}
    int size(int u)const{return f[0]->size(u);}
    // Commutative monoid over vertices. The aggregate order in ETT is not stable.
    S prod(int v)const{return Forest::weight_sum(Forest::root(f[0]->vertex_nodes[v]));}
    S get(int v)const{return f[0]->vertex_nodes[v]->weight;}
    void set(int v,S x){Node* p=f[0]->vertex_nodes[v];p->weight=std::move(x);f[0]->touch(p);}
    // Apply a monoid element to the value of vertex v (left-to-right op).
    void apply(int v,S x){Node* p=f[0]->vertex_nodes[v];p->weight=op(p->weight,x);f[0]->touch(p);}
    uintptr_t find(int u)const{return f[0]->find(u);}
    static int side(const Edge *e,int v){return e->u==v?0:1;}
    void add_extra(Edge&e){
        ++extra_count;int i=e.level;ensure(i);
        for(int t=0;t<2;t++){
            int x=t?e.v:e.u;
            Edge*&h=heads[i][x];
            if(!h)f[i]->bump(x,1);
            e.prev[t]=nullptr;e.next[t]=h;
            if(h)h->prev[side(h,x)]=&e;
            h=&e;
        }
    }
    void del_extra(Edge&e){
        --extra_count;int i=e.level;
        for(int t=0;t<2;t++){
            int x=t?e.v:e.u;
            Edge*&h=heads[i][x];
            Edge *p=e.prev[t],*q=e.next[t];
            if(p)p->next[side(p,x)]=q;else h=q;
            if(q)q->prev[side(q,x)]=p;
            if(!h)f[i]->bump(x,-1);
            e.prev[t]=e.next[t]=nullptr;
        }
    }
    // Markers store an edge index into the stable vector 'ids'.
    std::vector<uint64_t> ids;

    void unset_marker(Edge&e,int lev){auto*x=e.arcs[lev].first;x->mark=-1;f[lev]->touch(x);}
    void tree_link(Edge&e,int top,int id){ensure(top);e.tree=true;for(int j=int(e.arcs.size());j<=top;j++){Node*a,*b;f[j]->link(e.u,e.v,a,b,j==top?id:-1);e.arcs.push_back({a,b});}}
    bool add(int u,int v){if(u==v)return false;uint64_t k=key(u,v);if(es.find(k))return false;
        Undo candidate;
        for(unsigned j=0;j<32;j++){
            auto &p=recent[(recent_at-1-j)%32];
            if(p.key==k){candidate=p;p.key=0;break;}
        }
        int id=ids.size();ids.push_back(k);Edge* ptr=epool.make(u,v,id);es.insert(k,ptr);Edge&e=*ptr;e.label=(uint64_t(random32())<<32)|uint32_t(id+1);
        Edge *previous=candidate.rep;int previous_level=candidate.lev;
        bool restore=false;
        if(previous && previous->tree && previous->level==previous_level &&
            es.find(key(previous->u,previous->v))==previous &&
            f[previous_level]->same(u,v)){
            for(int j=0;j<=previous_level;j++)f[j]->cut(previous->arcs[j].first,previous->arcs[j].second);
            previous->arcs.clear();
            if(!f[previous_level]->same(u,v)) restore=true;
            else tree_link(*previous,previous_level,previous->id);
        }
        if(restore){
            previous->tree=false;previous->level=previous_level;
            e.level=previous_level;
            tree_link(e,previous_level,id);
            add_extra(*previous);
            f[0]->toggle_xor(previous->u,previous->label);
            f[0]->toggle_xor(previous->v,previous->label);
        }else if(same(u,v)){add_extra(e);f[0]->toggle_xor(u,e.label);f[0]->toggle_xor(v,e.label);}else tree_link(e,0,id);
        return true;
    }
    void replace_with(Edge& e,int lev){
        del_extra(e);f[0]->toggle_xor(e.u,e.label);f[0]->toggle_xor(e.v,e.label);
        tree_link(e,lev,e.id);
        remember(&e,lev);
    }
    bool erase(int u,int v){if(u==v)return false;uint64_t k=key(u,v);Edge* ptr=es.find(k);if(!ptr)return false;pending_key=k;Edge&e=*ptr;int lev=e.level;if(!e.tree){del_extra(e);f[0]->toggle_xor(u,e.label);f[0]->toggle_xor(v,e.label);es.erase(k);return true;}
        for(int j=0;j<=lev;j++)f[j]->cut(e.arcs[j].first,e.arcs[j].second);
        es.erase(k);
        if(extra_count==0)return true;
        // Cutset XOR. Internal edges cancel, so a unique crossing edge
        // can be recovered; use it only when its level is reached.
        // Replacing a level-i tree edge with a lower-level edge *before*
        // processing higher levels would break the HDT invariants.
        uint64_t cutset=f[0]->component_xor(u);
        uint32_t low=uint32_t(cutset);
        Edge* ce=(low && low<=ids.size()) ? es.find(ids[low-1]) : nullptr;
        if(!ce || ce->tree || f[0]->same(ce->u,ce->v))ce=nullptr;
        for(int i=lev;i>=0;i--){
            if(ce && ce->level==i){
                replace_with(*ce,i);return true;
            }
            Node* ra=Forest::root(f[i]->vertex_nodes[u]);
            Node* rb=Forest::root(f[i]->vertex_nodes[v]);
            Node* small=ra;
            if(ra->cnt>rb->cnt)small=rb;
            if(small->sum==0)continue;
            // Cheap probe of non-tree edges around the removed edge's endpoints.
            // A successful probe needs no HDT promotion; at most 12 probes per level.
            for(int x : {u,v}){
                Edge* t=heads[i][x];
                for(int j=0;t && j<6;j++){
                    Edge* next=t->next[side(t,x)];
                    if(!f[i]->same(t->u,t->v)){
                        replace_with(*t,i);
                        return true;
                    }
                    t=next;
                }
            }

            // Attempt the first available edge before costly tree-edge promotion.
            int first=f[i]->next_adj_root(small);
            if(first>=0){
                Edge* t=heads[i][first];
                if(!f[i]->same(t->u,t->v)){
                    replace_with(*t,i);return true;
                }
            }
            // Bounded search for a crossing edge before promotion.
            // A nonzero randomized XOR fingerprint indicates that the
            // deleted forest edge has a replacement in some level.
            if(cutset && DC_SCAN_BUDGET>0){
                int cnt=std::min(small->sum,DC_SCAN_BUDGET);
                uint32_t offset=random32();
                for(int rep=0;rep<cnt;rep++){
                    int k=(uint64_t(rep)*small->sum/cnt+offset)%small->sum;
                    int x=f[i]->kth_adj(small,k);
                    Edge* t=heads[i][x];
                    for(int j=0;t&&j<2;j++){
                        Edge* nxt=t->next[side(t,x)];
                        int other=t->u==x?t->v:t->u;
                        if(Forest::root(f[i]->vertex_nodes[other])!=small){
                            replace_with(*t,i);return true;
                        }
                        t=nxt;
                    }
                }
            }
            // Promote exact-level tree edges in smaller component.
            if(i+1<L){ensure(i+1);while(Node*x=f[i]->next_mark_root(small)){
                int eid=x->mark;uint64_t ek=ids[eid];Edge* ptr2=es.find(ek);assert(ptr2);Edge&ee=*ptr2;assert(ee.tree&&ee.level==i);unset_marker(ee,i);ee.level++;tree_link(ee,i+1,eid);
            }}
            while(true){
                int x=f[i]->next_adj_root(small);if(x==-1)break;
                Edge&ee=*heads[i][x];uint64_t ek=key(ee.u,ee.v);del_extra(ee);
                int other=ee.u==x?ee.v:ee.u;
                if(Forest::root(f[i]->vertex_nodes[other])==small){
                    assert(i+1<L);ee.level++;add_extra(ee);
                }else{
                    f[0]->toggle_xor(ee.u,ee.label);f[0]->toggle_xor(ee.v,ee.label);
                    ee.level=i;tree_link(ee,i,ee.id);remember(&ee,i);return true;
                }
            }
        }
        return true;
    }
};
template<class S, S (*op)(S,S), S (*e)()>
using DynamicConnectivity = DynamicConnectivityHDT<S,op,e>;


