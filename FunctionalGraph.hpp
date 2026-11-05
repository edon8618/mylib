#pragma once
#include <bits/stdc++.h>

// Static functional graph, 0 <= to[v] < n, n >= 1.
// The graph consists of directed cycles with in-trees.
// A long-path ladder + a single binary lift gives O(1) online jump.
struct FunctionalGraph {
    int n=0, lg=0;
    std::vector<int> to, depth, entry, cycle_id, cycle_pos;
    std::vector<int> cycle_begin, cycle_vertices;
    std::vector<std::vector<int>> up;
    std::vector<int> ladder, ladder_index;

    FunctionalGraph()=default;
    explicit FunctionalGraph(std::vector<int> successor) { build(std::move(successor)); }

    void build(std::vector<int> successor) {
        to=std::move(successor);
        n=(int)to.size();
        assert(n>=1);
        std::vector<int> deg(n,0);
        for(int v=0;v<n;++v) { assert(0<=to[v] && to[v]<n); ++deg[to[v]]; }
        std::vector<int> removed;
        removed.reserve(n);
        for(int v=0;v<n;++v) if(!deg[v]) removed.push_back(v);
        std::vector<int> heavy(n,-1), height(n,0);
        for(size_t qi=0;qi<removed.size();++qi) {
            int v=removed[qi], p=to[v];
            if(height[p]<height[v]+1) {
                height[p]=height[v]+1;
                heavy[p]=v;
            }
            if(!--deg[p])removed.push_back(p);
        }
        depth.assign(n,0);
        entry.resize(n);
        cycle_id.assign(n,-1);
        cycle_pos.assign(n,-1);
        cycle_begin.clear();cycle_vertices.clear();
        cycle_begin.push_back(0);
        for(int i=0;i<n;++i) {
            if(!deg[i] || cycle_pos[i]!=-1)continue;
            int c=(int)cycle_begin.size()-1;
            for(int v=i;cycle_pos[v]==-1;v=to[v]) {
                cycle_id[v]=c;
                cycle_pos[v]=(int)cycle_vertices.size()-cycle_begin.back();
                entry[v]=v;
                cycle_vertices.push_back(v);
            }
            cycle_begin.push_back((int)cycle_vertices.size());
        }
        for(int j=(int)removed.size()-1;j>=0;--j) {
            int v=removed[j],p=to[v];
            depth[v]=depth[p]+1;
            entry[v]=entry[p];
            cycle_id[v]=cycle_id[p];
        }
        int max_depth=0;
        for(int d:depth)max_depth=std::max(max_depth,d);
        lg=1;
        while((1ULL<<lg)<=(unsigned long long)max_depth)++lg;
        up.assign(lg,std::vector<int>(n));
        for(int v=0;v<n;++v)up[0][v]=depth[v]?to[v]:v;
        for(int j=1;j<lg;++j) {
            auto &cur=up[j];
            const auto &prev=up[j-1];
            for(int v=0;v<n;++v)cur[v]=prev[prev[v]];
        }
        // Long-path decomposition, extend each path upward by its own length.
        ladder.clear();
        ladder.reserve(2*n);
        ladder_index.assign(n,-1);
        std::vector<int> path;
        for(int v=0;v<n;++v) {
            if(depth[v]>0 && heavy[to[v]]==v)continue;
            path.clear();
            for(int x=v;x!=-1;x=heavy[x])path.push_back(x);
            int len=(int)path.size();
            for(int j=len-1;j>=0;--j) {
                int x=path[j];
                ladder_index[x]=(int)ladder.size();
                ladder.push_back(x);
            }
            int x=v;
            for(int j=0;j<len && depth[x]>0;++j) {
                x=to[x];
                ladder.push_back(x);
            }
        }
    }

    int cycles() const { return (int)cycle_begin.size()-1; }
    bool on_cycle(int v) const { return depth[v]==0; }
    int distance_to_cycle(int v) const { return depth[v]; }
    int entry_vertex(int v) const { return entry[v]; }
    int cycle_length(int v) const {
        int c=cycle_id[v];
        return cycle_begin[c+1]-cycle_begin[c];
    }

    // Only for 0 <= steps <= depth[v]. O(1).
    int ancestor(int v,int steps) const {
        assert(0<=v && v<n && 0<=steps && steps<=depth[v]);
        if(steps==0)return v;
        int j=31-__builtin_clz((unsigned)steps);
        int mid=up[j][v];
        int rem=steps-(1<<j);
        assert(ladder_index[mid]+rem<(int)ladder.size());
        return ladder[ladder_index[mid]+rem];
    }

    // f^t(v), including t == 0 and very large t. O(1).
    int jump(int v,unsigned long long t) const {
        assert(0<=v && v<n);
        if(t<(unsigned)depth[v])return ancestor(v,(int)t);
        int c=cycle_id[v];
        int len=cycle_begin[c+1]-cycle_begin[c];
        int p=cycle_pos[entry[v]];
        int r=(int)((t-(unsigned)depth[v])%(unsigned)len);
        int i=p+r;
        if(i>=len)i-=len;
        return cycle_vertices[cycle_begin[c]+i];
    }
};

// Vertex weight w[x] is consumed when leaving vertex x.
// sum(u,t) = w[u] + w[f(u)] + ... + w[f^(t-1)(u)] (t terms).
// T needs +, -, and multiplication by uint64_t, and default T{}=0.
// Accumulation must not overflow T.
template<class T>
struct FunctionalGraphSum {
    const FunctionalGraph* fg=nullptr;
    std::vector<T> tree_sum, cycle_prefix, cycle_total;
    FunctionalGraphSum()=default;
    FunctionalGraphSum(const FunctionalGraph& graph, const std::vector<T>& w) { build(graph,w); }
    void build(const FunctionalGraph& graph,const std::vector<T>& w) {
        fg=&graph;
        assert((int)w.size()==fg->n);
        int n=fg->n;
        tree_sum.assign(n,T{});
        // Any reverse topological order: sort by depth (increasing).
        // Rebuild O(N) tree order using depth buckets, without sorting.
        std::vector<int> count(n+2,0);
        for(int v=0;v<n;++v)++count[fg->depth[v]+1];
        for(int i=1;i<=n;++i)count[i]+=count[i-1];
        std::vector<int> order(n);
        for(int v=0;v<n;++v)order[count[fg->depth[v]]++]=v;
        for(int v:order)if(fg->depth[v])tree_sum[v]=w[v]+tree_sum[fg->to[v]];
        int C=(int)fg->cycle_vertices.size();
        cycle_prefix.assign(C+1,T{});
        for(int i=0;i<C;++i)cycle_prefix[i+1]=cycle_prefix[i]+w[fg->cycle_vertices[i]];
        cycle_total.resize(fg->cycles());
        for(int c=0;c<fg->cycles();++c) {
            int l=fg->cycle_begin[c],r=fg->cycle_begin[c+1];
            cycle_total[c]=cycle_prefix[r]-cycle_prefix[l];
        }
    }
    T sum(int v,unsigned long long t) const {
        assert(fg && 0<=v && v<fg->n);
        int d=fg->depth[v];
        if(t<=(unsigned)d) {
            int to=fg->ancestor(v,(int)t);
            return tree_sum[v]-tree_sum[to];
        }
        T ans=tree_sum[v];
        unsigned long long remain=t-(unsigned)d;
        int c=fg->cycle_id[v];
        int l=fg->cycle_begin[c], len=fg->cycle_begin[c+1]-l;
        ans=ans+cycle_total[c]*(remain/(unsigned)len);
        int r=(int)(remain%(unsigned)len);
        int p=fg->cycle_pos[fg->entry[v]];
        if(p+r<=len)ans=ans+(cycle_prefix[l+p+r]-cycle_prefix[l+p]);
        else ans=ans+(cycle_prefix[l+len]-cycle_prefix[l+p])+(cycle_prefix[l+p+r-len]-cycle_prefix[l]);
        return ans;
    }
    T fold(int v,unsigned long long t) const { return sum(v,t); }
    int jump(int v,unsigned long long t) const { return fg->jump(v,t); }
    std::pair<int,T> jump_fold(int v,unsigned long long t) const {
        return {fg->jump(v,t),sum(v,t)};
    }
};

// Static RMQ with 64-element micro blocks and a sparse table on block maxima.
// O(N) (word-RAM) space, O(1) query, supports max and min.
template<class T,bool IsMax=true>
struct FunctionalGraphRMQ {
    static constexpr int B=64;
    T identity{};
    std::vector<T> a, prefix, suffix;
    std::vector<uint64_t> masks;
    std::vector<std::vector<T>> sparse;
    FunctionalGraphRMQ()=default;
    FunctionalGraphRMQ(const std::vector<T>& data,T id) { build(data,id); }
    bool better_or_equal(const T& x,const T& y) const {
        if constexpr(IsMax)return x>=y;
        else return x<=y;
    }
    T merge(const T& x,const T& y) const {
        if constexpr(IsMax)return std::max(x,y);
        else return std::min(x,y);
    }
    void build(const std::vector<T>& data,T id) {
        identity=id;a=data;
        int n=(int)a.size();
        prefix.resize(n); suffix.resize(n);masks.resize(n);
        int nb=(n+B-1)/B;
        std::vector<T> block(nb);
        for(int i=0;i<n;++i) {
            int off=i%B;
            prefix[i]=(off==0)?a[i]:merge(prefix[i-1],a[i]);
            uint64_t m=(off==0)?0:masks[i-1];
            while(m) {
                int highest=63-__builtin_clzll(m);
                if(!better_or_equal(a[i],a[i-off+highest]))break;
                m^=1ULL<<highest;
            }
            masks[i]=m|(1ULL<<off);
        }
        for(int i=n-1;i>=0;--i) {
            suffix[i]=(i==n-1 || i/B!=(i+1)/B)?a[i]:merge(a[i],suffix[i+1]);
        }
        for(int b=0;b<nb;++b)block[b]=prefix[std::min((b+1)*B,n)-1];
        sparse.clear();
        if(nb) {
            sparse.push_back(block);
            for(int j=1;(1<<j)<=nb;++j) {
                int cnt=nb-(1<<j)+1;
                sparse.emplace_back(cnt);
                for(int i=0;i<cnt;++i)sparse[j][i]=merge(sparse[j-1][i],sparse[j-1][i+(1<<(j-1))]);
            }
        }
    }
    // Query [l,r). O(1).
    T query(int l,int r) const {
        assert(0<=l && l<=r && r<=(int)a.size());
        if(l==r)return identity;
        int lb=l/B,rb=(r-1)/B;
        if(lb==rb) {
            uint64_t m=masks[r-1] & (~0ULL << (l%B));
            int idx=__builtin_ctzll(m);
            return a[lb*B+idx];
        }
        T ans=merge(suffix[l],prefix[r-1]);
        int first=lb+1,last=rb;
        if(first<last) {
            int cnt=last-first;
            int j=31-__builtin_clz((unsigned)cnt);
            ans=merge(ans,merge(sparse[j][first],sparse[j][last-(1<<j)]));
        }
        return ans;
    }
};

// O(1) max/min over weighted walks.
// IsMax=true: max. IsMax=false: min.
// w[v] is included when leaving v; t=0 returns the supplied identity.
template<class T,bool IsMax=true>
struct FunctionalGraphExtreme {
    const FunctionalGraph* fg=nullptr;
    T identity;
    std::vector<std::vector<T>> block; // first 2^j vertices along the in-tree
    FunctionalGraphRMQ<T,IsMax> ladder_rmq, cycle_rmq;
    std::vector<T> cycle_extreme;
    FunctionalGraphExtreme(const FunctionalGraph& graph,const std::vector<T>& w,
        T id=IsMax?std::numeric_limits<T>::lowest():std::numeric_limits<T>::max()) { build(graph,w,id); }
    T merge(const T& a,const T& b) const {
        if constexpr(IsMax)return std::max(a,b);
        else return std::min(a,b);
    }
    void build(const FunctionalGraph& graph,const std::vector<T>& w,T id) {
        fg=&graph;identity=id;
        assert((int)w.size()==fg->n);
        block.assign(fg->lg,std::vector<T>(fg->n));
        block[0]=w;
        for(int j=1;j<fg->lg;++j)for(int v=0;v<fg->n;++v) {
            block[j][v]=merge(block[j-1][v],block[j-1][fg->up[j-1][v]]);
        }
        std::vector<T> lad;
        lad.reserve(fg->ladder.size());
        for(int v:fg->ladder)lad.push_back(w[v]);
        ladder_rmq.build(lad,id);
        std::vector<T> cyc;
        cyc.reserve(2*fg->cycle_vertices.size());
        cycle_extreme.assign(fg->cycles(),id);
        for(int c=0;c<fg->cycles();++c) {
            int l=fg->cycle_begin[c],r=fg->cycle_begin[c+1];
            for(int copies=0;copies<2;++copies)for(int i=l;i<r;++i)cyc.push_back(w[fg->cycle_vertices[i]]);
            for(int i=l;i<r;++i)cycle_extreme[c]=merge(cycle_extreme[c],w[fg->cycle_vertices[i]]);
        }
        cycle_rmq.build(cyc,id);
    }
    T tree_fold(int v,int t) const {
        assert(0<=t && t<=fg->depth[v]);
        if(t==0)return identity;
        int j=31-__builtin_clz((unsigned)t);
        int mid=fg->up[j][v];
        T ans=block[j][v];
        int rem=t-(1<<j);
        if(rem)ans=merge(ans,ladder_rmq.query(fg->ladder_index[mid],fg->ladder_index[mid]+rem));
        return ans;
    }
    T fold(int v,unsigned long long t) const {
        assert(fg && 0<=v && v<fg->n);
        if(t<=(unsigned)fg->depth[v])return tree_fold(v,(int)t);
        T ans=tree_fold(v,fg->depth[v]);
        unsigned long long rem=t-fg->depth[v];
        int c=fg->cycle_id[v],len=fg->cycle_length(v);
        if(rem>=(unsigned)len)ans=merge(ans,cycle_extreme[c]);
        int r=(int)(rem%(unsigned)len);
        if(r) {
            int off=2*fg->cycle_begin[c];
            int p=fg->cycle_pos[fg->entry[v]];
            ans=merge(ans,cycle_rmq.query(off+p,off+p+r));
        }
        return ans;
    }
    T extreme(int v,unsigned long long t) const {return fold(v,t);}
    T max(int v,unsigned long long t) const { static_assert(IsMax, "Use min() on FunctionalGraphMin"); return fold(v,t); }
    T min(int v,unsigned long long t) const { static_assert(!IsMax, "Use max() on FunctionalGraphMax"); return fold(v,t); }
    int jump(int v,unsigned long long t) const { return fg->jump(v,t); }
    std::pair<int,T> jump_fold(int v,unsigned long long t) const {
        return {fg->jump(v,t),fold(v,t)};
    }
};

template<class T>using FunctionalGraphMax=FunctionalGraphExtreme<T,true>;
template<class T>using FunctionalGraphMin=FunctionalGraphExtreme<T,false>;

// Generic associative fold, including noncommutative monoids.
// M must provide:
//   using S = ...;
//   static S op(const S&, const S&);
//   static S e();
// Optionally: static S repeat(const S&, uint64_t count); // fast cycle powers.
// Query: O(log N + log t), construction: O(N log N).
template<class M>
struct FunctionalGraphMonoid {
    using S=typename M::S;
    const FunctionalGraph* fg=nullptr;
    std::vector<std::vector<S>> block;
    int seg_size=1;
    std::vector<S> cycle_seg;
    FunctionalGraphMonoid()=default;
    FunctionalGraphMonoid(const FunctionalGraph& graph,const std::vector<S>& w) { build(graph,w); }
    void build(const FunctionalGraph& graph,const std::vector<S>& w) {
        fg=&graph;
        assert((int)w.size()==fg->n);
        block.assign(fg->lg,std::vector<S>(fg->n));
        block[0]=w;
        for(int j=1;j<fg->lg;++j)for(int v=0;v<fg->n;++v)
            block[j][v]=M::op(block[j-1][v],block[j-1][fg->up[j-1][v]]);
        int csz=(int)fg->cycle_vertices.size();
        seg_size=1;
        while(seg_size<2*csz)seg_size*=2;
        cycle_seg.assign(2*seg_size,M::e());
        for(int c=0;c<fg->cycles();++c) {
            int l=fg->cycle_begin[c],r=fg->cycle_begin[c+1],sz=r-l;
            for(int i=0;i<sz;++i) {
                S x=w[fg->cycle_vertices[l+i]];
                cycle_seg[seg_size+2*l+i]=x;
                cycle_seg[seg_size+2*l+sz+i]=x;
            }
        }
        for(int i=seg_size-1;i>0;--i)cycle_seg[i]=M::op(cycle_seg[i*2],cycle_seg[i*2+1]);
    }
    S cycle_range(int l,int r) const { // order-preserving fold of [l,r)
        S front=M::e(),back=M::e();
        for(l+=seg_size,r+=seg_size;l<r;l>>=1,r>>=1) {
            if(l&1)front=M::op(front,cycle_seg[l++]);
            if(r&1)back=M::op(cycle_seg[--r],back);
        }
        return M::op(front,back);
    }
    static S power(S x,unsigned long long n) {
        if constexpr(requires(const S& z, unsigned long long k){M::repeat(z,k);})
            return M::repeat(x,n);
        else {
            S res=M::e();
            while(n) {
                if(n&1)res=M::op(res,x);
                n>>=1;
                if(n)x=M::op(x,x);
            }
            return res;
        }
    }
    S tree_fold(int v,int k) const {
        assert(0<=k && k<=fg->depth[v]);
        S ans=M::e();
        for(int j=0;k;++j,k>>=1)if(k&1) {
            ans=M::op(ans,block[j][v]);
            v=fg->up[j][v];
        }
        return ans;
    }
    // Fold w[u], w[f(u)], ... w[f^(t-1)(u)] (t terms).
    S fold(int v,unsigned long long t) const {
        assert(fg && 0<=v && v<fg->n);
        int d=fg->depth[v];
        if(t<=(unsigned)d)return tree_fold(v,(int)t);
        S ans=tree_fold(v,d);
        unsigned long long rem=t-(unsigned)d;
        int c=fg->cycle_id[v], len=fg->cycle_length(v);
        int off=2*fg->cycle_begin[c]+fg->cycle_pos[fg->entry[v]];
        if(rem>=(unsigned)len) {
            S whole=cycle_range(off,off+len);
            ans=M::op(ans,power(whole,rem/(unsigned)len));
        }
        int r=(int)(rem%(unsigned)len);
        if(r)ans=M::op(ans,cycle_range(off,off+r));
        return ans;
    }
    int jump(int v,unsigned long long t) const { return fg->jump(v,t); }
    std::pair<int,S> jump_fold(int v,unsigned long long t) const {
        return {fg->jump(v,t),fold(v,t)};
    }
};
