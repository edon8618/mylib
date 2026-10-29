#pragma once
/** Dynamic Li Chao Tree: maximum queries, closed domain [min_x,max_x].
 * add_segment(m,c,l,r) uses a half-open interval [l,r).
 * Integer comparisons/evaluations use __int128; query throws if result cannot fit T.
 * Floating-point version subdivides at representable coordinates, without an EPS cutoff.
 * All floating inputs must be finite; evaluation accuracy follows T's arithmetic.
 * Nodes are owned by a vector pool. Copy/move have ordinary value semantics.
 */
template<typename T>
struct DynamicLiChaoTree {
    static_assert(std::is_floating_point_v<T> || (std::is_integral_v<T> && std::is_signed_v<T> && sizeof(T) <= 8));
    using Eval = std::conditional_t<std::is_integral_v<T>, __int128_t, T>;
    struct Line { T m, c; Eval eval(T x) const { return Eval(m) * Eval(x) + Eval(c); } };
    struct Node { Line line{0,0}; int left=-1, right=-1; bool has_line=false; };
    std::vector<Node> nodes;
    int root=-1;
    T min_x, max_x, INF_VAL;
    DynamicLiChaoTree(T lo, T hi, T identity) : min_x(lo), max_x(hi), INF_VAL(identity) {
        if (lo > hi) throw std::invalid_argument("invalid domain");
        check_finite(lo); check_finite(hi); check_finite(identity);
    }
    void add_line(T m,T c) {
        check_finite(m); check_finite(c);
        root=add_line_impl(root,{m,c},min_x,max_x);
    }
    void add_segment(T m,T c,T l,T r) {
        check_finite(m);check_finite(c);check_finite(l);check_finite(r);
        if(l>=r)return;
        root=add_segment_impl(root,{m,c},l,r,min_x,max_x);
    }
    T query(T x) const {
        check_finite(x);
        if(x<min_x||x>max_x)throw std::out_of_range("query outside domain");
        bool found=false; Eval ans{};
        query_impl(root,x,min_x,max_x,found,ans);
        if(!found)return INF_VAL;
        if constexpr(std::is_integral_v<T>) {
            if(ans < Eval(std::numeric_limits<T>::min()) || ans > Eval(std::numeric_limits<T>::max()))
                throw std::overflow_error("Li Chao answer does not fit T");
        }
        return T(ans);
    }
private:
    static void check_finite(T x) {
        if constexpr(std::is_floating_point_v<T>)
            if(!std::isfinite(x))throw std::invalid_argument("non-finite coordinate/coefficient");
    }
    static T midpoint(T l,T r) {
        if constexpr(std::is_integral_v<T>) return T(Eval(l)+(Eval(r)-Eval(l))/2);
        else {
            T m=l/2+r/2;
            if(m>=r)m=std::nextafter(r,l);
            if(m<l)m=l;
            return m;
        }
    }
    static T successor(T x) {
        if constexpr(std::is_integral_v<T>) return x+1;
        else return std::nextafter(x,std::numeric_limits<T>::infinity());
    }
    int make_node() { nodes.emplace_back();return int(nodes.size())-1; }
    int add_line_impl(int p,Line line,T l,T r) {
        if(p<0)p=make_node();
        if(!nodes[p].has_line){nodes[p].line=line;nodes[p].has_line=true;return p;}
        T mid=midpoint(l,r);
        bool left=line.eval(l)>nodes[p].line.eval(l);
        bool center=line.eval(mid)>nodes[p].line.eval(mid);
        if(center)std::swap(line,nodes[p].line);
        if(l==r)return p;
        if(left!=center){
            int child=add_line_impl(nodes[p].left,line,l,mid); nodes[p].left=child;
        }else if(line.eval(r)>nodes[p].line.eval(r)){
            int child=add_line_impl(nodes[p].right,line,successor(mid),r); nodes[p].right=child;
        }
        return p;
    }
    int add_segment_impl(int p,Line line,T a,T b,T l,T r) {
        if(r<a||b<=l)return p;
        if(a<=l&&r<b)return add_line_impl(p,line,l,r);
        if(p<0)p=make_node();
        if(l==r)return p;
        T mid=midpoint(l,r);
        int left=add_segment_impl(nodes[p].left,line,a,b,l,mid);nodes[p].left=left;
        int right=add_segment_impl(nodes[p].right,line,a,b,successor(mid),r);nodes[p].right=right;
        return p;
    }
    void query_impl(int p,T x,T l,T r,bool& found,Eval& ans) const {
        if(p<0)return;
        if(nodes[p].has_line){Eval val=nodes[p].line.eval(x);if(!found||val>ans){found=true;ans=val;}}
        if(l==r)return;
        T mid=midpoint(l,r);
        if(x<=mid)query_impl(nodes[p].left,x,l,mid,found,ans);
        else query_impl(nodes[p].right,x,successor(mid),r,found,ans);
    }
};
