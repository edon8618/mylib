#pragma once
// Linear Matroid Parity (unweighted) -- C++17, header-only.
// Field: Mint is a type over an ODD PRIME field with Mint::mod(), Mint(int),
//        +, -, *, /, ==, != (e.g. atcoder::modint998244353).
//
// INPUT: d-dimensional vector pairs (a_i,b_i), i = 0..m-1.
//        A chosen set of pairs is feasible iff its 2*|S| vectors are independent.
// OUTPUT: {maximum pair count, selected pair IDs if restore == true}.
//
// Usage:
//   auto ans = linear_matroid_parity<Mint>(d, m, get_pair, restore=false);
//   auto ans = linear_matroid_parity_sparse<Mint>(d, m, get_sparse_pair, restore=false);
// Dense callback: void get_pair(int i, vector<Mint>& a, vector<Mint>& b).
// Sparse callback: void get_sparse_pair(int i,
//          vector<pair<int,Mint>>& a, vector<pair<int,Mint>>& b).
// Each callback must produce the SAME vector pair each time it is called.
// Dense outputs are zero-initialized to d entries; sparse outputs are cleared.
// For sparse vectors, repeated indices are allowed (their values are summed).
//
// Optional: seed and trials parameters after restore (trials >= 1).
//   linear_matroid_parity<Mint>(d,m,get_pair,true,seed,3);
//
// MEMORY: O(d^2) auxiliary (input pairs never retained). Callback called multiple times.
// TIME: worst case per randomized trial O(m*d^2+d^3) for dense,
//       O(sum_i(nnz(a_i)*nnz(b_i)) + d^3) for sparse size-only;
//       restore adds O(m*d^2+d^3) in general, but sparse support can help.
//       Preliminary greedy / certificates add work, but can terminate early.
// PROBABILITY: randomized rank evaluation can underestimate OPT; never returns
//       an infeasible solution. The span certificates also have a one-sided
//       Monte Carlo error probability (roughly <= 1/p^2 per invalid test).
//       If restore==false, selected is always empty.
// NOTE: No guaranteed <500ms bound at d=500,m>=100000 with adversarial dense data.
// Optional fast paths use Mint::val() and Mint::raw() when available; the
// general fallback supports Mint without either method. GNU/Clang __int128
// enables faster dot products, with portable modint fallbacks otherwise.

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

struct LinearMatroidParityResult {
    int size = 0;
    std::vector<int> selected;
};

namespace linear_matroid_parity_detail {

inline std::uint64_t default_seed() {
    return static_cast<std::uint64_t>(
        std::chrono::steady_clock::now().time_since_epoch().count());
}

template<class T,class=void> struct HasVal:std::false_type {};
template<class T> struct HasVal<T,std::void_t<decltype(std::declval<T>().val())>>
    :std::true_type {};

template<class T,class=void> struct HasRaw:std::false_type {};
template<class T> struct HasRaw<T,std::void_t<decltype(T::raw(0u))>>
    :std::true_type {};

template<class Mint> Mint from_reduced_u32(std::uint32_t x){
    if constexpr(HasRaw<Mint>::value)return Mint::raw(x);
    else return Mint(x);
}

template<class Mint> using Vec = std::vector<Mint>;
template<class Mint> using SparseVec = std::vector<std::pair<int,Mint>>;

template<class Mint, class Get> struct DenseSource {
    int d;
    static constexpr bool is_sparse = false;
    Get& get;
    Vec<Mint> a, b;
    DenseSource(int d_, Get& g) : d(d_), get(g), a(d), b(d) {}
    void fetch(int id) {
        if ((int)a.size()!=d) a.resize(d);
        if ((int)b.size()!=d) b.resize(d);
        std::fill(a.begin(),a.end(),Mint(0));
        std::fill(b.begin(),b.end(),Mint(0));
        get(id,a,b);
        if ((int)a.size()!=d || (int)b.size()!=d)
            throw std::invalid_argument("linear_matroid_parity: wrong dense vector length");
    }
    void dense(Vec<Mint>& x, Vec<Mint>& y) const { x=a; y=b; }
    template<class Matrix> void add_wedge(Matrix& h, Mint coeff) const {
        for (int i=0;i<d;i++) {
            Mint ai=coeff*a[i], bi=coeff*b[i];
            Mint* row=h.data()+(std::size_t)i*d;
            for (int j=i+1;j<d;j++) row[j]+=ai*b[j]-bi*a[j];
        }
    }
    void projected(const std::vector<int>& idx, const std::vector<int>&,
                   Vec<Mint>& x, Vec<Mint>& y) const {
        for (int i=0;i<(int)idx.size();i++) { x[i]=a[idx[i]]; y[i]=b[idx[i]]; }
    }
};

template<class Mint, class Get> struct SparseSource {
    int d;
    static constexpr bool is_sparse = true;
    Get& get;
    SparseVec<Mint> a,b;
    SparseSource(int d_, Get& g):d(d_),get(g){}
    void fetch(int id) {
        a.clear(); b.clear();
        get(id,a,b);
        for (const auto& p:a) if(p.first<0 || p.first>=d)
            throw std::out_of_range("linear_matroid_parity: bad sparse coordinate");
        for (const auto& p:b) if(p.first<0 || p.first>=d)
            throw std::out_of_range("linear_matroid_parity: bad sparse coordinate");
    }
    void dense(Vec<Mint>& x,Vec<Mint>& y) const {
        x.assign(d,Mint(0));y.assign(d,Mint(0));
        for(const auto& [i,v]:a)x[i]+=v;
        for(const auto& [i,v]:b)y[i]+=v;
    }
    template<class Matrix> void add_wedge(Matrix& h, Mint coeff) const {
        for(const auto& [i,ai]:a) {
            Mint v=coeff*ai;
            for(const auto& [j,bj]:b) if(i!=j) {
                Mint t=v*bj;
                if(i<j)h.at(i,j)+=t;
                else h.at(j,i)-=t;
            }
        }
    }
    void projected(const std::vector<int>&,const std::vector<int>& pos,
                   Vec<Mint>& x,Vec<Mint>& y) const {
        std::fill(x.begin(),x.end(),Mint(0));
        std::fill(y.begin(),y.end(),Mint(0));
        for(const auto& [i,v]:a)if(pos[i]>=0)x[pos[i]]+=v;
        for(const auto& [i,v]:b)if(pos[i]>=0)y[pos[i]]+=v;
    }
    bool singleton_columns(int m){
        for(int i=0;i<m;i++){
            fetch(i);
            if(a.size()>1 || b.size()>1)return false;
        }
        return true;
    }
};

template<class Mint> struct SkewMatrix {
    int n;
    Vec<Mint> v;
    explicit SkewMatrix(int n_):n(n_),v((std::size_t)n_*n_,Mint(0)){}
    Mint& at(int i,int j){return v[(std::size_t)i*n+j];}
    const Mint& at(int i,int j)const{return v[(std::size_t)i*n+j];}
    Mint* data(){return v.data();}
};

template<class Mint> void swap_skew_index(SkewMatrix<Mint>& h,int p,int q){
    if(p==q)return;
    int n=h.n;
    auto get=[&](int i,int j)->Mint{
        return i==j?Mint(0):(i<j?h.at(i,j):-h.at(j,i));
    };
    auto set=[&](int i,int j,Mint value){
        if(i<j)h.at(i,j)=value;
        else if(i>j)h.at(j,i)=-value;
    };
    Vec<Mint> a(n),b(n);
    for(int j=0;j<n;j++)if(j!=p && j!=q){a[j]=get(p,j);b[j]=get(q,j);}
    Mint cross=get(p,q);
    for(int j=0;j<n;j++)if(j!=p && j!=q){set(p,j,b[j]);set(q,j,a[j]);}
    set(p,q,-cross);
}

template<class Mint> std::vector<int> nonsingular_principal_minor(SkewMatrix<Mint> h){
    int d=h.n,r=0;
    std::vector<int> perm(d);
    std::iota(perm.begin(),perm.end(),0);
    for(int k=0;k+1<d;k+=2){
        int p=-1,q=-1;
        if(h.at(k,k+1)!=Mint(0)){p=k;q=k+1;}
        else {
            for(int i=k;i<d && p<0;i++)for(int j=i+1;j<d;j++)
                if(h.at(i,j)!=Mint(0)){p=i;q=j;break;}
            if(p==-1)break;
        }
        if(p!=k){swap_skew_index(h,k,p);std::swap(perm[k],perm[p]);}
        if(q==k)q=p;
        if(q!=k+1){swap_skew_index(h,k+1,q);std::swap(perm[k+1],perm[q]);}
        Mint inv=Mint(1)/h.at(k,k+1);
        if constexpr(HasVal<Mint>::value){
            if(Mint::mod()<=1000000000ULL){
                const std::uint64_t mod=Mint::mod(),offset=mod*mod;
                for(int i=k+2;i<d;i++){
                    Mint ai=h.at(k,i)*inv,bi=h.at(k+1,i)*inv;
                    std::uint64_t av=ai.val(),bv=bi.val();
                    for(int j=i+1;j<d;j++){
                        std::uint64_t z=offset+h.at(i,j).val()
                                     +bv*h.at(k,j).val()-av*h.at(k+1,j).val();
                        h.at(i,j)=from_reduced_u32<Mint>((std::uint32_t)(z%mod));
                    }
                }
            }else{
                for(int i=k+2;i<d;i++){
                    Mint ai=h.at(k,i)*inv,bi=h.at(k+1,i)*inv;
                    for(int j=i+1;j<d;j++)
                        h.at(i,j)+=bi*h.at(k,j)-ai*h.at(k+1,j);
                }
            }
        }else{
            for(int i=k+2;i<d;i++){
                Mint ai=h.at(k,i)*inv,bi=h.at(k+1,i)*inv;
                for(int j=i+1;j<d;j++)
                    h.at(i,j)+=bi*h.at(k,j)-ai*h.at(k+1,j);
            }
        }
        r+=2;
    }
    perm.resize(r);
    return perm;
}

template<class Mint> bool invert(SkewMatrix<Mint>& x){
    const int n=x.n;
    SkewMatrix<Mint> iv(n);
    for(int i=0;i<n;i++)iv.at(i,i)=Mint(1);
    for(int k=0;k<n;k++){
        int p=k;
        while(p<n && x.at(p,k)==Mint(0))p++;
        if(p==n)return false;
        if(p!=k)for(int j=0;j<n;j++){
            std::swap(x.at(p,j),x.at(k,j));
            std::swap(iv.at(p,j),iv.at(k,j));
        }
        Mint mul=Mint(1)/x.at(k,k);
        for(int j=0;j<n;j++){x.at(k,j)*=mul;iv.at(k,j)*=mul;}
        for(int i=0;i<n;i++)if(i!=k && x.at(i,k)!=Mint(0)){
            Mint t=x.at(i,k);
            x.at(i,k)=Mint(0);
            for(int j=k+1;j<n;j++)x.at(i,j)-=t*x.at(k,j);
            for(int j=0;j<n;j++)iv.at(i,j)-=t*iv.at(k,j);
        }
    }
    x=std::move(iv);
    return true;
}

// Inversion of a nonsingular alternating matrix in the pivot order
// supplied by nonsingular_principal_minor. O(n^3) but ~1/6 n^3
// pairwise operations instead of full Gauss-Jordan on [A | I].
// Return false if this order has a zero 2x2 pivot (caller can fall back).
template<class Mint> bool invert_skew_in_pivot_order(SkewMatrix<Mint>& a){
    int n=a.n;
    if(n%2)return false;
    std::vector<Mint> fac((std::size_t)n*n,Mint(0)), piv(n/2);
    for(int k=0;k<n;k+=2){
        Mint t=a.at(k,k+1);
        if(t==Mint(0))return false;
        Mint invt=Mint(1)/t;
        piv[k/2]=invt;
        Mint* f0=fac.data()+(std::size_t)k*n;
        Mint* f1=fac.data()+(std::size_t)(k+1)*n;
        for(int j=k+2;j<n;j++){
            f0[j]=-a.at(k+1,j)*invt;
            f1[j]= a.at(k,j)*invt;
        }
        if constexpr(HasVal<Mint>::value){
            if(Mint::mod()<=1000000000ULL){
                const std::uint64_t mod=Mint::mod(),offset=2*mod*mod;
                for(int i=k+2;i<n;i++){
                    const std::uint64_t x=f0[i].val(),y=f1[i].val();
                    for(int j=i+1;j<n;j++){
                        std::uint64_t z=offset+a.at(i,j).val()
                             -x*a.at(k,j).val()-y*a.at(k+1,j).val();
                        a.at(i,j)=from_reduced_u32<Mint>((std::uint32_t)(z%mod));
                    }
                }
            }else{
                for(int i=k+2;i<n;i++){
                    const Mint x=f0[i],y=f1[i];
                    for(int j=i+1;j<n;j++)
                        a.at(i,j)-=x*a.at(k,j)+y*a.at(k+1,j);
                }
            }
        }else{
            for(int i=k+2;i<n;i++){
                const Mint x=f0[i],y=f1[i];
                for(int j=i+1;j<n;j++)
                    a.at(i,j)-=x*a.at(k,j)+y*a.at(k+1,j);
            }
        }
    }
    SkewMatrix<Mint> iv(n);
    std::vector<Mint> w0(n),w1(n);
    for(int k=n-2;k>=0;k-=2){
        const Mint* f0=fac.data()+(std::size_t)k*n;
        const Mint* f1=fac.data()+(std::size_t)(k+1)*n;
        Mint cross=-piv[k/2];
        for(int j=k+2;j<n;j++){
            Mint s0=0,s1=0;
#if defined(__SIZEOF_INT128__)
            if constexpr(HasVal<Mint>::value){
                if(Mint::mod()<=1000000000ULL){
                    unsigned __int128 acc0=0,acc1=0;
                    for(int i=k+2;i<n;i++){
                        std::uint64_t z=iv.at(i,j).val();
                        acc0+=(std::uint64_t)f0[i].val()*z;
                        acc1+=(std::uint64_t)f1[i].val()*z;
                    }
                    s0=from_reduced_u32<Mint>((std::uint32_t)(acc0%Mint::mod()));
                    s1=from_reduced_u32<Mint>((std::uint32_t)(acc1%Mint::mod()));
                }else{
                    for(int i=k+2;i<n;i++){
                        Mint zij=iv.at(i,j);
                        s0+=f0[i]*zij;s1+=f1[i]*zij;
                    }
                }
            }else{
                for(int i=k+2;i<n;i++){
                    Mint zij=iv.at(i,j);
                    s0+=f0[i]*zij;s1+=f1[i]*zij;
                }
            }
#else
            for(int i=k+2;i<n;i++){
                Mint zij=iv.at(i,j);
                s0+=f0[i]*zij;s1+=f1[i]*zij;
            }
#endif
            w0[j]=s0;w1[j]=s1;
            iv.at(k,j)=-s0;iv.at(j,k)=s0;
            iv.at(k+1,j)=-s1;iv.at(j,k+1)=s1;
        }
#if defined(__SIZEOF_INT128__)
        if constexpr(HasVal<Mint>::value){
            if(Mint::mod()<=1000000000ULL){
                unsigned __int128 acc=0;
                for(int j=k+2;j<n;j++)
                    acc+=(std::uint64_t)w0[j].val()*f1[j].val();
                cross+=from_reduced_u32<Mint>((std::uint32_t)(acc%Mint::mod()));
            }else{
                for(int j=k+2;j<n;j++)cross+=w0[j]*f1[j];
            }
        }else{
            for(int j=k+2;j<n;j++)cross+=w0[j]*f1[j];
        }
#else
        for(int j=k+2;j<n;j++)cross+=w0[j]*f1[j];
#endif
        iv.at(k,k+1)=cross;
        iv.at(k+1,k)=-cross;
    }
    a=std::move(iv);
    return true;
}

// A linear basis, preserving its own state on a failed insertion.
template<class Mint> struct Basis {
    int d, rank=0;
    Vec<Mint> rows;
    std::vector<unsigned char> occupied;
    Vec<Mint> work;
    explicit Basis(int d_):d(d_),rows((std::size_t)d_*d_,Mint(0)),occupied(d_,0),work(d_){}
    int insert(const Vec<Mint>& v){
        work=v;
        for(int p=0;p<d;p++){
            Mint c=work[p];
            if(c==Mint(0))continue;
            Mint* row=rows.data()+(std::size_t)p*d;
            if(!occupied[p]){
                Mint inv=Mint(1)/c;
                row[p]=Mint(1);
                for(int j=p+1;j<d;j++)row[j]=work[j]*inv;
                occupied[p]=1;
                ++rank;
                return p;
            }
            if constexpr(HasVal<Mint>::value){
                if(Mint::mod()<=1000000000ULL){
                    const std::uint64_t mod=Mint::mod(),offset=mod*mod;
                    const std::uint64_t cv=c.val();
                    for(int j=p+1;j<d;j++){
                        std::uint64_t t=offset+work[j].val()-cv*row[j].val();
                        work[j]=from_reduced_u32<Mint>((std::uint32_t)(t%mod));
                    }
                }else{
                    for(int j=p+1;j<d;j++)work[j]-=c*row[j];
                }
            }else{
                for(int j=p+1;j<d;j++)work[j]-=c*row[j];
            }
        }
        return -1;
    }
    void undo(int pivot){
        if(pivot>=0){occupied[pivot]=0;--rank;}
    }
};

// Independent pairs and a deterministic lower bound on OPT.
template<class Mint,class Source>
std::vector<int> greedy_pairs(int d,int m,Source& src,Basis<Mint>& basis){
    Vec<Mint> a(d),b(d);
    std::vector<int> chosen;
    chosen.reserve(d/2);
    // Dense dependent pairs may cost O(d^2) each. The randomized general
    // algorithm handles the tail, so do not let this prepass become worse
    // than its own purpose.
    const int probe_limit=std::min(m,std::max(64,d));
    for(int id=0;id<probe_limit;id++){
        src.fetch(id);
        src.dense(a,b);
        int p=basis.insert(a);
        if(p<0)continue;
        int q=basis.insert(b);
        if(q<0){basis.undo(p);continue;}
        chosen.push_back(id);
        if((int)chosen.size()==d/2)break;
    }
    return chosen;
}

// Safe Monte Carlo upper bound: if every first endpoint lies in the span
// of the first endpoints of the chosen pairs, OPT <= chosen.size().
// Also tested for second endpoints and a+b, a-b.
// Two independently chosen annihilating functionals make a false certificate
// probability <= p^-2 for any counterexample (assuming a prime field).
template<class Mint,class Source>
bool endpoint_span_certificate(int d,int m,Source& src,
    const std::vector<int>& chosen,std::mt19937_64& rng){
    int k=(int)chosen.size();
    if(k==0 && m==0)return true;
    if(k==d/2)return true;
    if(d<=0)return true;
    Vec<Mint> x(d),y(d),z(d);
    for(int flavor=0;flavor<4;flavor++){
        Basis<Mint> basis(d);
        for(int id:chosen){
            src.fetch(id);src.dense(x,y);
            if(flavor==0)z=x;
            else if(flavor==1)z=y;
            else {
                for(int j=0;j<d;j++)z[j]=(flavor==2?x[j]+y[j]:x[j]-y[j]);
            }
            basis.insert(z);
        }
        if(basis.rank!=k)continue;
        Vec<Mint> probe1(d),probe2(d);
        for(int trial=0;trial<2;trial++){
            Vec<Mint>& probe=trial==0?probe1:probe2;
            // Choose a random vector annihilating every row in the basis.
            for(int j=0;j<d;j++)probe[j]=Mint((std::uint32_t)(rng()%Mint::mod()));
            for(int p=d-1;p>=0;p--)if(basis.occupied[p]){
                Mint t=0;
                const Mint* row=basis.rows.data()+(std::size_t)p*d;
                for(int j=p+1;j<d;j++)t+=row[j]*probe[j];
                probe[p]=-t;
            }
        }
        bool good=true;
        for(int id=0;id<m;id++){
            src.fetch(id);
            Mint t1=0,t2=0;
            if constexpr(Source::is_sparse){
#if defined(__SIZEOF_INT128__)
                if constexpr(HasVal<Mint>::value){
                    if(Mint::mod()<=1000000000ULL){
                        unsigned __int128 a1=0,a2=0,b1=0,b2=0;
                        if(flavor!=1)for(const auto& [j,v]:src.a){
                            std::uint64_t x=v.val();
                            a1+=x*probe1[j].val();a2+=x*probe2[j].val();
                        }
                        if(flavor!=0)for(const auto& [j,v]:src.b){
                            std::uint64_t x=v.val();
                            b1+=x*probe1[j].val();b2+=x*probe2[j].val();
                        }
                        auto mod=Mint::mod();
                        if(flavor==0){ t1=Mint((std::uint64_t)(a1%mod));t2=Mint((std::uint64_t)(a2%mod)); }
                        else if(flavor==1){ t1=Mint((std::uint64_t)(b1%mod));t2=Mint((std::uint64_t)(b2%mod)); }
                        else {
                            auto u1=a1%mod,u2=a2%mod,w1=b1%mod,w2=b2%mod;
                            if(flavor==2){ t1=Mint((u1+w1)%mod);t2=Mint((u2+w2)%mod); }
                            if(flavor==3){ t1=Mint((u1+mod-w1)%mod);t2=Mint((u2+mod-w2)%mod); }
                        }
                    } else {
                        if(flavor!=1)for(const auto& [j,v]:src.a){t1+=v*probe1[j];t2+=v*probe2[j];}
                        if(flavor!=0)for(const auto& [j,v]:src.b){
                            Mint vv=flavor==3?-v:v;t1+=vv*probe1[j];t2+=vv*probe2[j];
                        }
                    }
                } else {
                    if(flavor!=1)for(const auto& [j,v]:src.a){t1+=v*probe1[j];t2+=v*probe2[j];}
                    if(flavor!=0)for(const auto& [j,v]:src.b){
                        Mint vv=flavor==3?-v:v;t1+=vv*probe1[j];t2+=vv*probe2[j];
                    }
                }
#else
                if(flavor!=1)for(const auto& [j,v]:src.a){t1+=v*probe1[j];t2+=v*probe2[j];}
                if(flavor!=0)for(const auto& [j,v]:src.b){
                    Mint vv=flavor==3?-v:v;t1+=vv*probe1[j];t2+=vv*probe2[j];
                }
#endif
            }else{
                src.dense(x,y);
                if(flavor==0)z=x;
                else if(flavor==1)z=y;
                else for(int j=0;j<d;j++)z[j]=(flavor==2?x[j]+y[j]:x[j]-y[j]);
                for(int j=0;j<d;j++){
                    t1+=z[j]*probe1[j];t2+=z[j]*probe2[j];
                }
            }
            if(t1!=Mint(0)||t2!=Mint(0)){good=false;break;}
        }
        if(good)return true;
    }
    return false;
}


// Delayed modular reductions. This is safe for the common ACL 998244353
// field: 16 products of residues fit into unsigned 64-bit.
template<class Mint> void matvec(const SkewMatrix<Mint>& a,
             const Vec<Mint>& x,Vec<Mint>& y,std::vector<std::uint64_t>& scratch){
    const int n=a.n;
#if defined(__SIZEOF_INT128__)
    if constexpr(HasVal<Mint>::value){
        if(Mint::mod()<=1000000000ULL){
            if((int)scratch.size()!=n)scratch.resize(n);
            for(int j=0;j<n;j++)scratch[j]=x[j].val();
            const std::uint64_t mod=Mint::mod();
            for(int i=0;i<n;i++){
                const Mint* row=a.v.data()+(std::size_t)i*n;
                unsigned __int128 acc=0;
                for(int j=0;j<n;j++)
                    acc+=(std::uint64_t)row[j].val()*scratch[j];
                y[i]=Mint((std::uint32_t)(acc%mod));
            }
            return;
        }
    }
#endif
    for(int i=0;i<n;i++){
        Mint t=0;
        const Mint* row=a.v.data()+(std::size_t)i*n;
        for(int j=0;j<n;j++)t+=row[j]*x[j];
        y[i]=t;
    }
}

template<class Mint,class Source>
std::vector<int> recover_singleton_sparse(int d,int m,Source& src,
                        const std::vector<int>& idx,SkewMatrix<Mint> iv){
    std::vector<int> active=idx,pos(d,-1),chosen;
    for(int i=0;i<(int)active.size();i++)pos[active[i]]=i;
    chosen.reserve(active.size()/2);
    for(int id=0;id<m && !active.empty();id++){
        src.fetch(id);
        if(src.a.empty() || src.b.empty() ||
           src.a[0].second==Mint(0) || src.b[0].second==Mint(0))continue;
        int p=pos[src.a[0].first],q=pos[src.b[0].first];
        if(p<0 || q<0 || p==q || iv.at(p,q)==Mint(0))continue;
        chosen.push_back(id);
        int n=(int)active.size();
        Mint denom=Mint(1)/iv.at(p,q);
        SkewMatrix<Mint> next(n-2);
        std::vector<int> rem;
        rem.reserve(n-2);
        for(int i=0;i<n;i++)if(i!=p && i!=q)rem.push_back(i);
        for(int i=0;i<n-2;i++)for(int j=i+1;j<n-2;j++){
            int x=rem[i],y=rem[j];
            Mint val=iv.at(x,y)+(iv.at(x,q)*iv.at(y,p)
                                     -iv.at(x,p)*iv.at(y,q))*denom;
            next.at(i,j)=val;
            next.at(j,i)=-val;
        }
        pos[active[p]]=pos[active[q]]=-1;
        std::vector<int> next_active;
        next_active.reserve(n-2);
        for(int i:rem)next_active.push_back(active[i]);
        active=std::move(next_active);
        for(int i=0;i<(int)active.size();i++)pos[active[i]]=i;
        iv=std::move(next);
    }
    return chosen;
}

template<class Mint,class Source>
std::vector<int> recover(int d,int m,Source& src,const SkewMatrix<Mint>& h,
                         const std::vector<int>& idx){
    int r=(int)idx.size();
    if(r==0)return {};
    SkewMatrix<Mint> iv(r);
    for(int i=0;i<r;i++)for(int j=0;j<r;j++){
        int a=idx[i],b=idx[j];
        iv.at(i,j)=a==b?Mint(0):(a<b?h.at(a,b):-h.at(b,a));
    }
    // The principal minor index order already includes all skew pivots.
    if(!invert_skew_in_pivot_order(iv)){
        // A zero block pivot requires a general inverse. Rebuild from the
        // original skew matrix, since the fast attempt modified iv.
        for(int i=0;i<r;i++)for(int j=0;j<r;j++){
            int a=idx[i],b=idx[j];
            iv.at(i,j)=a==b?Mint(0):(a<b?h.at(a,b):-h.at(b,a));
        }
        if(!invert(iv))throw std::runtime_error("linear_matroid_parity: rank/inverse mismatch");
    }
    if constexpr(Source::is_sparse){
        if(src.singleton_columns(m))
            return recover_singleton_sparse<Mint>(d,m,src,idx,std::move(iv));
    }
    std::vector<int> pos(d,-1);
    for(int i=0;i<r;i++)pos[idx[i]]=i;
    Vec<Mint> a(r),b(r),u(r),v(r);
    std::vector<std::uint64_t> scratch(r);
    std::vector<std::pair<int,Mint>> sparse_terms;
    sparse_terms.reserve(32);
    auto sparse_matvec=[&](const SparseVec<Mint>& terms,const Vec<Mint>& dense,Vec<Mint>& result){
        // Sparse vectors are particularly common in partition matroid
        // representations. A matrix-vector product needs only O(r*nnz),
        // not O(r^2), and duplicate coordinates are naturally summed.
        if(terms.size()*8 >= (std::size_t)r){
            matvec(iv,dense,result,scratch);
            return;
        }
        sparse_terms.clear();
        for(const auto& [j,val]:terms){
            int p=pos[j];
            if(p>=0 && val!=Mint(0))sparse_terms.push_back({p,val});
        }
        if constexpr(HasVal<Mint>::value){
            if(Mint::mod()<=1000000000ULL && sparse_terms.size()<=16){
                const std::uint64_t mod=Mint::mod();
                for(int i=0;i<r;i++){
                    const Mint* row=iv.v.data()+(std::size_t)i*r;
                    std::uint64_t sum=0;
                    for(const auto& [p,val]:sparse_terms)
                        sum+=(std::uint64_t)row[p].val()*val.val();
                    result[i]=from_reduced_u32<Mint>((std::uint32_t)(sum%mod));
                }
                return;
            }
        }
        for(int i=0;i<r;i++){
            const Mint* row=iv.v.data()+(std::size_t)i*r;
            Mint sum=0;
            for(const auto& [p,val]:sparse_terms)sum+=row[p]*val;
            result[i]=sum;
        }
    };
    std::vector<int> chosen;
    chosen.reserve(r/2);
    for(int id=0;id<m && (int)chosen.size()<r/2;id++){
        src.fetch(id);
        src.projected(idx,pos,a,b);
        if constexpr(Source::is_sparse)sparse_matvec(src.b,b,v);
        else matvec(iv,b,v,scratch);
        Mint dot=Mint(0);
        for(int i=0;i<r;i++)dot+=a[i]*v[i];
        if(dot==Mint(0))continue;
        if constexpr(Source::is_sparse)sparse_matvec(src.a,a,u);
        else matvec(iv,a,u,scratch);
        Mint invdot=Mint(1)/dot;
        for(int i=0;i<r;i++)u[i]*=invdot;
        // Conditioning rank-2 update. One modular reduction per matrix
        // entry instead of two multiplications + three modint operations.
        if constexpr(HasVal<Mint>::value){
            if(Mint::mod()<=1000000000ULL){
                const std::uint64_t mod=Mint::mod(),mod2=mod*mod;
                std::vector<std::uint32_t> uv(r),vv(r);
                for(int i=0;i<r;i++){uv[i]=u[i].val();vv[i]=v[i].val();}
                for(int i=0;i<r;i++){
                    const std::uint64_t ui=uv[i],vi=vv[i];
                    Mint* row=iv.v.data()+(std::size_t)i*r;
                    for(int j=i+1;j<r;j++){
                        std::uint64_t z=vi*uv[j]+mod2-ui*vv[j]+row[j].val();
                        Mint val=from_reduced_u32<Mint>((std::uint32_t)(z%mod));
                        row[j]=val;
                        iv.v[(std::size_t)j*r+i]=-val;
                    }
                }
            }else{
                for(int i=0;i<r;i++)for(int j=i+1;j<r;j++){
                    Mint z=iv.at(i,j)+v[i]*u[j]-u[i]*v[j];
                    iv.at(i,j)=z;
                    iv.at(j,i)=-z;
                }
            }
        }else{
            for(int i=0;i<r;i++)for(int j=i+1;j<r;j++){
                Mint z=iv.at(i,j)+v[i]*u[j]-u[i]*v[j];
                iv.at(i,j)=z;
                iv.at(j,i)=-z;
            }
        }
        chosen.push_back(id);
    }
    return chosen;
}

// If the selected 2k columns span EVERY input column, no (k+1)-pair
// solution can exist. Two random annihilators check this in O(m*d),
// regardless of the rank of the entire input matrix.
template<class Mint,class Source>
bool full_span_certificate(int d,int m,Source& src,
                           const Basis<Mint>& basis,std::mt19937_64& rng){
    int k=basis.rank/2;
    if(2*k==d)return true;
    Vec<Mint> a(d),b(d);
    Vec<Mint> z1(d),z2(d);
    auto make_probe=[&](Vec<Mint>& z){
        for(int i=0;i<d;i++)z[i]=Mint((std::uint32_t)(rng()%Mint::mod()));
        for(int i=d-1;i>=0;i--)if(basis.occupied[i]){
            Mint val=0;
            const Mint* row=basis.rows.data()+(std::size_t)i*d;
            for(int j=i+1;j<d;j++)val+=row[j]*z[j];
            z[i]=-val;
        }
    };
    make_probe(z1);make_probe(z2);
    // One reduction per dot product (two independent random annihilators).
    // d*(mod-1)^2 fits in 128 bits for ordinary feasible vector lengths.
    std::vector<std::uint64_t> zz1(d),zz2(d);
    if constexpr(HasVal<Mint>::value) {
        if(Mint::mod()<=1000000000ULL)
            for(int j=0;j<d;j++){zz1[j]=z1[j].val();zz2[j]=z2[j].val();}
    }
    auto ann_dense=[&](const Vec<Mint>& vec)->bool{
#if defined(__SIZEOF_INT128__)
        if constexpr(HasVal<Mint>::value){
            if(Mint::mod()<=1000000000ULL){
                unsigned __int128 s1=0,s2=0;
                for(int j=0;j<d;j++){
                    std::uint64_t v=vec[j].val();
                    s1+=v*zz1[j];s2+=v*zz2[j];
                }
                auto mod=Mint::mod();
                return (s1%mod)==0 && (s2%mod)==0;
            }
        }
#endif
        Mint s1=0,s2=0;
        for(int j=0;j<d;j++){s1+=vec[j]*z1[j];s2+=vec[j]*z2[j];}
        return s1==Mint(0)&&s2==Mint(0);
    };
    auto ann_sparse=[&](const SparseVec<Mint>& vec)->bool{
#if defined(__SIZEOF_INT128__)
        if constexpr(HasVal<Mint>::value){
            if(Mint::mod()<=1000000000ULL){
                unsigned __int128 s1=0,s2=0;
                for(const auto& [j,vv]:vec){
                    std::uint64_t v=vv.val();
                    s1+=v*zz1[j];s2+=v*zz2[j];
                }
                auto mod=Mint::mod();
                return (s1%mod)==0 && (s2%mod)==0;
            }
        }
#endif
        Mint s1=0,s2=0;
        for(const auto& [j,v]:vec){s1+=v*z1[j];s2+=v*z2[j];}
        return s1==Mint(0)&&s2==Mint(0);
    };
    for(int id=0;id<m;id++){
        src.fetch(id);
        if constexpr(Source::is_sparse){
            if(!ann_sparse(src.a)||!ann_sparse(src.b))return false;
        }else{
            if(!ann_dense(src.a)||!ann_dense(src.b))return false;
        }
    }
    return true;
}

// k=1 certificate: if every independently usable pair plane contains
// a common line, no two pairs can be jointly independent. The common line
// is reconstructed by intersecting the first pair plane with another one.
// This catches arbitrary changes of basis WITHIN each two-vector pair.
template<class Mint,class Source>
bool common_line_certificate(int d,int m,Source& src,int chosen_id,
                             std::mt19937_64& rng){
    if constexpr(Source::is_sparse) return false;
    else {
        Vec<Mint> a(d),b(d),v(d),u(d),w(d);
        src.fetch(chosen_id);src.dense(a,b);
        bool found=false;
        // Rejection-only pilot includes late and widely separated pairs.
        const int samples=std::min(m,160);
        for(int t=0;t<samples;t++){
            int id=t<32?m-1-t:(int)((std::uint64_t)(t-32)*m/128);
            if(id==chosen_id)continue;
            src.fetch(id);src.dense(u,w);
            std::array<std::array<Mint,4>,4> piv{};
            std::array<unsigned char,4> used{};
            int rank=0;
            for(int j=0;j<d;j++){
                std::array<Mint,4> row{a[j],b[j],-u[j],-w[j]};
                for(int k=0;k<4;k++){
                    if(row[k]==Mint(0))continue;
                    if(!used[k]){
                        Mint inv=Mint(1)/row[k];
                        for(int s=k+1;s<4;s++)row[s]*=inv;
                        row[k]=Mint(1);piv[k]=row;used[k]=1;rank++;
                        break;
                    }
                    Mint factor=row[k];
                    for(int s=k+1;s<4;s++)row[s]-=factor*piv[k][s];
                }
                if(rank==4)break;
            }
            // Both pairs are independent; disjoint planes allow 2 pairs.
            if(rank==4)return false;
            if(rank!=3)continue;
            std::array<Mint,4> coef{};
            for(int j=0;j<4;j++)if(!used[j]){coef[j]=Mint(1);break;}
            for(int j=3;j>=0;j--)if(used[j]){
                Mint val=0;
                for(int k=j+1;k<4;k++)val+=piv[j][k]*coef[k];
                coef[j]=-val;
            }
            bool nonzero=false;
            for(int j=0;j<d;j++){
                v[j]=coef[0]*a[j]+coef[1]*b[j];
                if(v[j]!=Mint(0))nonzero=true;
            }
            if(nonzero){found=true;break;}
        }
        if(!found)return false;
        Vec<Mint> z1(d),z2(d);
        for(int j=0;j<d;j++){
            z1[j]=Mint((std::uint32_t)(rng()%Mint::mod()));
            z2[j]=Mint((std::uint32_t)(rng()%Mint::mod()));
        }
        // Both random probes must certify membership of v in each plane.
        std::vector<std::uint64_t> zz1(d),zz2(d);
        if constexpr(HasVal<Mint>::value){
            if(Mint::mod()<=1000000000ULL)
                for(int j=0;j<d;j++){
                    zz1[j]=z1[j].val();zz2[j]=z2[j].val();
                }
        }
        auto dot2=[&](const Vec<Mint>& x)->std::pair<Mint,Mint>{
#if defined(__SIZEOF_INT128__)
            if constexpr(HasVal<Mint>::value){
                if(Mint::mod()<=1000000000ULL){
                    unsigned __int128 s1=0,s2=0;
                    for(int j=0;j<d;j++){
                        std::uint64_t xj=x[j].val();
                        s1+=xj*zz1[j];s2+=xj*zz2[j];
                    }
                    return {Mint((std::uint64_t)(s1%Mint::mod())),
                            Mint((std::uint64_t)(s2%Mint::mod()))};
                }
            }
#endif
            Mint s1=0,s2=0;
            for(int j=0;j<d;j++){s1+=x[j]*z1[j];s2+=x[j]*z2[j];}
            return {s1,s2};
        };
        auto [vz1,vz2]=dot2(v);
        auto verify_pair=[&](int id)->bool{
            src.fetch(id);src.dense(u,w);
            int p=-1,q=-1;
            for(int j=0;j<d;j++)if(u[j]!=Mint(0)){p=j;break;}
            if(p<0)return true; // zero vector => unusable pair
            Mint determinant=0;
            for(int j=0;j<d;j++)if(j!=p){
                determinant=u[p]*w[j]-u[j]*w[p];
                if(determinant!=Mint(0)){q=j;break;}
            }
            if(q<0)return true; // dependent pair => unusable
            Mint alpha=v[p]*w[q]-v[q]*w[p];
            Mint beta=u[p]*v[q]-u[q]*v[p];
            auto [ua,ub]=dot2(u);
            auto [wa,wb]=dot2(w);
            return alpha*ua+beta*wa==determinant*vz1
                && alpha*ub+beta*wb==determinant*vz2;
        };
        // Quick rejection checks prevent an exceptional late pair from
        // forcing a full scan for an invalid certificate.
        for(int t=0;t<samples;t++){
            int id=t<32?m-1-t:(int)((std::uint64_t)(t-32)*m/128);
            if(!verify_pair(id))return false;
        }
        for(int id=0;id<m;id++)if(!verify_pair(id))return false;
        return true;
    }
}

// Accumulate wedge products over 16 pairs before modular reduction.
// The positive and negative parts are separate: each fits in uint64_t
// when modulus <= 1e9 and the batch has <=16 terms.
template<class Mint, class Source>
void build_skew_matrix(SkewMatrix<Mint>& h,int d,int m,Source& src,
                       std::mt19937_64& rng){
    if constexpr(!Source::is_sparse && HasVal<Mint>::value){
        if(Mint::mod()<=1000000000ULL && d>=24){
            struct Sum {std::uint64_t plus=0,minus=0;};
            std::vector<Sum> sums((std::size_t)d*(d-1)/2);
            std::vector<std::size_t> start(d);
            for(int i=0;i<d;i++)
                start[i]=(std::size_t)i*(2*(std::size_t)d-i-1)/2;
            std::vector<std::uint32_t> a(d),b(d),ca(d),cb(d);
            auto mod=Mint::mod();
            for(int base=0;base<m;base+=16){
                std::fill(sums.begin(),sums.end(),Sum{});
                for(int id=base;id<std::min(m,base+16);id++){
                    src.fetch(id);
                    Mint x=Mint((std::uint32_t)(rng()%(Mint::mod()-1)+1));
                    for(int j=0;j<d;j++){
                        a[j]=src.a[j].val();b[j]=src.b[j].val();
                        ca[j]=(x*src.a[j]).val();cb[j]=(x*src.b[j]).val();
                    }
                    for(int i=0;i<d;i++){
                        Sum* row=sums.data()+start[i];
                        std::uint64_t ai=ca[i],bi=cb[i];
                        for(int j=i+1;j<d;j++){
                            auto& cell=row[j-i-1];
                            cell.plus+=ai*b[j];
                            cell.minus+=bi*a[j];
                        }
                    }
                }
                for(int i=0;i<d;i++){
                    const Sum* row=sums.data()+start[i];
                    for(int j=i+1;j<d;j++){
                        auto &cell=row[j-i-1];
                        auto plus=cell.plus%mod,minus=cell.minus%mod;
                        Mint delta=Mint(plus>=minus?plus-minus:plus+mod-minus);
                        h.at(i,j)+=delta;
                    }
                }
            }
            return;
        }
    }
    for(int id=0;id<m;id++){
        src.fetch(id);
        Mint x=Mint((std::uint32_t)(rng()%(Mint::mod()-1)+1));
        src.add_wedge(h,x);
    }
}

template<class Mint,class Source>
LinearMatroidParityResult solve(int d,int m,Source& src,bool restore,
                               std::uint64_t seed,int trials){
    if(d<0 || m<0 || trials<=0 || Mint::mod()<=2)
        throw std::invalid_argument("linear_matroid_parity: invalid d/m/trials/field");
    LinearMatroidParityResult best;
    if(d<2 || m==0)return best;
    std::mt19937_64 rng(seed);
    // A cheap deterministic lower bound and randomized upper-bound certificates.
    Basis<Mint> greedy_basis(d);
    auto greedy=greedy_pairs<Mint>(d,m,src,greedy_basis);
    best.size=(int)greedy.size();
    if(restore)best.selected=greedy;
    if(best.size==d/2 || full_span_certificate<Mint>(d,m,src,greedy_basis,rng)
       || (best.size==1 && common_line_certificate<Mint>(d,m,src,greedy.front(),rng))
       || endpoint_span_certificate<Mint>(d,m,src,greedy,rng))return best;
    for(int tr=0;tr<trials;tr++){
        SkewMatrix<Mint> h(d);
        build_skew_matrix(h,d,m,src,rng);
        auto idx=nonsingular_principal_minor(h);
        if((int)idx.size()/2 <= best.size)continue;
        if(!restore){best.size=(int)idx.size()/2;}
        else {
            auto solution=recover<Mint>(d,m,src,h,idx);
            if((int)solution.size()>best.size){
                best.size=(int)solution.size();
                best.selected=std::move(solution);
            }
        }
        if(best.size==d/2)break;
    }
    return best;
}
} // namespace linear_matroid_parity_detail

// Dense callback version: get(i, vector<Mint>& a, vector<Mint>& b)
template<class Mint,class Get>
LinearMatroidParityResult linear_matroid_parity(
    int d,int m,Get&& get,bool restore=false,
    std::uint64_t seed=linear_matroid_parity_detail::default_seed(),int trials=1){
    if(d<0 || m<0 || trials<=0)throw std::invalid_argument("linear_matroid_parity: invalid dimensions/trials");
    using G=std::remove_reference_t<Get>;
    linear_matroid_parity_detail::DenseSource<Mint,G> src(d,get);
    return linear_matroid_parity_detail::solve<Mint>(d,m,src,restore,seed,trials);
}

// Sparse callback version: get(i, vector<pair<int,Mint>>& a,
//                                vector<pair<int,Mint>>& b)
template<class Mint,class Get>
LinearMatroidParityResult linear_matroid_parity_sparse(
    int d,int m,Get&& get,bool restore=false,
    std::uint64_t seed=linear_matroid_parity_detail::default_seed(),int trials=1){
    if(d<0 || m<0 || trials<=0)throw std::invalid_argument("linear_matroid_parity: invalid dimensions/trials");
    using G=std::remove_reference_t<Get>;
    linear_matroid_parity_detail::SparseSource<Mint,G> src(d,get);
    return linear_matroid_parity_detail::solve<Mint>(d,m,src,restore,seed,trials);
}
