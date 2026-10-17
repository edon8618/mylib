#pragma once
#include <bits/stdc++.h>

// Integer geometry: all predicates and squared lengths are exact.
// Coordinates / resulting coordinates must fit in long long.
// Products and area2 are computed with __int128 (assuming they fit there).
// Every function is in namespace geo_int; the header is standalone.
namespace geo_int {
using ll = long long;
using i128 = __int128_t;
using ld = long double;

// Useful for exact area / rational answers.
i128 abs128(i128 x) { return x < 0 ? -x : x; }
i128 gcd128(i128 a, i128 b) {
    a=abs128(a); b=abs128(b);
    while (b) { i128 r=a%b; a=b; b=r; }
    return a;
}
std::string toString(i128 x) {
    if (x==0) return "0";
    bool neg=x<0;
    // uint128 avoids overflow when x is the smallest signed 128-bit value.
    __uint128_t u=neg ? (__uint128_t)(-(x+1))+1 : (__uint128_t)x;
    std::string s;
    while(u){s.push_back(char('0'+u%10));u/=10;}
    if(neg)s.push_back('-');
    std::reverse(s.begin(),s.end());
    return s;
}

struct P {
    ll x = 0, y = 0;
    P() = default;
    P(ll x_, ll y_) : x(x_), y(y_) {}
    P operator+(P b) const { return {x + b.x, y + b.y}; }
    P operator-(P b) const { return {x - b.x, y - b.y}; }
    P operator-() const { return {-x, -y}; }
    P operator*(ll k) const { return {x * k, y * k}; }
    friend P operator*(ll k, P a) { return a * k; }
    bool operator==(P b) const { return x == b.x && y == b.y; }
    bool operator!=(P b) const { return !(*this == b); }
    bool operator<(P b) const { return x != b.x ? x < b.x : y < b.y; }
};
std::istream& operator>>(std::istream& in, P& p) { return in >> p.x >> p.y; }
std::ostream& operator<<(std::ostream& out, P p) { return out << p.x << ' ' << p.y; }

int sgn(i128 x) { return (x > 0) - (x < 0); }
i128 dot(P a, P b) { return (i128)a.x * b.x + (i128)a.y * b.y; }
i128 cross(P a, P b) { return (i128)a.x * b.y - (i128)a.y * b.x; }
i128 cross(P a, P b, P c) { return cross(b - a, c - a); }
i128 norm2(P a) { return dot(a, a); }
i128 dist2(P a, P b) { return norm2(a - b); }
// +1: ccw, 0: collinear, -1: cw
int orient(P a, P b, P c) { return sgn(cross(a,b,c)); }

bool onSegment(P a, P b, P p) { // closed segment, a==b supported
    return orient(a,b,p) == 0 && dot(p-a,p-b) <= 0;
}
bool segmentsIntersect(P a, P b, P c, P d) { // closed segments, includes touching
    int x = orient(a,b,c), y = orient(a,b,d);
    int z = orient(c,d,a), w = orient(c,d,b);
    if (x == 0 && onSegment(a,b,c)) return true;
    if (y == 0 && onSegment(a,b,d)) return true;
    if (z == 0 && onSegment(c,d,a)) return true;
    if (w == 0 && onSegment(c,d,b)) return true;
    return x * y < 0 && z * w < 0;
}

using Poly = std::vector<P>;
// CCW, with no redundant collinear vertices by default.
// 0/1 distinct point -> 0/1 vertex; all collinear -> two endpoints.
// keepCollinear=true includes all input boundary points (all-collinear: sorted points).
Poly convexHull(Poly pts, bool keepCollinear = false) { // O(n log n)
    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
    int n = (int)pts.size();
    if (n <= 2) return pts;
    bool allCollinear = true;
    for (int i = 2; i < n; ++i)
        if (orient(pts[0],pts[1],pts[i])) { allCollinear = false; break; }
    if (allCollinear) return keepCollinear ? pts : Poly{pts.front(), pts.back()};
    Poly lo, hi;
    for (P p : pts) {
        while (lo.size() >= 2 && (keepCollinear ? orient(lo[lo.size()-2],lo.back(),p) < 0
                                                   : orient(lo[lo.size()-2],lo.back(),p) <= 0)) lo.pop_back();
        lo.push_back(p);
    }
    for (int i = n-1; i >= 0; --i) {
        P p = pts[i];
        while (hi.size() >= 2 && (keepCollinear ? orient(hi[hi.size()-2],hi.back(),p) < 0
                                                   : orient(hi[hi.size()-2],hi.back(),p) <= 0)) hi.pop_back();
        hi.push_back(p);
    }
    lo.pop_back(); hi.pop_back();
    lo.insert(lo.end(), hi.begin(), hi.end());
    return lo;
}

i128 area2(const Poly& poly) { // signed twice-area; CCW positive
    i128 s = 0;
    for (int i=0,n=(int)poly.size(); i<n; ++i)
        s += cross(poly[i],poly[(i+1)%n]);
    return s;
}
ld perimeter(const Poly& poly) {
    ld s = 0;
    for (int i=0,n=(int)poly.size(); i<n; ++i)
        s += sqrtl((ld)dist2(poly[i], poly[(i+1)%n]));
    return s;
}
// Any simple polygon, arbitrary orientation. 0: outside, 1: boundary, 2: inside.
int pointInPolygon(const Poly& poly, P p) { // O(n)
    int n = (int)poly.size();
    if (n == 1) return p == poly[0] ? 1 : 0;
    if (n == 2) return onSegment(poly[0],poly[1],p) ? 1 : 0;
    bool inside = false;
    for (int i=0; i<n; ++i) {
        P a=poly[i], b=poly[(i+1)%n];
        if (onSegment(a,b,p)) return 1;
        if ((a.y<=p.y && p.y<b.y && orient(a,b,p)>0) ||
            (b.y<=p.y && p.y<a.y && orient(a,b,p)<0)) inside = !inside;
    }
    return inside ? 2 : 0;
}
// Strictly convex, CCW polygon (e.g. convexHull(pts)), no collinear vertices.
// 0: outside, 1: boundary, 2: inside.
int pointInConvex(const Poly& h, P p) { // O(log n)
    int n = (int)h.size();
    if (n <= 2) return pointInPolygon(h,p);
    int a=orient(h[0],h[1],p), b=orient(h[0],h[n-1],p);
    if (a < 0 || b > 0) return 0;
    if (a == 0) return onSegment(h[0],h[1],p) ? 1 : 0;
    if (b == 0) return onSegment(h[0],h[n-1],p) ? 1 : 0;
    int lo=1, hi=n-1;
    while (hi-lo>1) {
        int mid=(lo+hi)/2;
        if (orient(h[0],h[mid],p)>=0) lo=mid;
        else hi=mid;
    }
    int t=orient(h[lo],h[lo+1],p);
    if (t<0) return 0;
    if (t==0) return onSegment(h[lo],h[lo+1],p) ? 1 : 0;
    return 2;
}

// Square of the diameter of a strictly convex CCW polygon (convexHull result).
i128 diameter2(const Poly& h) { // O(n)
    int n=(int)h.size();
    if (n<2) return 0;
    if (n==2) return dist2(h[0],h[1]);
    i128 ans=0;
    int j=1;
    for (int i=0; i<n; ++i) {
        int ni=(i+1)%n;
        while (cross(h[ni]-h[i],h[(j+1)%n]-h[j])>0) j=(j+1)%n;
        ans=std::max({ans,dist2(h[i],h[j]),dist2(h[ni],h[j])});
    }
    return ans;
}

// Convex Minkowski sum. Inputs: strict CCW convex polygons, possibly 1 or 2 points.
// Output: strict CCW polygon. Runs in O(|a|+|b|). Empty input -> empty output.
// For arbitrary point sets, first call convexHull on each.
Poly minkowskiSum(Poly a, Poly b) {
    if (a.empty() || b.empty()) return {};
    if (a.size()==1) { for (P &p:b) p=p+a[0]; return b; }
    if (b.size()==1) { for (P &p:a) p=p+b[0]; return a; }
    auto reorder=[](Poly& v) {
        int k=0;
        for (int i=1;i<(int)v.size();++i)
            if (v[i].y<v[k].y || (v[i].y==v[k].y && v[i].x<v[k].x)) k=i;
        std::rotate(v.begin(),v.begin()+k,v.end());
    };
    reorder(a); reorder(b);
    int n=a.size(),m=b.size(),i=0,j=0;
    Poly result;
    result.reserve(n+m);
    while (i<n || j<m) {
        result.push_back(a[i%n]+b[j%m]);
        if (i==n) { ++j; continue; }
        if (j==m) { ++i; continue; }
        i128 z=cross(a[(i+1)%n]-a[i%n],b[(j+1)%m]-b[j%m]);
        if (z>=0) ++i;
        if (z<=0) ++j;
    }
    // Remove collinear intermediate vertices / duplicate cyclic endpoints in linear time.
    Poly out;
    for (P p: result) {
        if (!out.empty() && p==out.back()) continue;
        while (out.size()>=2 && orient(out[out.size()-2],out.back(),p)==0) out.pop_back();
        out.push_back(p);
    }
    while (out.size()>=3 && orient(out[out.size()-2],out.back(),out[0])==0) out.pop_back();
    size_t start=0;
    while (out.size()-start>=3 && orient(out.back(),out[start],out[start+1])==0) ++start;
    if (start) out.erase(out.begin(),out.begin()+start);
    if (out.size()>1 && out.front()==out.back()) out.pop_back();
    if (out.size()>=3 && area2(out)==0) {
        auto [mn,mx]=std::minmax_element(out.begin(),out.end());
        return {*mn,*mx};
    }
    return out;
}

// Closest pair of points: exact squared distance. Returns nullopt for <2 points.
std::optional<i128> closestPair2(Poly pts) { // O(n log n)
    if (pts.size()<2) return std::nullopt;
    std::sort(pts.begin(),pts.end());
    Poly buf(pts.size());
    const i128 INF = (i128(1)<<125); // requires squared distances < this bound
    auto rec=[&](auto&& self,int l,int r)->i128 {
        if (r-l<=3) {
            i128 d=INF;
            for (int i=l;i<r;++i) for (int j=i+1;j<r;++j)
                d=std::min(d,dist2(pts[i],pts[j]));
            std::sort(pts.begin()+l,pts.begin()+r,[](P a,P b){
                return a.y!=b.y ? a.y<b.y : a.x<b.x;
            });
            return d;
        }
        int mid=(l+r)/2;
        ll mx=pts[mid].x;
        i128 d=std::min(self(self,l,mid),self(self,mid,r));
        std::merge(pts.begin()+l,pts.begin()+mid,pts.begin()+mid,pts.begin()+r,
                   buf.begin()+l,[](P a,P b){return a.y!=b.y ? a.y<b.y : a.x<b.x;});
        std::copy(buf.begin()+l,buf.begin()+r,pts.begin()+l);
        Poly strip;
        for (int i=l;i<r;++i) {
            P p=pts[i];
            if ((i128)(p.x-mx)*(p.x-mx)>=d) continue;
            for (int j=(int)strip.size()-1;j>=0;--j) {
                if ((i128)(p.y-strip[j].y)*(p.y-strip[j].y)>=d) break;
                d=std::min(d,dist2(p,strip[j]));
            }
            strip.push_back(p);
        }
        return d;
    };
    return rec(rec,0,(int)pts.size());
}
} // namespace geo_int
