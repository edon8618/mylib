#pragma once
#include <bits/stdc++.h>

// Floating point geometry: long double for constructions involving division,
// circles and distances. EPS is absolute; adjust for the problem's scale.
// Standalone header; separate namespace from geo_int.hpp.
namespace geo_real {
using ld = long double;
constexpr ld EPS = 1e-12L;
constexpr ld PI = acosl(-1.L);
int sgn(ld x) { return (x > EPS) - (x < -EPS); }
bool eq(ld a, ld b) { return fabsl(a-b)<=EPS; }

struct P {
    ld x=0, y=0;
    P()=default;
    P(ld x_,ld y_): x(x_),y(y_) {}
    P operator+(P b) const {return {x+b.x,y+b.y};}
    P operator-(P b) const {return {x-b.x,y-b.y};}
    P operator-() const {return {-x,-y};}
    P operator*(ld k) const {return {x*k,y*k};}
    P operator/(ld k) const {return {x/k,y/k};}
    friend P operator*(ld k,P a){return a*k;}
    bool operator==(P b) const {return x==b.x&&y==b.y;} // exact, for sorting/unique only
    bool operator<(P b) const {return x!=b.x?x<b.x:y<b.y;} // exact strict weak order
};
std::istream& operator>>(std::istream& in,P& p){return in>>p.x>>p.y;}
std::ostream& operator<<(std::ostream& out,P p){return out<<p.x<<' '<<p.y;}
ld dot(P a,P b){return a.x*b.x+a.y*b.y;}
ld cross(P a,P b){return a.x*b.y-a.y*b.x;}
ld cross(P a,P b,P c){return cross(b-a,c-a);}
ld norm2(P p){return dot(p,p);}
ld len(P p){return sqrtl(norm2(p));}
ld dist(P a,P b){return len(a-b);}
int orient(P a,P b,P c){return sgn(cross(a,b,c));}
bool samePoint(P a,P b){return len(a-b)<=EPS;}
P rotate90(P p){return {-p.y,p.x};}
P rotate(P p,ld theta){ld c=cosl(theta),s=sinl(theta);return {c*p.x-s*p.y,s*p.x+c*p.y};}

struct Line {P a,b;}; // infinite line through distinct points
struct Segment {P a,b;}; // closed segment, a==b permitted
struct Circle {P c; ld r;}; // r>=0
using Poly=std::vector<P>;

P projection(Line l,P p){P d=l.b-l.a;return l.a+d*(dot(p-l.a,d)/norm2(d));}
P reflection(Line l,P p){return projection(l,p)*2-p;}
// Unique intersection only: nullopt if parallel OR coincident.
std::optional<P> lineIntersection(Line a,Line b){
    P u=a.b-a.a,v=b.b-b.a;
    ld den=cross(u,v);
    if(fabsl(den)<=EPS*len(u)*len(v)) return std::nullopt;
    return a.a+u*(cross(b.a-a.a,v)/den);
}
bool onSegment(P a,P b,P p){
    return orient(a,b,p)==0 && sgn(dot(p-a,p-b))<=0;
}
bool segmentsIntersect(Segment s,Segment t){
    int a=orient(s.a,s.b,t.a),b=orient(s.a,s.b,t.b);
    int c=orient(t.a,t.b,s.a),d=orient(t.a,t.b,s.b);
    if(a==0&&onSegment(s.a,s.b,t.a)) return true;
    if(b==0&&onSegment(s.a,s.b,t.b)) return true;
    if(c==0&&onSegment(t.a,t.b,s.a)) return true;
    if(d==0&&onSegment(t.a,t.b,s.b)) return true;
    return a*b<0&&c*d<0;
}
ld distance(Line l,P p){return fabsl(cross(l.b-l.a,p-l.a))/dist(l.a,l.b);}
ld distance(Segment s,P p){
    P d=s.b-s.a;
    if(sgn(norm2(d))==0)return dist(s.a,p);
    ld t=std::clamp(dot(p-s.a,d)/norm2(d),0.L,1.L);
    return dist(s.a+d*t,p);
}
ld distance(Segment a,Segment b){
    if(segmentsIntersect(a,b))return 0;
    return std::min({distance(a,b.a),distance(a,b.b),distance(b,a.a),distance(b,a.b)});
}

// Convex hull, CCW, no collinear vertices (all collinear => endpoints).
Poly convexHull(Poly pts){ // O(n log n)
    std::sort(pts.begin(),pts.end());
    pts.erase(std::unique(pts.begin(),pts.end()),pts.end());
    if(pts.size()<=2) return pts;
    Poly lo,hi;
    for(P p:pts){
        while(lo.size()>=2&&orient(lo[lo.size()-2],lo.back(),p)<=0)lo.pop_back();
        lo.push_back(p);
    }
    for(int i=(int)pts.size()-1;i>=0;--i){
        P p=pts[i];
        while(hi.size()>=2&&orient(hi[hi.size()-2],hi.back(),p)<=0)hi.pop_back();
        hi.push_back(p);
    }
    lo.pop_back();hi.pop_back();
    lo.insert(lo.end(),hi.begin(),hi.end());
    return lo;
}
ld signedArea(const Poly& p){ // signed, positive for CCW
    ld s=0;
    for(int i=0,n=p.size();i<n;++i)s+=cross(p[i],p[(i+1)%n]);
    return s/2;
}
ld area(const Poly& p){return fabsl(signedArea(p));}
// Arbitrary simple polygon. 0: outside, 1: boundary, 2: inside.
int pointInPolygon(const Poly& poly,P p){
    int n=poly.size();
    if(n==1)return samePoint(poly[0],p)?1:0;
    if(n==2)return onSegment(poly[0],poly[1],p)?1:0;
    bool inside=false;
    for(int i=0;i<n;++i){
        P a=poly[i],b=poly[(i+1)%n];
        if(onSegment(a,b,p))return 1;
        if((a.y<=p.y&&p.y<b.y&&orient(a,b,p)>0)||
           (b.y<=p.y&&p.y<a.y&&orient(a,b,p)<0))inside=!inside;
    }
    return inside?2:0;
}
// Strict convex CCW input only (e.g. convexHull). 0/1/2 same meaning.
int pointInConvex(const Poly& h,P p){ // O(log n)
    int n=h.size();
    if(n<=2)return pointInPolygon(h,p);
    int a=orient(h[0],h[1],p),b=orient(h[0],h[n-1],p);
    if(a<0||b>0)return 0;
    if(a==0)return onSegment(h[0],h[1],p)?1:0;
    if(b==0)return onSegment(h[0],h[n-1],p)?1:0;
    int lo=1,hi=n-1;
    while(hi-lo>1){int m=(lo+hi)/2;if(orient(h[0],h[m],p)>=0)lo=m;else hi=m;}
    int c=orient(h[lo],h[lo+1],p);
    if(c<0)return 0;
    if(c==0)return onSegment(h[lo],h[lo+1],p)?1:0;
    return 2;
}
// Clip a CCW convex polygon by the LEFT closed half-plane of the directed line.
Poly convexCut(const Poly& poly,Line l){ // O(n)
    Poly ans;
    int n=poly.size();
    for(int i=0;i<n;++i){
        P a=poly[i],b=poly[(i+1)%n];
        ld da=cross(l.b-l.a,a-l.a),db=cross(l.b-l.a,b-l.a);
        if(sgn(da)>=0) ans.push_back(a);
        if(sgn(da)*sgn(db)<0) ans.push_back(a+(b-a)*(da/(da-db)));
    }
    Poly dedup;
    for(P p:ans)if(dedup.empty()||!samePoint(dedup.back(),p))dedup.push_back(p);
    if(dedup.size()>1&&samePoint(dedup.front(),dedup.back()))dedup.pop_back();
    return dedup;
}
// Both polygons must be CCW convex and nondegenerate (>=3 vertices).
Poly convexIntersection(Poly a,const Poly& b){ // O(|a|*|b|)
    if(a.empty()||b.empty())return {};
    for(int i=0,n=b.size();i<n&&!a.empty();++i)
        a=convexCut(a,{b[i],b[(i+1)%n]});
    return a;
}

// Circumcircle of 2 points, or of 3 points (fallback: diameter for collinear).
Circle circle2(P a,P b){P c=(a+b)/2;return {c,dist(c,a)};}
Circle circle3(P a,P b,P c){
    P u=b-a,v=c-a;
    ld d=2*cross(u,v);
    if(sgn(d)==0){
        if(norm2(a-b)>=norm2(a-c)&&norm2(a-b)>=norm2(b-c))return circle2(a,b);
        if(norm2(a-c)>=norm2(b-c))return circle2(a,c);
        return circle2(b,c);
    }
    P center=a+P{(norm2(u)*v.y-norm2(v)*u.y)/d,
                  (u.x*norm2(v)-v.x*norm2(u))/d};
    return {center,dist(center,a)};
}
bool contains(Circle c,P p){return dist(c.c,p)<=c.r+EPS;}
// Randomized incremental smallest enclosing circle; expected O(n), O(n^3) worst.
Circle minimumEnclosingCircle(std::vector<P> p){
    if(p.empty())return {{0,0},0};
    std::mt19937_64 rng(712367821);
    std::shuffle(p.begin(),p.end(),rng);
    Circle c={p[0],0};
    for(int i=1;i<(int)p.size();++i)if(!contains(c,p[i])){
        c={p[i],0};
        for(int j=0;j<i;++j)if(!contains(c,p[j])){
            c=circle2(p[i],p[j]);
            for(int k=0;k<j;++k)if(!contains(c,p[k]))
                c=circle3(p[i],p[j],p[k]);
        }
    }
    return c;
}
// Finite intersection points, 0-2. Coincident circles (infinitely many) -> {}.
std::vector<P> circleLineIntersection(Circle c,Line l){
    P q=projection(l,c.c);
    ld d2=norm2(q-c.c),h2=c.r*c.r-d2;
    if(h2 < -EPS)return {};
    if(h2 <= EPS)return {q};
    P unit=(l.b-l.a)/dist(l.a,l.b);
    ld h=sqrtl(std::max(0.L,h2));
    return {q+unit*h,q-unit*h};
}
std::vector<P> circleCircleIntersection(Circle a,Circle b){
    ld d=dist(a.c,b.c);
    if(d<=EPS)return {}; // concentric: none or infinitely many
    if(d>a.r+b.r+EPS||d<fabsl(a.r-b.r)-EPS)return {};
    ld x=(a.r*a.r-b.r*b.r+d*d)/(2*d);
    ld h2=a.r*a.r-x*x;
    if(h2< -EPS)return {};
    P foot=a.c+(b.c-a.c)*(x/d);
    if(h2<=EPS)return {foot};
    P h=rotate90((b.c-a.c)/d)*sqrtl(std::max(0.L,h2));
    return {foot+h,foot-h};
}
// Tangent contact points on circle, 0/1/2.
std::vector<P> tangentPoints(Circle c,P p){
    P u=p-c.c;
    ld d2=norm2(u),r2=c.r*c.r;
    if(d2<r2-EPS)return {};
    if(d2<=r2+EPS)return {p};
    P base=c.c+u*(r2/d2);
    P h=rotate90(u)*(c.r*sqrtl(std::max(0.L,d2-r2))/d2);
    return {base+h,base-h};
}
// Common tangents as valid Lines. First point is on circle a;
// second is the contact point on b, except at tangency where it is first + direction.
// Coincident circles have infinitely many tangents -> empty.
std::vector<Line> commonTangents(Circle a,Circle b){
    P d=b.c-a.c;
    ld z=norm2(d);
    if(z<=EPS*EPS)return {};
    std::vector<Line> ans;
    for(int s:{-1,1}){
        ld r=a.r-s*b.r;
        ld h2=z-r*r;
        if(h2< -EPS)continue;
        if(h2<0)h2=0;
        for(int sign:{-1,1}){
            P v=(d*r+rotate90(d)*(sqrtl(h2)*sign))/z;
            P pa = a.c + v * a.r;
            P pb = b.c + v * (s * b.r);
            if (h2 <= EPS) pb = pa + rotate90(v); // contact point + tangent direction
            ans.push_back({pa, pb});
            if(h2<=EPS)break;
        }
    }
    return ans;
}
// Area of intersection between simple polygon (CW or CCW) and circle, O(n).
ld circlePolygonIntersectionArea(const Poly& poly,Circle c){
    ld total=0;
    int n=poly.size();
    for(int i=0;i<n;++i){
        P a=poly[i]-c.c,b=poly[(i+1)%n]-c.c,d=b-a;
        ld aa=norm2(d);
        if(aa<=EPS*EPS)continue;
        std::vector<ld> ts={0,1};
        ld bb=dot(a,d),cc=norm2(a)-c.r*c.r;
        ld det=bb*bb-aa*cc;
        if(det>EPS){
            ld sq=sqrtl(det);
            ld t1=(-bb-sq)/aa,t2=(-bb+sq)/aa;
            if(EPS<t1&&t1<1-EPS)ts.push_back(t1);
            if(EPS<t2&&t2<1-EPS)ts.push_back(t2);
        }
        std::sort(ts.begin(),ts.end());
        for(int j=0;j+1<(int)ts.size();++j){
            if(ts[j+1]-ts[j]<=EPS)continue;
            P u=a+d*ts[j],v=a+d*ts[j+1];
            P mid=(u+v)/2;
            if(norm2(mid)<=c.r*c.r+EPS)total+=cross(u,v)/2;
            else total+=c.r*c.r*atan2l(cross(u,v),dot(u,v))/2;
        }
    }
    return fabsl(total);
}
} // namespace geo_real
