#pragma once
/** ll-only ax+by=c, count solutions with lx<=x<=rx and ly<=y<=ry.
 * Use Diophantine::count(a,b,c,lx,rx,ly,ry).
 * Internal arithmetic uses __int128; throws if solution count exceeds ll.
 */
struct Diophantine {
private:
    using W=__int128_t;
    static W floor_div(W a,W b){W q=a/b,r=a%b;return q-(r!=0&&((r>0)!=(b>0)));}
    static W ceil_div(W a,W b){W q=a/b,r=a%b;return q+(r!=0&&((r>0)==(b>0)));}
    static W egcd(W a,W b,W& x,W& y){
        if(b==0){x=1;y=0;return a;}
        W x1,y1;W g=egcd(b,a%b,x1,y1);x=y1;y=x1-(a/b)*y1;return g;
    }
    static long long narrow(W x){
        if(x>std::numeric_limits<long long>::max())throw std::overflow_error("solution count exceeds ll");
        return (long long)x;
    }
public:
    static long long count(long long a0,long long b0,long long c0,long long lx0,long long rx0,long long ly0,long long ry0){
        W a=a0,b=b0,c=c0,lx=lx0,rx=rx0,ly=ly0,ry=ry0;
        if(lx>rx||ly>ry)return 0;
        if(a==0&&b==0){
            if(c!=0)return 0;
            W product;
            if(__builtin_mul_overflow(rx-lx+1,ry-ly+1,&product))throw std::overflow_error("solution count exceeds ll");
            return narrow(product);
        }
        if(a==0)return c%b==0&&ly<=c/b&&c/b<=ry?narrow(rx-lx+1):0;
        if(b==0)return c%a==0&&lx<=c/a&&c/a<=rx?narrow(ry-ly+1):0;
        W x,y,g=egcd(a<0?-a:a,b<0?-b:b,x,y);
        if(c%g)return 0;
        if(a<0)x=-x;if(b<0)y=-y;x*=c/g;y*=c/g;
        W lo=std::numeric_limits<W>::min(),hi=std::numeric_limits<W>::max();
        auto bound=[&](W v,W coeff,W L,W R){
            if(coeff>0){lo=std::max(lo,ceil_div(L-v,coeff));hi=std::min(hi,floor_div(R-v,coeff));}
            else {lo=std::max(lo,ceil_div(R-v,coeff));hi=std::min(hi,floor_div(L-v,coeff));}
        };
        bound(x,b/g,lx,rx);bound(y,-a/g,ly,ry);
        if(lo>hi)return 0;
        return narrow(hi-lo+1);
    }
};
