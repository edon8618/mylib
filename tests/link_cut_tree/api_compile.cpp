#include "../../LinkCutTree.hpp"
#include <cassert>
#include <type_traits>
long long sum(long long a,long long b){return a+b;} long long zero(){return 0;}
using Tree=LinkCutTree<long long,sum,zero>;
template<class T> concept HasThreeArgumentEndpointSet = requires {
    static_cast<void(T::*)(int,int,long long)>(&T::set_edge);
};
template<class T> concept HasThreeArgumentLink = requires(T& t) { t.link_edge(0,1,1LL); };
static_assert(!HasThreeArgumentEndpointSet<Tree>);
static_assert(!HasThreeArgumentLink<Tree>);
static_assert(std::is_same_v<decltype(static_cast<void(Tree::*)(int,long long,long long)>(&Tree::set_edge)),void(Tree::*)(int,long long,long long)>);
// A monoid without a default constructor or equality operator.
struct NoDefault { int x; NoDefault()=delete; explicit NoDefault(int x):x(x){}; };
NoDefault op(NoDefault a,NoDefault b){return NoDefault(a.x+b.x);}
NoDefault e(){return NoDefault(0);}
int main(){
    LinkCutTree<NoDefault,op,e> t(2);
    t.link(0,1); int id=t.get_edge_id(0,1);
    t.set_edge(id,NoDefault(2),NoDefault(3));
    assert(t.prod(0,1).x==2 && t.prod(1,0).x==3);
}
