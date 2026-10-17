# 競プロ用・幾何ライブラリ（個人用・コピペ優先）

単体で貼れる2本。`#include <bits/stdc++.h>` を含み、C++17以降で動作。`using namespace std;` は不要。両方同時に読み込んでも衝突しない。

- **geo_int.hpp**（最初に使う）：整数座標の向き判定・交差・凸包・内外判定・面積・直径・Minkowski 和・最近点対。`geo_int::P` は `long long` 座標で、外積・内積・面積の2倍・距離の2乗は `__int128`。厳密判定。
- **geo_real.hpp**（割り算、円、接線など）：`geo_real::P` は `long double` 座標。`EPS=1e-12L`。近接判定には絶対誤差があるため、問題のスケールに合わせて調整。交点が1つに決まらない場合は一部 `nullopt` や空配列を返す。

## 使い始める順番

1. 入力が整数なら、できるだけ整数版を選ぶ（浮動小数点誤差を避ける）。
2. `P{3, 4}` で点を作り、`p.x`, `p.y` で座標を見る。`p+q`, `p-q`, `p*k` はベクトル演算。
3. `cross(b-a,c-a)` または `orient(a,b,c)` で向きを取る（反時計回りなら正）。
4. 点集合なら `convexHull(pts)` を呼ぶ。返り値は反時計回りで、一直線上の余分な点は取り除かれる。
5. 角度比較はなるべく `atan2` でなく、外積の符号を使う。

## geo_int.hpp 関数一覧

| 関数 | 使い方・意味 | 計算量 |
|---|---|---|
| `cross(a,b)`, `cross(a,b,c)` | 外積、`cross(b-a,c-a)` | O(1) |
| `dot(a,b)`, `dist2(a,b)` | 内積、距離の2乗（`i128`） | O(1) |
| `orient(a,b,c)` | 反時計回り +1 / 一直線 0 / 時計回り -1 | O(1) |
| `onSegment(a,b,p)` | `p` が端点を含む線分 `ab` 上か | O(1) |
| `segmentsIntersect(a,b,c,d)` | 閉線分同士が交わるか（接触含む） | O(1) |
| `convexHull(pts)` | 凸包。CCW、辺上の中間点を省く | O(n log n) |
| `area2(poly)` | 多角形の符号付き **2倍面積** (`i128`) | O(n) |
| `perimeter(poly)` | 周長 (`long double`) | O(n) |
| `pointInPolygon(poly,p)` | 単純多角形、0:外 / 1:境界 / 2:内 | O(n) |
| `pointInConvex(h,p)` | **CCW・余分な共線点なし** の凸多角形、同じ返り値 | O(log n) |
| `diameter2(h)` | 同上の凸多角形の直径の2乗 | O(n) |
| `minkowskiSum(a,b)` | 同上の凸多角形同士の Minkowski 和 | O(a+b) |
| `closestPair2(pts)` | 最近点対の距離の2乗 (`optional<i128>`) | O(n log n) |
| `toString(i128)`, `gcd128(a,b)` | 128bit 値の出力、約分用の最大公約数 | 桁数に比例 |

`minkowskiSum` は `convexHull` した結果を渡す。`a.size()==1` や `a.size()==2` も許容。

## geo_real.hpp 関数一覧

| 関数 | 用途 |
|---|---|
| `projection(line,p)`, `reflection(line,p)` | 射影、反射 |
| `lineIntersection(l1,l2)` | 一意な直線交点 (`optional<P>`、平行・一致は `nullopt`) |
| `onSegment(a,b,p)`, `segmentsIntersect(s,t)` | 線分上判定、線分交差 |
| `distance(line,p)`, `distance(segment,p)`, `distance(segment,segment)` | 距離 |
| `convexHull`, `signedArea`, `area`, `pointInPolygon`, `pointInConvex` | 基本多角形操作 |
| `convexCut(poly,line)` | CCW 凸多角形を有向直線の左半平面で切る |
| `convexIntersection(a,b)` | CCW 凸多角形同士の共通部分 |
| `circle2`, `circle3`, `minimumEnclosingCircle` | 2/3点通過円、最小包含円（期待値 O(n)） |
| `circleLineIntersection(c,l)` | 円・直線の交点 0/1/2 個 |
| `circleCircleIntersection(a,b)` | 円・円の交点 0/1/2 個。同一円は空扱い |
| `tangentPoints(c,p)` | 点から円への接点 0/1/2 個 |
| `commonTangents(a,b)` | 2円の共通接線（異なる2点による `Line`。接触時は接点＋方向点） |
| `circlePolygonIntersectionArea(poly,c)` | 円と単純多角形の共通面積 O(n) |

## コピペ例

```cpp
#include "geo_int.hpp"
using namespace std;
namespace gi = geo_int;
int main() {
    vector<gi::P> ps={{0,0},{2,0},{0,2},{1,1}};
    auto h=gi::convexHull(ps);
    auto twiceArea=gi::area2(h); // 4 (__int128)
    cout << gi::toString(twiceArea) << '\n';
    cout << gi::pointInConvex(h,{1,1}) << '\n'; // 1 (斜辺上)
}
```

有理数の面積は `area2` を分子として `2 * 拡大率^2` を分母に置き、`gcd` で既約化する。中間計算が `__int128` になるケースも、`geo_int::gcd128` と `geo_int::toString` でそのまま約分・出力できる。

## 注意事項

- 整数座標の差や和自体は `long long`。外積だけ128bitにしても、入力座標が64bit限界付近なら `p-q` の段階で溢れる。各問題の上界を確認する。
- `geo_int::closestPair2` は2乗距離が `2^125` 未満であることが必要（座標が通常の競プロ範囲なら十分）。
- 浮動小数点の凸包は誤差で破綻する条件がある。整数入力なら整数版を優先。
- `lineIntersection` で一致する2直線は一意な交点が存在しないので `nullopt`。
- `circleCircleIntersection` で一致する円も交点が無数にあるので空配列。
- `pointInConvex`, `diameter2`, `minkowskiSum` は `convexHull` の結果（CCW、共線な中間点なし）を渡す。
- 実数版の接線関数は、半径が正の円を前提に使う。
- **円・線分や円・多角形の交点で特別な退化条件** は特に EPS に注意。

- `commonTangents` の返却直線は常に異なる2点で表す。接触時は第2点を「接点＋接線方向」とし、円Bの接点そのものにはしない。
- 旧 `geometoric.txt` は廃止。独自機能の移植は今回の修正には含まない。
