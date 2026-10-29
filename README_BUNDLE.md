# oj-bundle で使う

ライブラリの .txt 43個を同名の .hpp に変更し、先頭に #pragma once を追加しました。既存の .hpp は変更していません。tests/VERIFIED_RESULTS.txt は検証記録なのでそのままです。

## a.cpp の例

元のライブラリと同様、標準ヘッダー、using namespace std、ll、all などは呼び出し側で先に定義してください。ACL が必要な場合は、展開対象になるダブルクォートで読み込みます。

```cpp
#include <bits/stdc++.h>
#include "atcoder/all"
using namespace std;
using namespace atcoder;
using ll = long long;
using mint = modint998244353;
#define all(x) (x).begin(), (x).end()

#include "mylib/WaveletMatrix.hpp"
#include "mylib/range_set.hpp"

int main() {
    range_set<int> rs;
}
```

たとえば ~/lib/mylib/WaveletMatrix.hpp に配置した場合、a.cpp のあるディレクトリで次のように実行します。

```bash
oj-bundle a.cpp -I "$HOME/lib" -I "$HOME/lib/ac-library-master" > a_submit.cpp
g++ -std=c++20 -O2 a_submit.cpp -o a.out
```

oj-bundle はファイルを1本に展開するコマンドです。その後 g++ でコンパイルします。a_submit.cpp を提出してください。

ファイル名に空白や日本語を含むものも、同じ名前でダブルクォート内に指定できます。

```cpp
#include "mylib/Weighted Union Find.hpp"
#include "mylib/LCA ダブリング.hpp"
```

## 元のZIPに含まれていないファイル

tests/ 内の一部テストは FPS/ と Tree/ のファイルを参照しますが、元のZIPにこれらのフォルダーはありません。それらのテストは不足分を用意するまで実行できません。データ構造同士の名前の衝突など、元からある制約もそのままです。

Zone.Identifier の付随ファイルと .git の管理データはこの配布用ZIPに含めていません。
