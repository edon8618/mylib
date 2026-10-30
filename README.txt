qandori FPS用追加ライブラリ（3本）

product_of_polynomials.hpp
stirling_number.hpp
sparse_pow_mul.hpp
各使用説明は docs/ の同名txt。

qandoriのソースは同梱していません。
https://github.com/qandori/cp-library
検証版: e2a1b76404b44706026ab3fdf739c72d7fd77f5a
各ヘッダーを読み込む前にqandoriのヘッダーをincludeしてください。
qandoriのmintと、tmp.cppの using mint=atcoder::modint998244353 は共存不可。
qandoriを使うソースではmintの別名定義を重ねないでください。

tests/test.cpp: 乱数seed固定の愚直比較・DP比較・境界条件検証。
実行:
g++ -std=c++20 -O2 -I/path/to/cp-library -I. tests/test.cpp -o test
./test

tests/example_*.cpp: 各txtに掲載した独立した使用例。
ソースコードを変更せず同じファイルを複数回includeしてもよいよう、
#pragma once と inline を使用しています。

以前の実装との意図的な差:
・product_of_polynomialsは通常呼出しで入力配列を変更しない。
・stirling_numberはqandoriのpow・taylor_shiftに合わせた。
・sparse_pow_mulは0^0=1、負指数の底は定数項非零。
  ローラン級数の相殺は扱わない。不正入力で零を返す挙動は廃止。
・LLONG_MINの符号反転オーバーフローを回避。
