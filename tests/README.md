# 検証手順

mylibディレクトリで実行。GCC/ClangとC++17、ACL、Boostヘッダを用意する。
実行可能ファイルは本体ディレクトリ外へ出す。

```bash
ACLPATH=/path/to/ac-library
BOOSTPATH=/path/to/boost/include
OUTDIR=$(mktemp -d)
g++ -std=c++17 -O2 -D_GLIBCXX_ASSERTIONS tests/persistent_test.cpp -o "$OUTDIR/persistent"
"$OUTDIR/persistent"
g++ -std=c++17 -O2 -D_GLIBCXX_ASSERTIONS -I "$BOOSTPATH" tests/revision_test.cpp -o "$OUTDIR/revision"
"$OUTDIR/revision"
g++ -std=c++17 -O2 -D_GLIBCXX_ASSERTIONS -I "$ACLPATH" tests/fps_test.cpp -o "$OUTDIR/fps"
"$OUTDIR/fps"
"$OUTDIR/fps" negative
g++ -std=c++17 -O2 -D_GLIBCXX_ASSERTIONS -I "$ACLPATH" tests/fps_extended_test.cpp -o "$OUTDIR/fps_extended"
"$OUTDIR/fps_extended"
g++ -std=c++17 -O2 -DNDEBUG -I "$ACLPATH" tests/fps_test.cpp -o "$OUTDIR/release"
"$OUTDIR/release" negative
```

sanitizerは `-O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all` を付けて再コンパイルする。この実行環境ではLeakSanitizerが利用できないため、`ASAN_OPTIONS=detect_leaks=0` で実行した。

```bash
g++ -std=c++17 -O2 -DBINOMIAL_BENCH_COMPARE tests/binomial_bench.cpp -o "$OUTDIR/bench"
"$OUTDIR/bench"
g++ -std=c++17 -O2 tests/modmul_bench.cpp -o "$OUTDIR/modmul"
"$OUTDIR/modmul" 998244353
```

`revision_test.cpp` はBoostがある場合cpp_intのテストも実行する。Boostがない場合はその部分を省略して、binomialはmodで計算したPascal三角形と照合する。今回の検証ではBoostありで全分岐を実行。

通常のcorrectnessテストはassertを使うためNDEBUGを付けない。NDEBUGのテストは負pow拒否の専用経路を実行する。
