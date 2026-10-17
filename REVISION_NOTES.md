# mylib 修正記録 — 2026-10-15

ベース：直前の `mylib_revision_20261015.zip`。元の `mylib(2).zip` とも比較した。
以前の修正版で行われたB01・FPSの次数API変更・永続lazyセグ木削除は維持。今回の指定を追加し、変更箇所を実行検証した。

## FPSの次数API

`inv/log/exp/pow/sqrt/composition/compositional_inverse` と `dense_*/sparse_*` の引数 `deg` は **保持する最高次数**。`deg=D` の結果は原則 `D+1` 係数。

```cpp
FPS f{1, 1};
auto a = f.inv(3);       // x^0 ... x^3 の4項
assert(a.size() == 4);
auto b = f.inv(3).pow(2, 3); // f(x)^(-2)、4項
```

- 省略時は入力の `size()` 項を保持（既定最高次数は `size()-1`）。`deg=-1` は省略扱いの番兵。
- `pre(sz)` は係数**個数**の指定のまま。`pre(4)` は4項。
- `sparse_pow_mul(..., D)` もD+1項。`LIF(f,g,n,D)` はk=0..DのD+1個の係数 `[x^n] f*g^k`。
- 他のFPSファイルからの呼び出しも次数指定に合わせた。`stirling_number`等のnやmax_n、標本点個数mなどは従来の数学的意味を維持。
- 既存コードで「n項」のつもりで `f.inv(n)` 等を呼んでいた箇所は `f.inv(n-1)` 等へ変更が必要。
- 内部の `_size` 実装は非公開で係数個数を扱う。

## B01〜B20の対応

| ID | 対応 |
|---|---|
| B01 | 直前版の修正を維持し再検証。永続segtreeのmax_right/min_leftを境界・非可換順序・過去versionと照合。 |
| B02 | BitVectorの64bit境界では余計なblockを読まずprefix sumだけ返す。 |
| B03 | WM/Trieの検索上限がbit範囲以上なら全件、下限以下なら0。完全一致の範囲外値も0。bit幅の静的検査とTrieの入力条件assertも追加。 |
| B04 | AVL countのIterator/optional不一致を修正。 |
| B05 | FPS powは負指数をinvalid_argumentで拒否。NDEBUGでも有効。負乗は逆元が存在する場合 `f.inv(D).pow(r,D)` を使用（rは非負）。 |
| B06 | 旧実装で係数数0のinvがres[0]へ書き込んだ問題を解消。新inv(0)は定数項1個。平方根の出力範囲に先頭項が現れない場合はゼロ列を返す。合成逆関数の次数0も0を返す。 |
| B07 | sparse_pow_mulで指数0の因子を最初に1へ正規化し、0^0=1。従来の専用関数の符号付き指数対応は維持。 |
| B08 | 空多項式の等比多点評価は、一般多点評価に合わせてm個の0を返す。 |
| B09 | `geometoric/geometoric.txt` を削除。旧版独自機能の移植は行っていない。 |
| B10 | 接する2円の共通接線は接点＋接線方向の別点でLineを構成。一般ケースでは従来どおり各円の接点を使う。 |
| B11 | modの素数冪がUINT32_MAX以下なら従来のBarrett、超える場合はunsigned __int128の乗算・剰余を使用。正mod／非負n_maxの検査も追加。 |
| B12 | CRTは入力Tを維持して計算。ll/128bit等の固定幅型は中間の加減乗算のoverflowを検出してoverflow_error。cpp_intは任意精度計算。extend_gcdファイルの追加貼付は不要になった。 |
| B13 | floor/ceilを商と余りの符号で計算。ゼロ除数とLLONG_MIN/-1は例外。 |
| B14 | LCA2版ともk<0またはk>depth[v]なら-1。 |
| B15 | 動的Li Chaoをvector poolで再実装。定義域は閉区間、線分追加は半開区間。整数比較は128bit、実数は表現可能な隣接座標まで分割。 |
| B16 | AVL/Treapのコピーを禁止し、move constructor/assignmentを実装。move元は空になり再利用できる。 |
| B17 | `Diophantine::count(...)` のll専用APIへ変更。内部は128bitで計算し、解数がllに入らなければ例外。`Diophantine<ll>` という旧呼び方は使わない。 |
| B18 | 保留。直前版から変更していない。短いfを自動ゼロ拡張する修正は未実施。 |
| B19 | TreeのLCA加算距離とXOR前処理を分離。実数重みでbuild_lca可能。XORには整数型と明示的build_xorが必要。LCA再構築時のtourもclear。 |
| B20 | pow_modの初期値を1%modに修正。mod=1, exponent=0なら0。 |

## 重要なAPI例

```cpp
// 加算距離とXOR距離を別々に準備する
Tree<long long> t(n);
// ... add_edge ...
t.build_lca(root);
t.build_xor(root);       // get_xor_distを使う場合だけ必要

// CRT: 全引数の型を揃える。型は呼び出し側が選ぶ
using I = __int128_t;
auto solution = CRT(I(r1), I(m1), I(r2), I(m2));
// T=llでoverflow_errorが出たら、必要に応じてIやcpp_intで再計算。

DynamicLiChaoTree<long long> lc(-10, 10, LLONG_MIN); // [-10,10]、右端を含む
lc.add_segment(2, 3, -5, 5);                       // [-5,5)、右端を含まない
```

CRTは自動で大きい型へ昇格しない。最終解がTに収まっても中間積が収まらない場合は例外になる。固定128bitを強制する方法や、全件cpp_int化する方法を採用していない。

動的Li Chaoの整数query結果がTに収まらなければoverflow_error。実数の座標・係数は有限値を要求し、評価の精度はTの算術に依存する。実数版は座標のbit構造によって再帰の深さが異なり、一律のEPS近似で打ち切らない。

## 実行検証

環境：x86_64、GCC 13.3.0、C++17。**実際のACLヘッダとBoost 1.83を使用した**。

| テスト | 内容／結果 |
|---|---|
| persistent_test | 長さ0..80、ランダム永続更新と全l/rの探索、和・文字列連結の非可換順序を照合。PASS。 |
| revision_test | WM長さ0..130（64境界含む）、WM/Trie検索上限、AVL5,000操作、AVL/Treap move、CRTのll/128bit/cpp_intとoverflow、整数除算極値、Diophantine愚直解、LCA2版、Tree<double>とXOR、binomial、接線、整数／実数Li Chaoを照合。PASS。 |
| fps_test | deg=0を含む基本演算、疎密経路、合成逆関数・合成、呼び出し元の長さ。PASS。 |
| fps_extended_test | 1,006ケース。0^0、空等比評価、sqrt出力範囲、負pow、LIF、一般評価、疎pow積、dense FPS（実際のNTT含む）、スターリング数を愚直解等と照合。PASS。 |
| NDEBUG実行 | 負powが例外で拒否されることを実行確認。PASS。 |
| ASan+UBSan | 上の4 correctness suiteを同じソースで実行。PASS。LeakSanitizerのみ実行環境の制約で無効化。 |

これは今回変更した機能に対する検証。変更対象外を含むライブラリ全機能や全入力条件の保証ではない。

## B11の計測

200,000回の二項係数クエリ、n_max=200,000、同じ乱数列、3回の中央値。前処理時間は含めていない。全32bit mod variantのchecksumは一致。

| 方式 | mod | 計測1 (ms) | 計測2 (ms) |
|---|---:|---:|---:|
| 旧Barrett | 998244353 | 59.718 | 56.505 |
| 128bit乗算を常用 | 998244353 | 59.327 | 55.372 |
| 今回の切替版 | 998244353 | 52.157 | 48.660 |
| 今回の切替版・大mod | 10460353203 (=3^21) | 183.589 | 160.501 |

大modは素数冪の計算段数が違うので、上の大mod行と小mod行の差を128bit化の遅延と解釈しない。環境差・コンパイラ最適化もあるため、切替版の数値を一般的な速度向上として保証しない。

mod=998244353の乗算だけ2,000,000回では、Barrett/128bit剰余は計測1が4.073/8.009ms、計測2が6.655/9.715ms。乗算単体は128bitが約1.5〜2倍遅い。そこで**小modの従来経路を残し、大modは128bitで正確に扱う**。32bit上限への制限は設けなかった。

二項係数の正しさはn<=150、mod=1,2,4,12,99991,998244353,3^21,10000000019をBoost cpp_intのPascal三角形と照合した。大きなmodの因数分解の速度や前計算上限は従来の制約が残る。

## ファイル構成

- 本体のパス／名前は削除対象を除いて維持。
- `CHANGES_FROM_AUDIT.diff`：添付mylib(2).zipとの差分。
- `CHANGES_FROM_PREVIOUS.diff`：直前の修正版との差分。
- `tests/`：実行テスト・ベンチ・再現手順。bench_referenceは速度比較専用で、本体ライブラリとして使わない。
- 差分はCRLF/LFを正規化して作成した。既存の関数雛形や今回対象外の改善候補は維持。
