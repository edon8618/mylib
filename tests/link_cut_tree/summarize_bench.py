#!/usr/bin/env python3
import csv,statistics,pathlib
base=pathlib.Path(__file__).resolve().parent
labels=['頂点のみ・ランダム木','無向辺のみ','方向付き辺のみ','頂点＋方向付き辺','頂点のみ・長いパス','頂点のみ・link/cut']
parts=['# LinkCutTree 検証・性能報告\n',
'''対象の旧版: `c64718ccd2586e5391177c174c448c0067be566e`。

## 変更内容とAPI

添付の確定APIを実装。詳細と全シグネチャは `docs/LinkCutTree.txt`。
端点指定の辺更新は4引数のみ。ID指定は `(id,x)` と `(id,f,g)`。
`link_edge(u,v,x)`、端点用 `set_edge(u,v,x)`、flip等の追加APIは作っていない。
ただし `S=int` の3整数呼び出しはID指定と同じ型であり、C++から意図を識別できない。
これを拒否するにはIDの型またはID指定側の関数名の変更が必要。

通常linkは直接接続と辺メタデータのみを使い、初めてset_edgeする時だけ辺ノードを追加。
ID・基準方向・現在の根を維持したまま展開する。値付き辺は順逆の値と集約、
基準方向との差を示すフラグで管理。辺ノードを再利用し、辺IDは再利用しない。
Splay作業スタックを再利用し、kthも反復処理。try_prodはmake_root/accessを1組だけ行う。
辺ノード未展開の辺がある場合、kth_edgeは2つのkth_vertexと端点検索で求める。

## 互換性・計算量

旧split/get_prod/get_sz/get_kthを維持。サイズ・kthは内部辺ノードを除いた元の頂点基準。
旧linkの同成分時、旧cutの辺なし時のno-opも維持。根の具体的な選択や補助木の内部状態を
固定する保証はない。全操作は有効な頂点・辺ID等の前提条件の下で使用する。

経路操作は償却O(log N)、辺ID検索はハッシュの期待O(1)。空間はO(N+A)、Aは累積辺追加数。
辺IDの履歴が残るため、削除しても履歴領域は減らない。演算とSのコピーがO(1)の場合。

## テスト結果

- 整数加算（S=int）とmod998244353の非可換アフィン合成、各100seed×10000反復: 計200万反復PASS。
- BFSの愚直な森と、全要求API、両方向値、頂点・辺混在、切断と再接続、旧APIを比較。
- 係数0、逆順端点更新→ID更新、削除→同じ端点の逆向き再追加、普通のlinkへの後付け値を含む。
- 20万頂点のパスで繰り返し順逆照会・kth・切断再接続: PASS。
- 1万頂点の非可換・方向付きパスの順逆反復照会、N=0/1、展開時の根とID維持: PASS。
- ASan+UBSan: 20万ランダム反復＋上記重点ケースPASS。LeakSanitizerはptrace環境制約で無効化。
- 3引数端点オーバーロードの不在（S=long long）、デフォルト構築・等値比較を持たないS: コンパイル・実行PASS。
- docsの6使用例: C++17でコンパイル・実行PASS。
- 通貨変換例: ACL modintでコンパイル、逆順更新・非連結等の固定入出力が一致。
- 公式オンラインジャッジへの提出は行っていない。

## 計測条件と読み方

GNU g++ 14.2.0、`-std=c++20 -O2 -DNDEBUG`。CPU情報: Intel Xeon Platinum 8573C。
N=Q、固定seed。3回の独立プロセスの中央値。新旧の順序を交互に入れ替える。
構築とQ操作を別計測。25%更新＋75%経路照会、link/cut行は25%切断・再追加ペア。
RSSはプロセス全体の最大常駐メモリ（入力用配列やランタイムを含む）、MiB=1024 KiB。
比較可能な全行で新旧のチェックサム一致を確認。

旧版の頂点用途は元の実装そのまま。無向辺用途は利用者側でN-1個の辺ノードを確保。
旧版は方向付き辺に直接対応しないので、その比較欄は「—」。
共有環境で時間に大きなばらつきがあり、比率はこの計測の中央値の比で、一般的な速度保証ではない。
プロセスCPU時間と経過時間の両方を保存し、全生データと再現用スクリプトを同梱。
''']
for n in [200000,500000]:
    filename=base/f'benchmark_{n}.csv'
    if not filename.exists():continue
    rows=list(csv.DictReader(filename.open()))
    assert len(rows)==60,(filename,len(rows))
    for mono in ['sum','affine']:
        parts.append(f'\n## N=Q={n:,} / '+('整数和' if mono=='sum' else 'アフィン合成')+'\n')
        parts.append('|ケース|CPU秒 旧→新|新/旧|経過秒 旧→新|構築CPU秒 旧→新|最大RSS MiB 旧→新|\n|---|---:|---:|---:|---:|---:|')
        for mode in range(6):
            groups={a:[r for r in rows if r['implementation']==a and r['monoid']==mono and int(r['mode'])==mode] for a in ['old','new']}
            if groups['old']:
                assert {r['checksum'] for r in groups['old']}=={r['checksum'] for r in groups['new']}
            def median(a,key):return statistics.median(float(r[key]) for r in groups[a]) if groups[a] else None
            def pair(key,scale=1):return ' → '.join(f'{median(a,key)/scale:.3f}' if groups[a] else '—' for a in ['old','new'])
            ratio=f'{median("new","operations_cpu_s")/median("old","operations_cpu_s"):.2f}' if groups['old'] else '—'
            parts.append(f'|{labels[mode]}|{pair("operations_cpu_s")}|{ratio}|{pair("operations_s")}|{pair("build_cpu_s")}|{pair("peak_rss_kib",1024)}|')
    parts.append(f'\n生データ: `benchmark_{n}.csv`。\n')
parts.append('''
## 性能上の変更と制約

最初の案では普通のlinkにも辺ノードを作っていたが、頂点用途への負担を減らすため遅延展開へ変更。
最終版は従来の頂点用途で補助木のノード数を増やさない。Splayごとの一時vector確保も削除。
コンストラクタは初めに容量を確保してから頂点を構築し、初期配列の再確保・全コピーを避ける。
一方、辺ID用のハッシュ・端点履歴と、方向値・個数を保持するノード領域は追加で必要。
link/cutではハッシュ更新と履歴追加のコストが増える。全用途で旧版より高速になる保証はない。
未展開辺があるとkth_edgeは経路抽出を2回行うが、償却O(log N)は維持する。

通貨変換問題は正式な入力形式が未提供。`examples/currency_network.cpp` は先頭コメントに
明記した仮の入力形式で検証したAPI使用例であり、そのまま提出できるという保証はない。
GitHubへのpushは実施していない。
''')
(base/'REPORT.md').write_text('\n'.join(parts))
print(base/'REPORT.md')
