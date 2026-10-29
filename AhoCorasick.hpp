#pragma once
//https://atcoder.jp/contests/abc419/submissions/me
// char_size: 文字の種類数 (例: 英小文字なら26, 数字なら10)
// base: 基準となる文字 (例: 'a', 'A', '0')
// stepは失敗辺も含めてうまくたどってくれる
// AhoCorasickとして使うならbuild()をちゃんと呼ぶ
// 逆にbuild()を呼ばないことでtrieとして使える
template <int char_size = 26, char base = 'a'>
struct AhoCorasick {
    struct Node {
        vector<int> next;   // 遷移先ノード
        int fail;           // 失敗辺（次にマッチする可能性のある最長接尾辞のノード）
        vector<int> accept; // このノードで完全にマッチする文字列のID
        
        Node() : next(char_size, -1), fail(-1) {}
    };

    vector<Node> nodes;

    AhoCorasick() {
        nodes.emplace_back(); // 初期状態（根）となるノード0を追加
    }

    // 文字列を辞書に追加する
    void insert(const string& s, int word_id) {
        int now = 0;
        for (char c : s) {
            int v = c - base;
            if (nodes[now].next[v] == -1) {
                nodes[now].next[v] = nodes.size();
                nodes.emplace_back();
            }
            now = nodes[now].next[v];
        }
        nodes[now].accept.push_back(word_id);
    }

    // すべての文字列を追加し終わった後に、失敗辺(fail)と遷移先(next)を構築する
    void build() {
        queue<int> que;
        
        // 根(0)の直接の子の処理
        for (int i = 0; i < char_size; i++) {
            if (nodes[0].next[i] != -1) {
                nodes[nodes[0].next[i]].fail = 0;
                que.push(nodes[0].next[i]);
            } else {
                nodes[0].next[i] = 0; // 根から遷移できない文字は根自身へ
            }
        }

        // BFSで失敗辺を配っていく
        while (!que.empty()) {
            int now = que.front();
            que.pop();

            for (int i = 0; i < char_size; i++) {
                if (nodes[now].next[i] != -1) {
                    int next_node = nodes[now].next[i];
                    
                    // 子の失敗辺は、「親の失敗辺の先にある、同じ文字の遷移先」
                    nodes[next_node].fail = nodes[nodes[now].fail].next[i];
                    
                    // 失敗辺の先にある受理状態（マッチング情報）も引き継ぐ
                    auto& fail_accept = nodes[nodes[nodes[now].fail].next[i]].accept;
                    nodes[next_node].accept.insert(
                        nodes[next_node].accept.end(),
                        fail_accept.begin(),
                        fail_accept.end()
                    );
                    
                    que.push(next_node);
                } else {
                    // 遷移先がない場合は、失敗辺の先の遷移先をショートカットとして登録 (オートマトン化)
                    nodes[now].next[i] = nodes[nodes[now].fail].next[i];
                }
            }
        }
    }
    
    // 現在のノード(状態)から文字 c で進んだ次のノード(状態)を返す (O(1))
    int step(int now, char c) const {
        return nodes[now].next[c - base];
    }

    // 現在のノード(状態)におけるマッチング数を返すなどの便利関数
    int count_matches(int now) const {
        return nodes[now].accept.size();
    }
};
