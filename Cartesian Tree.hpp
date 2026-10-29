#pragma once
/* https://judge.yosupo.jp/submission/347989
 * Cartesian Tree (O(N) Construction)
 * ----------------------------------------------------------------
 * 数列から「順序はBST、値はHeap」の性質を持つ木を構築する。
 *
 * [メンバ変数]
 * root: 根のインデックス
 * P[i]: 頂点 i の親 (なければ -1)
 * L[i]: 頂点 i の左の子 (なければ -1)
 * R[i]: 頂点 i の右の子 (なければ -1)
 *
 * [使い方]
 * 1. Min-Heap (デフォルト: 値が小さい方が親)
 * CartesianTree<long long> ct(A);
 *
 * 2. Max-Heap (値が大きい方が親)
 * CartesianTree<long long> ct(A, [](long long a, long long b) { return a > b; });
 *
 * 3. カスタム (例: popcountが多い方が親 / 同じなら値が小さい方が親)
 * auto comp = [](int a, int b) {
 *      int pa = __builtin_popcount(a);
 *      int pb = __builtin_popcount(b);
 *      if (pa != pb) return pa > pb; // popcountが多い方が強い
 *      return a < b;                 // 同じなら値が小さい方が強い
 * };
 * CartesianTree<int> ct(A, comp);
 */

template<typename T>
struct CartesianTree {
    int root;
    std::vector<int> P, L, R;

    template<typename Compare = std::less<T>>
    CartesianTree(const std::vector<T>& a, Compare comp = Compare()) {
        int n = a.size();
        P.assign(n, -1);
        L.assign(n, -1);
        R.assign(n, -1);

        std::stack<int> st;
        for (int i = 0; i < n; ++i) {
            int last_popped = -1;
            
            // 新しい要素 a[i] の方が、スタックトップより「親になる力が強い」場合
            // スタックトップを pop して、a[i] の左の子にする
            while (!st.empty() && comp(a[i], a[st.top()])) {
                last_popped = st.top();
                st.pop();
            }
            
            // pop された最後の要素（a[i]より弱かった中で一番右のもの）が左の子
            if (last_popped != -1) {
                P[last_popped] = i;
                L[i] = last_popped;
            }
            
            // スタックに残っているもの（a[i]より強いもの）がいれば、a[i] はその右の子
            if (!st.empty()) {
                P[i] = st.top();
                R[st.top()] = i;
            }
            
            st.push(i);
        }
        
        // 根を探す (親が -1 の場所)
        root = 0;
        if (n > 0) {
            while (P[root] != -1) root = P[root];
        }
    }
};
