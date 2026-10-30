#pragma once
#include <queue>
#include <vector>
#include <utility>

// Prerequisite: qandori fps/fps.hpp (fps, mint).
// Takes ownership of the copy; pass std::move(polys) to avoid copying.
inline fps product_of_polynomials(std::vector<fps> polys) {
    if (polys.empty()) return fps{1};
    for (const auto& f : polys) if (f.empty()) return {};
    auto cmp = [&](int a, int b) { return polys[a].size() > polys[b].size(); };
    std::priority_queue<int, std::vector<int>, decltype(cmp)> q(cmp);
    for (int i = 0; i < (int)polys.size(); ++i) q.push(i);
    while (q.size() > 1) {
        int a = q.top(); q.pop();
        int b = q.top(); q.pop();
        polys[a] *= polys[b];
        fps().swap(polys[b]);
        q.push(a);
    }
    return std::move(polys[q.top()]);
}
