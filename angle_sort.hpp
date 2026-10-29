#pragma once
// {0,0}は取り除いておかないと壊れる
struct Angle {
    ll x, y;

    // コンストラクタ
    Angle(ll x = 0, ll y = 0) : x(x), y(y) {}

    // ソート用の比較演算子 (<)
    // 偏角ソート (反時計回り: 0°〜360°) の順序を定義
    bool operator<(const Angle &another) const {
        // 場所（上半分か下半分か）を取得
        int place_self = location();
        int place_other = another.location();

        // 異なる領域なら、領域番号が小さい方が先
        if (place_self != place_other) {
            return place_self < place_other;
        }

        // 同じ領域なら、外積を使って角度を比較
        // オーバーフロー防止のために__int128を使うのが安全ですが、
        // 速度優先や制約内なら ll のままでも可
        return another.x * y < x * another.y;
    }

    // 等価演算子 (==)
    bool operator==(const Angle &another) const {
        return !(*this < another) && !(another < *this);
    }

private:
    // 旧 ok関数: 外部からは見えないように private に配置
    // 0: 上半分 (y>0) および 正のx軸 (y=0, x>0) -> [0, pi)
    // 1: 下半分 (y<0) および 負のx軸 (y=0, x<0) -> [pi, 2pi)
    int location() const {
        if (y < 0) return 1;
        if (y > 0) return 0;
        if (x > 0) return 0; // y == 0, x > 0
        if (x < 0) return 1; // y == 0, x < 0
        return 0;            // (0, 0)
    }
};