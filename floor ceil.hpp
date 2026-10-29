#pragma once
// Signed ll division, rounded mathematically. Throws if divisor is zero or result cannot fit.
ll ceil(ll a, ll b) {
    if (b == 0) throw std::invalid_argument("division by zero");
    if (a == std::numeric_limits<ll>::min() && b == -1) throw std::overflow_error("division overflow");
    ll q = a / b, r = a % b;
    return q + (r != 0 && ((r > 0) == (b > 0)));
}
ll floor(ll a, ll b) {
    if (b == 0) throw std::invalid_argument("division by zero");
    if (a == std::numeric_limits<ll>::min() && b == -1) throw std::overflow_error("division overflow");
    ll q = a / b, r = a % b;
    return q - (r != 0 && ((r > 0) != (b > 0)));
}
