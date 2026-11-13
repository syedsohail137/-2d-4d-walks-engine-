#ifndef BIG_RATIONAL_HPP
#define BIG_RATIONAL_HPP

#include <algorithm>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace spectral {

// Arbitrary-precision signed integer (base 10^9) to prevent denominator explosion
struct BigInt {
    static constexpr uint32_t BASE = 1000000000U;
    int sign = 0; // -1, 0, or 1
    std::vector<uint32_t> limbs; // little-endian

    BigInt(long long v = 0) {
        if (v == 0) {
            sign = 0;
        } else if (v < 0) {
            sign = -1;
            unsigned long long uv = static_cast<unsigned long long>(-(v + 1)) + 1ULL;
            while (uv > 0) {
                limbs.push_back(static_cast<uint32_t>(uv % BASE));
                uv /= BASE;
            }
        } else {
            sign = 1;
            unsigned long long uv = static_cast<unsigned long long>(v);
            while (uv > 0) {
                limbs.push_back(static_cast<uint32_t>(uv % BASE));
                uv /= BASE;
            }
        }
    }

    void trim() {
        while (!limbs.empty() && limbs.back() == 0) limbs.pop_back();
        if (limbs.empty()) sign = 0;
    }

    bool is_zero() const { return sign == 0; }
    bool is_one() const { return sign == 1 && limbs.size() == 1 && limbs[0] == 1; }

    BigInt abs() const {
        BigInt r = *this;
        if (r.sign < 0) r.sign = 1;
        return r;
    }

    BigInt operator-() const {
        BigInt r = *this;
        r.sign = -r.sign;
        return r;
    }

    static int cmp_abs(const BigInt& a, const BigInt& b) {
        if (a.limbs.size() != b.limbs.size())
            return (a.limbs.size() < b.limbs.size()) ? -1 : 1;
        for (std::size_t i = a.limbs.size(); i-- > 0;) {
            if (a.limbs[i] != b.limbs[i])
                return (a.limbs[i] < b.limbs[i]) ? -1 : 1;
        }
        return 0;
    }

    bool operator==(const BigInt& o) const { return sign == o.sign && limbs == o.limbs; }
    bool operator!=(const BigInt& o) const { return !(*this == o); }
    bool operator<(const BigInt& o) const {
        if (sign != o.sign) return sign < o.sign;
        if (sign == 0) return false;
        int c = cmp_abs(*this, o);
        return (sign > 0) ? (c < 0) : (c > 0);
    }

    static BigInt add_abs(const BigInt& a, const BigInt& b) {
        BigInt r;
        r.sign = 1;
        std::size_t n = std::max(a.limbs.size(), b.limbs.size());
        r.limbs.resize(n + 1, 0);
        uint64_t carry = 0;
        for (std::size_t i = 0; i < n; ++i) {
            uint64_t sum = carry;
            if (i < a.limbs.size()) sum += a.limbs[i];
            if (i < b.limbs.size()) sum += b.limbs[i];
            r.limbs[i] = static_cast<uint32_t>(sum % BASE);
            carry = sum / BASE;
        }
        r.limbs[n] = static_cast<uint32_t>(carry);
        r.trim();
        return r;
    }

    // Assumes |a| >= |b|
    static BigInt sub_abs(const BigInt& a, const BigInt& b) {
        BigInt r;
        r.sign = 1;
        r.limbs.resize(a.limbs.size(), 0);
        int64_t borrow = 0;
        for (std::size_t i = 0; i < a.limbs.size(); ++i) {
            int64_t diff = static_cast<int64_t>(a.limbs[i]) - borrow;
            if (i < b.limbs.size()) diff -= static_cast<int64_t>(b.limbs[i]);
            if (diff < 0) {
                diff += BASE;
                borrow = 1;
            } else {
                borrow = 0;
            }
            r.limbs[i] = static_cast<uint32_t>(diff);
        }
        r.trim();
        return r;
    }

    BigInt operator+(const BigInt& o) const {
        if (sign == 0) return o;
        if (o.sign == 0) return *this;
        if (sign == o.sign) {
            BigInt r = add_abs(*this, o);
            r.sign = sign;
            return r;
        }
        int c = cmp_abs(*this, o);
        if (c == 0) return BigInt(0);
        if (c > 0) {
            BigInt r = sub_abs(*this, o);
            r.sign = sign;
            return r;
        } else {
            BigInt r = sub_abs(o, *this);
            r.sign = o.sign;
            return r;
        }
    }

    BigInt operator-(const BigInt& o) const { return *this + (-o); }

    BigInt operator*(const BigInt& o) const {
        if (sign == 0 || o.sign == 0) return BigInt(0);
        BigInt r;
        r.sign = sign * o.sign;
        r.limbs.assign(limbs.size() + o.limbs.size(), 0);
        for (std::size_t i = 0; i < limbs.size(); ++i) {
            uint64_t carry = 0;
            for (std::size_t j = 0; j < o.limbs.size(); ++j) {
                uint64_t cur = r.limbs[i + j] +
                               static_cast<uint64_t>(limbs[i]) * o.limbs[j] + carry;
                r.limbs[i + j] = static_cast<uint32_t>(cur % BASE);
                carry = cur / BASE;
            }
            r.limbs[i + o.limbs.size()] += static_cast<uint32_t>(carry);
        }
        r.trim();
        return r;
    }

    static std::pair<BigInt, BigInt> divmod_abs(const BigInt& a, const BigInt& b) {
        if (b.is_zero()) throw std::domain_error("BigInt division by zero");
        BigInt ua = a.abs(), ub = b.abs();
        if (cmp_abs(ua, ub) < 0) return {BigInt(0), ua};

        BigInt q, r;
        q.sign = 1;
        q.limbs.assign(ua.limbs.size(), 0);

        for (std::size_t i = ua.limbs.size(); i-- > 0;) {
            r.limbs.insert(r.limbs.begin(), ua.limbs[i]);
            r.sign = 1;
            r.trim();

            uint32_t lo = 0, hi = BASE - 1, best = 0;
            while (lo <= hi) {
                uint32_t mid = lo + (hi - lo) / 2;
                BigInt prod = ub * BigInt(static_cast<long long>(mid));
                if (cmp_abs(prod, r) <= 0) {
                    best = mid;
                    lo = mid + 1;
                } else {
                    if (mid == 0) break;
                    hi = mid - 1;
                }
            }
            q.limbs[i] = best;
            if (best > 0) {
                r = sub_abs(r, ub * BigInt(static_cast<long long>(best)));
            }
        }
        q.trim();
        r.trim();
        return {q, r};
    }

    BigInt operator/(const BigInt& o) const {
        auto [q, r] = divmod_abs(*this, o);
        if (!q.is_zero()) q.sign = sign * o.sign;
        return q;
    }

    BigInt operator%(const BigInt& o) const {
        auto [q, r] = divmod_abs(*this, o);
        if (!r.is_zero()) r.sign = sign;
        return r;
    }

    static BigInt gcd(BigInt a, BigInt b) {
        a = a.abs();
        b = b.abs();
        while (!b.is_zero()) {
            BigInt t = a % b;
            a = std::move(b);
            b = std::move(t);
        }
        return a;
    }

    std::string str() const {
        if (sign == 0 || limbs.empty()) return "0";
        std::string s = (sign < 0 ? "-" : "") + std::to_string(limbs.back());
        for (std::size_t i = limbs.size() - 1; i-- > 0;) {
            std::string part = std::to_string(limbs[i]);
            s += std::string(9 - part.size(), '0') + part;
        }
        return s;
    }
};

struct BigRational {
    BigInt num;
    BigInt den;

    BigRational(BigInt n = BigInt(0), BigInt d = BigInt(1)) {
        if (d.is_zero()) {
            num = BigInt(0);
            den = BigInt(1);
            return;
        }
        if (d.sign < 0) {
            n = -n;
            d = -d;
        }
        if (n.is_zero()) {
            num = BigInt(0);
            den = BigInt(1);
            return;
        }
        BigInt g = BigInt::gcd(n, d);
        num = n / g;
        den = d / g;
    }

    BigRational operator+(const BigRational& o) const {
        return BigRational(num * o.den + o.num * den, den * o.den);
    }

    BigRational operator-(const BigRational& o) const {
        return BigRational(num * o.den - o.num * den, den * o.den);
    }

    BigRational operator*(const BigRational& o) const {
        return BigRational(num * o.num, den * o.den);
    }

    std::string str() const {
        if (den.is_one()) return num.str();
        return num.str() + "/" + den.str();
    }

    std::string compact_str() const {
        std::string full = str();
        if (full.size() <= 32) return full;
        return "[" + std::to_string(num.str().size()) + "d/" + std::to_string(den.str().size()) + "d]";
    }
};

inline std::ostream& operator<<(std::ostream& os, const BigRational& r) {
    return os << r.str();
}

} // namespace spectral

#endif // BIG_RATIONAL_HPP
