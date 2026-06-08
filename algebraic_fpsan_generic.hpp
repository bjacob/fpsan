// algebraic_fpsan_generic.hpp
//
// Copyright 2026 AMD
// SPDX-License-Identifier: Apache-2.0
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ----------------------------------------------------------------------------
// A prototype of the *algebraic* alternative to FPSan studied in
// algebraic-fpsan.md.  Where FPSan encodes a float by a scrambling bijection
// and works in the FREE ring Z/2^w, this encodes a float by the genuine ring
// homomorphism
//
//     phi_n : Z[1/2] --> Z/nZ ,   phi_n(m * 2^e) = (m mod n) * (2^-1)^(-e)
//
// so the residue carried is phi_n of the expression's EXACT (un-rounded) value.
// Two compact 32-bit variants (see algebraic-fpsan.md for the design space):
//
//   * single-prime (1b):  n = p, a prime just below 2^32.  Z/n is a FIELD:
//     division is total, x/x == 1 always.  No exp homomorphism (prime field).
//
//   * CRT (p*d):  n = p*d with p,d distinct odd primes and p == 1 (mod d),
//     n < 2^32.  Z/n ~ F_p x F_d is a product of fields (not a field: it has
//     zero-divisors).  It DOES support an exact exp homomorphism
//     exp(v) = g^(v mod d) with g of order d, so exp(a+b) == exp(a)*exp(b).
//
// Non-finite values use the projective line plus an absorbing NaN
//   Z/nZ  u  {Inf}  u  {NaN}                         (option C in inf-nan notes)
// with Inf the single (unsigned) pole 1/0 and NaN the absorbing bottom that
// totalizes the indeterminate forms (Inf+Inf, 0*Inf, Inf/Inf, 0/0).
//
// Order / min / max / abs / comparisons are NOT modeled: a finite field has no
// compatible order.  Overflow-to-Inf is NOT modeled either: the model is exact
// (un-rounded), so Inf arises only from division by a true zero.
//
// Requires C++20 and the __int128 extension (clang/gcc).  No dependencies
// beyond the standard library.
// ----------------------------------------------------------------------------
#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>

namespace algebraic_fpsan {

// ===========================================================================
//  Modular arithmetic mod n, with n < 2^32 (so products fit in uint64_t).
// ===========================================================================
inline uint64_t addmod(uint64_t a, uint64_t b, uint64_t n) {
  uint64_t s = a + b;
  return s >= n ? s - n : s;
}
inline uint64_t submod(uint64_t a, uint64_t b, uint64_t n) {
  return a >= b ? a - b : a + n - b;
}
inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t n) {
  return (a * b) % n; // a,b < n < 2^32  =>  a*b < 2^64
}
inline uint64_t powmod(uint64_t base, uint64_t e, uint64_t n) {
  uint64_t r = 1 % n;
  base %= n;
  while (e) {
    if (e & 1)
      r = mulmod(r, base, n);
    base = mulmod(base, base, n);
    e >>= 1;
  }
  return r;
}
// Extended Euclid inverse mod n (n need not be prime).  Returns false if a is
// not a unit (gcd(a,n) != 1) -- i.e. a zero-divisor.
inline bool invmod(uint64_t a, uint64_t n, uint64_t &out) {
  int64_t t = 0, newt = 1;
  int64_t r = (int64_t)n, newr = (int64_t)(a % n);
  while (newr != 0) {
    int64_t q = r / newr;
    int64_t tmp = t - q * newt;
    t = newt;
    newt = tmp;
    tmp = r - q * newr;
    r = newr;
    newr = tmp;
  }
  if (r != 1)
    return false; // not invertible
  if (t < 0)
    t += (int64_t)n;
  out = (uint64_t)t;
  return true;
}

// Trivial primality (n < 2^32 here; trial division to sqrt is plenty).
inline bool is_prime(uint64_t n) {
  if (n < 2)
    return false;
  if (n % 2 == 0)
    return n == 2;
  for (uint64_t i = 3; i * i <= n; i += 2)
    if (n % i == 0)
      return false;
  return true;
}

// ===========================================================================
//  Exact dyadic rationals  m * 2^e  (the elements of Z[1/2]).  Used as the
//  ground-truth value type: every finite float decodes to one of these, and
//  we can do exact arithmetic on them to test the homomorphism.
// ===========================================================================
struct Dyadic {
  __int128 m = 0; // signed significand
  int e = 0;      // power of two

  static Dyadic from_int(long long v) { return {(__int128)v, 0}; }

  Dyadic operator-() const { return {-m, e}; }

  Dyadic operator+(const Dyadic &o) const {
    if (m == 0)
      return o;
    if (o.m == 0)
      return *this;
    int lo = e < o.e ? e : o.e;
    __int128 a = m << (e - lo);   // shifts are by >= 0
    __int128 b = o.m << (o.e - lo);
    return {a + b, lo};
  }
  Dyadic operator-(const Dyadic &o) const { return *this + (-o); }
  Dyadic operator*(const Dyadic &o) const { return {m * o.m, e + o.e}; }
};

// ===========================================================================
//  Config: the modulus and (optionally) the exp generator.
//
//  Residues live in [0, n).  The on-the-wire encoding of a Value is a single
//  uint32_t `code`: finite residues are themselves, and two reserved codes
//  above n encode Inf and NaN.  n + 1 < 2^32 is required (true for both
//  variants), so the whole thing is a compact 32-bit word.
// ===========================================================================
struct Config {
  std::string name;
  uint64_t n = 0;    // modulus; residues in [0, n)
  uint64_t p = 0;    // field prime  (n == p for single-prime variant)
  uint64_t d = 0;    // exp exponent modulus (0 if no exp homomorphism)
  uint64_t inv2 = 0; // inverse of 2 mod n
  uint64_t g = 0;    // element of order d in (Z/n)^*  (0 if no exp)
  bool has_exp = false;

  uint32_t INF() const { return (uint32_t)n; }
  uint32_t NAN_() const { return (uint32_t)(n + 1); }
};

// Single-prime variant 1b: n = p, a prime just below 2^32.  Field; no exp.
inline Config make_single_prime(uint64_t p) {
  assert(is_prime(p) && (p & 1) && p < (uint64_t{1} << 32) - 1);
  Config c;
  c.name = "single-prime(p=" + std::to_string(p) + ")";
  c.n = c.p = p;
  c.d = 0;
  c.inv2 = (p + 1) / 2; // inverse of 2 mod odd p
  c.has_exp = false;
  return c;
}

// CRT variant: n = p*d, p == 1 (mod d), g of order d => exp homomorphism.
inline Config make_crt(uint64_t p, uint64_t d, uint64_t g) {
  assert(is_prime(p) && is_prime(d) && p != d);
  assert((p - 1) % d == 0 && "need d | p-1 for an order-d element");
  uint64_t n = p * d;
  assert(n + 1 < (uint64_t{1} << 32) && "n must fit a compact 32-bit code");
  assert((n & 1) && "n must be odd so 2 is invertible");
  Config c;
  c.name = "CRT(p=" + std::to_string(p) + ",d=" + std::to_string(d) + ")";
  c.n = n;
  c.p = p;
  c.d = d;
  c.inv2 = (n + 1) / 2;
  c.g = g;
  c.has_exp = true;
  // sanity: g must have order exactly d mod n.
  assert(powmod(g, d, n) == 1 && powmod(g, 1, n) != 1);
  return c;
}

// Build CRT parameters automatically: pick a Sophie-Germain-style pair
// (d prime, p = 2d+1 prime) near 2^15 so n = p*d ~ 2^31 < 2^32, then an
// order-d generator g.
inline Config make_crt_auto() {
  for (uint64_t d = 32749; d >= 3; --d) {
    if (!is_prime(d))
      continue;
    uint64_t p = 2 * d + 1;
    if (!is_prime(p))
      continue;
    uint64_t n = p * d;
    if (n + 1 >= (uint64_t{1} << 32))
      continue;
    // order-d element h in F_p^*: a^((p-1)/d) for some non-d-th-power a.
    uint64_t h = 0;
    for (uint64_t a = 2; a < p; ++a) {
      uint64_t cand = powmod(a, (p - 1) / d, p);
      if (cand != 1) {
        h = cand;
        break;
      }
    }
    assert(h && powmod(h, d, p) == 1);
    // CRT-lift to g == h (mod p), g == 1 (mod d): order = lcm(d,1) = d.
    uint64_t pInvModD;
    bool ok = invmod(p % d, d, pInvModD);
    (void)ok;
    assert(ok);
    uint64_t k = mulmod(submod(1 % d, h % d, d), pInvModD, d);
    uint64_t g = h + p * k; // < p*d = n
    return make_crt(p, d, g);
  }
  assert(false && "no CRT params found");
  return {};
}

// ===========================================================================
//  Value: a point of  Z/nZ u {Inf} u {NaN}, the algebraic-fpsan datum.
//
//  The real encoding is the single uint32_t `code_`; `cfg_` is ambient
//  configuration (shared, not part of the per-value footprint), carried here
//  only for ergonomic operators -- exactly as fpsan_generic.hpp carries the
//  FPFormat in each FPSanFloat.
// ===========================================================================
class Value {
public:
  enum Kind { Finite, Inf, NaN };

  Value(const Config &cfg, uint32_t code) : cfg_(&cfg), code_(code) {}

  static Value finite(const Config &cfg, uint64_t residue) {
    return Value(cfg, (uint32_t)(residue % cfg.n));
  }
  static Value inf(const Config &cfg) { return Value(cfg, cfg.INF()); }
  static Value nan(const Config &cfg) { return Value(cfg, cfg.NAN_()); }

  // phi_n of an exact dyadic value.
  static Value encode(const Config &cfg, const Dyadic &v) {
    if (v.m == 0)
      return finite(cfg, 0);
    __int128 mr = v.m % (__int128)cfg.n;
    if (mr < 0)
      mr += cfg.n;
    uint64_t r = (uint64_t)mr;
    uint64_t pw = v.e >= 0 ? powmod(2, (uint64_t)v.e, cfg.n)
                           : powmod(cfg.inv2, (uint64_t)(-v.e), cfg.n);
    return finite(cfg, mulmod(r, pw, cfg.n));
  }

  // phi_n of an actual IEEE double (exact -- a double *is* a dyadic rational).
  static Value encode(const Config &cfg, double x) {
    if (std::isnan(x))
      return nan(cfg);
    if (std::isinf(x))
      return inf(cfg);
    if (x == 0.0)
      return finite(cfg, 0);
    int exp;
    double frac = std::frexp(x, &exp);    // x = frac * 2^exp, frac in [0.5,1)
    __int128 m = (__int128)std::ldexp(frac, 53); // exact 53-bit integer
    return encode(cfg, Dyadic{m, exp - 53});
  }

  Kind kind() const {
    if (code_ == cfg_->INF())
      return Inf;
    if (code_ == cfg_->NAN_())
      return NaN;
    return Finite;
  }
  bool is_finite() const { return kind() == Finite; }
  uint32_t code() const { return code_; }
  uint64_t residue() const { return code_; } // valid when finite
  const Config &config() const { return *cfg_; }

  Value operator-() const {
    if (kind() == Finite)
      return finite(*cfg_, submod(0, code_, cfg_->n));
    return *this; // -Inf == Inf (unsigned), -NaN == NaN
  }

  Value operator+(const Value &o) const {
    if (anyNaN(o))
      return nan(*cfg_);
    Kind a = kind(), b = o.kind();
    if (a == Inf || b == Inf) {
      if (a == Inf && b == Inf)
        return nan(*cfg_); // Inf + Inf : indeterminate (no signs)  -> NaN
      return inf(*cfg_);   // Inf + finite
    }
    return finite(*cfg_, addmod(code_, o.code_, cfg_->n));
  }
  Value operator-(const Value &o) const { return *this + (-o); }

  Value operator*(const Value &o) const {
    if (anyNaN(o))
      return nan(*cfg_);
    Kind a = kind(), b = o.kind();
    bool aZero = (a == Finite && code_ == 0);
    bool bZero = (b == Finite && o.code_ == 0);
    if (a == Inf || b == Inf) {
      if (aZero || bZero)
        return nan(*cfg_); // 0 * Inf -> NaN
      return inf(*cfg_);   // Inf * (nonzero / Inf)
    }
    return finite(*cfg_, mulmod(code_, o.code_, cfg_->n));
  }

  Value operator/(const Value &o) const {
    if (anyNaN(o))
      return nan(*cfg_);
    Kind a = kind(), b = o.kind();
    if (a == Inf && b == Inf)
      return nan(*cfg_); // Inf / Inf
    if (a == Inf)
      return inf(*cfg_); // Inf / finite
    if (b == Inf)
      return finite(*cfg_, 0); // finite / Inf -> 0
    // both finite
    if (o.code_ == 0)
      return code_ == 0 ? nan(*cfg_) : inf(*cfg_); // 0/0 -> NaN ; x/0 -> Inf
    uint64_t inv;
    if (!invmod(o.code_, cfg_->n, inv))
      return nan(*cfg_); // divisor is a zero-divisor (CRT variant) -> poison
    return finite(*cfg_, mulmod(code_, inv, cfg_->n));
  }

  // exp.  CRT variant: the genuine homomorphism g^(v mod d).  Single-prime
  // variant: a deterministic hash token (no exp homomorphism is possible in a
  // prime field), tagged so different transcendentals don't collide.
  Value exp() const {
    if (kind() == NaN)
      return nan(*cfg_);
    if (kind() == Inf)
      return nan(*cfg_); // e^(+inf)=inf, e^(-inf)=0; unsigned -> indeterminate
    if (cfg_->has_exp)
      return finite(*cfg_, powmod(cfg_->g, code_ % cfg_->d, cfg_->n));
    return transcendental(/*tag=*/0x6578'70ULL); // "exp"
  }

  // Structure-less transcendentals: a deterministic, op-tagged hash token
  // (a fresh free generator), per algebraic-fpsan.md section 8.
  Value transcendental(uint64_t tag) const {
    if (kind() != Finite)
      return nan(*cfg_);
    uint64_t h = splitmix(tag ^ (0x9E37'79B9'7F4A'7C15ULL * (code_ + 1)));
    return finite(*cfg_, h % cfg_->n);
  }
  Value log() const { return transcendental(0x6C6F'67ULL); }
  Value sqrt() const { return transcendental(0x7371'7274ULL); }
  Value sin() const { return transcendental(0x73'696EULL); }
  Value cos() const { return transcendental(0x636F'73ULL); }

  // Fingerprint equality (NaN == NaN here: we compare model fingerprints, not
  // IEEE values -- two computations that both blow up should match).
  bool operator==(const Value &o) const {
    return cfg_ == o.cfg_ && code_ == o.code_;
  }
  bool operator!=(const Value &o) const { return !(*this == o); }

  std::string str() const {
    switch (kind()) {
    case Inf:
      return "Inf";
    case NaN:
      return "NaN";
    default:
      return "<" + std::to_string(code_) + ">";
    }
  }

private:
  bool anyNaN(const Value &o) const {
    return kind() == NaN || o.kind() == NaN;
  }
  static uint64_t splitmix(uint64_t z) {
    z += 0x9E37'79B9'7F4A'7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58'476D'1CE4'E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D0'49BB'1331'11EBULL;
    return z ^ (z >> 31);
  }

  const Config *cfg_;
  uint32_t code_;
};

} // namespace algebraic_fpsan
