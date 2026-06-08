// algebraic_demo.cpp
//
// Copyright 2026 AMD
// SPDX-License-Identifier: Apache-2.0
//
// Licensed under the Apache License, Version 2.0 (the "License").  See
// algebraic_fpsan_generic.hpp for the full notice.
// ----------------------------------------------------------------------------
// Exercises and stress-tests the algebraic-fpsan prototype on both compact
// variants:
//   * single-prime (1b):  n = 2^32 - 5  (a field; no exp homomorphism)
//   * CRT (p*d):           n ~ 2^31      (exp homomorphism; zero-divisors)
//
// Asserts the must-hold properties (ring homomorphism, ring laws, the Inf/NaN
// table, and -- CRT only -- the exp homomorphism), and reports the statistical
// ones (collision rates, exp image size, zero-divisor rate).  The tests are
// deliberately heavy: they double as a fuzz/stress harness for the prototype.
//
//   c++ -std=c++20 -O2 algebraic_demo.cpp -o algebraic_demo && ./algebraic_demo
// ----------------------------------------------------------------------------
#include "algebraic_fpsan_generic.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <unordered_set>
#include <vector>

using namespace algebraic_fpsan;

// ---- tiny test harness -----------------------------------------------------
static long g_pass = 0, g_fail = 0;
static void check(bool ok, const char *msg) {
  if (ok)
    ++g_pass;
  else {
    ++g_fail;
    std::printf("    FAIL: %s\n", msg);
  }
}

// ---- deterministic RNG -----------------------------------------------------
struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed) {}
  uint64_t next() {
    uint64_t z = (s += 0x9E37'79B9'7F4A'7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58'476D'1CE4'E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D0'49BB'1331'11EBULL;
    return z ^ (z >> 31);
  }
};

static Dyadic rand_dyadic(Rng &rng) {
  // m in [-2^39, 2^39),  e in [-20, 20]: products/sums stay well inside int128.
  int64_t m = (int64_t)(rng.next() & ((uint64_t{1} << 40) - 1)) -
              (int64_t)(uint64_t{1} << 39);
  int e = (int)(rng.next() % 41) - 20;
  return Dyadic{(__int128)m, e};
}

// ===========================================================================
static void test_encoding_basics(const Config &c) {
  std::printf("  encoding basics:\n");
  check(Value::encode(c, Dyadic::from_int(0)).code() == 0, "phi(0)=0");
  check(Value::encode(c, Dyadic::from_int(1)).code() == 1, "phi(1)=1");
  check(Value::encode(c, Dyadic::from_int(2)).code() == 2 % c.n, "phi(2)=2");
  check(Value::encode(c, Dyadic::from_int(-1)).code() == c.n - 1, "phi(-1)=-1");
  check(Value::encode(c, Dyadic{1, -1}).code() == c.inv2, "phi(1/2)=inv2");
  // 0.1 etc. are exact dyadics in double; just ensure they encode finite.
  check(Value::encode(c, 0.1).is_finite(), "phi(0.1) finite");
}

// The headline property: phi is a ring homomorphism on exact dyadics.
static void test_homomorphism(const Config &c, int N) {
  Rng rng(0xC0FFEEULL);
  long add_ok = 0, sub_ok = 0, mul_ok = 0;
  for (int i = 0; i < N; ++i) {
    Dyadic a = rand_dyadic(rng), b = rand_dyadic(rng);
    Value va = Value::encode(c, a), vb = Value::encode(c, b);
    add_ok += (Value::encode(c, a + b) == va + vb);
    sub_ok += (Value::encode(c, a - b) == va - vb);
    mul_ok += (Value::encode(c, a * b) == va * vb);
  }
  check(add_ok == N, "phi(a+b) == phi(a)+phi(b)  for all sampled pairs");
  check(sub_ok == N, "phi(a-b) == phi(a)-phi(b)  for all sampled pairs");
  check(mul_ok == N, "phi(a*b) == phi(a)*phi(b)  for all sampled pairs");
  std::printf("  homomorphism over %d random dyadic pairs: +%ld -%ld *%ld\n", N,
              add_ok, sub_ok, mul_ok);
}

// Ring laws + value-faithfulness (constant folding etc.).
static void test_ring_laws(const Config &c, int N) {
  Rng rng(0xBEEF1ULL);
  long assoc = 0, distrib = 0, comm = 0;
  for (int i = 0; i < N; ++i) {
    Value a = Value::finite(c, rng.next() % c.n);
    Value b = Value::finite(c, rng.next() % c.n);
    Value d = Value::finite(c, rng.next() % c.n);
    assoc += ((a + b) + d == a + (b + d));
    distrib += (a * (b + d) == a * b + a * d);
    comm += (a * b == b * a);
  }
  check(assoc == N, "associativity");
  check(distrib == N, "distributivity");
  check(comm == N, "commutativity");
  // value-faithful coincidences the FREE model would miss:
  auto E = [&](long long k) { return Value::encode(c, Dyadic::from_int(k)); };
  check(E(2) + E(2) == E(4), "2+2 == 4 (constant folding)");
  check(E(3) * E(3) == E(9), "3*3 == 9");
  Value x = Value::encode(c, Dyadic{12345, -3});
  check(x + x == E(2) * x, "x+x == 2*x");
}

// Field / division behaviour.
static void test_division(const Config &c, int N) {
  Rng rng(0xD1D1DEULL);
  long unit_id = 0, units = 0, zerodiv = 0;
  for (int i = 0; i < N; ++i) {
    uint64_t r = 1 + rng.next() % (c.n - 1); // nonzero residue
    Value x = Value::finite(c, r);
    Value q = x / x;
    uint64_t dummy;
    if (invmod(r, c.n, dummy)) {
      ++units;
      unit_id += (q == Value::encode(c, Dyadic::from_int(1)));
    } else {
      ++zerodiv;
      check(q.kind() == Value::NaN, "x/x for a zero-divisor poisons to NaN");
    }
  }
  check(unit_id == units, "x/x == 1 for every unit");
  std::printf("  division: %ld units (x/x==1 all), %ld zero-divisors (NaN)\n",
              units, zerodiv);
}

// The projective-line + absorbing-NaN table (option C).
static void test_inf_nan(const Config &c) {
  std::printf("  Inf / NaN table:\n");
  Value one = Value::encode(c, Dyadic::from_int(1));
  Value two = Value::encode(c, Dyadic::from_int(2));
  Value zero = Value::finite(c, 0);
  Value INF = Value::inf(c);
  Value NAN_ = Value::nan(c);
  check((one / zero).kind() == Value::Inf, "1/0 == Inf");
  check((one / INF) == zero, "1/Inf == 0");
  check((INF + one).kind() == Value::Inf, "Inf + finite == Inf");
  check((INF + INF).kind() == Value::NaN, "Inf + Inf == NaN");
  check((INF - INF).kind() == Value::NaN, "Inf - Inf == NaN");
  check((INF * two).kind() == Value::Inf, "Inf * 2 == Inf");
  check((zero * INF).kind() == Value::NaN, "0 * Inf == NaN");
  check((INF * INF).kind() == Value::Inf, "Inf * Inf == Inf");
  check((INF / INF).kind() == Value::NaN, "Inf / Inf == NaN");
  check((zero / zero).kind() == Value::NaN, "0/0 == NaN");
  check((NAN_ + one).kind() == Value::NaN, "NaN absorbs +");
  check((NAN_ * INF).kind() == Value::NaN, "NaN absorbs *");
  check((INF / two).kind() == Value::Inf, "Inf / 2 == Inf");
}

// exp: must be a homomorphism on the CRT variant; cannot be on the prime field.
static void test_exp(const Config &c, int N) {
  Rng rng(0xE7E7ULL);
  if (c.has_exp)
    check(Value::finite(c, 0).exp() == Value::encode(c, Dyadic::from_int(1)),
          "exp(0) == 1");
  long hom = 0;
  for (int i = 0; i < N; ++i) {
    Value a = Value::finite(c, rng.next() % c.n);
    Value b = Value::finite(c, rng.next() % c.n);
    hom += (a + b).exp() == a.exp() * b.exp();
  }
  if (c.has_exp) {
    check(hom == N, "exp(a+b) == exp(a)*exp(b)  (CRT: exact homomorphism)");
    std::printf("  exp homomorphism holds %ld/%d (exact)\n", hom, N);
  } else {
    std::printf("  exp homomorphism holds %ld/%d  (prime field: ~0 expected, "
                "no exp)\n",
                hom, N);
  }
}

// ---- statistical stress: collisions, exp image size, zero-divisor rate -----
static std::vector<float> distinct_floats(int N) {
  Rng rng(0x5EED5ULL);
  std::unordered_set<uint32_t> seen;
  std::vector<float> out;
  while ((int)out.size() < N) {
    uint32_t bits = (uint32_t)rng.next();
    float f;
    std::memcpy(&f, &bits, 4);
    if (!std::isfinite(f))
      continue;
    if (seen.insert(bits).second)
      out.push_back(f);
  }
  return out;
}

static void test_collisions(const Config &c) {
  const int N = 300000;
  std::vector<float> fs = distinct_floats(N);

  // (a) leaf collisions: distinct floats sharing a residue.
  std::vector<uint64_t> res;
  res.reserve(N);
  long zerodiv = 0;
  for (float f : fs) {
    Value v = Value::encode(c, (double)f);
    res.push_back(v.code());
    uint64_t dummy;
    if (v.code() != 0 && !invmod(v.code(), c.n, dummy))
      ++zerodiv;
  }
  std::sort(res.begin(), res.end());
  long distinct = std::unique(res.begin(), res.end()) - res.begin();
  double expected = (double)N * (N - 1) / 2.0 / (double)c.n;
  std::printf("  leaf collisions: %ld inputs, %ld collided  (expected ~%.1f, "
              "n=%llu)\n",
              (long)N, (long)(N - distinct), expected,
              (unsigned long long)c.n);

  // (b) zero-divisor rate among nonzero floats (CRT only sees these).
  double zexp = (double)N * (1.0 / (double)c.p + (c.d ? 1.0 / (double)c.d : 0));
  std::printf("  zero-divisor inputs: %ld / %d  (expected ~%.1f)\n", zerodiv, N,
              zexp);

  // (c) exp image size: exp ranges over at most d values (CRT) -> tiny image.
  if (c.has_exp) {
    std::unordered_set<uint32_t> img;
    Rng rng(0x1234ULL);
    int M = 300000;
    for (int i = 0; i < M; ++i)
      img.insert(Value::finite(c, rng.next() % c.n).exp().code());
    std::printf("  exp image: %zu distinct outputs from %d inputs  (bounded by "
                "d=%llu)\n",
                img.size(), M, (unsigned long long)c.d);
  }
}

// ===========================================================================
static void run(const Config &c) {
  std::printf("\n================ %s  (n=%llu)%s ================\n",
              c.name.c_str(), (unsigned long long)c.n,
              c.has_exp ? "  [exp homomorphism]" : "  [field, no exp]");
  test_encoding_basics(c);
  test_homomorphism(c, 100000);
  test_ring_laws(c, 100000);
  test_division(c, 100000);
  test_inf_nan(c);
  test_exp(c, 50000);
  test_collisions(c);
}

int main() {
  Config single = make_single_prime(4294967291ULL); // 2^32 - 5, prime
  Config crt = make_crt_auto();

  std::printf("variant 1b : %s\n", single.name.c_str());
  std::printf("variant CRT: %s  (p=%llu, d=%llu, g=%llu)\n", crt.name.c_str(),
              (unsigned long long)crt.p, (unsigned long long)crt.d,
              (unsigned long long)crt.g);

  run(single);
  run(crt);

  std::printf("\n---------------------------------------------\n");
  std::printf("checks passed: %ld   failed: %ld\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}
