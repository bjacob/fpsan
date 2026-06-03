// demo.cpp
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
// Toy examples confirming, by hand, the properties the FPSan blog post claims
// for the embedding of float32 into FPSanF32.
//
// Build:   c++ -std=c++20 -O2 demo.cpp -o demo
// Run:     ./demo
// ----------------------------------------------------------------------------
#include "fpsan_f32.hpp"

#include <cstdint>
#include <cstdio>
#include <vector>

using fpsan::FPSanF32;
using fpsan::phi;

static int g_failures = 0;

static void check(const char *name, bool ok) {
  printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
  if (!ok)
    ++g_failures;
}

static const char *hex(uint32_t p) {
  static char buf[16];
  snprintf(buf, sizeof buf, "0x%08X", p);
  return buf;
}

int main() {
  // A handful of "interesting" floats to exercise the laws on.
  std::vector<float> xs = {0.0f,     1.0f,      -1.0f,   2.0f, 0.5f,
                           3.14159f, -2.71828f, 1e-3f,   1e6f, -42.0f,
                           0.1f,     123.456f,  -0.333f, 7.0f, -9.0f};

  // ------------------------------------------------------------------
  // 1. Stable fixed points of the embedding phi.
  // ------------------------------------------------------------------
  printf("== Fixed points of phi ==\n");
  printf("  phi(+0.0) = %s\n", hex(phi(0.0f).payload()));
  printf("  phi(+1.0) = %s\n", hex(phi(1.0f).payload()));
  printf("  phi(-1.0) = %s\n", hex(phi(-1.0f).payload()));
  check("phi(+0.0) == 0", phi(0.0f).payload() == 0x00000000u);
  check("phi(+1.0) == 1", phi(1.0f).payload() == 0x00000001u);
  check("phi(-1.0) == 0xFFFFFFFF", phi(-1.0f).payload() == 0xFFFFFFFFu);

  // ------------------------------------------------------------------
  // 2. phi is a bijection: phi^{-1}(phi(x)) recovers x exactly.
  // ------------------------------------------------------------------
  printf("\n== Round trip phi^{-1}(phi(x)) == x ==\n");
  {
    bool ok = true;
    for (float x : xs) {
      float back = phi(x).to_float();
      // exact bit-for-bit round trip
      if (fpsan::bits_of(back) != fpsan::bits_of(x)) {
        ok = false;
        printf("    mismatch: x=%g back=%g\n", x, back);
      }
    }
    check("round trips bit-exactly for all sample floats", ok);
  }

  // ------------------------------------------------------------------
  // 3. phi(-x) == -phi(x)  for nonzero x.
  // ------------------------------------------------------------------
  printf("\n== phi(-x) == -phi(x) (nonzero x) ==\n");
  {
    bool ok = true;
    for (float x : xs) {
      if (x == 0.0f)
        continue;
      uint32_t lhs = phi(-x).payload();
      uint32_t rhs = (0u - phi(x).payload()); // negation mod 2^32
      if (lhs != rhs) {
        ok = false;
        printf("    x=%g  phi(-x)=%s  -phi(x)=%s\n", x, hex(lhs), hex(rhs));
      }
    }
    check("phi(-x) == -phi(x) for all nonzero samples", ok);
  }

  // ------------------------------------------------------------------
  // 4. THE headline property: fast-math-style reassociation is EXACT.
  //    Over IEEE floats, a+(b+c) != (a+b)+c in general.  Over FPSan it
  //    is the same modular integer, so the payloads match exactly.
  // ------------------------------------------------------------------
  printf("\n== Associativity of + : (a+b)+c == a+(b+c) ==\n");
  {
    float a = 1e8f, b = -1e8f, c = 1.0f; // classic IEEE catastrophe
    float ieee_left = (a + b) + c;
    float ieee_right = a + (b + c);
    printf("  IEEE float:   (a+b)+c = %g, a+(b+c) = %g  -> %s\n", ieee_left,
           ieee_right, ieee_left == ieee_right ? "equal" : "DIFFER");

    FPSanF32 A = phi(a), B = phi(b), C = phi(c);
    FPSanF32 fs_left = (A + B) + C;
    FPSanF32 fs_right = A + (B + C);
    printf("  FPSan payload:(a+b)+c = %s, a+(b+c) = %s  -> %s\n",
           hex(fs_left.payload()), hex(fs_right.payload()),
           fs_left == fs_right ? "EQUAL" : "differ");
    check("FPSan + is exactly associative on the catastrophe example",
          fs_left == fs_right);
  }
  // ...and exhaustively over the sample set.
  {
    bool ok = true;
    for (float a : xs)
      for (float b : xs)
        for (float c : xs) {
          FPSanF32 A = phi(a), B = phi(b), C = phi(c);
          if ((A + B) + C != A + (B + C))
            ok = false;
        }
    check("FPSan + associative over all sample triples", ok);
  }

  // ------------------------------------------------------------------
  // 5. Commutativity, identities, ring laws.
  // ------------------------------------------------------------------
  printf("\n== Ring identities ==\n");
  {
    bool addComm = true, mulComm = true, mulAssoc = true, distrib = true;
    bool addId = true, mulId = true, negInv = true, subSelf = true;
    FPSanF32 ZERO = FPSanF32::zero(), ONE = FPSanF32::one();
    for (float a : xs) {
      FPSanF32 A = phi(a);
      if (A + ZERO != A)
        addId = false; // x + 0 = x
      if (A * ONE != A)
        mulId = false; // x * 1 = x
      if (A + (-A) != ZERO)
        negInv = false; // x + (-x) = 0
      if (A - A != ZERO)
        subSelf = false; // x - x = 0
      for (float b : xs) {
        FPSanF32 B = phi(b);
        if (A + B != B + A)
          addComm = false; // commutativity of +
        if (A * B != B * A)
          mulComm = false; // commutativity of *
        for (float c : xs) {
          FPSanF32 C = phi(c);
          if ((A * B) * C != A * (B * C))
            mulAssoc = false; // * assoc
          if (A * (B + C) != (A * B) + (A * C))
            distrib = false; // distrib
        }
      }
    }
    check("x + 0 == x", addId);
    check("x * 1 == x", mulId);
    check("x + (-x) == 0", negInv);
    check("x - x == 0", subSelf);
    check("a + b == b + a", addComm);
    check("a * b == b * a", mulComm);
    check("(a*b)*c == a*(b*c)", mulAssoc);
    check("a*(b+c) == a*b + a*c (distributivity)", distrib);
  }

  // ------------------------------------------------------------------
  // 6. Division laws.  x/1 == x, 1/(1/x) == x; for odd payloads y/y == 1.
  // ------------------------------------------------------------------
  printf("\n== Division ==\n");
  {
    bool divId = true, doubleRecip = true;
    FPSanF32 ONE = FPSanF32::one();
    for (float x : xs) {
      FPSanF32 X = phi(x);
      if (X / ONE != X)
        divId = false;
      if (ONE / (ONE / X) != X)
        doubleRecip = false;
    }
    check("x / 1 == x", divId);
    check("1 / (1 / x) == x", doubleRecip);
    // x/x == 1 only for odd payloads (documented caveat). Show one of each.
    for (float x : {2.0f, 3.0f}) {
      FPSanF32 X = phi(x);
      bool isOne = (X / X == FPSanF32::one());
      printf("  x=%g payload=%s (parity %s): x/x == 1 ? %s\n", x,
             hex(X.payload()), (X.payload() & 1) ? "odd" : "even",
             isOne ? "yes" : "no");
    }
  }

  // ------------------------------------------------------------------
  // 7. The exponential homomorphism: exp(x+y) == exp(x)*exp(y),
  //    exp(0) == 1, exp(-x) == 1/exp(x).
  // ------------------------------------------------------------------
  printf("\n== Exponential homomorphism ==\n");
  {
    bool hom = true, exp0 = (fpsan::exp(FPSanF32::zero()) == FPSanF32::one());
    bool expNeg = true;
    for (float a : xs)
      for (float b : xs) {
        FPSanF32 A = phi(a), B = phi(b);
        if (fpsan::exp(A + B) != fpsan::exp(A) * fpsan::exp(B))
          hom = false;
      }
    for (float a : xs) {
      FPSanF32 A = phi(a);
      if (fpsan::exp(-A) != FPSanF32::one() / fpsan::exp(A))
        expNeg = false;
    }
    check("exp(x + y) == exp(x) * exp(y)", hom);
    check("exp(0) == 1", exp0);
    check("exp(-x) == 1 / exp(x)", expNeg);

    // Same homomorphism for exp2.
    bool hom2 = true;
    for (float a : xs)
      for (float b : xs) {
        FPSanF32 A = phi(a), B = phi(b);
        if (fpsan::exp2(A + B) != fpsan::exp2(A) * fpsan::exp2(B))
          hom2 = false;
      }
    check("exp2(x + y) == exp2(x) * exp2(y)", hom2);
  }

  // ------------------------------------------------------------------
  // 8. Worked example: two algebraically-equal-but-IEEE-different kernels.
  //    f(x) = (x*x - 1)            vs   g(x) = (x-1)*(x+1)
  //    These are equal as polynomials, so FPSan payloads must match,
  //    even though IEEE rounding can make them differ.
  // ------------------------------------------------------------------
  printf("\n== Worked example: (x^2 - 1) vs (x-1)(x+1) ==\n");
  {
    bool ok = true;
    FPSanF32 ONE = FPSanF32::one();
    for (float x : xs) {
      FPSanF32 X = phi(x);
      FPSanF32 f = X * X - ONE;
      FPSanF32 g = (X - ONE) * (X + ONE);
      if (f != g) {
        ok = false;
        printf("    x=%g  f=%s  g=%s\n", x, hex(f.payload()), hex(g.payload()));
      }
    }
    check("(x*x - 1) == (x-1)*(x+1) for all samples", ok);
  }

  printf("\n%s (%d failure%s)\n", g_failures == 0 ? "ALL PASSED" : "FAILURES",
         g_failures, g_failures == 1 ? "" : "s");
  return g_failures == 0 ? 0 : 1;
}
