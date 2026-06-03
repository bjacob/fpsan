// narrow_demo.cpp
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
// Works out the FPSan integer payloads on narrow floating-point types, and
// verifies the algebraic properties by exhaustive enumeration.  Every table in
// NARROW_FLOATS.md is produced by this program, so the document stays honest.
//
// Build:   c++ -std=c++20 -O2 narrow_demo.cpp -o narrow_demo
// Run:     ./narrow_demo
// ----------------------------------------------------------------------------
#include "fpsan_f32.hpp"
#include "fpsan_generic.hpp"

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

using namespace fpsan_generic;

static int g_failures = 0;
static void check(const char *name, bool ok) {
  printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
  if (!ok)
    ++g_failures;
}

static std::string bin(uint32_t v, unsigned bw) {
  std::string s;
  for (int i = int(bw) - 1; i >= 0; --i)
    s += ((v >> i) & 1u) ? '1' : '0';
  return s;
}

// Split a code word into sign|exp|mantissa with spaces, e.g. "0 11 1".
static std::string bin_fields(const FPFormat &fmt, uint32_t v) {
  std::string s;
  s += ((v >> (fmt.bitWidth - 1)) & 1u) ? '1' : '0';
  s += ' ';
  for (int i = int(fmt.expBits) - 1; i >= 0; --i)
    s += ((v >> (fmt.mantBits + i)) & 1u) ? '1' : '0';
  s += ' ';
  for (int i = int(fmt.mantBits) - 1; i >= 0; --i)
    s += ((v >> i) & 1u) ? '1' : '0';
  return s;
}

// negate a code word: flip the sign bit.
static uint32_t neg_bits(const FPFormat &fmt, uint32_t bits) {
  return bits ^ (1u << (fmt.bitWidth - 1));
}

static void print_config(const FPFormat &fmt) {
  const MixConfig &c = mix_config_for(fmt);
  printf(
      "format %s: bitWidth=%u (sign 1, exp %u, mant %u), bias=%d, inf/nan=%s\n",
      fmt.name, fmt.bitWidth, fmt.expBits, fmt.mantBits, fmt.bias,
      fmt.hasInfNan ? "yes" : "no");
  printf("  oneBits=0x%X  signMask=0x%llX  magMask=0x%llX  shift=%u\n",
         fmt.oneBits(), (unsigned long long)c.signMask,
         (unsigned long long)c.magMask, c.shift);
  printf("  mulA=%llu mulAInv=%llu  mulBPos=%llu mulBPosInv=%llu  "
         "mulBNeg=%llu mulBNegInv=%llu\n",
         (unsigned long long)c.mulA, (unsigned long long)c.mulAInv,
         (unsigned long long)c.mulBPos, (unsigned long long)c.mulBPosInv,
         (unsigned long long)c.mulBNeg, (unsigned long long)c.mulBNegInv);
}

// Exhaustively verify the FPSan algebra over a format whose every code word we
// can enumerate (bitWidth small).
static void verify_exhaustive(const FPFormat &fmt) {
  const MixConfig &c = mix_config_for(fmt);
  uint32_t N = static_cast<uint32_t>(c.fullMask) + 1u;

  // phi is a bijection that round-trips.
  std::vector<int> hit(N, 0);
  bool roundtrip = true;
  for (uint32_t b = 0; b < N; ++b) {
    uint32_t p = mix_float_bits_to_int(c, b);
    hit[p]++;
    if (unmix_int_to_float_bits(c, p) != b)
      roundtrip = false;
  }
  bool bijection = true;
  for (uint32_t p = 0; p < N; ++p)
    if (hit[p] != 1)
      bijection = false;
  check("phi is a bijection on all code words", bijection);
  check("phi^-1(phi(x)) == x for all code words", roundtrip);

  // Fixed points.
  check("phi(+0) == 0", mix_float_bits_to_int(c, 0u) == 0u);
  check("phi(+1) == 1", mix_float_bits_to_int(c, fmt.oneBits()) == 1u);
  check("phi(-1) == all-ones",
        mix_float_bits_to_int(c, neg_bits(fmt, fmt.oneBits())) ==
            static_cast<uint32_t>(c.fullMask));

  // phi(-x) == -phi(x) for nonzero magnitudes.
  bool negsym = true;
  for (uint32_t b = 0; b < N; ++b) {
    if ((b & c.magMask) == 0)
      continue; // skip +0/-0
    uint32_t lhs = mix_float_bits_to_int(c, neg_bits(fmt, b));
    uint32_t rhs =
        static_cast<uint32_t>((0u - mix_float_bits_to_int(c, b)) & c.fullMask);
    if (lhs != rhs)
      negsym = false;
  }
  check("phi(-x) == -phi(x) for nonzero x", negsym);

  // Ring laws over all payloads (payloads ARE Z/2^bitWidth, so these are the
  // ordinary modular-integer laws; we still check to be concrete).
  bool addAssoc = true, mulAssoc = true, distrib = true, addComm = true,
       mulComm = true, addId = true, mulId = true, negInv = true;
  FPSanFloat ZERO(fmt, 0u), ONE(fmt, 1u);
  for (uint32_t i = 0; i < N; ++i) {
    FPSanFloat A(fmt, i);
    if (A + ZERO != A)
      addId = false;
    if (A * ONE != A)
      mulId = false;
    if (A + (-A) != ZERO)
      negInv = false;
    for (uint32_t j = 0; j < N; ++j) {
      FPSanFloat B(fmt, j);
      if (A + B != B + A)
        addComm = false;
      if (A * B != B * A)
        mulComm = false;
      for (uint32_t k = 0; k < N; ++k) {
        FPSanFloat C(fmt, k);
        if ((A + B) + C != A + (B + C))
          addAssoc = false;
        if ((A * B) * C != A * (B * C))
          mulAssoc = false;
        if (A * (B + C) != (A * B) + (A * C))
          distrib = false;
      }
    }
  }
  check("x + 0 == x", addId);
  check("x * 1 == x", mulId);
  check("x + (-x) == 0", negInv);
  check("+ commutative", addComm);
  check("* commutative", mulComm);
  check("+ associative over ALL payload triples", addAssoc);
  check("* associative over ALL payload triples", mulAssoc);
  check("distributive over ALL payload triples", distrib);

  // exp2 homomorphism over all payloads.
  bool hom = true;
  for (uint32_t i = 0; i < N; ++i)
    for (uint32_t j = 0; j < N; ++j) {
      FPSanFloat A(fmt, i), B(fmt, j);
      if ((A + B).exp2() != A.exp2() * B.exp2())
        hom = false;
    }
  check("exp2(x+y) == exp2(x)*exp2(y) over ALL payloads", hom);
}

static void section(const char *title) {
  printf("\n========================================================\n");
  printf("%s\n", title);
  printf("========================================================\n");
}

int main() {
  section("FP4 (E2M1): the embedding table -- all 16 values");
  const FPFormat &fp4 = formats::E2M1;
  print_config(fp4);
  printf("\n bits   s e m  | value | payload  dec | -phi  | phi(-x)\n");
  printf(" ----------------------------------------------------------\n");
  for (uint32_t b = 0; b < 16; ++b) {
    uint32_t p = mix_float_bits_to_int(mix_config_for(fp4), b);
    uint32_t negp = (0u - p) & 0xFu;
    uint32_t phineg =
        mix_float_bits_to_int(mix_config_for(fp4), neg_bits(fp4, b));
    printf(" %s  %s | %5s | %s  %3u | %4u  | %4u\n", bin(b, 4).c_str(),
           bin_fields(fp4, b).c_str(), value_name(fp4, b).c_str(),
           bin(p, 4).c_str(), p, negp, phineg);
  }

  section("FP4 (E2M1): exhaustive verification");
  verify_exhaustive(fp4);

  section("FP4 (E2M1): the reassociation that breaks IEEE but not FPSan");
  {
    // a=6, b=-6, c=1. In real FP4: (a+b)+c = 1; a+(b+c) = 6+(-5 rounds to
    // -4)=2.
    uint32_t aB = encode_nearest(fp4, 6.0), bB = encode_nearest(fp4, -6.0),
             cB = encode_nearest(fp4, 1.0);
    printf("  a=6 (bits %s), b=-6 (bits %s), c=1 (bits %s)\n",
           bin(aB, 4).c_str(), bin(bB, 4).c_str(), bin(cB, 4).c_str());
    // "real" FP4 with round-to-nearest-even at each step:
    auto fp4_add = [&](uint32_t x, uint32_t y) {
      return encode_nearest(fp4, decode(fp4, x).value + decode(fp4, y).value);
    };
    uint32_t ieeeL = fp4_add(fp4_add(aB, bB), cB);
    uint32_t ieeeR = fp4_add(aB, fp4_add(bB, cB));
    printf("  real FP4:  (a+b)+c = %s,  a+(b+c) = %s   -> %s\n",
           value_name(fp4, ieeeL).c_str(), value_name(fp4, ieeeR).c_str(),
           ieeeL == ieeeR ? "equal" : "DIFFER");
    FPSanFloat A = FPSanFloat::embed(fp4, aB), B = FPSanFloat::embed(fp4, bB),
               C = FPSanFloat::embed(fp4, cB);
    FPSanFloat L = (A + B) + C, R = A + (B + C);
    printf("  FPSan:     (a+b)+c payload=%u,  a+(b+c) payload=%u   -> %s\n",
           L.payload(), R.payload(), L == R ? "EQUAL" : "differ");
    check("FPSan reassociation is exact where IEEE FP4 is not", L == R);
  }

  section("FP4 (E2M1): division -- x/x == 1 only for ODD payloads");
  {
    const MixConfig &c = mix_config_for(fp4);
    printf("  value | payload | parity | x/x == 1 ?\n");
    for (uint32_t b = 0; b < 16; ++b) {
      if ((b & c.magMask) == 0)
        continue; // skip zeros
      FPSanFloat X = FPSanFloat::embed(fp4, b);
      bool one = (X / X == FPSanFloat(fp4, 1u));
      printf("  %5s | %5u   | %-4s   | %s\n", value_name(fp4, b).c_str(),
             X.payload(), (X.payload() & 1u) ? "odd" : "even",
             one ? "yes" : "no");
    }
  }

  section("FP8 (E5M2): structure, special values, and selected payloads");
  const FPFormat &fp8 = formats::E5M2;
  print_config(fp8);
  printf("\n  Special / landmark code words:\n");
  printf("   %-8s %-10s %-12s payload\n", "bits", "fields", "value");
  uint32_t landmarks[] = {0x00, 0x01, 0x3C, 0x40, 0x38, 0x3D,
                          0x7B, 0x7C, 0x7D, 0xBC, 0xFC, 0x80};
  for (uint32_t b : landmarks) {
    uint32_t p = mix_float_bits_to_int(mix_config_for(fp8), b);
    printf("   0x%02X     %-10s %-12s %u\n", b, bin_fields(fp8, b).c_str(),
           value_name(fp8, b).c_str(), p);
  }

  section("FP8 (E5M2): exhaustive verification over all 256 code words");
  verify_exhaustive(fp8);

  section("FP8 (E5M2): a numerically catastrophic reassociation");
  {
    // 256 + (-256 + 1): -256+1 = -255 -> rounds in E5M2 (mantissa step at that
    // magnitude is 64) to -256; so a+(b+c) = 0, but (a+b)+c = 1.
    uint32_t aB = encode_nearest(fp8, 256.0), bB = encode_nearest(fp8, -256.0),
             cB = encode_nearest(fp8, 1.0);
    auto fp8_add = [&](uint32_t x, uint32_t y) {
      return encode_nearest(fp8, decode(fp8, x).value + decode(fp8, y).value);
    };
    uint32_t ieeeL = fp8_add(fp8_add(aB, bB), cB);
    uint32_t ieeeR = fp8_add(aB, fp8_add(bB, cB));
    printf("  a=256, b=-256, c=1\n");
    printf("  real FP8:  (a+b)+c = %s,  a+(b+c) = %s   -> %s\n",
           value_name(fp8, ieeeL).c_str(), value_name(fp8, ieeeR).c_str(),
           ieeeL == ieeeR ? "equal" : "DIFFER");
    FPSanFloat A = FPSanFloat::embed(fp8, aB), B = FPSanFloat::embed(fp8, bB),
               C = FPSanFloat::embed(fp8, cB);
    check("FPSan E5M2 reassociation is exact", (A + B) + C == A + (B + C));
  }

  section("Cross-check: generic F32 path matches fpsan_f32.hpp");
  {
    bool ok = true;
    float xs[] = {0.0f, 1.0f, -1.0f, 3.14159f, -2.5f, 1e6f, 0.1f};
    for (float x : xs) {
      uint32_t generic =
          FPSanFloat::embed(formats::F32, fpsan::bits_of(x)).payload();
      uint32_t concrete = fpsan::FPSanF32::phi(x).payload();
      if (generic != concrete)
        ok = false;
    }
    check("generic FPSanFloat(F32) payloads == FPSanF32 payloads", ok);
  }

  printf("\n%s (%d failure%s)\n", g_failures == 0 ? "ALL PASSED" : "FAILURES",
         g_failures, g_failures == 1 ? "" : "s");
  return g_failures == 0 ? 0 : 1;
}
