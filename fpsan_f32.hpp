// fpsan_f32.hpp
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
// A stand-alone, header-only C++ re-implementation of the core ideas behind
// Triton's "FPSan" floating-point sanitizer, restricted to float32.
//
// FPSan rewrites floating-point arithmetic into *purely integer* arithmetic
// modulo 2^32.  It is not meant to compute anything numerically useful; the
// numerical outputs are scrambled.  Its purpose is to give floating-point a set
// of EXACT algebraic laws: associativity, commutativity, distributivity, and an
// exponential homomorphism exp(x+y) = exp(x)*exp(y).  These laws do not hold
// for real IEEE arithmetic, so two kernels that are "the same up to
// fast-math-style rewrites" (e.g. a+(b+c) vs (a+b)+c) disagree under IEEE but
// agree EXACTLY under FPSan.
//
// The construction (see the blog post linked from Triton's docs,
// docs/programming-guide/chapter-3/fpsan.rst):
//
//   * A bijection  phi : float32  <->  Z/2^32   ("embed" / "unembed").
//       phi(+0.0) = 0,  phi(+1.0) = 1,  phi(-1.0) = 0xFFFFFFFF (= -1 mod 2^32),
//       phi(-x) = -phi(x)  for nonzero x.
//   * FP ops become the corresponding operations on the integer "payload":
//       x + y   ->  phi^{-1}( phi(x) + phi(y) )        (mod 2^32)
//       x - y   ->  phi^{-1}( phi(x) - phi(y) )
//       -x      ->  phi^{-1}( -phi(x) )
//       x * y   ->  phi^{-1}( phi(x) * phi(y) )
//       x / y   ->  phi^{-1}( phi(x) * inv(phi(y)) )
//       exp2(x) ->  phi^{-1}( C ^ phi(x) )             (modular exponentiation)
//       exp(x)  ->  exp2( x scaled by 1/ln(2) in payload space )
//
// The exact mixing/operation constants below are taken verbatim from the Triton
// sources so that this class produces the same payloads Triton would:
//   lib/Conversion/TritonInstrumentToLLVM/FpSanToLLVM.cpp   (phi / phi^{-1})
//   lib/Dialect/TritonInstrument/Transforms/FpSanitizer.cpp (div/exp/exp2)
//
// Requires C++20 (uses std::bit_cast); no other dependencies.
// ----------------------------------------------------------------------------
#pragma once

#include <bit>
#include <cstdint>

namespace fpsan {

// ---- bit_cast helpers ------------------------------------------------------
inline uint32_t bits_of(float f) { return std::bit_cast<uint32_t>(f); }
inline float float_of(uint32_t u) { return std::bit_cast<float>(u); }

// ---- modular inverse of an odd 64-bit number (Newton/Hensel lifting) -------
// Verbatim port of invOddU64 from the Triton sources.  For odd `a`, returns the
// inverse of `a` modulo 2^64.
inline uint64_t inv_odd_u64(uint64_t a) {
  // a must be odd.
  uint64_t x = 2 - a;
  for (unsigned correctBits = 2; correctBits < 64; correctBits *= 2)
    x *= 2 - a * x;
  return x;
}

// ---------------------------------------------------------------------------
//  The payload-mixing configuration ("phi") for float32.
//
//  Reproduces getPayloadMixConfig() from FpSanToLLVM.cpp for a 32-bit float.
//  All magnitude arithmetic happens modulo 2^31 (the sign bit is carried
//  separately), so we keep the low 31 bits and mask with magMask = 2^31 - 1.
// ---------------------------------------------------------------------------
struct MixConfig {
  static constexpr unsigned kBitWidth = 32;
  static constexpr uint32_t kSignMask = 1u << (kBitWidth - 1); // 0x80000000
  static constexpr uint32_t kMagMask = kSignMask - 1;          // 0x7FFFFFFF

  // 1.0f has IEEE-754 bit pattern 0x3F800000; its number of trailing zero bits
  // is the mantissa width (23), used as the xorshift amount.
  static constexpr uint32_t kOneBits = 0x3F800000u;
  static constexpr unsigned kShift = 23; // = countr_zero(0x3F800000)

  uint32_t mulA;    // odd multiplier, mixes low->high bits, invertible mod 2^31
  uint32_t mulAInv; // inverse of mulA mod 2^31
  uint32_t
      mulBPos; // second multiplier (non-negative branch); maps mixed 1.0 -> 1
  uint32_t mulBNeg; // second multiplier (negative branch) = -mulBPos mod 2^31
  uint32_t mulBPosInv; // inverse of mulBPos mod 2^31
  uint32_t mulBNegInv; // inverse of mulBNeg mod 2^31

  static const MixConfig &get() {
    static const MixConfig cfg = build();
    return cfg;
  }

private:
  static MixConfig build() {
    MixConfig c;
    // First multiplier: an arbitrary odd constant (922291) reduced mod 2^31.
    c.mulA = 922291u & kMagMask;

    // Mix 1.0 through the first stage so we can choose the second multiplier
    // such that phi(1.0) == 1.
    uint64_t oneMixed = (static_cast<uint64_t>(kOneBits) * c.mulA) & kMagMask;
    oneMixed ^= oneMixed >> kShift; // must be odd by construction

    c.mulBPos = static_cast<uint32_t>(inv_odd_u64(oneMixed)) & kMagMask;
    // magMask == -1 mod 2^31, so mulBNeg == -mulBPos mod 2^31.  This is what
    // makes phi(-x) = -phi(x).
    c.mulBNeg = static_cast<uint32_t>(
        (static_cast<uint64_t>(c.mulBPos) * static_cast<uint64_t>(kMagMask)) &
        kMagMask);

    c.mulAInv = static_cast<uint32_t>(inv_odd_u64(c.mulA)) & kMagMask;
    c.mulBPosInv = static_cast<uint32_t>(inv_odd_u64(c.mulBPos)) & kMagMask;
    c.mulBNegInv = static_cast<uint32_t>(inv_odd_u64(c.mulBNeg)) & kMagMask;
    return c;
  }
};

// xorshift-right step: v ^= v >> shift.  For 31-bit magnitudes with shift 23,
// this is an involution (its own inverse), since (v>>23)>>23 == 0.
inline uint32_t xorshift_right(uint32_t v, unsigned shift) {
  return v ^ (v >> shift);
}

// ---------------------------------------------------------------------------
//  phi : float32 -> Z/2^32   (the "embed" / mixFloatToInt operation)
//
//  Port of mixFloatToInt() from FpSanToLLVM.cpp.  Input is the raw IEEE bits.
// ---------------------------------------------------------------------------
inline uint32_t mix_float_bits_to_int(uint32_t u) {
  const MixConfig &cfg = MixConfig::get();
  // signFlip selects whether to flip the sign bit: nonzero only when the input
  // is negative.  XORing with it clears the sign bit going in and restores it
  // coming out, so the sign bit is carried verbatim while the magnitude is
  // mixed.
  uint32_t signFlip = (u & MixConfig::kSignMask) ? MixConfig::kSignMask : 0u;
  uint32_t x = u ^ signFlip; // magnitude bits, sign cleared (top bit 0)

  uint32_t y = static_cast<uint32_t>((static_cast<uint64_t>(x) * cfg.mulA) &
                                     MixConfig::kMagMask);
  uint32_t z = xorshift_right(y, MixConfig::kShift);

  uint32_t mulB = (u & MixConfig::kSignMask) ? cfg.mulBNeg : cfg.mulBPos;
  uint32_t w = static_cast<uint32_t>((static_cast<uint64_t>(z) * mulB) &
                                     MixConfig::kMagMask);

  return w ^ signFlip; // restore the original sign bit
}

// ---------------------------------------------------------------------------
//  phi^{-1} : Z/2^32 -> float32   (the "unembed" / unmixIntToFloat operation)
//
//  Port of unmixIntToFloat() from FpSanToLLVM.cpp.  Returns raw IEEE bits.
// ---------------------------------------------------------------------------
inline uint32_t unmix_int_to_float_bits(uint32_t v) {
  const MixConfig &cfg = MixConfig::get();
  uint32_t signFlip = (v & MixConfig::kSignMask) ? MixConfig::kSignMask : 0u;
  uint32_t w = v ^ signFlip;

  uint32_t mulBInv =
      (v & MixConfig::kSignMask) ? cfg.mulBNegInv : cfg.mulBPosInv;
  uint32_t z = static_cast<uint32_t>((static_cast<uint64_t>(w) * mulBInv) &
                                     MixConfig::kMagMask);

  // inverse of the xorshift: for shift = cfg.shift; shift < bitWidth; shift*=2.
  // With shift = 23 and bitWidth = 32, only the single shift-23 step runs, and
  // it is its own inverse.
  uint32_t y = z;
  for (unsigned shift = MixConfig::kShift; shift < MixConfig::kBitWidth;
       shift *= 2)
    y = xorshift_right(y, shift);

  uint32_t x = static_cast<uint32_t>((static_cast<uint64_t>(y) * cfg.mulAInv) &
                                     MixConfig::kMagMask);

  return x ^ signFlip;
}

// ---------------------------------------------------------------------------
//  Division helper: inv() on payloads.
//  Port of fpsanIntInv() from FpSanitizer.cpp.
//    * odd payloads -> their true modular inverse mod 2^32
//    * even payloads -> a parity-preserving involution
// ---------------------------------------------------------------------------
inline uint32_t fpsan_int_inv(uint32_t u) {
  uint32_t a = u | 1u; // force odd
  uint32_t x = 2u - a; // Newton iteration seed (mod 2^32)
  for (unsigned correctBits = 2; correctBits < 32; correctBits *= 2)
    x = x * (2u - a * x);
  uint32_t evenPart = x & 0xFFFFFFFEu; // drop the low bit
  uint32_t originalParity = u & 1u;    // ...and put back u's parity
  return evenPart | originalParity;
}

// ---------------------------------------------------------------------------
//  Modular exponentiation in payload space.
//  Port of fpsanExp2FromInt() from FpSanitizer.cpp.
//  Computes C ^ xI mod 2^32 (left-to-right square-and-multiply), with the fixed
//  odd generator C = 0xA343836D (note C == 5 mod 8).  This is what gives
//      exp2(x + y) = exp2(x) * exp2(y).
// ---------------------------------------------------------------------------
inline uint32_t fpsan_exp2_payload(uint32_t xI) {
  constexpr uint32_t C = 0xA343836Du;
  uint32_t y = 1u;
  for (int i = 0; i < 32; ++i) {
    y = y * y;
    int bitIndex = 31 - i;
    uint32_t bit = 1u << bitIndex;
    uint32_t factor = (xI & bit) ? C : 1u;
    y = y * factor;
  }
  return y;
}

// exp(x) = exp2(x / ln 2).  In payload space the "/ ln 2" rescaling is the
// modular multiply by 0x236EE9BF, taken verbatim from fpsanExp() in
// FpSanitizer.cpp.
inline uint32_t fpsan_exp_payload(uint32_t xI) {
  constexpr uint32_t kRcpLog2 = 0x236EE9BFu;
  return fpsan_exp2_payload(xI * kRcpLog2);
}

// ===========================================================================
//  FPSanF32: a float32 whose arithmetic is the FPSan integer-payload algebra.
//
//  The whole state is the 32-bit payload.  Every operator is exact modular
//  integer arithmetic, so the algebraic identities are EXACT.
// ===========================================================================
class FPSanF32 {
public:
  FPSanF32() : payload_(0) {}

  // --- the embedding phi: take an ordinary float and produce its FPSan value.
  static FPSanF32 phi(float x) {
    return FPSanF32::from_payload(mix_float_bits_to_int(bits_of(x)));
  }

  // --- the inverse embedding phi^{-1}: read the payload back out as a float.
  // (Numerically meaningless on its own; useful to inspect/compare results.)
  float to_float() const { return float_of(unmix_int_to_float_bits(payload_)); }

  // Construct directly from a raw payload (e.g. for a constant or test).
  static FPSanF32 from_payload(uint32_t p) {
    FPSanF32 v;
    v.payload_ = p;
    return v;
  }
  uint32_t payload() const { return payload_; }

  // --- ring operations: pure modular arithmetic on the payloads -------------
  FPSanF32 operator+(const FPSanF32 &o) const {
    return from_payload(payload_ + o.payload_);
  }
  FPSanF32 operator-(const FPSanF32 &o) const {
    return from_payload(payload_ - o.payload_);
  }
  FPSanF32 operator*(const FPSanF32 &o) const {
    return from_payload(payload_ * o.payload_);
  }
  FPSanF32 operator/(const FPSanF32 &o) const {
    return from_payload(payload_ * fpsan_int_inv(o.payload_));
  }
  FPSanF32 operator-() const { return from_payload(0u - payload_); }

  FPSanF32 &operator+=(const FPSanF32 &o) { return *this = *this + o; }
  FPSanF32 &operator-=(const FPSanF32 &o) { return *this = *this - o; }
  FPSanF32 &operator*=(const FPSanF32 &o) { return *this = *this * o; }
  FPSanF32 &operator/=(const FPSanF32 &o) { return *this = *this / o; }

  // --- transcendental ops ---------------------------------------------------
  FPSanF32 exp() const { return from_payload(fpsan_exp_payload(payload_)); }
  FPSanF32 exp2() const { return from_payload(fpsan_exp2_payload(payload_)); }

  // --- exact payload equality (the thing FPSan checks compare) --------------
  bool operator==(const FPSanF32 &o) const { return payload_ == o.payload_; }
  bool operator!=(const FPSanF32 &o) const { return payload_ != o.payload_; }

  // Handy named constants matching FPSan's stable fixed points.
  static FPSanF32 zero() { return from_payload(0u); }               // phi(+0.0)
  static FPSanF32 one() { return from_payload(1u); }                // phi(+1.0)
  static FPSanF32 minus_one() { return from_payload(0xFFFFFFFFu); } // phi(-1.0)

private:
  uint32_t payload_;
};

// Free-function forms so generic code reads naturally.
inline FPSanF32 exp(const FPSanF32 &x) { return x.exp(); }
inline FPSanF32 exp2(const FPSanF32 &x) { return x.exp2(); }
inline FPSanF32 phi(float x) { return FPSanF32::phi(x); }

} // namespace fpsan
