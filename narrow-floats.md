# FPSan by hand: the integer payloads on narrow floats

FPSan is Triton's; this document just works its construction out by hand on very
small float types. See the [README](README.md) for credits and links.

Every table below is produced and exhaustively verified by
[`narrow_demo.cpp`](narrow_demo.cpp) (built on
[`fpsan_generic.hpp`](fpsan_generic.hpp)); rebuild it to regenerate everything:
`c++ -std=c++20 -O2 narrow_demo.cpp -o narrow_demo && ./narrow_demo`.

## Why narrow floats

FPSan replaces floating-point arithmetic with integer arithmetic modulo `2^w`,
through a bijection `φ : float ↔ ℤ/2^w` (the [README](README.md) explains the
construction for `float32`). The appeal of FPSan is *algebraic*: the integer
image obeys **exact** ring laws and an exponential homomorphism, so two programs
that differ only by fast-math rewrites such as `a+(b+c) → (a+b)+c` become
*bit-identical* even though those rewrites are invalid over IEEE floats.

For `float32` that bijection has 2³² entries — impossible to eyeball. So this
document applies the **exact same recipe** to tiny types where you can see the
whole thing:

- **FP4 (E2M1)** — 4 bits, **16 values**, no Inf/NaN. Small enough to print the
  entire embedding table and check every property by eye.
- **FP8 (E5M2)** — 8 bits, 256 values, a *real* IEEE-754-style type with
  exponent, mantissa, subnormals, and Inf/NaN. Shows how FPSan treats a
  realistic format (including the special values it does **not** model).

> These narrow types go through the **same width-generic construction** that
> Triton uses (`getPayloadMixConfig` in
> [`FpSanToLLVM.cpp`](https://github.com/triton-lang/triton/blob/main/lib/Conversion/TritonInstrumentToLLVM/FpSanToLLVM.cpp)
> is parameterized by bit width). We are not inventing a new scheme; we are
> instantiating Triton's at `w = 4` and `w = 8` for legibility.

---

## Part 1 — FP4 (E2M1)

### The type

A code word is 4 bits: `[sign][exp:2][mantissa:1]`, exponent bias 1, **no Inf or
NaN** (all 16 codes are finite numbers). Decoding:

- normal (`exp ≠ 0`): `(-1)^s · 2^(exp − 1) · (1 + mantissa/2)`
- subnormal (`exp = 0`): `(-1)^s · 2^0 · (mantissa/2)`

So the non-negative values are `{0, 0.5, 1, 1.5, 2, 3, 4, 6}`, plus their
negatives.

### The mixing constants at width 4

Running the generic FPSan recipe at `bitWidth = 4`:

```
oneBits = 0x2   signMask = 0x8   magMask = 0x7   shift = 1
mulA = 3   mulAInv = 3
mulBPos = 5   mulBPosInv = 5
mulBNeg = 3   mulBNegInv = 3
```

The magnitude (low 3 bits) is mixed modulo `2^3 = 8`; the sign bit is carried
separately; the full payload lives in `ℤ/16`.

### The whole embedding table φ

| bits | s e m | value | payload (bin) | payload | −φ(x) mod 16 | φ(−x) |
|------|-------|------:|---------------|--------:|-------------:|------:|
| 0000 | 0 00 0 | 0    | 0000 | 0  | 0  | 8  |
| 0001 | 0 00 1 | 0.5  | 0010 | 2  | 14 | 14 |
| 0010 | 0 01 0 | **1**| 0001 | **1** | 15 | 15 |
| 0011 | 0 01 1 | 1.5  | 0101 | 5  | 11 | 11 |
| 0100 | 0 10 0 | 2    | 0110 | 6  | 10 | 10 |
| 0101 | 0 10 1 | 3    | 0100 | 4  | 12 | 12 |
| 0110 | 0 11 0 | 4    | 0111 | 7  | 9  | 9  |
| 0111 | 0 11 1 | 6    | 0011 | 3  | 13 | 13 |
| 1000 | 1 00 0 | −0   | 1000 | 8  | 8  | 0  |
| 1001 | 1 00 1 | −0.5 | 1110 | 14 | 2  | 2  |
| 1010 | 1 01 0 | **−1**| 1111 | **15** | 1 | 1 |
| 1011 | 1 01 1 | −1.5 | 1011 | 11 | 5  | 5  |
| 1100 | 1 10 0 | −2   | 1010 | 10 | 6  | 6  |
| 1101 | 1 10 1 | −3   | 1100 | 12 | 4  | 4  |
| 1110 | 1 11 0 | −4   | 1001 | 9  | 7  | 7  |
| 1111 | 1 11 1 | −6   | 1101 | 13 | 3  | 3  |

Read straight off the table:

- **It is a bijection.** The `payload` column is a permutation of `0..15`.
- **`φ(+0) = 0`, `φ(+1) = 1`, `φ(−1) = 15` (all-ones)** — the three fixed points
  FPSan pins down so that `x+0=x`, `x·1=x`, and negation behave naturally.
- **`φ(−x) = −φ(x)`** for every nonzero `x`: the last two columns match on every
  row except the two zeros (`φ(−0)=8 ≠ −φ(+0)=0` — the documented "nonzero"
  caveat; `−0` gets the lone payload `8 = 2³`).

### Computing φ by hand

Take `x = 1.5`, code word `0011`:

1. sign bit is 0, so `signFlip = 0`; magnitude `x = 0011₂ = 3`.
2. multiply by `mulA` mod 8: `3 · 3 = 9 ≡ 1 (mod 8)`.
3. xorshift right by `shift = 1`: `1 ^ (1 >> 1) = 1 ^ 0 = 1`.
4. multiply by `mulBPos = 5` mod 8: `1 · 5 = 5`.
5. re-apply sign: `5 ^ 0 = 5`. ✓ (matches the table: `1.5 → 5`).

Now `x = −1.5`, code word `1011`:

1. sign bit is 1, so `signFlip = 8`; magnitude `x = 1011₂ ^ 1000₂ = 0011₂ = 3`.
2. `3 · 3 ≡ 1 (mod 8)`; 3. xorshift: `1`; 4. multiply by `mulBNeg = 3` mod 8:
   `1 · 3 = 3`; 5. re-apply sign: `3 ^ 8 = 11`. ✓

And `mulBNeg = 3 = −5 (mod 8)`, which is exactly why step 4 negates the
magnitude and the whole thing satisfies `φ(−1.5) = 11 = −5 = −φ(1.5) (mod 16)`.

Inverting `φ⁻¹(5)`: `signFlip=0`, `w=5`; `z = 5 · mulBPosInv = 5·5 = 25 ≡ 1`;
inverse-xorshift of `1` is `1`; `y · mulAInv = 1·3 = 3 = 0011₂ = 1.5`. ✓

### Arithmetic = ℤ/16, so reassociation is exact

Once embedded, **addition is just `(p + q) mod 16`** and **multiplication is
`(p · q) mod 16`** — ordinary integer arithmetic, which is associative,
commutative, and distributive *on the nose*. The demo verifies this over **all
16³ = 4096 triples**.

The payoff is the rewrite that wrecks IEEE but not FPSan. Take `a = 6`,
`b = −6`, `c = 1` in FP4:

| | real FP4 (round-to-nearest each step) | FPSan payload |
|---|---|---|
| `(a+b)+c` | `(6 + −6) + 1 = 0 + 1 = 1` | `(φ6 + φ(−6)) + φ1 = (3+13)+1 = 0+1 = `**`1`** |
| `a+(b+c)` | `6 + (−6 + 1) = 6 + (−5→−4) = 2` | `φ6 + (φ(−6)+φ1) = 3 + (13+1) = 3+14 = `**`1`** |

In real FP4 the two orders give **1 vs 2** — because `−5` is not representable
and rounds to `−4`. In FPSan both are payload **1**, exactly. (`16` wraps to
`0`, so `3+13=16≡0` and `13+1=14`, `3+14=17≡1`.)

### Multiplication is *scrambled* but *consistent*

Individual products are numerically meaningless. E.g. `2 · 3`: payloads
`φ(2)=6`, `φ(3)=4`, product `24 ≡ 8 (mod 16)`, and `φ⁻¹(8) = −0`. So FPSan says
"`2 · 3 = −0`" — nonsense as a number. What survives is the *algebra*:
`x·1 = x` (because `φ(1)=1`), and **distributivity holds for all 4096 triples**.
FPSan preserves the structure, not the values.

### Division: `x/x = 1` only for odd payloads

Division is `φ⁻¹(φ(x) · inv(φ(y)))`, where `inv` is the true modular inverse for
**odd** payloads and a parity-preserving involution for **even** ones. So
`x / x = 1` exactly when `φ(x)` is odd:

| value | payload | parity | `x/x == 1`? |
|------:|--------:|:------:|:-----------:|
| 0.5 | 2  | even | no  |
| 1   | 1  | odd  | yes |
| 1.5 | 5  | odd  | yes |
| 2   | 6  | even | no  |
| 3   | 4  | even | no  |
| 4   | 7  | odd  | yes |
| 6   | 3  | odd  | yes |
| −1  | 15 | odd  | yes |
| −2  | 10 | even | no  |

This is the documented caveat, made concrete: half the FP4 values are their own
"divisible-to-1" elements, half are not.

### The exponential

`exp2` is modular exponentiation `C^p mod 16` with `C = 0xA343836D & 0xF = 13`.
Note `13 ≡ 5 (mod 8)`, the property FPSan requires of its generator (truncating
the 32-bit constant preserves the low bits, so this holds at every width). The
homomorphism `exp2(x+y) = exp2(x)·exp2(y)` is verified over all 256 payload
pairs.

Honest caveat: on a ring this tiny the exponential is *degenerate* — `13` has
multiplicative order only 4 mod 16 (`13, 9, 5, 1, …`), so `exp2` is far from
injective. The exponential only becomes interesting at full width, where the
generator's order is huge. What FP4 shows is just that the *identity* holds
exactly; it is not meant to look like a real `exp`.

---

## Part 2 — FP8 (E5M2)

### The type

`[sign][exp:5][mantissa:2]`, bias 15, **with Inf/NaN**. This is the
"IEEE-754-style" FP8: same philosophy as binary16/binary32, just narrower.
`exp = 11111₂` is special (Inf if mantissa 0, else NaN); `exp = 0` is
zero/subnormal; otherwise normal `(-1)^s · 2^(exp−15) · (1 + mantissa/4)`.

### The mixing constants at width 8

```
oneBits = 0x3C   signMask = 0x80   magMask = 0x7F   shift = 2
mulA = 51   mulAInv = 123
mulBPos = 89   mulBPosInv = 105
mulBNeg = 39   mulBNegInv = 23
```

256 values is a lot to print, so here are the landmark code words (the full set
is enumerated and verified by the demo):

| bits | fields | value | payload |
|------|--------|------:|--------:|
| 0x00 | 0 00000 00 | 0 | 0 |
| 0x01 | 0 00000 01 | 1.526e−05 (smallest subnormal) | 103 |
| 0x38 | 0 01110 00 | 0.5 | 82 |
| 0x3C | 0 01111 00 | **1** | **1** |
| 0x3D | 0 01111 01 | 1.25 | 126 |
| 0x40 | 0 10000 00 | 2 | 80 |
| 0x7B | 0 11110 11 | 57344 (largest finite) | 89 |
| 0x7C | 0 11111 00 | **+Inf** | 81 |
| 0x7D | 0 11111 01 | **NaN** | 78 |
| 0x80 | 1 00000 00 | −0 | 128 |
| 0xBC | 1 01111 00 | **−1** | 255 |
| 0xFC | 1 11111 00 | **−Inf** | 175 |

The same three fixed points hold (`φ(0)=0`, `φ(1)=1`, `φ(−1)=255`), and the
exhaustive check confirms bijection, `φ(−x)=−φ(x)`, the full ring laws over all
**256³ ≈ 16.7M triples**, and the `exp2` homomorphism over all 256² pairs.

### Inf and NaN are *just bit patterns*

FPSan does **not** model Inf/NaN semantics — it mixes their bit patterns like
any others. So `+Inf → 81` and `−Inf → 175`, and notice `175 = −81 (mod 256)`:
even infinities obey `φ(−x) = −φ(x)`, because at the bit level negation is just
a sign-bit flip. Arithmetic on these payloads is perfectly well-defined integer
arithmetic, and perfectly meaningless as floating-point. The rule from the docs
stands: **FPSan results are only meaningful compared against other FPSan
results.**

### The same catastrophe, at width 8

`a = 256`, `b = −256`, `c = 1`:

- real FP8: `(a+b)+c = 1`, but `a+(b+c) = 256 + (−255 → −256) = 0` — different,
  because the mantissa step near 256 is 64, so `−255` rounds to `−256`.
- FPSan: `(a+b)+c == a+(b+c)` exactly (verified).

---

## Takeaways

1. FPSan's embedding is one **width-generic recipe**; FP4 just runs it at
   `w = 4`, where the entire bijection fits in a 16-row table you can verify by
   eye.
2. The three pinned fixed points (`φ(0)=0`, `φ(1)=1`, `φ(−1)=` all-ones) plus
   `φ(−x)=−φ(x)` are exactly what make `x+0=x`, `x·1=x`, and negation behave.
3. After embedding, arithmetic is **plain modular-integer arithmetic**, so the
   ring laws — associativity, commutativity, distributivity — are *exact*. That
   is the whole point: numerically-invalid fast-math rewrites become exact
   integer identities.
4. The values are scrambled; only the **algebra** is preserved. Division is
   exact-inverse only for odd payloads; Inf/NaN are unmodeled bit patterns; the
   exponential relies (at real widths) on Schanuel's conjecture.

To regenerate every number here, build and run
[`narrow_demo.cpp`](narrow_demo.cpp).
