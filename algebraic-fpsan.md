# An algebraic alternative to FPSan: the residue homomorphism $\mathbb{Z}[1/2] \to \mathbb{F}_p$

*Working notes. Audience: us. Builds on [understanding-fpsan.md](understanding-fpsan.md).*

## 0. The one-sentence difference from FPSan

FPSan encodes a float by a **scrambling bijection** $\varphi$ of its bit pattern
and then does arithmetic in the payload ring $\mathbb{Z}/2^{w}$. That map is *not*
a ring homomorphism: $\varphi(x{+}y) \neq \varphi(x){+}\varphi(y)$. The ring
laws hold because the *payload operations* are the ring operations and the
*leaves* behave like free generators — FPSan realizes (an approximation of) the
**free** commutative ring/exponential-ring on the program's inputs. Scrambling
is what makes the leaves look free: the only coincidences are the ones forced by
the axioms.

The idea here is the opposite. Make the **encoding itself** a ring homomorphism,
so the model computes the *true value* of the expression, reduced into a finite
field. Instead of the free algebra on the inputs, we get the actual element of
$\mathbb{Q}$ the program denotes, fingerprinted mod $p$.

That single change — *free model* $\to$ *value model* — is the whole story, and
every property below is a consequence of it.

## 1. The map

Finite binary floats are dyadic rationals: a finite `float32` is $m\cdot 2^{e}$
with $m,e\in\mathbb{Z}$ and $|m| < 2^{24}$ (significand incl. implicit bit; sign
folded into $m$), $e \in [-149, 104]$. They live in

$$\mathbb{Z}[1/2] \;=\; \{\, a/2^k : a\in\mathbb{Z},\,k\ge 0 \,\} \;=\; \mathbb{Z}\big[\tfrac12\big],$$

the localization of $\mathbb{Z}$ at the powers of $2$ — "$\mathbb{Z}$ with $2$
inverted." For any **odd** prime $p$, $2$ is a unit mod $p$, so there is a
**unique** ring homomorphism

$$\varphi_p : \mathbb{Z}[1/2] \longrightarrow \mathbb{F}_p, \qquad
\varphi_p\!\left(\tfrac{a}{2^k}\right) = \bar a \cdot (\bar 2^{-1})^{k} \bmod p,$$

forced by $1\mapsto 1$ and $\tfrac12 \mapsto \overline{2}^{-1}$. It is a genuine
homomorphism: it respects $+,-,\times$ on the nose.

It is cleaner to name the largest domain. Let $\mathbb{Z}_{(p)} = \{a/b : p\nmid b\}$
be the localization of $\mathbb{Z}$ at the prime $p$; it is a local ring with
maximal ideal $p\mathbb{Z}_{(p)}$ and residue field $\mathbb{Z}_{(p)}/p\mathbb{Z}_{(p)} \cong \mathbb{F}_p$.
Then

$$\mathbb{Z}[1/2] \;\subset\; \mathbb{Z}_{(p)} \;\xrightarrow{\ \text{residue}\ }\; \mathbb{F}_p$$

and $\varphi_p$ is just the restriction of the residue map. The point of the
larger domain: **division** by any float whose residue is nonzero lands back in
$\mathbb{Z}_{(p)}$, even when the exact quotient (e.g. $1/3$) is not dyadic.
$\mathbb{F}_p$ being a *field* is what makes division total.

**The model.** Carry, for each value, its residue $r=\varphi_p(v)\in\mathbb{F}_p$.
Evaluate $+,-,\times,\div$ as field operations on residues. The residue of an
expression is then $\varphi_p$ of its *exact* (un-rounded) rational value, so two
expressions with the same exact value get the same residue — automatically, for
**every** identity true over $\mathbb{Q}$, not just the ones we thought to build
in. That is strictly more than FPSan's free model honors.

We restrict attention to **compact 32-bit encodings** throughout. Two variants
are developed in §8: a single prime $p$ just below $2^{32}$, and a composite
$n=pd$ (still $<2^{32}$) that buys an exact exponential. A working prototype of
both lives in [`algebraic_fpsan_generic.hpp`](algebraic_fpsan_generic.hpp) +
[`algebraic_demo.cpp`](algebraic_demo.cpp).

## 2. The central dichotomy: free model vs. value model

This deserves to be made sharp, because it decides everything and it cuts **both
ways**.

- **FPSan = free.** Inputs are scrambled into algebraically-independent
  generators; the model is the free (exponential) ring on them. Two expressions
  match iff they are provably equal *from the axioms alone*. Consequences:
  $2{+}2$ does **not** match $4$; $x{\cdot}(1/x)$ does **not** match $1$ in
  general; $a{+}a$ does **not** match $2a$. The model is blind to the actual
  numerical values of the leaves.
- **$\varphi_p$ = value.** Leaves are their true residues; the model is the
  actual element of $\mathbb{Q}$. Two expressions match iff they have the **same
  exact value** (mod $p$). So $2{+}2 = 4$, $x/x = 1$, $a{+}a = 2a$,
  $(a{+}b)^2 = a^2{+}2ab{+}b^2$, constant-folding, common-subexpression — all
  match.

Neither dominates universally:

- On the **algebraic fragment** (polynomials/rational functions of the inputs),
  the true values have *many* relations ($2{+}2=4$ really is $4$), so the value
  model is the *faithful* one and the free model over-discriminates (it would
  flag value-preserving rewrites as differences).
- On **transcendentals** the situation flips — see §6 — because transcendental
  values are algebraically *independent*, so freeness becomes the faithful
  choice there.

So the honest framing of this whole note: $\varphi_p$ is the right tool for the
*algebraic* core of a float program; FPSan-style freeness is the right tool for
the *transcendental* parts; and the interesting object is the hybrid (§8).

## 3. Injectivity and collisions

We give up injectivity ($p < 2^{32}$). Three separate questions hide under
"how bad are collisions."

**(a) Collision with zero — the one that would hurt — does not happen.**
$\varphi_n(m\cdot 2^e)=0 \iff n\mid m \iff m=0$, because $|m|<2^{24}$ and the
modulus $n>2^{24}$ (both variants: $n\approx 2^{32}$ for 1b, $\approx 2^{31}$ for
CRT). So the *only* finite float with residue $0$ is $\pm 0.0$ — the crucial
design constraint that keeps $0$ honest and the zero/infinity story (§7) clean.
Nothing nonzero ever silently becomes zero. (This is about residue-$0$; the CRT
variant separately has nonzero *zero-divisors* — §6, §8 — a distinct,
division-only issue.)

**(b) Leaf collisions (two distinct floats, same residue) are common but
mostly harmless.** With $N\approx 2^{32}$ finite floats dropped into $p\approx
2^{32}$ residues, occupancy is Poisson($\approx 1$): a constant fraction of
floats share their residue with another float. (The image
$\{\bar m\,\bar 2^{\,e}\}$ is a union of $\sim 250$ scaled copies of the
significand interval — structured, but spread across $\mathbb{F}_p$.) This only
matters for a *literal-substitution* bug: replacing input float $x$ by a
different $y$ with $\varphi_p(x)=\varphi_p(y)$. For the main use case —
comparing **re-orderings / refactors of the same computation on the same
inputs** — both sides share the same leaves, so leaf collisions are irrelevant.

**(c) Circuit-level false matches are fingerprinting at a random prime.**
Two algebraically *different* expressions on shared inputs collide iff
$p \mid \operatorname{num}(v_A - v_B)$. The difference $v_A-v_B = N/2^k$ has
$|N| \lesssim 2^{278}$ (significand $<2^{24}$, exponent span $\sim 254$), so $N$
has at most $\sim 278/32 \approx 9$ prime factors above $2^{32}$. With $\sim
2^{26.5}$ primes in $[2^{31},2^{32}]$, a **random** such prime divides $N$ with
probability $\lesssim 9/2^{26.5} \approx 2^{-23.5}$.

This is exactly Freivalds / Schwartz–Zippel / Rabin-fingerprint territory:
checking an identity by evaluating at a random point (here, reducing mod a
random prime). Consequences:

- Per comparison, false-match $\approx 2^{-24}$ to $2^{-32}$ with a 32-bit
  modulus. Staying compact (32-bit), the lever for more margin is **repetition**:
  $k$ runs with $k$ different primes multiply the rates (independent
  fingerprints / CRT). Two runs $\approx 2^{-48}$–$2^{-64}$. The sanitizer is a
  test tool, so spending *time* (extra runs) rather than *space* (a wider word)
  is the right trade — and it keeps the instrumented program's footprint
  unchanged (§8).
- Unlike FPSan's *fixed* scramble — which has a fixed set of colliding pairs that
  collide on every run — randomizing $p$ makes collisions **non-repeatable and
  adversarially robust**: there is no bad pair that always fools you.

The prototype confirms the rate: over 300k distinct random floats, observed leaf
collisions track $\binom{N}{2}/n$ to within noise (16 vs ~10.5 at $n\approx2^{32}$,
21 vs ~21.1 at $n\approx2^{31}$).

Verdict: the loss of injectivity is essentially cosmetic for the intended use.
The substantive cost is a slightly higher per-prime collision rate, bought back
by repeating with another prime.

## 4. Scrambling: a hard trade-off, not a dial

FPSan scrambles to make leaves free. Can we keep $\varphi_p$'s homomorphism
*and* re-introduce scrambling on top, $v \mapsto \sigma(\varphi_p(v))$? Only if
$\sigma$ is a ring automorphism of $\mathbb{F}_p$ — otherwise the homomorphism
(the entire point) is destroyed. But $\operatorname{Aut}(\mathbb{F}_p)=\{\mathrm{id}\}$
(the Frobenius $x\mapsto x^p$ is the identity on the prime field). **So you
cannot scramble inside $\mathbb{F}_p$ at all.** Homomorphism and value-scrambling
are mutually exclusive; this is not a tunable trade-off but a dichotomy.

What plays the role of the hash seed is instead the **choice of $p$**. A random
prime is the randomness; collisions depend on the seed exactly as in
fingerprinting (§3c). So we don't lose randomization — we relocate it from
"scramble the values" to "randomize the evaluation point," which is the more
principled place for it and is what gives the adversarial robustness above.

What we *do* lose by not scrambling: diffusion as a debugging aid (FPSan's
payloads look random, so a single wrong bit anywhere produces a wildly different
payload that is easy to eyeball). $\varphi_p$ residues of related values are
related; small structured errors can produce small residue changes. For a
*detector* this is fine (we compare for equality); for *forensics* FPSan's
avalanche is nicer.

## 5. $+,-,\times,\div$ and casts

This is where $\varphi_p$ shines.

- **$+,-,\times$:** exact homomorphic images; every ring identity holds for free.
- **$\div$:** when the modulus is **prime** (variant 1b), $\mathbb{F}_p$ is a
  field, so division is **total except by residue $0$**, and $x/x = 1$ holds
  **always** (any $x$ with nonzero residue — i.e. any nonzero float, by §3a).
  Compare FPSan, whose $\mathbb{Z}/2^{w}$ has zero-divisors and only a
  *parity-preserving involution* for "inverse": $x/x=1$ there holds only for odd
  payloads. The field is a clean win. Division by residue $0$ occurs iff
  dividing by a genuine $\pm 0.0$ → honest infinity (§7). (The composite CRT
  variant gives this up for a small zero-divisor rate — §8.)
- **Casts / mixed precision:** all of `f16, bf16, f32, f64` embed in
  $\mathbb{Z}[1/2]$ via the *same* $\varphi_p$. A value-preserving cast
  (widening, or any exact narrowing) does not change the rational value, so it is
  **literally the identity on residues.** Mixed precision is unified: a residue
  doesn't know or care what width produced it. Contrast FPSan, which needs a
  per-width $\varphi_w$ and models casts as a non-value-preserving signed resize
  of the payload (see [understanding-fpsan.md](understanding-fpsan.md), iter. 3).
  Only *lossy* narrowing rounds, which — like all rounding — this model
  abstracts away (§9).
- **Subnormals:** also dyadic, handled by the same map with no special case
  (FPSan punts on them).

What is **lost**, and it is fundamental: **order, magnitude, `abs`,
`min`/`max`, comparisons.** $\mathbb{F}_p$ is not an ordered field; the residue
discards size and sign-as-order entirely. FPSan loses order too (its $\varphi$
is a hash), so this is not a differentiator, but here it is *structural*: there
is no ordered-field quotient to hope for. A reduction like `max(x₁,…,xₙ)` is
simply **not representable** (it needs the total order); FPSan "solves" min/max
by imposing an order on the *scrambled payload* — associative and
reassociation-invariant, but value-*infaithful* (it may select the wrong
element). For us, min/max/`abs`/comparisons are **out of scope**: they need a
side channel (carry the actual float, or tag).

**Cost.** Mod-$n$ arithmetic is costlier than FPSan's free mod-$2^{w}$
wraparound (a reduction per op). With a 32-bit modulus the reduction is a single
multiply-high (Barrett/Montgomery); a pseudo-Mersenne $p$ makes it a couple of
adds. We stay 32-bit on purpose (§8): footprint fidelity matters more than the
per-op cost for a sanitizer.

## 6. Transcendentals: exp is recoverable (CRT); the rest stay free

`exp, log, sin, cos, sqrt, …` are not rational functions, so the homomorphism
says nothing directly. The question is whether we can *build* an `exp` with the
right identity $\exp(a{+}b)=\exp(a)\exp(b)$, the way FPSan builds $g^{p}$.

**As a function of the mod-$p$ residue, no.** Such an `exp` would be a
homomorphism $(\mathbb{F}_p,+)\to(\mathbb{F}_p^\times,\times)$; since
$\gcd(p,\,p-1)=1$ the only one is trivial. So **variant 1b** (single prime)
cannot honor the exponential law — its `exp` is just a tagged hash token, and
the prototype duly finds the law holding $0/50000$ times (i.e. never). FPSan
escapes this only because $\mathbb{Z}/2^{w}$'s additive and unit groups share
the prime $2$.

**As a function of a second, mod-$d$ residue, yes.** The value-level map
$\exp:(\mathbb{Z}[1/2],+)\to(\mathbb{F}_p^\times,\times)$ *does* exist: pick $g$
of **odd order** $d$, set $\exp(1)=g$; the odd-order $2^{n}$-th roots are unique,
giving the closed form $\exp(v)=g^{\,v\bmod d}$, which satisfies the law exactly.
The catch: it depends on $v\bmod d$, **not** $v\bmod p$ — it does not factor
through the mod-$p$ residue (else $h(pw)=h(w)^{p}=h(w)$ forces triviality). So we
must *carry* $v\bmod d$.

**CRT keeps this compact (variant CRT).** Carrying $(v\bmod p,\,v\bmod d)$ is, by
CRT, carrying $v\bmod n$ with $n=pd$ — one residue, one 32-bit word. With
$p\equiv1\pmod d$ and $g$ of order $d$ in $(\mathbb{Z}/n)^\times$, the exponent
$v\bmod d = r\bmod d$ is read straight off the single residue, and
$\exp(v)=g^{\,r\bmod d}$. The prototype's CRT variant honors
$\exp(a{+}b)=\exp(a)\exp(b)$ **exactly** ($50000/50000$). So honoring exp costs
neither footprint nor a second word — it costs (see §8): a *composite* modulus,
hence zero-divisors at rate $\sim 1/p+1/d$ (observed $19/300\mathrm{k}$), and an
`exp` whose image is only the order-$d$ subgroup, so exp results collide at
$\sim 1/d$ (observed: 300k inputs → only $\approx d$ distinct exp outputs).

`sin`/`cos` get the analogous order-$d$ rotation (the Euler construction).
`log`, `sqrt`, `erf`, … have no usable identity — a faithful `log` would be
$\exp$'s inverse, the discrete logarithm, which is expensive and again not a
function of the residue — so they stay **hash tokens**. And that is *correct*,
not a fallback: by Lindemann–Weierstrass / Schanuel, transcendental values at
distinct algebraic arguments are **algebraically independent**, so a fresh free
generator per call is the faithful model.

| fragment | true-value structure | model here |
|---|---|---|
| algebraic ($+,-,\times,\div$) | richly related ($2{+}2=4$) | **value**: $\varphi_n$ residue |
| $\exp$ (and $\sin,\cos$) | one relation: $\exp(a{+}b)=\exp(a)\exp(b)$ | **structured**: $g^{\,v\bmod d}$ (CRT only) |
| $\log,\sqrt,\dots$ | algebraically independent | **free**: hash token |

## 7. Zero, infinity, NaN, signs

Finite floats are in $\mathbb{Z}[1/2]$; $\pm\infty$ and NaN are not. Extend the
codomain to the projective line plus an absorbing symbol:

$$\widehat{\mathbb{F}_p} \;=\; \mathbb{P}^1(\mathbb{F}_p)\ \sqcup\ \{\mathrm{NaN}\}
\;=\; \mathbb{F}_p \cup \{\infty\} \cup \{\mathrm{NaN}\}.$$

- finite nonzero $\mapsto$ its residue in $\mathbb{F}_p^\times$;  $\pm 0.0 \mapsto 0$;
  $\pm\infty \mapsto \infty$;  NaN $\mapsto$ NaN.
- $\mathbb{P}^1(\mathbb{F}_p)$ gives $x/0=\infty$ for $x\neq 0$ for free, and the
  Möbius/IEEE-shaped rules $\infty+\text{finite}=\infty$, $\infty\cdot x=\infty$
  ($x\neq0$). NaN is absorbing, and seeds the indeterminate forms
  $\infty+\infty,\ \infty-\infty,\ 0\cdot\infty,\ 0/0,\ \infty/\infty \mapsto
  \mathrm{NaN}$. With *unsigned* $\infty$, even $\infty+\infty$ is indeterminate
  — IEEE resolves it via the sign we dropped — so the model emits NaN slightly
  more freely than IEEE, but only at the already-non-finite fringe.

The §3a fact ($0$ is hit only by true zero) is what makes this principled:
$1/x=\infty$ exactly when $x$ is a genuine zero, never because some nonzero float
collided onto $0$. So $0$ and $\infty$ are honest. This is a **better** Inf/NaN
story than FPSan, which models none of it.

Not modeled (out of scope — these compact variants carry no sign bits):

- **Signed zero / signed infinity.** $+0.0,-0.0\mapsto 0$ and $\pm\infty\mapsto$
  one $\infty$; $\mathbb{F}_p$ has one zero, $\mathbb{P}^1$ one pole. (Negation
  $\varphi_n(-v)=-\varphi_n(v)$ is exact on *finite* values — only the $\pm0$ and
  $\pm\infty$ poles merge.)
- **(CRT variant)** division by a nonzero *zero-divisor* poisons to NaN, not
  $\infty$ — a units issue (§8), separate from the honest $1/0=\infty$.

## 8. The two compact variants

Everything above instantiates as one of two 32-bit designs. Both carry a single
residue plus the projective+NaN extension of §7; they differ only in the modulus.

**Variant 1b — single prime.** $n=p$ just below $2^{32}$ (the prototype uses
$2^{32}-5$). $\mathbb{Z}/n=\mathbb{F}_p$ is a **field**: division total, $x/x=1$
always, no zero-divisors; cleanest algebra and best collision margin
($\sim2^{-32}$). No exponential law — `exp` and all transcendentals are hash
tokens.

**Variant CRT — composite $n=pd$.** $p,d$ distinct odd primes, $p\equiv1\pmod d$,
$n<2^{32}$ (prototype: a Sophie-Germain pair near $2^{15}$, $n\approx2^{31}$).
$\mathbb{Z}/n\cong\mathbb{F}_p\times\mathbb{F}_d$ — a product of fields, *not* a
field. Buys the **exact exp homomorphism** $\exp(v)=g^{\,v\bmod d}$ (§6). Costs:

- **zero-divisors** at rate $\sim 1/p+1/d$: division by one is undefined and
  poisons to NaN ($x/x\ne1$ there). Observed $19/300\mathrm{k}$.
- **exp collisions** at $\sim 1/d$: `exp`'s image is only the order-$d$ subgroup.
  Observed: 300k inputs → $\approx d$ distinct exp outputs.
- a slightly smaller algebraic modulus → marginally more leaf collisions.

Prefer $n=pd$ (two distinct primes — a product of fields) over $\mathbb{Z}/p^2$:
same exp trick (via the order-$p$ principal units) but with **nilpotents**, which
are algebraically uglier and less value-faithful.

| | variant 1b (prime) | variant CRT ($pd$) |
|---|---|---|
| modulus | $p\approx2^{32}$ | $n=pd\approx2^{31}$ |
| algebra | field; division total, $x/x=1$ | product of fields; zero-divisors $\sim1/p{+}1/d$ |
| $\exp(a{+}b)=\exp(a)\exp(b)$ | ✗ (hash token) | ✓ ($g^{\,v\bmod d}$) |
| exp collisions | — | $\sim 1/d$ |
| leaf collisions | $\sim2^{-32}$ | $\sim2^{-31}$ |

Both share the **hybrid transcendental scheme**: carry $\varphi_n(\text{exact
value})$ for $+,-,\times,\div$ and casts (exact, value-faithful,
precision-agnostic, honest $0/\infty/$NaN); at a transcendental call, either
apply the structured map ($\exp$, and $\sin/\cos$, in the CRT variant) or mint a
fresh free generator $t=H(f,r_x)$ (`log`, `sqrt`, … — and *all* transcendentals
in 1b). Faithful by §6.

A prototype of both — [`algebraic_fpsan_generic.hpp`](algebraic_fpsan_generic.hpp)
and [`algebraic_demo.cpp`](algebraic_demo.cpp) — passes the ring-homomorphism,
ring-law, Inf/NaN, and (CRT) exp-homomorphism checks and reports the
collision / zero-divisor / exp-image statistics quoted above.

### Resolved and remaining

- **$\exp(a{+}b)=\exp(a)\exp(b)$ — resolved.** The worry that the prime field
  forbids it was about exp *as a function of the mod-$p$ residue*; carrying the
  mod-$d$ residue (folded in by CRT) delivers it exactly, at the costs tabled
  above. Extra collision margin, if needed, comes from **repeating with different
  prime pairs** (time, not space — §3).
- **Order / `min`/`max` / `abs` / comparisons — open, and fundamental.**
  $\mathbb{F}_p$ is unordered; no sign or exponent trick recovers a total order.
  This is the real limitation of the approach, shared with FPSan (which fakes it
  with payload-order, sacrificing value-faithfulness). A side channel (the actual
  float, or a tag) is the only route.
- **Rounding — out of scope by design (§9).**

## 9. A note on rounding (shared boundary with FPSan)

Both systems carry *exact* model values and never round, so both answer "are
these the same up to exact algebra?" and are deliberately **invariant to
reassociation and to rounding-order**. The sanitizer flags rewrites that change
the *exact* value, not ones that merely reround. This is a feature (it is what
lets a reduction in two orders compare equal), and it is the same boundary FPSan
draws. $\varphi_n$ does not change it.

A corollary worth stating, since it surprises: **overflow to $\infty$ is not
modeled.** The model is exact and carries no magnitude ($\mathbb{F}_p$ is
unordered), so a sum of finite values is always a finite residue — never
$\infty$. The *only* source of $\infty$ is division by a true zero. Thus
`1/(sum of huge values)` is an ordinary finite inverse here, where IEEE would
overflow to $\infty$ and then give $0$; both behaviors are deliberate
consequences of abstracting rounding and magnitude away.

## 10. Scorecard vs. FPSan

| property | FPSan ($\mathbb{Z}/2^{w}$, scrambled) | variant 1b (prime) | variant CRT ($pd$) |
|---|---|---|---|
| encoding is a ring homomorphism | no (free model) | **yes (value)** | **yes (value)** |
| respects all rational-function identities | no (only axioms) | **yes** | **yes** |
| $x/x=1$, division total | partial (odd only) | **yes (field)** | mostly (zero-div $\sim1/p{+}1/d\to$NaN) |
| mixed precision / casts | per-width $\varphi_w$, resize | **one map, exact casts = id** | **one map, exact casts = id** |
| subnormals | punted | **uniform** | **uniform** |
| Inf / NaN | not modeled | **$\mathbb{P}^1\sqcup\{$NaN$\}$, honest $0$** | **$\mathbb{P}^1\sqcup\{$NaN$\}$, honest $0$** |
| $\exp(a{+}b)=\exp(a)\exp(b)$ | yes ($g^{p}$) | no (hash token) | **yes ($g^{\,v\bmod d}$)** |
| other transcendentals | tagged scramble (free) | hash token (free) | hash token (free) |
| order / `abs` / `min` / `max` | faked (payload-order, infaithful) | out of scope | out of scope |
| footprint | 32-bit | 32-bit | 32-bit |
| leaf-collision rate | $\sim2^{-w}$ | $\sim2^{-32}$ | $\sim2^{-31}$ |
| exp-collision rate | $\sim2^{-w}$ | — (no exp) | $\sim1/d\;(\approx2^{-15})$ |
| randomization | fixed scramble (fixed bad pairs) | **random prime** | **random prime pair** |
| per-op cost | free wraparound | mod-$p$ reduction | mod-$n$ reduction |

**Bottom line.** On the *algebraic* fragment both variants are a strictly better
foundation than FPSan — the honest value model, clean division and mixed
precision, a principled Inf/NaN story — at the cost of mod-$n$ arithmetic and a
few bits of collision margin (bought back by repeating with another prime).
**Variant 1b** is the clean field: pick it when the exponential law isn't needed.
**Variant CRT** adds the exact exp homomorphism (and $\sin/\cos$), paying with
zero-divisors and a small exp image. The genuine remaining limitation, shared
with FPSan and not fixable by any sign/exponent trick, is **order / `min` /
`max`** — a finite field has no compatible order.
