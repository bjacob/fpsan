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

Choose $p$ slightly below $2^{32}$ (or see §5 on $2^{61}-1$).

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
$\varphi_p(m\cdot 2^e)=0 \iff p\mid m \iff m=0$, because $|m|<2^{24} < p$. So,
provided $p>2^{24}$, the *only* finite float with residue $0$ is $\pm 0.0$. This
is the crucial design constraint: it keeps $0$ honest, which is what makes
division and the zero/infinity story (§7) clean. Nothing nonzero ever silently
becomes zero.

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

- Per comparison, false-match $\approx 2^{-24}$ with a 32-bit prime — a few bits
  weaker than FPSan's $\approx 2^{-32}$, but **tunable**: a 64-bit prime gives
  $\approx 2^{-53}$; $k$ independent primes (CRT) multiply the rates.
- Unlike FPSan's *fixed* scramble — which has a fixed set of colliding pairs that
  collide on every run — randomizing $p$ makes collisions **non-repeatable and
  adversarially robust**: there is no bad pair that always fools you. This is a
  real qualitative gain in the randomized regime.

Verdict: the loss of injectivity is essentially cosmetic for the intended use.
The substantive cost is a slightly higher per-prime collision rate, bought back
with a bigger or repeated prime.

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
- **$\div$:** $\mathbb{F}_p$ is a field, so division is **total except by
  residue $0$**, and $x/x = 1$ holds **always** (any $x$ with nonzero residue —
  i.e. any nonzero float, by §3a). Compare FPSan, whose $\mathbb{Z}/2^{w}$ has
  zero-divisors and only a *parity-preserving involution* for "inverse":
  $x/x=1$ there holds only for odd payloads. The field is a clean win.
  Division by residue $0$ occurs iff dividing by a genuine $\pm 0.0$ → honest
  infinity (§7).
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
is no ordered-field quotient to hope for. Comparisons need a side channel (carry
the actual float alongside, or tag) — same as FPSan.

**Cost.** Mod-$p$ arithmetic is costlier than FPSan's free mod-$2^{w}$
wraparound (a reduction per op). Mitigation: a pseudo-Mersenne or Mersenne prime
for cheap reduction. $p = 2^{61}-1$ is a sweet spot — prime, fits in 64 bits,
trivial reduction, collision $\approx 2^{-53}$, and $2$ has large order so the
$\bar2^{\,e}$ are well spread. (A *fixed* prime forfeits the adversarial
robustness of §3c; use a random prime when that matters.)

## 6. Transcendentals: where the prime field bites

`exp, log, sin, cos, sqrt, …` are not rational functions, so the homomorphism
says nothing: $\varphi_p(\exp x)$ has no relation to $\exp(\varphi_p x)$
($\exp$ isn't even defined on $\mathbb{F}_p$). The question is whether we can
*build* an `exp` on residues with the right identity, the way FPSan builds
$\exp(p)=g^{p}$ to get $\exp(a{+}b)=\exp(a)\exp(b)$.

**In a prime field, we cannot.** Such an `exp` would be a homomorphism
$(\mathbb{F}_p,+)\to(\mathbb{F}_p^\times,\times)$, i.e. from a group of order
$p$ to one of order $p-1$. Since $\gcd(p,\,p-1)=1$, the only such homomorphism is
**trivial**. So the single most useful transcendental identity — the one FPSan
gets — is *unavailable* here. FPSan's trick works precisely because it lives in
$\mathbb{Z}/2^{w}$, where the additive group (order $2^{w}$) and the unit group
(order $2^{w-1}$) share the prime $2$, so $g^{p}$ is a nontrivial homomorphism.
The prime field's coprimality is exactly what forbids it.

Partial escapes:

- **$\mathbb{Z}/p^2$** (still a homomorphic image of $\mathbb{Z}[1/2]$): the
  principal units $1+p\mathbb{Z}/p^2$ form a cyclic group of order $p$, so
  $\exp(r)=(1+p)^{r}$ *is* a nontrivial additive→multiplicative homomorphism.
  But it has period $p$ (kernel of order $p$), so it is wildly non-injective, and
  $\mathbb{Z}/p^2$ is no longer a field (we lose §5's clean division). A real
  trade.
- More generally, only $\mathbb{Z}/n\mathbb{Z}$ with $\gcd(n,\varphi_{\text{Euler}}(n))>1$
  supports a nontrivial exp; primes are the worst case for this and the best case
  for division. Tension.

Now the **conceptual payoff**, which I think is the most interesting point in
these notes. By Lindemann–Weierstrass / Schanuel, the values of $\exp$ (and the
other transcendentals) at distinct algebraic arguments are **algebraically
independent**. So the *faithful* value model for the transcendental fragment is
the one in which each transcendental result is a **fresh, free generator** —
which is exactly the *free* model, i.e. exactly what FPSan's scramble produces!

So the dichotomy of §2 resolves with a clean division of labour:

| fragment | true values | faithful model |
|---|---|---|
| algebraic ($+,-,\times,\div$, rational fns) | richly related ($2+2=4$) | **value** ($\varphi_p$) |
| transcendental ($\exp,\log,\sin,\dots$) | algebraically independent | **free** (FPSan scramble) |

FPSan applies the free model everywhere (paying for it on the algebraic part by
over-discriminating: $2{+}2\neq4$). The pure-$\varphi_p$ system would apply the
value model everywhere (paying for it on the transcendental part by having no
honest construction at all, *worse* than FPSan there). The right object uses each
where it is faithful.

## 7. Zero, infinity, NaN, signs

Finite floats are in $\mathbb{Z}[1/2]$; $\pm\infty$ and NaN are not. Extend the
codomain to the projective line plus an absorbing symbol:

$$\widehat{\mathbb{F}_p} \;=\; \mathbb{P}^1(\mathbb{F}_p)\ \sqcup\ \{\mathrm{NaN}\}
\;=\; \mathbb{F}_p \cup \{\infty\} \cup \{\mathrm{NaN}\}.$$

- finite nonzero $\mapsto$ its residue in $\mathbb{F}_p^\times$;  $\pm 0.0 \mapsto 0$;
  $\pm\infty \mapsto \infty$;  NaN $\mapsto$ NaN.
- $\mathbb{P}^1(\mathbb{F}_p)$ gives $x/0=\infty$ for $x\neq 0$ for free, and the
  Möbius/IEEE-shaped rules $\infty+\text{finite}=\infty$, $\infty\cdot x=\infty$
  ($x\neq0$). NaN is absorbing, and seeds the genuinely-undefined forms
  $\infty-\infty,\ 0\cdot\infty,\ 0/0,\ \infty/\infty \mapsto \mathrm{NaN}$,
  matching IEEE's structure.

The §3a fact ($0$ is hit only by true zero) is what makes this principled:
$1/x=\infty$ exactly when $x$ is a genuine zero, never because some nonzero float
collided onto $0$. So $0$ and $\infty$ are honest. This is a **better** Inf/NaN
story than FPSan, which models none of it.

Caveats, both minor and both about signs the residue can't see:

- **Signed zero** is lost ($+0.0,-0.0\mapsto 0$); $\mathbb{F}_p$ has one zero.
- **Signed infinity** is lost ($\pm\infty\mapsto$ one $\infty$); $\mathbb{P}^1$
  has one point at infinity. (Negation $\varphi_p(-v)=-\varphi_p(v)$ works
  perfectly on *finite* values — sign of finite numbers is fine; it is only the
  $\pm0$ and $\pm\infty$ poles that merge.) If a use case needs them, carry a
  sign bit beside the residue.

## 8. Synthesis: the hybrid worth building

The pieces compose into one design:

> Carry a residue $r=\varphi_p(\text{exact value})\in\widehat{\mathbb{F}_p}$.
> Do $+,-,\times,\div$ and casts as field/projective operations — exact,
> total, value-faithful, precision-agnostic, with honest $0/\infty/$NaN.
> At each transcendental call $f(x)$, where no rational residue exists, mint a
> **fresh free generator**: a value $t$ determined deterministically by
> $(f, r_x)$ (a keyed scramble / lookup, as FPSan does for its tagged ops), and
> continue algebraically with $t$. This is faithful because transcendental
> outputs *are* free (§6).

This dominates both parents: the algebraic core gets $\varphi_p$'s superior
properties (field division, unified mixed precision, true value-equivalence,
clean Inf/NaN); the transcendental parts get FPSan-style freeness, which is the
correct model there.

Genuinely open design questions:

1. **Honoring $\exp(a{+}b)=\exp(a)\exp(b)$ in the hybrid.** A fresh generator
   for each $\exp$ call forgets this one true relation among $\exp$-values.
   Recovering it requires stepping outside $\mathbb{F}_p$ — e.g. a side channel
   in the order-$p$ part of $(\mathbb{Z}/p^2)^\times$ (§6), or carrying an extra
   "exponent" coordinate. Is there a clean structure that is a field for the
   algebraic part *and* carries the one exponential relation? (Schanuel says
   that relation, plus its angle-addition cousins for $\sin/\cos$, is essentially
   *all* the structure there is — so this is a bounded amount to capture.)
2. **Comparisons / order / `abs`.** Unavoidably a side channel; how cheaply?
3. **Rounding.** Like FPSan, the model abstracts rounding entirely (§9). Is there
   any value-faithful way to make rounding visible without leaving the algebraic
   world? (Likely no — rounding is not algebraic — but worth stating as the
   boundary of the approach.)

## 9. A note on rounding (shared boundary with FPSan)

Both systems carry *exact* model values and never round, so both answer "are
these the same up to exact algebra?" and are deliberately **invariant to
reassociation and to rounding-order**. The sanitizer flags rewrites that change
the *exact* value, not ones that merely reround. This is a feature (it is what
lets a reduction in two orders compare equal), and it is the same boundary FPSan
draws. $\varphi_p$ does not change it.

## 10. Scorecard vs. FPSan

| property | FPSan ($\mathbb{Z}/2^{w}$, scrambled) | $\varphi_p$ ($\mathbb{Z}[1/2]\to\mathbb{F}_p$) |
|---|---|---|
| encoding is a ring homomorphism | no (free model) | **yes (value model)** |
| respects all rational-function identities | no (only axioms) | **yes** |
| $x/x=1$, division total | partial (odd only) | **yes (field)** |
| mixed precision / casts | per-width $\varphi_w$, resize | **one map, exact casts = id** |
| subnormals | punted | **uniform** |
| Inf / NaN | not modeled | **$\mathbb{P}^1\sqcup\{$NaN$\}$, honest $0$** |
| $\exp(a{+}b)=\exp(a)\exp(b)$ | **yes ($g^{p}$)** | no (prime-field obstruction) |
| transcendentals | tagged scramble (faithful: free) | needs the same fallback |
| order / `abs` / compare | lost | lost (structural) |
| injectivity | bijection | $\sim$ injective; $0$ safe; tunable collisions |
| randomization | fixed scramble (fixed bad pairs) | **random prime (no fixed bad pairs)** |
| collision rate / comparison | $\sim 2^{-w}$ | $\sim 2^{-24}$ (32-bit) … $2^{-53}$ ($2^{61}{-}1$) |
| per-op cost | free wraparound | mod-$p$ reduction (Mersenne helps) |

**Bottom line.** The residue homomorphism is a strictly better foundation for
the *algebraic* fragment of float programs — it is the honest value model, it
makes division and mixed precision clean, and it gives a principled Inf/NaN
story — at the cost of mod-$p$ arithmetic and a few bits of collision margin
(recoverable). It is *not* a drop-in win, because the prime field structurally
forbids FPSan's one genuinely-nice transcendental trick. The synthesis in §8 —
$\varphi_p$ for algebra, free generators for transcendentals — is the design
these notes point to, and the open question in §8.1 is the crux of whether it
can also keep the exponential law.
