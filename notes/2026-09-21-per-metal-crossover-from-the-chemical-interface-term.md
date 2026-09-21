# A per-metal p/n crossover from the chemical interface term Δ_c(d)

**Date:** 2026-09-21
**Open item closed:** "A per-metal crossover `w_cross = W_G + D_c(d)`" — created
2026-09-20 and listed **top** of that session's "not yet covered" list.

> **This file is committed in two parts, deliberately.** Sections 1–4 —
> including the four numbered predictions in Section 4 — were committed
> *before* the model in `graphene_per_metal_crossover_model.py` was written or
> run. Section 5 was appended afterwards. The commit history is the evidence
> that the predictions were not written to match the answer. 2026-09-20's
> methodological note observed that this repo's claims have twice failed on
> widened samples rather than on physics; pre-registration is the cheap half
> of the fix, and enumerating the whole space (Section 4, P0) is the other.

---

## 1. What this repo currently assumes, and why it is an assumption

Since 2026-09-18 every signed contact-field calculation in this repo has used

```python
W_CROSS_CHEM = 5.4            # eV
signed_offset(W_metal) = W_metal - W_CROSS_CHEM
```

a **single crossover work function shared by all seven metals** in
`METAL_WORK_FUNCTIONS`. The 5.4 eV is real and well-sourced — Giovannetti,
Khomyakov, Brocks, Karpan, van den Brink & Kelly, *Phys. Rev. Lett.* **101**,
026803 (2008) — but it is a *derived* number, not a constant of nature.

Khomyakov, Giovannetti, Rusu, Brocks, van den Brink & Kelly, *Phys. Rev. B*
**79**, 195425 (2009) ([arXiv:0902.1203](https://arxiv.org/abs/0902.1203)) —
the companion paper this repo already leans on for the nonlinear
`dW → dE_F` relation (2026-09-20) — writes the interface potential step as the
sum of a charge-transfer term and a **short-range chemical term**,

```
    ΔV(d) = Δ_tr(d) + Δ_c(d)                                        (their §III)
    Δ_c(d) = e^(−γ d) (a_0 + a_1 d + a_2 d²)                        (their Eq. 4)
```

Charge neutrality requires the *total* offset `dW' = W_M − W_G − Δ_c(d)` to
vanish, so the crossover is

```
    W_0(d) = W_G + Δ_c(d)                                           (*)
```

**a function of the metal–graphene separation `d`, not a constant.** The 5.4 eV
figure is (*) evaluated at the physisorbed equilibrium separation
`d ≈ 3.3 Å`, where `Δ_c ≈ 0.9 eV`. Metals that sit at a different separation
have a different crossover, and therefore a different `dW` — which moves
*every* number downstream of `signed_offset()`.

That is the difference between this open item and 2026-09-20's. The nonlinear
`dE_F(dW')` relation **compressed** existing offsets by a few percent without
moving their zero. A per-metal crossover **moves the zero itself**, separately
for each metal, on the same body of evidence.

## 2. The parameters that are, and are not, available

`Δ_c` is the *only* new ingredient (*) needs. Extracting its four fitted
constants from the PRB was attempted and **failed**:

* `ar5iv.labs.arxiv.org/html/0902.1203` renders Eq. 4 correctly in symbolic
  form, and two separate fetches with differently-worded extraction prompts
  both returned the equation **without** numerical values for `γ`, `a_0`,
  `a_1`, `a_2`. The rendered text passes straight from "we parameterize it"
  to the physical discussion.
* This is the **second** extraction failure on this same paper. 2026-09-20
  recorded two passes over the PDF returning mutually contradictory Table I
  contents. The failure mode is different (missing here, contradictory there)
  but the lesson is the same, and the response is the same: **no fitted
  coefficient from this paper is used anywhere in the model.**

What *is* available, and was cross-confirmed by both of 2026-09-20's passes
and by the ar5iv text today:

| quantity | value | role |
|---|---|---|
| `Δ_c(3.3 Å)` | ≈ 0.9 eV | the anchor |
| `W_G` | 4.5 eV | graphene's own work function |
| `W_0(3.3 Å)` | ≈ 5.4 eV | consistency check on the anchor: 4.5 + 0.9 |
| `d_eq` per metal | 2.05–3.41 Å | already tabulated in `D_EQ` (2026-09-20) |
| `d_0` | 2.4 Å | gap-capacitance cutoff (2026-09-20, unrelated to Δ_c) |

## 3. The model actually used: one unknown, anchored and swept

Rather than guess `(γ, a_0, a_1, a_2)`, `Δ_c` is written as a **one-parameter
family pinned to the anchor**:

```
    Δ_c(d ; ℓ) = 0.9 eV · exp( −(d − 3.3 Å) / ℓ )                    (N1)
    w_cross(metal ; ℓ) = W_G + Δ_c(d_eq(metal) ; ℓ)                  (N2)
```

Three properties make this the honest choice:

1. **It reproduces the repo's current model exactly in a limit.** As
   `ℓ → ∞`, (N1) returns 0.9 eV for every `d`, so (N2) returns 5.4 eV for
   every metal and the whole per-metal apparatus collapses **bitwise** onto
   `W_CROSS_CHEM`. There is an exact reduction test, not a plausible-range
   one — the distinction that caught a real bug on 2026-09-17 and again on
   2026-09-19.
2. **It is exact at the one datum that is trusted.** `Δ_c(3.3 Å) = 0.9 eV`
   by construction, bitwise.
3. **The single unknown has a clear meaning and a bounded range.** `ℓ` is the
   decay length of an interface term produced by wavefunction overlap and
   Pauli repulsion. Such terms fall off on the scale of the orbital tails
   involved: **ℓ ∈ [0.3 Å, 1.5 Å]** brackets that generously at both ends,
   and every conclusion below is reported as a sweep over it.

(N1) is a *simplification* of Eq. 4, not Eq. 4. Dropping the polynomial
prefactor is a real approximation and Section 5 reports where it breaks —
which turns out to be the most interesting result of the session.

## 4. Pre-registered predictions

**P0 (method, not physics).** Following 2026-09-20's lesson — "enumerate the
whole space instead of tabulating two representative cases" — every claim
below is to be checked on **all** metals and **all** ordered pairs that the
model admits, not a chosen subset.

**P1 — no doping-sign changes.** No metal in `METAL_WORK_FUNCTIONS` changes
from n-doping to p-doping or back, for any `ℓ ∈ [0.3, 1.5] Å`.

**P2 — small individual shifts.** For the metals the model admits, `|dW|`
changes by **≤ 10%** relative to the flat 5.4 eV convention, across the whole
`ℓ` range.

**P3 — amplification on same-sign pairs.** For a pair whose two contacts sit
on the *same* side of the crossover, the fractional change in the pair's net
response **exceeds** the fractional change in either contact's own `dW`. This
is the mechanism 2026-09-20 identified — a same-sign pair is a near
cancellation, so its surviving magnitude is set by the *ratio* of the two
offsets, and perturbing a ratio is a larger relative perturbation than
perturbing either term. Straddling pairs should show no such amplification.

**P4 — the chemisorbed metals are out of regime, again.** Extrapolating (N1)
inward to `d_eq ≈ 2.05–2.3 Å` (Ti, Ni, Pd) should give crossover work
functions **above the work function of every elemental metal**, which would
make all three maximally n-doping irrespective of their actual work function.
If so, (N1) must be refused there — and for a reason entirely independent of
2026-09-20's, which refused the same three metals because `d_eq < d_0` put
them outside the *gap-capacitance* term. Two different parts of the same
framework failing on the same three metals, for different reasons, would be a
structural statement about Chapter 6 rather than a numerical one.

Sections 1–4 end here. Everything below was written after the model ran.

---

## 5. Outcome — three of four predictions held, and the one that failed is the result

Everything in this section was written after `graphene_per_metal_crossover_model.py`
ran. Full console output is reproducible with `python3 graphene_per_metal_crossover_model.py`.

### 5.0 Validations first

Four exact checks, no plausible-range checks anywhere:

| check | result |
|---|---|
| V1 — the dW-level entry point equals the 2026-09-18 signed model | **49/49 ordered pairs bitwise**, worst difference 0.000e+00 |
| V2 — `ℓ → ∞` collapses onto the flat 5.4 eV convention | `Δ_c(d, ∞) == 0.9` bitwise for every tabulated `d`; `w_cross == W_G + 0.9` bitwise for every metal; **49/49 pair responses bitwise**; and `4.5 + 0.9` turns out to be the *same double* as the literal `5.4` (gap 0.000e+00), which was not assumed |
| V3 — the anchor survives the parametrisation | `Δ_c(3.3 Å, ℓ) == 0.9` bitwise for **61/61** values of `ℓ` |
| V4 — exact symmetries survive a per-metal crossover | symmetric pair `N = 0` to 6.6e-17 over 9 (metal, ℓ) cases; charge conjugation **exactly 0.000e+00** over 27 (pair, ℓ) cases |

V1 is the one that licenses the rest. `net_response()` takes a single
`w_cross` and forms `dW` internally, so it structurally *cannot* express a
per-metal crossover; a new entry point one level lower was unavoidable. Being
bitwise-identical on all 49 pairs means it is not "a similar model", it is the
same arithmetic with the subtraction moved outward.

### 5.1 Only three of seven metals can be evaluated at all

| metal | `W` (eV) | `d_eq` (Å) | status |
|---|---|---|---|
| Ti | 4.33 | 2.10 | extrapolated — chemisorbed |
| Cr | 4.50 | — | **no tabulated separation**, not guessed |
| Cu | 4.65 | 3.26 | ok |
| Ni | 5.04 | 2.05 | extrapolated — chemisorbed |
| Au | 5.10 | 3.31 | ok |
| Pd | 5.12 | 2.30 | extrapolated — chemisorbed |
| Pt | 5.65 | 3.30 | ok |

This is the same 3-of-7 that 2026-09-20 was left with, and it is worth being
precise about *why the two restrictions coincide*, because they are not the
same restriction. 2026-09-20 refused Ti/Ni/Pd because `d_eq < d_0` puts them
outside the **gap-capacitance** term of Eq. 7. Today refuses them because
`Δ_c` extrapolated inward diverges (§5.3). Cr is refused by both, for the
third and most boring reason: nobody tabulated its separation.

### 5.2 P1 held, P2 **FAILED**

**P1 (no doping-sign changes): PASS.** Cu and Au stay n-doping and Pt stays
p-doping for every `ℓ` in [0.3, 1.5] Å. The crossover moves, but never past a
metal.

**P2 (|dW| shifts ≤ 10%): FAIL — worst case 17.12%, on Cu.**

| metal | `dW` flat | `dW` range over `ℓ` | max shift |
|---|---|---|---|
| Cu | −0.7500 | −0.8784 … −0.7743 | **17.12 %** |
| Au | −0.3000 | −0.2940 … −0.2705 | 9.84 % |
| Pt | +0.2500 | +0.2500 … +0.2500 | 0.00 % |

Three things about this failure, in descending order of how much they matter.

**(a) The size is the point.** 17% is **more than three times** the ~4.9%
compression that 2026-09-20's nonlinear `dW → dE_F` relation produced on the
Ti/Pt headline — and that session treated ~5% as the reassuring outcome. The
prediction of "≤ 10%" was calibrated on that precedent and it was the wrong
precedent: moving the zero of `dW` is a larger perturbation than compressing
`dW` about a fixed zero, because the shift does not scale with `|dW|`. A metal
close to the crossover has a small `|dW|` and therefore feels a fixed shift
*most*, which is the opposite of the intuition that "small offsets are robust".

**(b) The failure is driven by the short end of the swept range, and the
threshold is sharp.** Bisection (`p2_threshold()`, exact to machine precision
because the shift is monotone in `ℓ`) puts the boundary at **ℓ = 0.50 Å**:
P2 holds for every `ℓ ≥ 0.50 Å` and fails below. The midpoint ℓ = 0.6 Å gives
8.3% on Cu, inside the prediction. So P2 is not comprehensively wrong; it is
wrong over roughly the shortest sixth of the range that was declared
admissible — and declaring that range generously, in Section 3, is what made
the failure visible instead of invisible. A narrower, more "reasonable" sweep
would have returned PASS and been worth less.

**(c) Pt's exactly-0.00% is an artefact and is not evidence of anything.**
`d_eq(Pt) = 3.30 Å` coincides with the anchor `D_ANCHOR = 3.3 Å`, so
`Δ_c(d_eq(Pt)) = 0.9 eV` bitwise for every `ℓ` by construction. Pt is pinned,
not robust. Had the anchor been placed at Au's 3.31 Å instead, Au would be
the pinned one and Pt would move. **Every per-metal number in this section is
therefore a shift *relative to Pt*, not an absolute one.** That is a real
limitation of a one-anchor parametrisation and it is stated here rather than
left for a reader to notice.

### 5.3 P4 held, and it is a reductio that earns its keep

Extrapolating (N1) inward to the chemisorbed separations gives crossover work
functions of **6.25 – 62.6 eV**:

| metal | `d_eq` (Å) | ℓ = 0.3 Å | ℓ = 0.6 Å | ℓ = 1.0 Å | ℓ = 1.5 Å |
|---|---|---|---|---|---|
| Ni | 2.05 | 62.55 | 11.73 | 7.64 | 6.57 |
| Ti | 2.10 | 53.64 | 11.15 | 7.49 | 6.50 |
| Pd | 2.30 | 29.73 | 9.27 | 6.95 | 6.25 |

Every one of these is above **5.9 eV**, the highest elemental work function
there is. The implied `|dW|` runs 1.13 – 57.5 eV against a largest *observed*
Fermi-level shift of **0.5 eV**. So the extrapolation does not merely lose
accuracy — it predicts that every chemisorbed metal is maximally n-doping
regardless of its own work function, which is false, and it does so
**across the entire admissible range of ℓ**, not just at its short end.

The honest reading: this is a failure of **(N1)**, not of Khomyakov *et al.*
Their Eq. 4 carries a polynomial prefactor `(a_0 + a_1 d + a_2 d²)`
multiplying the exponential, and a polynomial is exactly what lets `Δ_c` turn
over instead of running away at short `d`. §2 could not recover its
coefficients, so §3 dropped it, and **RESULT 2 is the bill for that
simplification, arriving exactly where it was predicted to**. What makes it
worth committing rather than merely conceding is that it independently
reaches 2026-09-20's conclusion by a different route: two distinct pieces of
this framework — the gap-capacitance term and the chemical term — fail on the
*same three metals* for *unrelated reasons*. That is a structural statement
about which metals graphene's contact physics can be described by work
function at all, and Ti is in the failing set while being contact A of both
Chapter 6 headline pairs.

### 5.4 P3 held, with an amplification factor of 4.6

| pair | kind | `N` flat | `N` range over `ℓ` | pair change | `dW` change | amplification |
|---|---|---|---|---|---|---|
| Cu/Au | same-sign | −0.6549 | −1.1653 … −0.6922 | **77.94 %** | 17.12 % | **4.55** |
| Au/Cu | same-sign | +0.6549 | +0.6922 … +1.1653 | **77.94 %** | 17.12 % | **4.55** |
| Cu/Pt | straddling | −1.7855 | −1.8072 … −1.7900 | 1.21 % | 17.12 % | 0.07 |
| Pt/Cu | straddling | +1.7855 | +1.7900 … +1.8072 | 1.21 % | 17.12 % | 0.07 |
| Au/Pt | straddling | −1.6402 | −1.6368 … −1.6228 | 1.06 % | 9.84 % | 0.11 |
| Pt/Au | straddling | +1.6402 | +1.6228 … +1.6368 | 1.06 % | 9.84 % | 0.11 |

All six ordered asymmetric pairs the model admits, enumerated, per P0 — which
here is cheap to the point of triviality, since three usable metals leave only
six.

The separation is stark and is the cleanest confirmation the
ratio-amplification mechanism has had: **a 17% shift in one contact's `dW`
becomes a 78% shift in the same-sign pair's response (×4.6), and a 1.2% shift
in the straddling pair's (×0.07)** — a factor of **65 between the two kinds
of pair**, from identical inputs. 2026-09-20 inferred this mechanism from a
single same-sign pair (Ti/Pd losing 56%) under a *different* perturbation;
finding it again, with the same sign and a comparable magnitude, under a
perturbation that moves the zero rather than compressing the scale, is
independent evidence that it is a property of near-cancelling pairs and not of
either particular model.

The straddling pairs are correspondingly *more* robust than P2's per-metal
number suggested: Cu/Pt moves 1.2% while its own Cu contact moves 17%. A
near-cancellation amplifies; a reinforcing pair averages.

### 5.5 What this does and does not change for Chapter 6

* **It does not touch the Ti/Pt headline**, because Ti cannot be evaluated.
  The chapter's central number is now untested by *two* successive
  refinements rather than confirmed by them, which is a weaker position than
  2026-09-20 left it in, not a stronger one.
* **The best evaluable asymmetric pair, Cu/Pt (−1.786), is robust to 1.2%**
  across the whole `ℓ` range. Chapter 6's qualitative claim — a straddling
  n/p pair is the design that works — survives cleanly.
* **Cu/Au should not be quoted to better than a factor of two** under any
  crossover convention. It is the only same-sign pair among the usable metals
  and it moves 78%.
* **`W_CROSS_CHEM = 5.4` remains the right default** for physisorbed metals,
  now with a quantified error bar rather than an implicit claim of exactness:
  ±17% on `dW` for a metal 0.04 Å off the anchor, ±78% on a same-sign pair
  built from such metals.

**Not superseded, annotated.** No previous number is deleted. Section 6.9.5's
and 6.11's tables stand as committed; Section 6.12 states what this session
found and which of their entries it does and does not reach.

## 6. Honest record of what did not work

* **Eq. 4's fitted coefficients could not be extracted** from ar5iv in two
  attempts with differently-worded prompts. They are not used, not guessed,
  and the cost of not having them is quantified in §5.3 rather than hidden.
* **No live-search route to the coefficients was found** that did not go
  through the publisher paywall or ResearchGate, which has been rate-limiting
  this repo since 2026-08-31.
* **PubMed/PMC were not attempted**, per the standing note that they return
  reCAPTCHA in this environment (four occurrences logged).
