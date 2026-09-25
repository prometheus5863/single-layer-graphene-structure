# How far can the p/n crossover move before Chapter 6's conclusions change?
## A condition-number enumeration over all 21 asymmetric pairs

**Date:** 2026-09-22
**Status of this file at the moment it is committed:** the literature section
and the four predictions in Section 4 are written **before**
`graphene_crossover_sensitivity_model.py` exists. Git history is the evidence.
This is the second session to pre-register; 2026-09-21 established the
practice and the practice immediately earned its keep by falsifying its own
prediction at 17%.

---

## 1. The open item this addresses

2026-09-21's log put this at the **top** of "Not yet covered":

> **Propagate the per-metal crossover back through Sections 6.9-6.11.** Those
> sections' numbers are still at the flat 5.4 eV convention, and today's Cu/Au
> result says a same-sign pair under the flat convention can be wrong by 78%.
> The 2026-09-20 retraction table (21 asymmetric pairs) is the obvious first
> target, since its worst violators were same-sign pairs.

A direct propagation is **blocked**, and the blockage is worth stating before
anything else. The per-metal crossover of Section 6.12 is
`w_cross(metal) = W_G + Delta_c(d_eq(metal))`, and `d_eq` is available for
only three of this repo's seven metals: Cu, Au and Pt. Ti, Ni and Pd are
chemisorbed (`d_eq < d0`), where Section 6.12's Result 2 showed the anchored
exponential diverges; Cr has no tabulated separation at all. But **every one
of the five pairs that the 2026-09-20 retraction table makes its headline of**
-- Au/Pd, Ti/Cr, Ni/Au, Cr/Cu, Ni/Pd -- contains at least one metal outside
that set. The per-metal machinery cannot confirm or refute a single one of
them.

So the propagation that was asked for cannot be done with the tool that
prompted the asking. This session does the propagation that **can** be done,
and which reaches all seven metals: treat the crossover as a scalar with an
unknown offset, and measure how hard each pair amplifies that offset.

## 2. Why a scalar offset is the right perturbation, and how big it is

Both sources of the 5.4 eV number write it with a tilde.

- Giovannetti, Khomyakov, Brocks, Karpan, van den Brink, Kelly, *Doping
  Graphene with Metal Contacts*, **Phys. Rev. Lett. 101, 026803 (2008)**
  (arXiv:0802.2267): *"the crossover from p-type to n-type doping occurs for a
  metal work function of ~5.4 eV"*, against *"the graphene work function of
  4.5 eV"*. Fetched 2026-09-22 via arxiv.org/abs/0802.2267; the abstract gives
  the number with no error bar and no significant-figure claim.
- Khomyakov, Giovannetti, Rusu, Brocks, van den Brink, Kelly, **Phys. Rev. B
  79, 195425 (2009)** (arXiv:0902.1203): same sentence, *"a metal with a work
  function ~5.4 eV"*. Fetched 2026-09-22. The abstract states the physisorbed
  Fermi-level shift as *"up to 0.5 eV"* and says chemisorption *"reduces
  considerably"* the graphene work function **without giving a number** --
  the third distinct extraction failure on this pair of papers, after
  2026-09-20's contradictory Table I and 2026-09-21's twice-missing Eq. 4
  coefficients. Recorded as a failure, not worked around.

**This repo has been reading a tilde as three significant figures since
2026-09-18.** Every `w_cross=5.4` in `graphene_photodetector_signed_carrier_model.py`
and every table in Sections 6.8-6.11 inherits it.

A scalar offset `delta` on `w_cross` is exactly the right shape for that
uncertainty, for three reasons:

1. It is what a tilde means. A tilde on a single reported number is an
   unknown additive offset, not a per-metal structure.
2. It reaches all seven metals, because it needs no `d_eq`. Section 6.12's
   per-metal form reaches three.
3. It is **algebraically identical** to a common shift of all seven work
   functions, since `dW_i = W_i - w_cross`. That gives a free exact test
   (Validation 1 below) and it means the same sweep simultaneously bounds a
   systematic error in `METAL_WORK_FUNCTIONS` -- which is a real worry, given
   that Ni's entry already carries a 4.9-5.35 eV spread in its own comment.

How big should `delta` be? Three independent handles, none of them
authoritative, all quoted rather than averaged:

- Khomyakov *et al.*'s own `Delta_c(3.3 A) ~ 0.9 eV`, the quantity whose
  variation *is* the crossover's variation (Section 6.12, Eq. *): if `Delta_c`
  is known to ~20%, `delta ~ 0.2 eV`.
- The spread within `METAL_WORK_FUNCTIONS` itself: Ni is entered as 5.04 eV
  from a range commonly quoted as 4.9-5.35 eV, i.e. +-0.2 eV on one entry.
- The one experimental anchor this repo has flagged since 2026-09-18 and never
  reconciled: Mueller *et al.* (arXiv:0902.1479) measure a **0.12 eV**
  potential step where `METAL_WORK_FUNCTIONS` implies 0.25-1.07 eV.

**`|delta| <= 0.2 eV` is adopted as the nominal band and 0.5 eV as the outer
sweep.** Nothing below rests on the exact choice; every result is reported as
a function of `delta`.

## 3. The quantity being measured

For an ordered pair (A, B) define the **crossover sensitivity**

    S_pair  =  d ln|N| / d(delta)          [units: 1/eV]

evaluated by central difference at `delta = 0`, and the **2026-09-21-style
amplification**

    A_pair  =  ( d|N| / |N| ) / ( d<|dW|> / <|dW|> )

where `<|dW|>` is the mean of the two contacts' offset magnitudes. `A_pair` is
dimensionless and is directly comparable with 2026-09-21's numbers: that
session measured 4.55 for the same-sign pair Cu/Au and 0.07-0.11 for the two
straddling pairs Cu/Pt and Au/Pt, a 65-fold gap, **from a sample of three
pairs under a per-metal perturbation**. Whether that gap is a property of
near-cancelling pairs or a property of that particular perturbation is exactly
what a 21-pair enumeration under a *different* perturbation can decide.

There are 7 metals, so 21 unordered pairs, of which 21 are asymmetric
(A != B); ordered direction only flips the sign of `N`, so unordered is the
right count. "Straddling" means one contact above 5.4 eV and one below, i.e.
`dW_A * dW_B < 0`; "same-sign" means both on one side.

## 4. PRE-REGISTERED PREDICTIONS

These are committed before the model is written. Each is falsifiable and each
states a number.

**P1 -- the dichotomy generalises, cleanly.** Over all 21 asymmetric pairs,
`A_pair` separates without overlap: every straddling pair has `A_pair < 1.0`
and every same-sign pair has `A_pair > 2.0`, and
`min(same-sign) / max(straddling) > 10`. *Rationale:* 2026-09-21 argued the
amplification belongs to near-cancellation, which is a property of the pair's
sign structure, not of the perturbation. If that argument is right the
separation should be clean. If it is only a property of that perturbation, a
different perturbation should scramble the two groups.

**P2 -- the extremum is Au/Pd.** The single largest `A_pair` at
`delta = +0.1 eV` is Au/Pd, whose contacts are 0.02 eV apart -- the closest
pair in the table. *Rationale:* if amplification is near-cancellation, the
nearest-cancelling pair should be the worst. This is falsifiable in a way P1
is not: near-cancellation could be non-monotone in separation once the
kernel's spatial structure enters, in which case some other pair wins.

> **2026-09-25.** That is what happened, and this paragraph is the reason the
> failure could be scored rather than merely noticed. The verdict **HELD**
> recorded below was an artefact of a finite-difference step 3.5 decades too
> coarse; converged, Ti/Cu (0.32 eV apart) wins at 2.828, Cr/Au (0.60 eV) is
> second, Ni/Pd (0.08 eV) third, and Au/Pd fourth at 2.016. The two closest
> pairs in the table place third and fourth. Near-cancellation is **not**
> monotone in contact separation at finite `delta`.

**P3 -- Ti/Pt survives.** Chapter 6's headline pair changes `|N|` by **less
than 5%** over the whole nominal band `|delta| <= 0.2 eV`. *Rationale:* Ti/Pt
straddles, and the two previous refinements moved it 4.9% (2026-09-20) and
1.2-1.8% (2026-09-21).

**P4 -- the retraction table loses a member.** Section 6.11.2's claim that
**five** pairs gain more than 2x from a perfect mask changes membership
somewhere inside `|delta| <= 0.2 eV`: at least one of Au/Pd, Ti/Cr, Ni/Au,
Cr/Cu, Ni/Pd drops below 2x, or a sixth pair rises above it. *Rationale:* all
five are same-sign pairs, and P1 says same-sign pairs are exactly the
amplifying ones. If P1 is right, P4 should follow -- so P4 failing while P1
holds would be informative about the mask ceiling specifically.

**D5 -- a derivation, not a prediction, to be confirmed numerically.** `N`
changes sign for a pair only when `w_cross` moves *between* that pair's two
work functions, since that is the only way `sign(dW_A * dW_B)` can change.
Hence the reachable sign flips, in order of `|delta|`:

| pair | W_A, W_B | delta window for a flip |
|---|---|---|
| Au/Pd | 5.10, 5.12 | -0.30 to -0.28 |
| Ni/Pd | 5.04, 5.12 | -0.36 to -0.28 |
| Ni/Au | 5.04, 5.10 | -0.36 to -0.30 |
| Cr/Cu | 4.50, 4.65 | -0.90 to -0.75 |
| Ti/Cr | 4.33, 4.50 | -1.07 to -0.90 |

**The nearest sign flip is therefore at `delta ~ -0.29 eV`: outside the
nominal +-0.2 eV band, but comfortably inside what a tilde on 5.4 eV can
mean.** To be confirmed to 1e-4 eV by bisection. This is stated as a
derivation precisely so that it is not later mistaken for a successful
prediction -- the boundary follows from arithmetic on
`METAL_WORK_FUNCTIONS`, not from the model.

## 5. Validations planned, all against exactly known values

The 2026-09-17 and 2026-09-19 sessions each had a real bug survive a
plausible-range check and die on an exact one. Four exact checks:

1. **Shift invariance.** `N(W_A, W_B; w_cross = 5.4 + delta)` must equal
   `N(W_A + delta, W_B + delta; w_cross = 5.4)` -- the same algebra that makes
   a crossover offset and a common work-function offset the same
   perturbation. *This identity is exact in real arithmetic but is NOT
   expected to be bitwise in floating point*, because `4.33 - 5.5` and
   `(4.33 + 0.1) - 5.4` round differently. Predicted agreement: `<= 1e-15`
   absolute on `N`. Stating the expected failure mode in advance is the point;
   a bitwise claim here would be wrong.
2. **`delta = 0` reduction.** At `delta = 0` every pair must reproduce the
   2026-09-18 signed model **bitwise**, since it is the identical call.
3. **Symmetric cancellation at every `delta`.** Identical contacts at zero
   bias give `N == 0` for any offset whatever; must hold across the full
   sweep, not only at `delta = 0`.
4. **Charge conjugation at every `delta`.** `N(-dW) == -N(dW)`, implemented as
   `W -> 2*(5.4 + delta) - W`. 2026-09-18 got exactly `0.000e+00` on 27 cases
   at `delta = 0`; it must stay exact as the crossover moves.

Plus one identity carried over rather than re-derived: the Validation-5
kernel ceiling `|N[g]| <= max|k|` of Section 6.9.2, re-checked at every
`delta` in the mask sweep.

## 6. What this session will not claim

- It will not produce a *corrected* crossover. It measures sensitivity to an
  unknown offset; it does not determine the offset.
- A scalar offset is not the per-metal structure of Section 6.12. The two are
  complementary: 6.12 has the right shape for three metals, this has the wrong
  shape for all seven but the right reach. Neither supersedes the other, and
  the log should say so.
- Nothing here touches `lambda = 250 nm`, which is still one value for all
  metals, or the photo-thermoelectric term, which is still absent.

---

# OUTCOME (appended after the model was written and run)

Scoring is mechanical: `check_predictions()` in
`graphene_crossover_sensitivity_model.py` prints the verdicts below.

| # | prediction | verdict |
|---|---|---|
| P1 | clean dichotomy in `A_pair`, ratio > 10 | **FALSIFIED as written** |
| P1' | the same dichotomy in `S = dln\|N\|/d(delta)` | **HELD, 30.4x** — *re-measured 2026-09-25: still held, at 26.7x* |
| P2 | Au/Pd is the most sensitive pair | **HELD** — **RETRACTED 2026-09-25: falsified; the converged answer at `delta = +0.1 eV` is Ti/Cu. See `notes/2026-09-25-numeric-default-scale-audit.md` and thesis Section 6.13.11.** |
| P3 | Ti/Pt moves < 5% over the nominal band | **HELD, 0.45%** |
| P4 | the ">2x mask gain" membership changes | **HELD, twice** |
| D5 | 21 sign flips, nearest at -0.29 eV | **FALSIFIED: there are none** |

## 1. P1 falsified as written, and the reason is an exact identity

`A_pair` came out **infinite for all six straddling pairs**. Not large --
infinite, and exactly so. For a pair that straddles,

    <|dW>| = ((w_cross - W_A) + (W_B - w_cross)) / 2 = (W_B - W_A) / 2

which does not contain `w_cross` at all. A uniform crossover offset changes
the pair's response (`S` is −0.008 to −0.022 eV⁻¹, small but nonzero) while
changing the input measure by **exactly nothing**. `A_pair` divides by zero.
Checked to `0.000e+00` over 126 (pair, `delta`) combinations and promoted to
**Validation 5**, which was not planned; it was found while trying to score
P1 and turned out to be exact.

**This is the session's main methodological result.** The amplification ratio
2026-09-21 reported is not a property of a pair. It is a property of the
*(pair, perturbation)* couple, and it is undefined for half of this repo's
pairs under a perfectly reasonable perturbation. Quoting "65x amplification
for near-cancelling pairs" as if it were a device property, which the
2026-09-21 log came close to doing, would have been wrong.

## 2. P1's substance survives, in the right statistic

`S = dln|N|/d(delta)` is well defined for every pair:

| group | range of \|S\| (eV⁻¹) |
|---|---|
| straddling (6 pairs) | 0.0077 – 0.0218 |
| same-sign (15 pairs) | 0.6625 – 2.5628 |

`min(same-sign) / max(straddling) = 30.4`, with **no overlap**. So the
dichotomy 2026-09-21 saw on three pairs is real and generalises to all 21
under a *different* perturbation — which is the independent confirmation that
was wanted — but it survives only when stated in a statistic that does not
divide by the input. **A 65x gap measured on three pairs is a 30x gap on
twenty-one.** The direction of that revision is the one this repo has now
seen three times: widening the sample shrinks the claim.

## 3. D5 falsified, and the bug that nearly hid it

D5's premise — that `sign(N)` is set by `sign(dW_A * dW_B)` — is false. Write
`dW_A = m - s`, `dW_B = m + s` with `s = (W_B - W_A)/2`. A uniform offset
moves `m` and leaves `s` **exactly** fixed, so

    E(x) = m [f(L-x) - f(x)]  +  s [f(x) + f(L-x)]

splits into an antisymmetric part carrying `m`, which contributes zero to `N`
by the Validation-3 argument, and a symmetric part carrying `s`, which the
perturbation cannot touch. `sign(N)` is therefore fixed by which contact has
the higher work function, and **no scalar offset can reverse it**. A scan of
`delta` over ±3 eV, 601 points × 21 pairs = 12621 evaluations, finds **zero**
sign changes; the smallest `|N|` anywhere is 5.1e-03 (Au/Pd at +3 eV),
approached and never crossed. `m` only rescales `|N|`, vanishing as
`|m| → ∞` — and that limit *is* Section 6.11.2's near-cancellation.

**The first draft of `sign_flip_table()` reported 21 roots.** It bisected
between `delta = 0` and the midpoint of D5's window without ever checking
that a root was bracketed; with no sign change, the loop walks its lower
bound up to the upper bound and returns the **endpoint**. Every printed root
equalled `(W_A + W_B)/2 - 5.4`, every one looked physical, and the table's
headline — "nearest sign flip: Pd/Pt at −0.0150 eV" — is the kind of number
that goes straight into a thesis. It was caught by hand-checking one row
against what the field actually does at `dW_A = -dW_B`, **not by any of the
five exact validations**, all of which passed.

That is a new failure mode for this repo's log. 2026-09-17 and 2026-09-19
were caught *by* exact checks. 2026-09-21 established that exact checks
cannot reach a claim about the *size* of an effect. Today adds: exact checks
cannot reach a **root-finder that was never asked whether a root exists**,
because the validations test the model and the bug was in the analysis
wrapped around it. A bracket assertion is now the first line of that
function.

## 4. P3 held, and Ti/Pt is now the most robust thing in Chapter 6

0.45% over the whole nominal band, against 4.9% (2026-09-20) and 1.2–1.8%
(2026-09-21). Three independent refinements have now failed to move it. The
`w_cross` uncertainty is **not** where Ti/Pt's risk lives.

## 5. P4 held, twice, and it cuts both ways

| pair | δ=−0.20 | δ=−0.10 | δ=0 | δ=+0.10 | δ=+0.20 |
|---|---|---|---|---|---|
| Au/Pd | 2.985 | 5.741 | **8.421** | 11.055 | 13.663 |
| Ti/Cr | 2.735 | 3.066 | 3.396 | 3.727 | 4.056 |
| Ni/Au | 1.534 | **2.442** | 3.357 | 4.277 | 5.195 |
| Cr/Cu | 2.496 | 2.868 | 3.242 | 3.614 | 3.989 |
| Ni/Pd | 1.122 | 1.801 | **2.486** | 3.176 | 3.869 |
| Ti/Cu | 1.337 | 1.513 | 1.689 | 1.865 | **2.041** |

The count above 2x is **3, 4, 5, 5, 6** across the band. Section 6.11.2's
"five pairs gain more than 2x" is true only at `delta = 0`; at −0.2 eV it is
three and at +0.2 eV it is six. The *qualitative* claim that replaced
"masks are for symmetric devices only" — that masks help whenever
`|N_uniform|` is small — survives at every offset, and so does the design
recommendation, since Au/Pd's best masked response is 0.377 per incident
photon at `delta = 0` against unmasked Ti/Pt's 1.832. **What does not survive
is the number five.** It is an artefact of evaluating a threshold at one
point of a parameter the literature gives with a tilde.

The Section 6.9.2 ceiling `|N[g]| <= max|k|` held at every offset: 0
violations in 210 checks.

## 6. One prediction was quantitatively too tight, and it is recorded

The note predicted shift invariance would agree to `<= 1e-15` absolute. The
measured worst case is **1.776e-15** (Cu/Pd at `delta = +0.100`), with
215/231 combinations bitwise identical. The prediction is wrong by a factor
of 1.8. The right statement is *a few ulp of `N`*, not a fixed absolute
bound: 1.78e-15 is 8 ulp of 0.593. Small, but it is a pre-registered number
that missed, and the practice is worth nothing if only the comfortable
misses get recorded.

## 7. Not yet covered — updated

- **A second anchor for `Delta_c`** (from 2026-09-21) is now sharper: today
  shows the *scalar* part of the crossover uncertainty is nearly harmless for
  straddling pairs and dominant for same-sign ones. What is untested is the
  **differential** part, `w_cross(A) != w_cross(B)`, which is the only part
  that can move `s` — and `s` is what sets the sign. **A per-metal crossover
  difference of even 0.02 eV could flip Au/Pd, where a 3 eV scalar offset
  cannot.** That is the sharpest form the propagation question has taken, and
  it is the new top item.
- **Sections 6.9-6.11's numbers still stand at `delta = 0`** and now carry a
  stated band; the *five*-pair count is retracted as a count.
- A description of Ti, Ni and Pd that does not go through work function.
- Re-check other "for every ..." claims in the repo against widened samples;
  Chapter 4's per-metal `Rc` recalibration and Chapter 5's liner scenarios.
- Ask whether Chapter 4's `Rc` and Chapter 5's resistivity are
  near-cancellations — today gives the test a well-defined statistic to use
  (`S`, not `A`).
- Chapter 7 now has a **sixth** thread: three distinct classes of claim
  failure are on record (implementation error, caught by exact checks;
  unrepresentative sample, caught by enumeration; analysis-layer error,
  caught by neither).
- Chapters 2-3 remain undrafted.
- Mueller *et al.* 0.12 eV step vs the 0.25-1.07 eV `METAL_WORK_FUNCTIONS`
  offsets — and today makes it quantitative: a common error there is the same
  perturbation as `delta`, hence bounded by the `S` table.
- Photo-thermoelectric term; Shimomura comb electrodes; plasmonic/contact
  integration; Section 4.7 negative residual; Ti and Cr `Rc` recalibration;
  second edge-contact dataset.
