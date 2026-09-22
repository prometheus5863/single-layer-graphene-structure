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
