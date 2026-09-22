# Are Chapter 4's contact resistance and Chapter 5's resistivity near-cancellations?

**Date:** 2026-09-22
**Repo:** single-layer-graphene-structure
**Open item addressed:** "Ask whether Chapter 4's `Rc` and Chapter 5's
resistivity are near-cancellations" — created 2026-09-21 and described there
as *the sharpest form the synthesis question has taken*. Secondary: "Re-check
whether other 'for every ...' claims in the repo rest on small samples"
(2026-09-20) and "Isolating the root cause of the Section 4.7
negative-residual result" (open since 2026-08-31).

---

## 1. Why this question exists

On 2026-09-21 the per-metal crossover session measured a **65-fold
difference** in how two classes of quantity respond to the *same* 17.12%
input perturbation:

| class | example | output shift | amplification |
|---|---|---|---|
| same-sign (near-cancelling) pair | Cu/Au | 77.94 % | 4.55 |
| straddling (reinforcing) pair | Cu/Pt | 1.21 % | 0.07 |

The physics of the two cases is identical; what differs is that a same-sign
pair's response is a *difference of two comparable numbers*, so the
fractional error in the difference is the fractional error in the inputs
multiplied by the ratio of input magnitude to output magnitude.

That is not a fact about photodetectors. It is a fact about arithmetic, and
it applies to **every** quantity this thesis computes. Chapters 4 and 5 have
never been asked which kind they are. This session asks.

## 2. The diagnostic, stated before it is run

For a computed quantity `Q(x_1, ..., x_n)`, the **signed logarithmic
sensitivity** of `Q` to input `x_i` is

    S_i = (x_i / Q) * (dQ / dx_i)                                      (C1)

read as: *a 1% change in `x_i` produces an `S_i` % change in `Q`.* The
**condition number** of `Q` is

    kappa(Q) = max_i |S_i|                                             (C2)

The classification this session will use throughout:

* `kappa <= 1` — **reinforcing**. Input errors are damped or passed through.
  A sum of same-signed positive terms is always in this class.
* `1 < kappa < 3` — **mildly ill-conditioned**.
* `kappa >= 3` — **near-cancelling**. An input error is amplified. Chapter
  6's same-sign pairs sit here at 4.55.

`kappa` is a property of the *formula*, not of the physics, and it is
computable without knowing whether any input is right.

## 3. PRE-REGISTERED PREDICTIONS

These are committed to git **before `graphene_sensitivity_audit.py` exists**,
following the procedure adopted on 2026-09-21. The commit order is the
evidence.

**Disclosure on prediction integrity.** Two of the quantities below can be
condition-numbered by hand from numbers this thesis has already published,
and it would be dishonest to present those as blind. Each prediction is
therefore labelled `[blind]` or `[hand-derivable]`, and for the
hand-derivable ones the hand estimate is written down here so that the
machine result can be checked against it rather than compared to nothing.

**P1 — exact identity, must hold.** Every decomposition audited here is a
homogeneous function of degree 1 in the inputs it is decomposed into (a sum
or a difference of them). Euler's homogeneous-function theorem then forces

    sum_i S_i = 1        exactly                                       (C3)

This must hold to machine precision for **every** quantity audited, in both
classes. If it fails anywhere, the implementation is wrong and no result
below may be believed. This is the session's exact validation, and it is
deliberately one that cannot be satisfied by an approximately-correct
derivative: a finite-difference `dQ/dx` that is slightly off will not sum to
1.000000000.

**P2 — `[blind]`.** Chapter 5's graphene-nanoribbon resistivity
`resistivity_vs_width()` is a Matthiessen sum of three strictly positive
terms, so audited **against its own arguments** it must be reinforcing:
`kappa <= 1` for every input at every width in 18-52 nm. But
`lambda_impurity` is not an argument — it is *solved* from a calibration
point by `_solve_lambda_impurity_nm()`, and that solve contains a
subtraction. Prediction: audited **against the calibration inputs**
(`rho_calibration`, `rho_bulk`, `p_default`, `W_calibration`, `lambda_bulk`),
`kappa > 1` somewhere in 18-52 nm — i.e. *the chapter with no manifest
subtraction in its headline formula still conceals a near-cancellation, one
level down.*

*Hand estimate, disclosed:* at the calibration point the residual is
`3.0 - 1.0 - 1.848 = 0.152`, so the solve itself has `kappa ~ 20`; the
impurity term carries only ~8% of the total ratio at the wide end, which
should damp it to order 1.5 in `rho` itself. The prediction `kappa > 1` is
therefore *expected* to pass; what is blind is **by how much, and whether the
damping is enough to keep it under 3** (the near-cancelling threshold).
Stated sharply: **`kappa(rho_GNR)` stays below 3 across the whole width
range.**

**P3 — `[blind]`.** Rank order of the four audited families by worst
`kappa`, largest first:

    Ch5 calibration solve  >  Ch4 R_transmission  >  Ch5 liner  >  Ch5 rho_GNR

The first two places are hand-derivable; **the last two are fully blind** and
are what P3 is actually betting on.

**P4 — `[hand-derivable]`.** Chapter 4's `R_transmission = Rc_measured -
R_extra` is near-cancelling, `kappa >= 3`, for at least three of the four
recalibrated metals, and the metal with the **smallest** residual (Pd,
+51.0 of 584) has the **largest** `kappa`.

*Hand estimate, disclosed:* Pd gives `584 / 50.9 = 11.5`; Au `609.2 / 90.2 =
6.8`. Both exceed Chapter 6's 4.55. This is arithmetic on Section 4.7's
published table, not a prediction, and is recorded so that a machine result
disagreeing with it exposes an implementation bug.

**P5 — `[blind]`.** Applying the *same general machinery* to Chapter 6's
2026-09-21 pairs reproduces the published amplifications (4.55 for Cu/Au,
0.07-0.11 for the straddling pairs) to within **10%**. This is the
cross-chapter check that the diagnostic is measuring the thing 2026-09-21
measured, rather than a differently-normalised cousin of it.

**P6 — `[blind]`.** At least one number **currently published** in Chapter 4
or Chapter 5 carries a `>100%` error bar under a 5% perturbation of a single
input — i.e. the audit finds a printed figure whose sign or order of
magnitude is not established by the model that produced it.

## 4. What this session will *not* do

* It will not change any input value, recalibrate anything, or "fix" a large
  `kappa`. A large condition number is a property of a formula; it is
  reported, annotated in place, and carried as an error bar.
* It will not delete or overwrite Section 4.7's table or Chapter 5's
  resistivity numbers. Per the repo's standing rule, superseded numbers are
  annotated, not removed.
* It will not treat `kappa` as an error bar on its own. `kappa` multiplies
  an input uncertainty; where no input uncertainty is available, the audit
  reports `kappa` and the *hypothetical* 5% bar, labelled as hypothetical.
