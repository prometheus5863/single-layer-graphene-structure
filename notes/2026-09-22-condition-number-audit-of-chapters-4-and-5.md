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

---

# OUTCOME (written after `graphene_sensitivity_audit.py` was run)

Full machine output: `audit_output.txt`. Figure: `sensitivity_audit.png`.

## 5. Verdict on the six predictions

| | prediction | status | measured |
|---|---|---|---|
| P1 | Euler sum rule `= 1` exactly | **PASS** | worst `|ΣS − 1| = 1.8e-15` over 15 decompositions |
| P2 | `ρ_GNR` `κ > 1` when the calibration is propagated | **PASS** | `κ = 1.551` at 52 nm |
| P2′ | …and stays below 3 | **PASS** | max 1.551 |
| P3 | rank order Ch5-calib > Ch4-`R_trans` > Ch5-liner > Ch5-`ρ` | **PASS** | 19.71 > 11.46 > 1.50 > 0.66 |
| P4 | ≥3 of 4 metals near-cancelling, Pd worst | **FAIL** | 2/4; Pd worst (correct) |
| P5 | Ch.6 amplifications reproduced within 10% | **PASS** | worst 1.9% off |
| P6 | a published number with a >100% bar at 5% input error | **FAIL** | largest bar 57% (Pd) |

**Four passed, two failed, and the two failures are the informative ones.**

## 6. Why P4 failed, and why the failure matters

P4 was `[hand-derivable]` and the hand arithmetic was *correct on the two
metals it was done for* — Pd 11.46 and Au 6.76, both above the
near-cancelling threshold of 3, both above Chapter 6's 4.55. What the hand
estimate did not do was **compute the other two**, and they come in at Ni
2.30 and Cu 1.23: *mildly* ill-conditioned, not near-cancelling.

This is the third consecutive session in which a claim generalised from a
partial enumeration failed on the full one (2026-09-20 on widened samples,
2026-09-21 on the pre-registered 10% bound, today on "at least three of
four"). The enumeration rule adopted on 2026-09-20 was applied to the
*model* this session but not to the *prediction* about it.

## 7. Why P6 failed, and what replaces it

No published Chapter 4 or Chapter 5 number carries a >100% bar under a 5%
single-input perturbation. The worst is Pd's residual at 57%. So the audit
did **not** find a printed figure whose order of magnitude is unestablished.

The reason is worth stating, because it is the opposite of what the
prediction assumed: `κ` is large exactly where the *output* is small, and a
57% bar on a number reported as `+51.0 Ω·µm` is not the same kind of problem
as a 57% bar on a headline. What P6 should have asked — and what Section 9
below asks instead — is whether the **sign** of a published number is
established, not its magnitude.

## 8. RESULT — the Chapter 5 calibration is the worst-conditioned step in the thesis

`κ = 19.71`, larger than anything in Chapter 4 (11.46) or Chapter 6 (4.55).

    residual = target_ratio − bulk_term − edge_term
             = 3.0000 − 1.0000 − 1.8478  =  0.1522
    λ_impurity = λ_bulk / residual = 361.4 nm

Three terms of order unity produce a result of order 0.15. The parameter-level
sensitivities are correspondingly severe: `d ln λ_imp / d ln ρ_calibration =
−19.71`, `d ln λ_imp / d ln λ_bulk = +13.14`, `d ln λ_imp / d ln W_calibration
= −12.14`. **A 5% error in the single calibration datum `ρ = 3.6 µΩ·cm at
22 nm` moves the extracted impurity mean free path by 99%.**

This is a real finding about Chapter 5 and it was invisible in the chapter's
own formula, which is a manifestly reinforcing Matthiessen sum (term-level
`κ ≤ 0.66` at every width — P2's first half, confirmed). **The
near-cancellation is one level down, in the calibration, and it propagates
into `ρ` damped but not removed: `κ(ρ_GNR)` rises from 0.88 at 18 nm to 1.55
at 52 nm.** A quantity can be well-conditioned in its own arguments and
ill-conditioned in the things those arguments were derived from; auditing
only the visible formula would have missed this entirely.

**An exact structural result found on the way** (Validation 4, not designed
in). At `W = W_calibration = 22 nm` the model reproduces `ρ_calibration`
**bitwise**, and `S(ρ_bulk) = S(p) = S(λ_bulk) = 0.0` **exactly** — three
independent exact zeros. At the calibration width the impurity term absorbs
whatever the others do, so the prediction there carries *no information* from
`ρ_bulk`, `p` or `λ_bulk`. `S(ρ_bulk)` then **changes sign through that
width**: `+0.120` at 18 nm, `0.000` at 22 nm, `−0.551` at 52 nm. A
one-calibration-point model inherits a sensitivity structure that pivots
about its calibration point, and Chapter 5 has never said so.

## 9. RESULT — Section 4.7's negative residual is *strengthened* by the audit

This is the session's most consequential finding and it closes off one branch
of an item open since 2026-08-31.

`κ` says how an input error is amplified. The operationally useful inverse is:
**how wrong would `R_extra` have to be to flip the *sign* of the residual?**
For `Q = R_c − R_extra` that fraction is `f = −Q / R_extra`:

| metal | residual (Ω·µm) | `κ` | `R_extra` must move by | sign is |
|---|---|---|---|---|
| Pd | **+50.9** | **11.46** | **−9.6 %** | **fragile** |
| Au | −90.2 | 6.76 | +14.8 % | intermediate |
| Ni | −360.5 | 2.30 | +43.4 % | intermediate |
| Cu | −787.5 | 1.23 | **+81.1 %** | **robust** |

**`κ` and sign-robustness run in exactly opposite order.** Pd — the *only*
positive residual, and the one Section 4.7 described as "a small, plausible
positive residual", i.e. the single data point consistent with the additive
decomposition surviving — is the one whose sign is *least* established: a
9.6% error in `R_extra` erases it. Cu, the largest and most apparently
unphysical violation, would need `R_extra` to be wrong by 81%.

Two consequences, both against the grain of the earlier reading:

1. **The additive decomposition `R_c = R_extra + R_transmission` fails
   robustly, not marginally.** Section 4.7 hedged that only Pd was
   "plausible"; the audit says Pd is the *weakest* evidence in the table,
   not the strongest.
2. **Amplified input error is eliminated as the cause of the negative
   residuals.** If the negatives were an artefact of a large `κ` acting on a
   mis-specified input, the *worst* negatives would sit at the *largest* `κ`.
   They sit at the smallest. The cause therefore lies where Section 4.7's
   first hypothesis put it — TLM double-counting the same near-contact
   region — or in `λ_decay`, and **not** in amplification. That is one
   candidate removed by measurement rather than by argument.

## 10. RESULT — the Cu liner model contains an unremarked subtraction

`W_eff = W − 2t` is a difference, and its `κ` rises from 1.130 at 52 nm to
**1.500 at 18 nm**, with `|S_t(ρ_eff)| = 1.167` there: at the narrow end of
the modelled range a 1% error in the 3 nm liner thickness is a 1.17% error in
the reported effective resistivity, and it grows without bound as `W → 2t`.
Chapter 5 states `t = 3.0 nm` from two literature sources that themselves
differ (3 nm vs a "2–3 nm functional floor"), i.e. a ~17% spread — which the
audit converts into a **~20% bar on `ρ_eff` at 18 nm**, previously unstated.

## 11. What the audit does not claim

1. It does **not** claim any number in Chapters 4 or 5 is wrong. `κ` is a
   property of a formula; it converts an input error into an output error and
   is silent on whether any input error exists.
2. The 5% input perturbation is **hypothetical** wherever the literature
   supplies no uncertainty, which is everywhere except the liner thickness.
   Every percentage above that derives from it is labelled accordingly.
3. `κ(ρ_GNR)` is audited over 18–52 nm only, the range Chapter 5 plots. It
   grows monotonically with `W` in the propagated case, so a wider wire would
   be worse conditioned, not better — untested beyond 52 nm.
4. Validation 3(c) (scale invariance of `κ`) is **3/5 bitwise, not 5/5**,
   with worst relative deviation `6.2e-16`. The invariance is exact in exact
   arithmetic; the two misses are the rounding of a division after rescaling
   by `1e±9`, a floating-point statement and not a modelling one. Recorded
   rather than rounded away.
5. No literature was fetched this session and none was needed: every input
   used was already cited in Chapters 4 and 5.

**Web search availability:** WebSearch/WebFetch **not used** — the session's
question was entirely internal to the existing model and citations. This is
recorded rather than left ambiguous.
