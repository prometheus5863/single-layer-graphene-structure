# Pre-registration — `log_sensitivity`'s step and its convergence estimate

**Date:** 2026-09-27. **Committed before the audit module was written**, per the
practice established 2026-09-21 and continued 09-22/09-23/09-25/09-26.

## Why this item

`graphene_sensitivity_audit.log_sensitivity(f, p, rel_step=1e-5)` is the single
estimator behind **every parameter-level sensitivity published in Chapter 4
Section 4.10 and Chapter 5 Section 5.5**, including the `κ = 19.71` calibration
result that Section 5.5 calls the worst-conditioned step in this thesis.

It reports, with every value, `conv = |S(h) − S(2h)|` — described in its own
docstring as "a convergence estimate that is reported with every value rather
than assumed small."

**That is step-doubling.** 2026-09-25 showed a step-doubling convergence
estimate answering *yes* while the answer was wrong by 44%, and blind by five
orders of magnitude. 2026-09-26 then showed that refining one step converges to
the derivative of whatever other discretisation was held fixed. This estimator
has never been checked against an anchored criterion, and it has been the top
open numerical item since 2026-09-25.

## Pre-registered questions, with numbers attached

**Q1 — the estimator is not the estimator its docstring names.** The two
evaluation points are `p(1 ± h)`. In the variable the derivative is taken with
respect to — `L = ln p` — those sit at `L + ln(1+h)` and `L + ln(1−h)`, and
`ln(1+h) ≠ −ln(1−h)`. **Prediction:** the estimator is a secant over an
*asymmetric* log interval whose midpoint is displaced from `L` by
`½ln(1−h²) ≈ −h²/2`, so its leading truncation error contains a term
`−(h²/2)·g''` — proportional to `dS/d ln p` — that a true centred difference
does **not** have. Predicted leading coefficient `c₂ = −g''/2 + g'''/6`
against a centred difference's `g'''/6`. **Predicted magnitude at the shipped
`h = 1e-5`: below `1e-9` relative on every Chapter 4/5 call site** — i.e. this
is a claim about the estimator's *shape*, not about any published number. Class
call: this is 2026-09-26 item 5's fault class (an asymmetric secant returns the
derivative *somewhere else*, correctly computed) on a **default** path rather
than on a guard path that never fired.

**Q2 — `conv` is conservative in the asymptotic regime, by exactly 3×.** If
`E(h) = c₂h² + O(h⁴)` then `conv = |E(h) − E(2h)| = 3|c₂|h² + O(h⁴)` while the
true error is `|c₂|h²`. **Prediction: measured `conv/|error| = 3.00 ± 0.02`**
for a function with a known exact `S`, at an `h` large enough that truncation
dominates cancellation (`h = 1e-2`).

**Q3 — `conv` has exact blind spots, and a step-size sweep is the procedure
that finds them.** `conv` vanishes wherever `E(h) = E(2h)`, which for
`E = c₂h² + c₄h⁴` happens at `h*² = −c₂/(4c₄)`, i.e. whenever `c₂` and `c₄`
have opposite signs. **Prediction: such an `h*` is constructible in closed
form, `conv(h*)` is zero to machine precision there, and the relative error in
`S` at that step is exactly `(3/8)h*²` — independent of the function's
coefficients.** So the *ratio* by which `conv` understates the error is
unbounded (it is 0/nonzero) while the *absolute* error hidden is of ordinary
`O(h²)` size. Both halves of that are predicted; the second is the honest limit
on how bad this can be.

**Q4 — whether the real call sites have such a blind spot.** **Prediction: NO**
— for Family A (`R_transmission = Rc − R_extra`, Chapter 4 §4.7/§4.10) and
Family B (the `λ_impurity` calibration, Chapter 5 §5.5) no `h*` exists in
`h ∈ [1e-8, 1e-1]`, because those functions are smooth rationals in which
`g''` dominates `g'''`, giving `c₂` and `c₄` the same sign and `conv/|error| ≥ 3`
throughout. Stated as a falsifiable *no*, with the range named.

**Q5 — whether the published numbers move.** **Prediction: every `S` published
in §4.10 and §5.5 moves by less than `1e-6` relative** under both (i) 16×
step refinement and (ii) the corrected symmetric estimator. If Q1 is right and
Q5 is right, the finding is structural and no chapter number changes.

**Q6 — Family B is the near-cancellation, so it should be the worst.**
`κ = 19.71` there against `11.46` in Chapter 4. **Prediction: Family B's step
error exceeds Family A's by more than 10×**, because the `1/residual` in
`λ_imp = λ_bulk/residual` puts a near-zero in a denominator and near-zeros in
denominators are where log-space curvature lives.

**Q7 — `conv` is never compared with anything.** **Prediction: no call site in
the repository asserts on `conv`; every one of them prints it.** If so, the
docstring's "reported with every value rather than assumed small" is true of
*reported* and false of everything that follows from it — a convergence
estimate with no criterion attached is a number, not a check.

## Validations to be built, and what each can and cannot discriminate

- **[V1] Pure power law `f = A·p^n` → `S = n` exactly, for every `h`.** Must be
  bitwise. **This check cannot discriminate** the shipped estimator from the
  corrected one, nor `h = 1e-5` from `h = 0.5`: for a power law the log-space
  integrand is linear and every secant over it is exact. It is retained *because*
  it cannot discriminate, and labelled as such — continuing 09-25 item 11 and
  09-26 item 10. **Validation 2 of `graphene_sensitivity_audit.py` is entirely
  built out of this case**, which is the point.
- **[V2] The discriminating exact check.** For `ln|f| = A(ln p)²` the shipped
  estimator's error is `A·ln(1−h²)` **exactly, to all orders**, while a true
  centred difference is **exactly zero** on the same function (a centred
  difference of a quadratic has no error). One function separates the two
  estimators by a closed form at every `h`. Predicted agreement with the closed
  form: `≤ 1e-15` relative.
- **[V3] A symmetry that must give exactly zero.** `S[1/f] = −S[f]` bitwise,
  since `ln|1/f| = −ln|f|` and the estimator is linear in `g`. Predicted
  `31/31` bitwise over a step sweep.
- **[V4] `conv` for a power law must be exactly `0.0` bitwise at every `h`,
  including an absurd `h = 0.5`** — the non-discriminating property of [V1] and
  of `conv` itself, made explicit as a measurement.
- **[V5] The blind spot, exactly.** For a quartic `g(L)`, the secant is exact
  in closed form, so `E(h)` is known with no truncation at all; `h*` is found by
  a **bracketed** root solve (09-24's rule) and `E(h*)` compared against the
  closed form. Predicted `≤ 1e-15` relative.

**Scoring rule, as always:** predictions are scored on class call *and*
magnitude band separately, and failures stay in the record.
