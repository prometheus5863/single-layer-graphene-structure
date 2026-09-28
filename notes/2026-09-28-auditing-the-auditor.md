# Auditing the auditor: what `graphene_default_scale_audit`'s own probe measures

**2026-09-28.** Pre-registered in
`notes/2026-09-28-covariance-probe-preregistration.md`, committed before the
audit module existed.  Code and full output:
`graphene_covariance_probe_audit.py`, `default_scale_probe_audit_output.txt`,
`default_scale_probe_audit.png`.  Scored: **Q1, Q4 PASS; Q2, Q3, Q5 FAIL**;
5/5 validations against derived bounds, with V5's first (bitwise) form failing
and kept in the source.

No published number moves.  Every finding below is about a **prediction's
stated reason**, a **validation's discriminating power**, or a **probe's
choice of direction**.

---

## 1. The logged top item closes as a null, and the null is load-bearing

2026-09-27 item 7 found that minimising `conv` chooses a step a median **877×**
below the one that minimises the error.  The obvious follow-up was whether this
repository's own `rel_step` sweeps had been guided by `conv`.  They had not.

The census is by AST rather than by eye: every `log_sensitivity` call in
`graphene_default_scale_audit.py` (two, both inside `_logsens_invariance`)
binds the second return value to `_` — discarded at the call site. The `H_DIFF`
and RESULT 7 sweeps go through `sensitivity()`, which uses an **absolute** step
and never computes `conv` at all. Six regex patterns for
`conv`-minimisation match nothing.

A null from a detector is worth only as much as the detector's demonstrated
sensitivity, so the census is given a **positive control**: the same census run
on a snippet that *does* pick a step by minimising `conv`, which it detects
(one call site, one pattern hit). Without that control the null would be
indistinguishable from a census that cannot see anything. This is the same move
as mutation-testing a self-checking suite, applied to a *static* census.

**Verdict: none of 2026-09-25's conclusions inherit 09-27 item 7. Item closed.**

## 2. D2 passes identically after the mechanism it names was deleted

`graphene_default_scale_audit`'s pre-registered prediction D2 reads, verbatim:

> `log_sensitivity(rel_step=1e-5)` passes covariance EXACTLY (bitwise zero
> deviation), **because `rel_step` multiplies the parameter.**

On 2026-09-27 the default path was changed so that `rel_step`
**exponentiates** the parameter: the samples are `p·exp(±h)`, not `p·(1 ± h)`.
The mechanism D2 names is no longer on the default path. D2 still returns
**bitwise `0.0`**, and it returns bitwise `0.0` on the old path too — the same
verdict with and without its stated cause.

A test whose verdict is unchanged by the removal of its explanation was never
testing that explanation. **D2's number stands; D2's reason is falsified.**

## 3. What the check actually tests — and it is not about graphene

This was the session's **unpredicted** result, and it subsumes two of its own
failed predictions.

`_logsens_invariance(rel_step, lam)` compares `S` for `f(p)` at `p = 1` against
`S` for `f(p/λ)` at `p = λ`. The scaled branch reaches the function through
`(λ·a)/λ` where `a = exp(±h)`. **When that round trip is exact, the scaled and
unscaled branches are literally the same floating-point expression**, so the
deviation is bitwise zero — for any function, any estimator, any step. Nothing
about the derivative is being tested.

Checked as an **exact set equality**, with no tolerance anywhere:

| estimator | test function | λ with non-zero deviation | λ with inexact round trip | equal |
|---|---|---|---|---|
| symmetric | power law `3p^2.5` | {1e7} | {1e7} | yes |
| symmetric | log-quadratic | {1e7} | {1e7} | yes |
| old | power law `3p^2.5` | {1e4} | {1e4} | yes |
| old | log-quadratic | {1e4} | {1e4} | yes |

over λ = 1e1 … 1e15. The two estimators fail at **different** decades — `1e7`
and `1e4` — which is the cleanest available evidence that the exactness is
arithmetic and not structural: nothing about the mathematics distinguishes
those decades, and nothing about the estimators distinguishes them either
except which sample points they happen to evaluate.

Three of this session's own scored items fall out of this one mechanism:

- **Q2 FAILED on magnitude** (non-zero at 7% of decades, predicted ≥ 20%) —
  the round trip is inexact only rarely. Class call correct.
- **Q3 FAILED outright.** I predicted a function with curvature in `ln p` would
  break the exactness, on the strength of 09-27 item 11's power-law finding.
  It does not, because the test never reaches the estimator. **The power-law
  criticism of D2 is true and is not the operative one** — a sharper criticism
  was available and the pre-registration reached for the familiar one.
- **V5 FAILED as first written**, one level up: `covariance_deviation` divides
  by `λ` as well. An exactly-proportional solver — covariant by construction —
  gives a worst deviation of **1.29e-16**, not zero, against a derived one-ulp
  bound of `eps = 2.22e-16`. So the probe's docstring claim that the deviation
  is *"identically zero for a scale-free procedure"* is **false**; it is zero
  to one ulp, and Validation 1's bitwise zero survives only because λ = 1
  makes the round trip trivial.

## 4. The probe's dynamic range is entirely in the direction it never runs

Every conviction in the audited module's RESULT 2 is measured at **λ = 1e3**:
the variable's unit made three decades smaller. Q5 predicted the probe is
one-sided, and predicted that running **downward** would be the *lenient*
direction. **That was backwards.**

| λ | `H_DIFF = 1e-3` (absolute) | relative step (scale-free) |
|---|---|---|
| 1e-6 | 2.12e+12 | ≤ 1 ulp |
| 1e-3 | 2.47e+01 | ≤ 1 ulp |
| 1e-2 | 7.02e-02 | ≤ 1 ulp |
| 1e+1 | 1.85e-03 | ≤ 1 ulp |
| **1e+3** | **1.327e-03  ← the conviction** | ≤ 1 ulp |
| 1e+6 | 1.327e-03 | ≤ 1 ulp |

Downward is **1.86e4× more severe**, not more lenient. And the structure is not
one curve with two halves:

- **Upward the deviation saturates.** λ = 1e3 and λ = 1e6 agree to 5e-5
  relative. The probe has a **ceiling**, so its magnitude cannot express
  severity: it would report the same 1.3e-3 for this default at λ = 1e15.
- **Downward it diverges**, and steepens faster than any power of λ — the
  difference quotient breaking down rather than a law. Stated as a breakdown,
  not fitted: two or three points are not a law (the discipline 09-27 applied
  to the divisor-invariance question in the other repository).

The convictions are **correct** — `H_DIFF` does carry dimensions, and the
scale-free alternative sits at one ulp in both directions, so the asymmetry
belongs to the default and not to the probe's arithmetic. What is wrong is
reading the reported magnitude as *how badly* a default encodes a scale.

## 5. A true claim, withdrawn by the wrong oracle, then recovered

The saturated ceiling looked like it should *be* the shipped difference's own
relative truncation error in the original units. That is a checkable claim, so
it was checked — and **the first check rejected it.**

A Richardson reference built from `(h, h/2)` implied a relative error of
**1.53e-2**, missing the saturated 1.327e-3 by **11.6×**, and the module duly
printed *"the identification DOES NOT HOLD, and the claim is withdrawn."*
Direct refinement to `h = 1e-6` gives **1.3288e-3** — the same number to
**0.14%**.

Richardson is the invalid one here, and the evidence is in the refinement
sequence itself:

| h | `1/S` | relative vs h = 1e-6 |
|---|---|---|
| 1e-2 | −0.3888439521 | 2.131e-03 |
| 3e-3 | −0.3889975904 | 1.737e-03 |
| 1e-3 | −0.3901922819 | 1.329e-03 |
| 3e-4 | −0.3878410937 | 4.705e-03 ← **the error grew as h shrank** |
| 1e-4 | −0.3894691862 | 5.268e-04 |
| 3e-5 | −0.3896560278 | 4.731e-05 |

Richardson extrapolation assumes a smooth `h²` error expansion. This quotient
does not have one, so extrapolating from two points of a non-monotone sequence
produces a confident wrong reference.

This is the **mirror image** of the last three sessions' theme. 09-25: a
procedure asked whether it had converged answered *yes* while 44% wrong. 09-27:
`conv` reported 1e-14 where `S` was 14% wrong. Today: an oracle reported
**failure where the claim was right**. Both first forms are left in the source,
per 09-24, because the failed oracle is the only record of how the right answer
was nearly discarded.

## 6. What this says, stated as narrowly as the evidence allows

Five weeks of these sessions have accumulated one methodological result per day
about *instruments*. Today's is about **explanations attached to instruments**.

A passing check carries two things: a **verdict** and, usually in a comment or
a docstring, a **reason it passes**. The verdict is tested every time the suite
runs. **The reason is tested never.** D2's reason survived the deletion of the
code it referred to, in a repository that re-runs its audits daily, because
nothing in a green suite is sensitive to why it is green.

The concrete defence found today is the one used in §3: **name the mechanism as
a separate, exactly-checkable proposition** — here, "the deviation is non-zero
exactly where `(λa)/λ ≠ a`" — and check *that*, as a set equality, with no
tolerance. That turns an explanation from prose into a test. It is the same
structural move as 09-27's "report the measured value beside a derived bound
and let the ratio be the verdict", one level up: **report the mechanism beside
a prediction it must reproduce, and let the set equality be the verdict.**

Second, smaller, and pointed at this session: **Q3 reached for the most recent
finding rather than the sharpest available one.** 09-27's power-law lesson was
one day old and it supplied a criticism of D2 that is true but not operative.
The pre-registration is what made that visible — without it the note would have
recorded the power-law criticism as the finding and never looked for the round
trip.

---

## Open, and sharpened by today

- **Whether other `== 0.0` exactness checks in this repository are round-trip
  tautologies.** Created today. The mechanism generalises to any check that
  compares a quantity against a rescaled recomputation of itself, and the
  repository has several.
- **Whether RESULT 2's other convictions are also read from the saturated
  branch** — created today. The conviction of `tol=1e-14` on a metre-scale root
  is the one 09-24 measured independently, so it has an anchor; the others do
  not.
- **A probe that reports both directions of λ**, since the whole dynamic range
  is downward. Cheap, and today's table is most of it.
- **Migrating every remaining validation in the repo to measured-value-beside-
  derived-bound form** — created 09-27 item 12, and today's V5 is one more
  instance of the same fault class caught by the same defence.
