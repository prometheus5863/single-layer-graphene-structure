# Pre-registration — the covariance probe's own blind spots (2026-09-28)

Written and committed **before** the audit module exists, continuing the
practice established 2026-09-21 and kept every session since.  Five questions
Q1–Q5, five validations V1–V5.  Scored verbatim in
`AUTOMATION_LOG.md` whether they pass or fail.

## Why this target

2026-09-27 left as its **top open numerical item**: *whether
`graphene_default_scale_audit`'s own `rel_step` sweeps are affected by the
finding that minimising `conv` walks a median 877× away from the right step.*
Reading the module to scope that question surfaced something sharper, so this
session answers the logged item **and** the thing found while scoping it.

`graphene_default_scale_audit.py` is the module that judges every numeric
default in this repository.  Two of its own components are now suspect for
reasons the repository itself established in the last three sessions:

1. **Pre-registered prediction D2** reads, verbatim:

   > `log_sensitivity(rel_step=1e-5)` passes covariance EXACTLY (bitwise zero
   > deviation), **because `rel_step` multiplies the parameter.**

   It is scored by `d2 = _logsens_invariance(1e-5, 1e3)[0] == 0.0`.  Two
   things have happened to that line since it was written:

   - 2026-09-27 changed `log_sensitivity`'s default path so that `rel_step`
     **no longer multiplies the parameter**.  It exponentiates it: the samples
     are `p·exp(±h)`, not `p·(1 ± h)`.  D2's stated *mechanism* has been
     removed from the code.  If D2 still passes, its explanation is wrong.
   - 2026-09-27 item 10 established that **an exact-zero check written as
     `v == 0.0` cannot distinguish structure from floating-point luck** — that
     is exactly how D2 is written — and item 11 established that a **power
     law is non-discriminating** for anything about a log-derivative
     estimator, because it makes `ln|f|` linear in `ln p` so every secant is
     exact.  D2's test function is `f(p) = 3p^2.5`.  A pure power law.

2. **The probe's scale factor `lam = 1e3` is a one-sided choice.**
   `covariance_deviation(solve, lam)` measures the deviation at a single λ, and
   every conviction in RESULT 2 rests on λ = 1e3, i.e. on rescaling the
   variable **up** by three decades.  Nothing in the module argues that the
   direction is irrelevant.  Chapter 7 §7.5's theme is unequal scrutiny; this
   would be an instance inside the auditing machinery itself.

## Questions

**Q1 — does D2 still pass, and if so is its stated reason still true?**
Prediction: **D2 still returns bitwise `0.0`** under the corrected
`symmetric=True` estimator, and therefore its stated reason ("because
`rel_step` multiplies the parameter") is **falsified by its own continued
success** — a test that passes identically after its named mechanism is
removed was never testing that mechanism.  Predicted class: PASS-for-the-
wrong-reason, not a numerical failure.

**Q2 — is D2's bitwise zero a property of the estimator, or of λ = 1e3?**
Prediction: it is a **floating-point round-trip coincidence**.  For
`f2(p) = 3(p/λ)^2.5` evaluated at `p = λ`, the sample points reach the power
law through `(λ·exp(±h))/λ`, and multiply-then-divide by the same constant is
*usually* but **not always** exact in IEEE-754.  Prediction: sweeping λ over
the decades `1e1 … 1e15` at the shipped `h = 1e-5`, **at least one decade
gives a non-zero deviation**, and the failing fraction is **≥ 20%**.  Stated
so it can fail: if all fifteen decades give bitwise zero, Q2 is wrong and
something structural is protecting the test.

**Q3 — is D2's test function non-discriminating in the 09-27 item 11 sense?**
Prediction: **yes, and demonstrably.** Replacing `3p^2.5` with a function of
the same `S` magnitude but genuine curvature in `ln p` (`ln|f| = A(ln p)^2`
with `A` chosen so `S(1) = 2.5`) makes the covariance deviation **non-zero,
≥ 1e-13**, for *both* the old and the corrected estimator, while the exact `S`
is unchanged by the rescaling in both cases.  If so, D2's word "EXACTLY" is a
statement about the test function and not about `rel_step`.

**Q4 — the logged top item: does any conclusion in
`graphene_default_scale_audit` rest on a `conv`-guided step choice?**
Prediction: **NO — a null result, stated so it can fail.** Census of every
step choice in the module: `log_sensitivity` is called at exactly one site
(`_logsens_invariance`) and its returned `conv` is **discarded into `_`**;
`_hdiff_solver` and RESULT 7's sweep use `sensitivity()` with an *absolute*
step and never compute `conv` at all.  Predicted: zero sites select a step by
minimising `conv`, so none of 2026-09-25's conclusions inherit 09-27 item 7,
and the logged item closes as a genuine null.  **Prediction of a null result,
and it fails if any site is found.**

**Q5 — is the probe's verdict one-sided in λ?**
Prediction: **yes, and it reverses.** `H_DIFF = 1e-3` is the default D3
convicted as scale-encoding, at λ = 1e3.  An absolute step on a variable in eV
has a relative error that scales like the step over the variable's scale, so
rescaling the variable **down** shrinks the deviation.  Prediction: the same
probe run at **λ = 1e-3** reports a deviation at least **1e4× smaller** than
at λ = 1e3, small enough that the module's own conviction threshold would
**acquit** the default it convicted.  The probe measures `|ln λ|` in one
direction only.

## Validations (each against an exactly known value, per 09-27 item 12)

Every validation below reports its **measured value beside a derived bound**
and lets the ratio be the verdict.  No round absolute tolerances: three of
five validations failed that way on 09-27, inside the session auditing that
very fault class, and the conclusion drawn there was that the defence has to
be structural rather than exhortative.

- **V1 — identity.** λ = 1 is the identity map, so the covariance deviation
  must be **bitwise 0.0** for every probe target, good default or bad.  A
  non-zero here would mean the harness varies, not the default.  Exact.
- **V2 — the power law's exact `S`.** For `f = 3p^2.5`, `S = 2.5` exactly, at
  every `h`, for both estimators.  Measured `|S − 2.5|` is compared against
  the derived rounding floor `4ε·|ln f|/(2h)`, not against a round number.
- **V3 — the corrected estimator is exact on a log-quadratic.** A central
  difference of a quadratic has **zero** truncation error, so for
  `ln|f| = A(ln p)^2` the `symmetric=True` estimator must return `2A ln p`
  with only rounding error, while `symmetric=False` must be wrong by
  `A·ln(1−h²)` **exactly**.  Two exact statements, one function, and they
  disagree — so this validation discriminates, which V2's power law cannot.
- **V4 — the round-trip mechanism is demonstrated, not asserted.** Q2's
  mechanism requires that `(λ·a)/λ != a` for some `(λ, a)` in double
  precision.  Validation: search the λ decades and the shipped `h` and
  **count** the non-round-trip pairs; the count must be **> 0** for Q2's
  mechanism to be the explanation, and the count is reported whatever it is.
  If the count is 0 while deviations are non-zero, the mechanism is wrong and
  that is recorded as the finding.
- **V5 — odd symmetry.** The deviation functional evaluated on a solver whose
  output is exactly proportional to its scale must vanish identically for
  every λ — an exact zero that does not depend on rounding, giving V1 a
  companion that varies λ rather than fixing it.

## What this session will NOT claim

No published number is expected to move.  Q1–Q3 concern a **prediction's
stated reason** and a **validation's discriminating power**, not any physical
result; Q4 is predicted to close as a null; Q5 concerns the **direction** of a
probe, not the correctness of the default it convicted (`H_DIFF`'s plateau
margin was measured independently in RESULT 3 and is not in question here).
Superseded text will be annotated in place, per the standing rule since
2026-09-24, and not deleted.
