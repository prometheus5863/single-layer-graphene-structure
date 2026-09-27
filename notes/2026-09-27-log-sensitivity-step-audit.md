# The step behind Chapter 5's conditioning numbers, and what its convergence estimate cannot see

**Date:** 2026-09-27. Pre-registration (committed first, commit `e2eed30`):
`notes/2026-09-27-log-sensitivity-preregistration.md`. Code:
`graphene_log_sensitivity_step_audit.py`. Output:
`log_sensitivity_step_audit_output.txt`. Figure:
`log_sensitivity_step_audit.png`. Live web search **not used** — this is an
internal audit of this repository's own numerics and needed none.

## The item, and why it was the top one

`graphene_sensitivity_audit.log_sensitivity(f, p, rel_step=1e-5)` is the only
finite-difference estimator in this thesis's conditioning work. Every
parameter-level sensitivity in **Chapter 5 Section 5.5** comes out of it,
including the `κ = 19.71` λ_impurity calibration the chapter calls the
worst-conditioned step in the thesis. It reports with every value
`conv = |S(h) − S(2h)|`, described in its own docstring as "a convergence
estimate that is reported with every value rather than assumed small."

That is step-doubling — the pattern 2026-09-25 caught answering *yes* while the
answer was 44% wrong and the estimate was blind by five orders of magnitude.
It has been the top open numerical item since.

**A correction to the pre-registration's own framing, recorded because it
matters for scope.** Q1 and Q4 named "Chapter 4/5 call sites". There are no
Chapter 4 call sites. Section 4.10's numbers are all term-level and closed
form (`S_i = T_i/Q`, no finite differences, no tolerance); this estimator's
entire exposure is Chapter 5. The pre-registration asserted a call site it had
not checked for, which is a smaller version of 09-26 item 6's finding that a
census read signatures instead of call sites.

## Finding 1 — the estimator is not the estimator its docstring names

The two evaluation points are `p(1 ± h)`. In `L = ln p` they sit at
`L + ln(1+h)` and `L + ln(1−h)`, and `ln(1+h) ≠ −ln(1−h)`. So it is a **secant
over an asymmetric log interval**, and a secant equals the derivative at its
interval's **midpoint**:

    L + ½ln(1−h²) = L − h²/2 + O(h⁴)

Expanding with `a = ln(1+h)`, `b = ln(1−h)`, `s = a+b`, `d = a−b`:

    S_affine(h) = g′ − (h²/2)g″ + (h²/6)g‴ + O(h⁴)        c₂ = −g″/2 + g‴/6
    S_geom(h)   = g′ +            (h²/6)g‴ + O(h⁴)        c₂ =        g‴/6

The `−(h²/2)g″` term is a **displacement of the evaluation point**, not a
truncation error, and it is proportional to `dS/d ln p` — a quantity a central
difference's error does not involve at all. **This is 2026-09-26 item 5's fault
class** ("an asymmetric secant is a second-order estimate at that interval's
midpoint, so the caller who asked for the derivative at `L` receives the
derivative somewhere else, correctly computed"), found there on a guard path
that had never fired and found here on the **default** path.

**Proved, not asserted.** For `ln|f| = A(ln p)²` the old estimator's error is

    E_affine(h) = A·ln(1−h²)     EXACTLY, at every order

while a central difference on the same function is **exactly zero**, because a
central difference of a quadratic has no error. One function separates the two
estimators by a closed form at every step. Measured agreement with that closed
form: at the derived cancellation bound throughout (worst measured/bound 2.08).

**Magnitude: 9.47e-10 relative, worst over 23 call sites.** No published number
moves. The finding is about the estimator's shape and about what any future
step audit should expect its error to scale with.

Fixed in `log_sensitivity(..., symmetric=True)`, now the default; the old path
is kept as `symmetric=False` because 09-24 established that deleting a
superseded implementation destroys the only oracle for judging its replacement.
Running `graphene_sensitivity_audit.py` before and after: **every value printed
to four decimals is identical**, the sole exception being one `0.0000`
becoming `-0.0000`.

## Finding 2 — the fix costs one of Validation 4's three exact zeros, and that is informative

Validation 4(c) has claimed since 09-22 that at the calibration width
`S(ρ_bulk) = S(p) = S(λ_bulk) = 0` exactly — "three independent exact zeros,
the strongest form of check this repo uses." Under the corrected stepping
`S(ρ_bulk) = −1.11e-11`, so the **bitwise** count is 2/3.

`1.11e-11` is the derived cancellation floor `4ε|ln ρ|/(2h) = 5.7e-11` at
`h = 1e-5`. All three zeros survive as mathematics; one has lost its *bitwise*
exactness, and that exactness turns out to have depended on the old
estimator's particular evaluation points happening to cancel. **An exact-zero
check written as `v == 0.0` cannot distinguish structure from floating-point
luck.** Both counts are now printed, the criterion is the derived floor, and
the annotation is in the source. No structural claim changes: at the
calibration width the prediction still carries no information from `ρ_bulk`,
`p` or `λ_bulk`.

## Finding 3 — THE RESULT: `conv` reports machine precision where `S` is 14% wrong

`conv` is exactly 3× conservative where truncation dominates — measured
3.0000, 3.0001, 3.0006, 3.0054 at `h = 1e-3 … 3e-2`, which is Q2's prediction
confirmed to four digits. It vanishes wherever the `h²` and `h⁴` error terms
cancel between `h` and `2h`, i.e. wherever `c₂` and `c₄` have opposite signs.

Q4 predicted **no** such step at any real call site, reasoning that these are
smooth rationals in which `g″` dominates `g‴`. **Q4 is falsified.** Four zeros
sit on the truncation branch — located by bracketed bisection (09-24's rule)
against the derived rounding floor rather than by eyeballing a sweep — and all
four are in **Family B, the λ_impurity calibration**:

| call site | h at conv = 0 | conv there | error in S | relative |
|---|---|---|---|---|
| `λ_imp / ρ_calibration` | 2.928e-2 | 4.97e-14 | 2.76 | **14.0 %** |
| `λ_imp / ρ_bulk` | 3.084e-2 | 1.42e-14 | 2.62 | **13.3 %** |
| `λ_imp / W_calibration` | 4.391e-2 | 1.24e-14 | 1.83 | **15.1 %** |
| `λ_imp / λ_bulk` | 4.752e-2 | 1.95e-14 | 1.69 | **12.9 %** |

Understatement factor up to **1.8e14**. At `h ≈ 3%` the only instrument the
code offers reports convergence to machine precision while the sensitivity is
wrong by 14%. That is 09-25's finding reproduced on this estimator, at 14%
instead of 44%, on the step this thesis calls its worst-conditioned.

The reason Q4's argument failed is the reason Family B was the right place to
look: `λ_imp = λ_bulk/residual` with a near-zero residual puts large higher
derivatives into `ln|f|`, so `c₂` and `c₄` acquire opposite signs. **The
near-cancellation that makes `κ = 19.71` is the same near-cancellation that
creates the blind spot** — one mechanism, two symptoms, and Section 5.5 had
measured only the first.

**What is NOT claimed.** Nothing published moves. Every `S` in Sections 4.10
and 5.5 is stable to **1.3e-8** relative under 16× refinement, against the
corrected estimator, and against a Richardson reference (Q5 PASS). The nearest
zero is **≈3000× above** the shipped `h = 1e-5`, not near it; the shipped step
is safe by a wide margin and this session does not claim otherwise. What is
unsafe is the **procedure**: a sweep over the natural range `h ∈ [1e-8, 1e-1]`
passes straight through all four points, and the repo ran exactly such a sweep
on 09-25 and 09-26.

**The exactly-solvable half, which is what licenses the above.** For
`ln|f| = A L² + B L⁴` the secant is exact in closed form, so a blind spot can
be constructed rather than found: with `A = −2B + δ` at `L = 1`, `conv` vanishes
near `h*² = 3δ/(20B)`, and — independently of `A`, `B` and `δ` — the relative
error there is `(3/8)h*²`. Measured ratio to that prediction **1.0668** (Q3
FAIL on its magnitude band; the `h⁶` term contributes 6.7% at `h* = 0.0179`;
class call correct). So the *ratio* by which `conv` understates is unbounded
while the *absolute* error it hides is of ordinary `O(h²)` size — which is the
honest limit on the mechanism in the constructed case, and precisely what the
real Family B result exceeds, because there `h*` is 3% rather than 1.8%.

## Finding 4 — `conv` is minimised where the answer is worst (unpredicted)

If a step-size sweep minimised `conv`, which `h` would it pick?

- median `h` at min error / `h` at min `conv` = **877**
- minimising `conv` gives a **worse** answer than the shipped `h = 1e-5` at
  **18 of 23** call sites, worst penalty **3.5e4×**

Far below the cancellation knee `conv` differences two noise samples and is
small for that reason. **A sweep steered by `conv` walks away from the answer
while reporting success.** This is the operationally important finding: the
blind spots of Finding 3 are four isolated points, but this is the whole
low-`h` half of the range.

## Finding 5 — for two of the three families, `conv` measures rounding

`conv` at the shipped step, divided by the derived rounding floor
`4ε|ln Q|/h`:

| family | min | median | max |
|---|---|---|---|
| B (λ_impurity calibration) | 10.2 | **431** | 1461 |
| C (ρ_GNR(W), λ re-solved) | 0.0 | **0.4** | 7.3 |
| D (Cu liner, `W_eff = W − 2t`) | 0.5 | **0.5** | 0.6 |

Family B's `conv` is a real truncation estimate. **Families C and D sit at the
floor**, so Section 5.5's "worst convergence estimate over the table:
2.7e-10" and "worst convergence estimate: 1.9e-10" lines are measurements of
double precision, and would print roughly the same number for a model with any
amount of curvature. The reassurance scales with `|ln Q|` and `h`, not with the
quality of the answer.

## Finding 6 — no caller ever compares `conv` with anything

Eleven `log_sensitivity` call sites; **zero** compare the returned `conv`
against a threshold. Every one prints it. Q7's class call is right and its
scope was wrong: one *unrelated* module (`graphene_crossover_sensitivity_model`)
does compare its own convergence quantity to a tolerance, so the repository
knows how and simply did not here. A convergence estimate with no criterion
attached is a number, not a check — and Findings 3–5 are what that number was
hiding.

## Validations, and the fault class this session committed three times

| | what it checks | result |
|---|---|---|
| **V1/V4** | power law: `S = n`, `conv = 0`, every `h` | **retained BECAUSE it cannot discriminate** |
| **V2** | `E_affine = A·ln(1−h²)` exactly | worst measured/derived-bound **2.08** |
| **V3** | `S[1/f] = −S[f]` | 27/31 bitwise, worst/bound **0.22** |
| **V5** | the constructed blind spot, closed form | measured/bound **0.49** |

**Three of the five failed in their first form, and all three failed the same
way: a round absolute tolerance — 1e-13, 1e-14, 1e-13 — that silently encoded
a step size.** V2 required 1e-13 and measured 5.10e-08, which *is* the derived
floor at `h = 1e-3`. V5 required 1e-13 and measured its own floor. V3 required
1e-14 and was wrong in *shape*, not only magnitude: `ln|1/f| = −ln|f|` holds in
exact arithmetic, but `1/f` is a separately rounded number, so `log(1/x)` is
not `−log(x)` bitwise. All three are rebuilt on derived bounds; all three first
forms are left in the record.

**This is the fault class documented on 2026-09-25 (item 10) and reproduced on
2026-09-26 (item 9), committed three more times by the session auditing it.**
Two sessions ago the conclusion was "knowing the class demonstrably does not
prevent reproducing it, and what caught it both times was a test that reports a
number instead of asserting a verdict." Today is the third consecutive
occurrence and the third consecutive rescue by the same mechanism. The
inference is no longer about carelessness: **an absolute tolerance is the
default way a person writes a numerical assertion, and the only reliable
defence found so far is structural — report the measured value and the derived
bound, and let the ratio be the verdict.**

Q4 deserves the same treatment. It was falsified twice over in one session:
first as written (243 sign changes over `h ∈ [1e-8, 1e-1]`, all cancellation
noise, because the prediction named a range it had not thought about), then on
its class call (4 real zeros on the truncation branch). Two boundary
definitions for "truncation branch" were tried and discarded before one that is
derived — `|conv| > 30·4ε|ln Q|/h` — worked. `10 × argmin|conv|` is two to
three decades too low, for the reason Finding 4 gives; a log-log slope test
reads noise as slope +2 whenever two adjacent noise samples happen to rise.

## Scorecard

**PASS:** Q1 (class and magnitude), Q2, Q5, Q7 (class).
**FAIL:** Q3 (magnitude band, class right), Q4 (as written *and* on class),
Q6 (Family B worse by 4.35×, not >10×), Q7 (as written — over-scoped).
**Unpredicted:** Q8 (Finding 4), Q9 (Finding 3's real zeros), the
Chapter-4-has-no-call-sites correction, and Finding 2.

## Not covered here

- The 14% error at `h ≈ 3%` was measured against a Richardson reference, not
  against a closed form, because Family B has none. Its own accuracy is
  `O(h⁴)` in truncation and `~1e-12` in cancellation, which is ample at the
  1e-1 level being measured but is stated rather than assumed.
- Whether the same blind-spot analysis applies to `graphene_default_scale_audit`'s
  own `rel_step` sweeps, which used this estimator on 09-25.
- Whether any *other* near-cancellation in this repo (Chapter 4 §4.7's
  residual, Chapter 6's same-sign pairs) would show the same `c₂`/`c₄` sign
  flip if it were differentiated rather than evaluated.
