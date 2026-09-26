# The two entangled defaults behind g_ds, and the first clean exoneration

**Date:** 2026-09-26
**Code:** `graphene_gds_quadrature_audit.py`,
`gds_quadrature_audit_output.txt`, `gds_quadrature_audit.png`
**Touches:** `rf_small_signal_model.py` (a fix), `thesis_draft/04-…` §4.6/§4.6.1,
`graphene_default_scale_audit.py` and `notes/2026-09-25-…` (in-place annotations)
**Score:** 4/6 pre-registered predictions, 5/5 exact validations (two of them
after a first form failed)

---

## 1. What was carried in, and why it was carried rather than done

2026-09-25's default-scale census found
`rf_small_signal_model.output_conductance(dVds=1e-3)` and logged it as the top
numerical item without measuring it. The reasoning recorded that day was
explicit: `dVds` is an absolute step **2% of its variable** — the identical
construction to `H_DIFF = 1e-3`, which the same session had convicted of
moving eleven of Chapter 6.13's sensitivities by up to **44%** and of
falsifying a prediction carried as HELD since 09-22 — but `g_ds` feeds
`f_max`, a Chapter 4 published quantity, so moving it needed its own
before/after comparison rather than a footnote.

That deferral was correct and this note is the comparison. The headline is a
**null result**, and the most useful thing in the session is the structure
that makes the null result meaningful instead of merely reassuring.

## 2. Why this case is structurally worse than `H_DIFF`

`H_DIFF` differentiates a closed-form expression. `dVds` differentiates
`graphene_fet_model.transfer_characteristic()`, which is **itself a
discretised object**:

```python
V_channel_profile = np.linspace(0, Vds, n_segments)   # n_segments = 50
R_channel_local   = np.mean([channel_resistance(V_g, V_ch)
                             for V_ch in V_channel_profile])
Id[i] = Vds / (R_channel_local + Rc_total)
```

The quadrature grid **is set by `Vds`**. So differencing in `Vds` differences
the quadrature error along with the physics, and two defaults are entangled
and not independent: a step size `dVds = 1e-3` and a resolution
`n_segments = 50`.

This makes a sharper failure mode available than anything in the 09-22…09-25
series. 09-25's anchored step criterion — accept a step only if the derivative
is unchanged at `h/10` *and* `h/100` — tests convergence **in the step**. A
step-refinement study of a discretised function converges to the derivative of
**the discretisation you fixed**, not to the derivative of the model. It can
therefore pass at any tolerance while sitting on a bias it cannot see. The
anchored criterion was the fix for the *previous* failure class; here is a
class it does not reach.

## 3. Result: both defaults are innocent, and the item closes

Five `V_g` points, `Vds = 0.05 V` (the signature default) and `0.1 V` (the
call site — see §6).

| what | measured | prediction |
|---|---|---|
| step `dVds = 1e-3` vs. anchored-plateau step | **1.25 × 10⁻⁶** | Q1 **PASS** |
| `n_segments = 50` vs. `n → ∞` (Richardson in 1/(n−1)) | **8.4 × 10⁻⁷** | Q2 **FAIL** (band) |
| order of the quadrature error in 1/n | ratios 2.031 → 2.001 ⇒ **first order** | Q3 **PASS** |
| propagated to peak `f_max` | **1.1 × 10⁻⁷** (0.00001%) | Q5 **PASS** |
| `V_g` grid (400 pts) behind `g_m`, hence peak `f_T` | **1.0 × 10⁻⁵** under 16× refinement | Q6 **PASS** |

`f_T` contains no `g_ds` at all and does not move. **Q2 failed on its
magnitude band (predicted ≥ 0.1%) and not on its class call** — the mechanism
it described is exactly what was observed, four decades smaller than
predicted. The item is **closed**, not carried: this is the first clean
exoneration in the audit series that began 09-22.

One number in that table is worth pausing on. The `g_ds` term carries
**0.9999** of the `f_max` denominator at peak `f_max`, so the dilution should
have been the square root's 0.5. It measured **0.229**. The remainder comes
from where the two curves peak: `f_max` is evaluated at its own argmax, and a
near-flat maximum converts a shift in the curve into a smaller shift in the
peak value. A weight of 0.9999 in the denominator does not imply a weight of
0.9999 in the reported figure of merit, and that gap is about the reporting
statistic rather than the physics.

## 4. The mechanism is real, and it is not small in general

The exoneration is a fact about graphene's channel resistance over a 50 mV
drop, not a general licence. Two measurements separate those.

**(a) The bias is provably invisible to the step criterion.** On the real
model at `V_g = 2.0 V`, the `n=50`-vs-`n=6400` gap is

| `dVds` | gap |
|---|---|
| 1e-2 | +3.089 × 10⁻⁷ |
| 1e-3 | +3.039 × 10⁻⁷ |
| 1e-5 | +3.038 × 10⁻⁷ |
| 1e-7 | +3.038 × 10⁻⁷ |

drifting **2.3 × 10⁻⁴ of itself across four decades of step refinement**.
Refining the step removes none of it, exactly as §2 predicts.

**(b) In an exactly-solvable case the same mechanism is 1.32%.** Take
`R(V_ch) = a + b·V_ch²`. Then both averages are closed forms:

- discrete mean of `V_ch²` over an endpoint-inclusive `linspace(0,V,n)`
  = `V²(2n−1) / (6(n−1))`
- continuum mean over `[0,V]` = `V²/3`
- difference = **`V²/(6(n−1))`**, exact

so `Id(V) = V/(A + b c V²)` and `g(V) = (A − b c V²)/(A + b c V²)²` with
`c = c_n` for the discretisation and `c = 1/3` for the model. Measured
`c_n = 0.3367346938775511` against the exact `0.3367346938775510`. The
central difference of the discretised model converges to the **discretised**
closed form at second order in `h` (order test 3.9995), and its residual
against the **continuum** closed form converges to `1.3220563663 × 10⁻²`
where the exact quadrature gap is `1.3220563698 × 10⁻²` — i.e. to the gap,
not to zero.

So the criterion "the derivative is unchanged at `h/10` and `h/100`" is
satisfied to any tolerance one likes while the answer is 1.32% wrong, and the
only reason it is 8 × 10⁻⁷ wrong here is that `R(V_ch)` is nearly linear
across this device's drain drop. Validation C makes that precise and
oracle-free: the endpoint-inclusive mean is **exact** for a constant and
**exact** for a linear integrand (mean of an arithmetic sequence, worst
residual 1.6 × 10⁻¹⁶ over n ∈ {2,3,50,501}), and fails only from curvature
onward. A higher-`V_ds` model, or a shorter channel with real pinch-off
curvature, does not inherit this chapter's exoneration.

## 5. A real latent bug, exactly

`output_conductance` read

```python
Id_minus = transfer_characteristic(Vg, Vds=max(Vds - dVds, 1e-4))
gds = (Id_plus - Id_minus) / (2 * dVds)
```

The guard moves the **interval** without changing the **divisor**. When it
fires, the returned value is the true secant slope times exactly
`W/(2·dVds)` with `W = Vds + dVds − 1e-4` — an algebraic factor with no
approximation in it, which is why validation E reproduces it to `rel err 0.0`:

```
measured shipped/true = 0.9158333333333334
exact     W/(2*dVds)  = 0.9158333333333334      (Vds=0.05, dVds=0.06)
```

Against the anchored reference at `V_g = 2.0 V`: **7.98160187e-04 S vs
8.71919960e-04 S, −8.46%**, silently. Q4 **FAILED on its magnitude band**
(predicted > 10%) and not its class call.

It is reachable for `Vds ≤ dVds + 1e-4` — 1.1 mV with the default step — and
**from any upward step-size sweep**, which is precisely what a step-size
audit does. The bug was one careless sweep away from contaminating an audit
of itself. It never fired at either operating point in use, so no published
number was ever affected.

The fix **shrinks the step symmetrically** rather than widening the interval,
and that choice is the substantive part. Widening is the smaller edit, but an
asymmetric secant over `[1e-4, Vds+dVds]` is a second-order estimate of
`dId/dVds` at that interval's **midpoint**, not at `Vds` — the caller who
asked for `g_ds` at `Vds` would have received `g_ds` somewhere else,
correctly computed. Shrinking keeps the difference centred where it was
asked for: same point, **+0.0017%** instead of −8.46%. A shrink now warns
(09-24: a correct guard that silently degrades a correct number is worse than
one that fails loudly) and `Vds ≤ 1e-4` raises instead of returning a number
with no defensible meaning. The non-firing path is unchanged expression by
expression and verified **bitwise** — not to a tolerance — over
`Vds ∈ {0.05, 0.1, 0.2} × dVds ∈ {1e-4, 1e-3, 1e-2}` on the standard
400-point sweep, with warnings promoted to errors.

## 6. Unpredicted: the census read a signature, not a call site

This was not pre-registered and is recorded as such. This module's peak `f_T`
came out **10.140 GHz** against Chapter 4's published "**≈20 GHz**". Neither
is wrong: `output_conductance`'s *signature* default is `Vds = 0.05 V`, while
`plot_fT_fmax()` — the path that produces `rf_figures_of_merit.png` and the
chapter's numbers — passes `Vds = 0.1 V`, where the model gives
**20.279 / 18.731 GHz**. One model, two bias points, and §4.6 never named
which.

The consequence reaches back into yesterday's instrument. The census scored
this default as "a step **2%** of its variable" by reading the signature. At
the call site it is **1%**. The census was low **by exactly the factor between
the signature default and the call site** — here 2×, and in general unbounded,
since nothing constrains how far a call site sits from a signature.

**A default-scale census must read call sites, not signatures.** Both the
census code and the 09-25 note are annotated in place; the 2e-2 stays, right
about the signature and wrong about the claim it was aimed at. §4.6 now names
its operating point.

Re-measured at the call site: step error 7.3 × 10⁻⁷, quadrature bias
1.4 × 10⁻⁶ — **1.66×** the `Vds = 0.05` value, not the 4× a pure `V²`
curvature scaling would give, because the local curvature of `R(V_ch)` is
itself bias-dependent. The guard remains latent at both points.

## 7. The auditor committed the audited error, for the second time in two days

Two of five validations failed in their first form, both for the same reason,
and both failures are recorded in the source rather than quietly corrected.

**[B]** required the residual below an absolute `1e-12` and measured
`6.08e-10` at `h = 1e-9`. The arithmetic was right and the assertion wrong:
`6.08e-10` **is** the cancellation floor `eps·|Id|/(2h|g|)` at that step, so
the constant `1e-12` silently encoded `h ≳ 1e-7`. That is 09-25 item 10's
fault class — an absolute tolerance is a hidden scale claim — reappearing
inside the auditor **one day later, in a fresh module, written by someone who
had just documented it.**

**[A]** asserted that the continuum residual would exceed the discretised one
by "100×" and measured 67×. Same family: a round number standing in for a
derived one.

Both were rebuilt on derived quantities — an oracle-free order test and the
exact algebraic gap — and both now pass. **[B] is retained and labelled
non-discriminating** (09-25 item 11): it is exact for every `(n, h)` tried
and passes for a good default and an absurd one alike, so it proves the
harness and nothing about the default. **[D]** states its cancellation bound
as a derivation and reports measured/bound ratios (0.62, 1.06) rather than
clearing a constant.

The lesson is not "be more careful." An absolute tolerance is the natural
thing to write, it passes in the case you wrote it for, and knowing about the
failure class demonstrably does not prevent writing it again. What caught it
both times was **a test that reports a number instead of asserting a
verdict** — the same property that made 09-25's convergence accessor useful.

## 8. Methodological note, continuing the series

09-20: exact validation does not protect against an unrepresentative sample.
09-21: pre-registration reaches what exact validation cannot.
09-22: pre-registration does not reach the analysis layer.
09-23: the other side of the ledger — how claims consolidate.
09-24: a correct guard silently degraded a correct number because its default
encoded the scale of its original caller.
09-25: a procedure asked whether it has converged can answer yes and be wrong
by 44%; only an anchored comparison distinguishes them.

**09-26: the anchored comparison is anchored in one variable.** Refining a
step converges to the derivative of whatever *other* discretisation was held
fixed, and no amount of step refinement reveals that. When two numerical
defaults are entangled — and they are entangled whenever one sets the other's
grid — convergence in one is not evidence about the pair.

And a second, quieter point, which is the one this session earned the hard
way: **a null result is worth a session only if it is a null result about a
mechanism you have shown can be large.** "We refined the step and nothing
moved" is compatible with the criterion being blind. "We refined the step,
nothing moved, the bias is provably invariant under that refinement, and the
same mechanism measures 1.32% in a case with real curvature" says something
about this device instead of about the absence of evidence. The first
formulation would have closed the item too, and would have closed it for the
wrong reason.

## 9. Left open

- **Every remaining `h`-like step and tolerance in `rf_small_signal_model.py`
  and the photodetector modules** — open since 09-25, and today's §6 sharpens
  the method: score them at **call sites**. The census's signature-reading is
  a systematic error in that whole list, not a one-off.
- **`log_sensitivity`'s convergence estimate is step-doubling** — open since
  09-25 and now the top numerical item by default. Chapter 4's
  `R_transmission` decomposition and Chapter 5's liner scenarios rest on the
  pattern 09-25 showed blind by five orders of magnitude, never checked
  against an anchored criterion.
- **Is any Chapter 4 or 5 *ranking* a step artefact?** Open since 09-25 and
  untouched. Today adds nothing for or against it.
- **`n_segments = 50` elsewhere.** Today exonerated it for `g_ds` at two bias
  points. It also sets every `Id` in Chapter 4, and the exoneration was
  curvature-dependent — so the natural follow-up is whether any *other*
  Chapter 4 quantity differentiates through that quadrature at a bias where
  the curvature is larger. Created today.
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open
  since 09-21 and still the top *physics* item, untouched today for the sixth
  consecutive session.
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural
  weakness.
- **Chapter 7 (discussion/outlook)** — now ten threads, material entirely in
  hand, displaced again.
- **Chapters 2–3 remain undrafted** despite complete computational results.
- **`__pycache__` is tracked in this repo** — open since 09-24; it dirties the
  working tree on every run and cost real time on 09-25.
