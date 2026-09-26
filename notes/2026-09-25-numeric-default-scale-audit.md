# The numeric-default scale audit: a detector without an oracle, and the published claim it overturned

**2026-09-25.** Closes the item the 2026-09-24 session created in the act of
closing the previous one: *"audit every other numeric DEFAULT in the repo for
the scale assumption it encodes — the bracket class is closed; this class has
one measured member and **no detector**."*

Code: `graphene_default_scale_audit.py`. Recorded output:
`default_scale_audit_output.txt`. Figure: `default_scale_audit.png`. Follow-on
code: `sensitivity_converged()` in
`graphene_crossover_sensitivity_model.py`, with its recorded output in
`sensitivity_converged_output.txt`.

---

## 1. Why "no detector" was the operative half of that sentence

On 2026-09-24 a correct guard was moved to a new call site and silently cost
ten significant digits, because its default tolerance — an *absolute* interval
width of `1e-14`, chosen where the variable was of order 1 eV — was applied to
a decay length of order `5e-11` m. Nothing raised. The number was simply wrong
in its fifth figure.

It was found only because the unguarded loop it replaced had been returning the
right answer for three days and could be used as an **oracle**. That is not a
method. Most defaults in this repo have no older correct implementation
standing beside them, so a technique that needs one detects nothing. The open
item was not "check the other tolerances by hand"; it was "find something that
can check them."

## 2. The detector: a unit-covariance probe

A numerical default parameterises a procedure `P` applied to a variable `x`.
Change the variable's units, `x → λx`, and restate **the same mathematical
problem** in the new units. Nothing about the mathematics has moved, so the
answer must transform covariantly:

```
P_λ(f_λ)  ==  λ · P(f)          exactly, up to rounding
```

A default that is a pure ratio — a relative tolerance, a relative step, a
halving count — satisfies this identically. A default carrying the dimensions
of `x` does not, and the deviation is what the probe returns. No oracle, no
reference implementation, no correct answer: the probe tests a procedure
against its own units.

Applied to the five tolerance-class defaults in the repo, all five verdicts
matched pre-registration:

| default | worst deviation | verdict |
|---|---|---|
| `bracketed_bisect(tol=1e-14)` | 4.9e-02 | **SCALE-BOUND** (the known member) |
| `bracketed_bisect(tol=0, rtol=eps)` | 0.0 | scale-free |
| `bracketed_bisect(tol=0, rtol=0)` — 200 halvings | 0.0 | scale-free |
| `H_DIFF = 1e-3` (central-difference step, eV) | 2.5e+01 | **SCALE-BOUND** |
| `h = 0.5% of DELTA_NOMINAL` (the relative alternative) | 1.4e-16 | scale-free |

The probe reproduces 2026-09-24's fault with nothing borrowed from it.

**One property of the probe, learned by a failed prediction.** D1 predicted the
deviation for `tol=1e-14` would land in `1e-6 … 1e-2`; it is 4.9e-02 at
`λ = 1e-3` and 4.8e-05 at `λ = 1e+3`. The band was written from the `λ > 1`
direction alone, and **the probe is asymmetric in λ**: shrinking the variable
makes an absolute tolerance coarse relative to the interval (one or two
halvings — catastrophic), while growing it makes the tolerance fine, and the
residual deviation is then dominated by the error of the `λ = 1` *baseline*
rather than of the rescaled run. So a deviation from this probe is a reliable
yes/no and **not** a calibrated error estimate. The same asymmetry
`bracketed_root.py` already records for the unguarded loop, met again in a
different guise.

## 3. Three classes, not one — the 09-24 item conflated two of them

That item listed "grid counts" alongside "`tol` everywhere". They need
different instruments, and saying so is most of what makes the audit finite:

| class | n | instrument |
|---|---|---|
| TOLERANCE-like — carries the dimensions of a variable | 5 | the covariance probe |
| COUNT-like — a dimensionless discretisation count | 4 | a convergence study; **the probe says nothing** |
| PHYSICAL — a statement about graphene, not about numerics | 10 | out of scope |

A count is dimensionless, so it passes the probe trivially; what a count claims
is that the integrand has no structure finer than the grid, and only a
convergence study reaches that. Ten of the nineteen defaults are modelling
choices — `T = 300 K`, `hopping = 2.8 eV`, `sigma = 20 nm`, `L_junction = 1 µm`
— and probing those would report covariance failures that mean nothing. The
PHYSICAL rows are enumerated anyway, so the audit's boundary is auditable
rather than asserted. The sharpest instance: `p2_threshold(tol=0.10)` *looks*
like a numerical tolerance and is not — 10% is the definition of the predicate
P2, and tightening it would change the question, not the accuracy.

The two counts that had never been studied were studied. `n_points = 400` in
`junction_extra_resistance` is converged to 6.6e-06 relative against a 25 600-node
reference (D6 held). `max_iter = 200` needs 61 halvings on the widest interval
the repo bisects, a margin of 3.3× (D5 held).

## 4. The screen is not the verdict

Failing covariance means a default **encodes** a scale. It does not mean the
number is wrong. So the audit is two-stage and both stages are reported:

* **(a) the screen** — is this default scale-free? Cheap, needs nothing.
* **(b) the margin** — for the ones that encode a scale, how far is the shipped
  value from where it stops working?

Reporting (a) as the verdict would condemn correct numbers, which is the mirror
image of the 09-24 fault and the more expensive mistake, because it moves
published tables. Both defaults the screen flagged went to stage (b). One was
acquitted — `tol=1e-14` is never *used* at a bad scale any more, since the one
metre-scale caller passes `rtol`. The other was convicted.

## 5. H_DIFF is 3.5 decades above its own plateau

`H_DIFF = 1e-3 eV` is the central-difference step behind every logarithmic
sensitivity in thesis Section 6.13. The plateau of that derivative — the range
of `h` over which `S(h)` is stationary in `h` — was measured for the first time
today and runs from about `1e-10` to `3e-7 eV`. The shipped step sits **3.5
decades above its top.** Prediction D4, written as a null result precisely so
it could fail, failed.

The damage is uneven and the unevenness is the story. All six straddling pairs
agree to better than 1 part in 10⁵. Four same-sign pairs move under 0.4%. The
other eleven move 5.9% to **44.4%** — and not all in the same direction: Cr/Ni
and Cu/Au were *under*-estimated, the other nine over-estimated, so no uniform
correction factor would have helped.

**The dichotomy of Section 6.13.4 — that section's only load-bearing result —
survives**, at 26.7× rather than 30.4×, with no overlap. That is the outcome to
lead with, because it is the one a reader of the chapter needs.

## 6. Prediction P2 is falsified, and it is falsified in the way it said it might be

P2, pre-registered in the 2026-09-22 note and recorded there as **HELD**:
the largest `|S|` at `δ = +0.1 eV` belongs to **Au/Pd**, on the rationale that
amplification is near-cancellation and Au/Pd is the nearest-cancelling pair
(0.02 eV apart). Converged:

| rank | pair | `|S|` (eV⁻¹) | separation |
|---|---|---|---|
| 1 | **Ti/Cu** | 2.828 | 0.32 eV |
| 2 | Cr/Au | 2.518 | 0.60 eV |
| 3 | Ni/Pd | 2.291 | 0.08 eV |
| 4 | Au/Pd | 2.016 | 0.02 eV |

identical at every step inside the plateau, over four decades. The two closest
pairs place third and fourth; a pair 16× further apart wins. **Near-cancellation
is not monotone in contact separation at finite `δ`** — which is exactly the
escape hatch the 2026-09-22 note wrote into P2's rationale: *"near-cancellation
could be non-monotone in separation once the kernel's spatial structure enters,
in which case some other pair wins."* Some other pair wins. Writing the
falsifying mechanism into the prediction is what converts this from a surprise
into a scored result, and it is the strongest evidence so far for the
pre-registration habit adopted on 2026-09-21.

The mechanism of the artefact generalises past this table. Au/Pd's `|S|` drifts
**2%** across five decades of step size; Ti/Cu's collapses **69%**. Au/Pd did
not lead at `h = 1e-3` by being the most sensitive pair — it led by being the
pair whose derivative that step happened to get right, while its competitors'
were destroyed. **A ranking taken at a step outside the plateau ranks its
entries by how well the step suits them.** No amount of care about the physics
finds that; only a sweep over `h` does.

## 7. The fifth failure class: a procedure can certify its own convergence and be wrong by 44%

The obvious fix was already in the repo. `log_sensitivity` has returned
`|S(h) − S(2h)|` alongside every value since 2026-09-22 — "a convergence
estimate that is reported with every value rather than assumed small." The first
version of `sensitivity_converged()` adopted that pattern, and it is **blind**:

| pair | true error at `h = 1e-3` | step-doubling estimate | ratio |
|---|---|---|---|
| Ti/Ni | 44.38% | 2.33e-06 | **190 268×** |
| Cu/Ni | 10.16% | 5.48e-05 | 1 854× |
| Cr/Pd | 31.99% | 3.96e-04 | 808× |

Ti/Ni certifies itself converged to seven digits while sitting 44% from its
limit. The mechanism is exact rather than mysterious: step-doubling measures
`dE/d(log h)` of the error curve `E(h) = S(h) − S(0)`, so it reads zero wherever
that curve is **stationary**, and `S(h)` for Ti/Ni is flat to four digits from
`1e-3` to `1e-2`. A shipped default has no reason to avoid a stationary point of
its own error curve, and if it lands on one, every local self-consistency test
agrees with it.

The criterion is therefore **anchored** and not self-consistent: accept a step
only if the derivative is unchanged at `h/10` **and** `h/100`. Walking downward
fires; wiggling locally does not. Over all 21 pairs, at `H_DIFF` it catches
11 of 11 pairs that moved more than 5%, zero escapes, Ti/Ni included; at
`H_DIFF_CONVERGED = 1e-8` all 21 pass, worst residual 9.6e-05. The two
populations are separated by **2.82 decades** (9.60e-05 … 6.36e-02) and
`rtol = 1e-3` is placed inside that measured gap rather than chosen for
roundness — an odd thing to do on the day the repo learned that a threshold is
a default and a default is a claim.

`H_DIFF` itself is **not** changed. Section 6.13's numbers stay reproducible and
comparable, annotated in place, with the converged values beside them.

## 8. This audit's own validation committed the error the audit exists to find

Validation 4 asserted that a central difference of an exactly linear function
recovers the slope to better than `1e-14` relative — no truncation error at any
step — and it **measured 4.8e-11 and failed.** The arithmetic was right and the
assertion was wrong: cancellation leaves a rounding error of order `eps·|c₀|` in
the numerator, so the relative error carries a derived bound
`eps·|c₀| / (2h|c₁|)`, which is 7.7e-15 at `h = 1e-3` and 7.7e-10 at `h = 1e-8`.
The threshold `1e-14` was an **absolute constant that silently encoded
`h ≈ 1e-3`** — the audited fault class, one level up, inside the auditor. It was
caught by the test failing rather than by review, and it is recorded rather than
quietly corrected.

One validation is retained **because it cannot discriminate**. Bisecting an odd
function on a symmetric interval returns bitwise `0.0` after one evaluation at
any tolerance, for any interval width, with or without `rtol`: the answer is
exact and the test passes for the good default and the bad one alike. It proves
the harness and proves nothing about the default. Recording a
non-discriminating test *as* non-discriminating is cheaper than rediscovering
that it was — this repo has twice mistaken one for the other.

## 9. Where this leaves the series

* **09-20** — exact validation does not protect against an unrepresentative sample.
* **09-21** — pre-registration reaches what exact validation cannot.
* **09-22** — pre-registration does not reach the analysis layer.
* **09-23** — the other side of the ledger: how claims consolidate.
* **09-24** — a correct guard, correctly applied, silently degraded a correct
  number because its default tolerance encoded the scale of its original caller.
* **09-25** — **a procedure asked whether it has converged can answer yes and be
  wrong by 44%.** Self-consistency near a stationary point of the error curve is
  not evidence. Only an anchored comparison — against a step known independently
  to be in the plateau — distinguishes them.

Today's is the first class in the series where the **fix for the previous class
was itself the fault**: 09-24's conclusion was "report your convergence with
every value," which is what the first version of `sensitivity_converged()` did,
and which Ti/Ni defeats by a factor of 190 000.

There is also an ordering worth stating, because it was not obvious in advance.
The screen came first and cost nothing; the margin measurement came second and
cost seconds; the published-claim check came last and is the only one that
produced a correction. The screen alone would have flagged `H_DIFF` without
knowing whether to care. The margin study alone would have shown a bad step
without knowing which claim it reached. Neither alone is an audit.

## 10. Left open

* **`output_conductance(dVds=1e-3)` at `Vds = 0.05 V`** — the identical
  construction, a step 2% of its variable, found by the census rather than the
  probe. Deliberately **not** measured today: `g_ds` feeds `f_max`, a published
  number, and moving it needs its own before/after comparison. This is the
  direct successor to today's work and the top numerical item.

  > **ANNOTATED 2026-09-26 — this item is now CLOSED, and the "2%" above is
  > the ratio at the function's *signature*, not at the call site.** Measured
  > this session (`graphene_gds_quadrature_audit.py`,
  > `notes/2026-09-26-gds-step-quadrature-audit.md`): `dVds = 1e-3` sits
  > inside a broad anchored plateau (worst error 1.25e-6), the entangled
  > `n_segments = 50` quadrature bias is 8.4e-7, and peak f_max moves 1.1e-7
  > — four to five decades below anything that could touch a published
  > number. That audit's prediction Q2, which expected ≥ 0.1%, **FAILED on
  > its magnitude band.** Separately, Chapter 4's RF numbers are computed at
  > `Vds = 0.1 V` (`plot_fT_fmax`), where the ratio is 1e-2 and not the 2e-2
  > this census recorded from the signature default of 0.05 V: **a
  > default-scale census must read call sites, not signatures.** The 2e-2 is
  > left in place above — it is right about the signature and wrong about the
  > claim it was aimed at.
* **Every remaining `h`-like step and tolerance in `rf_small_signal_model.py`
  and the photodetector modules**, now that the instrument exists and takes
  minutes per default.
* **Whether any Chapter 4 or 5 ranking is a step artefact of the same kind.**
  The Section 6.13 ranking was, and the failure mode — a ranking ordering
  entries by how well a shared numerical choice suits each — is not specific to
  derivatives or to Chapter 6.
* **`log_sensitivity`'s own convergence estimate is step-doubling**, i.e. the
  pattern Section 7 shows to be blind. Its callers in
  `graphene_sensitivity_audit.py` (Chapter 4's `R_transmission` decomposition,
  Chapter 5's liner scenarios) have never been checked against an anchored
  criterion. That is a concrete, bounded, and now clearly motivated job.
