# Appendix B: A taxonomy of how this thesis's claims failed

Moved from Chapter 7 (§§7.5–7.7) on 2026-11-16 with section numbers kept as anchors. It is also the raw material for a separate methods paper; Chapter 7 keeps a one-paragraph summary.

<a id="s7-5"></a>

## 7.5 A taxonomy of how this thesis's claims failed

Five failure classes are on record, each with a *distinct* detector. The list
matters because the detectors are not interchangeable: the check that catches
one class is systematically blind to the others.

<a id="s7-5-1"></a>

### 7.5.1 Class (i) — implementation error inside the model

**Detector: an exactly known value.**

The 2026-09-17 transit-time integral produced `inf − inf` in a limit where the
answer was analytically known. Caught immediately by a validation against that
known value. This is the class every codebase checks for, and Chapters 4, 5 and
6 all check for it: the repository now carries dozens of exact validations,
several of them bitwise.

**But an exact validation can be structurally incapable of discriminating.**
2026-09-27 established this in its sharpest form. `validate_power_laws` in the
conditioning audit consists of four sub-checks, all of which are pure power laws
or closed forms. For a power law `ln|f|` is *linear* in `ln p`, so every secant
over it is exact: the validation passes identically for an estimator that was
provably wrong, for its corrected replacement, and for an absurd step size of
25 %. Its own output proves the point — switching to the demonstrably *more*
accurate estimator moved its "worst error" from 1.66e-11 to 3.79e-11. **A
validation whose number gets worse when the estimator gets better is measuring
rounding.** Three such non-discriminating checks are now retained deliberately,
labelled as such, because a labelled blind check is more useful than a deleted
one.

<a id="s7-5-2"></a>

### 7.5.2 Class (ii) — a claim quantified on an unrepresentative sample

**Detector: enumerating the full population. Nothing else works.**

Section 6.9's `<1.02` ceiling (§7.4.2) is the canonical case: every exact
validation of the model passed, because the model was right and the *claim about
it* was not. Section 6.11.2's count of five `>2×` mask-gain pairs is another —
the count runs 3-4-5-5-6 across Section 6.13's nominal uncertainty band, so the
number was never a property of the physics.

Pre-registration addresses this class and essentially only this class. It has
worked: on 2026-09-21 four predictions were committed to git before the model
was written, three held, and one was falsified at **17.12 %** against a
pre-registered ≤10 %. The falsification was informative *precisely because* a
number had been attached in advance. Had the prediction been "Cu will move
somewhat", nothing would have been learned.

**The conditioning audit is where this class was finally measured rather than
avoided.** Sections 4.10 and 5.5 asked whether Chapters 4 and 5 contain
near-cancellations, and the answer is yes, with two results that belong in this
chapter:

* **The worst-conditioned step in this thesis is in Chapter 5, not Chapter 6.**
  The λ_impurity calibration has `κ = 19.71`, against Chapter 4's worst of
  11.46 and Chapter 6's 4.55. See §7.7.
* **`κ` alone is the wrong instrument.** In Section 4.7 it ran *opposite* to the
  robustness of the conclusions drawn from it. Pd's positive residual — the only
  row consistent with the additive decomposition surviving — flips sign on a
  **9.6 %** error in `R_extra`, while Cu's large negative needs **81 %**. Large
  `κ` occurs where the *output* is small, so it threatens **signs, not orders of
  magnitude**, and the useful quantity is the sign-flip margin. The consequence
  is a genuine strengthening of a negative result: **the additive decomposition
  of Section 4.7 fails robustly**, and amplified input error is eliminated as the
  cause of the negative residuals.

<a id="s7-5-3"></a>

### 7.5.3 Class (iii) — an error in the analysis layer wrapped around a validated model

**Detector: hand-checking one row against the physics.**

On 2026-09-22 an unbracketed bisection returned its own endpoint and printed 21
plausible fake roots, while *all five* of that module's exact validations
passed — because they validated the model, and the fault was in the analysis
wrapped around it. The fix (2026-09-24) is a root-finder that raises on an
unbracketed interval; it fired on **21/21** deliberately root-free intervals
where the previous loop had returned a plausible number on every one.

The general lesson is that **a validated model is not a validated result**, and
the boundary between them is where this repository has been weakest.

<a id="s7-5-4"></a>

### 7.5.4 Class (iv) — created by a fix

**Detector: keeping the superseded implementation as an oracle.**

The 2026-09-24 guard was *correct*. Its default tolerance was not: `tol=1e-14`
was chosen for a variable of order 1 eV, where it is machine precision, and was
then applied to a decay length of order 5e-11 m, where the same absolute
constant stops the bisection after 16 halvings at a relative precision of
2e-4 — wrong in the fifth significant figure, with nothing raised.

> **A default tolerance is a claim about the scale of the caller's variable, and
> a guard moved to a new call site carries its defaults with it.**

Reviewing the guard for correctness would never have found this. It was caught
only because the unguarded loop had been returning the right answer for three
days and could be used as an oracle — which is why **deleting a superseded
implementation before validating against it destroys the only evidence that the
replacement is worse**, and why this repository now keeps old code paths behind
a flag rather than removing them.

<a id="s7-5-5"></a>

### 7.5.5 Class (v) — a procedure that certifies itself

**Detector: none that is reliable. Structural mitigation only.**

This is the class discovered most recently and the one that should worry a
reader most.

* **2026-09-25.** A convergence estimate built on step-doubling answered *yes*
  while the answer was wrong by **44 %**, blind by five orders of magnitude.
* **2026-09-26.** Refining one step converges to the derivative of whatever
  *other* discretisation was held fixed. When two numerical defaults are
  entangled — and they are entangled whenever one sets the other's grid —
  convergence in one is not evidence about the pair.
* **2026-09-27.** The conditioning work's own sensitivity estimator carries a
  step-doubling estimate `conv = |S(h) − S(2h)|` that is exactly 3× conservative
  where truncation dominates (measured 3.0000) and has **exact blind spots**
  where the `h²` and `h⁴` error terms cancel between `h` and `2h`. Four such
  points lie on the truncation branch of the Chapter 5 calibration — the
  `κ = 19.71` step — and at `h = 2.93e-2` the estimate reports **4.97e-14**
  while the sensitivity is wrong by **14.0 %**. Understatement factor 1.8e14.
  Worse, and more general: minimising that estimate chooses a step a median
  **877×** below the one that minimises the error, and gives a *worse* answer
  than the step in use at **18 of 23** call sites. **The instrument is smallest
  where the answer is worst.**

  The same session found that the estimator was never the central difference its
  own documentation claimed — it sampled `p(1 ± h)`, whose logarithms are not
  symmetric, making it a secant over an asymmetric log interval that returns the
  derivative at `ln p − h²/2`. This is the *same* fault class as the 09-26 guard
  (an asymmetric secant returns the derivative somewhere else, correctly
  computed), found there on a path that never fired and here on the default path
  behind every published sensitivity.

  **No published number in this thesis moves.** Every sensitivity in Sections
  4.10 and 5.5 is stable to 1.3e-8 relative under 16× refinement, against the
  corrected estimator, and against an independent reference; the nearest blind
  spot is ~3000× above the step in use. What is unsafe is the *procedure*, and
  the repository ran exactly that procedure twice in the preceding two days.

<a id="s7-5-6"></a>

### 7.5.6 The recurring class, and the only defence found

There is a sixth entry that is not a separate class so much as a demonstration
that knowing a class does not prevent committing it.

**Three consecutive sessions — 09-25, 09-26, 09-27 — each asserted a numerical
criterion as a round absolute constant that silently encoded a step size.** On
09-27 it happened **three times in one session**, in the session auditing
exactly that fault: three of five validations required deviations below 1e-13,
1e-14 and 1e-13 respectively, and all three measured their own derived
cancellation floor.

On 09-26 the conclusion drawn was "knowing the class demonstrably does not
prevent reproducing it, and what caught it both times was a test that reports a
number instead of asserting a verdict." Three more instances support something
stronger and less comfortable:

> **An absolute tolerance is the default way a numerical assertion gets
> written. Exhortation does not fix it. The only defence found so far is
> structural — report the measured value beside a *derived* bound and let the
> ratio be the verdict.**

Every validation in the repository is being migrated to that form. The failing
first forms are left in the source, because a record of what was asserted is
part of the result.

<a id="s7-5-7"></a>

### 7.5.7 The allocation argument

The five classes have five different detectors, and the detectors do not
overlap:

| class | detector | do Chapters 4–5 have it? |
|---|---|---|
| (i) model implementation | exact known value | **yes** |
| (ii) unrepresentative sample | full enumeration + pre-registration | **no** |
| (iii) analysis layer | hand-check a row against physics | **no** |
| (iv) created by a fix | keep the superseded path as oracle | now yes |
| (v) self-certifying procedure | derived bound beside the measurement | partly, as of 09-27 |

> **Verification effort has to be allocated by failure class, not uniformly.**

Sections 4.10 and 5.5 say *where* Chapters 4 and 5 are numerically fragile — the
calibration, `κ = 19.71`. The taxonomy says *which kind of check* would find the
fragility. Those are different questions and this thesis needed both.

---

<a id="s7-6"></a>

## 7.6 The other side of the ledger: how claims consolidate

Retractions dominate the record because they are dramatic. The opposite
movement happened twice and is more valuable, because a consolidated result is
cheaper to trust than a surveyed one.

**Consolidation 1 — a parity that absorbed three earlier results (2026-09-23).**
The identity "`N` is exactly odd in the differential offset and exactly even in
the common offset" (measured 2.3e-15) was *unplanned*. It:

* implies Section 6.13.5's sign-invariance result, upgrading a 12 621-point
  numerical scan over ±3 eV from evidence-over-an-interval to a statement valid
  without one;
* implies `τ* = W_B − W_A`, which had been derived separately;
* **excludes** two of that session's own pre-registered predictions rather than
  merely missing them — odd-in-`s` forces `|N|` to be symmetric about the flip,
  so the predictions were not wrong by a margin, they were impossible.

A prediction that a later symmetry *forbids* is a better outcome than one that
merely fails, and this is the only instance of it in the thesis.

**Consolidation 2 — an assertion converted to a proof (2026-09-24).** The
monotonicity of the `|ΔW|` shift in the screening length ℓ had been a docstring
claim since 09-21, and a bracketed bisection buys nothing without it: a
bisection on a non-monotone function can bracket a root and return the wrong one
of several. It is provable term by term. The sharper test is an exact
consequence rather than the proof: **Pt's equilibrium separation equals the
anchor separation exactly (3.30 Å), so Pt's term is bitwise 0.0 at all 801 grid
points and the maximum is really over Cu and Au alone.**

**What this implies for how the thesis should be read.**

> **A claim checked against a wide sample and a claim derived from a symmetry
> are not the same kind of knowledge.** This thesis has one and a half of the
> second kind and roughly a dozen of the first, and the ratio is the honest
> summary of its epistemic state.

---

<a id="s7-7"></a>

## 7.7 Unequal scrutiny

Three of Chapter 6's sessions overturned the previous session's headline, always
by widening a sample rather than by finding a physics error. Read naively, that
makes Chapter 6 the unreliable chapter.

The conditioning audit refutes the naive reading. The worst-conditioned step in
the thesis is Chapter 5's calibration at `κ = 19.71` — larger than anything in
Chapter 4 (11.46) or Chapter 6 (4.55). And today's audit found that the *only*
call sites in the repository where the convergence instrument has real blind
spots are, again, that same Chapter 5 calibration: four of them, at 13–15 %
error apiece.

> **The chapter that has failed most often publicly is not the chapter whose
> arithmetic is most fragile. It is the chapter that has been *looked at*
> most.**

Chapter 6 has had six sections' worth of samples widened, conventions removed,
and error bars attached — `±17 %` per metal and `±78 %` per same-sign pair from
the crossover assumption alone, stated explicitly at Section 6.8.7. Chapters 4
and 5 state general design rules from tabulated subsets and, before Sections
4.10 and 5.5, carried no error bars at all.

**A worked example of what scrutiny does, from Chapter 4.** Section 4.6 reported
"peak `f_T` ≈ 20 GHz" for four weeks. On 2026-09-26 the same quantity came out
at **10.140 GHz** — and *neither number was wrong*. The function's signature
default is `V_ds = 0.05 V`; the plotting path that produced the chapter's figure
passes `V_ds = 0.1 V`, giving 20.279 GHz. The chapter had simply never named
which. The consequence was systematic: an audit that had scored numerical
defaults by reading *signatures* had mis-scored one of them by exactly the
factor between the signature default and the call site — 2× in that case and in
general unbounded. **An unstated operating point is not a small omission; it is
a number that cannot be reproduced, and it survived a month of use.**

So the reliability claim this chapter makes is:

> **A result's reliability tracks (a) how wide a sample it was checked against,
> (b) whether its error bar was stated in advance, and (c) which of its terms
> cancel — not how carefully its physics was derived.**

All three of Chapter 6's headline retractions were detected by (a); the 17.12 %
falsification was informative because of (b); the robustness of Section 4.7's
negative result was established by (c).

---
