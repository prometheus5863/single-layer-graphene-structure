# Chapter 7: Discussion and Outlook — Where Single-Layer Graphene Realistically Fits, and What This Thesis Actually Established

> **Drafting note (2026-09-27).** This chapter's material has been complete for
> some weeks; it was planned and displaced on seven consecutive working
> sessions, each time by a numerical result that was worth the displacement.
> The ten synthesis threads it was handed are listed in Chapter 1's status
> table. This draft discharges them, and adds two that arose on the day it was
> finally written.

---

## 7.1 The question this chapter has to answer

Chapters 2 and 3 establish what single-layer graphene *is*: a gapless
two-dimensional semimetal with a linear dispersion near **K**, a density of
states vanishing linearly at the Dirac point, a quantum Hall sequence at
half-integer filling, and an optical absorption of `πα ≈ 2.3 %` per layer that
is independent of frequency and of every material parameter.

Chapters 4, 5 and 6 then ask what that *buys* in three device contexts a
semiconductor company would recognise: a field-effect transistor, an
interconnect, and a photodetector.

The honest answer that emerges is not a ranking of the three. It is that **the
three chapters are not asking questions of the same shape**, and that most of
the apparent contradictions between them dissolve once that is noticed — while
the one that does not dissolve is the most useful result in the thesis.

This chapter argues four things:

1. Chapters 4 and 5 optimise a **single interface**; Chapter 6's figure of
   merit is a **difference between two** interfaces. That is a structural
   split, not a disagreement, and it has a concrete industrial consequence
   (§7.2, §7.3).
2. A design rule expressed as a **bound** or a **parity** is a different kind
   of knowledge from one expressed as a **ranking over a tabulated set**, and
   this thesis spent most of its effort acquiring the second while believing it
   was acquiring the first (§7.4).
3. The claims in this thesis failed in a small number of recurring ways, each
   with its own detector, and **verification effort should be allocated by
   failure class rather than uniformly** (§7.5). The counterpart — how claims
   *consolidate* — is a shorter but more valuable list (§7.6).
4. The chapter that has publicly overturned its own headlines most often is
   **not** the chapter whose arithmetic is most fragile. It is the chapter that
   has been looked at most (§7.7). This is a claim about method, and it is the
   one a reader should carry away.

§7.8–§7.10 then state the application-space conclusions and what would change
them.

---

## 7.2 The durable tension: one junction versus the difference between two

### 7.2.1 How the tension was first stated, and why that version is retracted

On 2026-09-17 the two-contact collection model (Section 6.7) appeared to
produce a clean contradiction with Chapter 4: the metals Chapter 4 likes for
low contact resistance — Ni, Au, Pd, all high work function — were the *worst*
for two-terminal photoresponse. Chapter 1 recorded that as this chapter's
synthesis task.

**That version does not survive.** The signed, carrier-resolved model of
Section 6.8 places the n/p crossover at its physical value
`w_cross = W_G + Δ_c(d) ≈ 5.4 eV` rather than at graphene's own 4.5 eV, and at
the physical crossover Pd sits near the *top* of the symmetric ranking, not the
bottom. The symmetric ranking itself is crossover-dependent and Section 6.8.5
withdraws it. The metal-list contradiction was an artefact of a convention.

This chapter does not restate the retracted version, and Section 6.8.6
annotates it in place rather than deleting its numbers — the standing practice
of this work.

### 7.2.2 What survives, and it is sharper

The durable statement is structural:

> **Chapter 4 optimises a single junction. The two-terminal zero-bias
> photoresponse depends on the *difference* between two junctions, and is to
> first order blind to either junction's own quality.**

The evidence is quantitative. Writing `N` for the net zero-bias response of a
two-contact device, Section 6.8.4 finds `|N| = 1.832` for the asymmetric pair
Ti/Pt against **at best 0.556** for *any* symmetric pair — a factor of 3.3 in
favour of deliberately mismatching the contacts. The nonlinear
work-function-to-doping relation of Section 6.11 compresses this to 1.742
(−4.9 %) and does not touch its sign or its ordering.

Section 6.14 then makes the blindness exact rather than approximate. Decomposing
the two contacts' work-function offsets into a common part `m` and a
differential part `s`, with `ΔW_{A,B} = m ∓ s`:

    N is exactly ODD in s and exactly EVEN in m      (measured 2.3e-15)

`s` is what sets `sign(N)`. `m` — the part that a single-junction figure of
merit like contact resistance is sensitive to — enters only evenly, so it cannot
reverse the response and cannot appear at first order in the response's
magnitude about a symmetric pair. **A contact metal chosen for its own low
resistance is being chosen on a criterion the photoresponse is, by a parity, not
merely approximately but exactly insensitive to at leading order.**

### 7.2.3 The consequence an industrial reader acts on

Section 6.14 states the operational form, and it belongs here verbatim because
it is the kind of statement that changes a process decision:

> **At zero bias, a differential contact-chemistry error of unknown direction
> leaves the signal *magnitude* predictable and the signal *polarity*
> unpredictable.**

Since `N` is odd in `s`, `|N|` is symmetric about the flip: an unknown
differential error of size `|Δs|` gives a magnitude that is known to second
order and a sign that is a coin toss once `|Δs|` approaches the threshold. And
the threshold is exact — Section 6.14.2:

    τ* = W_B − W_A

The differential crossover offset needed to reverse a pair's photoresponse is
*exactly that pair's work-function gap*, independent of the common offset, of
the screening length λ, of the channel length L, and of every transport
parameter in the model.

For a differential readout, a wrong-sign channel is worse than a weak one: it
subtracts. **No single-junction figure of merit anywhere in Chapters 4 or 5 can
see this failure mode**, because all of them are computed for one interface at a
time. That is the concrete thing this thesis has to say to a process engineer
about two-terminal graphene photodetectors, and it is a consequence of a parity
rather than of a survey.

### 7.2.4 Where the risk is worst, and why the model cannot reach it

Section 6.14.5 records an anti-correlation that should be stated plainly rather
than buried in a caveat list. The four smallest reversal margins in the whole
table —

| pair | τ* |
|---|---|
| Au/Pd | **0.020 eV** |
| Ni/Au | 0.060 eV |
| Ni/Pd | 0.080 eV |
| Cr/Cu | 0.150 eV |

— belong to **exactly the four pairs the per-metal crossover model of Section
6.12 refuses**. Ti, Ni and Pd are chemisorbed and sit ~1.2 Å below the
physisorbed anchor, which is both why Section 6.12.4's anchored exponential
diverges on them (predicting `w_cross` of 6.25–62.6 eV, above the highest
elemental work function anywhere) and why their crossovers ought to differ most
from one another.

> **The model's reach and the thesis's risk are anti-correlated by
> construction.**

Au/Pd's 0.020 eV is therefore reported as *a threshold no available model can
reach*, not as a prediction of reversal. This is the single most important
caveat in Chapter 6 and it is a physics gap, not a numerical one: **a
description of Ti, Ni and Pd that does not go through their work function** has
been Chapter 6's central open problem since 2026-09-20 and remains open.

---

## 7.3 Two design rules pointing opposite ways, and Chapter 5's inherited framing

Section 6.8.6 produces two rules:

* **Symmetric device:** choose the metal that dopes graphene *least* — sit near
  the crossover.
* **Asymmetric device:** *maximise* `|ΔW_A − ΔW_B|`.

They point in opposite directions, and the asymmetric one wins outright (1.83
against at best 0.56), so the symmetric rule is largely a statement about what
to avoid. But the two rules are not in conflict once §7.2's split is granted:
the first minimises a quantity that enters `N` evenly, the second maximises the
quantity that enters oddly. They are rules about different variables that happen
to be written in the same units.

**Chapter 5 inherits Chapter 4's framing, and this has not been tested.** The
interconnect argument of Chapter 5 is a single-conductor argument throughout:
`ρ(W)` for one wire, compared against one scaled-copper baseline. Nothing in
Chapter 5 asks what happens when two graphene interconnects of nominally
identical specification are placed in a circuit that measures their
*difference* — a differential pair, a matched load, a current mirror. Given that
the analogous question in Chapter 6 produced the thesis's sharpest result, this
is the most obvious unexplored direction the synthesis exposes, and it is
recorded here as a prediction rather than a result:

> **Prediction (untested).** Any interconnect figure of merit that is a
> *difference* between two nominally identical graphene wires will be dominated
> by edge-scattering variance rather than by the mean resistivity that Chapter 5
> computes, because the edge term `((1−p)/(1+p))(λ_bulk/W)` carries the process
> variation and the bulk term does not.

---

## 7.4 Bounds and parities versus rankings over tabulated sets

### 7.4.1 The distinction

Chapter 6 produced exactly one result of each of three kinds, and the
distinction between them is the methodological core of this thesis:

| kind | example | what it rests on |
|---|---|---|
| **ranking** | "Ti/Pt is the best pair of the 21" | the seven-metal table |
| **bound** | `|N[g]| ≤ max|k|` for *every* illumination pattern at fixed photon number (§6.9.1) | the structure of the collection kernel |
| **parity** | `N` odd in `s`, even in `m`; hence `τ* = W_B − W_A` (§6.14.4) | a symmetry of the contact field |

A ranking over a tabulated set is contingent on the set. A bound holds for
every member of a continuum. A parity holds identically and *implies* results
that had previously been established by survey.

**Chapter 6 has six survey-derived results and one symmetry-derived one**, and
until Section 6.14 it had been accumulating the first while assuming it was
acquiring the second. Chapters 4 and 5 state general design rules — which
contact metal, which liner, which linewidth — from tabulated subsets, with no
error bars stated in advance and no bound or parity anywhere.

### 7.4.2 A bound is not safe either, if its sample was narrow

The cautionary case is Section 6.9's own ceiling. It stated
`max|k|/|N_uniform| < 1.02` for every asymmetric pair — a bound in *form*, and
therefore read as the strongest kind of claim. Section 6.11 found that **14 of
21 pairs violate it under Section 6.9's own linear model**, worst Au/Pd at
22.4×. The bound was never a bound; it was a ranking over the pairs that had
been enumerated, wearing a bound's clothing.

So the distinction of §7.4.1 is not "bounds are reliable". It is:

> **A claim is as strong as the population it was checked against, and the
> grammatical form of the claim is not evidence about that population.**

The design recommendation of Section 6.9 nevertheless survived the retraction,
on different grounds: a perfect mask gains 8.4× on Au/Pd, but 8.4× of a small
number is still 4.9× worse than unmasked Ti/Pt. **A recommendation can outlive
the argument that produced it, and saying so is more useful than quietly
re-deriving it.**

### 7.4.3 The synthesis question this poses to Chapters 4 and 5

> Which of Chapters 4 and 5's design rules could be restated as **parities or
> bounds** rather than as rankings over a tabulated set?

Three candidates are visible and none has been attempted:

* **Chapter 4's contact-resistance crossover.** The channel length at which
  contact resistance overtakes channel resistance is currently a per-metal
  number read off a table. It is a ratio of two terms with different
  `L`-scaling, so it should admit a statement of the form "for every metal with
  `R_c` in this range, the crossover lies in this interval" — a bound, not a
  list.
* **Chapter 5's liner scenarios.** The conclusion currently ranges from "no
  crossover" (bare 3 nm TaN/Co, thin-Ru liner) to "≈3 nm" (thin-Co liner), i.e.
  it is **materially liner-choice-dependent**. That spread is presented as an
  uncertainty; it is better read as evidence that the wrong quantity is being
  ranked, and that the robust statement is a *parity* in whether the liner
  conducts in parallel or not.
* **Chapter 4 §4.7's additive decomposition.** Its failure is the cleanest
  candidate for a bound: see §7.5.2 below, where the conditioning result turns a
  negative result into a robust one.

**Today's audit adds a fourth, and it is a caution rather than a candidate.**
The `log_sensitivity` result of 2026-09-27 shows that the *same*
near-cancellation which makes Chapter 5's calibration ill-conditioned also
creates blind spots in the instrument used to measure that ill-conditioning. A
bound derived from a near-cancelling quantity inherits the near-cancellation. So
"restate it as a bound" is not a way around conditioning; it is a way of making
the conditioning explicit.

---

## 7.5 A taxonomy of how this thesis's claims failed

Five failure classes are on record, each with a *distinct* detector. The list
matters because the detectors are not interchangeable: the check that catches
one class is systematically blind to the others.

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

## 7.8 Graphene's realistic near-term application space

With those caveats stated, the three device cases separate cleanly.

### 7.8.1 The logic transistor: not a near-term application, for a reason that is not fixable

Graphene has no bandgap, so a GFET has no off state: Chapter 4's transfer
characteristics show the ambipolar V-shape with on/off ratios of order 10, not
10⁶. Nothing in Chapters 4–6 changes this, and the standard mitigations
(nanoribbon confinement, bilayer displacement fields, substrate-induced gaps)
buy a gap at the cost of the mobility that motivated graphene in the first
place.

**The honest conclusion is that single-layer graphene is not a candidate for
digital logic and this thesis offers no argument that it is.** Chapter 4's value
is not a transistor recommendation; it is the contact physics (§7.8.2), which is
a prerequisite for every other application here.

### 7.8.1a RF/analog: the one application this thesis cannot speak to, and the reason is structural

Section 7.8.1 rules out digital logic on the bandgap. The usual next
sentence — and the one the GFET literature reaches for — is that graphene's
mobility makes it a candidate for *analog and RF* rather than digital, where
no off state is required and f_T and f_max are the figures of merit.
Chapter 4 Section 4.6 computes both. **This thesis is nevertheless not
entitled to an RF conclusion, and Section 4.6.2 is where that became
clear.**

The reason is not a tuning problem. Decomposing the f_max denominator
(2026-10-04) showed that 99.78 % of it is the term g_ds(R_g+R_s), that
g_ds ≈ 1/R_total to within 0.3–4.9 % across V_ds = 0.05–1 V, and therefore
that

```
  f_max/f_T  ≈  (1/2) sqrt( R_total / (R_g + R_s) )
```

to 0.36 %. The model's f_max is a restatement of the ratio between its
access resistance and its channel resistance. That is not a *bad estimate
of* f_max; it is a different quantity wearing f_max's name. The underlying
cause is that `transfer_characteristic()` is a resistor —
I_d = V_ds/(R_channel + R_c) — with no saturation mechanism, so the device
has no output resistance, and f_max is primarily a measurement of output
resistance.

**Three consequences, in increasing order of how much they cost this
thesis.**

First, the R_g work was aimed at a non-binding constraint. Multi-finger
gate-resistance reduction, the N² scaling, and the 40× delivery bug found
on 2026-10-03 all act on a term that cannot close the gap: with R_g = 0
*exactly* — a perfect gate, not a better one — f_max/f_T is 0.935276, still
28 % below the bottom of Feijoo *et al.*'s 1.3–1.4 band. Sixty-four fingers
reaches 0.927229, within 1 % of that ceiling, so the mechanism is
exhausted at about sixteen. This is a clean instance of Section 7.5's
taxonomy applied to *effort* rather than to claims: the work was correct and
the target was wrong, and nothing in the audit record could have said so,
because every check asked whether R_g was computed correctly and none asked
whether R_g mattered.

Second, the apparent agreement in Section 4.6 was never evidence. Peak f_T
scales ×19.8 over a ×20 drain-bias range (10.140 GHz at V_ds = 0.05 V to
200.440 GHz at 1 V), because g_m ∝ V_ds in a resistor model. The ≈20 GHz
that matched "the literature range for non-exotic gate lengths" is the value
at the bias `plot_fT_fmax()` passes. Asked at V_ds = 1 V the same model
matches record 200 GHz devices. **A model with a free parameter that moves
the answer by ×20 and a literature range spanning ×20 will agree with that
literature at some bias, and did.** This is the 2026-10-02 lesson —
agreement between two artefacts is silent about the world — in its most
expensive form so far, because here the second artefact was the published
literature rather than another file in this repository.

Third, and the only constructive part: the gap is now a number. The g_ds
that would put this device in the Feijoo band at its own R_g, R_s and C_gd
is 0.1965 mS/µm against the model's 0.9527 mS/µm — a factor of **4.8483**.
So "add current saturation" is no longer a direction, it is a target with a
tolerance, and a saturation-velocity term in `transfer_characteristic()`
can be accepted or rejected against it. Equally, the *sign* of the trend is
a free diagnostic: f_max/f_T falls with drain bias in this model
(0.6442 → 0.6218) where a saturating device's rises, so the first thing any
such term has to fix is a trend, not a magnitude.

**The judgement, then, differs from the one Section 7.8.1 reached about
logic, and the difference matters.** Digital logic is ruled out by physics
graphene does not have — a bandgap — and no model improvement changes that.
RF is not ruled out; it is *unaddressed*, because the instrument this thesis
built for it measures something else. Those are not the same verdict, and
the honest statement of the second is that **this thesis contributes nothing
to the RF case for graphene, for or against.** The contact physics of
§7.8.2 is a prerequisite for an RF device too — R_s enters the same
denominator and is 47 % of it here — so the Chapter 4 work is not wasted on
this application. But the figures of merit are not results, and Section
4.6's numbers should be read as what they are: a resistance ratio, computed
correctly, reported under the wrong name for five weeks.

#### 7.8.1b Annotation (2026-10-05): the saturation term was built, and it did not change the verdict — it sharpened the reason

Section 7.8.1a left the RF case *unaddressed* rather than refuted, and named
the one repair that could change that: a saturation term worth ×4.8483 in
g_ds. **That term has now been built and the verdict stands — but the reason
it stands is no longer "the model has no output resistance".**

What changed. With velocity saturation present and v_sat taken from
optical-phonon physics rather than fitted, the trend §7.8.1a called "the first
thing any such term has to fix" **is fixed**: f_max/f_T now rises with drain
bias (0.6632 → 1.3844 over V_ds = 0.05–1 V) where it previously fell. The
model is no longer a pure resistor, and the sign of its central RF diagnostic
is now right.

What did not change, and why it is a sharper statement than before. The same
term delivers only **1.1361** of the required 4.8483 — 23.4 % — and the
shortfall is not a weakness in the saturation law. At the literature-scale
geometry, R_c,total is 15.0 Ω of a ≈ 26 Ω device, so **roughly half of g_ds is
a contact resistance that no channel mechanism can reach.** Inverting the
question: criterion A would need v_sat = 6.50 × 10⁶ cm/s, *below* Dorgan
*et al.*'s measured 1–3 × 10⁷ cm/s band on SiO₂, and no phonon energy in the
physically available range (0.059–0.196 eV) gets past a factor of 1.2491.

**Three levers, three non-binding results.** The gate resistance is exhausted
(R_g = 0 exactly still falls 28 % short). The feedback capacitance is
negligible (C_gd deleted entirely moves the answer 0.11 %). Velocity
saturation is diluted by the contacts. Each was proposed in turn as the
missing piece of the f_max shortfall, each was built or bounded, and **all
three point at the access resistance.**

This converts §7.8.1a's verdict from a confession into a finding. The earlier
statement was that this thesis contributes nothing to the RF case because its
instrument measured the wrong quantity — a fact about the instrument. The
statement now available is a fact about the device: **in this model, at this
geometry, the RF figures of merit are contact-limited, and the channel physics
that the GFET literature treats as the interesting part is not where the
constraint sits.** That is a claim this thesis is entitled to make, because
it was reached by building the mechanism and measuring it rather than by
assuming it.

It also makes Section 7.8.2 load-bearing for an application it did not
previously claim. §7.8.1a already noted that R_s enters the f_max denominator
and carries 47 % of it, so the contact work was "not wasted" on RF. The
stronger version is now available: the contact work is the *only* lever left
standing for RF in this model, which means §7.8.2's negative result — that
work function alone does not describe Ti, Ni and Pd — is a limitation on the
RF case as much as on the contact case.

**What would overturn this annotation.** A single counterfactual, and it has
not been run: f_max/f_T at R_c = 0 *exactly*, in the saturated model, in the
way §4.6.2 ran R_g = 0 and R_s = 0. If a perfect contact also falls short of
1.3, then all four levers are exhausted and §7.8.1a's structural verdict
becomes final rather than provisional. If a perfect contact reaches the band,
then the RF case for graphene in this model is a contact-engineering problem
with a quantified target, which is a far more useful conclusion than either
of the two this chapter currently carries. **That counterfactual is the top
open item created by this annotation**, and it is deliberately not asserted
either way here.

See Chapter 4 §4.6.3 and
`notes/2026-10-05-saturation-arrived-and-the-contacts-ate-it.md`.

#### 7.8.1c Annotation (2026-10-06): the counterfactual was run, and it answered both branches at once

§7.8.1b closed by naming one counterfactual as the thing that would overturn
it — `f_max/f_T` at `R_c = 0` *exactly*, in the saturated model — and stated
its two branches in advance. **It has now been run, and the answer is not
either branch. It is both, because the question was asked about a ratio.**

At the 40 µm / 8-finger geometry and `V_ds = 0.1 V`, with the contact removed
exactly and coherently:

| quantity | baseline | `R_c = 0` | change |
|---|---|---|---|
| `f_max/f_T` | 0.683186 | 0.759363 | +11.15 % — **58.4 % of the 1.3 floor** |
| `f_T` | 20.740 GHz | 74.650 GHz | **×3.5994** |
| `f_max` | 14.169 GHz | 56.687 GHz | **×4.0007** |

**On the ratio, §7.8.1b's first branch holds and the verdict is now final.** A
*perfect* contact — not a better one — leaves `f_max/f_T` 41.6 % below the
bottom of Feijoo *et al.*'s band. Four levers have now been built or bounded
and all four are non-binding: `R_g = 0` exactly gives 0.935276, deleting `C_gd`
entirely moves the answer 0.11 %, velocity saturation from optical-phonon
physics delivers ×1.1361 of a required ×4.8483, and `R_c = 0` exactly gives
0.759363. §7.8.1a's structural verdict on the *ratio* is no longer
provisional. A residual `g_ds` requirement of **×2.9955** survives the perfect
contact, and §7.8.1b already showed no phonon energy in the physically
available 0.059–0.196 eV reaches even ×1.2491 of it.

**On `f_max` itself, §7.8.1b's second branch also holds, and this is the part
an industrial reader acts on.** The same counterfactual is worth a factor of
**four** in `f_max`, and most of it is reachable: at the lowest contact
resistance reported in the literature (65 Ω·µm per contact, against this
model's 300) `f_max` rises ×2.79, which is **59.6 %** of the perfect-contact
gain. So the RF case in this model *is* a contact-engineering problem with a
quantified target — ×4.0 available, ×2.8 reachable — and that statement
coexists with the final verdict on the ratio rather than replacing it.

**Why both branches could be true, and it is not a drafting failure.**
`f_max/f_T` divides out exactly the quantity the contact dominates: the
contact limits `I_d`, `I_d` sets `g_m`, `g_m` sets `f_T`, and the same contact
sets `R_s` in the `f_max` denominator. Removing it moves numerator and
denominator together and the ratio keeps only the residue. **A criterion
written on a ratio is silent about any mechanism that scales both of its
arguments.** This is checkable rather than rhetorical: `f_max/f_T` is
**non-monotonic** in `R_c`, bottoming at 0.680047 at 470 Ω·µm and rising again
to 0.698615 at 4000 Ω·µm, while `f_max` falls **37×** across those same two
rows. There is a region of this design space where satisfying a ratio
criterion means losing on the figure of merit.

**One objection is now closed.** This model's contact is 300 Ω·µm per contact,
**1.8× the 165 Ω·µm** Feijoo *et al.* report for the very devices whose
1.3–1.4 band Chapter 4 is measured against. The whole `f_max` thread could
therefore have been nothing but a worse contact. Given Feijoo's own contact
resistance, this model reaches **0.695514**, not 1.3.

**Two findings consolidate, and one claim is narrowed.** §7.8.1b's dilution
*explanation* — that the contacts absorb the saturation — made an untested
prediction, that the `g_ds` saturation factor at `R_c = 0` must be the
undiluted 1.699. Measured: **1.70661**, i.e. 100.45 %, from a derivative of a
solved terminal current rather than from the channel integral that produced
the 1.699. And §7.8.1a's resistance-ratio identity acquires a sign test: the
*superseded* resistor model is made **worse** by the same perfect contact
(0.642901 → 0.581707), because that identity carries `R_total` in the
numerator and therefore rewards a worse contact. A figure of merit that
improves when the device gets worse is not measuring the device, which is what
§7.8.1a asserted and now has evidence for.

**What this annotation does NOT claim, and the reason is a literature result
that cuts against §7.8.1b's framing.** Feijoo, Pasadas, Bonmann *et al.*,
*Nanoscale Advances* **2** (2020) — a paper whose title asks "Does carrier
velocity saturation help to enhance `f_max` in graphene field-effect
transistors?" — answer their own question in the negative: "the largest
`f_max` are located at biases close to the onset of bipolar conduction and
**far from the saturated velocity regime**", with drift at only ~45 % of
`v_sat` and the **diffusion** contribution to the current comparable to drift
at the peak-`f_max` bias, concluding that "our results do not support that
operating in the regime of velocity saturation results in the highest
`f_max`". Three things follow, and none of them is comfortable:

1. **§7.8.1b's criterion A was constructed from this model's algebra, not from
   device physics.** Its arithmetic stands — ×4.8483 required, ×1.1361
   delivered — but it licenses a statement about *this model* and not about
   real GFETs, and that narrowing is recorded here rather than left implied.
2. **This model has no diffusion current at all.** Eq. (4) of the saturated
   model is drift-only. If diffusion is ~40 % of the current near the drain at
   the peak-`f_max` bias in a measured device, a drift-only model cannot be
   asked about that bias, and the `V_ds` ladders of §7.8.1b and of this
   annotation are a drift-only slice of a two-mechanism problem. **This is now
   the top open item for the RF thread**, ahead of anything else in it.
3. **This model has no self-heating**, a mechanism the same paper measures at
   65 → 40 GHz in `f_max`. Today's ×4.0007 is a ratio between two runs of one
   model and is largely insulated from this; the absolute 56.687 GHz is not.

So the RF verdict of §7.8.1a — *unaddressed*, not refuted — survives in a
third form. The ratio is settled and contact-limited only mildly; `f_max`
itself is contact-limited by a factor of four; and the reason this thesis
still cannot speak to the RF case is no longer the absence of saturation, but
the absence of diffusion and of self-heating in the model that was built to
replace the resistor.

See `graphene_perfect_contact_counterfactual.py`,
`perfect_contact_counterfactual_output.txt` and
`notes/2026-10-06-the-perfect-contact-and-the-ratio-that-divided-out-the-prize.md`.

#### 7.8.1d Annotation (2026-10-07): the diffusion term was built, and the item's premise did not survive it

§7.9 item 7, created 2026-10-06 and named that session's top item, said that
Eq. (4) is drift-only and that every `V_ds` ladder in §7.8.1b and §7.8.1c is
therefore a drift-only slice. The term has been derived and built
(§4.6.5, `graphene_diffusion_current_model.py`, 30 checks, 30 passed). Three
parts of the item are now settled and they do not settle the same way.

**It is worth ±0.50 % of `I_d`, and the sign flips between branches.**
+0.5018 % at the saturated-peak bias and −0.5038 % at the resistor-peak bias.
`f_max` gains +0.2941 % while `f_max/f_T` *falls* 0.2514 %, from 0.683186 to
0.681468. §7.8.1a's structural verdict is untouched — 52.42 % of 1.3 either
way — and the ratio moving against `f_max` is §7.8.1c's finding in a third
mechanism.

**The premise is the part that moved.** Pasadas and Jiménez
(*IEEE TED* 2016, [arXiv:1605.08235](https://arxiv.org/pdf/1605.08235)), from
the same group as the Feijoo *et al.* paper that created this item, write the
GFET channel current as `I = −W Q_tot v` with `v = μF`, `F = −dV/dx`, and
state explicitly that **"V(x) is the quasi-Fermi level along the graphene
channel"**. A drift expression whose driving potential is the quasi-Fermi
level is *already* the complete drift-diffusion current. Under that reading,
the term added in §4.6.5 is a double count, not missing physics.

`carrier_density(V_g, V_ch)` contains evidence for both readings and settles
neither. `V_ch` enters an electrostatic charge relation
`C_ox(V_g − V_ch − V_dirac)`, which reads electrostatic; it is then corrected
by a quantum-capacitance series factor, which is exactly what one applies when
`V_ch` is the quasi-Fermi level and the graphene drop `E_F/e` is a separate
voltage. **The two readings differ by exactly `λV_F`** — the term itself. So
±0.50 % is simultaneously the diffusion current under one reading and the size
of an ambiguity that was already present in every `I_d` this thesis has
committed.

This is §7.5's recurring class in a new location. 2026-10-04 recorded that a
quantity can be computed correctly under a name that claims more than it
measures; here the name is a **variable**, and what it fails to pin down is
which of two potentials it is. No value is wrong. The equation is not
determined.

**What this changes in this chapter.** The claim that §7.8.1b's and §7.8.1c's
ladders are *drift-only* slices is withdrawn: they are slices of a model whose
driving potential was never named, bounded at ±0.50 %. The repair is not the
one item 7 asked for — **it is to name the potential, not to add a term** —
and the recommendation is to declare `V_ch` the quasi-Fermi potential, keep
`λ = 0` as the production path, and keep §4.6.5 as the bound on the choice.

**And item 7's other claim is upheld for a different reason.** The
peak-`f_max` bias remains a bias this thesis cannot be asked about. Not for
the missing term: §4.6.5 shows the diffusion-dominated region near charge
neutrality is *absent* from this model, because the `n_puddle` regularisation
makes `d|n|/dV` exactly zero where Zebrev's `C_ox/C_Q` diverges. A
regularisation chosen to keep the conductivity minimum finite removed the one
regime the literature says matters, and no amount of care about the transport
term recovers it.

See §4.6.5 and
`notes/2026-10-07-the-diffusion-term-was-a-boundary-term-and-the-potential-was-never-named.md`.

#### 7.8.1e Annotation (2026-10-08): the potential is named, and the naming question was the smaller half

§7.8.1d recommended declaring `V_ch` the quasi-Fermi potential, keeping
`λ = 0` as the production path, and keeping §4.6.5 as the bound on the choice.
**All three have been done** (§4.6.6). The declaration lives in
`graphene_fet_model.py` as a machine-readable in-source marker, so the
statement of what this thesis's own variable means travels with the variable.

**What §7.8.1d got right.** The term is a double count, not missing physics:
+0.6965 % intrinsic and **+0.4279 % at the terminals** at `V_g` = 2.0 V,
`V_ds` = 0.1 V, bracketing §4.6.5's +0.5018 % at its own bias. The claim
withdrawn in §7.8.1d — that §7.8.1b's and §7.8.1c's `V_ds` ladders are
*drift-only* slices — can now be restated positively: **they are slices of the
complete drift-diffusion current in the quasi-Fermi variable**, with no
±0.50 % ambiguity attached. That ambiguity is resolved rather than bounded, and
§7.8.1a's structural verdict is untouched for the third consecutive session.

**What §7.8.1d did not see.** A declaration is a promise that the rest of the
model belongs to the declared variable, and that promise has a price. Under the
declaration the exact charge relation is
`V_g − V_dirac − V_ch = e·n/C_ox + E_F(n)/e`; §4.3 computes its linearisation.
Running Eq. (4) on the exact form moves `I_d` at the RF bias by **−0.6023 %**
— **1.4 times the term §7.9 items 7 and 8 were competing over.**

This is the uncomfortable part, and it belongs in this chapter rather than only
in Chapter 4. §7.5.6 identified the recurring class and §7.7 identified unequal
scrutiny as its cause. Here is a clean instance of the latter: **two sessions
(2026-10-06 and 2026-10-07) went into a transport term worth +0.43 %, while a
−0.60 % question about the same chapter's charge relation sat directly
underneath it, reachable by replacing one function and re-running.** Nothing
protected the charge model from scrutiny except that no item pointed at it, and
no item pointed at it because the attention had been captured by a literature
paper's framing (Feijoo *et al.*'s diffusion measurement) rather than by an
audit of this thesis's own equations. §7.7's thesis survives contact with a
third example.

**What it does not change.** Every number in Chapters 4 and 5 remains a
linearised-series-factor number, annotated in place and not deleted, per the
standing rule. The −0.60 % is a *measurement of a gap*, not a new result set,
and the two charge forms have been compared only to each other — never to
measurement — so nothing here licenses calling the committed model wrong, only
non-exact under its own declaration. §7.8.5's summary judgement is unchanged.

**And one ranking change that follows.** §7.9 item 9 (a signed-carrier FET
model, created 2026-10-07) is promoted from a sibling of item 7 to the
**prerequisite** of the new top item 12. Both charge forms carry `n` as a
magnitude, so re-deriving Chapter 4 on the exact relation before a
signed-carrier treatment exists would rebuild the branch-sign flaw in a new
equation rather than fix it. The hole branch is now blocking, not adjacent.

See §4.6.6 and
`notes/2026-10-08-the-declaration-and-the-question-underneath-it.md`.

#### 7.8.1f Annotation (2026-10-09): the prerequisite is cleared, and the off-state number this thread assumed was 2.07x off

§7.8.1e ended by promoting §7.9 item 9 to the prerequisite of the RF thread's
top item. **Item 9 is closed** (§4.6.7), item 12 is unblocked, and for the first
time in three sessions the top item has nothing in front of it. That is the
whole of this annotation's RF content: no `f_T` or `f_max` number moved today,
and the saturated-model queue — item 12, then the §4.6.5 rewiring, then
self-heating — is unchanged in order and one position shorter at the front.

**The non-RF finding is the one that bears on §7.8.1 rather than §7.8.1a.** This
thread's verdict that graphene is not a logic transistor rests on the off-state
floor, and §4.6.7 finds that Chapter 4's floor and Chapter 3's measured floor
disagree by **2.0678×**, with Chapter 4 on the more conductive side. Correcting
to Chapter 3's floor raises the on/off ratio from **1.2408 to 1.7820** at
`V_ds` = 0.1 V. **The verdict is undisturbed and this is worth saying plainly
rather than quietly**: an on/off ratio of 1.78 is no more a switch than one of
1.24, and §7.8.1's reason — that the floor is physics and not a process
limitation — is if anything strengthened by the discovery that the floor was
being *understated*. A finding that moves a number by a factor of two without
moving the conclusion it supports is the useful kind of robustness result, and
§7.6 is where it belongs.

But it is also an instance of §7.7. The floor comparison was a one-line
calculation flagged in §4.9 on 2026-09-29 and left for **eleven sessions**,
while two consecutive sessions went into a ±0.5 % transport term in the RF
thread. §7.7's unequal-scrutiny finding now has a fourth instance, and the
pattern across all four is the same: **the cheap check on the quantity the
conclusion actually rests on loses to the expensive check on the quantity the
current thread is about.**

#### 7.8.1g Annotation (2026-10-10): the fifth lever, and it was two levers that cancelled

§7.9 item 12 — re-deriving Chapter 4 on the exact charge relation — is closed
(§4.6.8). It is the **fifth** lever this thread has tried against the RF
verdict, after velocity saturation (10-05), the perfect contact (10-06), the
diffusion term (10-07/10-08) and the charge relation's prerequisite (10-09).
The verdict is unchanged and the *reason* it is unchanged has changed shape
twice.

**What moved.** `f_max/f_T` moves **+0.0228 %**, from 0.643007 to 0.643154 at
the 40 µm 8-finger geometry. Peak `f_T` moves −0.0594 %, peak `f_max`
−0.0366 %, extrinsic `f_T` −1.6720 %, extrinsic `f_max` −1.7091 %. The
structural statement of §7.8.1a — that this model's `f_max` is a restatement
of `R_total/(R_g+R_s)` and not an estimate of `f_max` — is therefore
untouched, which is the correct thing to have checked. An earlier draft of
that check tested a *frequency* instead (that `f_T` sat "two orders below the
100–300 GHz RF-frontend requirement") and was wrong twice: the gap is about
one order, and §7.8.1a never rested on a frequency threshold. §7.5's taxonomy
gains a small class (v) instance in the audit layer rather than in the model.

**What the small number conceals, and this is the new thing.** The −0.0594 %
in `f_T` is **a cancellation, not a smallness**. The charge change alone is
−1.9986 %; the gate-capacitance change that the same differential identity
forces is +1.9392 %; the two are within 3 % of each other in magnitude and
opposite in sign. Four of the five levers in this thread were reported as
single numbers. This one could not have been: reporting the charge alone would
have overstated the effect by 34×, and reporting only the total would have
called a pair of 2 % corrections negligible. **§7.4's distinction between a
bound and a ranking has a third member — a number that is small because two
larger numbers oppose — and nothing in this thesis's audit record would detect
it, because every check that has ever been written here asks about one
quantity at a time.**

**And the peak bias is where the sensitivity lives.** Peak `f_T` moves
0.0594 % in height and **6.76 %** in location, a ratio of 114. Both peaks stay
on the hole branch, so §4.6.7's open "does the peak-`f_T` bias change branch
between the two models" question is answered **NO** for the
(4.28)-vs-(4.29) pair, and stays open for the branched-vs-magnitude pair it
was actually asked about. The observation that this thread's figures of merit
are far more sensitive in *where* they peak than in *how high* they peak is new
and is not something any of the five levers was designed to reveal.

**What this does to the ledger.** The RF thread has now exhausted its
model-side levers; what remains (item 16, the `C_q` argument) is the first
open item with a cost larger than any lever tried, and it is a defect rather
than a missing physical effect. §7.8.1a's conclusion is unchanged and better
supported than it was, because the levers that could have moved it have been
tried and measured rather than assumed to be small.

### 7.8.2 Contacts: the real near-term contribution, and a negative result worth having

Contact resistance is the binding constraint on every graphene device in this
thesis, and Chapter 4's substantive results are about it:

* Contact resistance is **metal-specific for a structural reason**, not an
  empirical one: the metal's work function sets a doping profile that extends
  into the channel over a screening length, and the resulting junction
  resistance is what a two-probe measurement attributes to "the contact"
  (Section 4.5).
* **The additive decomposition `R_c = R_transmission + R_extra` fails for three
  of four metals, and it fails robustly.** The negative residuals are not
  amplified input error: Cu's would need an 81 % error in `R_extra` to flip
  sign (Section 4.10.3). This is a real modelling result and it says the
  two-term picture is incomplete, not mis-parameterised.
* **Edge contacts beat top contacts** and patterned (hole-array) contacts
  reproduce the qualitative direction and large-diameter trend of the published
  data — but not the small-diameter upturn, and not the full measured ~11×
  device-level reduction (Section 4.8).

The unresolved literature discrepancy should be stated here because it bounds
everything in Chapter 6: the *measured* potential step at a graphene/electrode
interface is ≈ **0.12 eV** with doping extending 0.2–0.3 μm into the channel,
against the **0.25–1.07 eV** offsets the metal table assumes. The measurement is
of the residual step in a *gated* device rather than the flat-band charge
transfer, which is the likely reconciliation, but it has not been done. Until it
is, every ΔW-derived number in Chapter 6 carries an unquantified systematic on
top of its stated ±17 %.

### 7.8.3 Interconnects: conditional, and the condition is a liner choice

Chapter 5's conclusion is genuinely conditional and should not be compressed. A
graphene interconnect beats scaled copper below some crossover linewidth, and
that crossover ranges from **"no crossover at all"** (bare 3 nm TaN/Co liner, or
a thin Ru liner) to **≈3 nm** (thin Co liner) depending on which copper baseline
is assumed.

> **The crossover is materially liner-choice-dependent, not a single number, and
> a paper quoting one figure has chosen a baseline.**

Two things temper this further. The `λ_impurity` calibration behind `ρ(W)` is
the worst-conditioned step in the thesis (`κ = 19.71`), and it propagates into
`ρ(W)` damped but not removed (`κ` from 0.88 to 1.55 over 18–52 nm). And the
liner model's `W_eff = W − 2t` contains an unremarked subtraction carrying a
~20 % bar on `ρ_eff` at 18 nm that Chapter 5 had not stated before Section 5.5.4.

The exact result of Section 5.5.3 is worth carrying because it is a
*structural* statement in a chapter otherwise short of them: **at the
calibration width, three of the five parameter sensitivities are exactly zero
and `S(ρ_bulk)` changes sign through it.** The model there carries no
information from `ρ_bulk`, `p` or `λ_bulk` at all — it is pinned to the datum.
Any claim about the interconnect that is evaluated near the calibration width is
a claim about the calibration, not about the physics.

*(Today's audit refines one detail of that result: the third of those three
zeros was bitwise-exact partly by floating-point luck, and is now stated as
"zero to the derived cancellation floor" rather than "exactly zero". The
structural fact is unaffected.)*

### 7.8.4 Photodetectors: the strongest case, for a reason specific to graphene

This is where graphene's properties line up rather than fight:

* The gaplessness that kills the logic transistor makes the detector
  **broadband** — absorption is `πα` independent of frequency.
* `2.3 %` absorption per layer is the bottleneck, and it is addressable by
  **plasmonic near-field enhancement** (Section 6.5 reproduces 25× and 8.5× for
  two published designs) rather than by changing the material.
* The gain–bandwidth tradeoff is explicit and quantified (Section 6.4), so the
  design space is a curve rather than a wish.
* **Asymmetric metallisation beats illumination engineering** by ≈4× per
  incident photon (≈8× against the shadow mask actually built in the
  literature), and the reason is structural: illumination enters only as a
  weight on a fixed collection kernel, so `max|k|` caps what *any* pattern can
  achieve (Section 6.9.1). Masks are worth using where metallisation cannot be
  changed.

The caveat that matters most is §7.2.4's: the pairs with the largest predicted
response are also the pairs whose reversal margins the model cannot compute.

### 7.8.5 Summary judgement

| application | verdict | binding constraint |
|---|---|---|
| Digital logic | **no** | no bandgap; not fixable without losing the mobility |
| RF / analogue | **unaddressed** (§7.8.1a–c) | the ratio `f_max/f_T` is settled below the literature band at `R_c = 0` *exactly*; `f_max` itself is contact-limited ×4.0, of which ×2.8 is reachable; the model lacks diffusion current and self-heating |
| Interconnect | **conditional on the copper baseline** | liner choice decides whether a crossover exists |
| Photodetector | **strongest near-term case** | 2.3 % absorption, addressable; polarity risk from contact chemistry |
| Contacts (as a topic) | **the enabling contribution** | the two-term decomposition is incomplete |

---

## 7.9 What would change these conclusions

Stated as falsifiable items, in order of how much they would move:

1. **A second anchor for `Δ_c` at any separation other than 3.3 Å.** Every
   "reachable" number in Section 6.12's margin table is one anchored
   exponential with a swept decay length. A second anchor would make Section
   6.12's per-metal crossover a measurement instead of an interpolation, and it
   is the top *physics* item in this thesis — open since 2026-09-21 and
   untouched for seven consecutive working sessions while the numerical thread
   ran.
2. **A description of Ti, Ni and Pd that does not go through work function.**
   Three independent failures are on record for these chemisorbed metals: the
   nonlinear doping relation excludes them (`d_eq < d_0`), the anchored
   exponential diverges on them, and they carry the four smallest reversal
   margins. Ti has the largest `|ΔW|` in the table (−1.07 eV) and is contact A
   of *both* Chapter 6 headline pairs. The chapter's central numbers rest on the
   one assumption the literature most explicitly disowns.
3. **Reconciling the 0.12 eV measured step with the 0.25–1.07 eV assumed
   offsets** (§7.8.2). Bounded in the common channel by Section 6.13's
   sensitivity table; **not** bounded in the differential channel, which is the
   one that sets polarity.
4. **A photo-thermoelectric term.** The parity of Section 6.14 is derived for a
   drift-only model with `μ_e = μ_h`. Whether it survives a
   photo-thermoelectric contribution is the single most consequential open
   question about the thesis's one symmetry-derived result — because if it does
   not, §7.2.3's polarity statement loses its exactness.
5. **Asking the differential question of Chapter 5** (§7.3's prediction).
6. **Re-running Chapter 4's contact-resistance results at the 5.4 eV
   crossover** rather than at graphene's 4.5 eV, for consistency with Chapter 6.
7. ~~**Naming the potential in Eq. (4)**~~ — **CLOSED 2026-10-08** (§4.6.6,
   §7.8.1e). `V_ch` is declared the **quasi-Fermi potential**, in the source
   where it is defined. The consequence is that 2026-10-07's diffusion term is
   a **double count** worth +0.4279 % of `I_d` at the terminals, not a missing
   0.50 %, so the committed `λ = 0` currents stand and `λ` is retained as the
   instrument that measures the difference between readings. The item's
   history — created 2026-10-06 as *a diffusion term in the saturated
   transport model*, rewritten 2026-10-07 when the term turned out to be the
   difference between two readings of a variable — is kept because the shape of
   it recurs. **It is replaced, not succeeded, by item 12.**
8. **A self-heating term.** The same paper measures 65 → 40 GHz in `f_max` from
   self-heating alone. This bounds every absolute frequency in Chapter 4; the
   *ratios* between model runs, which is what §7.8.1c reports, are much less
   exposed.
9. ~~**A signed-carrier FET model**~~ (added 2026-10-07) — **CLOSED
   2026-10-09**, `graphene_fet_signed_carrier_model.py`, 19/19 checks, §4.6.7.
   The sharper statement of why it blocked item 12 than 10-08 gave: both terms
   on the right of (4.28) are **odd** in `n`, so the relation has an odd root,
   and an odd root cannot come out of a solver whose first line is
   `dV = np.abs(dV)`. `carrier_density()` *computes* the sign — `n_eff` is odd
   in `dV` because `C_q` is even — and discards it in the quadrature floor.
   `n_net_exact()` carries it: odd **bitwise**, satisfying (4.28) to 7.690e-16
   relative. The decomposition `n_e,h = (S ± n_net)/2` reproduces the shipped
   magnitude as an identity, yields `n_e·n_h = n_puddle²/4` (so the August
   quadrature floor is algebraically a mass-action law), and makes graphene's
   ambipolar symmetry a bitwise check. **No number in Chapter 4 moved** — 3 ulp
   on `σ`. The two sub-questions this item carried are now separated: the
   hole-branch sign of Chapter 4's branch-asymmetric results is *expressible*
   but still unverified, and the peak-`f_T` branch-crossing question is
   untouched; both are retained in the backlog under their own entries rather
   than inside a closed item.
10. ~~**A census of every potential-like variable in this thesis**~~ —
   **CLOSED 2026-10-08** as a standing instrument,
   `graphene_potential_census_audit.py`. 86 names across 25 modules, each
   carrying an explicit adjudication; a name present in the repository and
   absent from the registry is reported as a fault, so the census cannot go
   stale silently — and it caught its own staleness on the very next file
   added. **One defect in 86 names.** All three siblings this item predicted
   were checked and none is one, for three different reasons: Chapter 5's
   interconnect model contains *no potential-like variable at all*; Chapter 6's
   two-contact `V_bias` is imported and printed and enters no computation in
   its own module; §4.4's `V_channel_shift` appears only in a docstring
   derivation the evaluated closed form never uses. The first of those is kept
   as a check that is *expected to fail*, because a pass there would mean a
   potential had appeared in Chapter 5. Seven names remain `UNADJUDICATED` and
   are recorded as open work; none is channel-swept, so none can carry the
   defect. **A predicted class is not a measured class**, and the measured rate
   is 1 in 86.
11. **Whether the `n_puddle` regularisation removes other regimes, not only
   the diffusion-dominated one** (added 2026-10-07). §4.6.5 found that
   `d|n|/dV` is *exactly* zero at `n_eff = 0`, which deletes the
   charge-neutrality physics that Zebrev's `C_ox/C_Q` and Feijoo *et al.*'s
   peak-`f_max` bias both live in. A floor added to keep one quantity finite
   silenced a different one. Which other near-Dirac results in Chapters 4
   and 6 rest on a derivative the floor has set to zero is not known, and this
   is a sibling of the §4.9 item on whether `n_puddle` reproduces §3.3's
   6.45 kΩ/sq floor — open since 2026-09-29 and untouched for ten sessions.
   **2026-10-08 adds a reason to prefer this item over item 13:** §4.6.6's
   exact charge relation (4.28) and the series factor (4.29) disagree by +68 %
   raw at `dV` = 0.01 V and by 0.001 % once the floor is applied, so the floor
   is now known to be the dominant term in the near-Dirac charge of *both*
   forms. Whatever it has silenced, it has silenced in the replacement model
   too.
   **2026-10-09 — PARTIALLY ANSWERED, and the premise of the question was
   wrong.** The question asked which *other near-Dirac regimes* the floor had
   removed. §4.6.7 measures that the floor is not a near-Dirac term at all: at
   `V_g` = 3.5 V, the largest overdrive swept anywhere in this thesis, it is
   **0.774** of the gate-induced net density and adds **+26.4 %** to the total,
   and its contribution falls below 1 % only at `dV` = **14.73 V**, 5.46× the
   maximum swept. It does not remove regimes; it adds a parallel conducting
   density to all of them, and `carrier_density()`'s docstring calling it a
   Dirac-point regularizer describes a regime this thesis never enters. It was
   found by a check written to PASS on that docstring's authority, which failed
   at −6.92 % in `I_on` — an allegedly off-state-only parameter moving the on
   state. **What remains open is the original `d|n|/dV` = 0 question**, which
   is about the derivative and is not answered by this. The sibling §4.9 floor
   item is **CLOSED** (see item 11a).
11a. ~~**Whether Chapter 4's `n_puddle` reproduces §3.3's 6.45 kΩ/sq floor**~~
   (§4.9, added 2026-09-29, **untouched for eleven sessions**) — **CLOSED
   2026-10-09**, §4.6.7 and §3.3.1. It does not. `σ_min = e·μ·n_puddle` =
   3.204353e-04 S/sq = **3.1208 kΩ/sq** against §3.3's measured `4q_e²/h` =
   **6.4532 kΩ/sq**: Chapter 4's channel at neutrality is **2.0678× more
   conductive** than the floor Chapter 3 calls measured, and the reconciling
   `n_puddle` is 2.418e11 cm⁻² against the shipped 5e11 cm⁻². It stayed open
   eleven sessions not because the calculation is hard — it is one line, as the
   item itself always said — but because `n_puddle` is a floor on a
   **magnitude**, and a magnitude cannot say whether it floors the net density,
   the total, or each species: three readings differing by factors of 1, 1 and 2
   in `σ_min`, the same order as the answer. Item 9 settles it as the total, and
   the per-species reading moves the discrepancy the **wrong way** (4.1357×), so
   the convention is not the explanation. Cost: §4.4's on/off ratio is
   understated by 1.7357× intrinsically and 1.4362× at the terminals, the
   contacts attenuating the error by a further 1.2085× — §4.6.6's
   negative-feedback structure on a different parameter, hence a property of the
   device topology. §7.2 is undisturbed: 1.78 is no more a switch than 1.24.
   **The generalisable finding is the eleven sessions**, not the factor: a
   one-line comparison was not blocked by its difficulty but by an undeclared
   signature in one of its two inputs.
12. **Re-deriving Chapter 4 on the exact charge relation (4.28)** (added
   2026-10-08, and **the new top item of the RF thread**). The declaration in
   item 7 makes
   `V_g − V_dirac − V_ch = e·n/C_ox + E_F(n)/e` the exact relation; §4.3
   computes its linearisation. Running Eq. (4) on the exact form moves `I_d` at
   the RF bias by **−0.6023 %**, which is **1.4 times** the transport term
   items 7 and 8 were competing over. Two sessions went into a +0.43 % term
   while this sat underneath it, reachable in one line. The item is to re-derive
   §§4.3–4.6 on (4.28) and re-report, keeping the (4.29) numbers annotated in
   place. **It is blocked by item 9**: (4.28) and (4.29) both carry `n` as a
   magnitude, so re-deriving before the signed-carrier treatment exists would
   rebuild the branch-sign flaw in a new equation. Item 9 is therefore promoted
   from a sibling of item 7 to the **prerequisite** of the top item.
   **2026-10-09: UNBLOCKED AND UNTOUCHED.** Item 9 is closed and
   `n_net_exact()` supplies the branched root (odd bitwise, residual 7.690e-16).
   This item has no prerequisite left and is the **top item with nothing in
   front of it** — the first time in three sessions that is true. It is also now
   the only remaining reason the §4.6.5 rewiring item is deferred, so it blocks
   work rather than merely ranking above it.
   **2026-10-10 — CLOSED**, §4.6.8,
   `graphene_exact_charge_rederivation.py`, 21/21 checks. The −0.6023 % is
   reproduced to four decimal places from an independently written patch path,
   and **the item's own description of its two equations was wrong**. (4.29) is
   not the linearisation of (4.28); the two do not share a leading order as
   `dV → 0` (exponents 1.99999854 and 1.00000000 at `dV` = 10⁻¹⁰ V, quadratic
   against linear). (4.29) is (4.28)'s **exact differential relation used
   algebraically** — a one-point rectangle rule for the integral of
   `dn/d(dV) = (C_ox/e)·C_Q/(C_Q+C_ox)` with `C_Q(n) = 2e√n/A_F`, verified an
   identity to 1 ulp and by integrating it back to the closed-form root at
   observed order 3.994 — compounded with a **second, independent** error:
   `C_q` evaluated at the overdrive rather than at `E_F/e`. Separated for the
   first time, the argument error (+4.7416 % of the charge at the RF bias) is
   **larger** than the quadrature error (+2.6396 %), and the two are not
   additive. §§4.3–4.4 move by at most 2.0859 % with every qualitative claim
   intact. **The item was also mis-scoped**: under the differential identity
   the gate capacitance *is* the charge relation's derivative, so `C_gs` moves
   too, and in `f_T` the charge half (−1.9986 %) and the capacitance half
   (+1.9392 %) nearly cancel to −0.0594 % — two coupled 2 % corrections, not a
   negligible effect. The three sentences in this thesis that call (4.29) a
   linearisation (§4.6.6, the §4.3 annotation, and this item) are marked
   superseded in place, not deleted. The §4.6.5 rewiring item is therefore
   **unblocked**, and the argument error this exposed becomes item 16.
13. **The `cosh` overflow in `quantum_capacitance()`** (added 2026-10-08).
   `log(2(1 + cosh η))` returns `inf`/`NaN` once `η = dV/(kT/e)` passes ~710,
   i.e. above `dV` = **18.3493 V** at 300 K, measured by bisection. The exact
   large-`η` limit is `|η| + log 2`. No committed number is affected — the
   maximum drive swept anywhere in this thesis is ~2.7 V — so this is a
   *latent* defect and is ranked last deliberately. It is recorded rather than
   patched because `graphene_fet_model.py` owns committed transcripts across
   Chapters 4, 5 and 6, and a numerics change there needs its own
   output-neutrality pass of the kind 2026-10-07 ran for the f-string fix.
   **2026-10-09 amends the "no committed number affected" clause, and the
   ranking stands.** On its first use from outside the module that owns it, the
   overflow produced a wrong number: above the onset `C_q` is `inf`, so
   `C_q/(C_q + C_ox)` is `inf/inf = nan`, `nan > 0` evaluates `False`, and the
   first draft of §4.6.7's 1 %-crossover bisection scored the `NaN` branch as
   "the floor dominates", walked away from the root and returned its own 200 V
   bracket ceiling as a result. The correct value is 14.73 V. Caught before
   commit because the printed figure was 74× the swept range; nothing reached a
   transcript. **Latent and harmless are not the same property** — the defect is
   still latent with respect to committed numbers, and is now known to be a trap
   for every new caller. The ranking is unchanged for the output-neutrality
   reason above, but the *mitigation* changed: callers outside this module must
   treat non-finite `C_q` as a hard stop rather than as a magnitude.
14. **Measuring what a declaration costs, as a standing habit** (added
   2026-10-08; the mechanisable form of §4.6.6's second half). A declaration
   that resolves an ambiguity is a promise that the rest of the model belongs
   to the declared variable, and the size of that promise is computable: here
   it was one monkey-patched function and one re-run. Nothing in this thesis
   asks it of the other declarations it has made. Items 2, 3 and 6 all turn on
   a named quantity, and none of them has been costed this way.
15. **An audit of every tolerance in this repository, recording what fixed it
   at that value** (added 2026-10-09; the mechanisable form of §4.6.7's last
   subsection). Three exactness checks failed as first written today, and all
   three were measurements of something other than what they were pointed at:
   one measured the floating-point behaviour of this thesis's `V_g` sweep
   parameterisation (the ambipolar symmetry is bitwise in `dV`, and
   `(V_dirac − x) − V_dirac ≠ −((V_dirac + x) − V_dirac)` to 2.22e-16 V, which
   the `C_ox/e` lever turns into 2.0 /m²); one measured its own denominator
   (PASS and FAIL 119× apart on a 1-ulp residual, because `n_net → 0` where `S`
   does not); one measured four significant figures against 53 bits. **Three
   repairs were available at the moment of failure — loosen the tolerance, fix
   the model, or find out what it measured — and only the third produced the
   findings.** The first would have silently accepted all three. The item is to
   sweep every numerical tolerance in the repository and flag any that was set
   *after* a failure rather than derived from the precision being claimed, which
   is the only way to tell a justified tolerance from a tuned one after the
   fact. Sibling of item 6's magnitude census and of 10-06's ratio-criterion
   item, and the natural instrument for both: a tolerance set after a failure
   and a criterion written as a ratio are the two ways a check stops being a
   check without anyone editing it.
   **2026-10-10 — a near-miss of exactly this kind, caught and recorded.**
   §4.6.8's first draft required "every number moves by under 2 %" and failed
   at 2.0859 %. Widening it to 2.5 % would have produced a passing criterion
   that measured nothing but the run that set it. The bound was **deleted**
   instead: the magnitudes are printed with no threshold and the check tests
   the qualitative claims. The same session replaced a guessed 10⁻¹⁰ error
   bound on a quadrature identity with the **observed convergence order and
   the absence of a plateau**, which the method supplies and no one has to
   choose. The item now has two patterns to look for and one worked repair for
   each: *delete the bound and check the claim*, and *let the method set the
   criterion*. It is still unbuilt, and now **five sessions overdue**.

16. **Fixing `quantum_capacitance()`'s argument** (added 2026-10-10, and the
   first item on this list that arrives already costed). `C_q` is evaluated at
   the gate overdrive `dV`; the dispersion it is derived from requires
   `E_F/e = A_F·√n`, which at the RF bias is **4.04 %** of the overdrive.
   §4.6.8 separates this from the quadrature error of (4.29) and sizes it at
   **+4.7416 %** of the raw charge at that bias — **larger than the quadrature
   error (+2.6396 %), and the largest cost attached to anything currently open
   on Chapter 4**. The compact-model review
   [arXiv:2209.00388](https://arxiv.org/pdf/2209.00388) evaluates its `C_q` at
   the chemical-potential shift and never at the overdrive, so this is a
   deviation from the literature's own convention rather than a modelling
   choice this thesis is entitled to. It is ranked here rather than acted on
   for the same process reason as item 13 — `graphene_fet_model.py` owns
   committed transcripts across Chapters 4, 5 and 6, and a numerics change
   there needs its own output-neutrality pass — but unlike item 13 it is **not
   latent**: it is active in every committed number in those chapters. The
   honest statement is that it has the largest known cost on this list and is
   ranked low only for process reasons. Acting on it would be the first change
   in this thesis to **move** committed Chapter 4 numbers rather than annotate
   them, and the puddle floor's 4.3× attenuation of charge errors (§4.6.8) is
   the reason that is survivable.

---

## 7.10 What this chapter does not claim

* It does not claim the retracted Section 6.7 contradiction. §7.2.1 states why,
  and the superseded numbers are annotated in place rather than deleted.
* It does not claim that any number in this thesis is *wrong* as a consequence
  of §7.5.5. Every published sensitivity is stable to 1.3e-8 under refinement.
  The claim is about procedures and instruments, not results.
* It does not claim the design rules of §7.8 are bounds. They are rankings over
  a seven-metal table and a handful of liner scenarios, and §7.4 is an argument
  for *converting* them, not a report that they have been converted.
* It does not claim that Chapter 6's error bars are complete. They quantify the
  crossover-convention uncertainty (±17 % / ±78 %). They do not include the
  0.12 eV discrepancy, the linear-vs-nonlinear doping relation beyond its
  measured 4.9 % compression, or anything about the chemisorbed metals, for
  which no bar can currently be computed at all.
* It does not claim that graphene will or will not be adopted industrially. It
  claims that *if* it is adopted first anywhere in the devices studied here,
  the photodetector is where the physics is least in tension with itself, and
  that contact engineering gates everything.
* The `f_T ≈ 20 GHz` figure is quoted at `V_ds = 0.1 V` and `V_g = 2.0 V`. It is
  10.140 GHz at the function's signature default of `V_ds = 0.05 V`. Neither is
  the device's property; both are the model's property at a stated operating
  point, and §7.7 records what happened when that was left unstated.
