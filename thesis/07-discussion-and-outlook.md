# Chapter 7: Discussion and Outlook — Where Single-Layer Graphene Realistically Fits, and What This Thesis Actually Established


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

## 7.3a The four strands are two process steps, and only one of the three couplings is a trade-off

### 7.3a.1 Why this is one problem and not four

The contact/photoresponse item has accumulated four strands since 2026-09-17,
each logged as it arrived and none of them yet synthesised. They are:

1. **Chapter 4 §4.7** — contact resistance is lowest for the high-work-function
   metals Ni, Au and Pd, and `R_c` is computed for one interface at a time.
2. **Chapter 6 §6.14** — the two-terminal zero-bias photoresponse `N` is
   *exactly* odd in the differential work-function offset `s` and *exactly*
   even in the common offset `m` (measured to 2.3×10⁻¹⁵), with a reversal
   threshold `τ* = W_B − W_A` independent of λ, of `L` and of every transport
   parameter.
3. **§4.11 and the 2026-10-30 result** — the same contact-induced
   work-function mismatch that pins the contact doping is the mechanism behind
   the contact-origin branch asymmetry `r₀`, and whether `r₀` survives at the
   operating density is a screening question governed by `x ≡ |n|/n_c`.
4. **The 2026-11-02 and 2026-11-03 results** — `n_c` is not a material
   constant but a charged-impurity density, `n_c = κ·N_imp`, with
   `N_imp ≥ N_A^(S) = C_ox·V_NP/e`; the dielectric that provides a strong gate
   also raises `N_imp` through the damage its deposition leaves behind.

Listed that way they look like four separate facts about graphene. They are
not. **Each is a consequence of one of exactly two process steps**, and the
synthesis this chapter owes is the statement of which step sets what, and of
which couplings between them are real constraints rather than apparent ones.

### 7.3a.2 The two steps, and the quantities each one sets

**Step A — contact metallisation.** Choosing the two contact metals fixes
`W_A` and `W_B`, and therefore both coordinates of the pair:

&nbsp;&nbsp;&nbsp;&nbsp;`m = ½(ΔW_A + ΔW_B)` — the common offset, which is what a
single-junction contact-resistance figure of merit responds to (strand 1);

&nbsp;&nbsp;&nbsp;&nbsp;`s = ½(ΔW_A − ΔW_B)` — the differential offset, which sets
`sign(N)` and `τ*` (strand 2);

and, through the contact-induced doping profile, the contact-origin component
of `r₀` (strand 3's source term).

**Step B — gate dielectric deposition.** Choosing thickness and permittivity
fixes `C_ox`, and the deposition leaves behind a charged-impurity density
whose *net* part is what the neutrality-point offset `V_NP` measures. Those
two together fix

&nbsp;&nbsp;&nbsp;&nbsp;`|n| = C_ox·ΔV/e` &nbsp;&nbsp;and&nbsp;&nbsp;
`n_c = κ·N_imp ≥ κ·C_ox·V_NP/e`,

and hence `x = |n|/n_c`, which is the only quantity that decides how much of
`r₀` survives (strands 3 and 4).

So there are two steps, four quantities of interest, and three couplings worth
asking about. Two of the three turn out not to be trade-offs at all, and the
third is a trade-off that cannot be resolved with the knob one would reach for.

### 7.3a.3 Coupling 1 — contact resistance against photoresponse is not a trade-off, and the reason is a parity

This is §7.2's result, read as a process statement. `m` and `s` are
*independent coordinates* on the metallisation choice, and §7.2.2 establishes
that `N` is exactly even in `m` and exactly odd in `s`. A process is therefore
free to optimise `R_c` through `m` and the photoresponse through `s` at the
same time: the two criteria are not competing for the same degree of freedom.

**The 2026-09-17 contradiction was never a trade-off.** It was a
single-junction criterion being applied to a quantity that depends on the
difference between two junctions, and §7.2.1 has already retracted the metal
list it produced. What the four strands add here is that the retraction is not
merely a correction of a ranking: it removes the apparent conflict between
strands 1 and 2 entirely.

The cost is real but it is a cost in *process complexity*, not in device
physics: `s ≠ 0` requires two metallisation steps where a symmetric device
needs one.

### 7.3a.4 Coupling 2 — the decoupling is bought with model risk, and the risk sits exactly where the model has no reach

§7.2.4 records the anti-correlation; as a process statement it is sharper.
Maximising `|s|` means choosing metals far apart in work function, and the
pairs with the smallest reversal margins — Au/Pd at 0.020 eV, Ni/Au at 0.060
eV, Ni/Pd at 0.080 eV, Cr/Cu at 0.150 eV — are **exactly the four pairs the
per-metal crossover model of §6.12 refuses**, because Ti, Ni and Pd are
chemisorbed and sit ~1.2 Å below the physisorbed anchor.

So Coupling 1's freedom is not free: it is exercised in the region where this
thesis's own model cannot predict the quantity that decides polarity. That is
a statement about the state of the modelling, not about the device, and
§7.2.4's wording — *the model's reach and the thesis's risk are
anti-correlated by construction* — is the correct summary of it.

### 7.3a.5 Coupling 3 — Step B cannot be used to fix Step A, because the obvious knob cancels exactly

This is the coupling the four strands were accumulated for, and it is the one
place a genuine constraint survives. It also **reverses the conclusion**
strand 3 was logged with, so the reversal is stated before the synthesis.

The 2026-10-30 result turned `n_c` into a crossover gate overdrive
`ΔV_c = e·n_c/C_ox` and evaluated it for two stacks:

| gate stack | `C_ox` (F m⁻²) | `ΔV_c`, one fixed `n_c` |
|---|---|---|
| 90 nm SiO₂ back gate | 3.84×10⁻⁴ | 12.53 V |
| 10 nm HfO₂ top gate | 2.21×10⁻² | 0.217 V |

and concluded that whether the asymmetry matters *"is a property of the gate
stack, not of the graphene"*, being screened away in the thin-oxide stacks a
semiconductor process would actually build. **That table uses one value of
`n_c` for both rows, and strand 4 says `n_c` is not shared between them.**

Substituting `N_imp = f·N_A^(S) = f·C_ox·V_NP/e` with `f ≥ 1`, the `C_ox` in
`n_c` cancels the `C_ox` in the denominator:

&nbsp;&nbsp;&nbsp;&nbsp;**`ΔV_c = κ·f·V_NP`** &nbsp;&nbsp;and&nbsp;&nbsp;
**`x = ΔV / (κ·f·V_NP)`**

The crossover overdrive is independent of oxide thickness and of oxide
permittivity alike. Under this reading the two stacks need the *same*
overdrive to reach `n_c`, and `graphene_screening_overdrive_invariance.py`
measures that difference as exactly `0.0` V rather than as small; the 57.69×
spread in the table above is exactly the inverse `C_ox` ratio, which is what
identifies it as a statement about the oxide rather than about the graphene.

**The device consequence is a reversal, not a rescaling.** With the bounded
`κ ≤ 2.0889` of §4.11.5a and Zion's `V_NP = 20 V`, `κ·V_NP = 41.78 V` — more
overdrive than *either* stack sees in operation. At 1 V of overdrive on the
10 nm HfO₂ top gate the fixed-`n_c` reading gives `r = 0.0180`, a 5.6×
screening; the scaled reading gives `r = 0.0977`, which is **97.7 % of the
unscreened `r₀ = 0.10`**, a factor of 192 in `x`. The branch asymmetry is
essentially *unscreened* in precisely the stack 10-30 said would screen it.

The superseded numbers are left in place in their own section and annotated
there, per the standing practice of this work; what is withdrawn is the
conclusion, not the arithmetic, which reproduces exactly (12.527 V and
0.2171 V) when the fixed-`n_c` premise is granted.

**And the same shape as §7.2.** It is worth naming, because this is now the
second time it has decided a device question in this thesis: in both cases the
knob a process engineer would reach for enters the quantity of interest in a
way that cancels. `m` enters `N` evenly, so contact quality cannot move the
photoresponse at leading order. `C_ox` enters `x` through both the numerator
and the denominator, so oxide scaling cannot move the screening at all. Two of
this thesis's three sharpest device statements are of the form *the obvious
lever is not a lever*, and both were found by asking what a quantity is
**exactly** even or invariant in, rather than by surveying values.

### 7.3a.6 The process statement, in the form a process engineer can act on

| quantity | set by | the lever that works | the lever that does **not** |
|---|---|---|---|
| contact resistance `R_c` | Step A | the common offset `m` — high-work-function metals | `s`: the pairing is irrelevant to one junction |
| `sign(N)`, margin `τ*` | Step A | the differential offset `s` — deliberately mismatch | `R_c` / `m`: exactly even, so no first-order effect |
| contact-origin `r₀` | Step A | contact chemistry and geometry | the gate stack: `r₀` is a source term, not a screened one |
| surviving asymmetry `r` | Step B | `V_NP` (how clean the deposition is) and the overdrive ratio `ΔV/V_NP` | `C_ox`: thickness and permittivity cancel exactly |

The last row is the new one and it is the practically expensive one. A process
that wants the electron/hole asymmetry screened away cannot buy that by
thinning the dielectric or by moving to a higher-permittivity one — the two
changes that a scaled device makes for every other reason. It can buy it only
by leaving less net charge behind (smaller `V_NP`) or by gating harder relative
to that charge. Those are deposition-quality and operating-point decisions, not
stack-geometry decisions.

### 7.3a.7 What this section does and does not claim

**It does not claim that `C_ox` is irrelevant to the device.** `C_ox` sets the
carrier density at a given gate voltage and therefore sets transconductance,
`f_T` and everything Chapter 4 computes from them. The claim is narrower and
exact: it does not appear in the *ratio* `x` that governs screening, because it
appears identically in `|n|` and in the lower bound on `n_c`.

**The cancellation rests on one identification, and it is named.** It holds
under `N_imp = f·C_ox·V_NP/e` with `f` not systematically different between
stacks — that is, under the assumption that the neutrality-point offset is a
faithful proxy for the charged-impurity density up to a stack-independent
factor. A dielectric whose deposition damage were uncorrelated with the net
charge it leaves would break it. That assumption is a literature question and
is logged as one; it is not established here.

**The direction of the inequality is the safe one.** `N_A^(S)` is a net charge
and therefore a lower bound on the impurity count, so `f ≥ 1`, so `x` is
bounded *above* and the surviving `r` is bounded *below*. The bound cannot
promise a device more screening than it gets. §7.3a.5's numbers are therefore
conservative in the direction that matters.

**`κ` is bounded, not known.** §4.11.5a excludes `κ = 2.4` and bounds
`κ ≤ 2.0889` from Zion's own two numbers. The 41.78 V inherits that bound and
would move with it. The *structural* result does not: the module checks that
`ΔV_c/(κ·V_NP) = 1` to two ulp over four decades of `κ`, so the cancellation —
and with it the withdrawal of the 57.69× spread — does not depend on resolving
the `κ` dispute at all.

**It does not resolve §4.11.5a's open question**, and that question cuts the
same way. If the neutral-defect channel is branch-asymmetric rather than
branch-symmetric, then at `s = 2r₀` the asymmetry does not screen at *any*
density, exactly. That would strengthen this section's conclusion rather than
weaken it, which is worth stating so that the open item is not read as a
caveat pointing the other way.

**It says nothing about Chapter 5.** §7.3's differential-interconnect
prediction remains untested and is unaffected by any of this.

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

## 7.5–7.7 How this thesis's claims failed, and how they consolidated (summary)

Five failure classes are on record, each with its own detector: (i) an implementation error inside a model, caught by an exactly known value; (ii) a claim quantified on an unrepresentative sample, caught only by full enumeration; (iii) an error in the analysis layer around a validated model; (iv) a defect created by a fix, caught by keeping the superseded path as an oracle; and (v) a procedure that certifies itself, caught by a derived bound reported beside the measurement. The detectors do not overlap, so verification effort has to be allocated by failure class, not uniformly; claims that survived did so by consolidating under independent checks, and scrutiny across the thesis has been unequal. The taxonomy, the ledger of consolidated claims and the scrutiny comparison are in [Appendix B](appendix-B-failure-taxonomy.md#s7-5).

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

Six dated annotations (2026-10-05 to 2026-10-10) record how each RF lever of Chapter 4 — saturation, the perfect contact, diffusion, the named potential, the signed carrier, the exact charge relation — was built and what it did to this verdict; none reverses it. They are in [Appendix A, §§7.8.1b–g](appendix-A-verification-record.md#s7-8-1b).

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

Stated as falsifiable items, in order of how much they would move. Each item's full history — what was tried, what it changed and how it was checked — is in [Appendix A, §7.9](appendix-A-verification-record.md#s7-9).

1. A second anchor for `Δ_c` at any separation other than 3.3 Å
2. A description of Ti, Ni and Pd that does not go through work function
3. Reconciling the 0.12 eV measured step with the 0.25–1.07 eV assumed offsets
4. A photo-thermoelectric term
5. Asking the differential question of Chapter 5
6. Re-running Chapter 4's contact-resistance results at the 5.4 eV crossover
7. Naming the potential in Eq. (4) — closed.
8. A self-heating term — closed.
9. A signed-carrier FET model — closed.
10. A census of every potential-like variable in this thesis — closed.
11. Whether the `n_puddle` regularisation removes other regimes, not only the diffusion-dominated one — closed.
12. Re-deriving Chapter 4 on the exact charge relation (4.28) — closed.
13. The `cosh` overflow in `quantum_capacitance()` — closed.
14. Measuring what a declaration costs, as a standing habit
15. An audit of every tolerance in this repository, recording what fixed it at that value
16. Fixing `quantum_capacitance()`'s argument
17. A two-knob check that reports the ratio, not smallness
18. Costing the double regularisation of `C_q` at neutrality — closed.
19. Making the rebind set a computed object rather than a remembered list
20. Converting the 20 inert knobs, and the four dead names behind them
21. The standing list is not a closed system, and two of its items asserted the opposite of what the tree contained
22. Which `v_F` this thesis is quoting, and the fact that it cannot currently say
23. The configurations §4.6.22 did not measure, and a verdict indexed by configuration
24. Does `n_c` differ between the electron and hole branches, and does the asymmetry therefore screen at all?
25. The `A` pairing, closed 2026-11-03, and what replaced it

## 7.10 What this chapter does not claim

* It does not claim the retracted Section 6.7 contradiction. §7.2.1 states why,
  and the superseded numbers are annotated in place rather than deleted.
* It does not claim that the error decomposition of §4.6.9(c) is the last word
  on (4.29). It separates two knobs at the model's own temperature against the
  exact integral of the model's own differential relation. The relation itself
  still uses a phenomenological thermal rounding rather than Fermi–Dirac
  integrals, so a third knob — the rounding form — has never been costed at
  all, and (4.28) remains the T = 0 limit of a standard literature form.
* It does not claim that the overflow-free spelling installed in §4.6.10 is
  more accurate than the one it replaces. Adjudicated at 60 digits it is
  strictly worse at 81 of 401 samples. The claim is that it is finite
  everywhere and that no committed number moves by a bit.
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
