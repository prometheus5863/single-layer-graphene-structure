# 2026-10-11 — The argument was a state, the item ranked last was the blocker, and the floor and the temperature are the same effect twice

Session target: Chapter 7 §7.9 **item 16**, "fixing `quantum_capacitance()`'s
argument" — created 2026-10-10 as the first item on that list to arrive
already costed (+4.7416 % of the raw charge at the RF bias), and the first
that would *move* committed numbers rather than annotate them. The previous
session's note said the next one should do this unless something larger
appeared. Nothing larger appeared; what appeared was inside the item.

Modules: `graphene_quantum_capacitance_argument.py` (20 checks, 20 passed),
`graphene_quantum_capacitance_overflow.py` (5/5). Thesis: §4.6.9, §4.6.10.
Figure: `quantum_capacitance_argument.png`.

---

## 1. The item understated the work, not the size

`quantum_capacitance(V)` takes a potential. The quantity the dispersion
requires it be evaluated at is `E_F/e = A_F·√n`, and `n` is what `C_q` is
being used to compute. **The correct argument is a state, not a voltage.** So
"fix the argument" is not an edit to a call site; it is a self-consistent
solve. That is why the largest open cost on Chapter 4 sat for a session as a
one-line description — not because it was hard to believe, but because the
sentence describing it made it sound like a one-liner.

At the RF bias `E_F/e = 0.0625 V`, 5.21 % of the 1.2 V overdrive, and the two
candidate arguments give `C_q` a factor **17.9** apart. That is only 2.30 % in
the series factor, because both are far above `C_ox` — i.e. **the error is
small because the device is not quantum-capacitance-limited here**, not
because the argument is nearly right. The regime where the argument would
matter is the one this chapter has least to say about.

## 2. The obvious fix deletes the temperature, and yesterday's number is its cost

`C_Q(n) = 2e√n/A_F` is a T = 0 expression. `quantum_capacitance()` is
thermally rounded for a stated reason — its own docstring says *"so C_q does
not diverge to zero exactly at the neutrality point"*. Substituting the
dispersion gives `C_q(0) = 0`, hence a series factor of zero and an induced
charge of **exactly** zero at the Dirac point.

The sharp form: with the dispersion, `n ≡ 0` solves
`dn/d(dV) = (C_ox/e)·f(C_Q(n))` through the origin, so the initial-value
problem is **degenerate** there. With the thermal form it is not. *Temperature
is what makes the charge relation well-posed at neutrality.*

2026-10-10 sized item 16 with that substitution. **So its +4.7416 % was the
cost of a fix that is wrong in a second way** — which is a different failure
from the ones this repo has been collecting. The previous ones were checks
that measured the wrong thing. This was a *cost estimate* that priced the
wrong repair.

The form that keeps the temperature — thermal `C_q` at `E_F/e`, solved
self-consistently — reproduces `C_q(0) = PREFACTOR·(kT/e)·ln 4` to **0 ulp**,
reduces to the dispersion as `T → 0`, and is what §4.6.9 measures.

## 3. Every number in yesterday's decomposition was inflated

(4.28) is a T = 0 relation; the model runs at 300 K. Differencing one against
the other puts a temperature difference inside every percentage. Re-measured
in a 2×2 at **one** temperature, against the exact integral of the model's own
differential relation (an IVP, validated to −1.5e-14 against (4.28)'s closed
form):

| | vs T = 0 (4.28) | vs 300 K reference | overstated |
|---|---|---|---|
| argument error | +4.7416 % | **+2.6408 %** | 1.80× |
| quadrature error | +2.6396 % | **+0.7603 %** | 3.47× |
| both | +5.2064 % | **+3.0963 %** | 1.68× |

The **ordering survives**: the argument is still the larger error, so item 16
was ranked correctly. The magnitudes do not, and the thesis quoted them in
four places. All annotated in place.

**The trap worth remembering.** Today's 300 K *argument* error (+2.6408 %)
agrees with yesterday's T = 0 *quadrature* error (+2.6396 %) to 0.0012 %. Two
different quantities against two different references. Read in sequence, the
agreement looks like a reproduction. 2026-10-02's rule — a coincidence between
two artefacts is not a validation — applied for the first time to a *numerical*
coincidence rather than to a re-run.

## 4. The item ranked last was gating the item ranked largest

The natural exact-limit check on the corrected form is `T → 0`. It could not
be run. `η = eV/kT` grows as `T` falls, `cosh(η)` overflows at
**η = 710.4758**, and the series factor returns `inf/inf = nan`. That is item
13, added 2026-10-08, judged *latent* and ranked **last** because η = 710 needs
18.35 V of overdrive at 300 K. At T = 1 K and the ordinary RF bias, η = 719.

> **Latency had been judged against the inputs the model takes, not against
> the inputs its validations take.**

A defect unreachable from a device's operating range is still reachable by a
check *on* that device's model. Ranking by reachability put the blocker last.
This is now a standing caution on any latency claim and part of item 15's
remit. Item 13 is closed in the same session.

## 5. The neutrality pass measured the wrong thing three times

The fix is an identity — `2(1+cosh x) = 4cosh²(x/2)` gives
`ln(2(1+cosh η)) = |η| + 2·log1p(exp(−|η|))` — so the output-neutrality pass
should have been trivial. It was not:

1. **Literal constants.** `e`, `ħ`, `k_B` re-typed by hand. `scipy.constants`
   has `ħ = h/2π = 1.05457181764615647e-34`; the literal `1.054571817e-34`
   differs at the tenth digit, squared in the prefactor ~1.2e-9. The pass
   reported **3983 of 4003** values differing by up to 1.1e7 ulp. Sixth
   instance of "find out what a failed check measured", and the first where
   the substituted quantity was a *physical constant*.
2. **Association.** Constants fixed, but written `prefactor·(kT/e)·log(…)`
   against the module's `(kT/e)·log(…)` then `prefactor·…`. Floating
   multiplication is not associative; the regrouping alone moved **118 of 400**
   values by up to 2 ulp.
3. **A silent no-op.** The patch meant to fix the association matched nothing
   and reported success anyway, so draft 2 ran again and was nearly read as a
   confirmed result.

> **A function that "restores the old spelling for comparison" is a
> re-implementation, and every detail of it not copied exactly — the constants,
> the grouping, the order of operations — becomes part of the difference it
> reports.**

The repair is not a better transcription. It is to make the thing that changed
**substitutable**: `_thermal_rounding()` now holds the single expression item
13 touched; the pass swaps that one function and leaves the rest of the call
path as shipped code. That form reproduces the direct file-against-file result
exactly — **4 of 400, 1 ulp** — where all three transcriptions did not.

And the honest result is not "output-neutral": `C_q` itself moves by 1 ulp at
4 of 400 biases within 0.093 V of the Dirac point; **no committed quantity
moves by a bit** over 3603 values.

## 6. Two self-corrections inside one session

**(a) A ≤ 1 ulp criterion I chose, which failed at 2 ulp.** Measured on 11
hand-picked values it was 1; on 401 it is 2. Not widened to 2 — 2026-10-10's
fourth repair applies, and for an *identity* a ulp count cannot answer the
question that matters. `mpmath` at 60 digits adjudicates, and **the answer is
not the expected one**:

| | worst rel. error | mean | strictly better | strictly worse |
|---|---|---|---|---|
| committed `ln(2(1+cosh η))` | 1.793e-16 | 4.643e-17 | — | — |
| installed `|η| + 2·log1p(e^−|η|)` | **1.527e-16** | 5.162e-17 | 45/401 | **81/401** |

The installed spelling is **worse at more samples than it is better**. All it
has is the smaller worst case. **Accuracy is not a reason to install it**;
being finite where the other is not is the only reason — which is the only
reason item 13 ever gave. Written into §4.6.10 rather than quietly dropped,
because today's first commit message states the ≤ 1 ulp figure and a pushed
commit message cannot be amended.

**(b) An arrival control that failed on a healthy run.** One counter,
incremented only in the capacitance patch, gating a section that exercises the
*charge* patch and never asks for `C_gs`. Zero capacitance calls was the
correct answer there, and the control called it a fault. 2026-10-10's rule — a
control must not read the thing it controls for — in its other failure mode:
**a control that counts a different arrival from the one it gates reports a
fault that does not exist.** One counter per patched name now, each checked
where that name is used.

## 7. What moves, and what does not

Raw charge **−2.2659 %**; the puddle floor attenuates it by **0.451** to
−1.0212 % in `n`(3.5 V) and the crossover lengths; `I_d` at the RF bias
−0.2530 %; `R_ch` at the Dirac point **0.0000 %** because the floor supplies
the whole density there.

Every qualitative claim of §§4.3–4.4 survives — on/off still between 1 and 2,
`L_x` at 300 Ω·μm still above the 200 nm channel (314.17 → 310.96 nm),
resistance ordering unchanged. **The correction moves numbers and no verdicts.**

`f_T`: charge alone −1.0911 %, the consistent `C_gs` **+1.8976 %**, total
+0.8065 %. Opposite signs again — the **second independent instance** of the
structure 2026-10-10 recorded, on one instance, as a property of the device
*topology* rather than of any one parameter. That claim now has what it needed.
`f_max/f_T` moves +0.0009 %: §7.8.1a untouched by a **sixth** lever, and the
binding constraint is still that `transfer_characteristic()` is a resistor
with no output resistance.

## 8. The two-knob check's first use fails in the direction it cannot see

Yesterday's new item: move each input separately, flag the case where the
combined move is *smaller* than either. First instance, opposite sign — the
knobs **reinforce**, combined/(sum of separate) = **1.109**, so a check testing
for smallness would not have flagged it. It must report the **ratio** and flag
any departure from 1.

And a requirement only running it could find: **the interaction measure
inherits the reference error of whatever it is differenced against.** The same
two knobs give a departure from additivity of 0.717 against the T = 0 relation
and 0.109 against a 300 K reference — **6.6×**. A two-knob check must be run at
the model's own operating point. The item is amended, not closed, and becomes
item 17.

## 9. Literature, and a new question it opens

Live web search available and used: one search, one fetch, both successful. No
PubMed/PMC attempted, per the standing note.

A **second, independent** compact-model lineage —
[the Southampton bilayer-GFET model](https://eprints.soton.ac.uk/366757/1/GFetIEEEtran.pdf)
— confirms §1 from a different source than 2026-10-10's
[Pasadas/Jiménez review](https://arxiv.org/pdf/2209.00388): its `C_q` is a
function of the **surface potential** `φ_s`, never of the gate overdrive, and
its layer potentials form a **coupled system solved together**. That is
precisely the claim that item 16 is a self-consistency problem, now held by two
independent groups rather than inferred from one.

But it regularises `C_q` at neutrality **differently, and the difference is the
interesting part.** Temperature never enters its `C_q` at all. It floors the
capacitance with the residual disorder density instead:
`C_q,min = e²√n₀/(√π·ħ·v_F)`. The puddle density does the job this thesis gives
to `kT`.

At this thesis's `n_puddle = 5e15 m⁻²` the two are **the same size**:

| | `C_q` at the Dirac point (F m⁻²) |
|---|---|
| this thesis, thermal rounding | 8.437399e-03 |
| literature form, at `n_puddle` | 9.710793e-03 |
| ratio | **1.1509** |

They would agree exactly at `n₀ = 3.77e15 m⁻²` (0.755 × `n_puddle`) or, holding
`n₀`, at **T = 345 K** (1.151 × 300 K). So the agreement is not a coincidence
of one parameter: SiO₂-supported graphene at room temperature sits right where
disorder and temperature regularise `C_q` equally.

**And this thesis applies both.** The thermal rounding is inside `C_q`, and
`n_puddle` is then added in quadrature to `n`. The literature form shows those
are two spellings of one physical effect at this `n₀`, so **using both is a
double count of the same regularisation, worth about a factor 2 in `C_q` at the
Dirac point, and it has never been costed.** New §7.9 **item 18**. It would not
be true for suspended graphene (`n₀` two orders lower) or at 4 K, where the two
conventions diverge — which is why neither is obviously wrong in the
literature and why the double count here is easy to miss.

This also puts a number on how little §4.6.9's measurement can say about
neutrality: the region where the argument error would be large is the region
where two overlapping regularisations are fighting over the same physics.

## 10. Environment

Device reachable and folder connected at the first firing of the day. `pip
install scipy` succeeded again (the device VM still lacks it); `mpmath` was
already present and is used as the 60-digit adjudicator. The whole session ran
on the device VM; nothing crossed the bridge. Each commit was pushed and
verified against the GitHub API before the next was written.

## 11. Not yet covered

- **COST THE DOUBLE REGULARISATION** — item 18, created today, and the first
  item to arrive with a literature source rather than from an internal audit.
  Ratio 1.1509 at the operating point, and a factor ~2 in `C_q(0)`.
- **FLIP THE DEFAULT TO `B_T`** — item 16 stays open. Now costed on the right
  form, with the blocker gone. It needs a session of its own for transcript
  regeneration across Chapters 4, 5 and 6.
- **BUILD ITEM 17 AS A SWEEP** — one worked instance and two requirements
  (report the ratio; run at the model's operating point), still no sweep.
- **AUDIT EVERY TOLERANCE** — item 15, now **six sessions overdue**, and
  today supplied a sixth worked instance: a ≤ 1 ulp bound I chose, which
  failed, and whose repair was to let 60-digit arithmetic set the criterion.
- **A THIRD KNOB: THE ROUNDING FORM.** §4.6.9's 2×2 holds the rounding fixed.
  Fermi–Dirac integrals have never been compared against the phenomenological
  `ln(2(1+cosh η))` at all, and item 18 is the same question from the other
  side.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — untouched for the third session;
  still the cheapest open physics item.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, nothing in front
  of it, untouched for a second session.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  **untouched for the fifth session**. Today produced a new instance after all:
  `f_max/f_T`, used again as the structural verdict's instrument.
- **MECHANISE the magnitude item** — created 2026-10-03, **six sessions
  overdue**; items 15 and 17 subsume parts of it.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, untouched for the third
  session.
- **WHETHER `n_puddle` SHOULD ACTUALLY CHANGE** — deliberately still not an
  item, and item 18 is now the reason to look at it.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** — still
  expressible and unverified, untouched for the third session.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, **untouched for
  TWENTY-ONE consecutive sessions**. Today was again a numerics-and-method
  session. §7.7 gained a sixth instance; this is the standing one.
