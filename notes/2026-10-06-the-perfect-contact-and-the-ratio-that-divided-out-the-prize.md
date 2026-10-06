# The perfect contact, and the ratio that divided out the prize

**Date:** 2026-10-06
**Modules:** `graphene_perfect_contact_counterfactual.py`,
`graphene_perfect_contact_mutation.py`
**Transcripts:** `perfect_contact_counterfactual_output.txt`,
`perfect_contact_mutation_output.txt`
**Figure:** `perfect_contact_counterfactual.png`
**Closes:** the top open item of 2026-10-05 (f_max/f_T at `R_c = 0` exactly, in
the saturated model)

---

## 1. The item, and why its dichotomy was not exhaustive

2026-10-05 created this item and wrote it as a clean two-branch decision:

> If a perfect contact also falls short of 1.3, four levers are exhausted and
> §7.8.1a's structural verdict is final; if it reaches the band, the RF case
> becomes a contact-engineering problem with a quantified target. **Either
> outcome is more useful than what this thesis currently carries.**

Both branches turned out to be true at once. At the literature-scale
40 µm / 8-finger geometry and `V_ds = 0.1 V`, with the contact removed exactly:

| quantity | baseline | `R_c = 0` | change |
|---|---|---|---|
| `f_max/f_T` | 0.683186 | 0.759363 | **+11.15 %**, still 41.6 % below 1.3 |
| `f_T` | 20.740 GHz | 74.650 GHz | **×3.5994** |
| `f_max` | 14.169 GHz | 56.687 GHz | **×4.0007** |

So the *ratio* falls short — the fourth lever is exhausted, §7.8.1a's verdict
on `f_max/f_T` is final — **and** the RF case becomes a contact-engineering
problem with a quantified target anyway, because the same counterfactual is
worth a factor of four in `f_max` itself.

The item could not see this, and the reason is mechanical rather than careless.
`f_max/f_T` divides out exactly the quantity the contact dominates. The contact
limits `I_d`, `I_d` sets `g_m`, `g_m` sets `f_T`; the contact also sets `R_s`,
which sits in the `f_max` denominator. Removing it moves numerator and
denominator together, and the ratio keeps only the residue. **A criterion
written on a ratio is silent about any mechanism that scales both of its
arguments.** That is this session's methodological finding and §6 states it as
the next item in the series.

## 2. The second result: the counterfactual's sign depends on its coherence

`R_c` enters this model in two places, and they are different code paths
reading different names:

1. as a series resistance inside `transfer_characteristic_saturated`'s
   self-consistent solve for the intrinsic channel drop, and
2. as `R_s = R_c/2` in the `f_max` denominator, via
   `rf_small_signal_model.source_access_resistance()`, which reads the
   **module global** `gfet.Rc_total`.

Passing `Rc_total=0.0` through `**kw` — the obvious way to write the
counterfactual, and the way `fT_fmax_saturated` would propagate it — zeroes
only (1). The result:

| counterfactual | `f_max/f_T` | change |
|---|---|---|
| baseline | 0.683186 | — |
| `R_c = 0` in the current path only | 0.552349 | **−19.15 %** |
| `R_c = 0` coherently (`R_s = 0` too) | 0.759363 | **+11.15 %** |

**The half-removed perfect contact answers the item with the opposite sign.**
Worse, the half-removed answer is the *plausible* one: a 19 % degradation is
what a careless reading of "the contacts are binding" predicts, and nothing
about the number looks wrong. Check C3 asserts the sign inversion so the trap
is in the transcript rather than in a comment.

This is a new location for a fault this repository has now hit repeatedly: a
quantity that enters a calculation through two names, where changing one is
indistinguishable, locally, from changing the thing. 2026-10-03 found it as a
frozen default in another module; today it is a counterfactual that arrives in
one of two places.

## 3. The third result: the two models disagree on the SIGN of the contact lever

The same coherent counterfactual, in the superseded resistor model:
**0.642901 → 0.581707, i.e. −9.52 %**, against the saturated model's
**+11.15 %**. Not a magnitude discrepancy — a sign disagreement, asserted by
check C4 and left unreconciled.

It has a diagnostic reading, and it is the strongest single sentence available
from this session. 2026-10-04 established that in the resistor model

```
  f_max/f_T  ≈  (1/2) sqrt( R_total / (R_g + R_s) )
```

to 0.36 %, and concluded that this is "a different quantity wearing f_max's
name". That identity carries `R_total` in the **numerator**, so any model
obeying it is *rewarded* for a worse contact. 2026-10-04 read that off the
algebra. C4 is its falsifiable consequence, and it holds. **A figure of merit
that improves when the device gets worse is not measuring the device.**

Section 6b then found that the saturated model is not exempt, only
better-behaved over a range (check Q3): `f_max/f_T` is **non-monotonic** in
`R_c`, bottoming out at 0.680047 at 470 Ω·µm per contact and rising again to
0.698615 at 4000 Ω·µm, while `f_max` itself falls 37× across those same two
rows. The ratio and the figure of merit it is supposed to summarise move in
**opposite directions** over part of the design space.

## 4. What 2026-10-05 predicted, and what was measured

2026-10-05's explanation of its own shortfall made an untested quantitative
prediction. It measured a `g_ds` saturation factor of 1.1361 against a required
4.8483 and attributed the gap to dilution: `µS/L = 0.3036` at that bias "alone
would divide the CHANNEL conductance by 1.699", and the contacts, being in
series, absorb the rest. If that is the mechanism, then **at `R_c = 0` the
factor must be 1.699.**

Measured today, from a finite difference of a self-consistent solve rather than
from a channel integral: **1.70661** at the item's own bias (100.45 % of the
prediction) and 1.69621 at the other reference bias (99.84 %). Check D1 asserts
it inside 2 %.

This is a stronger form of evidence than the agreement-between-two-artefacts
that 2026-10-02 warned about, because the two numbers come from different
computations of different objects — one an integral of `1/v_sat` along the
channel, the other a derivative of the solved terminal current — and the
prediction was recorded before the measurement existed.

## 5. The literature, and it complicates this thesis's framing rather than confirming it

### 5.1 The paper that asks today's question in its title

**P. C. Feijoo, F. Pasadas, M. Bonmann *et al.*, "Does carrier velocity
saturation help to enhance f_max in graphene field-effect transistors?",
*Nanoscale Advances* **2** (2020), DOI
[10.1039/c9na00733d](https://pubs.rsc.org/en/content/articlepdf/2020/na/c9na00733d).**

This is the most consequential literature finding of the session, because it
qualifies the premise of 2026-10-04's and 2026-10-05's entire thread. That
thread read: the model has no current saturation, therefore `g_ds` is the
channel conductance, therefore add velocity saturation and `f_max` should move.
Feijoo *et al.* 2020 report the opposite of the implied conclusion:

- "the largest `f_max` are located at biases close to the onset of bipolar
  conduction and **far from the saturated velocity regime**";
- at the bias where `f_max` peaks the drift contribution is only ~45 % of the
  saturated value and the **diffusion** contribution to the current is
  comparable to the drift contribution;
- "our results **do not support** that operating in the regime of velocity
  saturation results in the highest `f_max`".

Their device is `L_g = 500 nm` CVD graphene on 22 nm Al₂O₃ over 1 µm SiO₂ on
high-resistivity Si, with extrinsic `f_T,x = 34 GHz` and `f_max = 37 GHz`, and
they use `v_sat = 4 × 10⁷ cm/s`. They also report that **self-heating alone
degrades `f_max` from 65 to 40 GHz**.

Three consequences for this thesis, stated rather than absorbed:

1. **2026-10-05's criterion A may have been the wrong criterion**, not merely a
   missed one. It demanded that a saturation term divide `g_ds` by 4.8483; this
   paper says the devices that reach the 1.3–1.4 band are not operating in the
   saturated regime at all. The criterion was constructed from *this model's*
   algebra, and the agreement of its target with the literature band was never
   checked for mechanism — only for number. The criterion stands as a
   statement about this model; it does not stand as a statement about what real
   GFETs do.
2. **This model has no diffusion current.** Eq. (4) is a drift-only
   drift-diffusion form with the diffusion term dropped. If diffusion is ~40 %
   of the current near the drain at peak `f_max` in a real device, then a
   drift-only model cannot be asked about the peak-`f_max` bias at all, and the
   `V_ds` ladders of 2026-10-05 and today are reporting a drift-only slice of a
   two-mechanism problem. This is a new top-tier open item.
3. **This model has no self-heating.** A mechanism worth 65 → 40 GHz in `f_max`
   is absent, which bounds how much any of today's `f_max` numbers can mean in
   absolute terms. Today's ×4.0007 is a *ratio* between two runs of the same
   model and is less exposed to this than the absolute 56.687 GHz is.

### 5.2 Contact resistance: what is actually reachable

`RC_PER_WIDTH_LADDER` in the module is anchored to these:

- **165 Ω·µm** — Feijoo *et al.* 2020 (above) report `R_c W_g/2 = 165 Ω·µm`,
  with the metal/graphene contact resistivity alone ≈ 90 Ω·mm. This is the
  device family whose `f_max/f_T = 1.3–1.4` band Chapter 4 is measured against,
  which makes it the load-bearing row of Section 6b.
- **470 Ω·µm** — B. Khosravi Rad, A. H. Mehrfar, Z. Sadeghi Neisiani, M. Khaje
  and A. Eslami Majd, "Effect of fabrication process on contact resistance and
  channel in graphene field effect transistors", *Scientific Reports* **14**,
  9190 (2024), DOI
  [10.1038/s41598-024-58360-9](https://www.nature.com/articles/s41598-024-58360-9).
  Their "two-in-one" Ni process (metal deposited on graphene *before*
  photolithography, so photoresist residue never reaches the interface) gives
  470 Ω·µm, against 4 kΩ·µm for their hybrid method and 20.5 kΩ·µm for
  graphene-on-metal. They state that most of the literature exceeds 4 kΩ·µm.
- **65 Ω·µm** — the same paper's citation of Liu *et al.* (2019), electron-beam
  lithography with a bottom-contact configuration, as the lowest reported.
  **Recorded as a secondary citation: the primary was not fetched this
  session**, so the 65 Ω·µm row of Section 6b rests on Khosravi Rad *et al.*'s
  report of it and should not be quoted as independently verified.
- **300 Ω·µm** — `gfet.Rc_per_width_ohm_um`, this model's own value, **1.8× the
  165 Ω·µm of the devices it is compared against.** That ratio was not known
  until today and it is why Q2 exists.

### 5.3 Fetch failures, recorded

- `pubs.aip.org` returned **HTTP 403** for "Impact of contact resistance on the
  performances of graphene field-effect transistor through analytical study"
  (*AIP Advances* 11, 045220). Not retried. This is the fourth distinct
  publisher to 403 this repository (after ResearchGate rate-limiting, Wiley,
  and the standing PubMed/PMC reCAPTCHA note), and the pattern is now worth
  stating: **every literature number in this thesis that came from a 403'd
  source came from a secondary citation**, which is the 2026-10-01 fault (a
  citation is a claim) in its most ordinary form.
- PubMed/PMC were **not attempted**, per the standing note (three logged
  reCAPTCHA occurrences).
- `arxiv.org/pdf/1110.0978` (Rodriguez, Vaziri, Östling, Rusu, Alarcón, Lemme,
  "RF Performance Projections of Graphene FETs vs. Silicon MOSFETs") fetched
  successfully but is **not usable for this question**: it states explicitly
  that it includes "not a model for contact resistance, but only a contact
  resistance parameter", and it does not decompose `f_max`. Recorded so the
  next session does not re-fetch it for the same purpose.

## 6. The methodological item, continuing the series

- **09-28:** prose is a detector.
- **09-29:** a mutation that does not arrive is indistinguishable from a system
  that does not respond.
- **09-30:** a control has to sit where the failure enters.
- **10-01:** and name what it compares against in a way that cannot drift.
- **10-02:** when it agrees, that is a fact about two artefacts, not the world.
- **10-03:** a PASS/FAIL at zero is silent about magnitude.
- **10-04:** and a quantity can be computed correctly under the wrong name.
- **10-05:** and a criterion with a number in it can be missed in a way that
  locates the real constraint.
- **10-06: AND A CRITERION CAN BE DECIDED CORRECTLY ON THE QUANTITY IT NAMES
  AND STILL MISS THE FINDING, IF THAT QUANTITY IS A RATIO AND THE MECHANISM
  SCALES BOTH ITS ARGUMENTS.**

The 10-05 item was well-formed by every standard this series has accumulated:
it named a quantity, attached a number, and stated both branches of the
decision in advance. It was answered on exactly the quantity it named. And the
×4.0007 in `f_max` — the only part of the result an industrial reader would act
on — appears nowhere in either branch, because `f_max/f_T` is invariant under
the mechanism to within 11 %.

**The practical rule, and it is checkable:** for every criterion written on a
RATIO, ask what the mechanism does to the numerator and the denominator
separately, and if it does the same thing to both, the criterion is the wrong
instrument whatever number it carries. Q3 turns this from advice into a
measurement: over 470–4000 Ω·µm the ratio and `f_max` move in opposite
directions, so there is a region of this very design space where satisfying a
ratio criterion means *losing* on the figure of merit.

## 7. The mutation harness, and the three holes it found

Final: **7 of 7 killed, control A green, control B survived.** First run: 3 of
7. The three survivors across the runs were one fault in three places — a
quantity or a claim that nothing in the battery actually read.

- **M5** (one-sided `g_m` stencil) survived because `f_max/f_T` is nearly
  *blind* to `f_T` here: term B is 1.11 % of the `f_max` denominator, so a 1e-5
  relative shift in `f_T` moves the ratio by ~5e-8. Note that this is the same
  fact as §1's: the ratio divides out what matters. Pinning `f_T` against the
  committed transcript would **not** have repaired it, because the committed
  `f_T` is an `np.gradient` over a 0.015 V gate grid while this module
  differences at `dVg = 1e-4` — a different *estimator*, differing by 1.4e-5 on
  the unmutated module. That gap is the `O(h²)` discretisation of the committed
  `f_T` and had never been quantified. Repair: **X2b**, a Richardson ratio on
  the stencil order (measured 3.9999, against 4 for central and 2 for
  one-sided), plus **X2c**, which reports the estimator gap instead of
  asserting it away.
- **M6** (term B dropped from the residual-requirement formula) survived
  because R1 asserted only that the residual was finite and > 1. Repair:
  **R1b**, the same factor recomputed by bisection on `f_max/f_T`, sharing no
  algebra with the closed form.
- **M3** (`RC_ZERO` set to 1e-12) survived two runs because X1a was written
  against a **literal** `0.0` while the counterfactual used `RC_ZERO` — an
  exactness check pointed at a value the module under test did not use. The
  2026-09-30 arrival fault, inside an exactness check, which is a new location
  for it.

Battery 17/17 → 20/20 → 23/23 (Section 6b added Q1–Q3). This is the fifth
consecutive harness in this repository to improve the suite it was pointed at
rather than bless it.

## 8. What this note does not claim

- It does not claim the saturated model's absolute `f_max` numbers mean
  anything in the world. §5.1 names two mechanisms the model lacks (diffusion
  current, self-heating), one of which is worth 65 → 40 GHz in a real device.
- It does not claim the `V_ds = 0.5 V` row of §7 of the transcript, where the
  perfect contact enters the Feijoo band, is agreement with Feijoo *et al.*
  2026-10-04 and 2026-10-05 both recorded why: a ×20 bias knob against a ×20
  literature band will coincide somewhere.
- It does not retract 2026-10-05's criterion A. §5.1 states that it is a
  statement about this model rather than about real devices, which narrows what
  it licenses without changing its arithmetic. The numbers 4.8483 and 1.1361
  stand as computed.
- The 65 Ω·µm row rests on a secondary citation (§5.2).
