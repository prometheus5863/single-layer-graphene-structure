# 2026-10-05 — Saturation arrived, and the contacts ate it

*Opens the 2026-10-04 top item: a saturation term in `transfer_characteristic()`.
Module: `graphene_velocity_saturation_model.py` (26/26). Harness:
`graphene_velocity_saturation_mutation.py` (7 of 7, control A and control B).*

---

## 0. The item, and the fact that it came with numbers

The 2026-10-04 session closed the f_max decomposition with a diagnosis rather
than a direction: `transfer_characteristic()` is a bias-dependent resistor, so
its g_ds **is** its channel conductance, term A of the f_max denominator
collapses to a resistance ratio, and f_max is a measurement of output
resistance in a model that has none. It then created the first open item in
this series to arrive with a **measured acceptance criterion** rather than a
topic:

> **(A)** divide g_ds by ≈ 4.85 at the peak-f_T bias; **(B)** reverse the sign
> of d(f_max/f_T)/dV_ds, which currently runs the wrong way; **(C)** validate
> against the resistor model in the low-V_ds overlapping limit.

That is a better item than this repository usually writes, and it is worth
saying why before reporting against it: a criterion with a number in it can be
**missed**, and a missed criterion is information. A criterion phrased as "add
velocity saturation" could only have been *done*.

**Result: B met, C met, A missed by 4.3×** — and the reason A is missed is not
the saturation law.

---

## 1. The model, and why it stayed closed-form

Soft-saturation (Caughey–Thomas) local drift velocity at exponent β = 1, which
is the γ = 1 form Feijoo *et al.* (2019) use:

```
  v(E) = μE / (1 + μE/v_sat)
```

Current continuity fixes I_d along the channel, so the velocity the channel is
*required* to deliver at potential V is u(V) = I_d/(Weₙ(V)), and inverting for
the field that delivers it gives E = u/(μ(1 − u/v_sat)). Since E = dV/dx, the
channel-length integral separates:

```
  I_d = μWe·Q / (L + μ·S),    Q = ∫ n dV,    S = ∫ dV/v_sat
```

closed on the contacts by V_ds,ch = V_ds − I_d·R_c. **The separation survives a
position-dependent v_sat**, which matters because v_sat here is
density-dependent and the density varies along the channel. Setting S = 0
recovers drift-diffusion with no saturation, which is the overlapping limit
criterion (C) asks for — it is the *same code path*, not a second
implementation, so the comparison cannot drift.

Eq. (3) has a real ceiling: no solution exists once u ≥ v_sat. The model
**reports** it (max u/v_sat = 0.721 at V_ds = 1 V) rather than clipping it.
Section 5 below is about how nearly that distinction was lost.

### Why this is a new module and not an edit

`transfer_characteristic()` is read by five modules and quoted in committed
transcripts and in Chapter 4. Editing it in place would invalidate all of them
in one commit, and 2026-10-02's rule — when a rule changes, every number
derived from it is stale until re-derived — makes that a separate deliberate
act. The saturated model is added *alongside*, the resistor model is checked
bitwise unchanged (G4), and what the rewiring would cost is now measured
rather than guessed. Same discipline as the verification repository's
2026-10-04 session: a subclass, not an edit.

---

## 2. v_sat is not a free parameter, and that is the whole methodological point

Criterion A names a number. v_sat is the one knob that moves it. **Choosing
v_sat to make the criterion pass would be fitting the model to its own
acceptance test** — 2026-10-02's fault ("when it agrees, that is a fact about
two artefacts, not about the world") in its purest available form. So v_sat
comes from optical-phonon emission with no adjustable scale:

```
  v_sat(n) = (2/π)·Ω/√(πn),   ħΩ = 0.10 eV  (Feijoo et al. 2019, Table 1)
```

ħΩ is the only literature number entering, taken as published. At this
device's density Eq. (6) gives **6.63 × 10⁷ cm/s** — *above* Dorgan *et al.*'s
measured 1–3 × 10⁷ cm/s band on SiO₂, i.e. the **generous** end of the
physics, which is the direction that makes the shortfall below a bound rather
than a coincidence.

Then the question is inverted, which is what makes the gap falsifiable:

| | |
|---|---|
| v_sat that satisfies criterion A exactly | **6.50 × 10⁶ cm/s** |
| v_sat Eq. (6) gives at this bias | 6.63 × 10⁷ cm/s |
| Dorgan *et al.* measured on SiO₂ | ~1 – 3 × 10⁷ cm/s |

**Criterion A needs a v_sat below the measured band.** And a sweep over four
phonon energies spanning the physically available range reaches 1.2491 at
best:

| ħΩ (eV) | source | v_sat (cm/s) | g_ds factor |
|---|---|---|---|
| 0.059 | SiO₂ surface polar phonon, low mode | 3.91e7 | **1.2491** |
| 0.100 | Feijoo *et al.* 2019 fit | 6.63e7 | 1.1361 |
| 0.149 | SiO₂ surface polar phonon, high mode | 9.88e7 | 1.0893 |
| 0.196 | graphene intrinsic optical phonon | 1.30e8 | 1.0679 |

**No phonon energy in the physical range satisfies criterion A.** That table
also discharges the 2026-09-29 item "the remote-polar-phonon cap on SiO₂ is
asserted, not computed", which 2026-10-04 predicted this work would sharpen:
Eq. (6) is linear in Ω, so the phonon energy *is* the v_sat scale, and the
cap's effect on g_ds is now a measured 1.07–1.25 rather than an assertion.

β = 1 vs. Dorgan *et al.*'s own best-fit β = 2: at the largest field this
device reaches, β = 2 gives a velocity 24.6–26.0 % **higher**, i.e. it
saturates *less* hard. The β = 1 choice is the generous one in this direction
too, so the shortfall is a bound from both sides.

---

## 3. Criterion A: missed, and the reason is in Chapter 4's other chapter

| | resistor | saturated | factor |
|---|---|---|---|
| g_ds at the resistor's peak-f_T bias | 3.811e-2 S | 3.355e-2 S | **1.1361** |
| g_ds at the saturated model's peak | 3.767e-2 S | 3.354e-2 S | 1.1230 |
| required by criterion A | | | **4.8483** |

23.4 % of the requirement. But μS/L = 0.3036 at this bias, which alone would
divide the **channel** conductance by 1.699 — so most of even that 1.14 is
being eaten before it reaches g_ds. At the 40 µm / 8-finger geometry
R_c,total = 15.0 Ω of a ≈ 26 Ω device:

> **Roughly half of g_ds is a contact resistance that no saturation mechanism
> can touch. The contacts, not the channel, now cap g_ds.**

This is the useful part of a missed criterion. 2026-10-04 concluded that the
R_g thread was correct work aimed at a non-binding constraint and pointed the
next session at the drain-field physics. The drain-field physics is now in,
behaves correctly, and is **also** not the binding constraint. The binding
constraint is the access resistance — which is Chapter 4 §4.5, §4.7 and §4.8,
the part of this thesis with the most real content and the one place it has a
negative result worth having. The RF thread has now been chased through three
candidate levers (R_g, C_gd, v_sat) and all three have come back
non-binding, each pointing at the contacts.

A note on what is *not* claimed. f_max/f_T reaches 1.3844 at V_ds = 1 V, inside
Feijoo *et al.*'s 1.3–1.4 band. **That is not corroboration and is explicitly
not offered as any.** 2026-10-04 recorded the trap: peak f_T here scales ×19.8
over a ×20 drain-bias range and the literature band spans a comparable factor,
so a model with a knob that size agrees with it *somewhere*. The claim in
Section 4 below is the **sign**, which is bias-independent; where the ladder
happens to cross 1.3 is a fact about where the ladder was stopped.

---

## 4. Criterion B: met, and it is the one clean win

| V_ds | resistor | saturated |
|---|---|---|
| 0.05 | 0.644461 | 0.663187 |
| 0.10 | 0.642966 | 0.683186 |
| 0.20 | 0.640511 | 0.728061 |
| 0.50 | 0.633440 | 0.911727 |
| 1.00 | 0.621775 | 1.384374 |

Δ = −0.0227 → **+0.7212**, a reversal 31.8× the magnitude of the wrong-signed
slope it replaces. A saturating device's f_max/f_T rises with drain bias; this
one now does. The sign was the part of 2026-10-04's diagnosis that was a
statement about *mechanism* rather than about magnitude, and it is the part
that came out right.

---

## 5. Criterion C found something it was not written to find

The overlapping-limit check cannot be an assertion of bitwise agreement,
because the saturated model's S → 0 limit is the *arithmetic* average of n
while the resistor model takes a *harmonic* one. So four artefacts are made to
converge as V_ds → 0, which separates three mechanisms:

| mechanism | size at V_ds = 0.1 V | order in V_ds |
|---|---|---|
| quadrature rule (`np.mean` over 50 samples vs. trapezoid) | −4.3e-05 % | 2.04 |
| **profile domain (what the contacts drop)** | **+0.253555 %** | **1.00** |
| Jensen gap (⟨1/n⟩ vs. 1/⟨n⟩) | +0.000376 % | 2.00 |

**Two predictions of mine were wrong and the decomposition is what said so.**
Both are in the committed transcript rather than edited out.

**(i)** I predicted `np.mean` over `np.linspace(0, Vds, 50)` carried an O(1/N)
sample-mean bias. It does not — a mean of equally spaced samples including
both endpoints is exact for a linear integrand and second order for a smooth
one. `n_segments = 50`, open since 2026-09-26 and called *load-bearing* on
2026-10-04 because g_ds is a V_ds difference of that quadrature, is hereby
measured at 4.3e-05 % and is **not a problem at the committed biases**. A
nine-day item closes with a negative result.

**(ii)** I predicted the residual was the Jensen gap. It is 674× *not*. The
committed model profiles n(V_ch) over the **full V_ds**, but over half of V_ds
is dropped across the contacts and never appears across the channel at all, so
it evaluates the channel over a potential range about twice too wide. That is
first order in V_ds and it is the real content of criterion C.

No committed number is withdrawn. The total is one-signed, 0.2539 % at
V_ds = 0.1 V, and in the direction that makes I_d *larger* — so every shortfall
against literature in Chapter 4 is, if anything, understated. What changed is
that the bound is measured, its mechanisms are separated and ordered, and **the
larger one was not the one anybody here was looking for**. The three
mechanisms are separable precisely because their orders in V_ds differ (2, 1,
2); a single aggregate residual would have shown one number and attributed it
to whichever cause was in mind.

---

## 6. What the mutation harness changed

First run: **5 of 7**. Both survivors were real.

**M6 survived.** M6 multiplies v_sat by π/2 and does nothing else — no
structural change, no exception, no sign flip, a 57 % shift in the one knob
that decides whether criterion A was measured or fitted. Every check in the
battery as first committed was a sign, a zero or a ratio, and **every one of
them is invariant under a constant rescale of v_sat.** That is the standing top
methodological item (2026-10-03: a PASS/FAIL at zero is silent about
magnitude) landing on the most load-bearing quantity in the item being closed.

**M7 survived.** M7 clips the reported physical ceiling to 0.5 instead of
reporting it. X6 asserted only `max < 1.0`, which a guard that lies about its
own report satisfies exactly as well as one that does not — the same family as
2026-10-04's census renderer that detected and then crashed: a detection
strictly weaker than it looks from outside, invisible until something makes it
lie.

Three checks added: **X7**, a magnitude pin on v_sat against an anchor
computed outside the module at 30 significant digits (so the check is not
Eq. (6) agreeing with itself); **X6b**, the reported ceiling must equal an
independent recomputation; **X6c**, it must respond across the V_ds ladder.
Battery 23/23 → 26/26, harness 7 of 7, control B still survives.

X7's own first committed form used a 6-digit hand value of Ω and **failed at
rel. err 1.9e-04** — the check caught my arithmetic before it caught any
mutant. Recorded in the module.

**Read M6's killer column.** M6 is killed by **exactly one** check — X7 — and
X7 is the only check in the battery that retains a magnitude. That is the
identical shape 2026-10-04 reported for M4/R7 in the f_max decomposition, now
reproduced independently on a different module, a different quantity and a
different defect. **Two independent instances is no longer an anecdote**, and
the item should be mechanised rather than given a third worked example.

---

## 7. Methodological note, continuing the series

09-28: prose is a detector. 09-29: a mutation that does not arrive is
indistinguishable from a system that does not respond. 09-30: a control has to
sit where the failure enters. 10-01: and name what it compares against in a way
that cannot drift. 10-02: when it agrees, that is a fact about two artefacts,
not about the world. 10-03: a PASS/FAIL at zero is silent about magnitude.
10-04: and a quantity can be computed correctly under the wrong name.

**10-05: AND A CRITERION WITH A NUMBER IN IT CAN BE MISSED IN A WAY THAT
LOCATES THE REAL CONSTRAINT, WHICH A CRITERION WITHOUT ONE CANNOT.**

Had 2026-10-04 written the item as "add velocity saturation to
`transfer_characteristic()`", today would have added it, observed f_max/f_T
rise, noted agreement with Feijoo at V_ds = 1 V, and closed the item. Every one
of those steps would have been defensible and the conclusion would have been
wrong — because the thing that would have gone unmeasured is the 1.1361, and
the 1.1361 is what says the contacts are binding. **The number in the criterion
is what converted a successful addition into a located constraint.**

The corollary is about where such numbers come from. 4.8483 was not available
to the session that wrote the criterion until it had *decomposed* f_max rather
than computed it — it is the by-product of asking a denominator which of its
terms mattered. So the practice that makes items falsifiable is the same
practice that found the 10-04 result: attribution, not correctness. An item
that can only be *done* is an item written before the attribution was.

---

## 8. Not yet covered (candidates for future runs)

- **Rewire `transfer_characteristic()` to the saturated form, or decide not
  to** — created today, and the direct successor to the top item. The cost is
  now priced: five modules read it, their transcripts move, and the profile-
  domain correction of §5 moves them a further one-signed 0.25 %. The decision
  is whether Chapter 4's committed RF numbers should be re-derived under a
  model whose f_max/f_T has the right sign but whose g_ds is still
  contact-limited. **Recommended: yes, but as its own session**, because it is
  a rule change under the 2026-10-02 rule and nothing else should share the
  commit.
- **MECHANISE the magnitude item** — created 2026-10-03, worked example
  2026-10-04, **second independent instance today**, still not mechanised.
  Two different modules, two different quantities, same result: the only
  detector that catches a pure-magnitude defect is the only check that keeps a
  magnitude. The mechanisable form is a repository-wide census of every
  MUST_CHANGE that records its response size, and today is the last session at
  which "one more worked example" is a defensible answer.
- **Why does the peak-f_T bias MOVE BRANCH between the two models?** — created
  today. The resistor model peaks at V_g = −1.1128 V (hole side), the
  saturated model at +2.8571 V (electron side). Both were reported rather than
  one being chosen, but the branch flip is not explained, and a figure of merit
  whose optimum jumps 4 V between two models of the same device is either a
  real asymmetry or an artefact of where g_m peaks. Chapter 4's electron–hole
  asymmetry work (the `t'` item, 09-28) is the natural place to look.
- **The contacts now cap g_ds — so what is the f_max ceiling at R_c = 0?** —
  created today, and the obvious next counterfactual, exactly parallel to
  2026-10-04's R_g = 0 and R_s = 0 rows. If a perfect contact *also* falls
  short, the RF thread is closed by three independent exhausted levers rather
  than two, and §7.8.1a's structural verdict is final rather than provisional.
- **A β = 2 implementation** — created today. β = 2 destroys the separability
  of Eq. (4) and needs numerical quadrature on x rather than on V. Bounded
  today at 24.6–26.0 % in velocity at the largest field and shown to move the
  answer the *wrong* way, so this is a completeness item, not a lever.
- **Whether Chapter 4's `n_puddle` actually reproduces Section 3.3's measured
  6.45 kΩ/sq floor** — created 09-29, **untouched for seven sessions**, and
  today raises it again: n_puddle sets n at the Dirac point, n sets v_sat
  through Eq. (6), and the peak-f_T bias of the resistor model sits near the
  Dirac point.
- **A second anchor for `Delta_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, now **untouched for fifteen
  consecutive sessions**.
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural
  weakness.
- **The Section 3.6 Pauli edge is absent from Chapter 6's model** (09-29);
  **which other quantities here are computed correctly under a name that claims
  more than they measure** (10-04); **which checks are RATIOS of two quantities
  the same defect would scale** (10-03); **audit every remaining check for the
  10-02 §1 fault** (10-02); **do the four existing audits have references that
  are not content-pinned** (10-01); **is `D_CHEM_PHYS` the only duplicated
  literature VALUE** (10-01); **is `plot_fermi_surface` the only uncommitted
  figure** (10-02); **where does the 6.5430 % internal collection efficiency
  come from** (09-30); a finite-temperature optical conductivity (09-29);
  angular trigonal warping (09-28); finite-temperature `n(E_F, T)` (09-28);
  whether `t'` can be excluded quantitatively as the source of Chapter 4's
  electron–hole asymmetry (09-28); Ti and Cr per-metal `Rc` recalibration
  (ResearchGate rate-limiting); a second independent edge-contact dataset
  (Lee *et al.* 2022, Wiley 403'd).

---

## References

- V. E. Dorgan, M.-H. Bae and E. Pop, "Mobility and saturation velocity in
  graphene on SiO₂", *Appl. Phys. Lett.* **97**, 082112 (2010).
  <https://poplab.stanford.edu/pdfs/Dorgan-GrapheneVsat-apl10.pdf> — measured
  v_sat ≈ 3 × 10⁷ cm/s at low density on SiO₂, falling with increasing n;
  soft-saturation fit with β = 2.
- P. C. Feijoo, F. Pasadas, M. Bonmann, M. Asad, X. Yang, A. Generalov,
  … D. Jiménez, "Does carrier velocity saturation help to enhance f_max in
  graphene field-effect transistors?", *2D Materials* **6**, 035027 (2019).
  <https://arxiv.org/pdf/1910.08304> — Caughey–Thomas form at γ = 1,
  ħΩ = 0.10 eV, and the finding that g_sd does *not* minimise at the f_max
  peak. Same group as the Feijoo *et al.* 2016 band this chapter benchmarks
  against.
- P. C. Feijoo *et al.*, *Sci. Rep.* **6**, 35717 (2016) — the f_max/f_T =
  1.3–1.4 band used throughout Chapter 4.

**Fetch record.** Both PDFs above fetched successfully. PubMed/PMC not
attempted, per the standing note that they return reCAPTCHA in this
environment (now four occurrences). No fetch failures to report.
