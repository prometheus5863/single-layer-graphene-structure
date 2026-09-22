# Chapter 5: Graphene Interconnects — Resistivity vs. Linewidth

*Draft section — this chapter was "not started" as of the Chapter 1
status table; this is the first drafted content, written after the
corresponding research notes and computational model were completed (see
`notes/2026-08-23-interconnect-resistivity-vs-linewidth.md` and
`graphene_interconnect_model.py`). Subsequent sessions may extend this
chapter (e.g. with the graphene-all-around-metal liner/cap picture noted
as a follow-on item below).*

## 5.1 Motivation: the interconnect problem is a different problem from the transistor problem

Chapter 4 established that a graphene field-effect transistor's realistic
performance is set less by graphene's intrinsic properties (established
in Chapters 2-3) than by device-level, geometry-dependent effects —
contact resistance, quantum capacitance, absence of current saturation.
Interconnects present a structurally similar but physically distinct
version of the same lesson. Copper interconnect resistivity in a modern
back-end-of-line (BEOL) stack rises sharply as linewidth is scaled below
roughly 20 nm, driven by two effects: electron surface/grain-boundary
scattering (captured empirically by Fuchs-Sondheimer- and Mayadas-
Shatzkes-type thin-film models) and a diffusion-barrier liner whose
thickness does not scale down proportionally with the line, so it
consumes an increasing fraction of the available cross-section at small
linewidths. This "interconnect resistivity crisis" is one of the primary
motivations — independent of any transistor-level motivation — for
evaluating graphene and other 2D materials as alternative conductors,
liners, or caps in advanced BEOL wiring.

The question this chapter addresses is narrower and more specific than
"is graphene a better conductor than copper?" (a question best answered
by comparing bulk mobilities, where graphene wins trivially and
uninformatively). The right question, matching how a semiconductor
company would actually evaluate this, is: **at the linewidth a given
process node actually uses, and accounting for the fact that a real
graphene ribbon is patterned with lithographically rough (not
atomically smooth) edges, is graphene's resistivity lower than scaled
copper's?**

## 5.2 Physical model

Graphene ribbon resistivity is modeled here via a Matthiessen's-rule
combination of three scattering mechanisms limiting the carrier mean free
path λ: bulk (phonon/substrate) scattering, edge scattering, and residual
impurity/disorder scattering:

    1/λ_eff(W) = 1/λ_bulk + 1/λ_edge(W) + 1/λ_impurity

Edge scattering is modeled with a phenomenological specularity parameter
p ∈ [0, 1] (0 = fully diffuse edge scattering; 1 = fully specular,
non-resistive), following the approach used in the graphene-nanoribbon-
interconnect literature (adapted from the classical Fuchs-Sondheimer
thin-film treatment):

    1/λ_edge(W) = [(1 - p)/(1 + p)] · (1/W)

so that edge scattering grows in importance — and resistivity rises — as
W shrinks, with the rate of increase controlled by how rough (diffuse,
low p) or smooth (specular, high p) the ribbon edge is. This is directly
analogous to the specularity treatment of metal-graphene contact
resistance in Chapter 4 and the underlying notes: in both cases, an
idealized bulk/interior property (mobility here, quantum-limited contact
resistance there) is degraded by an interface- or edge-localized
scattering mechanism that a purely bulk material-property measurement
would miss entirely.

`graphene_interconnect_model.py` implements this model with two
literature-calibrated regimes, deliberately kept distinct rather than
collapsed into a single fit:

- A **best-achievable, edge-quality-limited** family of curves (varying
  p), calibrated against the *best individual* GNR resistivity reported
  by Murali et al. (~3× the phonon-scattering-limited intrinsic
  resistivity of 1.2 μΩ·cm) — representing what a well-patterned, low-
  impurity sample can achieve as a function of edge specularity alone.
- An **empirical reference band** showing the *average* measured
  resistivity cluster from the same study (15-25 μΩ·cm, 18 nm < W <
  52 nm) — reflecting the additional process-induced line-edge-roughness
  and impurity scattering present in typical (not best-case) lithographically
  patterned samples, which a clean specularity-only model does not
  capture.

Keeping these two regimes separate, rather than forcing a single p value
to reproduce the average measured data, is a deliberate modeling choice:
it makes explicit that the gap between "best demonstrated" and "typical"
graphene interconnect performance is currently a process/fabrication-
quality gap, not a fundamental physical limit — the same framing used in
Chapter 4 for contact resistance, where sample-to-sample variation (Pd
~110 Ω·µm vs. more typical values several times higher) was similarly
attributed primarily to interface cleanliness and process order rather
than an intrinsic material ceiling.

A simplified Fuchs-Sondheimer-style model for scaled copper (bulk
resistivity 1.68 μΩ·cm, mean free path 40 nm, specularity 0.6) is used as
one comparison baseline. As originally noted here, this surface-
scattering-only model is conservative in graphene's favor: it omits the
liner/barrier-thickness effect described in Section 5.1. That gap is now
closed. `graphene_interconnect_model.py`'s `cu_resistivity_with_liner()`
(added 2026-08-28; see
`notes/2026-08-28-copper-liner-barrier-thickness-effect.md`) adds a
second copper baseline that treats the diffusion-barrier/adhesion-liner
as non-conducting and consuming a fixed thickness (3 nm per side,
calibrated to a conventional Ta/TaN-barrier, Co-liner Cu stack — see the
Nanomaterials 12(10), 1760 (2022) resistance calculations, consistent
with the ~2-3 nm functional floor reported by Domenichini et al.,
arXiv:2406.09106) from *both* in-plane dimensions of the drawn wire. The
remaining Cu core of width `W_eff = W − 2·t_liner` is both narrower (and
so more surface-scattering-limited itself) and dilutes the effective
resistivity measured across the full drawn cross-section by a further
factor of `(W / W_eff)²`. Below `W = 2·t_liner = 6 nm` the liner consumes
the entire drawn cross-section and no Cu core can exist at all; the
model flags this explicitly rather than reporting a finite number there.
This construction is an original simplified analytic model built from
the two sources' reported physical picture and thickness figures, not a
formula copied from either — neither publishes a closed-form
liner-thickness correction, both instead using full resistance/Monte
Carlo simulation. Both copper baselines — surface-scattering-only and
surface-scattering-plus-liner — are now carried through this chapter side
by side, rather than replacing one with the other, since they represent
genuinely different (bare-wire vs. process-realistic) questions.

## 5.3 Results

Running the model (`graphene_interconnect_model.py`, reproducible via
`python graphene_interconnect_model.py`) gives two qualitatively
different outcomes depending on assumed edge quality, and — as of the
2026-08-28 update — a materially different picture once the copper
baseline itself is made more realistic:

- **Against the surface-scattering-only Cu baseline:** at realistic,
  diffuse edges (p = 0.15), graphene resistivity crosses below copper at
  approximately **W ≈ 130 nm** — i.e. graphene only wins at linewidths
  wide enough that they are not representative of aggressively-scaled
  local-interconnect tiers. At idealized, near-specular edges (p = 0.9),
  graphene stays below copper across the entire modeled range (W = 2-500
  nm), consistent with the more optimistic literature projections (e.g.
  Naeemi & Meindl's idealized single-layer-GNR crossover at very narrow,
  near-atomic linewidths).
- **Against the liner-aware Cu baseline (new):** at realistic, diffuse
  edges (p = 0.15), graphene now stays *below* the liner-aware copper
  model across the entire valid comparison range — no crossover at all.
  Representative points from the model: at W = 20 nm, graphene (p = 0.15)
  is 3.82 μΩ·cm versus liner-aware copper at 4.90 μΩ·cm; at W = 14 nm,
  4.87 μΩ·cm versus 9.00 μΩ·cm; at W = 12 nm, 5.45 μΩ·cm versus
  13.44 μΩ·cm. The gap widens sharply as W approaches the liner's
  6 nm cross-section-consumption floor, where the liner-aware Cu model
  diverges while graphene's edge-scattering model degrades only linearly.
  At idealized edges (p = 0.9), graphene stays below the liner-aware
  model as well, by an even wider margin.

**2026-08-29 update — parallel-conduction refinement and thinner (Ru/Co)
liner scenarios.** Section 5.2's liner-aware copper model treated the
liner as strictly non-conducting, a simplification flagged explicitly at
the time. `cu_resistivity_with_liner_parallel()` (added 2026-08-29; see
`notes/2026-08-29-liner-parallel-conduction-and-thin-liner-scenarios.md`)
replaces that assumption with an explicit core-plus-liner parallel-
conduction model, verified to reduce exactly to the non-conducting model
in the rho_liner → ∞ limit. Applied to the original 3 nm TaN/Co case, the
correction is nearly negligible — TaN's assumed effective resistivity is
still roughly two orders of magnitude above the Cu core's own resistivity
at relevant widths, so allowing it to conduct barely changes the result,
and the "graphene wins at every buildable linewidth" conclusion above is
unaffected. Applied to thinner, more conductive liner scenarios (Ru at
0.3 nm, Co at 1.0 nm — thicknesses from the same Nanomaterials review
cited in Section 5.2, with representative, order-of-magnitude liner
resistivities since this session's literature search located strong
qualitative but no quantitative source for the specific thin-film
values), the picture changes materially for the Co case specifically: a
thin, well-conducting liner keeps Cu resistivity comparatively low even
very close to the W = 2·t_liner cross-section-consumption floor, because
current can still flow through the liner as the Cu core area vanishes —
something the non-conducting model cannot represent at all, since it
forces resistance to diverge there regardless of what the liner is made
of. This pushes the realistic-edge (p = 0.15) graphene crossover down
from "no crossover" (bare 3 nm TaN/Co case) to approximately **W ≈ 3 nm**
against the thin-Co-liner Cu model — deep enough into the sub-nanometer
regime that it is no longer a practically relevant advantage at any
achievable edge quality. The thin-Ru-liner case, by contrast, tracks close
to the original liner-*free* baseline (crossover ≈ 123 nm, essentially
unchanged from the 130 nm liner-free figure), since a 0.3 nm liner
consumes too little cross-section to matter much either way. The
practical reading: whether graphene's process-realistic interconnect
advantage (main result above) survives once a specific liner material and
thickness are chosen is *not* a settled question — it depends materially
on which liner integration path the copper baseline assumes, with the
Co-liner scenario in particular substantially narrowing graphene's
advantage window. This is exactly the qualitative reason the
interconnect-materials literature gives for pursuing thin, low-resistance
liner materials in the first place (Section 1 of the 2026-08-29 notes),
now reproduced from this chapter's own model rather than only cited from
elsewhere.

This is a genuinely different qualitative conclusion from the 2026-08-23
draft's headline number (a 130 nm realistic-edge crossover), not merely a
refinement of it, and it is worth being precise about why: the original
130 nm figure was explicitly flagged at the time as resting on a
conservative (bare-wire, no-liner) copper baseline. Once the same
literature-reported liner/barrier consumption that Cu interconnects
actually require in a real damascene process is included, the realistic-
edge-quality graphene curve is lower than copper's at *every* linewidth
in the modeled range where a Cu wire could physically be built at all
(W ≳ 6 nm given a 3 nm liner). The practical conclusion is not that
graphene's interconnect advantage is now unconditional — this remains a
simplified analytic model on both sides (a phenomenological specularity
parameter for graphene; a non-conducting-liner, isotropic-cross-section
approximation for copper, both described in Section 5.2), not a full
transport or Monte Carlo simulation — but it does mean the dominant
open engineering question is narrower than "can graphene in principle
beat copper at process-realistic linewidths": with a realistic Cu
baseline, the model says it already does, and the remaining question is
whether lithographic edge quality can be kept near p ≈ 0.15 (not
degraded further) as linewidths shrink toward the liner-consumption
floor, and whether the non-conducting-liner approximation itself
(Section 5.2, Section 5.4) holds up against a more detailed treatment.

## 5.4 Planned follow-on work

- Add the graphene-all-around-metal (liner/cap on a conventional
  copper or cobalt core, rather than a pure graphene wire) picture as a
  second, distinct model — a 2024 result on a graphene-wrapped cobalt
  interconnect (27% resistance reduction vs. bare cobalt) suggests this
  hybrid integration path may be nearer-term-relevant than a pure-GNR
  wire, and is worth treating as a separate case rather than conflating
  with the results above. Still open — distinct from the liner-effect
  item below, which corrects the *pure-Cu* baseline rather than modeling
  graphene as a liner/cap itself.
- ~~Extend the copper comparison model to include the liner/barrier-
  thickness effect explicitly~~ — **done 2026-08-28**
  (`cu_resistivity_with_liner()`; see Section 5.2/5.3 above and
  `notes/2026-08-28-copper-liner-barrier-thickness-effect.md`).
- ~~Parallel-conduction (liner + core) refinement, and a thinner-liner
  (Ru/Co) scenario~~ — **done 2026-08-29**
  (`cu_resistivity_with_liner_parallel()`; see Section 5.3 above and
  `notes/2026-08-29-liner-parallel-conduction-and-thin-liner-scenarios.md`).
  Surfaced two narrower open items: (1) the specific liner-material
  resistivity values used (TaN 400, Ru 30, Co 20 μΩ·cm) are order-of-
  magnitude placeholders — literature search this session confirmed the
  qualitative resistivity ordering but could not access a quantitative
  thickness-resolved source (several candidates were paywalled or
  blocked by robots.txt); tightening these against a specific source is
  a direct follow-on. (2) The liner's current distribution is treated as
  spatially uniform across the whole liner "frame," not locally resolved
  near corners — a second-order refinement relative to item (1)'s
  uncertainty.
- Connect this chapter's edge-scattering specularity framework explicitly
  back to the contact-resistance specularity/mode-counting discussion in
  Chapter 4 as a unifying methodological point for the thesis discussion
  chapter (Chapter 7): several of graphene's device-relevant limitations
  (contact resistance, interconnect resistivity) trace to the same
  underlying physical picture — an atomically thin conductor's
  sensitivity to boundary/interface quality — rather than being
  unrelated, chapter-specific effects.

See `notes/2026-08-23-interconnect-resistivity-vs-linewidth.md` for full
citations and additional literature discussion underlying this chapter,
`notes/2026-08-28-copper-liner-barrier-thickness-effect.md` for the
liner/barrier-effect addition discussed in Sections 5.2-5.4 above, and
`notes/2026-08-29-liner-parallel-conduction-and-thin-liner-scenarios.md`
for the parallel-conduction refinement and thin-Ru/Co-liner scenarios
discussed in Section 5.3.


## 5.5 Conditioning: the calibration is the worst-conditioned step in this thesis

Section 4.10 introduces the condition-number diagnostic and its exact
validation, Eq. (4.22). This section applies it to Chapter 5. It changes no
number and recalibrates nothing; the results below are error bars that were
previously unstated, not corrections.

### 5.5.1 The formula is reinforcing; the calibration behind it is not

`resistivity_vs_width()` is a Matthiessen sum of three strictly positive
terms, and audited against its own arguments it behaves exactly as such:

| W (nm) | ρ (µΩ·cm) | S(bulk) | S(edge) | S(impurity) | κ |
|---|---|---|---|---|---|
| 18 | 4.093 | 0.293 | 0.662 | 0.045 | 0.662 |
| 22 | 3.600 | 0.333 | 0.616 | 0.051 | 0.616 |
| 30 | 3.009 | 0.399 | 0.541 | 0.061 | 0.541 |
| 40 | 2.602 | 0.461 | 0.469 | 0.070 | 0.469 |
| 52 | 2.321 | 0.517 | 0.404 | 0.079 | 0.517 |

Every `κ ≤ 0.66`: reinforcing at every width, as a sum of positive terms
must be. Read at that level, Chapter 5 contains no near-cancellation at all.

**But `λ_impurity` is not an argument of that formula — it is *solved* from a
single calibration point, and the solve is a subtraction:**

    residual = ρ_cal/ρ_bulk − 1 − [(1−p)/(1+p)]·(λ_bulk/W_cal)
             = 3.0000 − 1.0000 − 1.8478  =  0.1522                       (5.9)
    λ_impurity = λ_bulk / residual = 361.4 nm

Three terms of order unity producing a result of order 0.15 is the textbook
signature, and the audit measures

    κ(residual) = 19.71                                                  (5.10)

**larger than anything in Chapter 4 (11.46) or Chapter 6 (4.55). Equation
(5.9) is the worst-conditioned step in this thesis.** The parameter-level
sensitivities follow: `d ln λ_imp / d ln ρ_cal = −19.71`,
`d ln λ_imp / d ln λ_bulk = +13.14`, `d ln λ_imp / d ln W_cal = −12.14`,
`d ln λ_imp / d ln p = −3.73`. **A 5 % error in the single calibration datum
(ρ = 3.6 µΩ·cm at W = 22 nm, "best GNR" from Murali et al.) moves the
extracted impurity mean free path by 99 %.**

### 5.5.2 Result: the near-cancellation propagates, damped but not removed

Re-solving `λ_impurity` inside the audit — the chapter's *actual* dependency
chain, rather than the formula read in isolation — gives:

| W (nm) | S(ρ_cal) | S(ρ_bulk) | S(p) | S(W_cal) | S(λ_bulk) | κ |
|---|---|---|---|---|---|---|
| 18 | 0.880 | +0.120 | −0.037 | 0.542 | +0.120 | 0.880 |
| 22 | **1.000** | **0.000** | **0.000** | 0.616 | **0.000** | 1.000 |
| 30 | 1.197 | −0.197 | +0.060 | 0.737 | −0.197 | 1.197 |
| 40 | 1.384 | −0.384 | +0.118 | 0.852 | −0.384 | 1.384 |
| 52 | 1.551 | −0.551 | +0.169 | 0.956 | −0.551 | **1.551** |

The `κ ≈ 20` of (5.9) is damped by the small weight the impurity term
carries in the total ratio, but not removed: `κ(ρ)` crosses 1 at the
calibration width and reaches **1.551** at 52 nm. A prediction registered
before the audit — that `κ` would exceed 1 but stay below the
near-cancelling threshold of 3 — holds.

The methodological point is general and worth stating in Chapter 7's voice:
**a quantity can be well-conditioned in its own arguments and
ill-conditioned in the quantities those arguments were derived from.**
Auditing the visible formula alone would have returned `κ ≤ 0.66` and a
clean bill of health.

### 5.5.3 An exact result: the calibration point pins the model, and the
### sensitivity to ρ_bulk changes sign through it

The bold row above is not rounding. At `W = W_calibration = 22 nm` the model
reproduces `ρ_calibration = 3.6 µΩ·cm` **bitwise**, and

    S(ρ_bulk) = S(p) = S(λ_bulk) = 0.0        exactly                    (5.11)

— three independent exact zeros, verified as such and not merely small
(`conv = 0.0e+00`). The reason is structural: at the calibration width the
impurity term is defined as whatever makes the total come out to the
calibration datum, so it absorbs any change in the other three inputs
exactly. **At 22 nm this chapter's prediction carries no information from
`ρ_bulk`, `p` or `λ_bulk` whatsoever.**

`S(ρ_bulk)` then changes sign through that width — `+0.120` at 18 nm,
`0.000` at 22 nm, `−0.551` at 52 nm — because on the wide side an increase in
`ρ_bulk` raises the modelled bulk term but lowers the fitted impurity term by
more. A one-point-calibrated model inherits a sensitivity structure that
pivots about its calibration point, and **the further a prediction is made
from 22 nm, the more of its content comes from the single calibration datum
and the less from the physical inputs** — which is the opposite of how a
calibrated model is usually read. Section 5.3's headline comparisons are
drawn at the narrow end (18–22 nm), where this effect is weakest; the
wide-end numbers carry it.

### 5.5.4 Result: the Cu liner model contains an unremarked subtraction

`W_eff = W − 2t` (Section 5.2) is a difference of comparable quantities:

| W (nm) | W_eff (nm) | S(W) | S(2t) | κ(W_eff) | \|S_t(ρ_eff)\| |
|---|---|---|---|---|---|
| 18 | 12.0 | 1.500 | −0.500 | **1.500** | **1.167** |
| 22 | 16.0 | 1.375 | −0.375 | 1.375 | 0.852 |
| 30 | 24.0 | 1.250 | −0.250 | 1.250 | 0.550 |
| 40 | 34.0 | 1.176 | −0.176 | 1.176 | 0.379 |
| 52 | 46.0 | 1.130 | −0.130 | 1.130 | 0.276 |

`κ` grows as `W → 2t` and is mildly ill-conditioned across the whole
modelled range. At the narrow end a 1 % error in the liner thickness is a
1.17 % error in the reported effective resistivity. Section 5.2 takes
`t = 3.0 nm` from two sources that themselves differ — 3 nm for a Cu/TaN-Co
stack (*Nanomaterials* **12**(10), 1760 (2022)) against a "2–3 nm functional
floor" (Domenichini *et al.*, arXiv:2406.09106) — a spread of roughly 17 %,
which (5.12)

    δρ_eff/ρ_eff ≈ |S_t| · δt/t                                          (5.12)

converts into a **~20 % bar on ρ_eff at 18 nm**, previously unstated. This
is the one input in either chapter whose uncertainty is available from the
literature rather than hypothesised, and it is also the width range where
Chapter 5's graphene-vs-Cu crossover argument is made.

### 5.5.5 What this section does not claim

1. No number in Sections 5.2–5.4 is changed, retracted or recalibrated. The
   calibration datum is still `3.6 µΩ·cm at 22 nm` and `λ_impurity` is still
   361.4 nm.
2. `κ(ρ_GNR)` grows monotonically with `W` over the audited 18–52 nm, so a
   wider wire would be *worse* conditioned. Beyond 52 nm this is untested.
3. The parallel-conduction refinement (Section 5.3, thin Ru/Co liners) and
   the `cu_resistivity_vs_width()` Fuchs–Sondheimer baseline are **not**
   audited here. The baseline in particular contains its own `p_cu = 0.6`
   with no stated uncertainty.
4. `κ ≈ 20` on (5.9) is not evidence that `λ_impurity = 361 nm` is wrong. It
   is evidence that the value is determined by one datum to a precision
   twenty times worse than that datum's own — a statement about what the
   chapter is entitled to claim, not about what is true.

Figure: `sensitivity_audit.png`, panels (b) and (c). Machine output:
`audit_output.txt`. Pre-registration and outcome:
`notes/2026-09-22-condition-number-audit-of-chapters-4-and-5.md`.
