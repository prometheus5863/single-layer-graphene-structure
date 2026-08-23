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
the comparison baseline. This copper model is deliberately conservative
in graphene's favor: it omits the liner/barrier-thickness effect
described in Section 5.1, which pushes real scaled-copper resistivity
higher than a pure surface-scattering model predicts. The comparison in
this chapter should therefore be read as an upper bound on copper's
narrow-linewidth performance, making any graphene crossover point
reported below a conservative (if anything, optimistic-for-copper)
estimate.

## 5.3 Results

Running the model (`graphene_interconnect_model.py`, reproducible via
`python graphene_interconnect_model.py`) gives two qualitatively
different outcomes depending on assumed edge quality:

- **Realistic, diffuse edges (p = 0.15):** graphene resistivity crosses
  below the copper model at approximately **W ≈ 130 nm** — i.e. for
  lithographically-patterned edges of the roughness typical of current
  processing, graphene only becomes the lower-resistivity conductor at
  linewidths wide enough that they are not representative of the
  aggressively-scaled local-interconnect tiers where a resistivity
  advantage would matter most.
- **Idealized, near-specular edges (p = 0.9):** graphene resistivity
  stays below the copper model across the entire modeled range (W = 2-500
  nm) — no crossover is needed because graphene already wins at every
  width, consistent with the more optimistic literature projections
  (e.g. Naeemi & Meindl's idealized single-layer-GNR crossover at very
  narrow, near-atomic linewidths) that assume near-ideal edge
  termination.

The practical conclusion — consistent with the contact-resistance finding
in Chapter 4 — is that **graphene's interconnect advantage over copper is
conditional on edge/process quality, not automatic.** The two curves
bracket a wide range of possible real-world outcomes, and the dominant
open engineering question for graphene interconnects is not "can graphene
in principle beat copper" (yes, in the idealized limit) but "can edge
roughness be controlled well enough, at the linewidths that matter, to
realize that advantage in a manufacturable process."

## 5.4 Planned follow-on work

- Add the graphene-all-around-metal (liner/cap on a conventional
  copper or cobalt core, rather than a pure graphene wire) picture as a
  second, distinct model — a 2024 result on a graphene-wrapped cobalt
  interconnect (27% resistance reduction vs. bare cobalt) suggests this
  hybrid integration path may be nearer-term-relevant than a pure-GNR
  wire, and is worth treating as a separate case rather than conflating
  with the results above.
- Extend the copper comparison model to include the liner/barrier-
  thickness effect explicitly (currently omitted, and noted above as
  making the current copper baseline conservative), which would tighten
  the graphene-favorable crossover width estimate.
- Connect this chapter's edge-scattering specularity framework explicitly
  back to the contact-resistance specularity/mode-counting discussion in
  Chapter 4 as a unifying methodological point for the thesis discussion
  chapter (Chapter 7): several of graphene's device-relevant limitations
  (contact resistance, interconnect resistivity) trace to the same
  underlying physical picture — an atomically thin conductor's
  sensitivity to boundary/interface quality — rather than being
  unrelated, chapter-specific effects.

See `notes/2026-08-23-interconnect-resistivity-vs-linewidth.md` for full
citations and additional literature discussion underlying this chapter.
