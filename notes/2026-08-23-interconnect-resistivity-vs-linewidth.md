# Device Physics Notes: Graphene Interconnect Resistivity vs. Linewidth

**Date:** 2026-08-23
**Focus:** How graphene's resistivity scales as interconnect linewidth is
shrunk toward sub-20 nm dimensions, why edge scattering sets a
width-dependent mean free path, and how that compares with scaled copper
— the other main "device application" pillar for this thesis besides the
GFET contact-resistance/quantum-capacitance/RF work already covered
(see `notes/2026-08-21-*.md`, `notes/2026-08-22-*.md`).

## 1. Why this matters for device work

Copper interconnect resistivity rises sharply as linewidth scales below
~20 nm because of two compounding effects: electron surface/grain-boundary
scattering (the classic Fuchs-Sondheimer / Mayadas-Shatzkes thin-film
result) and the need for a diffusion-barrier liner that eats an
increasing fraction of the line's cross-section as W shrinks, since the
liner thickness does not scale down proportionally. This "interconnect
resistivity crisis" is one of the most-cited motivations (alongside pure
logic-device scaling) for looking at 2D materials — including graphene —
as an alternative or capping/liner material for back-end-of-line (BEOL)
wiring. This note focuses on graphene's own linewidth-dependent
resistivity, which is the natural counterpart to the GFET contact/channel
work already in this repo: contacts and channel set transistor
performance, interconnects set how fast signals move between transistors,
and both are evaluated by semiconductor companies using the same
figure of merit — resistance (or resistivity) as a function of a shrinking
critical dimension.

## 2. Physical origin of width-dependent graphene resistivity

Unlike bulk copper, a graphene ribbon narrow enough to be interesting as
a local-level interconnect (tens of nm down to single-digit nm) has a
non-negligible fraction of its carriers scattering off the *edges* of the
ribbon rather than only off bulk defects/phonons/impurities. This is
formally analogous to surface scattering in thin metal films, but the
relevant boundary is the ribbon edge rather than the film surface. Three
scattering channels compete, combined via Matthiessen's rule
(1/λ_eff = 1/λ_bulk + 1/λ_edge + 1/λ_impurity):

- **Bulk/phonon-limited scattering** — sets the highest achievable mean
  free path (MFP), largely width-independent, and is what the "best"
  measured GNR samples approach.
- **Line-edge-roughness (LER) / edge scattering** — the width-dependent
  term. A phenomenological *specularity parameter* p (0 = fully diffuse
  edge scattering, randomizing momentum; 1 = fully specular, i.e. mirror-
  like and non-resistive) is used to interpolate between these limits,
  following the same approach originally developed for thin-film surface
  scattering (Fuchs-Sondheimer) and adapted to graphene nanoribbons (GNRs)
  by Naeemi & Meindl and others. A commonly used approximate form for the
  edge-limited contribution to the mean free path is

      1/λ_edge ≈ [(1 - p) / (1 + p)] × (1 / W)

  so that λ_edge shrinks roughly linearly as the ribbon width W shrinks,
  and the edge contribution vanishes only in the fully specular (p → 1)
  limit. Real graphene edges produced by plasma/lithographic patterning
  are rough at the atomic scale (etched rather than exfoliated edges),
  so measured behavior corresponds to p closer to 0 (diffuse) than to 1.
- **Impurity/substrate scattering** — resist residue, substrate charged
  impurities, and dangling bonds/oxidation at the cut edge; empirically
  significant and, per the measurements below, comparable in magnitude
  to the pure LER contribution for lithographically patterned GNRs.

## 3. Representative experimental numbers

The most directly relevant experimental dataset is Murali, Brenner, Yang,
Beck & Meindl, "Resistivity of Graphene Nanoribbon Interconnects"
(Georgia Tech, IEEE EDL, arXiv:0906.0924):

- Few-layer GNRs with 18 nm < W < 52 nm, patterned by e-beam lithography
  + O2 plasma etch from exfoliated flakes on 300 nm SiO2, measured at
  300 K.
- Measured 3D resistivity clustered at **15-25 µΩ·cm** across that width
  range — about **2-3x higher than the phonon-scattering-limited
  intrinsic graphene resistivity of ~1.2 µΩ·cm** (at carrier density
  n = 5e12 cm^-2) and about **2-3x higher than the ITRS-2007-projected Cu
  resistivity** at the same linewidths.
- The *best* individual GNR (out of ten measured in parallel per device)
  at a given width had resistivity **comparable to Cu at that width** —
  i.e. sample-to-sample variation from LER/impurity differences is large
  enough that the best-case and average-case GNR performance differ by
  a factor of several.
- Scattering-mechanism decomposition (via Matthiessen's rule, comparing
  mobility before/after the plasma-etch step that converts a 2D flake
  into a narrow GNR) attributed the W=22 nm devices' effective mobility
  (4,000-8,000 cm^2/V·s) to a mix of impurity-limited mobility
  (2,500-19,000 cm^2/V·s, from an estimated impurity density of
  n_i = 2-19e11 cm^-2) and LER-limited mobility (6,000-9,000 cm^2/V·s) —
  i.e. at W ~20 nm the two mechanisms are comparably important, neither
  cleanly dominates.
- Breakdown current density was high (5-20e8 A/cm^2), noted as evidence
  of good electromigration robustness relative to Cu — a secondary but
  real interconnect-relevant advantage independent of resistivity.

Complementary theoretical/projection context (Naeemi & Meindl,
"Conductance modeling for graphene nanoribbon (GNR) interconnects," IEEE
EDL 28, 428 (2007), cited in the above): idealized single-layer GNRs are
projected to *undercut* 1:1-aspect-ratio Cu resistance-per-length only
once W is pushed below roughly **8 nm**, i.e. the crossover where
graphene's advantage becomes clearly favorable is a genuinely aggressive,
near-atomic linewidth — consistent with graphene's main projected
interconnect advantage being at the most advanced/narrowest tiers of the
BEOL stack (local-level wiring), not global/mid-level wiring where wider
lines are used and Cu's bulk conductivity advantage still wins.

A 2024 result on a different but related architecture — a "graphene-
all-around" Cobalt interconnect (graphene wrapped around a Co core,
rather than a pure graphene line) — reported a **27% resistance
reduction** versus a bare Co interconnect of the same geometry, in a
process compatible with standard BEOL integration (Nano Letters 2024 /
PMC10870778). This is a different, nearer-term commercialization path
(graphene as a low-resistance liner/cap on a conventional metal, rather
than graphene as the sole conductor) worth distinguishing from the
"pure GNR wire" picture above when discussing near-term vs. long-term
device relevance in the thesis.

## 4. Why linewidth scaling (not just intrinsic mobility) is the right lens

It would be easy to cite graphene's very high intrinsic carrier mobility
(demonstrated in the fundamental-properties scripts already in this
repo, e.g. `graphene_transport_properties.py`) and conclude graphene is
an obviously superior conductor. The interconnect literature makes clear
that this is the wrong comparison for a real BEOL wire: what matters is
resistivity (or resistance-per-length) *at the linewidth that a given
technology node actually uses*, and at those width scales edge/LER
scattering — not intrinsic mobility — is frequently the dominant, and
most process-sensitive, resistivity-limiting mechanism. This mirrors the
lesson from the contact-resistance notes (2026-08-21): a "textbook"
material-level number (mobility there, resistivity here) is not the
number semiconductor device/interconnect engineers actually use to
compare technologies; the geometry- and process-limited number is.

## 5. Planned follow-on work

- Implement a compact resistivity-vs-linewidth model in this repo
  combining a literature-calibrated bulk mean free path with the
  specularity-parameter edge-scattering term above (Matthiessen
  combination), and plot it against a simple Cu resistivity-vs-linewidth
  trend line for comparison — added alongside these notes, see
  `graphene_interconnect_model.py`.
- A future extension could add the graphene-all-around-metal (liner/cap)
  picture as a second, distinct model, since the 2024 Co result above
  suggests that may be the more near-term-relevant integration path than
  a pure-graphene wire.

## References

1. Murali, R., Brenner, K., Yang, Y., Beck, T., Meindl, J. D.
   "Resistivity of Graphene Nanoribbon Interconnects." arXiv:0906.0924.
   https://arxiv.org/pdf/0906.0924 (published as IEEE Electron Device
   Letters, 2009).
2. Naeemi, A., Meindl, J. D. "Conductance modeling for graphene
   nanoribbon (GNR) interconnects." IEEE Electron Device Letters 28,
   428-431 (2007).
3. "Modeling of Edge Scattering in Graphene Interconnects." IEEE Xplore
   document 8355485. https://ieeexplore.ieee.org/document/8355485
   (specularity-parameter edge-scattering model background; abstract-
   level access only in this session).
4. "Graphene-All-Around Cobalt Interconnect with a Back-End-of-Line
   Compatible Process." Nano Letters (2024).
   https://pubs.acs.org/doi/10.1021/acs.nanolett.3c04833 (open-access
   mirror: https://pmc.ncbi.nlm.nih.gov/articles/PMC10870778/).
5. "Prospects and challenges of compound conductors for advanced
   interconnect applications." Journal of Applied Physics 138, 090902.
   https://pubs.aip.org/aip/jap/article/138/9/090902/3361546
6. Semiconductor Industry Association, International Technology Roadmap
   for Semiconductors (ITRS) 2007 Cu resistivity projections, as cited
   in reference 1.

**Search availability note:** WebSearch was available this session and
was used for this research (queries on graphene interconnect resistivity
scaling, edge-scattering/specularity models, and recent 2024-2025
developments).
