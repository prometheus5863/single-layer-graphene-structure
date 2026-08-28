# Device Physics Notes: The Cu Liner/Barrier-Thickness Effect on Scaled-Copper Resistivity

**Date:** 2026-08-28
**Focus:** The one piece of the copper comparison baseline in
`graphene_interconnect_model.py` that Chapter 5 (Section 5.2, Section 5.4
follow-on list) has explicitly flagged since 2026-08-23 as omitted and
"conservative in graphene's favor": the diffusion-barrier/adhesion-liner
layer that a real damascene Cu interconnect requires, whose thickness
does not scale down with linewidth and therefore consumes a growing
fraction of the wire's cross-section as W shrinks. This note gathers the
literature basis for adding that effect; the corresponding code and plot
changes are in the same-day commit to `graphene_interconnect_model.py`,
and the thesis update is in `thesis_draft/05-graphene-interconnects.md`.

WebSearch and WebFetch were both available this session and used for the
sources below.

## 1. Why the existing Cu baseline needed this

`cu_resistivity_vs_width()` as written on 2026-08-23 models only
surface/grain-boundary scattering via a single-parameter Fuchs-Sondheimer
correction (bulk mean free path 40 nm, specularity 0.6) applied to bulk
Cu resistivity (1.68 μΩ·cm). That is a real and necessary effect, but it
is not the only one, and the Chapter 5 draft already said so explicitly:
"[the Cu model] omits the liner/barrier-thickness effect... which pushes
real scaled-copper resistivity higher than a pure surface-scattering
model predicts. The comparison in this chapter should therefore be read
as an upper bound on copper's narrow-linewidth performance." This note
and the accompanying code change replace that acknowledged gap with an
actual, literature-grounded model, rather than continuing to carry it as
a caveat.

## 2. Why a liner/barrier layer exists in the first place

Copper diffuses readily into silicon and into most low-k dielectrics at
BEOL processing and operating temperatures, and does not adhere well to
oxide/low-k surfaces on its own. Standard dual-damascene Cu
interconnects therefore require a thin diffusion-barrier layer (typically
Ta or TaN) plus an adhesion liner (historically Ta or TaN itself,
increasingly a thin Co or Ru layer in more advanced schemes) lining the
trench before the Cu fill. Both layers are far more resistive than Cu
(TaN in particular is essentially a poor conductor by interconnect
standards) and, in nearly all resistance models, are treated as
non-conducting relative to the Cu core — i.e. current is assumed to flow
only through the remaining Cu cross-section, not through the liner/barrier
itself.

## 3. Quantitative thickness figures used to calibrate the model

Two sources gave the concrete numbers used below (see full search
results in this session's tool transcript; summarized here):

- **Domenichini et al. (arXiv:2406.09106, "Selecting Alternative Metals
  for Advanced Interconnects")**: states that the combined barrier +
  liner thickness in a conventional Cu damascene stack "cannot be scaled
  below a combined thickness on the order of 2-3 nm without losing
  function" (diffusion-barrier integrity and adhesion). It also states
  the underlying mechanism directly in words: "for lines with reduced
  width, barrier and liner layers (with high resistivity) occupy an
  increasingly large volume fraction of the total metallization, leaving
  less and less space for Cu, while contributing little to the
  conductance of the wire" — exactly the effect Chapter 5 flagged as
  missing. The same paper's IRDS-based roadmap table gives projected
  minimum metal pitches of 23 nm (2024/25), 20 nm (2027/28), 18 nm
  (2029), 16 nm (2031), 14 nm (2033), and 12 nm (2035) — i.e. the regime
  where this effect matters is not a distant-future concern but the
  present and immediate-next few roadmap generations.
- **A companion nanoscale-interconnect-materials review (Nanomaterials
  12(10), 1760, 2022; doi:10.3390/nano12101760)**: gives concrete per-
  conductor liner thicknesses used in its own resistance calculations —
  **3 nm** for a conventional Cu/TaN-Co liner+barrier stack, versus much
  thinner liners for alternative conductors evaluated in the same study
  (Ru: 0.3 nm; Co: 1 nm; W: no liner required at all, "linerless
  deposition"). It scans total drawn linewidth (conductor + 2× liner)
  from 48 nm down to 10 nm and reports that "below about 20 nm, the
  superiority of Cu in resistance is significantly weakened" as the
  liner/barrier's fixed thickness becomes a large fraction of the total
  width. This 3 nm figure — for the same conventional Ta/TaN-barrier,
  Co-liner Cu stack this thesis's Cu baseline is meant to represent — is
  the value used below; it is also consistent with the 2-3 nm functional
  floor reported by Domenichini et al.
- Neither source gives a closed-form "effective resistivity including
  barrier thickness" equation; both compute it via full resistance/Monte
  Carlo simulation rather than a simple analytic correction. The model
  added here is therefore an original (if standard-in-spirit) simplified
  analytic construction built from their reported physical picture and
  thickness numbers, not a formula copied from either source — this is
  stated explicitly in the code docstring and below, to avoid overstating
  how directly it is "from the literature."

## 4. The simplified area-dilution model used here

Treat the barrier/liner as (a) non-conducting, in line with both sources'
qualitative framing, and (b) lining both sides of a wire whose two
in-plane dimensions (width and height) scale together — a standard
simplification for a first local-interconnect-tier estimate, since real
local-tier aspect ratios are close to unity to a few (IRDS local-tier
aspect ratios are typically in the 1.5-2.5 range; treating width and
height as scaling together is a reasonable single-parameter proxy for
this rather than an attempt to model a specific aspect ratio exactly).
Under those assumptions, a drawn linewidth W (matching the "linewidth
(w + 2t)" convention used in the Nanomaterials review) has an effective
conducting Cu core of width/height `W_eff = W - 2·t_liner`, and:

    rho_eff(W) = rho_Cu_surface_scattering(W_eff) × (W / W_eff)^2

The first factor is the *already-implemented* Fuchs-Sondheimer surface-
scattering resistivity (Section 5.2 of the thesis), now evaluated at the
narrower effective Cu core width rather than the full drawn width — the
remaining Cu is not only smaller in cross-section, it is also more
surface-scattering-limited than a same-sized non-liner-bounded wire would
be, because it is confined by the liner interface on all sides. The
second factor is the purely geometric "area dilution": the same
resistivity-limited core current path is now measured across the full
drawn cross-section (the metric a chip designer actually cares about,
since that is the area the process allocates to the wire), so resistance
scales up by the ratio of drawn area to conducting area. This
compounding of two effects (a smaller, more surface-scattering-limited
core, diluted further by the non-conducting rim) is qualitatively what
both sources describe, even though neither publishes this exact
factored form.

`t_liner_nm = 3.0` is used as the default (from the Nanomaterials
review's Cu/TaN-Co figure, Section 3 above), with the code structured to
take an override so a thinner-liner scenario (e.g. representing an Ru- or
Co-liner-compatible dual-damascene Cu process, a real integration option
discussed in both sources) can be explored without changing the core
model. Below `W = 2 × t_liner_nm = 6 nm`, the entire drawn cross-section
is consumed by barrier/liner and there is no remaining Cu core at all;
the code flags this regime explicitly (returns `NaN`, not a large finite
number) rather than reporting a physically meaningless resistivity, and
`summary_numbers()` states this floor explicitly rather than only showing
it as a break in the plotted curve.

## 5. Expected qualitative effect on the Chapter 5 crossover result

The 2026-08-23 crossover analysis found graphene beats the (liner-free)
Cu baseline only above W ≈ 130 nm at realistic (p = 0.15) edge quality.
Since the liner effect makes Cu resistivity rise *faster* than pure
surface scattering as W shrinks (the `(W/W_eff)^2` factor diverges as W
approaches 2×t_liner, well before surface scattering alone would look
severe), including it should push that crossover width down substantially
— i.e. make graphene's realistic-edge-quality advantage kick in at a
narrower, but more roadmap-relevant, linewidth than the previous
(deliberately conservative) baseline suggested. The exact new crossover
width is computed and reported in `graphene_interconnect_model.py`'s
`summary_numbers()` output and discussed quantitatively in the updated
Chapter 5 draft, rather than estimated here.

## 6. Open items after this note

- This model treats the liner/barrier as strictly non-conducting. In
  reality TaN and especially thin Co/Ru liners carry a small but non-zero
  fraction of current; a more complete treatment would use a parallel-
  conduction (liner + core) model rather than assuming zero liner
  conductance. Flagged as a refinement, not pursued today, since both
  source papers treat the liner as the *dominant* additional-resistance
  mechanism specifically because its resistivity is so much higher than
  Cu's that the non-conducting approximation is standard in this
  literature (both sources' own qualitative framing supports this
  simplification, even though their own full simulations do include
  liner conduction).
- A thinner-liner (e.g. Ru- or Co-liner-enabled, ~1 nm or less) Cu
  scenario is directly explorable with the new `t_liner_nm` parameter but
  not separately plotted/discussed today; worth a follow-on comparison
  since both source papers flag thin-liner-compatible Cu integration as
  an active industry direction distinct from switching the conductor
  material entirely.
- The graphene-all-around-metal (graphene-wrapped Co/Cu core) hybrid
  picture noted in Section 5.4 as a distinct follow-on item remains
  unaddressed by this note — that is a different device architecture
  (graphene as liner/cap on a conventional metal core) from the pure-Cu-
  liner-effect correction made here, and is still a candidate for a
  future session.

## References

- Domenichini, F. et al., "Selecting alternative metals for advanced
  interconnects," arXiv:2406.09106 (2024). https://arxiv.org/html/2406.09106v1
- "Mechanisms of Scaling Effect for Emerging Nanoscale Interconnect
  Materials," *Nanomaterials* 12(10), 1760 (2022).
  https://doi.org/10.3390/nano12101760 /
  https://www.mdpi.com/2079-4991/12/10/1760
- Murali et al., arXiv:0906.0924 and Naeemi & Meindl, IEEE EDL 28, 428
  (2007) — as previously cited in `notes/2026-08-23-*.md`, for the
  graphene-side model this note's Cu-side correction is compared against.
