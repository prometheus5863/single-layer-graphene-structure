# Device Physics Notes: Parallel-Conduction Liner Refinement and Thinner-Liner (Ru/Co) Scenarios

**Date:** 2026-08-29
**Focus:** The first open item flagged in `notes/2026-08-28-copper-liner-barrier-thickness-effect.md`
(Section 6): `cu_resistivity_with_liner()` treats the barrier/liner as
strictly non-conducting, which both cited sources describe as a standard
simplification but not a physically exact one, since TaN and especially
thin Co/Ru liners carry some non-zero current in reality. This note covers
(a) the parallel-conduction model that removes that simplification and (b)
using it to evaluate the thinner Ru- and Co-liner Cu scenarios also flagged
as open on 2026-08-28. Code and plot changes are in the same-day commit to
`graphene_interconnect_model.py`; the thesis update is in
`thesis_draft/05-graphene-interconnects.md`.

WebSearch and WebFetch were both available this session and used.

## 1. What search did and did not turn up

Search queries covering TaN, Ru, and Co thin-film liner resistivity
(individually and combined) consistently surfaced the *right kind* of
literature -- barrier/liner characterization papers (e.g. "Characteristics
of ultrathin Ta and TaN films," *J. Vac. Sci. Technol. B* 20, 2328 (2002);
"Co-W Barrier Layers for Metallization of Copper Interconnects," PMC
article PMC9144600; imec's public Ru/Co advanced-interconnect liner
program; "The Resistivity Size Effect in Epitaxial Ru(0001) and Co(0001)
Films") -- but every attempt to open one and read a quantitative
thickness-resolved resistivity table failed in this session: the PMC
article resolved to a Google reCAPTCHA challenge page rather than content,
ScienceDirect abstract pages were blocked by `robots.txt`, and a re-fetch
of the arXiv HTML version of Domenichini et al. (already cited in the
2026-08-28 notes) returned only the qualitative framing already
extracted then, not its Section II.3 ("ab initio screening of alternative
metals") numeric content. This is recorded here rather than silently
worked around, per this repo's established practice (e.g. the 2026-08-26
RF-parasitics notes recording two failed IEEE Xplore/ResearchGate
fetches).

What search *did* establish, reliably, across multiple independent
sources (imec's public Ru/Co liner press materials, the Domenichini et al.
abstract/framing text that *did* load, and general familiarity with why
Ru and Co are being pursued as barrier/liner replacements at all): TaN is
a comparatively poor conductor even by barrier-material standards, and Ru
and Co are being actively developed specifically *because* they combine
much lower resistivity than TaN with the ability to be deposited thinner
(down to sub-nm for Ru) while still providing adequate diffusion-barrier/
adhesion function -- i.e. the qualitative ordering rho_TaN >> rho_Ru,
rho_Co, and the associated thinner achievable thicknesses, is well
established, even without a single clean quantitative table from this
session's searches.

## 2. Representative values used (and their status)

Given (1), the code uses representative, order-of-magnitude effective
(few-nm-thickness) resistivities, exposed as adjustable function
parameters rather than hard-coded into the model structure, so a future
session can tighten them against a specific source without touching the
model itself:

| Material | Effective thin-film resistivity (this model) | Bulk resistivity (reference) | Default thickness (per side) |
|---|---|---|---|
| TaN barrier | 400 uOhm.cm | ~125-250 uOhm.cm (bulk, phase-dependent) | 3 nm (with Co, combined -- 2026-08-28 default) |
| Ru liner | 30 uOhm.cm | ~7.1 uOhm.cm | 0.3 nm (Nanomaterials 12(10),1760 (2022)) |
| Co liner | 20 uOhm.cm | ~6.2 uOhm.cm | 1.0 nm (Nanomaterials 12(10),1760 (2022)) |

The thickness values are the same Nanomaterials review figures already
cited in the 2026-08-28 notes (Section 3) for the Cu/TaN-Co, Ru, and Co
liner cases. The resistivity values are this session's representative
estimates and are the honest caveat of this note: they capture the
correct *qualitative* size-effect physics (thin-film resistivity is always
well above the bulk value for any of these metals, and TaN is much worse
than Ru/Co) but should be read as order-of-magnitude placeholders, not
precision literature figures, until a future session locates and reads a
quantitative source directly (a clear, scoped follow-on item -- see
Section 5).

## 3. The parallel-conduction model

`cu_resistivity_with_liner()` (2026-08-28) computes the effective drawn-
cross-section resistivity assuming the liner carries *no* current:

    rho_eff = rho_core(W_eff) * (W / W_eff)^2

This is the correct rho_liner -> infinity limit of a more general
parallel-conduction picture. Treating the Cu core and the liner "frame" as
two resistors in parallel (same length L, conductances add):

    G_total = A_core/(rho_core * L) + A_liner/(rho_liner * L)
    rho_eff = A_drawn / (G_total * L) = W^2 / (W_eff^2/rho_core + (W^2 - W_eff^2)/rho_liner)

implemented as `cu_resistivity_with_liner_parallel()`. Two limiting-case
checks (both run automatically at the top of the script's `__main__`
block, not just asserted in a docstring) verify the implementation:
taking rho_liner very large (1e12 uOhm.cm) reproduces the 2026-08-28
non-conducting-liner formula to 12 significant figures at a representative
width, and setting rho_liner exactly equal to the core's own resistivity
at that width collapses the model to the plain core resistivity
(independent of t_liner_nm, as expected when the "two conductors" are
electrically identical).

## 4. Results

Full numeric output is in the `summary_numbers()` section of
`graphene_interconnect_model.py`'s run output (committed alongside the
code and plot). Headline results:

- **TaN, parallel conduction vs. non-conducting approximation:** because
  TaN's assumed resistivity (400 uOhm.cm) is still ~100x the Cu core's
  own resistivity at relevant widths, allowing it to conduct changes the
  result only marginally -- the realistic-edge-quality graphene crossover
  conclusion is unchanged (graphene stays below the TaN-liner Cu model
  across the full 5-300 nm scan range either way). This confirms the
  qualitative point already made in the 2026-08-28 notes' Section 6 (the
  non-conducting approximation is standard in this literature specifically
  *because* TaN's resistivity is so much higher than Cu's that treating it
  as an open circuit barely changes the answer).
- **Thin Ru liner (0.3 nm, 30 uOhm.cm):** because the liner is so thin, it
  consumes little cross-sectional area even though it conducts reasonably
  well, so the result tracks close to the *no-liner-effect* surface-
  scattering-only baseline rather than either the TaN case or a
  liner-free ideal -- realistic-edge graphene crosses below it at
  W ~ 123 nm, essentially unchanged from the 130 nm liner-free-Cu
  crossover found on 2026-08-23.
- **Thin Co liner (1.0 nm, 20 uOhm.cm):** a qualitatively different,
  more interesting result. Because the liner conducts well *and* the
  parallel-conduction model does not force resistance to diverge as
  W_eff -> 0 the way the non-conducting model does (current can still
  flow through the liner even as the Cu core area vanishes), Cu's
  resistivity stays comparatively low all the way down to widths just
  above the W = 2*t_liner = 2 nm floor -- pushing the realistic-edge
  graphene crossover down to W ~ 3 nm, deep in the sub-nanometer-scale
  regime graphene would need aggressive edge quality (not the diffuse,
  p=0.15 assumption) to compete at all. This is a genuinely different
  qualitative story from either the TaN case or the original 2026-08-23
  liner-free baseline: a thin, well-conducting liner materially narrows
  graphene's interconnect advantage window at very small W, which is
  exactly the physical rationale the interconnect-materials literature
  gives for pursuing Ru/Co liner integration in the first place (Section
  1 above) -- this model reproduces that qualitative industry motivation
  from first principles rather than assuming it.
- Both idealized-edge (p = 0.9) graphene curves stay below every Cu
  variant modeled today across the full scan range, unchanged from prior
  sessions' findings.

## 5. Open items after this note

- **Tighten the representative liner resistivity values (Section 2)
  against an actual quantitative source.** This is the most direct
  follow-on: today's session confirmed the right qualitative ordering but
  could not read a specific thickness-resolved number for any of TaN, Ru,
  or Co due to access failures (Section 1). A future session with
  different search/fetch luck (or access to a specific open-access PDF,
  e.g. via an arXiv mirror rather than the publisher page) should retry
  this rather than treat today's placeholder table as final.
- The spatial current distribution *within* the liner frame is treated as
  uniform (a single lumped rho_liner over the whole frame area); a real
  liner's resistivity can vary with local thickness/microstructure near
  corners. Not pursued today -- a second-order refinement relative to the
  liner-resistivity-value uncertainty already flagged above.
- The graphene-all-around-metal (graphene-wrapped Co/Cu core) hybrid
  architecture noted in the 2026-08-23 and 2026-08-28 notes as a distinct
  follow-on item remains unaddressed -- a different device architecture
  from any Cu-liner correction, still a candidate for a future session.

## References

- Domenichini, F. et al., "Selecting alternative metals for advanced
  interconnects," arXiv:2406.09106 (2024); published version
  pubs.aip.org/aip/jap/article/136/17/171101 -- as cited in
  notes/2026-08-28-*.md; re-consulted today for the liner-material
  resistivity-ordering discussion (Section 1), qualitative content only.
- "Mechanisms of Scaling Effect for Emerging Nanoscale Interconnect
  Materials," *Nanomaterials* 12(10), 1760 (2022) -- as cited in
  notes/2026-08-28-*.md, source of the Ru (0.3 nm) and Co (1.0 nm) liner
  thickness figures used today.
- imec public materials on Ru/Co advanced-interconnect metallization
  (industry-motivation context; e.g. imec's "16nm Ru lines using
  semi-damascene integration approach" and "Alternative metals as a path
  to 2nm technology nodes" press pages) -- consulted for qualitative
  liner-material-choice rationale, not a quantitative source.
- Bulk resistivity reference values (Ta, TaN, Co, Ru, Cu) -- standard
  materials-science reference values, not independently re-derived today.
