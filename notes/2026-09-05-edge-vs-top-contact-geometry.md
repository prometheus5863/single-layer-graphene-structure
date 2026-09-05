# Contact geometry: edge contacts vs. top (surface) contacts to graphene

**Date:** 2026-09-05
**Status:** Research notes supporting `graphene_edge_contact_model.py`

## 0. Why this topic, now

Every "not yet covered" list in `AUTOMATION_LOG.md` since the 2026-08-31 run
has flagged the same gap: "Contact-geometry dependence (top vs. edge
contact) is not represented at all in `graphene_contact_doping_model.py`,
despite Section 4.7's own finding that it is a larger lever (~11x) than
metal choice (~3x) at fixed geometry." The "~11x" figure there was drawn
from a single same-paper, same-metal (Au) data point already cited in this
repo (Passi et al., arXiv:1807.04772: 519 vs. 45 Ω·µm). This session goes
back to that paper and one adjacent one to build an actual physical model of
*why* edge geometry helps, rather than only citing the top-line ratio.

## 1. The foundational result: true 1D edge contacts (Wang et al. 2013)

Wang, Meric, Huang, Gao, Gao, Tran, Taniguchi, Watanabe, Campos, Muller,
Guo, Kim, Hone, Shepard & Dean, "One-Dimensional Electrical Contact to a
Two-Dimensional Material," *Science* 342, 614-617 (2013),
[doi:10.1126/science.1244358](https://www.science.org/doi/10.1126/science.1244358).
(Full text was not directly fetchable this session -- science.org returned
HTTP 403 -- corroborated instead via
[a contemporaneous institutional/press summary](https://www.sciencedaily.com/releases/2013/10/131031142734.htm),
which quotes the paper's own reported number.)

Graphene, being sp2-bonded with delocalized pz orbitals but no out-of-plane
dangling bonds, "lacks the ability to make out-of-plane bonds, which makes
electrical contact through its surface difficult" -- this is the physical
reason top (surface) contacts are fundamentally compromised: current must
transfer via a weak van-der-Waals-coupled vertical tunneling path over the
metal-graphene overlap area. Wang et al.'s alternative: encapsulate the
graphene in hexagonal boron nitride, then etch a trench through the whole
stack and deposit metal into it, contacting graphene **only at its exposed
1D edge** -- "a 1D interface between the 2D active layer and 3D metal
electrode." Reported contact resistance: **~100 Ω per micron of contact
width**, stated by the authors as smaller than achievable with top
contacts, and (per standard follow-on literature discussion of this result)
essentially independent of contact length once true edge-only current
injection is achieved -- consistent with all current entering through a
1D line rather than being funneled across a finite-area 2D overlap.
This ~100 Ω·µm figure is of the same order as the *best* literature top
contacts already in this repo's `METAL_LITERATURE_RC` (Cu 184, this is not
a dramatically lower number by itself) -- the real advantage shows up in
the *combined* edge+patterning strategy below, and in contact
*reliability*/process independence rather than only the headline number.

## 2. Quantifying edge vs. top on the *same* device: Passi et al.

Passi, Gahoi, Marin, Cusati, Fortunelli, Iannaccone, Fiori & Lemme, "Ultra
Low Specific Contact Resistivity in Metal-Graphene Junctions via Atomic
Orbital Engineering," arXiv:1807.04772 -- **already partially cited** in
`graphene_contact_doping_model.py`'s `METAL_LITERATURE_RC["Au"]` for its
top-contact value (519 Ω·µm, back-gate biased on-state). This session
re-examined the same paper for its *edge*-contact data, which was not
pulled out in the 2026-08-31 session.

Passi et al.'s edge-contact approach is a *patterned* one: rather than a
clean encapsulated-edge geometry (Wang et al.'s route), they etch an array
of round holes through the graphene *underneath* the metal contact pad
before deposition, so the metal reaches graphene edges at every hole
perimeter in addition to the unpatterned top area outside the holes. Fixed
contact pad area 12 µm x 5 µm; hole diameters 50-1000 nm tested by TLM:

| Hole diameter | Rc at Dirac point (Ω·µm) | Rc at V_BG = -40 V (Ω·µm) |
|---|---|---|
| none (top contact) | 1372 | 519 |
| 50 nm | 620 | 212 |
| 100 nm | 732 | 352 |
| 200 nm | 456 | **45** |
| 500 nm | 1354 | 410 |
| 1000 nm | 1590 | 560 |

Two things stand out, both directly relevant to the model built this
session:

1. **Non-monotonic in hole diameter.** 200 nm holes give the *lowest*
   resistance (an 11.5x reduction vs. unpatterned, on-state); both smaller
   (50, 100 nm) and larger (500, 1000 nm) holes are worse, with 1000 nm
   holes actually *exceeding* the unpatterned baseline. The paper's own
   qualitative explanation: "contact resistivity does not just depend on
   the edge perimeter, but also on the remaining graphene" available for
   current transport, i.e. design should "maximize the perimeter to area
   ratio" -- but the paper gives no closed-form formula for this trade-off,
   only the guidance and the five-point table above.
2. **A DFT-computed microscopic reason edges inject better.** The paper's
   own DFT calculations found the metal-induced Fermi-level shift is
   **0.35 eV at a graphene edge vs. 0.14 eV at the flat surface**, for the
   same metal (Au) -- a factor of ~2.5x stronger doping-induced Fermi-level
   shift right at an edge. Physical origin: at an edge, carbon atoms can
   form sigma bonds directly with the metal ("stronger chemical binding,
   leading to higher transmission of carriers") rather than the purely
   van-der-Waals coupling of a flat-surface contact -- consistent with Wang
   et al.'s "no out-of-plane bonds" framing of *why* top contacts are
   compromised in the first place. This 0.35/0.14 eV pair is the one
   quantitative, metal-specific, first-principles number in either source
   that can be turned into a doping-density input for this repo's existing
   quantum-capacitance-based contact-doping machinery (Section 3 below).

## 3. How this connects to `graphene_contact_doping_model.py`

That module already computes a contact-edge carrier density n_contact from
a metal work function via a self-consistent quantum-capacitance /
interface-capacitance combination (`contact_edge_carrier_density()`), then
integrates the resulting spatial doping profile into an "extra junction
resistance" `R_extra` (Ω·µm) via `junction_extra_resistance()`. That
function assumes a *single* effective interface coupling (`C_interface`,
sub-nm vdW gap) appropriate to a **top** contact -- there was no edge-mode
counterpart.

Passi et al.'s DFT Fermi-shift numbers give a direct, metal-specific way to
add one: rather than re-deriving a new "edge interface capacitance" (which
would require guessing an edge-specific gap/dielectric that isn't reported
anywhere), the two reported Fermi-level shifts (0.35 eV edge, 0.14 eV
surface) can be converted *directly* to carrier densities via graphene's
own linear-dispersion relation,

  n(E_F) = sign(E_F) * E_F^2 / (pi * (hbar * v_F)^2)

(the same E_F <-> n relation implicit in `graphene_fet_model.py`'s quantum
capacitance derivation, using the same v_F = 1.0e6 m/s already used
throughout this repo), instead of going back through the work-function/
interface-capacitance self-consistency loop a second time. This sidesteps
having to invent an edge-specific interface capacitance value while still
reusing every other piece of the existing machinery (the doping_profile()
spatial falloff and the puddle-regularized sheet-resistance integral) via
the newly-factored-out `junction_extra_resistance_from_ncontact()`.
Au is p-type at both edge and surface sites here (Au's work function,
5.10 eV, is on the p-type side of the ~5.4 eV crossover reported by
Khomyakov et al. 2010 and used throughout this repo, so both Fermi shifts
are treated as downward/hole-doping, consistent sign).

## 4. Geometric model for patterned (hole-array) contacts

The 12x5 µm contact pad + circular-hole geometry gives a purely geometric
relation between hole diameter D, areal hole-fill fraction f (fraction of
the pad area etched away), and edge length generated: for holes on a square
lattice of pitch a with f = (pi/4)(D/a)^2, the ligament (remaining
graphene) width between adjacent hole edges is

  ligament(D, f) = a - D = D * (sqrt(pi/(4f)) - 1)

This session's model treats a location in the remaining graphene as
"edge-dominated" once it lies within the *existing* contact-doping decay
length (lambda_decay = 250 nm, Khomyakov et al. 2010, already used
elsewhere in this file) of some hole edge -- i.e. once neighboring holes'
doping fronts start to overlap across the ligament. This gives an
"edge-influenced area fraction"

  p_edge(D, f) = min(1, 2*lambda_decay / ligament(D, f))

with no new free length-scale parameter beyond the one this repo already
uses. The extra junction resistance is then a conductance-weighted mix of
the edge-mode and top-mode results, R_extra_eff = 1 / [(1-p_edge)/R_top +
p_edge/R_edge], with a final series "current-constriction" penalty
1/(1-f) applied on top, representing the reduced conducting cross-section
available once a fraction f of the contact-area graphene has been etched
away (the same area-dilution idiom already used for the liner-aware copper
model in `graphene_interconnect_model.py`, 2026-08-28).

**Important limitation, stated up front:** Passi et al.'s paper does not
report the actual hole areal density/fill fraction f used for each tested
diameter (only D is given in the fetched material) -- so this model cannot
be fit point-by-point against their five measured values. It is instead
evaluated as a swept function of D at a few illustrative, explicitly-labeled
assumed f values, and compared to their data *qualitatively*: does it
reproduce the existence of an optimal intermediate hole diameter, and the
right-hand (large-D, area-loss-dominated) branch of their trend? It is
**not** expected to reproduce the left-hand (small-D) upturn they observed,
since nothing in this geometric+doping picture predicts *very* small,
densely-packed holes to be worse than moderately-sized ones -- that upturn
more plausibly reflects a fabrication/lithography effect (proximity-effect
degradation of pattern fidelity, or edge-disorder-limited mobility in very
narrow remaining graphene ligaments) outside the scope of this compact
model, and is reported as such rather than reverse-engineered away with an
extra fitting parameter.

## 5. Citations

- L. Wang et al., "One-Dimensional Electrical Contact to a Two-Dimensional
  Material," *Science* 342, 614-617 (2013),
  [doi:10.1126/science.1244358](https://www.science.org/doi/10.1126/science.1244358)
  (full text 403'd this session; corroborated via
  [ScienceDaily's summary of the paper](https://www.sciencedaily.com/releases/2013/10/131031142734.htm),
  which directly quotes the ~100 Ω·µm figure and the "no out-of-plane
  bonds" mechanism).
- G. Passi et al., "Ultra Low Specific Contact Resistivity in
  Metal-Graphene Junctions via Atomic Orbital Engineering,"
  [arXiv:1807.04772](https://arxiv.org/pdf/1807.04772) (already cited in
  this repo for its top-contact Au value; this session additionally used
  its edge/hole-pattern TLM table and its DFT Fermi-shift comparison).
- P. Khomyakov et al., *Phys. Rev. B* 82, 115437 (2010), arXiv:0911.2027
  (already cited throughout this repo; reused here only for lambda_decay
  and the n/p crossover work function, no new claims from it this session).

**Search access notes:** WebSearch and WebFetch were both available and
used this session. `science.org` returned HTTP 403 on direct fetch (worked
around via a secondary source, not silently dropped). A Semantic Scholar
fetch of the Wang et al. paper returned an empty page body (no usable
content, also not silently substituted with an invented number). A direct
fetch of the Wiley abstract page for Lee et al. 2022 ("Contact Resistivity
in Edge-Contacted Graphene Field Effect Transistors," a second, independent
edge-contact paper identified via search) returned HTTP 403 and was not
pursued further this session -- it remains a candidate for a future
session's cross-check, listed below.

## 6. Not yet covered / candidates for future sessions

- A second independent edge-contact data point (Lee et al. 2022, Advanced
  Electronic Materials, `doi:10.1002/aelm.202101169`) was identified but
  not fetchable this session (Wiley 403) -- would let the edge-mode
  Fermi-shift/carrier-density estimate be cross-checked against a source
  other than Passi et al.
- The model's edge-mode carrier density is only calibrated for Au (the one
  metal with reported DFT edge-vs-surface Fermi shifts). Extending
  `METAL_WORK_FUNCTIONS`'s other metals to edge mode would need either
  their own DFT data (not found this session) or an assumption that the
  edge/surface Fermi-shift *ratio* (~2.5x for Au) generalizes across
  metals, which is untested here.
- The small-hole-diameter (50-100 nm) degradation in Passi et al.'s data is
  explicitly *not* reproduced by this session's geometric model (Section 4)
  and remains an open modeling gap.
