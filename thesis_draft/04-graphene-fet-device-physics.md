# Chapter 4: Device Physics of the Single-Layer Graphene FET

*Draft section — first Chapter 4 content (previously only referenced from
the Chapter 1 status table). This is a working draft, updated as the
accompanying computational studies and literature review progress; see
`notes/` for the underlying research and the individual scripts listed
below for the computational results referenced here.*

## 4.1 Overview

Chapters 2-3 established graphene's intrinsic electronic structure: a
linear, gapless dispersion near the K/K' points, a density of states that
vanishes at the Dirac point, and the resulting anomalous transport and
optical signatures. None of that intrinsic physics, on its own, determines
how a real graphene field-effect transistor (GFET) behaves — that requires
adding the device-level ingredients that connect the material to a usable
three-terminal structure: a gate dielectric, source/drain metal contacts,
and a finite channel geometry. This chapter develops those ingredients in
the order they were added to this thesis's codebase, and shows how each
one feeds into the next.

## 4.2 Quantum capacitance and the self-consistent channel charge

Because graphene's density of states is finite (rather than the
effectively infinite DOS of a bulk 3D metal gate electrode), the total gate
capacitance seen by an applied gate voltage is not just the geometric
oxide capacitance C_ox — it is the series combination of C_ox and a
gate-voltage-dependent **quantum capacitance** C_q(V_g):

C_total(V_g) = (C_ox^-1 + C_q(V_g)^-1)^-1

`graphene_fet_model.py`'s `quantum_capacitance()` implements a
thermally-broadened form of C_q (finite at the Dirac point rather than
diverging to zero, following the phenomenological treatment of
arXiv:1105.5827; see `notes/2026-08-21-contact-resistance-and-quantum-capacitance.md`,
Section 3). `carrier_density()` then combines this with the oxide term to
give the self-consistent channel sheet density n(V_g), regularized near
the Dirac point by a residual disorder-induced "puddle" density
(`n_puddle`) that reproduces the experimentally observed minimum-
conductivity plateau rather than a true zero. `quantum_capacitance.png`
shows C_q and C_total vs. gate overdrive for the 90 nm SiO2 back-gate
stack used throughout this thesis.

This is not a cosmetic correction: because C_q is comparable to (in fact,
for typical back-gate oxides, considerably larger than) C_ox except very
near the Dirac point, it measurably suppresses the induced charge relative
to a naive C_ox-only estimate, and — as shown in Section 4.5 below — it is
also the dominant charge-limiting element right at a metal contact, not
just under the gate.

## 4.3 GFET transfer characteristics and contact resistance

`transfer_characteristic()` computes Id-Vg at fixed Vds using a long-
channel drift approximation, dividing the channel into segments to capture
how the local charge-neutrality point shifts along the channel once Vds is
non-negligible. This reproduces the well-known **V-shaped, ambipolar**
GFET transfer curve (`gfet_transfer_characteristics.png`) — qualitatively
different from a conventional MOSFET's flat, saturated Id-Vg, because
graphene has no bandgap to shut the channel off on either side of the
Dirac point.

A literature-calibrated series contact resistance (`Rc_total`, per-width
values of ~110-500 Ω·µm; see
`notes/2026-08-21-contact-resistance-and-quantum-capacitance.md`, Section 2,
for representative measured values such as Pd ≈110 Ω·µm at 6 K and Ni
≈470 Ω·µm from a 2024 fabrication-process study) is added in series with
the channel resistance. Because graphene is atomically thin, current must
transfer into the channel through an interface rather than a bulk 3D
contact volume, making this contact term fundamentally different in origin
from a conventional ohmic contact — and, as `contact_resistance_crossover.py`
quantifies (Section 4.4), not a small correction at competitive channel
lengths.

## 4.4 Contact resistance vs. channel length: when does the contact dominate?

`contact_resistance_crossover.py` reuses the unmodified carrier-density and
sheet-conductivity functions from Section 4.2-4.3, sweeping channel length
L (rather than gate voltage) to compute total device resistance
(channel + Rc) and the contact-resistance fraction of that total, plus an
analytic crossover length L_x at which the two contributions are equal
(closed-form, since sheet conductivity is L-independent). Results
(`contact_resistance_crossover.png`): L_x = 115 nm for the best-case
literature contact (Pd, ~110 Ω·µm), 314 nm for a mid-range contact
(300 Ω·µm), and 524 nm for the worst case (500 Ω·µm) — confirming that even
the best literature graphene contacts are not negligible at the 200 nm
channel length used elsewhere in this thesis's GFET model, and that
contact engineering is not a secondary concern for any device competitive
with modern scaled nodes.

## 4.5 Why contact resistance is metal-specific: a spatially-resolved doping model

Sections 4.3-4.4 treat contact resistance as a single lumped number, taken
directly from measurements. That is sufficient for the resistance
bookkeeping above, but it does not explain *why* different contact metals
give different Rc, nor capture a second, distinct effect: the contact
metal also locally **dopes** the graphene beneath and beyond it, through
work-function-driven charge transfer, forming an in-plane p-n, p-p', or
n-n' junction that is separate from (and additional to) the interface
transfer resistance.

`notes/2026-08-26-contact-induced-doping-profile.md` reviews the
underlying physics: work-function mismatch between the metal and
graphene's own ~4.5 eV work function drives charge transfer at the
interface (Giovannetti et al., *Phys. Rev. Lett.* 101, 026803 (2008)), but
the resulting Fermi-level shift is limited by graphene's low density of
states — the same quantum-capacitance-limits-induced-charge physics as
Section 4.2, now applied at the (ungated) contact interface rather than
under the gate. Khomyakov et al. (*Phys. Rev. B* 82, 115437 (2010),
arXiv:0911.2027) further show that because graphene's screening is
anomalously weak, this induced doping is not confined under the metal: the
induced potential decays with distance x from the contact edge as x^-1 for
doped graphene, extending hundreds of nm into the exposed channel. The
same reference reports the n/p doping-type crossover work function at
≈5.4 eV — *above* graphene's own 4.5 eV work function — so most common
contact metals (Ti, Cr, Cu, Ag, Al) n-type dope graphene, while only the
highest-work-function metals (Au, Pt, and variably Pd) p-type dope it.

`graphene_contact_doping_model.py` turns this into a compact quantitative
model, in the same style as Sections 4.2-4.4:

- The contact-edge carrier density is computed with the *same*
  `quantum_capacitance()` series-capacitance approach as
  `carrier_density()`, but with the oxide capacitance replaced by a much
  larger effective **interface capacitance** (sub-nm effective separation
  appropriate for a direct/weakly-bonded metal-graphene contact, vs. the
  90 nm back-gate oxide) — reproducing the literature result that
  graphene's own quantum capacitance, not the interface geometry, is the
  dominant charge-limiting element at the contact.
- The spatial decay away from the contact edge is represented by a
  saturating profile f(x) = 1/(1 + x/λ), reducing to unity at the contact
  edge and falling off as ~λ/x for x ≫ λ, matching the Khomyakov et al.
  asymptotic; λ is set to a representative "few hundred nm" scale.
- The extra sheet resistance contributed specifically by this doping
  profile (beyond the already-tabulated lumped Rc) is obtained by
  integrating the local sheet resistance along the junction region and
  comparing to an equal length of bulk-doped channel, with both terms
  regularized by the same disorder-puddle density used in Section 4.2 so
  the integral stays finite across the p-n crossing.

Results (`contact_doping_profile.png`, computed against a representative
n-type gated bulk channel density of 2×10¹² cm⁻²): Ti (W=4.33 eV) induces
weak n-type contact doping and the smallest extra junction resistance
(~55 Ω·µm); Cu, Au, and Pd — all p-type contacts against the n-type bulk
channel, i.e. genuine p-n junctions — add several hundred Ω·µm; Pt
(W=5.65 eV, the most strongly p-type of the metals considered) induces
such a large contact-edge hole density that most of the junction region is
*better* doped, and hence *more* conductive, than the lightly gated bulk
channel, producing a net **negative** extra-resistance figure in this
model. That result is flagged explicitly at runtime rather than hidden:
it is a genuine feature of the classical drift-conductance picture used
here (a strongly overdoped access region conducts well), but the model
omits depletion-region and carrier-injection physics right at the
charge-neutrality crossing itself, which a full electrostatic/transport
solver would need to capture accurately. This is an explicit limitation to
revisit (see Section 4.6).

This section closes a gap flagged since the first automation run
(2026-08-21): the earlier note "currently only a lumped series resistance,
not a spatially resolved p-n junction model" no longer applies to the
*doping-profile* piece of contact physics, though the lumped Rc values
remain the primary quantity used in the resistance bookkeeping of Sections
4.3-4.4 (this model is diagnostic/explanatory, not yet a replacement
calibration for Rc itself).

## 4.6 RF figures of merit

`rf_small_signal_model.py` layers a standard hybrid-pi small-signal model
on top of the DC model above: transconductance g_m and output conductance
g_ds from numerical derivatives of `transfer_characteristic()`, gate-source
capacitance C_gs reusing the same quantum-capacitance-limited gate
capacitance from Section 4.2, C_gd as a fraction of C_gs, and a distributed
gate resistance estimate. It computes current-gain cutoff frequency
f_T = g_m/(2π C_gs) and maximum oscillation frequency
f_max = f_T / (2√(g_ds(R_g+R_s) + 2π f_T C_gd R_g))
(`rf_figures_of_merit.png`). Peak f_T ≈20 GHz for the 200 nm long-channel
device used throughout this thesis is consistent with the literature range
for non-exotic gate lengths (see
`notes/2026-08-22-rf-figures-of-merit-fT-fmax.md` for benchmarks up to
hundreds of GHz for aggressively scaled record devices).

The earlier version of this model (through 2026-08-22) omitted the
source/drain access resistance R_s from the f_max denominator and had no
way to separate an intrinsic (de-embedded-equivalent) estimate from an
extrinsic (as-measured, pad-parasitic-included) one — and, on that basis,
flagged f_max occasionally exceeding f_T as an unphysical artifact.
Revisiting this on 2026-08-26 against more specific literature (Feijoo
et al., *Sci. Rep.* 6, 35717 (2016)) showed that framing was wrong:
de-embedded graphene FETs with deliberately low, engineered gate
resistance routinely show f_max > f_T (that paper reports f_max/f_T
ratios of 1.3-1.4 at every gate length measured), because f_max ∝ 1/√R_g
while f_T is independent of R_g — the ratio is a design outcome, not an
intrinsic ceiling. The model was corrected accordingly rather than forced
to reproduce a rule that turned out not to hold in general: R_s
(= R_c,total/2, reusing the Section 4.3 contact-resistance calibration)
was added to the f_max denominator, and an explicit extrinsic estimate
(+15 fF GSG pad capacitance per pad, evaluated at a literature-scale
40 µm/8-finger device rather than this repo's normalized 1 µm width,
which would otherwise overstate the pad-capacitance effect by two orders
of magnitude) is now reported alongside the intrinsic one. At the
literature-scale device size, the extrinsic estimate degrades both f_T
and f_max to roughly 15-20% of their intrinsic values — the same
direction, if not the same magnitude, as the raw-vs-de-embedded gap
Feijoo et al. report (~60-70%). See
`notes/2026-08-26-fmax-parasitics-and-fT-fmax-ratio.md` for the full
derivation and literature citations.

## 4.7 Per-metal Rc recalibration: a genuine negative-residual result

Section 4.5's contact-doping model computes a work-function-dependent
"extra" resistance, R_extra, contributed by the extended in-plane doping
junction beyond a contact's edge. Section 4.5 used R_extra only
*diagnostically* (comparing metals to each other and to the model's own
qualitative predictions). This section attempts to go one step further —
recalibrate the model against real, individually-cited, per-metal
measured Rc values from the literature — and reports a genuine negative
result rather than a forced fit.

Four literature Rc values were sourced this session (see
`notes/2026-08-31-per-metal-contact-resistance-literature-and-
recalibration.md` for full citations and search-access notes; Ti and Cr
could not be recalibrated because no single clean per-metal figure was
obtained this session — two ResearchGate fetches were rate-limited),
all for top (surface) contacts at room temperature, gate-biased
(on-state): Cu 184 Ω·µm, Ni 470 Ω·µm, Au 519 Ω·µm, Pd 584 Ω·µm.

The natural hypothesis going in was an additive decomposition,
Rc,measured = R_extra (Section 4.5's doping-junction term) +
R_transmission (an unmodeled interface-tunneling/mode-limited-injection
term, concentrated at the contact itself). Computing R_extra at the same
on-state bulk channel density used throughout this chapter
(`n_bulk = 2.0e16 m^-2`) and subtracting it from each metal's measured
Rc (`graphene_contact_doping_model.recalibrate_metal_rc()`,
`rc_recalibration.png`):

| Metal | Rc, measured (Ω·µm) | R_extra, computed (Ω·µm) | Naive R_transmission (Ω·µm) |
|---|---|---|---|
| Cu | 184.0 | 971.5 | **−787.5** |
| Ni | 470.0 | 830.5 | **−360.5** |
| Au | 519.0 | 609.2 | **−90.2** |
| Pd | 584.0 | 533.1 | +51.0 (~9% of total) |

For three of the four metals, the computed doping-junction contribution
**alone exceeds the entire measured lumped contact resistance**, which
would require a negative "transmission resistance" to balance the
books — unphysical, since a resistance cannot be negative. Only Pd gives
a small, plausible positive residual. Rather than adjusting
`lambda_decay` or the bulk-density reference ad hoc until the four data
points come out positive — which would just be curve-fitting a
parameter to this session's four measurements — this is reported as-is:
**the simple additive decomposition does not hold**, most plausibly
because the transfer length method (TLM) used to extract all four
literature Rc values fits total device resistance vs. contact spacing
back to zero spacing, and its "contact resistance" term can already
partially absorb the same near-contact doping-gradient region this
model integrates separately — i.e. R_extra and Rc,measured likely
double-count part of the same physical region rather than summing as
independent series terms. A second, non-exclusive possibility is that
`lambda_decay` = 250 nm (Section 4.5, from Khomyakov et al. 2010's
doped-graphene asymptotic regime) over-integrates the extra-resistance
contribution for these specific literature devices' actual channel
geometry. Distinguishing between these needs either the source papers'
extracted transfer length L_T (not available from the excerpts obtained
this session) or an explicit simulation of the TLM extraction procedure
on top of the doping profile — both left as open items rather than
resolved here.

A secondary, more solid finding from the same literature search: within
a single paper and a single metal (Au, Passi et al., arXiv:1807.04772),
switching from a top contact to an edge (hole-patterned) contact reduces
Rc from 519 to 45 Ω·µm at the same gate bias — an ~11x reduction from
geometry alone, a substantially larger lever than the ~3x metal-to-metal
spread observed at fixed (top-contact) geometry in the table above. This
is not modeled quantitatively here (the current model has no notion of
contact geometry), but is worth flagging as context for why real
fabrication increasingly favors edge or quasi-edge contacts over simple
top contacts (see also the patterned-vs-normal Cu/Pd comparison already
in the Section 4.5 literature review).

## 4.8 Summary and open items

| Sub-topic | Status |
|---|---|
| Quantum capacitance / self-consistent channel charge | Complete (`graphene_fet_model.py`) |
| GFET DC transfer characteristics | Complete (`graphene_fet_model.py`) |
| Contact resistance (lumped, literature-calibrated) | Complete (`graphene_fet_model.py`) |
| Contact-resistance-vs-channel-length crossover | Complete (`contact_resistance_crossover.py`) |
| Spatially-resolved, work-function-dependent contact doping | Complete, diagnostic model (`graphene_contact_doping_model.py`) — see Section 4.5 for known limitations |
| RF figures of merit (f_T, f_max) | Complete, incl. access resistance + extrinsic pad-capacitance estimate (`rf_small_signal_model.py`, Section 4.6) |
| Using the doping-profile model to *recalibrate* Rc per metal (vs. using it only diagnostically) | Attempted (Section 4.7) — additive decomposition found NOT to hold for 3 of 4 metals (Cu, Ni, Au); only Pd gives a physically plausible residual. Root cause (TLM double-counting vs. `lambda_decay` mismatch) not yet isolated; Ti and Cr not recalibrated (no literature Rc sourced this session) |

Cross-references: Chapter 5 (interconnects) and Chapter 6 (photodetectors)
both trace performance-limiting effects to the same underlying physical
picture developed here — that an atomically thin conductor is unusually
sensitive to boundary/interface quality (edges for interconnects, metal
contacts for both FETs and photodetector charge collection). Chapter 6,
Section 6.5 in particular flags a "spatially resolved collection model" as
a prerequisite for improving the photodetector responsivity estimate;
today's Section 4.5 contact-doping model is a direct, reusable building
block toward that.
