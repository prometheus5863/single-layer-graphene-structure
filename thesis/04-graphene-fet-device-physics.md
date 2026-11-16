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

The ≈20 GHz is evaluated at V_ds = 0.1 V (peak f_T = 20.279 GHz, peak f_max = 18.731 GHz); see [Appendix A, §4.6 note a](appendix-A-verification-record.md#s4-6-note-a).

**Corrected 2026-10-03.** At the literature-scale 40 µm / 8-finger device the intrinsic peak f_max is 13.094 GHz (was 18.786 GHz), the extrinsic one 2.2196 GHz (was 3.1885 GHz) and f_max/f_T 0.647 (was 0.928); the normalized-device figures above are unchanged. The remaining shortfall is g_ds, not R_g (§4.6.2). Old and new values, cause and detecting check: [Appendix A, §4.6 note b](appendix-A-verification-record.md#s4-6-note-b) and the corrections table, [§A.1](appendix-A-verification-record.md#corrections).

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

### 4.6a Model refinements, and what each did to the numbers above

Sections 4.6.1–4.6.24 extend the RF model and audit the verification layer around it. They are kept under their original numbers in Appendix A. Their effect on this chapter's reported numbers is:

| § | refinement or audit | effect on a reported number |
|---|---|---|
| [4.6.1](appendix-A-verification-record.md#s4-6-1) | numerical provenance of g_ds | none; step/quadrature entanglement bounded |
| [4.6.2](appendix-A-verification-record.md#s4-6-2) | f_max shortfall decomposed | R_g·C_gd carries 0.2154 % of the denominator; f_max/f_T is 0.935276 at R_g = 0; the shortfall is g_ds |
| [4.6.3](appendix-A-verification-record.md#s4-6-3) | velocity-saturation term | two of three acceptance criteria met; the contacts cap g_ds |
| [4.6.4](appendix-A-verification-record.md#s4-6-4) | perfect-contact counterfactual | R_c must be removed from both places it enters |
| [4.6.5](appendix-A-verification-record.md#s4-6-5) | diffusion term | ±0.5 % boundary term; superseded by §4.6.6 |
| [4.6.6](appendix-A-verification-record.md#s4-6-6) | V_ch declared quasi-Fermi | diffusion term is a double count (+0.4279 %); exact charge relation costs −0.6023 % of I_d |
| [4.6.7](appendix-A-verification-record.md#s4-6-7) | signed-carrier model | no Chapter 4 number moved; disorder floor is 2.07× Chapter 3's |
| [4.6.8](appendix-A-verification-record.md#s4-6-8) | Chapter 4 on (4.28) | −0.60 % reproduced; (4.29) is not a linearisation |
| [4.6.9](appendix-A-verification-record.md#s4-6-9) | C_q argument | item 16's direction right, its three numbers wrong by factors of 1.68 to 3.47 |
| [4.6.10](appendix-A-verification-record.md#s4-6-10) | cosh overflow | no committed quantity moves |
| [4.6.11](appendix-A-verification-record.md#s4-6-11) | double regularisation | the item's premise was wrong twice over |
| [4.6.12](appendix-A-verification-record.md#s4-6-12) | item 16 blast radius | §4.7's R_extra moves; Pd's +51.0 becomes −203.276 under the correction |
| [4.6.13](appendix-A-verification-record.md#s4-6-13) – [4.6.24](appendix-A-verification-record.md#s4-6-24) | rebinding, capture and guard audits of the verification layer | no physics number |


## 4.7 Per-metal Rc recalibration: the negative residual was an artefact of the crossover

Section 4.5's contact-doping model computes a work-function-dependent
"extra" resistance, R_extra, contributed by the extended in-plane doping
region beyond a contact's edge. This section tests it against real,
individually cited, per-metal measured Rc values: does the measured
lumped Rc decompose as Rc,measured = R_extra + R_transmission, with
R_transmission an unmodelled interface-injection term concentrated at
the contact itself?

Four literature Rc values were sourced (see
`notes/2026-08-31-per-metal-contact-resistance-literature-and-
recalibration.md` for full citations and search-access notes; Ti and Cr
could not be recalibrated because no single clean per-metal figure was
obtained), all for top (surface) contacts at room temperature,
gate-biased (on-state): Cu 184 Ω·µm, Ni 470 Ω·µm, Au 519 Ω·µm,
Pd 584 Ω·µm.

**The model's scope.** The doping is computed from W_metal − W₀, where
W₀ = 5.4 eV is the n/p crossover work function of Giovannetti *et al.*
(*PRL* 101, 026803 (2008)) and Khomyakov *et al.* (*PRB* 79, 195425
(2009)) — graphene's own 4.5 eV plus the ≈0.9 eV interface step of the
short-range metal–graphene interaction. That relation holds for
**physisorbed** metals only. Of the four metals with a measured Rc, Cu
and Au are physisorbed; **Ni and Pd chemisorb on graphene and are outside
the model**, the same stance §6.11–§6.12 take. Their rows are computed
and shown only for continuity with the superseded table; no claim is
made from them.

Computing R_extra at the on-state bulk channel density used throughout
this chapter (`n_bulk = 2.0e16 m^-2`) and subtracting it from each
metal's measured Rc (`graphene_contact_doping_model.recalibrate_metal_rc()`,
`rc_recalibration.png`):

| Metal | regime | Rc, measured (Ω·µm) | R_extra, computed (Ω·µm) | Naive R_transmission (Ω·µm) |
|---|---|---|---|---|
| Cu | physisorbed | 184.0 | -445.04 | **+629.04** |
| Ni | *outside model (chemisorbed)* | 470.0 | -225.99 | +695.99 |
| Au | physisorbed | 519.0 | -161.14 | **+680.14** |
| Pd | *outside model (chemisorbed)* | 584.0 | -135.78 | +719.78 |

**What the corrected model shows.** Every tabulated metal with a
measured Rc has W < W₀, so each contact *n*-dopes graphene — in the same
sense as the n-type on-state channel, and more strongly than the gate
does. There is no p–n junction at any of these contacts. The
near-contact region therefore conducts *better* than the bulk channel,
R_extra is **negative**, and the naive residual is positive and larger
than the measured Rc itself. For the two metals the model covers, the
implied transmission resistance is +629.04 Ω·µm (Cu) and +680.14 Ω·µm (Au).

**The negative-residual result was an artefact.** Until 2026-11-16 the
code measured the doping from graphene's 4.5 eV rather than from W₀,
contradicting §4.5's own statement of the crossover. That p-doped Cu,
Ni, Au and Pd against the n-type channel, put a spurious p–n junction at
every contact, and inflated R_extra past the measured Rc for three of
the four metals (Cu −787.5, Ni −360.5, Au −90.2 Ω·µm). The conclusion
drawn from it — that the additive decomposition fails, robustly — does
not survive the correction and is withdrawn. The superseded table, its
two annotations and the correction are in
[Appendix A, §4.7 notes](appendix-A-verification-record.md#s4-7-note-a)
and the [corrections table](appendix-A-verification-record.md#corrections).

**What is and is not established.** A positive residual makes the
additive decomposition *sign-consistent* for Cu and Au; it does not
confirm it. For Au, |R_extra| is a minority of the measured Rc, so the
implied R_transmission is set mainly by the measurement. For Cu it is
not: |R_extra| exceeds Cu's entire measured Rc, so the decomposition
needs an interface resistance 3.4× the measured total, offset by
a channel-side reduction larger than that total. That is not
unphysical, but it is strained, and it says the model's near-contact
conductance enhancement is at least as uncertain as the residual it
produces. The positive signs are not fragile: the residuals are
well-conditioned (κ ≤ 0.81; Cu's sign would need a 141.3 % error in R_extra to
flip, `graphene_sensitivity_audit.py`), and with the `C_q`-argument correction
(form B_T) installed they stay positive (Cu 599.06, Au 647.16 Ω·µm,
`graphene_default_flip.py`). Two further caveats apply. First, the work
functions are tabulated polycrystalline values while W₀ comes from DFT
on (111) surfaces; those same calculations find Au(111) *p*-type at
their own W = 5.54 eV, whereas this table's W = 5.10 eV makes it
n-type, so Au's sign is not established here. Second, the TLM
extraction behind all four measured values can absorb part of the
near-contact region this model integrates separately, so R_extra and
Rc,measured may still overlap. Pt (W = 5.65 eV > W₀) is the one tabulated
physisorbed metal that p-dopes graphene and forms a genuine p–n junction
against the n-type channel; it has no measured Rc in this set.

A secondary, more solid finding from the same literature search: within
a single paper and a single metal (Au, Passi et al., arXiv:1807.04772),
switching from a top contact to an edge (hole-patterned) contact reduces
Rc from 519 to 45 Ω·µm at the same gate bias — an ~11x reduction from
geometry alone, a substantially larger lever than the ~3x metal-to-metal
spread observed at fixed (top-contact) geometry in the table above. This
was not modeled quantitatively in this section (the recalibration above
has no notion of contact geometry) — it is picked up in Section 4.8 below,
which finds that only part of this ~11x is attributable to the isolated
doping-density effect (see also the patterned-vs-normal Cu/Pd comparison
already in the Section 4.5 literature review).

## 4.8 Contact geometry: edge vs. top contacts, quantitatively

Section 4.7 flagged, but did not model, that switching from a top to an
edge (hole-patterned) contact reduces Au's measured Rc by ~11x (519 to
45 Ω·µm, Passi et al., arXiv:1807.04772) — a bigger lever than the ~3x
metal-to-metal spread found at fixed top-contact geometry. This section
closes that gap with `graphene_edge_contact_model.py` (full literature
review in `notes/2026-09-05-edge-vs-top-contact-geometry.md`).

**Why edges inject better.** Graphene is sp2-bonded with no out-of-plane
dangling bonds, so a top contact can only couple via a weak van-der-Waals
overlap (Wang et al., *Science* 342, 614 (2013), whose fully edge-contacted
devices — graphene contacted only at its 1D edge inside an hBN
encapsulation stack — reached ~100 Ω·µm, "smaller than what can be
achieved for contacts at the graphene top surface"). At an exposed edge,
carbon atoms can instead form direct sigma bonds with the metal. Passi et
al.'s own DFT calculations quantify this for Au: the metal-induced
Fermi-level shift is **0.35 eV at an edge vs. 0.14 eV at the flat
surface** — a ~2.5x stronger doping-induced shift right at the edge.

**Reusing Section 4.5's machinery.** Rather than inventing a new
edge-specific interface capacitance (nothing in the literature constrains
one), this session converts the two DFT Fermi shifts directly into contact
carrier densities via graphene's linear dispersion,
n(E_F) = E_F²/(π(ħv_F)²), then feeds them into the *same* doping-profile
integration Section 4.5 already uses (factored out this session as
`junction_extra_resistance_from_ncontact()` so both routes — work-function-
derived for top contacts, Fermi-shift-derived for edge contacts — share one
implementation). Result: R_extra = 982 Ω·µm pure top-mode vs. 379 Ω·µm pure
edge-mode for Au — edge-mode is lower, as expected, but by only ~2.6x, not
the full ~11x Passi et al. measured on their actual patterned device. That
gap is expected and explicitly not claimed to be closed here: the
measured 11x bundles in the *specific patterned geometry* (Section 4.8.1
below) on top of the pure doping-density effect isolated here, and R_extra
is already known (Section 4.7) to run higher than measured Rc in this
model's on-state convention, for the same TLM-double-counting reasons
discussed there.

### 4.8.1 Patterned (hole-array) contacts: geometry model

Passi et al.'s actual device etches an array of round holes through the
graphene under the contact metal before deposition, mixing edge-mode and
top-mode injection across the pad area. This session adds a purely
geometric model: for holes of diameter D at areal fill fraction f on a
square lattice, the ligament width between adjacent holes is
D·(√(π/4f) − 1); once that ligament shrinks below 2·λ_decay (the same
250 nm contact-doping decay length used throughout Section 4.5), doping
fronts from neighboring holes overlap and the local graphene is treated as
fully edge-dominated. The edge- and top-mode resistances are then mixed by
this area fraction (conductance-weighted) and a 1/(1−f) current-
constriction penalty is applied for the etched-away area.

Passi et al. do not report the actual fill fraction f used for their five
tested hole diameters (50–1000 nm), so this cannot be fit point-by-point;
it is instead swept at a few illustrative, explicitly-assumed f values and
compared qualitatively (`edge_vs_top_contact.png`). The model reproduces
the *shape* of the large-diameter branch of their data — resistance rising
again as holes grow past a few hundred nm and area loss starts to dominate
over the added edge length — but does **not** reproduce their measured
upturn at small diameters (50–100 nm), where the model instead predicts a
flat, saturated best case. That small-D upturn more plausibly reflects a
fabrication effect (lithographic proximity degradation, or edge-disorder-
limited mobility in very narrow graphene ligaments between closely packed
small holes) outside this compact model's scope, and is reported as an
open item rather than fitted away.

The ~11× edge-contact comparison was later found to have been made against the wrong column; the model is 14% off the right one ([Appendix A, §4.8.2](appendix-A-verification-record.md#s4-8-2)).

## 4.9 Summary and open items

| Sub-topic | Status |
|---|---|
| Quantum capacitance / self-consistent channel charge | Complete (`graphene_fet_model.py`) |
| GFET DC transfer characteristics | Complete (`graphene_fet_model.py`) |
| Contact resistance (lumped, literature-calibrated) | Complete (`graphene_fet_model.py`) |
| Contact-resistance-vs-channel-length crossover | Complete (`contact_resistance_crossover.py`) |
| Spatially-resolved, work-function-dependent contact doping | Complete, diagnostic model (`graphene_contact_doping_model.py`) — see Section 4.5 for known limitations |
| RF figures of merit (f_T, f_max) | Complete, incl. access resistance + extrinsic pad-capacitance estimate (`rf_small_signal_model.py`, Section 4.6); numerical provenance of g_ds audited 2026-09-26 and both entangled defaults exonerated to <10⁻⁶ (Section 4.6.1) |
| Using the doping-profile model to *recalibrate* Rc per metal (vs. using it only diagnostically) | Attempted (Section 4.7). Corrected 2026-11-16: with the doping measured from W₀ = 5.4 eV the naive residual is positive for both physisorbed metals with a measured Rc (Cu +629.04, Au +680.14 Ω·µm), so the additive decomposition is sign-consistent, not refuted; the earlier three-of-four failure was an artefact ([Appendix A](appendix-A-verification-record.md#s4-7-superseded)). Ni and Pd are outside the model; Au's doping sign is not established; Ti and Cr not recalibrated (no literature Rc sourced) |
| Contact-geometry dependence (edge vs. top, patterned contacts) | Complete, DFT-Fermi-shift-derived model (Section 4.8, `graphene_edge_contact_model.py`) — reproduces the qualitative direction (edge lower than top) and the large-hole-diameter branch of Passi et al.'s patterned-contact data; does not reproduce their small-diameter upturn or the full ~11x measured device-level reduction (only ~2.6x from the isolated doping-density effect)**— annotated 2026-10-01 (Section 4.8.2): the ~11x is an ON-STATE figure, and the model is gate-independent. Against Passi et al.'s Dirac-point column — read for the first time this session, having been a cited but dead constant for ten sessions — the same geometry gives 3.01x and the model's 2.59x is within 14% of it. The small-D upturn is confirmed to be a carrier-density effect (3.5x larger in the on state than at the Dirac point), so a purely geometric model cannot produce it** |

Cross-references: Chapter 5 (interconnects) and Chapter 6 (photodetectors)
both trace performance-limiting effects to the same underlying physical
picture developed here — that an atomically thin conductor is unusually
sensitive to boundary/interface quality (edges for interconnects, metal
contacts for both FETs and photodetector charge collection). Chapter 6,
Section 6.5 in particular flags a "spatially resolved collection model" as
a prerequisite for improving the photodetector responsivity estimate;
today's Section 4.5 contact-doping model is a direct, reusable building
block toward that.

## 4.10 Conditioning (moved)

A condition-number audit of §4.7 found that the sensitivity of its rows to input error runs opposite to their plausibility. It is in [Appendix A, §4.10](appendix-A-verification-record.md#s4-10).

## 4.11 Electron–hole mobility asymmetry: the floor is a geometric mean, and the Dirac point is not a doping measurement

Every number in this chapter so far has been computed with a single mobility.
`graphene_fet_model.sheet_conductivity()` takes `|n|` and multiplies by one
`μ`, so the electron and hole branches of every transfer curve in §4.3–§4.5
are identical by construction. That is not an approximation the chapter
declared; it is one the code made silently, and this section is about what it
costs.

### 4.11.1 Why the branches differ, and why the literature does not agree about it

Coulomb scattering of a Dirac fermion from a charged impurity depends, within
the Born approximation, on the *square* of the impurity charge. The sign drops
out, electrons and holes scatter identically, and electron–hole symmetry is
exact at that order. Any measured asymmetry therefore originates outside
leading-order Born, and two mechanisms compete:

| | channel origin | contact origin |
|---|---|---|
| source | Zion *et al.*, arXiv:1801.03565 | Huard *et al.*, arXiv:0804.2040 |
| geometry | Hall bar, probes 30 µm from the metal | probes crossing the strip vs. external arms |
| result | `μ_h < μ_e`, asymmetry shrinking with `\|n\|` | `ρ_e/ρ_h = 1 ± 0.03` away from the metal, ≈ 1.20 through it |
| mechanism | beyond-Born scattering off negatively charged acceptors (≈ 1.5 × 10¹² cm⁻² in the SiO₂) | the metal pins the graphene beneath it p-type, making the junction p–n for an n-type channel and p–p′ for a p-type one |
| scales with | channel length | metal/graphene **interface width** |
| worse branch | holes | electrons |

They disagree about sign as well as mechanism. This thesis does not adjudicate
that — it cannot, from a model — but §4.11.4 shows the two are distinguishable
in principle and says exactly what experiment would distinguish them.

Two of Huard *et al.*'s numbers land on this chapter's own calibration, which
is the reason to take the contact mechanism seriously here rather than treating
it as somebody else's samples: their specific contact resistance
`λ = 135 Ω·µm` sits at the **bottom of the 110–500 Ω·µm band** §4.2 calibrates
`Rc_total` against, and their mobility 4600 cm²/V·s is within 15 % of the
shipped 4000 cm²/V·s.

### 4.11.2 The model, and the closed form it turns on

Mobility is a transport parameter and does not enter the charge relation, so
the electrostatics of §4.3 are kept **untouched** and the asymmetry goes only
into the two-carrier Drude sum, evaluated on the signed species decomposition
of §4.6.7:

```
  n_e, n_h = (S ± n_net)/2,          S = sqrt(n_net² + n_p²)
  σ(n_net) = e(μ_e n_e + μ_h n_h)
           = e(μ̄ S + (Δμ/2) n_net)                                     (4.31)
```

with `μ̄ = (μ_e+μ_h)/2` and `Δμ = μ_e-μ_h`. Equation (4.31) is the useful form
because `S` is **even** in `n_net` while `n_net` is **odd** in the overdrive:
the chapter's existing conductivity is the even part, and the entire effect of
the asymmetry is an odd part linear in the mobility difference. Minimising it:

```
  dσ/dn_net = e(μ̄ n_net/S + Δμ/2) = 0
    ⇒  n_net* = -n_p r/sqrt(1-r²),   r ≡ (μ_e-μ_h)/(μ_e+μ_h)            (4.32)
    ⇒  σ_min  = e n_p μ̄ sqrt(1-r²) = e n_p sqrt(μ_e μ_h)                (4.33)
```

The last equality is an identity, not an expansion, because
`1-r² = 4μ_eμ_h/(μ_e+μ_h)²`. **The minimum sheet conductivity of a disordered
graphene channel is set by the geometric mean of the two branch mobilities.**

(4.32) is the part worth dwelling on: the conductivity minimum is **not at
charge neutrality**. It sits displaced onto the low-mobility branch, at the
density where the slow species' extra carriers exactly compensate the fast
species' lost ones. Evaluating (4.31) at that displaced point is what converts
an arithmetic mean into a geometric one.

### 4.11.3 Correction of record, and the result for a device

**§4.6.7's floor claim is correct where it was evaluated and overstated as
written.** `graphene_fet_signed_carrier_model.sigma_min_ch4()` returns
`e·μ·n_puddle` and calls it *exact*. On the `μ_e = μ_h` diagonal — where the
shipped model lives and where that function evaluates — it is exact, and (4.33)
reproduces it to **1 ulp**, checked rather than asserted. Off the diagonal:

1. `e·μ·n_puddle` **overestimates** the floor by `μ̄/sqrt(μ_eμ_h)`;
2. the accompanying phrase "at `dV = 0`" is no longer *where the minimum is*.

No number in §4.6.7 or §3.3.1 moves, and none is deleted; both are annotated
in place, and §3.3.2 carries the knock-on (the matched puddle density is an
underestimate by +0.88 % at `μ_e/μ_h = 1.3`).

The device-relevant result is the **difference in orders**:

```
  apparent Dirac shift :  ∝ r                        FIRST order
  floor change         :  sqrt(1-r²) = 1 - r²/2      SECOND order
```

Inverting (4.32) through this chapter's own series-factor charge relation at
the shipped 90 nm SiO₂ back gate (`C_ox = 3.8368 × 10⁻⁴ F/m²`):

| `μ_e/μ_h` | `r` | apparent `ΔV_dirac` | fraction of the shipped `V_dirac = 0.80 V` | misattributed density |
|---|---|---|---|---|
| 1.05 | 0.0244 | −0.0524 V | 0.065 | 1.25 × 10¹⁰ cm⁻² |
| 1.10 | 0.0476 | −0.1011 V | 0.126 | 2.42 × 10¹⁰ cm⁻² |
| 1.20 | 0.0909 | −0.1922 V | 0.240 | 4.60 × 10¹⁰ cm⁻² |
| 1.30 | 0.1304 | −0.2763 V | **0.345** | 6.62 × 10¹⁰ cm⁻² |
| 1.50 | 0.2000 | −0.4278 V | 0.535 | 1.02 × 10¹¹ cm⁻² |

A 30 % asymmetry — inside the range both papers above report — moves the
resistance maximum **0.276 V** to lower gate voltage, **34.5 %** of the 0.80 V
that §4.2 attributes entirely to built-in doping. The sign follows from (4.32)
with no sweep: `μ_e > μ_h` puts the minimum on the hole side, so the apparent
Dirac point moves **down** in `V_g`.

**The measurement consequence, stated plainly:** reading a doping level off the
position of a transfer-curve minimum assumes `r = 0`, and nothing in a
two-terminal transfer curve tests that assumption. This is a caveat on an
extraction method, not a defect in a model, and it is the practical reason this
section exists. It also explains a pattern in the literature — that reported
minimum conductivities cluster tightly while reported Dirac-point offsets
scatter widely — without needing to invoke sample quality.

### 4.11.4 Telling the two mechanisms apart, exactly

The odd part of the two-terminal resistance is

```
  R_odd(L) ≡ R(-ΔV) - R(+ΔV)
           = (L/W)(1/σ(-ΔV) - 1/σ(+ΔV)) + (Rc_p - Rc_n)                (4.34)
```

**exactly linear in channel length.** The slope carries the channel asymmetry;
the intercept carries the contact asymmetry. Neither can be absorbed into the
other at any length, because the contact term enters with no factor of `L`.
This is the structure Huard *et al.* report from the measurement side: their
odd resistance collapses when scaled by interface width, and a bulk-only
reading of the same data yields a resistivity asymmetry that **drifts with
probe spacing** — the signature of a nonzero intercept being misread as a
slope.

Because (4.34) is exactly linear, a planted truth must come back out of it, and
that is the strongest check on this section. Planting `μ_e/μ_h = 1.3` and
`ΔRc = 120 Ω` and solving two points of a 25-length sweep recovers the slope to
**2.2 × 10⁻¹⁶ relative** and the intercept to the printed digit, residuals
**1.4 × 10⁻¹² Ω**; all three scenarios (channel only, contact only, both) are
adjudicated correctly from the fitted pair alone.

**The cost, which is this section's honest limit.** At the shipped
`L = 200 nm`, channel-only `R_odd` is 80.93 Ω and contact-only `R_odd` is
120.00 Ω — a factor **1.48** apart. A single-length measurement is blind
between the two mechanisms, and the device this chapter models is a single
length. (4.34) is therefore a statement about an experiment this thesis does
**not** model, offered as a method rather than claimed as a result.

### 4.11.5 What this section does not claim

`μ_e` and `μ_h` are **constants** here. Zion *et al.*'s central observation is
that the asymmetry is density-dependent and dies off as the long-range
impurity potential screens (their hole crossover is at `p ≈ 3 × 10¹² cm⁻²`),
and a constant pair cannot represent that. A density-dependent pair is logged
as an open item rather than guessed at. Likewise the contact asymmetry here is
a constant resistance difference, not Huard *et al.*'s `n^{-1/2}` junction
law; §4.11.4's separation argument needs only that the contact term be
`L`-independent, which the `n^{-1/2}` law also is, so the argument survives the
simplification even though the fitted intercept would then be bias-dependent.

Nothing in §4.3–§4.10 is recomputed. The shipped model sits exactly on the
`r = 0` diagonal, so every published number in this chapter is the `r = 0` case
of the forms above. What this section changes is what those numbers are claims
*about*: the floor is a geometric mean evaluated on a diagonal, and the 0.80 V
`V_dirac` is a doping offset only under an assumption the chapter had never
stated and no committed check had ever tested.

Source: `graphene_branch_asymmetric_mobility.py`, 43/43 checks green. Figure:
`branch_asymmetric_mobility.png`. Machine output:
`branch_asymmetric_mobility_output.txt`. Derivation and literature:
`notes/2026-10-29-the-floor-is-a-geometric-mean.md`.

*Annotated 2026-11-03: the crossover p ≈ 3 × 10¹² cm⁻² quoted above is a hole-branch number, not a material constant, and is high by at least 13 %; see [Appendix A, §4.11.5a](appendix-A-verification-record.md#s4-11-5a).*
