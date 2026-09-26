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

> **Operating point, made explicit 2026-09-26.** The ≈20 GHz above is
> evaluated at **V_ds = 0.1 V** — the value `plot_fT_fmax()` passes, and
> therefore the value behind `rf_figures_of_merit.png`. Measured precisely:
> peak f_T = 20.279 GHz, peak f_max = 18.731 GHz. This sentence previously
> quoted the number without naming the bias, which invited exactly the
> confusion described in Section 4.6.1: `output_conductance`'s own signature
> default is V_ds = 0.05 V, at which the same model gives peak
> f_T = 10.140 GHz and peak f_max = 9.355 GHz. Neither number is wrong and
> neither supersedes the other; they are one model at two bias points, and
> the ≈20 GHz figure quoted throughout this chapter is the V_ds = 0.1 V
> one.

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

### 4.6.1 Numerical provenance of g_ds, and a caution about auditing a discretised model

The output conductance g_ds entering the f_max expression above is not a
closed-form derivative. It is a central difference in V_ds of
`transfer_characteristic()`, which is itself a 50-point quadrature of the
local channel resistance over the drain-bias drop
(`V_channel_profile = linspace(0, V_ds, 50)`). Two numerical defaults are
therefore entangled in every f_max number in this chapter: a step size
ΔV_ds = 10⁻³ V and a resolution n_segments = 50. They are not independent,
because the quadrature grid is itself set by V_ds — differencing in V_ds
differences the quadrature error too.

That entanglement was audited on 2026-09-26
(`graphene_gds_quadrature_audit.py`,
`notes/2026-09-26-gds-step-quadrature-audit.md`) because the analogous
absolute step in Chapter 6 (`H_DIFF = 10⁻³`) had been convicted eight days
earlier of moving eleven of Section 6.13's sensitivities by up to 44%.
**Both defaults here are innocent, and the numbers in this section stand
unchanged:**

| quantity | measured effect on g_ds | effect on peak f_max |
|---|---|---|
| step ΔV_ds = 10⁻³ V vs. anchored-plateau step | 1.25 × 10⁻⁶ | — |
| n_segments = 50 vs. n → ∞ (Richardson in 1/(n−1)) | 8.4 × 10⁻⁷ | — |
| both, propagated | — | 1.1 × 10⁻⁷ (0.00001%) |

f_T contains no g_ds and is untouched; the 400-point V_g grid behind
g_m = ∂I_d/∂V_g moves peak f_T by 0.001% under sixteenfold refinement.

The methodological point survives the null result, and is the reason this
subsection exists rather than a one-line footnote. **A step-refinement
study of a discretised function converges to the derivative of the
discretisation that was held fixed, not to the derivative of the model.**
The anchored step criterion adopted in Section 6.13 — accept a step only if
the derivative is unchanged at h/10 *and* h/100 — is therefore blind to a
quadrature bias by construction, at any tolerance. This was verified rather
than argued: the n = 50 bias measured 3.04 × 10⁻⁷ at ΔV_ds = 10⁻³ and
3.04 × 10⁻⁷ at ΔV_ds = 10⁻⁷, drifting 2.3 × 10⁻⁴ across four decades of
step refinement. And the mechanism is not small in general: in a test case
with R(V_ch) = a + b·V_ch², where both the discretised and the continuum
channel average are exact closed forms (the discrete mean of V_ch² over an
endpoint-inclusive grid is V²(2n−1)/(6(n−1)) against a continuum V²/3, so the
bias is exactly V²/(6(n−1))), the same mechanism reaches **1.32%**. Graphene's
channel resistance simply has too little curvature in V_ch across a 50 mV
drop for it to bite here. That is a property of this model, not a general
licence to ignore the effect — a higher-V_ds or a shorter-channel model with
stronger pinch-off curvature would not inherit this exoneration.

Two defects in the surrounding machinery were found by the same audit and
are recorded here because they bear on how the chapter's numbers should be
read. First, the guard `max(V_ds − ΔV_ds, 10⁻⁴)` in `output_conductance()`
moved the evaluation interval without changing the divisor, returning the
true secant slope times exactly (V_ds + ΔV_ds − 10⁻⁴)/(2ΔV_ds) — an 8.42%
silent under-report whenever it fired. It never fired at either operating
point in use (it is reachable only for V_ds ≤ 1.1 mV), so no number in this
chapter was ever affected; it is fixed, with the unaffected path verified
bitwise unchanged. Second, the 2026-09-25 default-scale census scored
ΔV_ds/V_ds as 2 × 10⁻² by reading the function's *signature* default, whereas
at the call site that actually produces this chapter's figures it is
1 × 10⁻². **A default-scale census must read call sites, not signatures**;
the discrepancy is the ratio between the two, here a factor of two and in
general unbounded.

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

> **Annotation added 2026-09-22 (Section 4.10).** No number in the table
> above is changed or retracted. A condition-number audit
> (`graphene_sensitivity_audit.py`) finds that the four rows' **sensitivity
> to input error runs in the opposite order to their apparent
> plausibility**: Pd's `+51.0` — the one row read below as "a small,
> plausible positive residual", i.e. the only evidence that the additive
> decomposition might survive — flips sign on a **9.6 %** error in
> `R_extra`, while Cu's `−787.5` would need **81.1 %**. The negative results
> are therefore the *robust* rows of this table and Pd is the fragile one.
> Section 4.10 gives the full margins and their consequence: the additive
> decomposition fails robustly, and amplified input error is eliminated as
> the cause of the negatives.

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

## 4.9 Summary and open items

| Sub-topic | Status |
|---|---|
| Quantum capacitance / self-consistent channel charge | Complete (`graphene_fet_model.py`) |
| GFET DC transfer characteristics | Complete (`graphene_fet_model.py`) |
| Contact resistance (lumped, literature-calibrated) | Complete (`graphene_fet_model.py`) |
| Contact-resistance-vs-channel-length crossover | Complete (`contact_resistance_crossover.py`) |
| Spatially-resolved, work-function-dependent contact doping | Complete, diagnostic model (`graphene_contact_doping_model.py`) — see Section 4.5 for known limitations |
| RF figures of merit (f_T, f_max) | Complete, incl. access resistance + extrinsic pad-capacitance estimate (`rf_small_signal_model.py`, Section 4.6); numerical provenance of g_ds audited 2026-09-26 and both entangled defaults exonerated to <10⁻⁶ (Section 4.6.1) |
| Using the doping-profile model to *recalibrate* Rc per metal (vs. using it only diagnostically) | Attempted (Section 4.7) — additive decomposition found NOT to hold for 3 of 4 metals (Cu, Ni, Au); only Pd gives a physically plausible residual. Root cause (TLM double-counting vs. `lambda_decay` mismatch) not yet isolated; Ti and Cr not recalibrated (no literature Rc sourced this session) |
| Contact-geometry dependence (edge vs. top, patterned contacts) | Complete, DFT-Fermi-shift-derived model (Section 4.8, `graphene_edge_contact_model.py`) — reproduces the qualitative direction (edge lower than top) and the large-hole-diameter branch of Passi et al.'s patterned-contact data; does not reproduce their small-diameter upturn or the full ~11x measured device-level reduction (only ~2.6x from the isolated doping-density effect) |

Cross-references: Chapter 5 (interconnects) and Chapter 6 (photodetectors)
both trace performance-limiting effects to the same underlying physical
picture developed here — that an atomically thin conductor is unusually
sensitive to boundary/interface quality (edges for interconnects, metal
contacts for both FETs and photodetector charge collection). Chapter 6,
Section 6.5 in particular flags a "spatially resolved collection model" as
a prerequisite for improving the photodetector responsivity estimate;
today's Section 4.5 contact-doping model is a direct, reusable building
block toward that.

## 4.10 Conditioning: which Chapter 4 numbers inherit their inputs' errors, and by how much

This section applies to Chapter 4 a diagnostic developed in Chapter 6. It
adds no physics and changes no input; it asks only how each of this
chapter's computed numbers would respond if one of its inputs were wrong.

### 4.10.1 Why the question is not rhetorical

Section 6.12.5 measured a **65-fold** difference in how two classes of
quantity respond to the *same* 17.12 % input perturbation: a same-sign
(near-cancelling) contact pair moved 77.94 %, a straddling (reinforcing)
pair 1.21 %. The physics of the two cases is identical. What differs is that
the first is a *difference of two comparable numbers*, so its fractional
error is the inputs' fractional error multiplied by the ratio of input
magnitude to output magnitude.

That is arithmetic, not photodetector physics, and Chapter 4's central
recalibration result is a difference of two comparable numbers:

    R_transmission(metal) = R_c,measured(metal) − R_extra(metal)          (4.21)

### 4.10.2 The diagnostic

For a quantity written as a sum of additive terms `Q = Σ_i T_i`, the signed
logarithmic sensitivity of `Q` to term `T_i` is exactly `S_i = T_i / Q`,
read as *a 1 % change in `T_i` produces an `S_i` % change in `Q`*, and

    Σ_i S_i = 1        exactly                                           (4.22)

by Euler's homogeneous-function theorem. The **condition number** is
`κ(Q) = max_i |S_i|`, and the classification used here is `κ ≤ 1`
reinforcing, `1 < κ < 3` mildly ill-conditioned, `κ ≥ 3` near-cancelling
(Chapter 6's same-sign pairs sit at 4.55).

Equation (4.22) is the implementation's exact validation: it held to
`1.8 × 10⁻¹⁵` over all fifteen decompositions audited across Chapters 4, 5
and 6. Three further exact checks are reported in Section 11 of the
2026-09-22 note, and a fourth — three independent exact zeros — is described
in Chapter 5, Section 5.5.

### 4.10.3 Result: Section 4.7's four metals, conditioned

| metal | residual (Ω·µm) | `S(R_c)` | `S(R_extra)` | `κ` | class |
|---|---|---|---|---|---|
| Pd | +50.9 | +11.462 | −10.462 | **11.46** | near-cancelling |
| Au | −90.2 | −5.756 | +6.756 | **6.76** | near-cancelling |
| Ni | −360.5 | −1.304 | +2.304 | 2.30 | mildly ill-cond. |
| Cu | −787.5 | −0.234 | +1.234 | 1.23 | mildly ill-cond. |

A prediction written into the note before the audit ran — that **at least
three** of the four would be near-cancelling — is **falsified**: only two
are. The hand arithmetic behind it had been done for Pd and Au and
generalised from those two; Ni and Cu were never computed. This is the third
consecutive session in which a claim generalised from a partial enumeration
failed on the full one, and the first in which the failure was in a
*prediction about* the model rather than in the model.

### 4.10.4 Result: κ and sign-robustness run in opposite order

`κ` says how an input error is amplified. The question Section 4.7 actually
needs answered is the inverse and sharper one: **how wrong would `R_extra`
have to be to flip the *sign* of the residual?** For (4.21) that fraction is

    f = (R_extra − R_c) / R_extra = −R_transmission / R_extra             (4.23)

| metal | residual | `κ` | `R_extra` must move by | sign is |
|---|---|---|---|---|
| Pd | +50.9 | 11.46 | **−9.6 %** | **fragile** |
| Au | −90.2 | 6.76 | +14.8 % | intermediate |
| Ni | −360.5 | 2.30 | +43.4 % | intermediate |
| Cu | −787.5 | 1.23 | **+81.1 %** | **robust** |

The two orderings are **exactly reversed**, and that reversal carries the
section's two conclusions.

**First, the additive decomposition fails robustly rather than marginally.**
Section 4.7 reported three failures and one survivor and hedged
accordingly. The survivor is the weakest row in the table: a 9.6 % error in
a quantity computed from a model with an admittedly uncertain `λ_decay`
erases Pd's positive residual entirely. Read with its conditioning, the
table does not contain one plausible case and three implausible ones; it
contains one case that establishes nothing and three that establish
something, the strongest being the one that looked worst.

**Second, amplified input error is eliminated as the cause of the negative
residuals** — one branch of an item open since 2026-08-31. Had the negatives
been an artefact of a large `κ` acting on a mis-specified input, the largest
negatives would sit at the largest `κ`. They sit at the smallest: Cu is both
the largest violation and the best-conditioned row. The cause must therefore
lie in the model's structure, which leaves Section 4.7's own two hypotheses
(TLM double-counting the near-contact region; `λ_decay = 250 nm`
over-integrating) standing and removes a third that had not been named. This
is a candidate eliminated by measurement rather than by argument, and it is
the first narrowing of that open item since it was opened.

### 4.10.5 What this section does not claim

1. It does **not** claim any Chapter 4 number is wrong. `κ` converts an
   input error into an output error and is silent on whether an input error
   exists.
2. The 5 % perturbation used for the illustrative error bars is
   **hypothetical**: none of the four literature `R_c` values is quoted with
   an uncertainty. The resulting bars (Pd 57 %, Au 34 %, Ni 12 %, Cu 6 %)
   inherit that status.
3. A second pre-registered prediction — that some *published* Chapter 4 or 5
   number would carry a >100 % bar under a 5 % input error — is also
   **falsified**. The worst is Pd's 57 %. `κ` is large exactly where the
   output is small, so large `κ` here threatens signs, not orders of
   magnitude; the sign-flip margin of (4.23), not the error bar, is the right
   instrument and it was not the one predicted.
4. Sections 4.5, 4.6 and 4.8 are **not** audited here. Only the Section 4.7
   decomposition and Chapter 5's two families were, and Section 4.8's
   geometry model in particular contains a ratio of two computed resistances
   that has never been conditioned.

Figure: `sensitivity_audit.png`, panel (a). Machine output:
`audit_output.txt`. Pre-registration and outcome:
`notes/2026-09-22-condition-number-audit-of-chapters-4-and-5.md`.
