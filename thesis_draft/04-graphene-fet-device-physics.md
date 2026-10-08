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

> **Correction, 2026-10-03: the literature-scale f_max numbers are revised
> downward; the normalized-device numbers above are NOT.** The
> distributed gate resistance `gate_resistance()` took its channel length and
> gate width as *function-parameter defaults* read out of
> `graphene_fet_model` at definition time. `compute_fT_fmax`'s
> literature-scale option rescales the device to the 40 µm / 8-finger
> geometry by *rebinding* `graphene_fet_model.W`, and R_g kept the 1 µm it
> had captured — so the gate resistance entering both f_max denominators at
> that geometry was too small by exactly W_RF/W = 40 (0.208 Ω instead of
> 8.33 Ω). The three numbers that change, with the superseded values kept:
>
> | quantity, 40 µm / 8 fingers, V_ds = 0.1 V | superseded | corrected |
> |---|---|---|
> | intrinsic peak f_max | 18.786 GHz | **13.094 GHz** |
> | extrinsic peak f_max (+15 fF/pad) | 3.1885 GHz | **2.2196 GHz** |
> | f_max / f_T at that geometry | 0.928 | **0.647** |
>
> Peak f_T is **bitwise** unchanged at every bias and both geometries,
> because R_g does not enter f_T = g_m/(2π C_gs) at all; and the
> **normalized 1 µm device** — the path behind the ≈20 GHz and 18.731 GHz
> figures quoted above and behind panels 1 and 2 of
> `rf_figures_of_merit.png` — is bitwise unchanged too, verified against the
> pinned pre-correction source rather than assumed. Only the literature-scale
> comparison moves.
>
> **Two consequences for the argument of this section, in opposite
> directions.** First, the extrinsic *degradation* claim below (f_T and f_max
> falling to roughly 15–20 % of their intrinsic values) survives unchanged:
> measured, it moves from 0.1697 to 0.1695. So does the intrinsic/extrinsic
> peak-f_max ratio, 5.892 → 5.899, a 0.13 % change. The error was
> multiplicative in both terms of every *ratio* this section uses, and
> divided out of all of them — which is exactly why it survived five weeks of
> audits that compare ratios against the literature. Second, and less
> comfortably, the corrected model puts **f_max/f_T = 0.65 rather than 0.93**
> at the literature-scale device, i.e. *further* from the f_max/f_T = 1.3–1.4
> that Feijoo *et al.* report for de-embedded devices, not closer. The
> discussion immediately below reads as though the model and that literature
> had been brought into rough agreement; with R_g correct, they have not
> been. The gap is now attributable rather than hidden: at 8 fingers this
> model's R_g = 8.33 Ω is already in the engineered-low range Feijoo *et al.*
> describe, so the remaining shortfall sits in g_ds and in the R_g·C_gd
> feedback term, not in the gate resistance. ~~Quantifying that split is
> left open and is recorded as such in `AUTOMATION_LOG.md`.~~
>
> **Quantified 2026-10-04 in Section 4.6.2, and the sentence immediately
> above is half wrong.** The R_g·C_gd feedback term carries **0.2154 %** of
> the f_max denominator, so it is not a contributing cause at all; and R_g
> is not the binding constraint either, since f_max/f_T is still 0.935276
> in the limit R_g = 0 exactly. The shortfall is g_ds, and g_ds is large
> because the DC model has no current saturation. See Section 4.6.2.
>
> See `graphene_cross_module_delivery_audit.py` (Sections 3, 4 and 7) and
> `notes/2026-10-03-a-frozen-default-in-another-module.md`.

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

### 4.6.2 The f_max shortfall, decomposed: one of the two named causes carries 0.22 %, and R_g was never the lever

The annotation above left the residual f_max/f_T shortfall attributed to
"g_ds and the R_g·C_gd feedback term". That attribution was written without
measuring it. `graphene_fmax_shortfall_decomposition.py` measures it, and
the measurement is separable exactly rather than perturbatively, because
the f_max denominator is a **sum**:

```
  f_max/f_T = 1 / (2 sqrt(denom)),   denom = g_ds(R_g+R_s) + 2π f_T C_gd R_g
                                             \_____ A ____/   \_____ B ____/
```

At the literature-scale 40 µm / 8-finger geometry and V_ds = 0.1 V:

| term | value | share of denominator |
|---|---|---|
| A = g_ds(R_g + R_s) | 6.033531×10⁻¹ | **99.7846 %** |
| — of which g_ds·R_g | 3.175543×10⁻¹ | 52.5182 % |
| — of which g_ds·R_s | 2.857988×10⁻¹ | 47.2664 % |
| B = 2π f_T C_gd R_g | 1.302534×10⁻³ | **0.2154 %** |

**Result 1: the feedback term is not a cause.** Deleting C_gd *entirely* —
not improving it, deleting it — moves f_max/f_T from 0.643007 to 0.643701,
an improvement of **0.108 %** against a shortfall of a factor of 2.2. Its
share rises with drain bias but only to 1.99 % at V_ds = 1 V. Section 4.6's
annotation named two residual causes and one of them is not one. The
superseded sentence is struck through in place above rather than rewritten,
per this thesis's standing practice.

**Result 2, which contradicts a thread this chapter has been pursuing for
weeks: R_g is not the binding constraint, even in the limit R_g = 0.** The
same module evaluates each counterfactual through the same code path with
one input deleted:

| counterfactual | f_max/f_T | vs. baseline |
|---|---|---|
| model as committed | 0.643007 | — |
| C_gd = 0 | 0.643701 | +0.108 % |
| 16 fingers (R_g/4) | 0.827025 | +28.618 % |
| 64 fingers (R_g/64) | 0.927229 | +44.202 % |
| **R_g = 0 (perfect gate)** | **0.935276** | +45.453 % |
| R_s = 0 (perfect contacts) | 0.885467 | +37.707 % |
| g_ds ÷ 4.8483 | 1.410000 | +119.282 % |

A zero-resistance gate still lands **28.06 % below the bottom of the
1.3–1.4 band**, and 64 fingers is already within 1 % of that limit, so the
multi-finger mechanism is exhausted at about 16 fingers. The 2026-10-03
session found a real 40× delivery bug in R_g and was right to; the belief
that R_g was the lever on f_max did not survive being measured.

**Why term A is so large: the DC model has no current saturation.**
`transfer_characteristic()` computes I_d = V_ds/(R_channel(V_g,V_ds) + R_c),
a bias-dependent *resistor* — no velocity saturation, no pinch-off, no
drain-field cutoff. Its output conductance is therefore the channel
conductance itself. Measured, g_ds·R_total lies in [1.0026, 1.0487] across
V_ds = 0.05–1 V, i.e. g_ds ≈ 1/R_total to between 0.3 % and 4.9 %. Term A
then collapses to a pure resistance ratio and

```
  f_max/f_T  ≈  (1/2) sqrt( R_total / (R_g + R_s) )
```

which reproduces the exact value to 0.36 % (0.645352 vs 0.643007). Reaching
1.41 requires R_total/(R_g+R_s) ≥ 4·1.41² = 7.95; this device has
26.38 Ω / 15.83 Ω = 1.67. **f_max is primarily a measurement of output
resistance, and this model has none to measure.** The intrinsic voltage gain
g_m/g_ds runs 0.00514–0.09648 over the same bias range, one to two orders of
magnitude below a real GFET's. The g_ds that would put this device at
f_max/f_T = 1.41 at its own R_g, R_s and C_gd is 7.859726×10⁻³ S
(0.1965 mS/µm) against the model's 3.810651×10⁻² S (0.9527 mS/µm): **the
missing physics is worth a factor of 4.8483.**

**Result 3, a caution that applies to f_T and therefore to the ≈20 GHz
figure quoted throughout this chapter.** Because g_m = d(V_ds/R)/dV_g is
nearly proportional to V_ds in a resistor model, peak f_T here scales
nearly linearly with drain bias:

| V_ds (V) | 0.05 | 0.1 | 0.3 | 0.6 | 1.0 |
|---|---|---|---|---|---|
| peak f_T (GHz) | 10.140 | 20.279 | 60.779 | 121.168 | 200.440 |
| f_max/f_T | 0.6442 | 0.6430 | 0.6385 | 0.6312 | 0.6218 |

That is ×19.8 over a ×20 bias range. A saturating device's f_T does not do
this. **So the agreement between this model's ≈20 GHz and the literature's
tens of GHz is a fact about the bias point `plot_fT_fmax()` happens to pass,
not independent corroboration** — the same model "agrees" with 200 GHz
record devices if asked at V_ds = 1 V. None of the numbers in Section 4.6
are withdrawn; what is withdrawn is reading their proximity to the
literature as evidence. Note also that f_max/f_T *falls* slightly with
drain bias here, where a real device's rises as it saturates: the sign of
that trend is itself a saturation diagnostic, and this model has it
backwards.

**What this section does not claim.** It does not claim the model is wrong
for the use Chapter 4 puts it to. A resistor-plus-contact model is the right
instrument for the contact-resistance questions of Sections 4.3–4.5, 4.7 and
4.8, which are its substantive results, and those are untouched. What it
claims is narrower and, for an RF reader, more important: **f_max from this
model is not a prediction of f_max.** It is a restatement of the device's
access-resistance-to-channel-resistance ratio. Making it a prediction needs
a saturation term in `transfer_characteristic()` — a drain-field-dependent
velocity, or an explicit saturation velocity v_sat with the usual
I_d = W·q·n·v_sat ceiling — and that is now the top open item for this
chapter in `AUTOMATION_LOG.md`, with a measured target (4.85× in g_ds)
rather than a direction.

Validations: 30/30 and 6 of 6 mutants killed, with control A and control B.
Five checks are tolerance-free, each naming the operation it is exact
under: bitwise reconstruction of the model's own f_max (exact under
summation *order* — re-associating the identical algebra breaks it, which
mutant M1 confirms), multiplication by zero, the 1/N² finger law at powers
of two, and the overlapping-limit identity where this module's
default-geometry path must agree with `rf_small_signal_model` bitwise. One
further result is worth reading off the mutation report: the 2026-10-03
frozen-default bug itself, reinstated as mutant M4, is killed by **exactly
one** check — and that check is the only one in the battery that retains a
*magnitude* rather than a sign. See
`notes/2026-10-04-fmax-was-a-resistance-ratio.md` and
`fmax_shortfall_decomposition.png`.

### 4.6.3 A saturation term, and the third exhausted lever: the contacts now cap g_ds

Section 4.6.2 ended with a target that had a tolerance rather than a
direction: the missing physics was worth **×4.8483** in g_ds, and velocity
saturation was the named route. That term has now been built
(`graphene_velocity_saturation_model.py`, 26/26; mutation report 7 of 7 with
control A and control B) and evaluated against all three of the acceptance
criteria §4.6.2 set. **Two are met. The third is missed by 4.3×, and the
reason it is missed is the result of this section.**

#### The model

Soft-saturation (Caughey–Thomas) local drift velocity at β = 1 — the γ = 1
form of Feijoo *et al.* (2019) — with current continuity fixing I_d along the
channel. The channel-length integral then separates, so the drain current
stays closed-form even though v_sat varies with position:

```
  I_d = μWe·Q / (L + μ·S),     Q = ∫ n dV,     S = ∫ dV/v_sat
```

closed on the contacts by V_ds,ch = V_ds − I_d·R_c. Setting S = 0 recovers
drift-diffusion with no saturation along the *same code path*, which is what
makes the overlapping-limit comparison below a comparison and not a second
implementation. The model is added **alongside** `transfer_characteristic()`,
which is left bitwise unchanged; rewiring the repository onto it is a
deliberate act under the 2026-10-02 rule and is not performed here.

#### v_sat is not fitted, and that is load-bearing

Criterion A names a number and v_sat is the only knob that moves it, so
choosing v_sat to satisfy the criterion would be fitting the model to its own
acceptance test. v_sat is therefore taken from optical-phonon emission,
v_sat(n) = (2/π)Ω/√(πn), with ħΩ = 0.10 eV as published by Feijoo *et al.*
(2019). At this device's density that gives 6.63 × 10⁷ cm/s — **above**
Dorgan *et al.*'s measured 1–3 × 10⁷ cm/s band on SiO₂, i.e. the generous
end of the physics.

Inverting the question makes the gap falsifiable rather than rhetorical:
**criterion A requires v_sat = 6.50 × 10⁶ cm/s, below the measured band.** A
sweep across the available phonon energies (0.059 eV and 0.149 eV, the SiO₂
surface polar modes; 0.100 eV, Feijoo's fit; 0.196 eV, graphene's intrinsic
optical phonon) reaches a g_ds factor of **1.2491 at best**. No phonon energy
in the physical range satisfies criterion A. This also discharges the standing
question of whether the remote-polar-phonon cap on SiO₂ matters here: it is
now a measured 1.07–1.25 in g_ds rather than an assertion.

Dorgan *et al.*'s own best fit is β = 2, which saturates *less* hard at these
fields (24.6–26.0 % higher velocity at the largest field reached), so the
β = 1 choice is generous in that direction too and the shortfall is bounded
from both sides.

#### Criterion B — met

| V_ds (V) | f_max/f_T, resistor | f_max/f_T, saturated |
|---|---|---|
| 0.05 | 0.644461 | 0.663187 |
| 0.10 | 0.642966 | 0.683186 |
| 0.20 | 0.640511 | 0.728061 |
| 0.50 | 0.633440 | 0.911727 |
| 1.00 | 0.621775 | 1.384374 |

The slope reverses from −0.0227 to **+0.7212**, 31.8× the magnitude of the
wrong-signed behaviour it replaces. A saturating device's f_max/f_T rises
with drain bias; this one now does, which was the part of §4.6.2's diagnosis
that made a claim about mechanism.

**What this table is not.** The V_ds = 1 V row lands inside Feijoo *et al.*'s
1.3–1.4 band, and that is **explicitly not offered as corroboration**. §4.6.2
recorded why: peak f_T here scales ×19.8 over a ×20 drain-bias range, and a
model with a knob that size will meet a comparably wide literature band
somewhere. The claim is the sign, which is bias-independent; where the ladder
crosses 1.3 is a fact about where the ladder was stopped.

#### Criterion A — missed, and the miss is the finding

| | resistor | saturated | factor |
|---|---|---|---|
| g_ds at the resistor's peak-f_T bias | 3.811e-2 S | 3.355e-2 S | **1.1361** |
| g_ds at the saturated model's peak | 3.767e-2 S | 3.354e-2 S | 1.1230 |
| **required** | | | **4.8483** |

23.4 % of the requirement. The saturation term is not weak: μS/L = 0.3036 at
this bias, which alone would divide the **channel** conductance by 1.699. It
is diluted because saturation acts only on the channel while the contacts sit
in series with it. At the 40 µm / 8-finger geometry R_c,total = 15.0 Ω of a
≈ 26 Ω device, so:

> **Roughly half of g_ds is a contact resistance that no saturation mechanism
> can reach. The contacts, not the channel, now cap g_ds.**

This is the third candidate lever on the f_max shortfall to come back
non-binding. §4.6.2 exhausted the gate resistance (R_g = 0 *exactly* still
falls 28 % short) and the feedback capacitance (C_gd deleted entirely moves
the answer 0.11 %). Velocity saturation now joins them. **All three point at
the access resistance** — which is §4.5, §4.7 and §4.8 of this chapter, the
part of this work with the most content and the one place it carries a
negative result worth having. The RF thread does not end in a tuning problem;
it ends in the contact.

#### Criterion C — met, and it corrects this chapter's own model

The overlapping limit cannot be checked as bitwise agreement, because the
saturated model's S → 0 limit takes the *arithmetic* average of n along the
channel while `transfer_characteristic()` takes a *harmonic* one. Making four
artefacts converge as V_ds → 0 separates three mechanisms, which are
separable precisely because their orders in V_ds differ:

| mechanism | size at V_ds = 0.1 V | order in V_ds |
|---|---|---|
| quadrature rule (50-sample mean vs. trapezoid) | −4.3 × 10⁻⁵ % | 2.04 |
| **profile domain** | **+0.253555 %** | **1.00** |
| Jensen gap (⟨1/n⟩ vs. 1/⟨n⟩) | +3.8 × 10⁻⁴ % | 2.00 |

Two things follow, and both were predicted wrongly before being measured.

First, the `n_segments = 50` discretisation that §4.6.1 flagged and that
§4.6.2 called load-bearing — because g_ds is a V_ds difference of exactly that
quadrature — is **measured at 4.3 × 10⁻⁵ % and is not a problem at the
committed biases.** A long-standing caution closes with a negative result.

Second, the dominant discrepancy is one nobody here was looking for.
`transfer_characteristic()` profiles the local channel resistance over the
**full V_ds**, but over half of V_ds is dropped across the contacts and never
appears across the channel at all, so the model evaluates n(V_ch) over a
potential range about twice too wide. The error is first order in V_ds and
674× the Jensen gap the check was written to find.

**No number in this chapter is withdrawn.** The total is one-signed,
0.2539 % at V_ds = 0.1 V, far below anything §4.6 concludes from, and in the
direction that makes I_d *larger* — so every shortfall against literature
reported in this chapter is, if anything, understated rather than overstated.
What has changed is that the bound is measured, its mechanisms are separated
and ordered, and the larger one was not the one that had been suspected.

#### One unexplained observation, recorded rather than resolved

The peak-f_T bias **moves branch** between the two models: V_g = −1.1128 V
(hole side) for the resistor model, +2.8571 V (electron side) for the
saturated one. Both readings are reported above rather than one being chosen.
A figure of merit whose optimum jumps 4 V between two models of the same
device is either a real electron–hole asymmetry — which this chapter has
independent reason to expect — or an artefact of where g_m peaks, and this
section does not distinguish them.

See `notes/2026-10-05-saturation-arrived-and-the-contacts-ate-it.md`,
`velocity_saturation_model.png` and `velocity_saturation_output.txt`.

### 4.6.4 The perfect contact: the fourth exhausted lever, and a factor of four the ratio hid

§4.6.3 ended by naming the contacts as what caps `g_ds` and left one
counterfactual unrun: `f_max/f_T` at `R_c = 0` *exactly*, in the saturated
model, in the way §4.6.2 ran `R_g = 0` and `R_s = 0`. This section runs it.

`R_c` has to be removed from **both** places it enters, and that is not a
formality. It is a series resistance inside the self-consistent solve for the
intrinsic channel drop, *and* it is `R_s = R_c/2` in the `f_max` denominator,
reached through a different module that reads `gfet.Rc_total` as a global.
Zeroing only the first — which is what the natural keyword-argument route
does — gives `f_max/f_T` = 0.552349, a **19 % degradation**. Zeroing both
gives 0.759363, an **11 % improvement**. The two answers have opposite signs,
and the misleading one is the plausible one, since a 19 % degradation is what
"the contacts are binding" predicts if read without care. Check C3 of
`graphene_perfect_contact_counterfactual.py` asserts the sign inversion so
that it is a recorded measurement rather than a caveat.

**The result, at the 40 µm / 8-finger geometry and `V_ds` = 0.1 V, at a fixed
gate bias:**

| quantity | baseline | `R_c = 0` coherent | change |
|---|---|---|---|
| `g_ds` [S] | 3.374491e-02 | 5.145076e-02 | ×1.5248 |
| `f_T` [Hz] | 2.07398e+10 | 7.46504e+10 | ×3.5994 |
| `f_max` [Hz] | 1.41692e+10 | 5.66868e+10 | **×4.0007** |
| `f_max/f_T` | 0.683186 | 0.759363 | +11.15 % |

Note that the comparison has to be at a **fixed** bias. With `R_c = 0` the
peak of `f_T(V_g)` leaves the sweep entirely and sits at the grid edge in both
models, because the contact is what produced the turn-over: `I_d` is capped by
`R_c` at large `|V_g − V_dirac|` while `C_gs` keeps growing, so
`f_T = g_m/2πC_gs` turns over. Remove the cap and it does not. Widening the
sweep moves the edge, not the physics. "At the peak-`f_T` bias" is therefore
not an available convention for this counterfactual, and check P1 asserts the
interior/boundary distinction rather than leaving it as a remark.

**The fourth lever is exhausted on the ratio.** 0.759363 is 58.4 % of the
bottom of Feijoo *et al.*'s 1.3–1.4 band. A residual `g_ds` requirement of
**×2.9955** survives the perfect contact, and §4.6.3 bounded what `v_sat` can
supply at ×1.2491 across the whole physically available phonon range. So all
four quantities that have been proposed as the `f_max` shortfall's cause —
`R_g`, `C_gd`, `v_sat`, `R_c` — have now been set to their most favourable
value and none closes the gap.

**And the ratio hid a factor of four.** The same counterfactual is worth
×4.0007 in `f_max` itself. Section 6b of the module sweeps contact resistances
that have been measured rather than assumed:

| `R_c` per contact [Ω·µm] | `f_max` [GHz] | `f_max/f_T` | share of the `R_c = 0` gain |
|---|---|---|---|
| 0 | 56.687 | 0.759363 | 100.0 % |
| 65 (lowest reported) | 39.526 | 0.722202 | 59.6 % |
| 165 (Feijoo *et al.*) | 24.308 | 0.695514 | 23.8 % |
| **300 (this model)** | **14.169** | **0.683186** | 0.0 % |
| 470 (practical Ni flow) | 8.268 | 0.680047 | −13.9 % |
| 4000 (typical literature) | 0.222 | 0.698615 | −32.8 % |

Two readings, and the second was not anticipated.

First, **this model's contact is 1.8× worse than that of the devices its
`f_max/f_T` is compared against**, which was not known before this section and
which closes an obvious objection: give this model Feijoo's own 165 Ω·µm and
it reaches 0.695514, not 1.3. The shortfall is not a contact-quality artefact.

Second, **`f_max/f_T` is non-monotonic in `R_c`.** It bottoms out at 0.680047
at 470 Ω·µm and rises again to 0.698615 at 4000 Ω·µm, while `f_max` falls
**37×** across those same two rows. This is §4.6.2's resistance-ratio identity
reasserting itself inside the *saturated* model: once `R_c` dominates
`R_total` faster than it dominates `R_g + R_s`, the ratio is rewarded for a
worse contact. There is a region of this design space where the ratio and the
figure of merit it is supposed to summarise **move in opposite directions**,
and no `f_max/f_T` criterion can be trusted inside it. §7.8.1c draws the
methodological consequence.

**One prediction of §4.6.3 is confirmed, from the other side.** That section
attributed the ×1.1361 `g_ds` factor to dilution by the contacts, computing
from a channel integral that the undiluted factor would be 1.699. That is a
testable prediction at `R_c = 0`, and the measured factor there is
**1.70661** — 100.45 % of it — obtained from a finite difference of a solved
terminal current rather than from the integral that produced the 1.699
(check D1).

**One claim is narrowed, by literature rather than by computation.** Feijoo,
Pasadas, Bonmann *et al.*, *Nanoscale Advances* **2** (2020), report that the
largest `f_max` in their measured GFETs occur "far from the saturated velocity
regime", at ~45 % of `v_sat`, with the **diffusion** contribution to the
current comparable to drift. Equation (4) of §4.6.3 is drift-only. §4.6.3's
criterion A therefore remains correct arithmetic about *this model* and is not
a statement about real devices, and the absence of a diffusion term — not the
absence of saturation — is now the RF thread's binding question. See §7.8.1c
and §7.9 item 7.

See `graphene_perfect_contact_counterfactual.py`,
`perfect_contact_counterfactual.png`,
`perfect_contact_counterfactual_output.txt` and
`notes/2026-10-06-the-perfect-contact-and-the-ratio-that-divided-out-the-prize.md`.

### 4.6.5 The diffusion term: a boundary term worth ±0.5 %, and a potential this chapter never named

§4.6.4 closed by naming the absence of a diffusion term in Eq. (4) as the RF
thread's binding question, on the strength of Feijoo *et al.*'s measurement
that diffusion is comparable to drift at the peak-`f_max` bias. The term has
now been derived and built (`graphene_diffusion_current_model.py`,
**30 checks, 30 passed**), and the result divides into three parts that do not
agree with each other.

**The derivation, and the `beta = 1` separability it preserves.** Zebrev's
generalised Einstein relation ([arXiv:1102.2348](https://arxiv.org/pdf/1102.2348),
Eq. 20) is `mu = e D/eps_D` with `eps_D = n/(dn/dmu_c)`, and graphene's linear
dispersion gives `eps_D = E_F/2` exactly in the degenerate limit, so
`D = mu E_F/(2e)`. Writing `V_F(n) = E_F/e = A_F sqrt(n)` with
`A_F = hbar v_F sqrt(pi)/e` — this thesis's own Dirac dispersion, with no new
parameter — the diffusion term collapses exactly:

```
    D dn/dx = mu n dV_F/dx
    I_d = W e mu n d(V - lambda V_F)/dx                                  (10)
```

so the drift-diffusion current is transport down the gradient of the
quasi-Fermi potential `Phi = V - V_F`, and the separability that makes Eq. (4)
closed-form survives intact:

```
    I_d = mu W e Q_D / (L + mu S_D)
    Q_D = Q + lambda (A_F/3) (n_s^{3/2} - n_d^{3/2})                     (11)
    S_D = S + int kappa(V)/v_sat(V) dV,     kappa = -lambda dV_F/dV      (12)
```

`lambda = 0` recovers Eq. (4) **bitwise** (check X1, 0.0 ULP at twelve
(`V_g`, `V_ds`) points), which is the overlapping-limit requirement of
§4.6.3's criterion C in its strongest available form.

**Eq. (11) is a boundary term, and that matters more than its size.**
`int n dV_F = A_F int n d(sqrt n) = (A_F/3)[n^{3/2}]` is an exact
antiderivative, so the diffusion contribution to the channel integral depends
only on the **source and drain densities** and not at all on the density
profile between them. §4.6.3's criterion C found that this chapter's
profile-domain handling carried a one-signed error 674× larger than the Jensen
gap it was written to look for; the diffusion term **cannot inherit that
error**. Check X6a asserts this structurally rather than arguing it — a
deliberately lopsided interior quadrature grid with the same endpoints leaves
the term bitwise unchanged — and X6c is the contrast, the quadrature route to
the same number moving by 0.4 % on the same regridding.

**The size, at this chapter's own biases, `V_ds = 0.1 V`:**

| bias | `I_d` (λ = 0) [A] | `I_d` (λ = 1) [A] | share | `Q_d/Q` |
|---|---|---|---|---|
| saturated-peak, `V_g = +2.857143 V` | 9.007791e-05 | 9.052996e-05 | **+0.5018 %** | +1.1621 % |
| resistor-peak, `V_g = −1.112782 V` | 8.910103e-05 | 8.865209e-05 | **−0.5038 %** | −1.1546 % |

**The correction is signed and its sign flips between branches.** It follows
`d|n|/dV`, and `carrier_density()` returns a magnitude
(`sqrt(n_eff² + n_puddle²)`), so on the hole branch `|n|` rises toward the
drain and the term subtracts. A model that has discarded the sign of the
carrier cannot settle whether that is the correct physical sign for holes;
this is recorded as an open item (§7.9 item 9) and not as a result. The share
of `I_d` is about half of `Q_d/Q` because `S_d > 0`: Eq. (12) **opposes**
Eq. (11), with `S_d/S = +1.1621 %`.

**The RF consequence.** At the literature 40 µm / 8-finger geometry:

| quantity | λ = 0 | λ = 1 | change |
|---|---|---|---|
| `f_T` | 2.073963e+10 | 2.085305e+10 | +0.5469 % |
| `f_max` | 1.416902e+10 | 1.421070e+10 | +0.2941 % |
| `f_max/f_T` | 0.683186 | 0.681468 | **−0.2514 %** |
| `g_ds` | 3.374491e-02 | 3.391519e-02 | +0.5046 % |

`f_max` rises while the ratio **falls** — §4.6.4's ratio finding in a third
mechanism, and this time the two move in opposite directions rather than
merely by different amounts. §7.8.1a's verdict is unaffected: 0.681468 is
52.42 % of 1.3 either way.

#### A prediction written down first, and refuted — with its mechanism measured

Zebrev's closed form for the local ratio, `kappa = C_ox/(C_Q + C_it)`
(his Eq. 62), **diverges** as `C_Q → 0` at charge neutrality, and Feijoo
*et al.* put the peak-`f_max` bias near the onset of bipolar conduction. Both
predict the diffusion share **peaks at the Dirac point**, and that prediction
was written into the module before its sweep was run.

Measured: the share **collapses** there — `−0.006925 %` at the nearest grid
point to the Dirac point against `0.5053 %` at the maximum, a factor of 73,
with the maximum two volts away at `V_g = −1.272 V`. The mechanism is this
chapter's own regularisation, measured rather than argued:
`n = sqrt(n_eff² + n_puddle²)` makes `d|n|/dV = (n_eff/n) dn_eff/dV`
**exactly 0.0** at `n_eff = 0`, while `C_Q` stays finite at `n_puddle`
(G4b: `0.000000e+00` against `−2.290589e+15 m⁻²/V` unfloored). **The puddle
floor flattens precisely the region where the capacitor ratio diverges.** The
diffusion-dominated regime is not small in this model — it is absent, and that
is a statement about the regularisation rather than about graphene. §4.6.4 was
right that the peak-`f_max` bias cannot be asked about; the reason is the
puddle floor, not the missing term, and it does not go away by adding terms.

#### Annotation (2026-10-07): §4.3's `C_q` is evaluated at the gate overdrive, not at `E_F/e`

`quantum_capacitance()` takes `(V_g − V_ch − V_dirac)` in the slot the
dispersion wants `E_F/e = V_F`. At the RF bias those are **2.057 V and
0.0977 V**, a ratio of **21.05**, so this chapter's `C_q` overstates the
dispersion-consistent `C_Q = (2e³/πℏ²v_F²)V_F` by that factor and would
understate `kappa` by it. §4.6.5's comparison with Zebrev therefore uses the
dispersion-consistent form. **No number in this chapter is withdrawn**: `C_q`
is used downstream only through the series factor `C_q/(C_q + C_ox)`, where
`C_q ≫ C_ox` makes the factor close to unity under either evaluation. The
figures and tables of §4.3 are annotated here rather than recomputed, per the
2026-10-02 rule that a changed rule makes derived numbers stale until
deliberately re-derived.

#### What this section does not claim

It does not claim that Eq. (4) was missing this term. §7.8.1d sets out why
that question is not decided by the measurement above: the standard GFET
compact model this chapter's Eq. (4) descends from writes its drift expression
in the **quasi-Fermi** potential, under which Eq. (4) is already a complete
drift-diffusion current and `lambda = 1` double-counts. `lambda = 1` is
therefore **not** wired into any other module of this thesis, and the
production path remains `lambda = 0` with ±0.50 % as a stated bound.

See `graphene_diffusion_current_model.py`, `graphene_diffusion_mutation.py`,
`diffusion_current_model.png`, `diffusion_current_output.txt`,
`diffusion_mutation_output.txt` and
`notes/2026-10-07-the-diffusion-term-was-a-boundary-term-and-the-potential-was-never-named.md`.

### 4.6.6 The potential is named — and the naming question was standing in front of a larger one

§4.6.5 ended with a variable that had two readings and a term that was exactly
the difference between them. This section names the variable, reports what that
does to §4.6.5's result, and reports the thing the naming question had been
hiding.

#### The declaration

> **`V_ch` is the quasi-Fermi (electrochemical) potential of the channel
> carriers, in volts, measured from the source.**

It is declared in `graphene_fet_model.py` itself, as a machine-readable
in-source marker above `carrier_density()`, so that the repository's statement
of what its own variable means travels with the variable rather than living in
a note. Three grounds, none of them new work:

1. The compact model Eq. (4) descends from defines exactly that variable.
   Pasadas and Jiménez (*IEEE TED* 63(7) 2016,
   [arXiv:1605.08235](https://arxiv.org/pdf/1605.08235)) write `v = μF` with
   `F = −dV/dx` and state that *"V(x) is the quasi-Fermi level along the
   graphene channel"*.
2. Classical long-channel theory uses the same variable and states what it
   buys: the channel variable `V(y)` is the quasi-Fermi potential, and the
   current-density equation it feeds is *"both drift and diffusion"*
   ([long-channel MOSFET notes](https://dunham.ece.uw.edu/ee531/Long_Channel_MOSFET.pdf)).
   In the quasi-Fermi variable, **one term is the complete current**.
3. `carrier_density()` already assumes it. The quantum-capacitance series
   factor `C_q/(C_q + C_ox)` is the correction one applies when the channel
   variable is the quasi-Fermi level and the graphene drop `E_F/e` is carried
   separately. Under a strictly electrostatic reading, that factor is itself a
   double count.

**What the literature does not say, recorded because it bears on the thesis.**
A review of GFET compact models (Lu *et al.*,
[arXiv:1703.09759](https://arxiv.org/pdf/1703.09759)) does *not* make this
identification: it calls the imref splitting "a.k.a. channel voltage" in one
section and describes `V(x)` only as "the voltage along the channel" in the
compact-model sections, never connecting them; its Eq. (12) carries a diffusion
term that its Eq. (15) writes with the same expression as drift; and its
compact model uses a drift form while stating only that drift-diffusion "is
assumed". **The ambiguity §4.6.5 found is inherited from the compact-model
literature rather than introduced here** — which does not excuse it in a
thesis, and does raise the value of declaring it explicitly.

#### §4.6.5's term is a double count, and its numbers stand

`graphene_diffusion_current_model.py` writes the current as `d(V − λV_F)/dx`
and names `V − V_F` the quasi-Fermi potential. Under the declaration `V` *is*
that potential, so `λ = 0` is already the complete drift-diffusion current.

**The term built in §4.6.5 is therefore a double count, not missing physics.**
Measured three ways at `V_g` = 2.0 V, `V_ds` = 0.1 V:

| Level | `I_d` (λ = 0) | `I_d` (λ = 1) | difference |
|---|---|---|---|
| intrinsic (fixed channel drop) | 1.432116×10⁻⁴ A | 1.442091×10⁻⁴ A | **+0.6965 %** |
| terminal (contact feedback on) | 8.190217×10⁻⁵ A | 8.225263×10⁻⁵ A | **+0.4279 %** |
| §4.6.5, at its saturated-peak bias | — | — | +0.5018 % |

The terminal figure must be the smaller of the first two, because the series
contact resistance is a negative feedback on `I_d`; that inequality is checked
rather than assumed, and its sign is fixed by the circuit, not by a fit.

**No number in §4.6.5 is deleted or changed.** The derivation is correct, the
boundary-term structure is correct, and the measured size is correct. What
changes is the direction: the production path is `λ = 0`, and `λ` is retained
as the instrument that *measures* the difference between the two readings. The
±0.50 % ambiguity §4.6.5 attached to every `I_d` in this chapter is therefore
**resolved rather than bounded**.

#### And the larger question the naming had been standing in front of

Declaring `V_ch` quasi-Fermi is not free: it promises that the charge relation
is the one that variable belongs to. Under the declaration the gate drive
divides exactly,

```
V_g − V_dirac − V_ch  =  e·n/C_ox  +  E_F(n)/e ,    E_F(n)/e = A_F·√n    (4.28)
```

which is a quadratic in `√n` with a closed-form root and no new parameter —
`A_F` is this thesis's own Dirac dispersion. §4.3's `carrier_density()` instead
uses the linearised series factor

```
n = (C_ox·dV/e) · C_q/(C_q + C_ox)                                      (4.29)
```

(4.29) is the linearisation of (4.28). Running Eq. (4) on (4.28) with nothing
else changed:

| | `I_d` at the RF bias |
|---|---|
| on (4.29), λ = 0 — **the committed model** | 8.190217×10⁻⁵ A |
| on (4.28), λ = 0 — the declaration taken exactly | 8.140885×10⁻⁵ A |
| **cost of the declaration** | **−0.6023 %** |
| the double count it removes | +0.4279 % |

**The charge-model question is 1.4 times larger than the transport term this
chapter spent 2026-10-06 and 2026-10-07 on.** That is the finding of this
section, and it is not a flattering one: two sessions went into a term worth
+0.43 % while a −0.60 % question about the chapter's own charge relation sat
unexamined directly underneath it, reachable in one line of code.

**A number this section declines to quote.** The raw (4.28)-vs-(4.29) charge
disagreement reaches **+68 %** at `dV` = 0.01 V. Quoting that would be wrong:
`n_puddle` = 5×10¹⁵ m⁻² swamps both forms there and no result in this chapter
lives in that regime. In what `carrier_density()` actually returns, the worst
disagreement across the drives this chapter sweeps is **+2.22 %**; in `I_d` it
is **−0.60 %**. The raw figure would have overstated the cost by a factor of
31. This is §4.10's lesson and 2026-10-06's rule applied to this section's own
result.

#### Annotation to §4.3 and §4.6.5, in place

§4.3's `carrier_density()` docstring is annotated in the source rather than
rewritten: the wording "electrostatic charge from `V_g − V_ch`" is retained and
marked loose, since under the declaration the exact relation is (4.28) and the
series factor is its linearisation rather than an additional physical effect.
**The numbers that function returns are unchanged — verified bitwise over a
41-point sweep — and what changed is what they are a model of.**

#### A numerical defect found on the way, measured and left open

`quantum_capacitance()` returns a finite value up to `dV` = **18.3493 V** and
`inf`/`NaN` above it, because `log(2(1 + cosh η))` overflows once
`η = dV/(kT/e)` exceeds about 710 (`kT/e` = 25.85 mV at 300 K). The exact
large-`η` limit is `|η| + log 2`. The onset is far outside the ~2.7 V maximum
drive this chapter sweeps, so **no committed number is affected**; it is
recorded as §7.9 item 13 rather than patched here, because that module owns
committed transcripts.

#### What this section does not claim

1. It does **not** claim (4.29) is wrong. (4.28) and (4.29) have been compared
   only to each other, never to measurement. Under the declaration (4.29) is
   non-exact; whether it is less accurate is not established here.
2. It does **not** re-derive this chapter. Every committed number remains a
   (4.29) number. The −0.60 % is a measurement of the gap, obtained by
   replacing one function and re-running, not a new result set.
3. It does **not** touch the hole branch. Both (4.28) and (4.29) carry `n` as a
   magnitude, which §4.6.5 already identified as leaving every branch-asymmetric
   result in this chapter unverified. The declaration does not help, and the
   signed-carrier item (§7.9 item 9) now **blocks** the (4.28) re-derivation
   rather than sitting beside it: re-deriving on a magnitude would rebuild the
   same flaw in a new equation.

Machine output: `potential_declaration_output.txt` (10 checks, 10 passed) and
`potential_census_output.txt` (16 checks, 15 passed, 1 expected failure).
Modules: `graphene_potential_declaration.py`,
`graphene_potential_census_audit.py`. Note:
`notes/2026-10-08-the-declaration-and-the-question-underneath-it.md`.

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

### 4.8.2 Annotation (2026-10-01): the ~11x was compared against the wrong column, and the model is 14% off the right one

**No number in 4.8 or 4.8.1 above is withdrawn.** What changes is which
measured number the model's 2.59x should have been held against, and the
reason is that Passi et al.'s table has two columns and this chapter had
only ever read one.

Their TLM table reports contact resistance at **two gate biases**: the
on-state (V_BG = −40 V) and the **Dirac point**. The repository ingested
both arrays in 2026-08 and read only the on-state one; the Dirac-point
array was a cited, unread constant until `graphene_dead_name_sweep.py`
found it (`PASSI_RC_DIRAC_OHM_UM`, dead for ten sessions). Reading it:

| hole D (nm) | Rc on-state (Ω·µm) | Rc at Dirac (Ω·µm) | ratio |
|---|---|---|---|
| 0 (unpatterned) | 519 | 1372 | 2.64 |
| 50 | 212 | 620 | 2.92 |
| 100 | 352 | 732 | 2.08 |
| **200** | **45** | **456** | **10.13** |
| 500 | 410 | 1354 | 3.30 |
| 1000 | 560 | 1590 | 2.84 |

The ratio is the **gate-tunable multiple** of the contact resistance. It
sits between 2.08 and 3.30 at five of the six diameters and jumps to
**10.13** at D = 200 nm — the one diameter whose 11x this chapter quotes,
and **3.7x the mean of the other five**.

**So the headline reduction is a gate-bias statement, not a geometry
statement.** D = 0 → D = 200 nm is a **11.53x** reduction in the on state
and a **3.01x** reduction at the Dirac point, for the same two devices and
the same etched geometry.

**And 3.01x is the comparison this model's own physics selects.**
`R_extra` is computed from the *contact-induced* carrier density alone:
the metal's Fermi-level shift sets it, and the gate does not enter. At the
Dirac point the gate contributes no channel carriers, so the measured
access resistance is dominated by exactly the contact-doping-limited term
the model isolates; in the on state a −40 V back gate floods the channel
with holes and the measured Rc reflects a different balance entirely.
Holding a gate-independent model against an on-state measurement is a
category error, and it is the error 4.8 above makes.

Against the right column the model's **2.59x** (982.1 → 379.0 Ω·µm) is
within **14%** of the measured **3.01x**, where against the on-state
column it was short by a factor of **4.5**.

**Two caveats, stated rather than absorbed.** (i) The model's 2.59x is a
*pure-mode* ratio — fully top-mode against fully edge-mode — while the
D = 200 nm device is a mixture at a fill fraction Passi et al. do not
report (4.8.1). Treating the two as comparable assumes the patterned
contact is essentially edge-dominated at that diameter, which the
ligament-overlap criterion of 4.8.1 supports but does not measure.
(ii) The agreement is a *ratio* agreement; Section 4.7's finding that
`R_extra` runs higher than measured Rc in absolute terms is untouched by
this and still open.

**Consequence for 4.8.1's open item.** The small-D upturn the geometry
model fails to reproduce was attributed above to a fabrication effect
"outside this compact model's scope". That attribution is now much better
supported, and more specifically: whatever produces the D = 200 nm optimum
is **3.5x more of an effect in the on state than at the Dirac point**
(4.71x against 1.36x, taking D = 50 nm as the reference), so it scales with
*carrier density*, not with interface area. A model built entirely from
interface geometry — which 4.8.1's is, edge length per unit area at fixed
fill fraction — cannot produce a minimum whose depth depends on the gate.
Its failure was not a missing refinement; it was the wrong class of model
for that feature, and the evidence was in the column the chapter cited and
did not read.

Source: `graphene_edge_contact_model.py`,
`gate_tunable_fraction_of_contact_resistance()` and `dirac_point_summary()`;
sweep `graphene_dead_name_sweep.py`, log `dead_name_sweep_output.txt`.

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
