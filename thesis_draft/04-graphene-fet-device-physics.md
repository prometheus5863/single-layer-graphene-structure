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

(4.29) is the linearisation of (4.28). [**SUPERSEDED 2026-10-10, §4.6.8:**
this sentence is wrong and is kept as written. (4.29) is (4.28)'s *exact
differential relation* used algebraically — a one-point rectangle rule for
∫(4.31) — compounded with `C_q` evaluated at the overdrive rather than at
`E_F/e`. The two forms do not share a leading order as `dV → 0`: (4.28) is
quadratic, (4.29) linear. **The −0.6023 % below is unaffected and was
reproduced independently.**] Running Eq. (4) on (4.28) with nothing
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
[**SUPERSEDED 2026-10-10, §4.6.8:** read "its exact differential relation,
used algebraically" for "its linearisation". The statement that it is not an
additional physical effect stands.]
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

### 4.6.7 The sign was computed and discarded — and the disorder floor is 2.07× Chapter 3's

**2026-10-09.** §4.6.6 declared `V_ch` and left the exact charge relation (4.28)
as §7.9's top item, blocked by item 9. Item 9 is now closed, and closing it
answered two other open items and falsified a docstring. Module:
`graphene_fet_signed_carrier_model.py` (19/19 checks). Transcript:
`fet_signed_carrier_output.txt`. Figure: `fet_signed_carrier_model.png`.
Session notes:
`notes/2026-10-09-the-sign-was-computed-then-thrown-away-and-the-floor-was-never-a-floor.md`.

#### Why (4.28) could not be solved before today

Both terms on the right of (4.28) are **odd** in `n`: `e·n/C_ox` trivially, and
`E_F(n)/e = A_F·sign(n)·√|n|` because the Fermi level of a hole gas lies below
the Dirac point. A relation with two odd sides has an odd root, and an odd root
cannot come out of a solver whose first line is `dV = np.abs(dV)` — which is the
first line of the 10-08 solver. That solver returns the correct *magnitude*;
what no object in this repository could return was **which branch the channel is
on at a given point**.

`carrier_density()` is the origin. It *computes* the sign — `n_eff` is odd in
`dV` because `C_q` is even — and then discards it in the quadrature floor
`√(n_eff² + n_puddle²)`. Everything downstream takes `abs(n)`, so for the
conductivity the discard is harmless; for (4.28) it removes the branch.

`n_net_exact()` restores it. Validated against exactly known values, not ranges:
the root is odd **bitwise** over 801 points, and substituting it back into
(4.28) leaves a residual of **7.690e-16** relative over `dV ∈ [−2, 2] V`.
**§7.9 item 12 is unblocked.**

#### The decomposition, and what it reveals `n_puddle` to be

    n_net = n_e − n_h  (electrostatics, signed),  S = √(n_net² + n_puddle²)
    n_e = (S + n_net)/2,    n_h = (S − n_net)/2

- `n_e + n_h` **is** `S`, i.e. exactly what `carrier_density()` returns. So the
  quantity this chapter has been flooring since August is the **total conducting
  density** — not the net density, not a per-species density.
- `n_e·n_h = n_puddle²/4` to 3.603e-16 relative: the August quadrature floor,
  chosen for smoothness, is algebraically a **mass-action law**.
- `n_e(+dV) = n_h(−dV)` bitwise: graphene's ambipolar symmetry, machine-checked.

**Contract: no number in this chapter moves.** The two-carrier conductivity
equals the shipped conductivity to 3 ulp (3.609e-16 relative, against the four
significant figures this chapter quotes), and the shipped sweep is re-run after
the module's monkey-patching with the constants verified restored.

[Wiedmann et al., Phys. Rev. B **84**, 115314 (2011), arXiv:1107.3929] splits
net-from-gate and total-from-disorder in the same way and argues from the
divergence of `1/R_H` at neutrality that `n + p` must remain finite —
independent, measurement-side support for the identification above. It imposes
**no** electron–hole product, so the mass-action law is a *stronger* assumption
than the literature makes and is recorded as a limitation, not a result.

#### §4.9's floor item, open since 2026-09-29 and untouched for eleven sessions

With `n_puddle` named as the total, `dV = 0` gives `n_e = n_h = n_puddle/2`
exactly and the floor is a one-line calculation:

| | value | sheet resistance |
|---|---|---|
| this chapter: `σ_min = e·μ·n_puddle` | 3.204353 × 10⁻⁴ S/sq | **3.1208 kΩ/sq** |
| Chapter 3 §3.3, measured `4q_e²/h` | 1.549618 × 10⁻⁴ S/sq | **6.4532 kΩ/sq** |

**Chapter 4 is 2.0678× more conductive at neutrality than Chapter 3.** The
`n_puddle` that reconciles them is 2.417989 × 10¹⁵ m⁻² = 2.418 × 10¹¹ cm⁻²,
against the shipped 5 × 10¹¹ cm⁻². (§3.3's 6.45 kΩ/sq is reproduced here from
CODATA to 0.0032 kΩ/sq rather than re-typed, so the comparison cannot drift.)

Quoted as a factor first, per 10-03: a PASS/FAIL on "do the two floors agree"
returns FAIL and says nothing about whether the gap is 2× or 200×.

**This is decidable only because item 9 is closed.** A quadrature floor on a
*magnitude* cannot say whether it floors the net density, the total, or each
species — three readings differing by factors of 1, 1 and 2 in `σ_min`, the same
order as the answer. And the convention is not the explanation: under the
per-species reading the total would be `2·n_puddle` and this chapter would be
**4.1357×** too conductive, so that reading moves the discrepancy the **wrong
way**.

**What is not claimed.** `4q_e²/h` is a quantum/ballistic minimum conductivity
observed experimentally; `e·μ·n_puddle` is a diffusive disorder floor. They are
not required to be equal. The defect is that this thesis quotes both, in
chapters that feed each other, and had never compared them — and that the
diffusive one lands **below** the measured one, i.e. this chapter's channel at
neutrality is less resistive than any measured graphene sheet.

#### What the factor costs, and the contacts attenuate it again

At `V_ds` = 0.1 V, on-state `V_g` = 3.5 V, off-state at `V_g = V_dirac`:

| | shipped | Ch.3-matched | on/off multiplier |
|---|---|---|---|
| terminal (`R_c` = 600 Ω) | on/off **1.2408** | 1.7820 | **1.4362×** |
| intrinsic (`R_c` = 0) | on/off **1.6146** | 2.8024 | **1.7357×** |

`I_off` falls 35.19 % at the terminals, 51.58 % intrinsically. The 2.0678× floor
error reaches the on/off ratio as 1.7357× intrinsically and only 1.4362× at the
terminals: **the contacts attenuate it by a further 1.2085×** — §4.6.6's
negative-feedback structure, on a completely different parameter, which makes it
a property of the device topology rather than of the charge relation.

**§4.4's on/off ratio is therefore understated**, by up to 1.44× as measured at
the terminals, *if* §3.3's measured floor is the right one for this device.
§4.4's numbers stand and are not deleted; this is the limitation attached to
them. It does not disturb §7.2's conclusion — an on/off ratio of 1.78 is no more
a logic transistor than one of 1.24.

#### Annotation to §4.3 and §4.4, in place: `n_puddle` is not a Dirac-point term

§7.9 item 11, created 10-07, asked whether `n_puddle` removes other regimes.
It does not remove regimes; **it adds a floor to all of them.** Measured:

| `V_g` (V) | \|n_net\| (m⁻²) | `n_puddle`/\|n_net\| | floor adds to total |
|---|---|---|---|
| 0.90 | 2.357 × 10¹⁴ | 21.22 | +2023.9 % |
| 2.00 | 2.870 × 10¹⁵ | 1.742 | +100.9 % |
| 3.00 | 5.265 × 10¹⁵ | 0.950 | +37.9 % |
| **3.50** | **6.462 × 10¹⁵** | **0.774** | **+26.4 %** |

At the largest overdrive this thesis sweeps, the *residual* puddle density is
**77 % of the gate-induced density** and still adds **+26.4 %** to the total.
Its contribution falls below 1 % only at `dV` = **14.73 V**, 5.46× anything swept
here. `carrier_density()`'s docstring — the floor "regularizes `n` at the Dirac
point" — describes a regime this thesis never enters. The docstring's wording is
annotated in place rather than deleted.

This was found by a **check written to pass**: §5 of the module asserted that
`I_on` is untouched by an off-state parameter, on the docstring's authority, and
it failed at −6.92 %. 09-28 established that prose is a detector; today adds
that so is a sanity check one expects to pass.

#### Three exactness checks failed as first written, and none was a code defect

Recorded because the *repair* is the transferable part.

1. The ambipolar symmetry failed "bitwise" by 2.0 /m². The symmetry is not
   broken — it is bitwise in `dV`. `(V_dirac − x) − V_dirac` is not
   `−((V_dirac + x) − V_dirac)` in binary floating point (worst 2.220e-16 V at
   `x` = −2 V), and the `C_ox/e` lever of 2.395 × 10¹⁵ /(V·m²) turns that into
   2.0 /m². **Every sweep in this thesis is a `linspace` over `V_g` followed by
   a subtraction of `V_dirac`**, so any bitwise symmetry claim must name its
   parameterisation. The residual is accounted for, not thresholded: feeding the
   reconstructed `dV` back into the `dV`-level model reproduces the `V_g`-route
   species bitwise.
2. "`n_e − n_h` recovers `n_net`" failed at 1.772e-14 and passes at 1.5e-16 —
   the same 1-ulp residual, **119× apart**, because `n_net → 0` at neutrality
   while `S` stays at the floor. The denominator decided PASS from FAIL: §7.9's
   10-06 ratio-criterion item, in a check written by the audit itself.
3. The `σ` identity failed "bitwise" at 3 ulp and was **not** repaired by
   widening a ulp budget — 09-29's α-anchor policy — but by stating the identity
   and naming the slack against the precision this chapter commits.

#### A number item 13 produced on its first use outside its own module

§4.6.6 recorded the `cosh` overflow above `dV` = 18.3493 V as item 13, latent,
"no committed number affected". The first draft of the 1 %-crossover bisection
above scored the resulting `NaN` as "the floor dominates" (`nan > 0` is `False`)
and returned its own 200 V bracket ceiling as the answer; the correct value is
14.73 V. The bracket is now capped below the onset and `NaN` is a hard stop.
**Latent and harmless are not the same property.** Item 13's ranking is
unchanged — it is a numerics change inside a module with committed transcripts
across three chapters — but "no committed number affected" is no longer the
whole of its description.

#### What this section does not claim

- It does not claim `n_puddle` should be changed. The 2.418 × 10¹⁵ m⁻² figure is
  what would reconcile two floors that are not required to agree; adopting it
  would move §4.4, §4.7 and §4.10 and is not done here.
- It does not claim the mass-action closure is physical. It is algebraically
  equivalent to a smoothness choice, and the literature it most resembles
  declines to impose it.
- It does not claim the signed model is a new physical model. Its contract is
  the opposite: every quantity it introduces is a decomposition of one already
  being computed, and the conductivity is unmoved to 3 ulp.
- It does not re-derive anything on (4.28). That is item 12, now unblocked and
  still open. [**CLOSED 2026-10-10, §4.6.8.**]

### 4.6.8 Chapter 4 re-derived on (4.28) — and (4.29) was never a linearisation

§4.6.6 declared `V_ch` the quasi-Fermi potential, made (4.28) the exact charge
relation, measured the gap to (4.29) at −0.6023 % of `I_d`, and left the
re-derivation as §7.9 item 12. §4.6.7 supplied the branched root the item was
blocked on. This section does the re-derivation. **The −0.60 % reproduces
exactly; the finding is that the chapter has been describing the pair of
equations wrongly, in three places, since 2026-10-08.**

Machine output: `exact_charge_rederivation_output.txt` (21 checks, 21 passed).
Module: `graphene_exact_charge_rederivation.py`. Figure:
`exact_charge_rederivation.png`.

#### (4.29) is not the linearisation of (4.28), and the test is one line

A linearisation of a function agrees with it to first order by construction.
These two do not share a *leading* order. In (4.28), `a·n + b·√n = dV` with
`a = e/C_ox` and `b = A_F`; as `dV → 0` the `√n` term dominates and
`n → (dV/b)²` — **quadratic**. In (4.29), `C_q(0)` is finite, so
`n → (C_ox/e)·f(C_q(0))·dV` — **linear**. Measured log-log exponents:

| `dV` (V) | d ln n/d ln dV, (4.28) | d ln n/d ln dV, (4.29) |
|---|---|---|
| 10⁻¹⁰ | **1.99999854** | **1.00000000** |
| 10⁻⁶ | 1.999387 | 1.000000 |
| 10⁻⁴ | 1.943766 | 1.000000 |
| 10⁻² | 1.274493 | 1.002201 |

#### What (4.29) actually is: the exact differential relation, used algebraically

Differentiating (4.28),

```
d(dV)/dn = e/C_ox + (1/e)·dE_F/dn = e/C_ox + e/C_Q(n),
C_Q(n) := e²·dn/dE_F = 2·e·√n / A_F                                  (4.30)

  =>  dn/d(dV) = (C_ox/e) · C_Q(n)/(C_Q(n) + C_ox)                   (4.31)
```

(4.31) **is** (4.29) with `dn/d(dV)` in place of `n/dV`. So the series factor
is not an approximation *of* the charge relation; it is the charge relation's
exact derivative, and (4.29) uses it as a **one-point rectangle rule** for
`∫₀^dV (4.31) dV'`. That is a quadrature error, not a truncated expansion, and
a rectangle rule is worst where the slope varies fastest — near the Dirac
point, which is exactly where §4.6.6's raw disagreement reached +68 %.

(4.30) and (4.31) are verified as identities, not asserted: `e·dn/d(dV)` from
(4.28)'s closed-form root equals the series combination of `C_ox` with `C_Q(n)`
to **1 ulp**, and numerically integrating (4.31) back to the closed-form root
converges at **observed order 3.994** against Simpson's 4, with the `h⁴`
extrapolation landing **1.5×10⁻¹³** relative from it. An error in the prefactor
of `C_Q` — a dropped factor of 2, say — would leave an O(1) plateau and drive
the observed ratio to 1 instead of 16; that is what the order, rather than a
threshold on the error, is testing.

#### Two independent errors, separated for the first time

(4.29) differs from (4.28) in two unrelated ways: the quadrature above, and
`quantum_capacitance()` evaluating `C_q` at the gate overdrive `dV` instead of
at `E_F/e = A_F·√n` (noticed in §4.6.5's 2026-10-07 annotation and never
costed). Running the 2×2 — {rectangle, integral} × {argument `dV`, argument
`E_F/e`} — on the raw charge, with the puddle floor off so that the floor does
none of the work:

| `dV` (V) | total error of (4.29) | quadrature alone | **argument alone** | sum |
|---|---|---|---|---|
| 0.01 | +68.2193 % | +32.1722 % | **+68.0927 %** | +100.2649 % |
| 0.50 | +8.0543 % | +4.1176 % | **+7.2144 %** | +11.3320 % |
| 1.20 (RF bias) | **+5.2064 %** | **+2.6396 %** | **+4.7416 %** | +7.3812 % |
| 2.70 | +3.4729 % | +1.7523 % | **+3.2193 %** | +4.9716 % |

> **SUPERSEDED 2026-10-11, §4.6.9(c). Numbers retained, not rewritten.**
> Every figure in this table is differenced against (4.28), which is a T = 0
> relation, while the model is a 300 K one — so each of them contains a
> temperature difference as well as the knob it names. Re-measured in a 2×2
> at one temperature, the RF-bias row reads **argument +2.6408 %, quadrature
> +0.7603 %, total +3.0963 %**: this table overstates them by 1.80×, 3.47×
> and 1.68× respectively. The *ordering* claimed below survives — the
> argument error is still the larger contributor — and the conclusion that
> they are not additive survives, though the departure from additivity is
> 1.109 rather than 1.42 once the reference is at the model's temperature.
> The mistake is kept visible because its shape (a decomposition measured
> against a reference at a different operating point) recurs, and because
> §4.6.9's 300 K argument error of +2.6408 % is numerically within 0.0012 %
> of this table's quadrature figure, which reads as a reproduction and is not
> one.

Three things follow. **The argument error is the larger single contributor**
(+4.74 % against +2.64 % at the RF bias) — which no earlier session could have
known, because none separated them. **Neither is the explanation on its own.**
And **they are not additive**: the sum overstates the total by 1.42× at the RF
bias and by 1.47× at `dV` = 0.01 V, because they compose multiplicatively in
the series factor. A decomposition reported as a sum would have been wrong in
the direction that flatters it.

The puddle floor then attenuates all of this by **4.3×** at the RF bias
(+5.21 % raw → +1.22 % in what `carrier_density()` returns). That attenuation
is a property of `n_puddle`, which §4.6.7 measured to be 2.0678× off Chapter
3's floor — so it is a reason the re-derivation moves little, not a reason for
confidence in either form.

#### §§4.3–4.4 re-derived

| quantity | on (4.29) | on (4.28) | change |
|---|---|---|---|
| `n` at `V_g` = 3.5 V (m⁻²) | 8.170466×10¹⁵ | 8.000036×10¹⁵ | **−2.0859 %** |
| `R_ch` at `V_g` = 3.5 V (Ω) | 381.9555 | 390.0926 | +2.1304 % |
| `R_ch` at the Dirac point (Ω) | 624.1509 | 624.1509 | **0.0000 %** |
| `I_d`, drift, RF bias (A) | 8.719345×10⁻⁵ | 8.671268×10⁻⁵ | −0.5514 % |
| `I_d`, saturated, λ = 0, RF bias (A) | 8.190217×10⁻⁵ | 8.140885×10⁻⁵ | **−0.6023 %** |
| `L_x`, Pd 110 Ω·µm (nm) | 115.20 | 112.79 | −2.0859 % |
| `L_x`, 300 Ω·µm (nm) | 314.17 | 307.62 | −2.0859 % |
| `L_x`, 500 Ω·µm (nm) | 523.62 | 512.70 | −2.0859 % |
| terminal on/off | 1.240756 | 1.230617 | −0.8171 % |
| intrinsic on/off (`R_c` = 0) | 1.614570 | 1.581261 | −2.0630 % |

The −0.6023 % line **reproduces §4.6.6's published figure to four decimal
places from an independently written patch path**. Per 2026-10-02 that is a
fact about two artefacts rather than a validation; what it establishes is that
the (4.29) column above is the committed model and not a re-run with drift in
it.

At the Dirac point the two forms agree to 0.0000 % because the floor supplies
the entire density there — the one place in this chapter where the charge
relation does not matter at all, and §4.6.7 already showed it is also the
place the chapter's off-state numbers come from.

**Every qualitative claim of §§4.3–4.4 survives.** The best literature contact
still crosses over below the 200 nm channel and the mid- and worst-case
contacts still do not; the on/off ratio is still between 1 and 2, i.e. still
not a switch; the channel is still more resistive at neutrality than in the on
state. §4.4's contact-engineering conclusion and §7.2's "no switch" are
untouched.

**On/off moves less at the terminals (−0.8171 %) than intrinsically
(−2.0630 %).** This is the third parameter in three sessions to show that
structure — the diffusion term (2026-10-08), `n_puddle` (2026-10-09), and now
the charge relation itself. It is therefore a property of the **device
topology**, not of any one parameter: the series contact resistance is a
negative feedback on `I_d` and attenuates whatever moves upstream of it.

#### §4.6 re-derived — and the item was mis-scoped

Under (4.31) the small-signal gate capacitance is **not an independent
ingredient**: `e·dn/d(dV)` *is* `(1/C_ox + 1/C_Q(n))⁻¹`. So a consistent
re-derivation moves `C_gs` as well as `I_d`, and `f_T = g_m/(2πC_gs)` sees
both. Applied separately, at the peak-`f_T` bias, `V_ds` = 0.1 V:

| | on (4.29) | (4.28) charge only | (4.28) charge + `C_Q` |
|---|---|---|---|
| peak-`f_T` bias `V_g` (V) | −1.1115 | −1.1992 | **−1.1867** |
| `f_T` (GHz) | 20.2787 | 19.8734 | 20.2666 |
| `f_max` (GHz) | 18.5566 | 18.1781 | 18.5495 |
| `C_gs` (fF) | 0.0767 | 0.0767 | 0.0752 |
| change in `f_T` | — | **−1.9986 %** | **−0.0594 %** |
| change in `f_max` | — | −2.0399 % | −0.0386 % |

**The combined −0.06 % is a cancellation, not a smallness.** The charge half
is −1.9986 % and the capacitance half +1.9392 % — opposite in sign and within
3 % of each other in magnitude. A re-derivation that changed only the charge
would have reported −2.0 % and been wrong by a factor of 34; one that reported
only the combined number would have called item 12 negligible in `f_T` when it
is in fact two coupled 2 % corrections that happen to oppose. **An item that
names one of two coupled quantities gets the other one whether or not it asks,
and the sum is not the way to find out which mattered.**

The capacitance half is small for a reason worth stating: `C_q` at the
overdrive is 4.500×10⁻¹ F/m² and `C_Q` at `E_F/e` is 1.820×10⁻² F/m², a factor
of **24.7** apart, but both are ≫ `C_ox` = 3.837×10⁻⁴ F/m², so a factor of 25
in `C_q` becomes −1.94 % in `C_gs`. `E_F/e` = 0.0773 V is 4.04 % of the 1.91 V
overdrive, which is the argument error of the previous subsection seen
directly.

**The peak bias is the sensitive quantity, not the peak value.** `f_T` moves
0.0594 % in height and 6.76 % in *location* (−1.1115 → −1.1867 V), a ratio of
**114**. Both peaks sit on the **hole branch** and neither crosses to the
electron branch, so the open "does the peak-`f_T` bias change branch between
the two models" question of §4.6.7 is **answered NO for the (4.28)-vs-(4.29)
pair**; it was asked about the branched-vs-magnitude pair and stays open there.

At the literature-scale geometry (40 µm, 8 fingers) that §4.6.2 and §7.8.1a
quote: peak `f_T` 20.2787 → 20.2666 GHz, `f_max` at that bias 13.0393 →
13.0346 GHz, extrinsic `f_T` 3.4423 → 3.3847 GHz, extrinsic `f_max` 2.2105 →
2.1727 GHz. **§7.8.1a is untouched, and for the right reason**: its verdict is
structural — that `f_max` is a restatement of `R_total/(R_g+R_s)` — so the
quantity to check is the ratio, and `f_max/f_T` moves **+0.0228 %**, from
0.643007 to 0.643154, still below the bottom of Feijoo *et al.*'s 1.3–1.4 band
that 2026-10-04 showed a *perfect* gate cannot reach. Five levers have now
been tried — velocity saturation, the perfect contact, the diffusion term, the
exact charge relation, and its capacitance — and the binding constraint is
still that `transfer_characteristic()` is a resistor with no output resistance.

#### Three exactness checks failed as first written, and two were wrong in their prose

Continuing 2026-10-09's series, with a new instance of it.

**(a) The quadratic-exponent check failed at `dV` = 10⁻⁴ V**, returning 1.9438
against a ±0.01 criterion. It did not measure the exponent. The two terms of
(4.28) are equal at `dV_x = 2·A_F²·C_ox/e` = **6.519 mV**, so 10⁻⁴ V is only
65× below the crossover, where the correction is O(√(dV/dV_x)) ≈ 12 %. **It
measured the width of the crossover**, and the repair is to evaluate the limit
where the limit is and to print `dV_x` so the criterion comes from the model
instead of from a round number.

**(b) The integral identity failed at 1.1×10⁻⁷ with order 3.49 on a uniform
grid — and it measured the same 6.519 mV.** The integrand of (4.31) turns over
inside a layer of width `dV_x` at the lower limit, and a uniform grid over
[0, 1.2 V] puts only 22 nodes in it even at 4001 points. A grid graded about
`dV_x` reaches order 3.994 with **12× fewer** points. **Two exactness checks
written for different purposes, in different sections, both failed by measuring
the same number** — which is the strongest evidence this thesis has that
2026-10-09's rule is a rule and not an anecdote.

**(c) "Every number moves by under 2 %" failed at 2.0859 %, and the bound was
removed rather than widened.** The 2 % was a round number chosen before the
run and derived from nothing, so widening it to 2.5 % would have been exactly
the tuned tolerance §7.9 item 15 exists to find. The magnitudes are now
printed with no threshold on them and the check tests the qualitative claims,
which are what a re-derivation can actually overturn. **Of 2026-10-09's three
available repairs — loosen, fix the model, find out what it measured — this
session used the third twice and the fourth option, delete the criterion, once.**

Two further checks were wrong in their *prose* and failed on it: one asserted
the capacitance half was the *larger* half (it is 0.970× the charge half), and
one asserted `f_T` was "two orders below the 100–300 GHz RF-frontend
requirement" — wrong twice, since 20.3 GHz is about one order below 100 GHz
and §7.8.1a does not rest on a frequency threshold at all. 2026-09-28
established that prose is a detector; here it detected the author's own
checks.

#### A delivery fault found on the way

`contact_resistance_crossover.py` opens with `from graphene_fet_model import
carrier_density`, so it holds its **own** reference, bound at import time.
Rebinding `graphene_fet_model.carrier_density` does not reach it, and a
re-derivation that patched only the owning module would have left §4.4 on
(4.29) while reporting it as re-derived. This is the 2026-10-03
frozen-default fault in its monkey-patch form. An **arrival control**
(2026-09-30's rule) now sits in the module ahead of the §4.4 numbers and fails
if the patch did not reach the module that computes them.

#### What this section does not claim

1. It does **not** claim (4.29) is less accurate than (4.28). Neither has been
   compared to measurement. (4.29) is non-exact under the declaration and is
   now additionally *named* — a quadrature of the exact differential — which
   is a statement about its derivation, not about its accuracy.
2. It does **not** claim (4.28) is exact physics. `E_F = ħv_F√(πn)` is the
   **T = 0** dispersion and carries no thermal broadening; (4.28) is exact
   with respect to the declaration and the Dirac dispersion only. The
   compact-model review [arXiv:2209.00388](https://arxiv.org/pdf/2209.00388)
   uses `C_q ≈ k·c₁·√(1 + (V_c/c₁)²)`, whose `c₁ → 0` limit is exactly (4.30),
   and `c₁` is the thermal rounding (4.28) drops. That regime is far inside
   the puddle floor and is stated as a limitation rather than costed.
3. It does **not** replace Chapter 4's numbers. Every committed number remains
   a (4.29) number, annotated with its (4.28) value above. The chapter's
   production path is unchanged and is now carrying a measured, decomposed
   statement of what that choice costs.
4. It does **not** touch the hole-branch sign of this chapter's
   branch-asymmetric results. §4.6.7 made that expressible; it remains
   unverified.
5. It does **not** change the floor. `n_puddle` is carried into (4.28)
   unchanged so the comparison is apples to apples, which means (4.28) + floor
   is not itself an exact model either.

#### Annotation to §4.6.6 and §4.3, in place

§4.6.6's sentence "(4.29) is the linearisation of (4.28)" and the matching
wording in `carrier_density()`'s 2026-10-08 annotation are **marked
superseded where they stand, not rewritten**: (4.29) is the exact differential
relation of (4.28) used algebraically, compounded with an argument error. No
number in §4.6.6 changes — the −0.6023 % it published is reproduced here — and
the superseded wording is kept visible because the shape of the mistake
(describing a quadrature as a truncation) is the kind that recurs.

### 4.6.9 `C_q`'s argument: the fix is a self-consistent solve, the obvious version of it deletes the temperature, and §4.6.8's decomposition was measured against the wrong reference

§4.6.8 named a second error alongside the quadrature error of (4.29) —
`quantum_capacitance()` is evaluated at the gate overdrive `dV` while the
dispersion it descends from requires `E_F/e = A_F·√n` — sized it at
**+4.7416 %** of the raw charge at the RF bias, called it the largest cost
attached to anything then open on this chapter, and made it §7.9 item 16.
This section acts on it. **The item's direction is right and all three of its
numbers are wrong, by factors of 1.68 to 3.47.**

Machine output: `quantum_capacitance_argument_output.txt` (18 checks, 18
passed, five against exactly known values). Module:
`graphene_quantum_capacitance_argument.py`. Figure:
`quantum_capacitance_argument.png`.

#### (a) The argument cannot be fixed by changing an argument

`quantum_capacitance(V)` takes a potential. The quantity it must be evaluated
at is `E_F/e = A_F·√n`, and `n` is what `C_q` is being used to compute. The
correct argument is a **state**, not a voltage, so the fix is a
self-consistent solve

    n = (C_ox·dV/e) · C_q(E_F(n)/e) / (C_q(E_F(n)/e) + C_ox)             (4.32)

and not an edit to a call site. That is the reason the largest open cost on
this chapter sat for a session as a one-line description: the item as written
understated the work, not the size.

At the RF bias, on the committed charge, `E_F/e = 0.062498 V` — **5.2082 %** of
the 1.2 V overdrive — and the two candidate arguments give

| argument | `C_q` (F m⁻²) | `C_q/C_ox` | series factor |
|---|---|---|---|
| `dV` = 1.2 V (committed) | 2.825142e-01 | 736.3 | 0.998644 |
| `E_F/e` = 0.0625 V | 1.575324e-02 | 41.06 | 0.976223 |

A **17.9×** difference in `C_q` is a **2.2966 %** difference in the series
factor. The reason the error is small is that both values are far above
`C_ox`, i.e. **the device is not quantum-capacitance-limited at this bias** —
not that the argument is nearly right. Where it would matter is exactly where
this chapter has least to say: near the Dirac point, where the puddle floor
already dominates.

#### (b) The obvious fix is a T = 0 expression, and it deletes the rounding the function exists for

The dispersion `C_Q(n) = 2e√n/A_F` is a zero-temperature result.
`quantum_capacitance()` is thermally rounded for a stated reason — its own
docstring says *"so C_q does not diverge to zero exactly at the neutrality
point"*. Substituting the dispersion removes that:

| | `C_q` at the Dirac point (F m⁻²) |
|---|---|
| committed thermal, at `dV` = 0 | 8.43739892e-03 |
| T = 0 dispersion, at `n` = 0 | **0** |
| thermal at `E_F/e`, at `n` = 0 | 8.43739892e-03 |

`C_q(0)` has a closed form, `PREFACTOR·(kT/e)·ln 4`, because `ln(2(1+cosh 0)) =
ln 4` exactly; the thermal-at-`E_F/e` form reproduces it to **0 ulp**. The
dispersion gives zero, hence a series factor of zero and an induced charge of
*exactly* zero at the Dirac point — the floor would supply 100 % of the
density there instead of its committed 99.9 %.

A sharper form of the same statement: with the dispersion, `n ≡ 0` solves
`dn/d(dV) = (C_ox/e)·f(C_Q(n))` through the origin, so the initial-value
problem is **degenerate** there. With the thermal form it is not. *Temperature
is what makes the charge relation well-posed at the Dirac point*, and it is
what the obvious fix throws away. §4.6.8 sized item 16 with that substitution,
so **its +4.7416 % is the cost of a fix that is wrong in a second way.**

The form used here, (4.32) with the thermal `C_q` at `E_F/e`, is called `B_T`.
It keeps the temperature, reduces to the dispersion as `T → 0` (checked, §(d)),
and keeps a finite `C_q` at neutrality.

#### (c) The decomposition, re-measured at one temperature

(4.28) is a T = 0 relation and the committed model is a 300 K one, so every
percentage in §4.6.8's 2×2 contains a temperature difference as well as the
knob it names. The repair is a 2×2 at **one** temperature — argument ∈ {`dV`,
`E_F/e`} × quadrature ∈ {rectangle, exact} — with the thermal `C_q`
throughout, and a reference `D_T` obtained by integrating the model's own
differential relation as an initial-value problem. That integrator is
validated against an exactly known value: with the dispersion in place of the
thermal form it must reproduce (4.28)'s closed-form root, and does, to
**−1.5e-14** relative.

Raw charge at `dV` = 1.2 V:

| | `n` (m⁻²) | vs `D_T` |
|---|---|---|
| A `C_q(dV)`, rectangle **(committed)** | 2.86980420e15 | +3.0963 % |
| `B_T` `C_q(E_F/e)`, rectangle | 2.80477739e15 | +0.7603 % |
| `C_T` `C_q(dV)`, exact quadrature | 2.85712585e15 | +2.6408 % |
| `D_T` `C_q(E_F/e)`, exact quadrature | 2.78361495e15 | 0 |

| | §4.6.8, vs T = 0 (4.28) | here, vs 300 K `D_T` | overstated by |
|---|---|---|---|
| argument error | +4.7416 % | **+2.6408 %** | 1.80× |
| quadrature error | +2.6396 % | **+0.7603 %** | 3.47× |
| both wrong | +5.2064 % | **+3.0963 %** | 1.68× |

§4.6.8's three figures are **marked superseded where they stand and not
deleted**, per this chapter's standing practice. What survives is the
*ordering*: the argument error is still the larger of the two, so item 16 is
still correctly ranked ahead of the quadrature error. What does not survive is
the headline number this chapter quotes in two places.

**A trap, recorded because it is the kind that is read as a confirmation.**
Today's argument error at 300 K (+2.6408 %) agrees with §4.6.8's *quadrature*
error measured at T = 0 (+2.6396 %) to 0.0012 %. They are different quantities
against different references. Anyone reading the two sections in sequence
would take the agreement for a reproduction. §4.6.6's rule — that a
coincidence between two artefacts is not a validation — is applied here for
the first time to a **numerical** coincidence rather than to a re-run.

#### (d) Item 13 was ranked last, and it is what blocked this section

The natural exact-limit check on `B_T` is `T → 0`: it must reproduce the
dispersion fixed point. It **could not be run on the committed code**.
`η = eV/kT` grows as `T` falls, `cosh(η)` overflows at **η = 710.4758**,
`ln(2(1+∞)) = ∞`, and the series factor returns `∞/∞ = nan`. That is §7.9 item
13, added 2026-10-08, judged *latent* and ranked **last** because η = 710 needs
**18.35 V** of overdrive at 300 K, which the device never sees. At T = 1 K and
the ordinary RF bias, η = 719.

> **Latency was judged against the inputs the model takes, not against the
> inputs its validations take.** A defect unreachable from a device's
> operating range can still be reached by a check on that device's model, and
> ranking by reachability missed it. The item ranked last gated the item
> ranked as the largest cost.

Item 13 is closed in the same session (§4.6.10). With the overflow-free
spelling in place the limit runs: at T = 10⁻³ K, `B_T` reproduces the
dispersion fixed point to **7e-16** relative, against `nan` before.

#### (e) What actually moves

`B_T` installed, with an arrival control on each patched name:

| | (4.29), `C_q(dV)` | `C_q(E_F/e)` | change |
|---|---|---|---|
| `n` at 3.5 V (m⁻²) | 8.170466e15 | 8.087032e15 | **−1.0212 %** |
| `R_ch` at the Dirac point (Ω) | 624.150907 | 624.150907 | **0.0000 %** |
| `L_x` Pd / 300 / 500 (nm) | 314.17 / 314.17 / 523.62 | 310.96 / 310.96 / 518.27 | −1.0212 % |
| terminal on/off | 1.240756 | 1.235801 | −0.3993 % |
| intrinsic on/off | 1.614570 | 1.598207 | −1.0135 % |
| `I_d` at the RF bias (A) | 8.190217e-05 | 8.169499e-05 | **−0.2530 %** |

The raw charge moves **−2.2659 %** and `n` at 3.5 V moves **−1.0212 %**: the
puddle floor attenuates by a factor **0.451**, which is why item 16 was
survivable as a change that moves committed numbers. `R_ch` at the Dirac point
is unchanged to the printed digit because the floor supplies the whole density
there.

Every qualitative claim of §§4.3–4.4 survives: the on/off ratio stays between
1 and 2, `L_x` at 300 Ω·μm stays above the 200 nm channel (314.17 → 310.96 nm),
and the resistance ordering by metal is unchanged. **The correction moves
numbers and no verdicts.**

RF, per-μm convention at `V_ds` = 0.1 V:

| | `C_q(dV)` | charge only | charge + `C_gs` |
|---|---|---|---|
| peak-`f_T` bias (V) | −1.1115 | −1.1617 | −1.1491 |
| `f_T` (GHz) | 20.2787 | 20.0574 | 20.4422 |
| `C_gs` (fF) | 0.0767 | 0.0767 | 0.0752 |

`f_T`: charge alone **−1.0911 %**, the consistent `C_gs` **+1.8976 %**, total
**+0.8065 %**. Opposite signs again. §4.6.8 found the same structure for item
12 and recorded it as a property of the **device topology** rather than of any
one parameter, on a single instance; it now has a second, independent one,
which is what that claim needed. `f_max/f_T` moves **+0.0009 %** and stays at
0.9151, so §7.8.1a is untouched by a **sixth** lever.

#### (f) The two-knob check's first use, and it fails in the direction the check cannot see

§4.6.8 created a backlog item for a check that can see a cancellation: move
each input separately, and flag the case where the combined move is *smaller*
than either. Its first instance is here, and it has the **opposite** sign.

| reference | argument only | quadrature only | sum | both | both/sum |
|---|---|---|---|---|---|
| 300 K, consistent | −2.2659 % | −0.4418 % | −2.7077 % | −3.0033 % | **1.109** |
| T = 0 (§4.6.8) | −2.4398 % | −0.4418 % | −2.8816 % | −4.9488 % | **1.717** |

The two knobs **reinforce**: the combined move is larger than the sum of the
separate ones, so a check that tests for smallness would not have flagged it.
The item is **amended, not closed** — it must report the ratio and flag any
departure from 1 in either direction.

And a requirement on the item that could only be found by running it: **the
interaction measure inherits the reference error of whatever it is differenced
against.** The T = 0 frame puts the departure from additivity at 0.717 where
the consistent frame puts it at 0.109 — an inflation of **6.6×**. A two-knob
check must be run at the model's own operating point.

#### (h) The literature regularises `C_q` with disorder, not with temperature — and here the two are the same size

Live search reached a **second, independent** compact-model lineage, the
[Southampton bilayer-GFET model](https://eprints.soton.ac.uk/366757/1/GFetIEEEtran.pdf).
It confirms §(a) from a different source than §4.6.8's
[Pasadas/Jiménez review](https://arxiv.org/pdf/2209.00388): its `C_q` is a
function of the **surface potential** `φ_s`, never of the gate overdrive, and
its layer potentials form a **coupled system solved together**. The claim that
item 16 is a self-consistency problem is now held by two independent groups
rather than inferred from one.

But it regularises `C_q` at neutrality **differently**, and the difference
opens a question this chapter has not asked. Temperature never enters its
`C_q` at all; it floors the capacitance with the residual disorder density,
`C_q,min = e²√n₀/(√π·ħ·v_F)`. The puddle density does the job this thesis
gives to `kT`.

| | `C_q` at the Dirac point (F m⁻²) |
|---|---|
| this thesis, thermal rounding | 8.437399e-03 |
| literature form, at `n_puddle` = 5e15 m⁻² | 9.710793e-03 |
| ratio | **1.1509** |

They agree exactly at `n₀ = 3.77e15 m⁻²` (0.755× `n_puddle`) or, holding `n₀`,
at **T = 345.28 K** (1.151× 300 K). The agreement is therefore not a
coincidence of one parameter: **SiO₂-supported graphene at room temperature
sits where disorder and temperature regularise `C_q` equally.**

And §4.3 applies **both**. The thermal rounding is inside `C_q`, and
`n_puddle` is then added in quadrature to `n`. The literature form shows those
are two spellings of one physical effect at this `n₀`, so **using both is a
double count of the same regularisation, worth about a factor 2 in `C_q` at the
Dirac point, and it has never been costed** (§7.9 item 18). It would not be
true for suspended graphene (`n₀` two orders lower) or at 4 K, where the two
conventions diverge — which is why neither is obviously wrong in the literature
and why the overlap here is easy to miss.

This also bounds what §(a)'s measurement can say about neutrality. The region
where the argument error would be *large* is precisely the region where two
overlapping regularisations are fighting over the same physics.

#### (g) What this section does not claim

1. It does **not** install `B_T` as the chapter's charge relation. The default
   remains (4.29); `B_T` is measured, committed and available. Flipping the
   default would move every number in Chapters 4, 5 and 6 and is its own
   session's work, with its own transcript regeneration.
   *[Annotated 2026-10-21, §4.6.17. Two corrections to this sentence, neither
   of which changes the claim that the default is unflipped. First, §4.6.12
   measured that Chapters 5 and 6 never call the charge path at all, so "every
   number in Chapters 4, 5 and 6" is wrong in that direction and "more of
   Chapter 4 than was measured" in the other. Second — and this is what
   §4.6.17 adds — "its own session's work" was read as meaning the work is
   the flip. It is not. The two halves of the flip are each correct alone and
   **do not compose**: the contact replacement reads the shipped
   `gfet.quantum_capacitance` at call time, which is the name the channel half
   rebinds, and installing both at once costs −51 % to −62 % in `n_contact`,
   4.86× to 10.77× the effect being installed. The numbers in this subsection
   stand; what was incomplete was the account of what flipping requires.]*
2. It does **not** claim `D_T` is exact physics. It is the exact integral of
   *this model's* differential relation at 300 K. The relation itself still
   uses a phenomenological rounding rather than Fermi–Dirac integrals, and
   (4.28) remains the T = 0 limit of a standard literature form.
3. It does **not** re-open §4.6.8's −0.6023 % in `I_d`. That number compared
   (4.29) against (4.28) and is reproduced unchanged; what is corrected is the
   *attribution* of the raw-charge error between the two knobs.
4. It does **not** claim the argument error is negligible where it matters.
   2.64 % at the RF bias is a statement about a bias where `C_q ≫ C_ox`. The
   quantum-capacitance-limited regime, where it would be large, is not studied
   in this chapter.

### 4.6.10 Item 13 closed, and a neutrality pass that measured the wrong thing three times

The rounding factor in `quantum_capacitance()` is now
`|η| + 2·log1p(exp(−|η|))`, an identity with `ln(2(1+cosh η))` because
`2(1+cosh x) = 4cosh²(x/2)`, and finite for every real η. The committed line
is retained as a comment above it.

Machine output: `quantum_capacitance_overflow_output.txt` (5 checks, 5 passed).
Module: `graphene_quantum_capacitance_overflow.py`.

**The result, stated as a measurement rather than as neutrality.** `C_q`
itself moves by **1 ulp at 4 of 400** sampled gate biases, all within 0.093 V
of the Dirac point where η is small and the two spellings round differently in
the last bit. **No committed quantity moves by a single bit**: `n`, `R_ch`,
`I_d` at both `V_ds` of record, `f_T`, `f_max`, `g_m`, `g_ds`, `C_gs` and the
crossover lengths are bit-identical over 3603 values. Saying it that way
rather than "output-neutral" is §4.6.6's prose-is-a-detector rule applied to
this section's own summary: a 1-ulp move in the function is a real difference,
and the claim that survives is about what it reaches.

**The pass failed three times before it measured the spelling.** Each draft
re-typed the old function body to compare against:

1. `e`, `ħ` and `k_B` entered as literals. `graphene_fet_model` takes them
   from `scipy.constants`, whose `ħ = h/2π = 1.05457181764615647e-34`; the
   literal `1.054571817e-34` differs at the tenth digit, squared in the
   prefactor that is ~1.2e-9, and the pass reported **3983 of 4003** values
   differing by up to 1.1e7 ulp. Sixth instance of §4.6.7's rule (find out
   what a failed check measured), and the first in which the substituted
   quantity was a *physical constant* rather than a model parameter.
2. Constants fixed, product regrouped as `prefactor·(kT/e)·log(…)` against the
   module's `smoothed_V = (kT/e)·log(…)` then `prefactor·smoothed_V`. Floating
   multiplication is not associative; the regrouping alone moved **118 of 400**
   `C_q` values by up to 2 ulp.
3. The patch meant to fix the association silently matched nothing, so draft 2
   ran again and was very nearly read as a confirmed result.

> **A function that "restores the old spelling for comparison" is a
> re-implementation, and every detail of it not copied exactly — the
> constants, the grouping, the order of operations — becomes part of the
> difference it reports.** The repair is not a better transcription. It is to
> make the thing that changed *substitutable*: `_thermal_rounding()` now holds
> the single expression item 13 touched, the pass swaps that one function, and
> the rest of the call path is shipped code running unmodified. That form
> reproduces the direct file-against-file result exactly — 4 of 400, 1 ulp —
> where all three transcriptions did not.

**And a correction to this section's own first draft.** It claimed the two
spellings agree "to ≤ 1 ulp over η ∈ [0, 700]". That was measured on eleven
hand-picked values; on 401 it is **2 ulp**. The ≤ 1 criterion was a round
number chosen before the run — precisely what §7.9 item 15 exists to find —
and it is **not widened to 2**. For an *identity*, a ulp count cannot answer
the question that matters, so `mpmath` at 60 digits adjudicates, and the answer
is not the expected one:

| | worst relative error | mean | strictly better | strictly worse |
|---|---|---|---|---|
| committed `ln(2(1+cosh η))` | 1.793e-16 | 4.643e-17 | — | — |
| installed `|η| + 2·log1p(e^−|η|)` | **1.527e-16** | 5.162e-17 | 45 / 401 | **81 / 401** |

The installed spelling is **worse at more samples than it is better**, with a
slightly larger mean error. All it has is the smaller worst case. **Accuracy is
therefore not a reason to install it**; being finite where the other is not is
the only reason, which is the only reason item 13 ever gave. The claim stops
there.

### 4.6.11 The double regularisation, costed: the two regularisers never meet in (4.29), and where they do they compose max-like

Section 7.9 **item 18** was created 2026-10-11 out of a literature reading,
not out of an internal audit, and it was ranked **first** on the next-run list
because it was the only open item that could change a *physical* claim rather
than a numeric one. Its premise was that this chapter regularises `C_q` at
neutrality with temperature *and* regularises `n` with the disorder density
`n_puddle`, that at this operating point the two are the same size, and that
applying both therefore double-counts one physical effect by about a factor 2.

Costing it closes it, and it closes by having been **wrong twice over**. The
measurement is `graphene_double_regularisation.py`; machine output
`double_regularisation_output.txt` (19 checks, 19 passed); figure
`double_regularisation.png`.

#### (a) In the committed relation the two regularisers are disjoint

The premise "this thesis applies both" is false for (4.29), and the reason is
structural rather than numerical.

`C_q`'s argument in (4.29) is the gate overdrive `dV`, and `dV` is never
floored. `n_puddle` enters only the final quadrature step
`n = √(n_eff² + n_puddle²)`, which happens *after* `C_q` has been evaluated.
So **the disorder regulariser cannot reach `C_q` at all** in the committed
form.

The converse is the sharper half. The induced charge is

```
n_eff = (C_ox·dV/e) · f(C_q),      f(C) = C/(C + C_ox) ∈ (0,1)
```

and `dV` appears as a **factor**. At `dV = 0` the product is zero whatever
`f` returns, because `f` is bounded. Measured on the shipped function with
four different `C_q` conventions substituted — the thermal rounding, the T = 0
dispersion, the Southampton disorder floor, and this thesis's own
`C_Q(n_puddle)` — `n_eff(0)` is **`0.0` bitwise in all four**, and
`carrier_density(V_dirac)` returns **exactly** `n_puddle`
(`5000000000000000.0`, difference `0.0`). So **the thermal regulariser cannot
reach the Dirac-point density at all.**

Both halves are checked against *exactly-known* values — `0.0` and
`n_puddle` — and not against a tolerance. That matters here for the reason
2026-09-17 established and §4.6.10 restated: a plausible-range check on these
quantities would have passed under either premise.

> **What the committed model does is apply exactly one regulariser to each of
> two different quantities. That is a reuse of one parameter for two jobs — a
> modelling choice to defend — and not a double count, which is an arithmetic
> error.**

Item 18 was created by comparing two *literature forms* without first asking
which quantity each term in **this** model reaches. The generalisable form:
a cross-model comparison establishes that two conventions exist; it does not
establish that one model contains both.

#### (b) The spelling that does double-count, and its size

There is one spelling in which both regularisers reach the same quantity.
Item 16's correction (form B_T, §4.6.9) makes `C_q` a function of `E_F(n)`;
place the floor **inside** that argument,

```
C_q = C_q( E_F( √(n² + n_puddle²) ) / e ),
```

and the capacitance now carries the thermal rounding and the disorder floor
at once. This is the most double-counting spelling available, and it is
arguably the *consistent* one: the channel's density of states should be
evaluated at the density the channel has.

Measured at the Dirac point, in form B_T:

| | thermal OFF (T = 0) | thermal ON |
|---|---|---|
| floor outside `C_q` | 0.000000000e+00 | 8.437398908e-03 |
| floor inside `C_q` | 1.942158579e-02 | **1.991221401e-02** |

| composition | value |
|---|---|
| both / max(alone) | **1.025262** |
| both / (sum of increments) | 0.714750 |
| thermal **on top of** the floor | **+2.5262 %** |
| floor **on top of** thermal | +135.9994 % |

The double count is **+2.5262 %** in `C_q(0)`. It is not a factor 2.

> **The "factor 2" was the ratio of the two candidate values. That ratio is
> the cost of *choosing wrongly* between the conventions. It is not the cost
> of *applying both*.**

Those are different quantities and only one of them was the item. Reading one
for the other is the 2026-10-06 backlog item — *census every criterion that is
a ratio* — in its **first instance where the ratio was not a criterion but an
estimate of a magnitude**. The item has been open five sessions as a question
about verdicts; it turns out to have been a question about magnitudes too.

And the reason the two do not add is worth stating as a rule, because it is a
property of regularisers and not of these two:

> **Two regularisers that enter through a monotone function of a quadrature
> sum compose max-like, not additively: the larger swallows the smaller.**

Here `max(8.44e-03, 1.94e-02) = 1.94e-02` and `both = 1.99e-02`, which is
2.53 % above pure max and 28.5 % *below* the additive null.

#### (c) What it costs downstream, at the model's own biases

Item 17's second requirement — run a two-knob check at the model's own
operating point, not at a reference that inflates it — applied here:

| `ΔV` | neither | thermal | disorder | both | both/max |
|---|---|---|---|---|---|
| 0 V | 5.000000e+15 | 5.000000e+15 | 5.000000e+15 | 5.000000e+15 | 1.000000 |
| 1.2 V | 5.730515e+15 | 5.732955e+15 | 5.741244e+15 | 5.741709e+15 | 1.000081 |
| 2.7 V | 8.085749e+15 | 8.087032e+15 | 8.095579e+15 | 8.096222e+15 | 1.000079 |

At the Dirac point **all four conventions agree bitwise**, because every one
of them is overwritten by the quadrature floor on the *result*. The quantity
item 18 was created to protect is the one quantity in the model that cannot
move. `σ_min` is identical to the last bit in all four
(`3.204353268e-04` S/sq): it is `n_puddle·e·μ` and nothing else, so **no
committed conductivity, resistance or on/off number in Chapters 4, 5 or 6 is
exposed to item 18 at neutrality.**

At the RF bias the whole double count is worth **+0.1527 %** in `n`, against
−1.0212 % for item 16 itself and −0.60 % for the exact charge relation. Item
18 is **real but small, and it is smaller at the operating point than at
neutrality** — the opposite of how it was ranked when it was created, and the
second time in two sessions that a cost estimate has been found to price
something other than the repair it named (§4.6.9's lesson, one level up).

#### (d) The floor's own spelling carries an exact factor 2

One finding sharpens the comparison rather than the cost, and it is the reason
yesterday's ratio should not be read as settled. The Southampton floor is

```
C_q,min = e²√n₀ / (√π·ħ·v_F)
```

while graphene's density of states, `D(E) = 2|E|/(π ħ² v_F²)`, carries the
spin–valley degeneracy `g = 4` and gives this thesis's own T = 0 capacitance

```
C_Q(n) = e² D(E_F) = 2e²√n / (√π·ħ·v_F) = 2e√n / A_F.
```

**The literature floor is exactly half this thesis's own `C_Q` at the same
density** — measured at **0.000 ulp from ½** with one consistent `ħ`. It is a
convention difference, not a numerical one. Whatever its origin — a
per-species `n₀`, a bilayer convention, a different definition of the residual
density — it is not a quantity that can be compared against this thesis's `kT`
before the factor is fixed.

The consequence for yesterday's number: the ratio **1.1509** is `kT` against
the *halved* spelling. On this thesis's own degeneracy convention it is
**2.3018**, and the claim "disorder and temperature regularise `C_q` equally
at this operating point" weakens from *the same size* to *the same order*.

**Yesterday's 1.1509 is not rewritten here, and §4.6.9(h) is annotated rather
than corrected.** Both figures are reported, because which one is right
depends on what `n₀` denotes in the source, and that is a question for the
literature rather than for arithmetic. It is now the only open item on this
list whose resolution requires reading rather than computing.

#### (e) The module's first draft failed its own substitutability check

Recorded because the check that caught it was written **before** the
measurement rather than during a neutrality pass afterwards, which is the
first time that ordering has held in this repository.

Draft 1 built the two switches on `graphene_diffusion_current_model`'s literal
`E_CHARGE` and `HBAR` while `quantum_capacitance()` uses `scipy`'s, and
reproduced the shipped `carrier_density(2.0)` only to **2375 ulp**. That is
the **eighth** instance of the 2026-10-09 rule and the second in which the
substituted quantity was a *physical constant*.

The repair is 2026-10-11's repair and not a better transcription: each switch
**substitutes one shipped name** — `gfet._thermal_rounding` for the thermal
regulariser, `gfet.n_puddle` for the disorder one — and leaves every other
operation in the shipped call path. Form A reached through the switches is
then the shipped function, and agrees with it at **0 ulp** by construction
rather than by care. One arrival counter per patched name, each checked at a
call that reads that name (§4.6.10's sixth repair).

One caveat is recorded in the output rather than smoothed over: the T → 0
validation reports a departure of *exactly* zero, which is weaker than
`1e-16` would be. At `η ~ 10⁷` the term `2·log1p(exp(−η))` underflows to
exactly `0.0`, so the two expressions reduce to the same floating-point
operations. The informative version of that check is at moderate `η`, where
§4.6.10 already adjudicates the identity against 60-digit arithmetic.

#### (f) The source was read to settle the factor, and it leaves `n₀` undefined — but its own composition law bounds the whole class of error

The one residue of (d) that arithmetic cannot reach is what `n₀` denotes in
the source. The paper was fetched and read for exactly that, and **it does not
say.** At zero surface potential it states only that "the channel has been
reported to have a charge density, `n₀`", cites two references for it, and
never calls it total, per-species, or a disorder parameter; in its resistance
section the same symbol is "the minimum charge density". **No degeneracy
factor appears anywhere in the paper, and it gives no density of states at
all** — it goes directly from `n` to `C_q`.

> **The comparison is blocked not by disagreement but by the source having
> left a symbol undefined and delegated it to a citation.** There is no stated
> convention to convert, so the factor 2 of (d) cannot be reconciled by
> reading harder; both 1.1509 and 2.3018 stay in the record and neither is
> installed.

Its fitted values are in the `10¹⁶ m⁻²` range — **1.2e16 m⁻²** for its Test A,
**2.4×** this thesis's `n_puddle`. That is a second literature anchor on the
standing `n_puddle` question, and it points the same way as §4.6.7's: this
thesis's floor may be low rather than high.

What the reading *does* settle is the composition law, and it settles it in
favour of (b). Their Eq. (3) is

```
C_q(φ_s) = (C_q,min² + 2h²) / √(C_q,min² + h²),      h = C_q,var/2
```

which returns `C_q,min` **bitwise** at `φ_s = 0`, tends to `C_q,var` to twelve
digits when that dominates, and is therefore a **smooth max and not a sum**.
Its excess over pure max is maximal exactly where the two terms are equal,
and there it is

```
sup  C_q / max(C_q,min, C_q,var)  =  1.5/√1.25  =  3/√5  =  1.341640786
```

measured at `C_var/C_min = 1.000000` over a 200 001-point sweep, against that
closed form.

Two consequences, both for item 18 and not for the model.

**The literature applies `n₀` twice itself.** Eq. (2) floors `C_q` with it and
its Drude relation `n = √(n₀² + n*²)` floors the density with it, with no
temperature anywhere. So the comparison item 18 rested on is not between a
model that double-counts and one that does not; it is between two
**allocations of regularisers to quantities**. This thesis spends `kT` on
`C_q` and `n₀` on `n`; the source spends `n₀` on both.

**And §4.6.9(h)'s intuition was right about where the worst case is and wrong
about how big it is.** A smooth-max composition is furthest from pure max
exactly when its two terms are *equal* — which is precisely the condition that
section identified as holding at this `n₀` and this `T`. That reasoning was
sound. But at that worst case the excess is `3/√5`, **+34.164 %**, and this
thesis sits at **+2.5262 %** because its two terms are 2.30:1 rather than 1:1.

> **A quadrature-composed regulariser pair cannot exceed the larger of its two
> terms by more than 34.2 %. The "factor 2" was unreachable in principle, not
> merely unmeasured.**

An upper bound on the entire class of error the item named, derived from the
source's own composition law rather than from this thesis's. It is quoted as
context for the item's estimate and **not** as a bound on this thesis's own
arithmetic: the nested spelling of (b) is a different function, and its
+2.5262 % is measured rather than bounded.

#### (g) What survives

Not a double count. What survives is the **reuse of one parameter for two
jobs**: `n₀` buys the density-of-states broadening in the literature and the
carrier floor here, and this thesis has never stated which of the two it means
`n_puddle` to be. That is the standing `n_puddle` question — §4.6.7 reached it
from the disorder floor being 2.07× Chapter 3's, and item 18 reaches it from
the opposite direction — and it is still deliberately not an item, because
changing `n_puddle` moves committed numbers in three chapters and the case for
changing it has not been made.

### 4.6.12 Item 16's blast radius, measured: the contact model was outside the harness, Chapters 5 and 6 were never inside it, and Pd's residual does not survive

Section 7.9 item 16 — "flip `quantum_capacitance()`'s argument from the gate
overdrive to `E_F/e`" — has been ranked first since 2026-10-12, carrying the
note that flipping it "moves every committed number in Chapters 4, 5 and 6 and
needs its own session". This is that session, and it does not flip the
default. It measures the blast radius first, and the blast radius is not what
that sentence says. It is **smaller across the thesis and larger inside
Chapter 4**, and inside Chapter 4 it reaches a verdict.

`graphene_contact_argument_blast_radius.py`, 24/24,
`contact_argument_blast_radius_output.txt`,
`contact_argument_blast_radius.png`.

#### (a) The harness that sized item 16 could not reach half of Chapter 4

`graphene_quantum_capacitance_argument.on_consistent_argument()` rebinds three
names: `gfet.carrier_density`, `gfet.quantum_capacitance`, and
`crx.carrier_density` — the third added 2026-10-10 after exactly this class of
fault. It does not rebind `graphene_contact_doping_model.quantum_capacitance`,
which that module bound at import time on its line 50

```python
from graphene_fet_model import quantum_capacitance, sheet_conductivity, ...
```

and which `contact_edge_carrier_density()` calls on its line 105 with

```python
C_q = quantum_capacitance(dV, T=T)        # dV = -(W_metal - W_graphene)
```

— **exactly the defect item 16 names**, in a second module, against a
work-function-induced potential instead of a gate overdrive.

Measured rather than read from the source: under `on_consistent_argument()`
that binding is reached **11 times**, and **exactly zero** of those calls
arrive at the corrected function. Every cost ever quoted for item 16 —
−2.2659 % in the raw charge, −1.0212 % in `n` at 3.5 V, −0.2530 % in `I_d`
(Section 4.6.9(e)) — was therefore measured with Sections 4.5, 4.7 and 4.8
**switched out of the circuit**.

> **Nothing errored and no number moved, which is the whole difficulty: a cost
> table built through this harness reports a model that was absent as a model
> that is insensitive.**

#### (b) And Chapters 5 and 6 were never inside it either

The same counters, on nine committed summary routines —
`graphene_interconnect_model.summary_numbers()` and
`self_check_parallel_limit()`, the collection model's `summary_numbers()` and
`modeled_eqe_by_metal()`, the signed-carrier model's three tables, the
two-contact model's two — record **0 charge-path calls between them**.
Chapter 5 imports nothing from the charge model. Chapter 6 imports `mu`,
`METAL_WORK_FUNCTIONS`, `W_GRAPHENE` and `lambda_decay`: parameters, not the
charge path.

> **Item 16 is a Chapter 4 operation that reaches two models, not a
> three-chapter one that reaches one.**

#### (c) The correction is 4.7× to 12.2× larger at the contact

The contact geometry is where the argument error is least forgiving, for two
compounding reasons.

`C_interface = ε₀/3.5 Å = 2.529768e-02 F m⁻²` is **65.93×** `C_ox`, and `C_q`
is *comparable* to it rather than small against it: the series factor
`C_q/(C_q + C_s)` runs **0.555–0.806** at the contact against **0.999** in the
channel at `V_g = 3.5 V`. The channel's near-unity factor is what attenuated
item 16 by 0.451 (Section 4.6.9, §7); the contact has no such attenuation.
And the potential is a work-function difference of 0.15–1.15 V rather than a
few hundred mV of overdrive, so the gap between `ΔV` and `E_F(n)/e` is wider
to begin with — `E_F/e` lands at 0.315 V where `ΔV` is 0.620 V for Pd.

| Metal | W_F (eV) | ΔV (V) | n_contact, form A (m⁻²) | form B_T (m⁻²) | change | E_F/e (V) |
|---|---|---|---|---|---|---|
| Ti | 4.33 | +0.170 | +1.644935e16 | +1.541997e16 | **−6.2578 %** | 0.1449 |
| Cr | 4.50 | −0.000 | −0.000000e+00 | −0.000000e+00 | — (ΔW_F = 0 exactly) | 0.0000 |
| Cu | 4.65 | −0.150 | −1.380517e16 | −1.314672e16 | **−4.7696 %** | 0.1338 |
| Ni | 5.04 | −0.540 | −7.111301e16 | −6.227844e16 | **−12.4233 %** | 0.2911 |
| Au | 5.10 | −0.600 | −8.034794e16 | −7.031424e16 | **−12.4878 %** | 0.3094 |
| Pd | 5.12 | −0.620 | −8.343500e16 | −7.300861e16 | **−12.4964 %** | 0.3152 |
| Pt | 5.65 | −1.150 | −1.660634e17 | −1.463454e17 | **−11.8738 %** | 0.4463 |

Against **−1.0212 %** for the same correction in the channel.

#### (d) And it changes a verdict — the one Section 4.10 said was fragile

Section 4.7's Table reports `R_transmission = Rc,measured − R_extra` for four
metals and finds three negative, i.e. unphysical, with **only Pd positive** at
+50.949 Ω·µm. Section 4.10's condition-number audit then found that this one
positive row — the only evidence the additive decomposition might survive —
**flips sign on a 9.6 % error in `R_extra`**, while Cu's would need 81.1 %.

The argument correction moves Pd's `R_extra` by **+47.6925 %**, five times
that margin.

| Metal | Rc, measured | R_extra, form A | R_extra, form B_T | change | R_tr, A | R_tr, B_T |
|---|---|---|---|---|---|---|
| Cu | 184.0 | 971.461 | 959.236 | −1.2583 % | −787.461 | −775.236 |
| Ni | 470.0 | 830.462 | 1003.490 | +20.8352 % | −360.462 | −533.490 |
| Au | 519.0 | 609.159 | 848.120 | +39.2281 % | −90.159 | −329.120 |
| Pd | 584.0 | 533.051 | 787.276 | **+47.6925 %** | **+50.949** | **−203.276** |

> **Pd's residual does not survive the correction. Under form B_T the additive
> decomposition fails on 4 of 4 metals rather than 3 of 4, and Section 4.10's
> 9.6 % margin — computed as a tolerance to INPUT error — is consumed by a
> MODELLING error that was already in the committed code.**

Two things that do *not* move, and both were checked rather than assumed: the
metal ordering by residual (Cu < Ni < Au < Pd under both forms) and the
n-type/p-type assignment of Section 4.5 for all seven metals.

Note the direction. `R_extra` goes **up** because `n_contact` goes **down**:
a less strongly doped contact edge is a more resistive one, and the
doping-profile integral of Section 4.5 integrates `1/σ(n(x))`. Cu is the
exception at −1.2583 % because Cu's contact density (−1.38e16 m⁻²) is below
the on-state bulk reference (2.0e16 m⁻²), so reducing it moves the profile
*toward* the bulk over part of the junction.

#### (e) Why the default is still not flipped

Nothing here is an argument for keeping form A. It is an argument that the
flip is a larger operation than the item said, and that two of its
preconditions were discovered today rather than before:

1. `contact_edge_carrier_density()` must be flipped **in the same commit** as
   `carrier_density()`, or Chapter 4 will hold a corrected channel against an
   uncorrected contact and Section 4.7 will mix the two.
2. Section 4.7's text and Section 4.10's annotation both have to be rewritten,
   not re-run: the prose reasons explicitly from "only Pd gives a small,
   plausible positive residual", and that sentence becomes false.

#### (f) What this section does not claim

It does not claim form B_T is right. `E_F/e` is the argument the dispersion
requires and the argument two independent literature lineages use (Sections
4.6.9(h) and 4.6.11(f)), but the contact model's *other* approximations —
a vacuum-like 3.5 Å gap with no interface dipole layer, a single
`lambda_decay` for every metal — are untouched here and are plausibly larger
than 12 %. What is claimed is narrower and is measured: **the two forms
differ by 4.7–12.5 % at the contact, the harness that sized item 16 could not
see that difference, and the difference is large enough to move a verdict
Section 4.10 had already flagged as the fragile one.**

It also does not claim the Pd flip is the *end* of Table 4.7's story. Under
form B_T all four rows are negative, which strengthens rather than weakens
Section 4.7's stated explanation — that TLM-extracted `Rc,measured` already
absorbs part of the same near-contact region this model integrates
separately. A four-of-four failure is more consistent with a double count
than a three-of-four one was.

#### (g) Two instrument defects found on the way, both kept in place

**(i) This section's own first draft counted the wrong arrival.** It counted
calls to the *name* `cd.quantum_capacitance` and asserted they would be zero
under the patch. They are 11, whether or not anything was substituted into
that name. The question is not whether the binding is reached but which
function is behind it, and the counter has to live inside the *replacement*.
**Ninth** instance of the 2026-10-10 rule, and the second in its "counts a
different arrival from the one it gates" direction.

> **A control that counts calls to a name cannot see whether the right
> function is bound to it.**

**(ii) Item 13's overflow is closed in the shipped function and open in the
harness.** `qca.cq_thermal` is still the committed `ln(2(1+cosh η))` spelling,
deliberately retained there, and still overflows at `η = 710.4758`. A fixed
point built on it takes the *electrostatic* density as its first iterate, and
that iterate grows with the series capacitance — so the solver returns `nan`
for any `C_s` above ≈ 100 F m⁻². `C_interface` is 2.53e-02 and the worst
first-iterate `η` in this session is **19.2**, so no number here is affected.
But it is item 13's lesson a third time: *a defect unreachable from the
device's operating range is reachable from a check on the device's model.*
Form B_T is built here on the **shipped** `gfet.quantum_capacitance`
accordingly — which also removes a **1.2254e-09** constants gap (prefactor
ratio `qca/gfet = 1.000000001225439`) from the A→B_T comparison, so the two
forms differ in the argument and in nothing else by construction.

**(iii) And one tolerance chosen before the run failed.** The
large-`C_s` limit check demanded `C_s`-independence to `1e-9` and measured
`5.254e-04`. The limit is approached at **first order** in `C_q/C_s`, so no
fixed tolerance states what it claims. Replaced by the rate — the departure
must fall by ten per decade of `C_s`, measured **9.9933** and **9.9993**.
**Item 15's eighth worked instance**, and the repair is 2026-10-11's: let the
structure set the criterion.

### 4.6.13 The rebind set is a runtime object: 25 stale sites, 1 that can fire, and a repair that blinded its own audit

§4.6.12 deferred item 16's flip for a stated reason — nobody knew what a
rebinding of the charge path reaches — and closed by noting that the
`cd.quantum_capacitance` miss found there had been found **by hand**, as had
the two before it. Three hand-found instances in ten days is a rate. This
section measures the exposure that rate implies, and then removes the
mechanism that produced it.

**(a) What a hand-maintained list leaves behind.** An AST pass over all 53
modules finds **149** module-level `from`-imports of repository names across
**73** distinct `(owner, name)` pairs. Against the repository's **14**
hand-listed rebinding sites (`alias.name = ...`), **25** of those bindings are
left stale. Three had been found by hand; **22 had not.** The identity classes
that matter are small and the two that item 16 needs are the smallest:

| name | sites a single `setattr` fails to reach |
|---|---|
| `gfet.carrier_density` | 1 (`crx.carrier_density`) |
| `gfet.quantum_capacitance` | 1 (`cd.quantum_capacitance`) |
| `gfet.mu` | 4 (`cd`, and Chapter 6's three photodetector modules) |
| `gfet.W` | 2 (`crx`, `cd`) |
| `gfet.n_puddle` | 2 (`cd.n_puddle`, `dr.N0_SHIPPED`) |
| `gfet.C_ox` | 2 (`cd.C_ox`, `dr.C_OX`) |
| `gfet.sheet_conductivity` | 2 (`crx`, `cd`) |

**(b) Which of them are USED, measured without entering any manager.**
Perturb **only** the holder's own copy and see whether the holder's own
numbers move. This is the decomposition that matters, and it needs no
manager — which is what makes it applicable to the ten inline `try/finally`
rebinding sites that cannot be entered at all.

| stale binding | max relative change in its holder's numbers | verdict |
|---|---|---|
| `cd.quantum_capacitance` | 9.8781 (988 %) | USED |
| `cd.n_puddle` | 35.764 % | USED |
| `crx.carrier_density` | 27.007 % | USED |
| `crx.W` | 27.007 % | USED |
| `pcm.mu` | 9.1743 % | USED |
| `psc.mu` | 5.3067 % | USED |
| `ptc.mu` | 3.9039 % | USED |
| `cd.mu`, `cd.C_ox`, `cd.W` | exactly 0 | imported, never reaches a number |

**(c) The repair, and why it is a computed class and not a longer list.** A
list is exactly what was wrong all three times, so the mechanism had to stop
being one. `from M import n` copies a **reference** at import time, so the set
of places to rebind is the set of names that at runtime hold the same
**object** — discoverable by identity, with no import graph, no AST pass and
no list. `graphene_binding_class.rebind()` computes it at entry and is now
installed at all three charge-path sites (`on_consistent_argument`,
`on_exact`, and `graphene_potential_declaration.main`, which had no second
line at all). Behaviour is unchanged, and that is **verified rather than
asserted**: all three modules' committed transcripts regenerate
**byte-identically**, and each was confirmed byte-reproducible *before* the
change so that the diff is a control and not a hope.

**(d) The class is necessary and not sufficient — item 16 needs two
replacements, not one rebinding.** `cd.quantum_capacitance` is deliberately
**not** rebound. It is reached at a contact where the series capacitance is
`C_interface = 2.529768e-02 F m⁻²`, **65.93×** `C_ox`, and the series factor
runs **0.555–0.806** against **0.999** in the channel (§4.6.12). The object
the owner holds after item 16's flip solves for `n` against `C_ox`; copying it
to the contact installs a **channel solver at the contact**, which is a
different wrong answer rather than a repair. The contact-correct form is the
one §4.6.12 measured separately (**−12.4964 %** in `n` for Pd). The exclusion
lives once, in code, with its reason attached, and `rebind()` **refuses an
exclusion without a stated reason** — because an unstated exclusion is how all
three 2026-10 faults were introduced.

**(e) And the static count over-reports, which retires a claim made earlier
the same day.** A miss needs the stale holder to be **imported in the patching
process**, and an import graph over the repository cannot see that. Measured
per rebinding site in a **fresh interpreter**, so that the audit's own imports
cannot supply the answer: of the 22 remaining static misses the holder is
**PRESENT for 1** and **ABSENT for 21**. `import
graphene_velocity_saturation_model` brings in three repository modules and
none of the four `mu` holders is among them. `graphene_potential_declaration`
does not import `contact_resistance_crossover` at all — four repository
modules enter its process and `crx` is not one of them — so the reading of
that site as "the 2026-10-10 fault unrepaired in a third module" was **wrong
about this repository while right about the shape**: there was no stale
binding to leave behind and no wrong number. The single survivor is
`blast_radius.series_capacitance` leaving `cd.C_ox` stale, and (b) already
measures `cd.C_ox` as UNUSED, so it is cosmetic too.

> **The 25 is retained in place as an upper bound rather than restated. The
> honest standing is three independent axes — PRESENT, USED, REACHED — and
> nothing currently in the repository clears all three. What the sweep found
> is a landmine count, not a defect count.**

**(f) Four instrument faults, all four kept in place.**

**(i) A repair that blinds its own audit is not a repair** — §7.5.4's "created
by a fix" class in its purest form yet. Retrofitting `rebind()` made all three
charge-path sites **invisible** to the sweep, whose only notion of a rebinding
site was `alias.name = ...`, and the validation that exists to keep
2026-10-13's finding on the record **failed while nothing was wrong**: the
manager it names had stopped being a rebinding site in the only sense the
instrument understood. The sweep now recognises both forms and reports a
class-complete site as covered by construction with its exclusions listed. The
validation is rewritten to assert the same fact in the form the repair left it
in — zero *silent* misses plus one exclusion named with a reason of substance —
and its original no-silent-miss clause is **kept**, not replaced.

**(ii) The planted-binding control failed on the instrument.** The module
scanner skipped every file whose name begins with an underscore, so the sweep
reported **zero** bindings in a file it had itself just created with one. An
underscore-prefixed module is as importable as any other; the filter was a
blind spot and not a convention. **Fourth** occurrence in this repository of a
control failing on the instrument rather than on the thing measured.

**(iii) The identity control found a blind spot in the AST rule, and the rule
was extended rather than the control dropped.** Twelve second-order
**renaming** aliases (`dr.C_OX is gfet.C_ox`, `dr.N0_SHIPPED is
gfet.n_puddle`, `psc.W_CROSS_VACUUM is cd.W_GRAPHENE`) are invisible to any
rule that reads import statements. Each is a real site. What identity yields
is not a directed edge but the **equivalence class** of names sharing one
object, and for the question "what must a patch rebind" the class is the safe
criterion and the import graph the unsafe one.

**(iv) The arrival counter carries a positive control, so that LATENT is a
measurement.** Under their own workloads both charge-path managers reach `cd`
**zero** times, which agrees with §4.6.12's finding that the contact model
sits outside the harness — but a counter that reads zero because it is broken
reads the same. The same counter reads **exactly 7** when `cd` is exercised
over its seven metals.

**(g) Validations against exactly known values, six, all passing.** The
2026-10-10 repair must show 0 misses (0, and both managers now cover it by
construction); 2026-10-13's finding must survive as a reasoned exclusion with
no silent miss (it does); a planted binding must be found exactly once by each
rule (1, 1, and 2 by identity — the import and the rename); identity-patch
neutrality over **98** wrapped sites at **0 ulp**; the usedness instrument's
own null, a perturbation by a factor of exactly 1.0, at **exactly 0**; and the
arrival counter's positive control at **7 of 7**.

What stays open is named and small: `graphene_contact_doping_model` needs its
own contact-correct replacement before item 16 can flip, and the eleven
surviving hand-listed sites outside the charge path are **not** converted,
because each is latent by import today and converting one changes nothing
until it is not.

### 4.6.14 The contact gets a replacement instead of an exclusion, and the identity class turns out to span two operating points

Section 7.9 item 16 — flip `quantum_capacitance()`'s argument from the gate
overdrive to `E_F/e` — has been the first-ranked open item on this chapter
since 2026-10-10. §4.6.13 reduced its blockers to exactly two, of which the
first was:

> `graphene_contact_doping_model` needs its **own** contact-correct
> replacement, because the contact cannot take the channel's object. The flip
> is two replacements, not one rebinding.

**That blocker is now discharged.** What follows is the discharge, and the
structural finding that came out of attempting it. Neither flips the default:
item 16 remains open, and what has been removed is the reason it could not be
flipped.

#### (a) A measurement harness is not a replacement

The contact-correct numbers were not unknown. §4.6.12 measured them a
fortnight before, on 2026-10-13: the correction at the contact edge runs
**−4.7696 %** (Cu) to **−12.4964 %** (Pd) in `n`, against **−1.0212 %** for
the same correction in the channel at `V_g = 3.5 V`. The contact is **4.7× to
12.2× less forgiving**, because its series factor `C_q/(C_q + C_s)` runs
**0.555–0.806** where the channel's is **0.999**, and it has a wider `dV` to
be wrong about.

Having the numbers was not having the object, for three reasons, and only the
third was a surprise:

1. **The object belonged to the audit.** It lived inside
   `graphene_contact_argument_blast_radius`, so nothing outside that module
   could install it. Item 16's flip had no second object to install even
   though its size was known to four decimal places.
2. **The harness rebound by hand.** A single assignment,
   `cd.quantum_capacitance = ...` — precisely the mechanism §4.6.13 retired
   for the channel after measuring that hand-listed rebinding leaves **25**
   stale sites at **14** listed sites, **22** of which no audit had ever
   found. The contact half was still on the retired mechanism.
3. **And the replacement mechanism could not express what the contact needed.**

#### (b) Identity answers *which names move together* and does not answer *to what*

`bc.rebind()` — §4.6.13's repair — computes a name's runtime **identity
class**, the set of names that hold the same object, and gives every member
**the same replacement**. Its one escape hatch is `CONTACT_OPERATING_POINT`:
leave this member alone, with a reason recorded. So its vocabulary is
`{same replacement, left alone}`, and it has no word for *rebind this member
to a **different** object* — which is exactly what the contact requires the
moment it has an object of its own.

Measured at call time, the class is:

| member | operating point | series capacitance | series factor |
|---|---|---|---|
| `graphene_fet_model.quantum_capacitance` | channel | `C_ox` | 0.999 |
| `graphene_contact_doping_model.quantum_capacitance` | contact | `C_interface` = 65.93 × `C_ox` | 0.555–0.806 |

> **Two names can hold one object and still sit at two operating points. An
> operating point is not something identity can see, so the class must be
> PARTITIONED by one — and the partition has to come from the caller.**

This is a limitation of §4.6.13's repair discovered by trying to use it, not a
defect in it. The repair's claim — that the set of names to rebind is
computable and needs no list — survives intact; what fails is the implicit
assumption that one class takes one replacement.

`bc.rebind_by_operating_point()` is that generalisation. `bc.rebind()` is
**left untouched**, and all three transcripts it owns regenerate
byte-identically. The new function takes a map from class member to
`(replacement, reason)`, where `replacement = None` expresses §4.6.13's
exclusion in the same vocabulary rather than as a separate mechanism, and the
`reason` is mandatory in both cases.

The refusal carries across in the form the same lesson takes here. §4.6.13
refuses an **exclusion without a stated reason**, because an unstated
exclusion is how all three of the 2026-10 delivery faults were introduced.
This refuses an **unassigned class member**:

> **Silently leaving a member at the shipped object is the identical fault,
> reached by omission from a map instead of omission from a list. The whole
> value of computing the class is that the omission is detectable, so it is
> detected: the manager raises and rebinds nothing.**

#### (c) The emptiness that discharges the precondition

`on_consistent_argument_everywhere()` installs the charge relation and **both**
capacitance operating points together, and its left-alone set is **empty** —
four sites installed, zero left at the shipped object:

```
graphene_fet_model.carrier_density                     channel
contact_resistance_crossover.carrier_density           channel
graphene_fet_model.quantum_capacitance                 channel  (C_ox)
graphene_contact_doping_model.quantum_capacitance      contact  (C_interface)
```

The emptiness is **asserted by the mechanism, not by this section**: the
manager refuses to open if any member of either class is unassigned, and a
negative control plants a missing member and requires the refusal to fire.
Without that control, an empty left-alone set and a mechanism that never
checks would read identically — the 2026-10-14 rule that a null result is a
measurement only when the instrument carries a control.

#### (d) Validations against exactly known values

Nine checks, 9/9, five of them against values known exactly rather than
plausible ranges:

| check | value required | result |
|---|---|---|
| `ΔWF = 0` → `n = 0` under **both** forms | exactly 0 | `−0.0` / `−0.0` |
| form `B_T` odd in `dV` | **0 ulp** | 0 ulp over 5 offsets |
| new machinery fed the **OLD** argument vs shipped | bitwise | **0 ulp**, all 7 metals |
| quantum limit `n = 4(dV/A_F)²`, error order in `1/C_s` | ratio exactly 100 | **100.000** at `C_s × 100` |
| today's owned object vs §4.6.12's independent construction | bitwise | **0 ulp**, all 7 metals |
| planted unassigned member | must raise | raises |
| three arrival counters | positive control on each | 1 / 19 / 1 |

The third is the one that makes the rest mean anything. Feeding the new code
path a `C_q` that ignores `n` and returns the shipped function at the **old**
argument reproduces the shipped `contact_edge_carrier_density()` **bitwise**
over all seven metals. That is what establishes that form A and form `B_T` at
the contact differ in `C_q`'s **argument** and in nothing else — no second
solver, no second set of constants, no altered series factor. Without it,
every per-metal number in §4.6.12 and here could be the residue of some other
difference between two implementations. The fourth is the honest form of a
closed-form check: at finite `C_s` the solve differs from the limit by
`O(C_q/C_s)`, so a single threshold would be a plausible-range check wearing a
closed form's clothes, and the exact statement is the **convergence order**.

#### (e) And the first draft of one check misattributed a real defect

The oddness check (second row) was first written through the shipped
`contact_edge_carrier_density()` rather than on the solver, and **failed at 3
ulp**. It then failed at **2–3 ulp under form A as well**, which is what
identified the cause:

> **The shipped function's own `dWF = work_function_metal - W_GRAPHENE` is not
> antisymmetric in double precision.** `(4.5 + 0.05) − 4.5 =
> 0.050000000000000266`, while `(4.5 − 0.05) − 4.5 = −0.04999999999999982`.
> The solver is exactly odd; the work-function subtraction is not.

So a pre-existing 2–3 ulp asymmetry in the shipped contact model was briefly
read as a defect introduced by the new solver. It is kept as its own
validation row rather than quietly dropped, because this is §4.6.6's rule for
the second time in two sessions: **a check that reads two knobs at once
attributes the effect to the wrong one.** The defect is cosmetic at 2–3 ulp
against a 12 % physical effect, and it is now on the record instead of being
unknown.

#### (f) What this does and does not settle

| | standing |
|---|---|
| Precondition (i) — the contact's own object | **DISCHARGED** |
| Precondition (ii) — §4.7's prose and §4.10's annotation | **OPEN**; writing, not plumbing |
| Item 16 itself — flipping the default | **OPEN**; nothing here flips it |

§4.7's prose becomes false under the flip — "only Pd gives a small, plausible
positive residual" cannot survive a 12.4964 % correction to Pd's `n_contact`
— and §4.10's conditioning annotation inherits that. Both are rewritten under
item 16 proper, not here, and §4.6.12's numbers stay annotated in place rather
than replaced.

Module: `graphene_contact_consistent_argument.py`. Machine output:
`contact_consistent_argument_output.txt`. Figure:
`contact_consistent_argument.png`. Mechanism:
`graphene_binding_class.rebind_by_operating_point()`.


### 4.6.15 The inert-knob census was counting a syntactic pattern for a semantic fault, and three instruments had been dark for twelve days

Sections 4.6.13 and 4.6.14 built the machinery item 16's flip needs: an
identity class computed at run time, and a partition of that class by
operating point. Both deliver a change the same way everything in this
repository does — by **rebinding a name**. That makes one question
load-bearing for every result in Chapters 4 to 6 that was ever established
by mutation: *which readers does a rebinding actually reach?*

Since 2026-09-28 the known answer has had one documented exception, the
**inert knob**:

```
K = 1.0
def f(x=K): ...          # x's default is bound when `def` executes
```

A rebinding of `K` never reaches `f`, so a mutation audit that rebinds `K`
and watches `f` measures nothing and reports a pass. Three separate
instruments census the pattern —
`graphene_figure_provenance_audit.census()` (Section C),
`graphene_mutation_arrival_probe.repo_wide_frozen_census()` (Section 7b),
and `graphene_cross_module_delivery_audit` for the cross-module subset.

#### (a) All three transcripts were written on the same day, and none of them has been read since

The committed transcripts read **0**, **0** (an empty census behind a
`[PASS]` on "no module-level constant is captured as a default anywhere")
and **8**. All three were last written **2026-10-04**. Run today they read
**21**, **20** and **16**.

> **A fault class that 2026-10-04 drove to zero regrew over twelve days in
> front of three censuses, none of which was being read — and the fourth
> instrument whose only purpose is to notice exactly that
> (`graphene_pristine_transcript_audit.py`, "A COMMITTED TRANSCRIPT IS A
> CLAIM ABOUT THE CODE BESIDE IT") was last run on the same day.** Asked
> today, unprompted, it reported four of its subjects STALE or ERROR.

This is §7.7's unequal-scrutiny habit one level out. The repository has
spent three weeks learning that a guard which is not run is not a guard;
what it had not noticed is that the rule applies to the guard that checks
the guards.

#### (b) 21 against 20 is not a rounding difference

Over one path set both censuses name the **same five files**, and **four of
the five agree exactly**. They differ by one occurrence, in
`graphene_velocity_saturation_model.py` at line 589 —
`_R_harmonic(V_g, V_hi, n_quad=N_QUAD)`, defined **inside the body of
another function**. Measured rather than asserted, in all three
configurations, against exactly known values:

| where the `def` is | when the `def` statement executes | captured default |
|---|---|---|
| module body | once, at import | **FROZEN** |
| class body | once, at import | **FROZEN** |
| function body | on **every call** of that function | **LIVE** |

With `K` rebound 1.0 → 99.0: the module-level arm returns exactly 1.0, the
class-body arm exactly 1.0, the function-nested arm exactly 99.0.

> **THE SYNTACTIC PATTERN `def f(x=K)` IS TWO CLASSES WITH OPPOSITE
> BEHAVIOUR UNDER THE ONE OPERATION THIS REPOSITORY USES TO DELIVER
> CHANGES. A FUNCTION-NESTED CAPTURE IS NOT AN INERT KNOB. IT IS LIVE.**

So neither census is simply right:

- **figure_provenance's 21** counts the syntax correctly and over-counts the
  **fault** by exactly the live captures. Its own historical figures — "49
  occurrences across 14 files, 14 of them inside the audit modules … 37
  remain, and the 13 in audit modules are the next run's top item" — are
  syntactic counts of a class that is not uniformly the fault, and have
  never been split into the two.
- **mutation_arrival_probe's 20** counts the fault correctly **for the wrong
  reason**. Its `classify()` confirmed each statically-found site with
  `getattr(module, fname)` — added correctly, to drop a function a later
  redefinition or decorator had replaced — and a function-nested function is
  not a module attribute, so the site was dropped **silently**: no row, no
  class, no reason. The right verdict by an accident of resolution.

The accident does not generalise. The same line also drops a **class-body**
capture, which is *frozen* — so the next capture of that shape would have
been **undercounted**, and undercounting is the dangerous direction here: an
instrument that silently omits a true site has the identical signature to the
fault being censused. The repository contains no class-body capture today;
that is luck, not a property of the mechanism.

#### (c) Reconciled to zero residue

```
21 syntactic  =  20 FROZEN (the genuine inert knobs)  +  1 LIVE
```

Four checks, because a reconciliation that leaves a remainder means a third
cause is present and unaccounted for: the arbiter reproduces
figure_provenance's total exactly (R1), mutation_arrival_probe's confirmed
total exactly (R2), accounts for the whole gap with nothing left over (R3),
and agrees **per `(file, constant)`** rather than only in total (R4 — added
because two totals can agree while the rows behind them disagree in
compensating directions; it is the only check that catches a single dropped
row).

Per file: `graphene_velocity_saturation_model` 9 = 8 + 1,
`graphene_diffusion_current_model` 5 = 5 + 0,
`graphene_exact_charge_rederivation` 3 = 3 + 0,
`graphene_perfect_contact_counterfactual` 2 = 2 + 0,
`graphene_second_anchor_bound` 2 = 2 + 0 — the last of these added by the
session immediately before this one.

#### (d) The repairs, and the one deliberately not made

`classify()` now takes the verdict from the **scope** instead of discarding
what it cannot resolve: function-enclosed → `live_by_scope`, otherwise →
`frozen_unresolved`, which still counts as frozen because the static evidence
stands and `getattr` merely cannot reach it. The fix is a verdict from the
scope and **not a widening of the resolver**, because resolvability was never
the question being asked. The amendment is **count-neutral** — the frozen
total is 20 before and after, and the suite's own output differs in exactly
the two lines that now carry `1 rescoped-live`. The `getattr` confirmation's
original purpose is preserved, not removed: a module attribute rebound after
import genuinely no longer carries the captured default (V6).

The pristine audit needed its own repair, created by this session's decision
rather than inherited. Every suite here ends `return 0 if failed == 0 else
1`, so **a red guard is a nonzero exit by design** — and `run_suite`
discarded stdout on any nonzero exit and reported `ERROR` without comparing
the transcript at all. Tolerable while every suite was green; not tolerable
once this repository began committing red transcripts in preference to stale
green ones, since two suites would then have had transcripts nothing checks,
permanently. ERROR is now reserved for a run producing no usable stdout, and
a new check requires a transcript's own `TOTAL: N passed, M failed` line to
agree with its suite's exit code — the two being independent claims about one
run, which contradicted each other for the whole of the twelve days.

**Not done here, and the reason matters more than the omission.** The 20
genuine inert knobs are not converted. Converting 21 captures as one class
would have rewritten a live capture into a `None`-sentinel for no reason and
recorded it as a fix; the work is ranked first for the next run and is now
**specified correctly for the first time** — 20, not 21, with the live one
named so it is not converted by mistake.

#### (e) A defect in this section's own instrument, kept on the record

The first draft of `graphene_frozen_default_semantics.py` classified by
**nesting depth**: depth 0 frozen, depth ≥ 1 live. Its own control planted a
method in a class body, called it depth 1, and asserted *in prose* that it
was therefore not an inert knob. A class body executes once at import, so
that is wrong — and the repository contains no class-body capture, so the
draft's count would have been right **by luck**, which is the same luck
`classify()` was running on. That is what made it worth recording rather than
quietly fixing. Controls C3/C4/C5 now pin all three cases, including a class
defined *inside* a function, where a depth rule and a class rule both get it
wrong.

#### (f) What this means for item 16

The flip's blast radius is established by regenerating committed transcripts
and reading what moved. The instrument of record for that is the pristine
transcript audit — which, at the start of this session, reported four of its
own subjects STALE or ERROR.

> **A flip measured against that baseline could not have distinguished its own
> blast radius from twelve days of accumulated drift. That is a precondition
> (iii) for item 16 that no entry had listed, and it is discharged here: 13
> IDENTICAL, 0 STALE, 0 NUMERIC_DRIFT, 0 ERROR.**

Module: `graphene_frozen_default_semantics.py` (**24/24**). Machine output:
`frozen_default_semantics_output.txt`. Amended:
`graphene_mutation_arrival_probe.py`,
`graphene_pristine_transcript_audit.py`. Mutation-tested outside the
repository, so the instrument and the measurement do not move together: the
depth rule restored fails C1–C4; every capture called frozen fails C1, C2,
C5, R2, R3, R4; one frozen row dropped fails only R2 and R4; `run_suite`
reverted fails C7; and the 2026-10-04 transcript restored verbatim over
today's red guard fails the new consistency check **by name**.



### 4.6.16 The list that ranks this work had never been audited, and two of its items asserted the opposite of the tree

Every session in this project chooses what to do from the "Not yet covered
(candidates for future runs)" list at the foot of the preceding
`AUTOMATION_LOG.md` entry. §4.6.12 audited a harness, §4.6.13 a rebind set,
§4.6.15 three censuses. Nothing had ever audited **the list**, and a claim
nothing checks is the standing subject of this chapter's instrument
sections — §4.6.15's finding was that three instruments had been dark for
twelve days, and the list is the artefact with the widest reach of all of
them, because it is upstream of which instrument gets read at all.

Two of the 2026-10-17 list's twenty-five items assert a state of this tree
that the tree denies.

**A2, item 16's precondition (ii).** Recorded 2026-10-17 as "§4.7's prose and
§4.10's annotation — is still untouched and is still writing, not plumbing",
and recorded in §7.9 on 2026-10-14 as "the only precondition still open".
Commit `7604c05` of **2026-10-13** is titled *"Section 4.6.12, and the two
sentences in Sections 4.7 and 4.10 that the measurement makes false, annotated
in place"*. Both annotations are in this chapter, above §4.7's table and below
it, carrying −1.2583 % (Cu), +20.8352 % (Ni), +39.2281 % (Au), +47.6925 %
(Pd), and `+51.0 → −203.276`. The §7.9 sentence was therefore written **one
day after** the work it describes as outstanding, and carried for four days
after that.

What is genuinely open is narrower than the item said and **cannot be done
ahead of the flip**, which is the part worth recording. Both annotations are
written in the subjunctive: §4.7's says the sentence it annotates is "true of
the committed model and false of the corrected one" and leaves it standing
"pending the item 16 flip". The writing that remains is the writing only the
flip makes possible — the point at which "pending the item 16 flip" is itself
stale, the committed table becomes the superseded column rather than the live
one, and §4.10's margin-against-input-error becomes a statement about a margin
already consumed. A session that believes (ii) is missing can discharge it by
writing a sentence that asserts a state of the code that does not exist: the
§7.5.4 "created by a fix" class, in prose. **Item 16's remaining content is
the flip and nothing else.**

**A1, and it is a revival rather than a staleness.** The 2026-10-17 list
carries *"CHAPTERS 2 AND 3 ARE STILL UNDRAFTED — their computational results
have been complete for weeks. This is the largest remaining block of pure
writing and has now gone three sessions without being touched"*, and closes it
*"Worth ranking against item 16 deliberately rather than by default"*. Chapter
2 was drafted 2026-09-28 (`5bbd5f3`, eleven sections) and Chapter 3 on
2026-09-29 (`c1d5abc`, ten sections); the tree holds 375 and 330 non-blank
lines and 11 and 9 `##` sections. So the item offered as a competitor for
first rank was describing writing committed three weeks earlier.

And the list had already got this right. Measured across every dated entry:

| wording | entries |
|---|---|
| "Chapters 2–3 remain/are still undrafted" | 2026-09-17 … 2026-09-27 (12), **and 2026-10-17** |
| "Chapters 2 and 3 **are drafted**, but §4.8.2's form of error is not audited for elsewhere" | 2026-10-01 … 2026-10-06 (6) |
| neither | 2026-10-07 … 2026-10-16 (10) |

Every "undrafted" entry up to 09-27 was **true when written**. The list then
tracked the drafting correctly for six entries, carried nothing for ten, and
the entry after those ten revived the retired wording. The two failure modes
have different cures and are worth separating: *staleness* is cured by
re-checking an item, *revival* by making retirement checkable. An item absent
long enough stops being remembered as retired and becomes available to be
re-derived.

**The cause is outside the repository.** These sessions are driven by a
scheduled brief carrying a frozen description of repository state dated
2026-09-17, which contains the A1 sentence almost verbatim — eleven entries
after the list retired it. No audit inside the tree can read that brief, so
the cause is recorded in `graphene_standing_list_audit.py` as a
`NOT_MECHANISABLE` row with its reason rather than silently omitted, and the
compensating rule is procedural: **a claim about the state of the tree that
arrives from outside the tree is re-verified against the tree before it ranks
anything.**

**The instrument, and why its retired rows are the load-bearing half.**
`graphene_standing_list_audit.py` parses the latest entry's list, matches each
baseline assertion, and evaluates a predicate over the tree — so a verdict is
gated on the mechanism it reports rather than on text similarity, which is
Design-Verification-Roadmap's 2026-10-16 rule. The trap this file had to avoid
is §4.6.15's: a suite that erases its own finding the moment the finding is
acted on. Once A1 and A2 are retired from the list, "absent" must not read as
a pass, so a **retired** row asserts two things — that the wording stays gone
(returning it reports `REVIVED`, which is a fault class, and is the check that
would have caught 2026-10-17 on the day) and that the tree-side fact which
justified the retirement still holds (`RETIREMENT-STALE` otherwise, so
truncating Chapter 2 later would re-open the item rather than leave a silent
retirement). `ABSENT` is also a fault class, so a list the parser cannot read
cannot pass as clean.

Two by-products. **The module's first draft failed its own plant, on its own
subject matter.** `p_item16_ii_untouched` tested `phrase in file.read()` for a
sentence this file *wraps*, and for numbers written here with a Unicode minus;
it therefore reported §4.10 "untouched" about the file carrying the amendment
— it **confirmed the item it was written to contradict**, by answering "are
these exact bytes present" instead of "is this claim present". That is the
same defect as §4.6.15's two docstrings, which passed an audit by spelling
while being the thing it audited for, and as the 2026-10-17 grep for a
`setattr` route over computed names. Three for three: this repository's
instruments fail first by answering a weaker question than the one asked, and
the fix each time is to normalise or to put the question to an oracle rather
than to a string.

**And the instrument already existed, in the other repository.**
Design-Verification-Roadmap built `tools/open_item_staleness_audit.py` for
this failure on 2026-10-13, and on 2026-10-17 used it to retire an item for
the identical reason — *"it names work finished before the item started being
carried"*. That repository caught its instance on the same day this one
revived its instance. The difference between the two days was one module; the
generalisation is that an instrument built in one place does not protect the
other, which is §7.5.2 across repositories rather than across modules.


### 4.6.17 The two halves of item 16 do not compose: a replacement that reads a name the same flip rebinds

§4.6.9 built form `B_T` in the channel. §4.6.14 gave the contact a
replacement of its own instead of an exclusion. §4.6.12 measured what each
costs. The log of 2026-10-20 recorded, on that basis, that item 16's
"remaining content is the flip and nothing else".

**It was not.** Composition was left as well, and composition fails.

#### (a) The measurement

`graphene_default_flip.py` installs both halves at once and asks what
`n_contact` becomes. Against the committed contact-only form `B_T`:

| Metal | `n_contact`, form A (m⁻²) | committed `B_T` (m⁻²) | both halves, unrepaired (m⁻²) | `B_T` vs A | unrepaired vs `B_T` | ratio |
|---|---|---|---|---|---|---|
| Ti | 1.644935e16 | 1.541997e16 | 7.281739e15 | −6.2578 % | −52.7772 % | 8.43× |
| Cr | −0.0 | −0.0 | −0.0 | `n = 0` | `n = 0` | — |
| Cu | −1.380517e16 | −1.314672e16 | −6.391282e15 | −4.7696 % | −51.3850 % | 10.77× |
| Ni | −7.111301e16 | −6.227844e16 | −2.466546e16 | −12.4233 % | −60.3949 % | 4.86× |
| Au | −8.034794e16 | −7.031424e16 | −2.760960e16 | −12.4878 % | −60.7340 % | 4.86× |
| Pd | −8.343500e16 | −7.300861e16 | −2.859729e16 | −12.4964 % | −60.8302 % | 4.87× |
| Pt | −1.660634e17 | −1.463454e17 | −5.573671e16 | −11.8738 % | −61.9143 % | 5.21× |

The unrepaired composition overshoots the effect item 16 exists to install by
**4.86× to 10.77×, metal by metal, and in the same direction** — so it would
not have presented as a sign error, and a session that landed the flip and
regenerated transcripts would have committed it as the answer.

> The per-metal ratio is the honest statement, and the first draft of this
> table did not give it. It divided the smallest-magnitude unrepaired
> deviation by the **largest**-magnitude committed effect and the other way
> round, quoting "4.1× to 13.0×" from two pairs of numbers belonging to
> different metals. That is this chapter's recorded dominant failure mode — a
> claim answering a weaker question than the one it was written to ask — and
> it was caught by *reading the transcript*, not by any check in the file.

#### (b) The cause is one attribute read

The contact replacement is deliberately built on the **shipped**
`gfet.quantum_capacitance`, for the two reasons §4.6.12 records: the
alternative carries its own constants and differs by 9.48e6 ulp, and the
committed `log(2(1+cosh η))` spelling overflows where the fixed point's first
iterate lands. It reads that name by attribute lookup at call time:

```python
return gfet.quantum_capacitance(qca.A_FERMI * np.sqrt(nn), T=T)
```

and the channel half of the flip rebinds exactly that name. With both halves
installed, the contact's self-consistency is solved against a `C_q` that has
itself already been solved — once, against the wrong series capacitance.

#### (c) The general form, which no mechanism in this thesis covers

§4.6.13 established that the set a patch must reach is the name's runtime
**identity class**, so that no binding is left stale. §4.6.14 established that
the class must be **partitioned** by operating point, so that no member is
silently left at the shipped object. Both are statements about *which names
move*.

**Neither asks whether a replacement READS a name the same flip rebinds.**

A replacement that does is neither stale nor unassigned. It is correctly
installed, computed by the right mechanism, reported as installed by that
mechanism's own report — and it reads a different object than the one it was
validated against. This is a third failure mode beside the two §4.6.13 and
§4.6.14 were built for, and it is invisible to both. The audit that would
catch it without running anything is a one-line question asked of each
replacement: *which module-level names does this function read at call time,
and does the flip rebind any of them?*

> **Annotated 2026-10-22 (§4.6.18). The question as worded above answers
> *no* for this very replacement.** It says "*the* flip" — the flip that
> installs the replacement. `blast.contact_on_consistent_argument()`
> installs `_cq_consistent_contact` and rebinds `cd.quantum_capacitance`;
> the replacement reads `gfet.quantum_capacitance`. The intersection is
> **empty**. The defect is **cross-flip**: the name it reads is rebound by a
> *different* manager, which is only installed at the same time because
> item 16 composes them. The audit of §4.6.18 measures both forms and
> reports 1-of-4 under the composition against 0-of-4 under the wording
> above. The sentence is left unrevised because it is the record of what was
> asked; the correct question is *every flip installable concurrently*.

#### (d) Two mechanisms, and they partially cancel

The failure has two parts, not one. **(i)** the inner `C_q` is the channel's
form-`B_T` solve rather than the shipped form-A `C_q`; **(ii)** that inner
solve reaches `qca.series_factor`, which reads `gfet.C_ox` at call time, and
the contact solve has `C_ox` substituted by `C_interface` for its duration —
so the inner channel solve runs against the **contact's** series capacitance
too. Isolating (i) by restoring `C_ox` for the inner call:

| Metal | (i) alone | (i)+(ii) | (ii) contributes |
|---|---|---|---|
| Ti | −53.8545 % | −52.7772 % | +1.0772 |
| Cu | −52.4054 % | −51.3850 % | +1.0204 |
| Ni | −62.1287 % | −60.3949 % | +1.7338 |
| Au | −62.5403 % | −60.7340 % | +1.8063 |
| Pd | −62.6595 % | −60.8302 % | +1.8293 |
| Pt | −64.2107 % | −61.9143 % | +2.2965 |

**(ii) offsets (i) for every metal.** Repairing (ii) alone — which is the more
obviously wrong-looking of the two, since a global being mutated across a
nested call is visible on inspection — would have made the error *larger*.
The repair has to be made at (i).

> **Annotated 2026-10-22 (§4.6.18). (i) and (ii) are not two defects.** A
> direct count of how often the channel replacement is reached during the
> contact census, and of the value `gfet.C_ox` holds at each arrival, gives
> **192 arrivals under the unrepaired flip, 181 of them with `C_ox` already
> substituted by `C_interface`** — and **0 arrivals under the repaired
> flip**. (ii) cannot occur without (i), because (ii) is a property of the
> context (i)'s arrivals land in. One cause, two consequences, 181/192 of
> the arrivals carrying both at once. The percentages in the table are
> unrevised and remain the correct decomposition *of the error*; what is
> corrected is the reading that they are two independent faults. It also
> explains why a single change at (i) repairs both, which the table above
> could not have predicted.

#### (e) The repair, and what it does not do

Bind the contact replacement to the shipped object **captured at manager
entry** rather than to the live name. `on_default_flip()` then installs the
whole flip through the two managers §4.6.14 named — `bc.rebind()` for the
charge, whose class is `{gfet, crx}` and both members in the channel, and
`bc.rebind_by_operating_point()` for the capacitance, whose class spans the
two operating points. This is the first caller to give that manager a
**non-`None` contact assignment**; until §4.6.14 the contact had no
replacement of its own and the exclusion was correct exactly as long as that
was true.

Validated against exactly known values, not plausible ranges: the identity
assignment (every member → its own shipped object) is **bitwise** neutral on
all nineteen committed numbers; the flip's contact half reproduces the
committed contact-only form `B_T` **bitwise** on all seven metals; its channel
half reproduces `on_consistent_argument()` **bitwise** on all eight census
entries; `n_contact(ΔW_F = 0)` is exactly `0.0`; the electron/hole antisymmetry
is bitwise; an incomplete assignment raises and rebinds **nothing**; and a
positive control asserts the repair is load-bearing, so that if the two ever
agree this section's premise is known to have expired. 18 passed, 0 failed;
transcript byte-identical across two consecutive runs.

**With both halves in the circuit for the first time, every cost quoted in
Chapter 7 item 16 is reproduced exactly** — `n` −1.0212 %, `I_d` −0.2530 %,
`n_contact` −4.7696 % to −12.4964 %, §4.7 `R_extra` −1.2583 % to +47.6925 %,
Pd's implied residual +50.9489 → −203.2764 Ω·µm, additive decomposition
4-of-4 negative. That reproduction is the evidence that the defect lives in
the composition and in neither half.

**This section does not flip the default.** `graphene_fet_model.py` is
unedited and the chapter's charge relation is still (4.29). What has changed
is that the flip is a single installable, validated, costed object, so what
remains is a decision about what this thesis claims rather than a build.

> **Annotated 2026-10-22 (§4.6.18). The repair is correct and its
> correctness is contingent.** Binding to the object captured at entry
> protects the **name** the cell holds. It does not protect the names the
> **held object** reads at call time. It is sufficient here only because the
> shipped `gfet.quantum_capacitance` reads nothing in the rebind set — its
> own live reads are `{gfet._thermal_rounding}`. Its neighbour
> `gfet.carrier_density` reads `gfet.quantum_capacitance` at call time, so
> the identical construction applied to *it* is **not** immune, which
> §4.6.18 demonstrates with a planted control rather than arguing. Had the
> contact needed the shipped charge relation instead of the shipped
> capacitance, the same repair would have been written, would have looked
> right, and would have left the defect in place. Separately, the live read
> is still **in** the function on the `live=True` branch, and
> `make_contact_cq` is public: this is a default, not a removal.

### 4.6.18 Auditing what a replacement reads: three readers, a cross-flip relation, and the depth of a capture

§4.6.17(c) named the audit that would have caught its own defect without
running anything, and the log of 2026-10-21 ranked it first precisely because
that defect had been found by accident. This section builds it. The result is
that the exposure §4.6.17 found is, as framed there, **complete** — and that
the instrument built to confirm that found two things §4.6.17 did not ask
about, one of which bears on its own repair.

#### (a) Three readers, chosen so that agreement is evidence

A replacement's dependence on a name is a dependence on the **live binding**,
so the question is which `(module, attribute)` pairs it resolves at call time.

**Reader A — static.** A bytecode walk of the transitive closure of live name
lookups: every `LOAD_GLOBAL` of a repository module followed by a `LOAD_ATTR`,
every repository function reached through such a name or through a global,
followed into recursively. Closure cells are resolved from the **real function
object**, so a *captured* object is distinguished from a *live name read* —
which is the entire content of §4.6.17's repair. Branch-**insensitive**, hence
an over-approximation; §(d) turns that into a result rather than a caveat.

**Reader B — sentinels.** An identity counting wrapper is installed at each
member of the flip's rebind set, the replacement is called, and the counters
that fired are the names it read. Object substitution, so it sees a read
whether it is spelled `mod.attr` or held as a module global by
`from mod import attr` — the spelling behind all three of the October delivery
faults.

**Reader C — module-type interception.** Each module's `__class__` is swapped
for a `ModuleType` subclass whose `__getattribute__` and `__setattr__` record
the access and delegate. It sees two things B cannot: reads of
**non-callables** (`gfet.C_ox` is a float and cannot be wrapped) and
**writes**.

B and C are blind in **opposite** directions — a read made from inside a
name's own module is a globals-dict lookup and never reaches the module type —
so their agreeing is a test and not a restatement. They agree exactly on all
five replacements. The planted controls confirm the blind spot is real: three
of them are flagged by B and invisible to C. Rebind sets are obtained by
entering each manager under reader C and collecting its own `setattr` calls,
so a manager that grows a name is reported with it.

Both instruments are bitwise neutral on the seven committed contact densities,
checked before any count is used.

#### (b) The answer, and the question's wording

Of the four replacements §4.6.17 names, **exactly one** reads a name that is
rebound while it is installed: `blast._cq_consistent_contact`. The accident of
2026-10-21 found the whole of the exposure as it framed it.

Asked in §4.6.17(c)'s own words, however — *does the flip **that installs
it*** rebind any of them — the answer for that same replacement is **empty**.
The defect is cross-flip, and the table makes the two readings explicit:

| Replacement | own flip | the composition |
|---|---|---|
| `qca._carrier_density_consistent` | — | — |
| `qca._quantum_capacitance_consistent` | — | — |
| `blast._cq_consistent_contact` | **—** | **`gfet.quantum_capacitance`** |
| `flip.make_contact_cq(live=False)` | — | — |
| `flip.make_contact_cq(live=True)` | — | `gfet.quantum_capacitance` |

A guard written to §4.6.17(c)'s wording would have passed the composed flip.

#### (c) What the question did not ask: a write/read relation

Reader C reports that **all five** replacements read `gfet.C_ox` at call time
and that **three** of them write it. `gfet.C_ox` is in nobody's identity
class, so it is invisible to §4.6.13's mechanism, to §4.6.14's, and to the
question of §4.6.17(c), all three of which are about names a flip *rebinds*.

The relation that governs composition is therefore not a set of names but a
**bipartite relation between what replacements write and what replacements
read**, and §4.6.17(d)'s mechanism (ii) is one edge of it rather than a quirk
of one function. Its reachability under the repaired flip is currently zero —
but that zero is a property of the capture, not of the write, and it returns
the moment the capture is undone.

#### (d) Capture is one level deep

A cell protects the **name** it holds. It does not protect the names the held
**object** reads at call time:

| Captured object | its own live reads | in the rebind set? |
|---|---|---|
| `gfet.quantum_capacitance` | `{gfet._thermal_rounding}` | **no** — safe to capture |
| `gfet.carrier_density` | `{gfet._thermal_rounding, gfet.quantum_capacitance}` | **yes** — not safe |

> **Annotated 2026-10-24 (§4.6.20). Every entry in this table is correct and
> the conclusion drawn in its third column is not.** "In the rebind set?" asks
> whether a flip over the name exists *today*; "safe to capture" was read off
> it as if it answered whether the name can move the result. §4.6.20 measures
> the second question directly and finds that a rebind of
> `gfet._thermal_rounding` — the single live read of the captured object, the
> `no` row — moves the repaired contact replacement by **+58.3421 %**, to
> bitwise the same value a rebind of `gfet.quantum_capacitance` moves the
> *unrepaired* one to. Not being in a rebind set is not safety; it is the
> absence of a mechanism that would notice. The numbers and the read sets
> above are unrevised.

Both captures are the same construction and are indistinguishable at the call
site. §4.6.17's repair is safe because of a contingent fact about the one
function it happens to capture. A planted control holding `gfet.carrier_
density` by the identical construction is still flagged, which demonstrates
the limit rather than asserting it.

Stated in the form the other members of this family take:

- §4.6.13 — a repair to one binding is not a repair to the class.
- §4.6.14 — a class can span two operating points, so it must be partitioned.
- §4.6.17 — a replacement that *reads* a rebound name is neither stale nor
  unassigned.
- **§4.6.18 — a capture that fixes one name does not fix the names the
  captured object reads. The safe criterion is the transitive live-read set
  of the captured object.**

#### (e) What this section does not claim

It does not claim the repository is clean. It audits the four replacements
§4.6.17 names plus the two variants of the contact replacement; the eleven
surviving hand-listed rebinding sites are not replacements and are not
covered. It does not flip the default. And reader A is an
over-approximation by construction: it reports the `live=True` branch for
both variants of `make_contact_cq`, which is correct — the read is in the
function — and is why the three readers are reported separately rather than
reduced to one verdict.

Machine output: `replacement_read_audit_output.txt` (30/30 validations,
transcript byte-identical across two runs). Figure:
`replacement_read_audit.png`. Notes:
`notes/2026-10-22-what-a-replacement-reads.md`.


### 4.6.19 The write/read relation as a guard: the two composition defects are two edges of one relation

Section 4.6.17 reported a replacement reading a name its own flip rebinds.
Section 4.6.18 §4 reported `gfet.C_ox`, written by three of the five
replacements and read by all five, in no identity class and therefore
invisible to `bc.rebind()`, to `bc.rebind_by_operating_point()` and to
4.6.18's own question — all three of which are about names a flip *rebinds*.
It closed by stating that what governs composition is not a set of names but
a bipartite relation between what replacements write and what they read.

**Those are two sources of writes into one relation, not two findings.** A
name is written either by a replacement or by the flip itself, which is every
member of the rebind set. Composing the two makes the composed column —
the one 4.6.18 identified as the column that answers the question —
computable instead of argued. With reads traced and each edge derived from
the read:

| composition | observed `(writer, name, reader)` edges |
|---|---|
| unrepaired (`live_contact=True`) | `THE FLIP --gfet.quantum_capacitance--> make_contact_cq`; `make_contact_cq --gfet.C_ox--> qca._quantum_capacitance_consistent`; `make_contact_cq --gfet.C_ox-->` itself |
| repaired (`live_contact=False`) | `make_contact_cq --gfet.C_ox-->` itself |

The first unrepaired edge **is** 4.6.17's defect; the second **is**
4.6.18 §§4–5. The survivor is the contact substituting `gfet.C_ox` for
`cd.C_interface` and reading it back inside its own call — the same operating
point, restored before return — so it is the substitution and not a leak.
**One change removes both cross-replacement edges**, which is 4.6.18 §5's
"one cause, two consequences" as a measurement rather than as a reading of a
coincidence. No number in 4.6.18 is revised; what is revised is the status of
the relation, from a census to a guard.

**The guard refuses, and it refuses in both directions.** Each edge must
carry a verdict from `{GOVERNED, IMPOSSIBLE, OCCURS}` and a mandatory reason,
in the vocabulary `bc.rebind_by_operating_point()` already uses for
assignments; an unadjudicated edge raises before anything runs. What is new
is that completeness is static but a verdict is a claim about
**reachability**, which no static set contains — so every verdict is checked
against an order probe. `GOVERNED` is refused if the name is in no rebind
set, `IMPOSSIBLE` if the probe observes the edge, and `OCCURS` **if the probe
never observes it**. The last is the half a census cannot have: a guard that
only catches a wrongly-permissive claim answers a weaker question than it was
written to ask, and an `OCCURS` reason describing an edge that has stopped
occurring stays true-looking indefinitely because nothing asks it the
opposite question.

The operational result is that the repaired composition's adjudication — the
same stated reasons, unchanged — is **accepted** for `live_contact=False` and
**refused** for `live_contact=True`, naming three disagreeing verdicts. The
distinction between the two compositions is now a raise rather than a table.

Validation is anchored to exactly-known values: the 4.6.18 arrival counts are
**reproduced** rather than cited (0 arrivals repaired, 192 unrepaired, exactly
181 of those with `C_ox` already substituted), the repaired composition has
**exactly zero** cross-replacement edges, and the deriver itself is put under
synthetic traces whose answers are known without it — a deriver that returns
few edges and a repository that has few are otherwise indistinguishable.

Two defects in the work that produced this section are recorded because both
generalise. The first probe **did not trace reads**: it closed an edge on
every *arrival* following a write, which asks whether a replacement *ran*
after a name moved rather than whether it *read* it, and reported twelve
edges of which four were on names the reader does not read at all. Same shape
as the dominant failure mode, reached through the choice of which events to
instrument. The second is sharper: a **planted control found a crash in a
refusal path that nothing else reaches**, so without it the guard would have
shipped with one of its three refusals unable to report itself while looking
correct, because the two refusals the live subject exercises both work. The
generalisation, and the fourth member of the 4.6.14 / 4.6.17 / 4.6.18 family:
**a refusal that nothing reaches is not a refusal**, so every verdict in a
vocabulary needs its own planted control and not only the ones the subject
happens to exercise.

Machine output: `write_read_guard_output.txt` (21/21). Module:
`graphene_write_read_guard.py`. Notes:
`notes/2026-10-23-two-findings-one-relation.md`.


### 4.6.20 A capture is adjudicated by a counterfactual, and the repaired contact's exact twin was open all along

§4.6.17 repaired a composition defect by binding the contact replacement to
the shipped `gfet.quantum_capacitance` **captured at manager entry**, so that
a flip over that name could not reach it. §4.6.18 (d) recorded the limit of
that construction — a cell protects the name it holds, not the names the held
object reads — and judged the repair safe on the ground that the held
object's own live reads, `{gfet._thermal_rounding}`, are in no rebind set.
§4.6.19 made the write/read relation a guard. This section asks what the
capture is actually worth, and answers it with a measurement rather than with
a read set.

#### (a) The instrument used to make the claim cannot check it

§4.6.17's argument and §4.6.18's table are both written in the vocabulary of
`rra.static_reads`. On the two objects that the entire repair is the
difference between, that reader returns:

| | static live-read set |
|---|---|
| `flip.make_contact_cq(shipped_cq, live=False)` — repaired | 8 names |
| `flip.make_contact_cq(shipped_cq, live=True)` — unrepaired | the same 8 |

with symmetric difference **zero**, both containing `gfet.quantum_capacitance`
as a *live* read. `cq_of_n` carries both branches in its bytecode and a
bytecode walk cannot know which runs, so the reader is sound and
over-approximate exactly as designed — and therefore **cannot distinguish the
repaired construction from the unrepaired one at all**. The correctness
argument was never checkable where it was made.

#### (b) The verdict is a counterfactual, and it is exact

For each name on the capture surface — the holder's own live reads together
with the live reads of everything it holds — the name is rebound to a
deliberately wrong stand-in scaled by ×1.5, the holder is called, and the two
results are compared **bitwise**. Nothing is compared to a tolerance: *does a
rebind of this name reach this object* has an exact answer. An **identity**
rebind of all nine names moves the result by exactly zero, checked before any
verdict is used.

Contact replacement at `V = 0.5 V`, base `C_q = 6.555495551420368 × 10⁻²`
F/m² (`0x1.0c835a6063492p-4`):

| Name on the capture surface | verdict | repaired | unrepaired |
|---|---|---|---|
| `gfet._thermal_rounding` | EXPOSED | **+58.3421 %** | +58.3421 % |
| `qca.A_FERMI` | EXPOSED | +58.3417 % | +58.3417 % |
| `qca.series_factor` | EXPOSED | +26.1720 % | +26.1720 % |
| `qca._fixed_point` | EXPOSED | +22.4740 % | +22.4740 % |
| `cd.C_interface` | EXPOSED | +18.1154 % | +18.1154 % |
| `gfet.T` | EXPOSED | +0.0034 % | +0.0034 % |
| `cd.quantum_capacitance` | UNREACHED | exact 0 | exact 0 |
| `gfet.C_ox` | PROTECTED | **exact 0** | exact 0 |
| `gfet.quantum_capacitance` | PROTECTED | **exact 0** | **+58.3421 %** |

The two replacements differ on **exactly one** of the nine names, and it is
the one §4.6.17 is about. That column is the repair, measured.

#### (c) The result: the repair has an exact twin it does not cover

A rebind of `gfet._thermal_rounding` moves the **repaired** replacement to
`0x1.a92b527d6852cp-4` — **bitwise** the value a rebind of
`gfet.quantum_capacitance` moves the **unrepaired** one to.
`gfet._thermal_rounding` is the *sole* live read of the captured object, so
the route is the capture itself. §4.6.17 closed one path to that number and
left an exactly equally powerful one open **through the object it captured**.

It is in no flip's rebind set and is written by no replacement, so
`bc.rebind_by_operating_point`, §4.6.19's write/read guard and §4.6.17's
argument all reach it equally: not at all. **A name outside every rebind set
is not thereby harmless — it is only unwatched**, and §4.6.18 (d) is
annotated in place to say so.

#### (d) A second result, in the opposite direction

Walking the holder's closure cells on **object identity** recovers not one
held name but the held object's whole identity class,
`{gfet.quantum_capacitance, cd.quantum_capacitance}`, reproducing
`bc.identity_class(gfet, 'quantum_capacitance')` exactly and independently —
two derivations of one set, which is this section's validation against the
existing implementation in the overlapping limit. The capture therefore
neutralises a rebind of **both** members, including
`blast.contact_on_consistent_argument()`'s, which is **stronger** than
§4.6.17 claimed. But only one member is ever looked up, so calling both
PROTECTED would credit the capture with work it did not do — §4.6.18 (d)'s
conflation with the sign reversed. `UNREACHED` keeps the two apart, and
separating them needs the one judgement `static_reads` **can** make soundly:
an over-approximate read set has reliable *negatives*, so *not on the read
surface* is trustworthy even though *on it* is not.

> **Annotated 2026-10-25 (§4.6.21). The sentence this paragraph ends on is
> false.** `rra.static_reads` records a `(module, attr)` pair only when a
> module object is loaded and an attribute taken from it, or when a bare
> global resolves to a repo *function*. A bare global resolving to a
> module-level **scalar**, read from inside the module that defines it —
> which is how every constant here is read by its own module's functions —
> is recorded **nowhere**, so a negative means *not read in a spelling the
> walker can see*, not *not read*. The demonstration is a single name:
> `qca.A_FERMI` is **off** the read surface of both channel replacements
> while a rebind of it moves one of them by **+29.8900 %**, and it is **on**
> the surface of both contact holders, where `blast` spells it qualified.
> §4.6.21 measures this. **Nothing in (d) above is withdrawn**: the identity
> class, the two independent derivations of it and the `PROTECTED`/`UNREACHED`
> distinction all stand, and the `UNREACHED` verdicts of §4.6.20 are
> re-measured in §4.6.21 and still hold. What is withdrawn is the *ground*
> offered for them here. `UNREACHED` is sound when it rests on a measured
> bitwise zero, and unsound when it rests on an absence.

#### (e) A device-physics consequence, not only a software one

`gfet.C_ox` is **exactly** protected for both variants. §4.6.19 left one
`OCCURS` edge standing — the contact replacement's write-then-read of
`gfet.C_ox` inside its own call — adjudicated *the substitution itself and
not a leak* on the argument that the write and the read are at the same
operating point. The counterfactual converts that argument into a measured
zero: the **external** value of `gfet.C_ox` reaches the contact
replacement's output not at all, because the replacement substitutes
`cd.C_interface` for the whole solve and restores the old value before
returning. **The contact self-consistency is referenced to `C_interface` and
is exactly independent of the oxide capacitance**, so no §4.6 contact number
inherits an oxide-thickness error through that path. The three names that
*do* carry the contact `C_q` are `gfet._thermal_rounding` and `qca.A_FERMI`
at ≈ +58 %, and `cd.C_interface` itself at +18 %, which is the sensitivity
ranking a reader of §4.5 should have.

> **Annotated 2026-10-25 (§4.6.21). The `gfet.C_ox` result stands; the
> sensitivity ranking in the last sentence is incomplete and its "three
> names" is wrong.** The exact zero on `gfet.C_ox` is re-measured in
> §4.6.21 on `blast._cq_consistent_contact`, which has **no capture at
> all**, so the zero is a property of the contact model and not of §4.6.17's
> repair — a strengthening of this paragraph, not a retraction. The ranking
> is another matter. Five further names reach the contact `C_q` and were
> invisible to the surface this paragraph was computed from:
>
> | name | leverage on contact `C_q` |
> |---|---|
> | `gfet.e` | **+278.4911 %** |
> | `gfet._thermal_rounding` | +58.3421 % |
> | `qca.A_FERMI` | +58.3417 % |
> | `gfet.hbar` | −63.5273 % |
> | `gfet.v_F` | −63.5273 % |
> | `qca.E_CHARGE` | −21.2518 % |
> | `cd.C_interface` | +18.1154 % |
>
> The true ranking is headed by `gfet.e`, and `gfet.hbar` and `gfet.v_F`
> are interleaved among the three names named above. The numbers quoted in
> this paragraph are all correct; the claim that they are *the* ranking is
> not, and §4.6.21 (c) draws the device-physics consequence for `v_F`.

#### (f) It is a guard, not a census

`graphene_capture_surface_guard.py` raises `UnadjudicatedCaptureSurface`
before measuring anything if a surface name carries no `(verdict, reason)`,
and `FalseCaptureVerdict` after measuring if a verdict disagrees with the
counterfactual **in either direction**: PROTECTED on a name that moves,
EXPOSED on a name that does not, UNREACHED on a name that *is* read,
PROTECTED on a name that is read by nothing, GOVERNED on a name in no rebind
set. One refusal is §4.6.18 (d)'s error verbatim — a PROTECTED verdict on a
held name obliges every live read of the object held there to carry its own
verdict, and any that moves must be adjudicated EXPOSED. The repaired
capture's adjudication is **refused** for the unrepaired replacement, which
(a) showed the static reader cannot tell apart.

29/29 checks, with a planted control for every verdict and every refusal
branch — the §4.6.19 lesson that a refusal nothing reaches is not a refusal,
applied inside the module that states it. Seven defects injected one at a
time into a copy of the repository are all detected (7/7); undoing §4.6.17's
repair drops the score to 23/29, and replacing the bitwise comparison with a
1 × 10⁻³ relative tolerance is caught too, so the result depends on exactness
rather than on a threshold.

Stated in the form the other members of this family take:

- §4.6.13 — a repair to one binding is not a repair to the class.
- §4.6.14 — a class can span two operating points, so it must be partitioned.
- §4.6.17 — a replacement that *reads* a rebound name is neither stale nor
  unassigned.
- §4.6.18 — a capture that fixes one name does not fix the names the captured
  object reads.
- §4.6.19 — a refusal that nothing reaches is not a refusal.
- §4.6.20 — a name outside every rebind set is not thereby harmless, only
  unwatched. Whether a name can move a result is a counterfactual, and no
  read set answers it.
- **§4.6.21 — a name's absence from a read surface is a fact about the
  surface's vocabulary, not about the name. Two names for one parameter are
  not one knob, and a precomputed copy is unreachable from the name it was
  computed from.**


**Annotated 2026-10-26 (§4.6.22).** This section's `gfet.C_ox` **PROTECTED**
verdict on `blast._cq_consistent_contact`, and the number behind it, are
**correct and are not withdrawn** — but every leverage figure in this section
and in §4.6.21 was measured with the replacement's **flip not entered**, and
§4.6.22 shows that for this one holder the two configurations disagree about
eleven of fifteen names. The stated reason here generalised from one
configuration; §4.6.22 §(c) states the scope and corrects the wording in the
code, keeping the superseded text. The §4.6.20 counterfactual instrument
itself is vindicated: it is bitwise identical to the as-installed instrument
for five of the six objects the flips install.

### 4.6.21 The negatives were never sound, and `v_F` was never one name

§4.6.20 adjudicated one capture and closed with a sentence it did not test:
that an over-approximate read set has reliable *negatives*. This section
tests it on the three remaining replacements, finds it false, and follows the
consequence into a statement about the Fermi velocity that is a device-physics
result rather than a software one.

#### (a) The blind spot, and one name that exhibits it

`rra.static_reads` records a name in two situations only: a module object
loaded and an attribute taken from it (`qca.A_FERMI`), or a bare global
resolving to a repo **function** (`cq_thermal`). A bare global resolving to a
module-level **scalar** — the spelling every constant in this repository is
read by inside its own module — matches neither, and is recorded nowhere.

The cleanest exhibit is one name read two ways:

| holder | on the §4.6.20 surface? | leverage of `qca.A_FERMI` |
|---|---|---|
| `qca._carrier_density_consistent` | **no** | +0.0121 % |
| `qca._quantum_capacitance_consistent` | **no** | **+29.8900 %** |
| `blast._cq_consistent_contact` | yes | +58.3417 % |
| `flip.make_contact_cq(live=False)` | yes | +58.3417 % |

One name, one value, one mechanism; `qca` reads it bare and `blast` spells it
`qca.A_FERMI`, and that spelling is the whole of the difference. This is an
exact binary prediction of the mechanism on a single name, and it is checked
as such.

#### (b) The answer to §4.6.20's own question, in both directions

§4.6.20 asked whether `gfet._thermal_rounding` and `qca.A_FERMI` carry its
≈ +58 % onto the other replacements, which would make them repository-wide
levers rather than contact facts. The two halves answer oppositely.
`gfet._thermal_rounding` is **bitwise exactly zero** on both channel
replacements — a contact fact, and the zero is a *measured counterfactual*,
not an absence from a surface, which after (a) would prove nothing.
`qca.A_FERMI` is repository-wide, and off-surface in both places where it
matters.

Completing the surface on the four holders adds **18 (holder, name) pairs**,
8 distinct names, and **every one moves the value**; the largest is `gfet.e`
at +278.4911 % and the largest derived model prefactor is `qca.CQ_PREFACTOR`
at +50.4145 %. None had ever carried a verdict, because a guard can only
demand one for a name its surface names — and the §4.6.20 guard is shown to
**accept** an adjudication of all three replacements with every one of them
omitted. Repository-wide, **75 of the 90 module-level scalars (84 %)** are in
the blind class; `gfet.t_ox` and `gfet.epsilon_0` are scalars outside it,
because they are read only at import time to build `C_ox`, which pins the
census from the other side.

#### (c) The device-physics consequence: `v_F` is two names, not one

`C_q` carries the Fermi velocity twice, as `CQ_PREFACTOR = 2e³/(π ħ² v_F²)`
and as `A_FERMI = ħ v_F √π / e`, so a consistent change of `v_F` by *k* must
move one as *k*⁻² and the other as *k*⁺¹. The two families do not share a
name for it:

| holder | `gfet.v_F` | `qca.CQ_PREFACTOR` | `qca.A_FERMI` |
|---|---|---|---|
| `qca._carrier_density_consistent` | **exact 0** | +0.000671 | +0.000247 |
| `qca._quantum_capacitance_consistent` | **exact 0** | +1.008232 | +0.519524 |
| `blast._cq_consistent_contact` | −2.323463 | **exact 0** | +1.161680 |

*(d ln value / d ln name; `exact 0` = a ×1.5 rebind leaves the value bitwise
unchanged.)* The contact replacement reads `gfet`'s **live** prefactor;
both channel replacements read `qca.CQ_PREFACTOR` and `qca.A_FERMI`, which
`graphene_exact_charge_rederivation` computed from `gfet.v_F` **at import
time** and copied in. So `gfet.v_F` is bitwise unreachable from every
Chapter 4 channel number, and `qca.CQ_PREFACTOR` is bitwise unreachable from
the contact number. **A `v_F` sensitivity study performed by editing
`gfet.v_F` would have reported channel insensitivity to the Fermi velocity
as a result.** This is the §4.6.13 delivery fault in a new place: the name is
neither stale nor unassigned but **precomputed**, and `bc.identity_class`
cannot see it because the two values are not the same object — they are a
float and a float derived from it.

#### (d) And the perturbation is the physical case

`v_F` is not a constant of nature. Hwang *et al.* measure, by ARPES,
1.15 ± 0.02, 1.49 ± 0.08 and 2.49 ± 0.30 (×10⁶ m s⁻¹) for graphene on
SiC(000-1) (ε = 7.26), BN (ε = 4.22) and quartz (ε = 1.80) respectively,
against a fully screened *ab initio* value of 0.85 × 10⁶ m s⁻¹, with a 1/ε
dependence at fixed carrier concentration [Sci. Rep. **2**, 590 (2012),
doi:10.1038/srep00590].

The value committed throughout this chapter, `v_F = 1.00 × 10⁶ m s⁻¹`, sits
just above the **bare** value. The 90 nm SiO₂ back gate modelled in §4.2 has
ε ≈ 3.9, between the BN and quartz cases, so the substrate-appropriate value
is near 1.5 × 10⁶ m s⁻¹ — which is exactly the ×1.5 the probe applies by
default. The **−63.5273 %** it reports for the contact `C_q` is therefore not
a synthetic perturbation: it is the SiO₂-screened Fermi velocity, and this
chapter quotes the contact `C_q` at a `v_F` its own substrate does not have.

This is left as an open question rather than resolved here, and §7.9 records
it as such. Using the bare value consistently is defensible — it is what most
compact device models do, and the renormalisation is itself carrier-density
dependent, so any single number is already an approximation. Quoting it as
*the* contact `C_q` while the same chapter models a specific SiO₂ stack is
weaker. What is not defensible is the present state, in which the question
cannot be *asked* by changing one name.

#### (e) Anchored to an exactly known symmetry

All of (c) and (d) rest on a probe, so the probe is anchored to a value known
in advance: `gfet.hbar` and `gfet.v_F` must have **bitwise identical**
leverage, because the prefactor depends on them only through `(ħ v_F)²`. They
do — both move the contact `C_q` to `0x1.87bc5bc0c26d0p-6`. A probe reaching
the value by any other route would not reproduce that coincidence. 36/36
checks, with 9 defects injected one at a time into a copy of the repository
and all 9 detected; the mutation that hides this section's two findings behind
the measuring module's own apparatus carve-out scored **31/31** against the
first draft of the check meant to prevent exactly that, which is recorded in
`notes/2026-10-25-the-negatives-were-never-sound.md` §8 rather than quietly
fixed.

**Annotated 2026-10-26 (§4.6.22).** Every number in (a)–(d) above is measured
on the replacement **in isolation**, with its own flip not entered, and
§4.6.22 shows that this is bitwise equivalent to the as-installed measurement
for both channel holders and for `flip.make_contact_cq(live=False)` — so
(a)–(c) stand unrevised, and the ×1.5 `v_F` figure of **−63.5273 %** holds in
*every* configuration on the holder the default flip actually installs. The
one exception is `blast._cq_consistent_contact`, on which `gfet.v_F` is
**bitwise dead** once either other flip is entered. The 18 added pairs and
the 84 % blind-class census are unaffected: both are statements about
spelling, not about configuration.


### 4.6.22 A verdict measured in isolation is a verdict about a configuration the replacement never runs in

§4.6.20 and §4.6.21 adjudicate four *holders* — the replacement functions the
flip managers of §4.6.16 install. Both measure them by calling them
**directly**, with no manager entered. A replacement only ever runs with its
manager entered, and all three instruments used here are relative to the
namespace at the moment of use: `csg.capture_surface` and `vrs.value_reads`
are static walks that resolve each name to the object *bound now*, and
`csg.leverage` calls the holder in the current namespace. So the question
§4.6.21 left untouched is whether isolation and installation agree. Until this
section, nothing had asked.

#### (a) The population, derived rather than listed

`rra.FLIPS` × `rra.rebind_set` gives **8 (flip, site) pairs and 6 distinct
installed objects**. §4.6.20 and §4.6.21 between them audited four. The two
that had never had a read surface computed at all are the arrival-counting
wrappers `flip._counting` builds — the objects `gfet.carrier_density`,
`gfet.quantum_capacitance` and `crx.carrier_density` *actually become* under
the default flip.

#### (b) The answer: the instruments are vindicated on five of six objects

Measured over all 6 × 3 cells, with the dependence of each (holder, name)
pair classified on a four-valued axis (`BITWISE_SAME`, `BOTH_LIVE`,
`DEAD_AS_INSTALLED`, `LIVE_ONLY_AS_INSTALLED`):

| installed object | under its own flip | under the other two |
|---|---|---|
| `qca._carrier_density_consistent` | 14/14 `BITWISE_SAME` | invariant |
| `qca._quantum_capacitance_consistent` | 11/11 | invariant |
| `_counting(_carrier_density_consistent)` | 15/15 | invariant |
| `_counting(_quantum_capacitance_consistent)` | 12/12 | invariant |
| `flip.make_contact_cq(live=False)` | 19/19 | invariant |
| **`blast._cq_consistent_contact`** | **15/15** | **not invariant under either** |

**Every one of the six is bitwise flip-invariant under its own flip** — every
name on its surface, every counterfactual, and the holder's own value. That
proposition is what licenses every leverage figure quoted in §4.6.20,
§4.6.21 and §4.6.18, and it had never been stated, let alone checked. It is
now an assertion rather than an assumption.

A **surface** change is also not a **verdict** change, and the measurement
separates them: four of the six objects gain names on their surface in some
configuration — the counting wrapper for `carrier_density` goes 15 → 17 under
`qca.on_consistent_argument()` — while every one of their counterfactuals
stays bitwise identical.

#### (c) The one row that moves, and what it does and does not revise

`blast._cq_consistent_contact` under `flip.on_default_flip(live_contact=False)`
returns `0x1.4ef081d146707p-7` instead of `0x1.0c835a6063492p-4`, a factor of
**6.4134**, and **eleven of its fifteen names cross the exact-zero
boundary**: five go bitwise dead (`gfet.e` +278.4911 %, `gfet.hbar` and
`gfet.v_F` −63.5273 %, `gfet.kB` +0.0242 %, and `gfet._thermal_rounding`
+58.3421 % — the whole of §4.6.20's finding) and six come alive from bitwise
zero (`qca.cq_thermal`, `qca.cq_thermal_of_n` and `qca.CQ_PREFACTOR` at
+52.4103 %, `qca.KB` +35.8863 %, `qca.n_arg_EF_thermal` +9.0738 %).

**The value change is not a new result.** It is §4.6.17's — *the two halves do
not compose* — which measured this exact composition and costed it at −51 %
to −62 % in `n_contact`. What is new is that the same non-composition appears
as a change in the **read surface**, which is the object §4.6.19's and
§4.6.20's guards demand adjudication over. §7.9's standing item on the
unadjudicated compositions is therefore partly discharged: the relation is
measured over all 6 × 3 cells and exactly one row is non-composing.

What it revises is one **reason**, not one number. §4.6.21's dictionary gives
`gfet.C_ox` on this holder the verdict PROTECTED with the reason *"the contact
solve closes against `C_interface`, so the oxide capacitance reaches this
value not at all … the zero is a property of the contact model and not of the
2026-10-21 repair."* The first clause is right. The generalisation is not: the
same holder gives `gfet.C_ox` **+7.8110 %** under the default flip, so the
exact zero is a property of the contact model **composed with nothing else**.
The reason is corrected in place in the code and the superseded wording is
kept there.

#### (d) And the vocabulary cannot express the state of that name

Exhaustively over the four verdicts, for `gfet.C_ox` on this holder:

| verdict | in isolation | as installed |
|---|---|---|
| `PROTECTED` | admissible | refused — it moves |
| `UNREACHED` | refused — it is read | refused — it moves |
| `EXPOSED` | refused — bitwise zero | admissible |
| `GOVERNED` | refused — in no rebind set | refused |

Exactly one verdict is admissible in each configuration **and they are
different ones**. The four-verdict vocabulary of §4.6.20 is therefore not
merely incomplete here; non-composition is **not dischargeable by
adjudication at all**, which is the structural reason that item has survived
three sessions of being noted. A fifth value — a verdict indexed by
configuration — is the repair, and it is not attempted here.

#### (e) A second finding: `flip._counting` is a capture

Each counting wrapper's read surface contains exactly one name the
replacement it wraps does not — that replacement's **own name** — and a
rebind of it leaves the wrapper's value **bitwise unchanged**, because
`_counting` holds the function object in a closure cell. This is a **third**
capture in the model code, after the §4.6.17 contact repair and the
`bc.rebind` family, and nothing had noticed it; its consequence is that the
channel replacement *name* is immune under the default flip. The verdict is
`UNREACHED` and not `PROTECTED`, by §4.6.20's own argument for the other
identity-class member: the wrapper never *looks the name up*, so the bitwise
zero is not a capture neutralising a read. This section's first draft said
`PROTECTED` and the §4.6.20 guard refused it — the fourth time a rule planted
in an earlier section has caught a claim written in a later one.

#### (f) Exact anchors, and two of this section's own checks

The result rests on probes, so the probes are anchored to values known in
advance. `gfet.hbar` and `gfet.v_F` move `flip.make_contact_cq(live=False)` to
the **bitwise identical** `0x1.87bc5bc0c26d0p-6` in all four configurations,
because the `C_q` prefactor depends on them only through `(ħ v_F)²`; the
counting wrapper returns a value bitwise equal to the replacement it wraps;
an identity control flip leaves value and surface bitwise unchanged; and —
the sharpest of the four — **the repaired and unrepaired contact replacements
are bitwise indistinguishable in isolation**, both `0x1.0c835a6063492p-4`, so
the entire content of the §4.6.17 repair appears *only* as installed. Five
sections measured both holders in the one configuration in which they cannot
be told apart.

Nine defects were injected one at a time into a copy of the repository, each
required to be caught **by a named check** rather than by a non-zero exit, and
**two of this section's own checks failed that bar on the first run**: the
positive control for the new PROTECTED-as-installed refusal planted its false
verdict on a holder whose adjudication was incomplete in that configuration,
so the missing-name refusal fired first and the branch under test was never
reached (26/26, undetected); and the flip-invariance check asserted two hex
values but never the invariance *flag* computed from them, so a relative
tolerance could report a factor of 6.4 as invariant (27/28, undetected). Both
are repaired. Final: **28/28 checks, 9/9 defects detected on their own terms,
control clean.**

#### (g) What this section does not claim

It measures three configurations, not all of them: `live_contact=True` and the
*simultaneous* composition of two managers are not covered. "Invariant under
its own flip" is a bitwise statement at one probe point with one perturbation,
not a theorem. And `vrs.value_reads` remains over-approximate with unsound
negatives (§4.6.21 (a)); nothing here repairs that.

Machine output: `installed_surface_audit_output.txt`,
`installed_surface_mutation_output.txt`. Working notes:
`notes/2026-10-26-the-configuration-the-replacement-runs-in.md`.


### 4.6.23 A repaired instrument's own repair was mis-predicted twice, and the second prediction was the informative one

Every adjudication in §4.6.16 through §4.6.22 is taken with two constants:
`csg.PROBE_V`, the gate voltage the holder is evaluated at, and
`csg.PERTURBATION`, the factor the counterfactual multiplies by. Until
2026-10-27 both were captured as function defaults rather than read at call
time, which put them in the delivery class `MIXED` — some readers live, some
frozen — and `MIXED` is not an intermediate state. A rebind of a `MIXED` name
reaches the live readers and not the frozen ones, so if the live readers are
an audit's expected-value expressions and the frozen ones are its measuring
functions, the audit reports a disagreement its own mutation manufactured.
These two names are not incidental to this chapter: **every counterfactual in
it is taken with them**, so the split was a defect in the instrument and not
in any model the instrument measures. Three further names were `MIXED` for the same
reason in `graphene_chord_decay_length.py`, the module behind §7.9 item 2's
2026-10-20 annotation — the one that shows `γ` cancels from the chord/local
decay-length gap and that `ℓ` is the local decay length at no metal in the
table. Eleven of the repository's
twelve capture sites are converted; the twelfth is discussed below.

**No number in this chapter moves.** The two affected modules' machine
transcripts are bitwise unchanged by the conversion, which is the standard the
2026-09-30 and 2026-10-01 conversions were held to and the only available
evidence that a repair to a measuring instrument did not also move the
measurement. What changes is what a *future* mutation of those names will do.

**The informative part is that the repair was mis-predicted, twice.** The
suite written to measure it predicted the textbook `MIXED` signature on
`D1_CHORD`: the live reader `L'(d)` moves under a rebind and the captured
reader `⟨L'⟩` freezes. It does not freeze. `⟨L'⟩` captured the chord endpoint
`d₁` as a default, but the function it calls to evaluate `P` reads `D1_CHORD`
*live*, to re-centre the polynomial at `d − d₁`. So under a rebind the
pre-repair code averaged `L'` over the **old** interval while evaluating a
polynomial centred on the **new** `d₁`. The resulting number corresponds to no
consistent choice of `d₁` at all. That is strictly worse than the frozen case
the prediction assumed, because a frozen number is at least one the model once
produced, and the quantity carrying it was Eq. (D) — the γ-free discrepancy
this chapter quotes.

Restating the check as *coherence* rather than *movement* then caught a second
wrong draft of the same kind. The reference value for "the old interval" was
first taken **with the rebind still in place**, which made the reference the
same incoherent quantity it was meant to be distinguished from; it agreed with
the measurement to `0.000000e+00` and the check failed for the right reason.
A reference value has to be taken in the state it claims to describe. Under a
rebind of `d₁` from 3.3 Å to 3.9 Å the pre-repair default path lands
`1.617099e-02` from the new-interval value and `3.838368e-02` from the
pre-rebind one — between the two and equal to neither. After the repair it
equals the new-interval value **bitwise**.

**The residual, and a conflict this chapter does not resolve.** One capture
site is left: `n_quad=N_QUAD` in the velocity-saturation module. Its `def` sits
inside a function *body*, so the default is re-evaluated on every call of the
enclosing function and a rebind does arrive — it is live by scope, not frozen,
and is measured as live. It is also the only *real* site the
frozen-default-semantics suite has for demonstrating that a static pass can
find a capture its confirmation step then silently discards. So the standing
guard "no module-level constant is captured as a default anywhere" cannot go
green while that demonstration keeps its site, while its sibling "and in
particular none is `MIXED`" now can, for the first time since it was written
on 2026-10-01. The conflict is in the guard's **vocabulary** — it counts
syntactic captures and names the class it cares about only in the second check
— and not in the tree. Resolving it means either splitting the guard, or
trading a real demonstration site for a synthetic one, which is the trade that
suite exists to argue against.

Machine output: `mixed_capture_repair_output.txt` (33/33),
`mutation_arrival_probe_output.txt` (Section 7b, 32/2 → 33/1). Notes:
`notes/2026-10-27-the-repair-that-refuted-its-own-prediction.md`.


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

> **Annotation added 2026-10-13 (Section 4.6.12).** No number in the table
> above is changed or retracted, and the four `R_extra` values remain the
> committed ones. They are computed through
> `contact_edge_carrier_density()`, which evaluates `C_q` at the
> work-function potential rather than at `E_F/e` — the Section 7.9 item 16
> defect, in a second module, and one that no measurement of item 16 had
> reached because the harness rebinds `gfet.quantum_capacitance` while this
> module holds its own import-time reference. Corrected, `R_extra` moves
> −1.2583 % (Cu), +20.8352 % (Ni), +39.2281 % (Au), **+47.6925 % (Pd)**, and
> **Pd's +51.0 becomes −203.276**: all four residuals are then negative. The
> sentence below, "Only Pd gives a small, plausible positive residual", is
> therefore **true of the committed model and false of the corrected one**,
> and is left standing with this annotation rather than rewritten, pending
> the item 16 flip. The ordering of the four metals and the n/p type
> assignments are unchanged. The section's *conclusion* — that the additive
> decomposition does not hold — is strengthened, not weakened: a
> four-of-four failure is more consistent with the TLM double-count
> explanation given below than a three-of-four one was.

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
>
> *Amended 2026-10-13 (Section 4.6.12). The 9.6 % figure stands as computed
> and is not revised. What this annotation did not know is that it is a
> margin against INPUT error only, and that a known MODELLING error in the
> same quantity already exceeds it by 5x: correcting `C_q`'s argument in
> `contact_edge_carrier_density()` moves Pd's `R_extra` by +47.6925 %. The
> sentence "Pd is the fragile one" was right; the fragility was realised by
> the model and not by its inputs.*

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

**Annotated 2026-10-23.** Section 4.6.18's closing claim — that what governs
composition is a bipartite write/read relation rather than a set of names —
is **correct and was a census**. It enumerated the relation; it did not
refuse a composition containing an edge. Section 4.6.19 makes it a guard and,
in doing so, shows that 4.6.17's defect and 4.6.18 §4's `gfet.C_ox` are **two
edges of the same relation** removed by **one** change. The numbers in
4.6.18 stand unrevised; its status does not.
