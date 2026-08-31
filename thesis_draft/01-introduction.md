# Chapter 1: Introduction

*Draft section — single-layer graphene: electronic properties and device
applications. This is a working draft, updated as the accompanying
computational studies and literature review progress; see `notes/` for the
underlying research and `README.md` / individual scripts for the
computational results referenced below.*

## 1.1 Motivation

Since its isolation by mechanical exfoliation in 2004 (Novoselov et al.,
*Science* 306, 666, 2004), single-layer graphene has served two somewhat
separate roles in condensed-matter physics and electrical engineering. In
the first, graphene is a model system for two-dimensional Dirac physics:
its honeycomb lattice produces a linear, gapless band dispersion near the
K and K' points of the Brillouin zone, giving rise to massless Dirac
fermions, a half-integer quantum Hall effect, and a density of states that
vanishes linearly at the charge-neutrality (Dirac) point. These properties
are well established and form the material-level foundation of this
thesis (Chapters 2-3; see `graphene_band_structure.py`,
`graphene_transport_properties.py`, and the associated plots for the
computational treatment of band structure, density of states, and the
anomalous quantum Hall sequence).

In the second, more application-driven role, graphene is evaluated as a
candidate channel, contact, or interconnect material for post-silicon
electronics. This is the role this thesis is primarily oriented toward.
The same features that make graphene scientifically interesting — zero
bandgap, extremely high carrier mobility, atomic thinness — are a mixed
blessing from a device engineering standpoint. Zero bandgap rules out
graphene as a direct silicon-MOSFET replacement for digital logic (no
usable on/off ratio), but the same properties make graphene genuinely
promising for three narrower, well-defined application spaces that
semiconductor and device companies actively evaluate 2D materials for:

1. **RF/analog electronics**, where high carrier velocity and
   transconductance-to-capacitance ratio matter more than on/off ratio
   (Chapter 4; see `graphene_fet_model.py` and `rf_small_signal_model.py`).
2. **Interconnects**, where graphene's high current-carrying capacity and
   potential for lower resistivity at aggressively scaled linewidths are
   of interest as copper interconnect scaling runs into electron-surface-
   scattering limits (planned Chapter 5).
3. **Photodetectors and other optoelectronic devices**, which exploit
   graphene's broadband, gate-tunable optical absorption rather than a
   fixed absorption edge set by a bandgap (planned Chapter 6).

This thesis is organized to reflect that split: Chapters 2-3 establish the
intrinsic material physics computationally (band structure, transport,
optical absorption — already implemented), and Chapters 4 onward build up
device-level, application-specific analysis on top of that foundation,
aimed at the kind of device physics questions a semiconductor/device
company evaluating graphene would actually ask.

## 1.2 Why device physics, not just material physics

A recurring theme that emerges once the intrinsic material properties are
combined with realistic device geometry is that **most of the practically
important performance limits in a real graphene device come from effects
that are entirely absent from a pure band-structure treatment.** Three
examples, developed in detail in this thesis:

- **Metal-graphene contact resistance.** Because graphene is atomically
  thin, current must transfer into the channel through an interface
  (rather than a bulk 3D contact volume), producing a contact resistance
  that is fundamentally different in origin from a conventional
  semiconductor ohmic contact, and that dominates total device resistance
  once channel length is scaled down into the sub-micron regime relevant
  to competitive logic/RF nodes (see
  `notes/2026-08-21-contact-resistance-and-quantum-capacitance.md`,
  Section 2, for the literature review; representative values in the
  ~110-500 Ω·µm range are used throughout this thesis's device models).

- **Quantum capacitance.** Graphene's finite (linear, vanishing-at-Dirac-
  point) density of states means the total gate capacitance is not simply
  the oxide capacitance C_ox but the series combination of C_ox and a
  gate-voltage-dependent quantum capacitance C_q. This has a direct,
  quantitative effect on every downstream device metric computed in this
  thesis: channel carrier density, transconductance, and (Chapter 4) the
  RF cutoff frequency f_T, since f_T = g_m/(2π·C_gs) and C_gs is itself
  set by this same C_q(V_g) physics (see
  `notes/2026-08-21-contact-resistance-and-quantum-capacitance.md`,
  Section 3, and `graphene_fet_model.py`).

- **Absence of current saturation.** Unlike a conventional MOSFET, an
  ambipolar, gapless GFET does not show a flat, well-saturated Id-Vg
  characteristic; the transfer curve is V-shaped, with the location and
  depth of the minimum controlled by the interplay of gate overdrive,
  drain bias, and the local charge-neutrality point along the channel
  (visible directly in `gfet_transfer_characteristics.png`). This
  qualitative difference from silicon is why RF figures of merit for
  graphene show a characteristic pattern not seen in III-V or Si RF
  transistors: cutoff frequency f_T can be very high (hundreds of GHz in
  record devices), while maximum oscillation frequency f_max is set by
  output conductance and parasitic gate/access resistance, both made
  worse by the lack of saturation. Whether f_max lags f_T or exceeds it
  turns out to be a design outcome rather than an intrinsic rule — early
  device generations with large, single-finger gate resistance showed
  f_max well below f_T, while later devices with engineered low gate
  resistance (multi-finger, T-shaped gates) report f_max/f_T ratios above
  1 (Chapter 4, Section 4.6; see
  `notes/2026-08-22-rf-figures-of-merit-fT-fmax.md` and
  `notes/2026-08-26-fmax-parasitics-and-fT-fmax-ratio.md` for the
  literature benchmarks and `rf_small_signal_model.py` for the
  corresponding computational estimate, which models both the intrinsic
  and extrinsic/pad-capacitance-limited cases).

## 1.3 Thesis structure (current status)

| Chapter | Topic | Status |
|---|---|---|
| 2 | Electronic band structure of single-layer graphene (tight-binding model, Dirac dispersion, density of states) | Computational results complete (`graphene_band_structure.py`) |
| 3 | Transport and optical properties (quantum Hall effect, universal optical absorption) | Computational results complete (`graphene_transport_properties.py`, optical absorption scripts) |
| 4 | Device physics: contact resistance, quantum capacitance, GFET transfer characteristics, RF figures of merit | In progress -- DC model, RF small-signal model (incl. access resistance + extrinsic pad-capacitance estimate), contact-resistance-vs-channel-length crossover analysis, and spatially-resolved work-function-dependent contact doping model implemented (`graphene_fet_model.py`, `rf_small_signal_model.py`, `contact_resistance_crossover.py`, `graphene_contact_doping_model.py`); literature review of contact resistance, quantum capacitance, f_T/f_max (incl. the 2026-08-26 access-resistance/pad-capacitance revision), and contact doping complete; draft written (`thesis_draft/04-graphene-fet-device-physics.md`); per-metal Rc recalibration against 4 literature-sourced values (Cu, Ni, Au, Pd) attempted 2026-08-31 -- additive decomposition found NOT to hold for 3 of 4 metals (Section 4.7), an open modeling question rather than a completed recalibration; Ti/Cr not yet recalibrated (no literature Rc sourced) |
| 5 | Graphene interconnects: resistivity vs. linewidth, comparison to scaled copper | In progress -- literature review and resistivity-vs-linewidth model complete (`notes/2026-08-23-interconnect-resistivity-vs-linewidth.md`, `graphene_interconnect_model.py`); copper baseline includes a literature-calibrated liner/barrier-thickness effect and a parallel-conduction refinement, plus thin-Ru/Co-liner scenarios (`notes/2026-08-28-copper-liner-barrier-thickness-effect.md`, `notes/2026-08-29-liner-parallel-conduction-and-thin-liner-scenarios.md`) -- the realistic-edge-quality crossover conclusion now ranges from "no crossover" (bare 3nm TaN/Co, or thin-Ru-liner) to "~3 nm" (thin-Co-liner), i.e. materially liner-choice-dependent rather than a single number; draft written and updated (`thesis_draft/05-graphene-interconnects.md`); still to add: graphene-all-around-metal liner/cap model (a distinct architecture from the Cu-liner corrections above), sourcing the specific liner-material resistivity values against a quantitative reference |
| 6 | Graphene photodetectors: responsivity, gate-tunable absorption | In progress -- literature review complete (`notes/2026-08-24-photodetector-responsivity.md`); quantitative responsivity/gain-bandwidth model implemented (`graphene_photodetector_model.py`, `photodetector_responsivity_gain_tradeoff.png`); draft written (`thesis_draft/06-graphene-photodetectors.md`); still to add: plasmonic-enhancement factor, spatially resolved collection model (building block now available in Chapter 4, Section 4.5 -- `graphene_contact_doping_model.py` -- not yet integrated into the photodetector model itself) |
| 7 | Discussion and outlook: graphene's realistic near-term application space in the semiconductor industry | Not started |

This draft will be expanded chapter-by-chapter as the corresponding
research notes and computational studies are completed; see
`AUTOMATION_LOG.md` for the dated record of what has been added and when.
