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

## 1.3 Thesis structure

The registered title of this thesis is *quantum transport in twisted bilayer
graphene–superconductor junctions*. Chapters 2–6 are its single-layer
foundation: the band structure, transport, contacts, interconnect and
photodetector physics of single-layer graphene that a twisted-bilayer junction
model inherits. Chapter 7 discusses what that foundation establishes, and
Chapter 8 outlines the step to twisted bilayer graphene–superconductor
junctions.

| Chapter | Topic | Role |
|---|---|---|
| 1 | Introduction | motivation and structure |
| 2 | Electronic structure of single-layer graphene: tight-binding model, Dirac dispersion, density of states, quantum capacitance | foundation |
| 3 | Transport and optical properties: mobility, minimum conductivity, quantum Hall effect, universal optical absorption, Pauli blocking | foundation |
| 4 | Device physics of the graphene FET: quantum capacitance, transfer characteristics, contact resistance and its metal dependence, RF figures of merit, contact geometry, electron–hole asymmetry | foundation |
| 5 | Graphene interconnects: resistivity versus linewidth, compared with scaled copper | foundation |
| 6 | Graphene photodetectors: responsivity and gain, plasmonic enhancement, metal-dependent photocarrier collection, the n/p crossover | foundation |
| 7 | Discussion and outlook: where single-layer graphene realistically fits, and what this thesis establishes | synthesis |
| 8 | Towards twisted bilayer graphene–superconductor junctions | outline |
| A | Verification record: audits of each model, audits of the audits, annotations, and the corrections-of-record table | record |
| B | A taxonomy of how this thesis's claims failed | record |

Each chapter keeps its model, the model's validation against exactly known
limits, its results and a short statement of its limitations. The record of how
each result was checked, corrected or superseded is in Appendix A under the
original section numbers; where that record changed a number the body reports,
the body says so in one sentence and links to it. The dated history of this
section's own status table is in [Appendix A, §1.3](appendix-A-verification-record.md#s1-3).
