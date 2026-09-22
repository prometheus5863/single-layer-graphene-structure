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
| 4 | Device physics: contact resistance, quantum capacitance, GFET transfer characteristics, RF figures of merit | In progress -- DC model, RF small-signal model (incl. access resistance + extrinsic pad-capacitance estimate), contact-resistance-vs-channel-length crossover analysis, spatially-resolved work-function-dependent contact doping model, and a DFT-Fermi-shift-derived edge-vs-top contact geometry model implemented (`graphene_fet_model.py`, `rf_small_signal_model.py`, `contact_resistance_crossover.py`, `graphene_contact_doping_model.py`, `graphene_edge_contact_model.py`); literature review of contact resistance, quantum capacitance, f_T/f_max (incl. the 2026-08-26 access-resistance/pad-capacitance revision), contact doping, and contact geometry (2026-09-05) complete; draft written (`thesis_draft/04-graphene-fet-device-physics.md`); per-metal Rc recalibration against 4 literature-sourced values (Cu, Ni, Au, Pd) attempted 2026-08-31 -- additive decomposition found NOT to hold for 3 of 4 metals (Section 4.7), an open modeling question rather than a completed recalibration; Ti/Cr not yet recalibrated (no literature Rc sourced); edge-vs-top geometry model (Section 4.8) reproduces the qualitative direction and large-hole-diameter trend of Passi et al.'s patterned-contact data but not their small-diameter upturn or the full measured ~11x device-level reduction **Conditioned 2026-09-22 (Section 4.10):** the Section 4.7 recalibration table audited for near-cancellation. `kappa` and sign-robustness run in *opposite* order -- Pd's positive residual, the only row consistent with the additive decomposition surviving, flips sign on a 9.6% error in `R_extra`, while Cu's large negative needs 81%. The additive decomposition therefore fails **robustly**, and amplified input error is eliminated as the cause of the negatives (one branch of an item open since 2026-08-31, narrowed by measurement for the first time). |
| 5 | Graphene interconnects: resistivity vs. linewidth, comparison to scaled copper | In progress -- literature review and resistivity-vs-linewidth model complete (`notes/2026-08-23-interconnect-resistivity-vs-linewidth.md`, `graphene_interconnect_model.py`); copper baseline includes a literature-calibrated liner/barrier-thickness effect and a parallel-conduction refinement, plus thin-Ru/Co-liner scenarios (`notes/2026-08-28-copper-liner-barrier-thickness-effect.md`, `notes/2026-08-29-liner-parallel-conduction-and-thin-liner-scenarios.md`) -- the realistic-edge-quality crossover conclusion now ranges from "no crossover" (bare 3nm TaN/Co, or thin-Ru-liner) to "~3 nm" (thin-Co-liner), i.e. materially liner-choice-dependent rather than a single number; draft written and updated (`thesis_draft/05-graphene-interconnects.md`); still to add: graphene-all-around-metal liner/cap model (a distinct architecture from the Cu-liner corrections above), sourcing the specific liner-material resistivity values against a quantitative reference **Conditioned 2026-09-22 (Section 5.5):** the resistivity formula is reinforcing at every width (`kappa` <= 0.66), but the `lambda_impurity` calibration *behind* it is a three-term subtraction with `kappa = 19.71` -- **the worst-conditioned step in this thesis**, larger than anything in Chapter 4 (11.46) or Chapter 6 (4.55). It propagates into `rho(W)` damped but not removed (`kappa` 0.88 -> 1.55 over 18-52 nm). At the calibration width three sensitivities are **exactly zero** and `S(rho_bulk)` changes sign through it. The liner model's `W_eff = W - 2t` carries a previously unstated ~20% bar on `rho_eff` at 18 nm. |
| 6 | Graphene photodetectors: responsivity, gate-tunable absorption | In progress -- literature review complete (`notes/2026-08-24-photodetector-responsivity.md`, `notes/2026-09-06-plasmonic-enhancement-graphene-photodetectors.md`, `notes/2026-09-07-spatial-photocarrier-collection-model.md`); quantitative responsivity/gain-bandwidth model implemented (`graphene_photodetector_model.py`, `photodetector_responsivity_gain_tradeoff.png`); literature-calibrated plasmonic near-field absorption-enhancement model added (`graphene_plasmonic_photodetector_model.py`, `plasmonic_photodetector_enhancement.png`: 25x/8.5x on-resonance enhancement for two literature designs); spatially resolved, metal-dependent photocarrier collection model added (`graphene_photodetector_collection_model.py`, `photodetector_collection_efficiency_by_metal.png`: 1.00x-1.47x collection-enhancement factor across this thesis's 7-metal table, single-reinforcing-contact scope), integrating the Chapter 4 Section 4.5 contact-doping-profile machinery into the photodetector model as planned since 2026-08-24; self-consistent TWO-CONTACT collection model added 2026-09-17 (`graphene_photodetector_two_contact_model.py`, `photodetector_two_contact_net_response.png`, `notes/2026-09-17-two-contact-self-consistent-collection.md`), which **reverses** the single-contact model's metal ranking for a symmetric device (Pt goes from best to worst, overstated 5.5x) and reproduces the literature's symmetric-cancellation result to machine precision; signed, carrier-resolved two-contact model added 2026-09-18 (`graphene_photodetector_signed_carrier_model.py`, `photodetector_signed_carrier_response.png`, `notes/2026-09-18-signed-carrier-resolved-contact-fields.md`), which places the n/p crossover at its physical 5.4 eV rather than graphene's 4.5 eV [Giovannetti et al., PRL 101, 026803 (2008)] and tracks electrons and holes separately: it **confirms and doubles** the asymmetric-pair result (best pair inverts from Cr/Pt at |N| = 0.917 to Ti/Pt at |N| = 1.832, insensitive to the crossover) and **contradicts** the 2026-09-17 symmetric ranking, which turns out to be crossover-dependent and is no longer claimed -- all three of its validations are against exactly known values, two passing bitwise; draft written (`thesis_draft/06-graphene-photodetectors.md`, Sections 6.1-6.9); still to add: non-uniform illumination, bolometric-mechanism model, and a reconciliation of the 0.12 eV measured contact potential step with the 0.25-1.07 eV offsets the metal table assumes; **Section 6.9 added 2026-09-19** removes the uniform-illumination simplification standing since Section 6.6, showing illumination enters only as a weight g(x) on a fixed collection kernel k(x), and produces the chapter's first *bound* rather than a ranking: |N[g]| <= max|k| for every illumination pattern at fixed photon number. Its practical verdict is that asymmetric metallisation beats illumination engineering by ~4x per incident photon (~8x against the shadow mask Shimomura et al. actually built), so masks are worth using only on symmetric devices. **Section 6.11 added 2026-09-20** removes the linear dW -> doping assumption standing since Section 6.6, implementing Khomyakov et al.'s (PRB 79, 195425 (2009)) self-consistent square-root relation (`graphene_contact_doping_nonlinear_model.py`, `contact_doping_nonlinear_relation.png`, `notes/2026-09-20-nonlinear-work-function-to-doping-relation.md`). It **confirms** Section 6.8's Ti/Pt headline (1.832 -> 1.742, a 4.9% compression, stable over alpha = 0-5 eV^-1) but **retracts two claims of Section 6.9**: the 'max|k|/|N_uniform| < 1.02 for every asymmetric pair' bound (14 of 21 pairs violate it under 6.9's *own* linear model, worst Au/Pd at 22.4x) and the conclusion that masks are for symmetric devices only (a perfect mask gains 8.4x on Au/Pd per incident photon). The design recommendation survives on different grounds: 8.4x of a small number is still 4.9x worse than unmasked Ti/Pt. It also establishes that **four of this thesis's seven metals -- Ti, Ni, Pd (chemisorbed, d_eq < d0) and Cr (untabulated) -- lie outside the regime of the relation the chapter depends on**, Ti most consequentially, so a description of the chemisorbed metals that does not go through their work function is now the chapter's central open problem  **Section 6.12 added 2026-09-21** removes the last of the chapter's three hidden conventions -- the single 5.4 eV p/n crossover shared by all seven metals (`graphene_per_metal_crossover_model.py`, `per_metal_crossover.png`, `notes/2026-09-21-per-metal-crossover-from-the-chemical-interface-term.md`). Khomyakov et al. put the crossover at w_cross(d) = W_G + Delta_c(d), so 5.4 eV is that expression at the *physisorbed* separation, not a constant. Four exact validations (49/49 ordered pairs bitwise against the Section 6.8 model; 49/49 bitwise in the ell -> infinity limit; the anchor bitwise for 61/61 values of ell; charge conjugation exactly zero over 27 cases). Results: no metal changes doping type, but **Cu's dW moves 17.12%** -- more than three times Section 6.11's 4.9% compression, **falsifying this session's own pre-registered <=10% prediction**, which had been calibrated on that precedent; and the same-sign pair Cu/Au moves **77.9%** against 1.06-1.21% for every straddling pair, a **65-fold difference from identical inputs** that independently reconfirms Section 6.11's ratio-amplification mechanism under a perturbation of a different kind. The anchored exponential also **diverges** on the chemisorbed metals (w_cross of 6.25-62.6 eV, above the highest elemental work function anywhere), refusing Ti, Ni and Pd for a reason wholly unrelated to Section 6.11's d_eq < d0 -- two distinct parts of one framework failing on the same three metals. Net effect on the chapter: the Cu/Pt straddling design rule is robust to 1.2%, the Ti/Pt headline remains **untested** by two successive refinements rather than confirmed by them, and Sections 6.8-6.11 now carry an explicit +-17% (per-metal) / +-78% (same-sign pair) error bar rather than an implicit claim of exactness |
| 7 | Discussion and outlook: graphene's realistic near-term application space in the semiconductor industry | Not started -- synthesis task **restated 2026-09-18**. As framed on 2026-09-17 (Chapter 4's low-contact-resistance metals Ni/Au/Pd being the worst for photoresponse) the contradiction does not survive Section 6.8's signed model: at the physical n/p crossover Pd sits near the top of the symmetric ranking, and the symmetric ranking itself is no longer claimed. The durable and more interesting tension is narrower: Chapter 4 optimises a **single** junction, whereas the two-terminal photoresponse depends on the **difference between two** junctions (Section 6.8.4: |N| = 1.832 for Ti/Pt vs at best 0.556 for any symmetric pair) and is to first order blind to either junction's own quality. A second, now sharper thread for this chapter: the symmetric-device design rule (pick the metal that dopes graphene least) and the asymmetric one (maximise |dW_A - dW_B|) point in opposite directions, and Chapter 5's interconnect argument shares Chapter 4's single-junction framing. **Third thread, added 2026-09-19:** Section 6.9 turns the design question into a bounded one -- max|k| caps what any illumination pattern can achieve -- while Chapters 4 and 5 still argue from unbounded single-junction optimisation. Whether an analogous ceiling exists for contact resistance and for interconnect resistivity is the sharpest form of the synthesis question this chapter has yet been handed. **Fourth thread, added 2026-09-20:** that third thread is now *itself* the example, because Section 6.11 retracted the ceiling it was built on. Three of Chapter 6's last four sessions have overturned the previous one's headline, always by widening a sample rather than by finding a physics error -- which is the methodological point Chapter 7 should make in its own voice: this thesis's chapters differ not in rigour but in how many cases each claim was checked against, and Chapters 4 and 5 have never had that widening done to them  **Fifth thread, added 2026-09-21:** Section 6.12 supplies the first clean counter-example to that pessimistic reading and sharpens it into something Chapter 7 can actually argue. The session pre-registered four predictions *before* writing the model and committed them to git first; three held and one was falsified at 17.12%, and the falsification was informative rather than embarrassing precisely because the prediction had been stated with a number attached. The methodological claim Chapter 7 should make is therefore not 'this thesis's claims keep failing' but the sharper and more useful one: **a claim's reliability tracks how wide a sample it was checked against and whether its error bar was stated in advance, not how carefully its physics was derived** -- Chapter 6 now has six sections' worth of evidence for that, and Chapters 4 and 5 state general design rules from tabulated subsets with no error bars at all. The concrete synthesis task: Section 6.12 shows a near-cancelling quantity amplifies an input uncertainty 65-fold relative to a reinforcing one. Chapter 4's contact resistance and Chapter 5's liner scenarios both contain differences of comparable quantities, and neither has been asked whether it is a near-cancellation **Sixth thread, added 2026-09-22 -- and the fifth thread's concrete synthesis task is now DONE.** The fifth thread asked whether Chapter 4's contact resistance and Chapter 5's liner scenarios are near-cancellations. They were audited (Sections 4.10, 5.5) and the answer is yes, with two findings Chapter 7 should carry. (i) **The worst-conditioned step in the thesis is not in Chapter 6 but in Chapter 5's calibration**, `kappa = 19.71` against Chapter 6's 4.55 -- so the chapter that has failed most often publicly is not the one whose arithmetic is most fragile; it is the one that has been *looked at* most. That reframes threads four and five from a story about Chapter 6's unreliability into a story about **unequal scrutiny**, which is a claim about method rather than about graphene. (ii) `kappa` alone is the wrong instrument: in Section 4.7 it ran *opposite* to the robustness of the conclusions drawn from it, and the useful quantity was the **sign-flip margin** (4.23). Chapter 7's methodological claim should accordingly be stated as: a result's reliability tracks how wide a sample it was checked against, whether its error bar was stated in advance, **and which of its terms cancel** -- with the honest rider that two of this session's own six pre-registered predictions failed, both by generalising from a partial enumeration, the same failure mode as 2026-09-20 and 2026-09-21. |

This draft will be expanded chapter-by-chapter as the corresponding
research notes and computational studies are completed; see
`AUTOMATION_LOG.md` for the dated record of what has been added and when.
