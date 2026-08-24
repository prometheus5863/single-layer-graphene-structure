# Chapter 6: Graphene Photodetectors — Responsivity, Gain, and the Absorption Limit

*Draft section -- this chapter was "not started" as of the Chapter 1
status table; this is the first drafted content, written after the
corresponding research notes were completed (see
`notes/2026-08-24-graphene-photodetector-responsivity.md`). Unlike
Chapters 4-5, this chapter does not yet have a dedicated code/analysis
contribution -- see Section 6.4 for the specific follow-on planned.*

## 6.1 Motivation: the same absorption number, a different consequence

Chapters 2-3 of this thesis established graphene's universal optical
absorption, alpha_abs ~= pi*alpha ~= 2.3% per layer, independent of
photon energy, as a direct consequence of graphene's linear
Dirac-fermion dispersion (`graphene_transport_properties.py`,
`optical_absorption.png`). In Chapters 2-3 this number was presented as a
distinguishing materials-physics result -- a flat, universal absorption
spectrum is unusual and diagnostic of the underlying band structure. In
this chapter the same number reappears, but now as a design constraint:
it directly upper-bounds the responsivity achievable from a bare,
ungated, single-pass graphene photodetector to well under 20 mA/W, since
responsivity R = (eta*q*lambda)/(h*c) is linear in the achievable quantum
efficiency eta and eta is itself capped by the ~2.3% absorption fraction
absent any enhancement.

This is the third instance in this thesis of a single graphene material
property showing up in two different device contexts with opposite
implications: zero bandgap enables ambipolar FET operation (Chapter 4)
but also removes current saturation and complicates RF power gain;
finite mean free path and boundary sensitivity limits interconnect
scaling (Chapter 5) but the same edge/contact physics also sets contact
resistance in Chapter 4; and now, fixed low absorption is a clean physics
result in Chapters 2-3 but a hard performance ceiling in this chapter.

## 6.2 Two multiplicative bottlenecks

Real graphene photodetectors underperform even the na\"ive 2.3%-absorption
estimate, because two bottlenecks stack multiplicatively rather than
independently:

1. **Absorption bottleneck.** At most ~2.3% of incident photons interact
   with the single atomic layer at all (Section 6.1).
2. **Collection bottleneck.** Of the photocarriers that are generated,
   only those created within a diffusion length of a built-in-field
   region (a metal contact or p-n junction, typically ~100-200 nm) are
   separated and collected before recombining. Graphene's photocarrier
   lifetime is extremely short (~1 ps), so carriers generated away from a
   junction contribute little to photocurrent (Mueller, Xia & Avouris,
   *Nature Photonics* 2010).

Reported external quantum efficiencies for bare metal-graphene-metal
photodetectors are consequently only ~0.1-0.2%, well below even the
2.3% absorption ceiling. This is directly analogous to the "boundary
physics dominates" theme of Chapters 4-5: just as contact resistance
(Chapter 4) and edge scattering (Chapter 5) are set by interface quality
rather than bulk graphene transport, photocurrent collection here is set
by the geometry and quality of the built-in-field region rather than by
bulk graphene absorption.

## 6.3 Breaking the limit: gain, and the tradeoff it creates

Because graphene has no bandgap, it has no conventional avalanche
multiplication mechanism, so essentially all reported photoconductive
gain in high-performance graphene photodetectors comes from *extrinsic*
engineering rather than the graphene itself:

- **Photogating**, where trap states (defects, adsorbates, or a coupled
  narrower-gap partner material) capture one carrier species and leave
  the other to recirculate through the external circuit many times
  before recombining, can push responsivity up by many orders of
  magnitude -- reported values range from ~10^3 A/W (fast, ~400 ns
  interfacial photogating) to ~10^10 A/W (slow, trap lifetimes extended
  to the millisecond-second range).
- **Heterostructure integration** with silicon (Schottky-junction
  self-powered detectors, ~510 mA/W; plasmonically enhanced variants,
  ~1.9 A/W at 1535 nm) or with TMDs such as WS2/MoS2 (quantum-confined
  photogating, up to ~4.4x10^6 A/W at very low incident power) uses
  graphene primarily as a fast, highly conductive transport/contact
  layer while a stronger-absorbing or carrier-trapping partner material
  supplies the gain.
- Because gain and response speed trade off directly against each other
  in a photogating picture (longer trap lifetime -> more gain but slower
  response), recent (2025) work has instead used device *geometry* --
  e.g. alternating electron- and hole-conduction channels -- to partially
  decouple the two, reporting ~1.7x10^7 mA/W responsivity with a
  3-4 us response time, and heterostructure-engineered zero-bias designs
  achieving 160 Gb/s data rates by attacking the tradeoff directly at the
  device-architecture level rather than through trap engineering.

The qualitative shape of this tradeoff -- gain bought at the cost of
speed unless geometry intervenes -- mirrors gain-bandwidth tradeoffs in
conventional avalanche and photoconductive detectors, but is unusually
sharp in graphene because there is no intrinsic gain mechanism to fall
back on; every reported high-gain graphene device relies on an added
material or engineered interface, not on graphene's own band structure.

## 6.4 Planned follow-on work

- Implement a compact responsivity model in this repo's codebase,
  combining the existing wavelength-independent 2.3% absorption result
  (`graphene_transport_properties.py::calculate_optical_conductivity`)
  with a parametrized photoconductive gain factor G, to reproduce the
  qualitative tradeoff curve (responsivity vs. response time) shown
  qualitatively in Section 6.3, in the same style as the quantitative
  models already built for Chapters 4-5.
- Extend that model with an optional plasmonic-absorption-enhancement
  factor, analogous to how contact resistance was added as a
  literature-calibrated correction term to the Chapter 4 DC model.
- Connect the collection-bottleneck picture (Section 6.2) to the
  quantum-capacitance-limited channel electrostatics already modeled in
  Chapter 4, since both ultimately trace back to the same finite-DOS,
  short-carrier-lifetime physics of graphene near the Dirac point.

## References

See `notes/2026-08-24-graphene-photodetector-responsivity.md` for the
full literature review and citation list supporting this chapter.
