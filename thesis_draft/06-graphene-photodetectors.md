# Chapter 6: Graphene Photodetectors — Responsivity, Gain, and the Absorption Limit

*Draft section -- first drafted 2026-08-24 after the corresponding
research notes were completed (see
`notes/2026-08-24-graphene-photodetector-responsivity.md`). Section 6.4
below adds the quantitative responsivity/gain model that was originally
only planned as follow-on work, closing out this chapter's code/analysis
contribution alongside Chapters 4-5 (`graphene_fet_model.py`,
`graphene_interconnect_model.py`).*

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

## 6.4 A compact responsivity/gain-bandwidth model

`graphene_photodetector_model.py` implements the tradeoff described
qualitatively in Section 6.3 as a quantitative model, following the same
compact-analytic-model philosophy as the Chapter 4 GFET model
(`graphene_fet_model.py`) and the Chapter 5 interconnect model
(`graphene_interconnect_model.py`): not a TCAD-grade solver, but a model
built from the same handful of literature-anchored parameters used
elsewhere in this thesis, intended to reproduce the right qualitative
and order-of-magnitude behavior.

**Bare-device responsivity.** Using the net external quantum efficiency
reported for bare (no-gain) MGM graphene photodetectors, EQE ~= 0.15%
(midpoint of the 0.1-0.2% range in Section 6.2 -- already the product of
the absorption bottleneck and the collection bottleneck), the standard
responsivity relation R = EQE*e*lambda/(hc) gives R_bare ~= 0.665 mA/W at
lambda = 550 nm -- comfortably under the ~20 mA/W ceiling implied by the
bare 2.3% absorption limit (Section 6.1), consistent with bare MGM
devices being reported in the sub-mA/W-to-mA/W range.

**Photoconductive gain.** Because graphene has no intrinsic (bandgap/
avalanche) gain mechanism, gain is modeled as extrinsic photogating in
its simplest form: a trap-captured carrier lifetime tau_trap versus a
carrier transit time tau_transit, giving the classic photoconductor gain
G = tau_trap/tau_transit and 3 dB bandwidth f_3dB = 1/(2*pi*tau_trap).
Reusing this thesis's own Chapter 4 device parameters (L = 200 nm
channel, mu = 0.4 m^2/(V.s), V_ds = 0.1 V -- the same GFET parameters
used throughout `graphene_fet_model.py`) rather than introducing new
free parameters gives tau_transit = L^2/(mu*V_ds) = 1.0 ps, which is
also consistent with graphene's independently-reported ~1 ps
photocarrier lifetime (Section 6.2), so the bare (tau_trap =
tau_transit) limit of this model reproduces the "fast, gain-free,
GHz-class" bare-device regime described in the literature review by
construction.

**Model vs. literature.** Sweeping tau_trap from 1 ps to 10 s and
comparing R(tau_trap) = R_bare*G(tau_trap) against the three literature
anchor points from `notes/2026-08-24-photodetector-responsivity.md`
(see `photodetector_responsivity_gain_tradeoff.png`, left panel):

| Literature point | tau_trap | R_literature | R_model | Ratio |
|---|---|---|---|---|
| Interfacial photogating | 400 ns | 1.0e3 A/W | 2.7e2 A/W | 3.8x |
| Alternating-channel (2025) | 3.5 us | 1.7e4 A/W | 2.3e3 A/W | 7.3x |
| Extended-trap-lifetime photogating | 1 s | 1.0e10 A/W | 6.7e8 A/W | 15.0x |

The model reproduces the correct *slope* of the tradeoff across six
decades of tau_trap (R scales linearly with tau_trap by construction,
and the three literature points -- measured independently, in unrelated
devices -- track that line to within a roughly constant offset), while
under-predicting absolute responsivity by a factor that itself grows
mildly with tau_trap (3.8x -> 7.3x -> 15.0x). The natural reading is the
gain-bandwidth invariant implied by this simple model,
GBP_model = 1/(2*pi*tau_transit) ~= 159 GHz (right panel of the same
figure), is a lower bound set purely by transit-time physics; real
photogated devices exceed it because trap capture cross-section, local
field concentration, and channel/contact geometry -- none of which are
represented by a single lumped tau_transit -- add gain beyond the bare
transit-time picture. Notably, the 2025 alternating-channel device,
explicitly engineered to beat the naive gain/speed tradeoff by geometry
rather than trap lifetime, shows the same kind of super-invariant GBP as
the trap-based photogating points rather than a qualitatively different
signature -- consistent with Section 6.3's reading that its improvement
comes from decoupling gain from speed, not from a fundamentally
different physical gain mechanism.

**Caveat.** As with the RF f_max model's known artifact (flagged in
`notes/2026-08-22-rf-figures-of-merit-fT-fmax.md`), this model's
absolute-magnitude offset is a known limitation, not a hidden one: a
single lumped tau_transit cannot capture device-specific trap physics,
and closing that gap (rather than just the tradeoff slope) is exactly
the kind of refinement Section 6.5 below scopes out.

## 6.5 Further follow-on work

- Extend `graphene_photodetector_model.py` with an optional plasmonic-
  absorption-enhancement factor, analogous to how contact resistance was
  added as a literature-calibrated correction term to the Chapter 4 DC
  model, to close part of the offset identified in Section 6.4.
- Connect the collection-bottleneck picture (Section 6.2) to the
  quantum-capacitance-limited channel electrostatics already modeled in
  Chapter 4, since both ultimately trace back to the same finite-DOS,
  short-carrier-lifetime physics of graphene near the Dirac point.
- Replace the lumped EQE_bare = 0.15% parameter with a spatially
  resolved diffusion-length collection model once the Chapter 4
  spatially-resolved contact-doping-profile follow-on (still open as of
  2026-08-24) is implemented, since both rely on the same underlying
  built-in-field-region geometry near a metal contact.

## References

See `notes/2026-08-24-photodetector-responsivity.md` for the
full literature review and citation list supporting this chapter.
