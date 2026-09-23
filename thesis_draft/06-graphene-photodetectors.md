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

## 6.5 Plasmonic near-field absorption enhancement

Section 6.4 flagged its own model's under-prediction of absolute
responsivity as a limitation to be closed by "an optional plasmonic-
absorption-enhancement factor" (former Section 6.5, now folded into this
one -- see `notes/2026-09-06-plasmonic-enhancement-graphene-
photodetectors.md` for the full literature review). A resonant metal
nanostructure on or near graphene locally concentrates the incident
optical near-field; since a 2D sheet's absorbed power density scales
with local field *intensity* (|E_local|^2), this multiplies graphene's
effective absorption -- and, assuming unchanged downstream collection
efficiency, its EQE and responsivity -- in a narrow band around the
structure's plasmon resonance.

`graphene_plasmonic_photodetector_model.py` models this as a Lorentzian
intensity-enhancement factor F(lambda) applied multiplicatively to
`EQE_bare`, reusing `responsivity_bare()` from Section 6.4's model
unmodified with this wavelength-dependent EQE. Two designs are
calibrated against literature near-field/absorption enhancement numbers
(not device-to-device responsivity comparisons, which would bundle in
unrelated contact/bias/collection differences -- see the notes file's
discussion of why a third candidate design, Fang et al.'s bowtie
nano-antenna, is deliberately excluded from the fit for exactly this
reason):

- **Echtermeyer et al. (Nat. Commun. 2011)** Ti/Au finger-grating
  designs at a graphene p-n junction: on-resonance intensity enhancement
  F_max = 25x (from the paper's reported ~5x near-field *amplitude*
  enhancement, squared, since absorption is an intensity effect), at two
  measured resonances -- 514 nm (110 nm finger) and 633 nm (130 nm
  finger). Applied to this thesis's own EQE_bare = 0.15%, this raises
  on-resonance responsivity from 0.62-0.77 mA/W to 15.6-19.2 mA/W at the
  two resonances (`plasmonic_photodetector_enhancement.png`, left panel).
- **Bowtie-antenna array on a Si waveguide, telecom C-band
  (arXiv:1808.10823):** a directly-simulated single-element absorption
  enhancement of 8.5x at ~1550 nm, raising responsivity from 1.88 mA/W
  to 15.9 mA/W on resonance (same figure, right panel).

Both designs' resonance linewidth (quality factor Q) is an **assumed**
value (Q = 8 visible, Q = 6 telecom), not measured or fitted, since
neither source paper reports FWHM in the material reviewed this
session -- stated explicitly rather than absorbed silently into the
Lorentzian shape, consistent with this thesis's practice elsewhere (e.g.
the Chapter 4 Section 4.8.1 hole-array fill-fraction assumptions) of
flagging assumed vs. measured/fitted parameters.

`graphene_plasmonic_photodetector_model.py::summary_numbers()` also
shows an illustrative combination of the plasmonic EQE enhancement with
Section 6.4's photoconductive gain at tau_trap = 1 s -- explicitly
labeled as a demonstration that the two mechanisms multiply in this
compact picture, not a claim about any real device, since the telecom
bowtie/waveguide device's actual reported mechanism is photo-bolometric
(resistance change from absorbed-light heating under bias), not
photogating, and this thesis does not yet model the bolometric
mechanism at all (Section 6.6).

This closes part, not all, of Section 6.4's offset: the plasmonic factor
explains how a *spectrally narrow* enhancement can bridge much of the
gap between the bare-transit-time model and literature devices at their
specific operating wavelength, but the three literature anchor points in
Section 6.4's table were not reported as plasmonically-enhanced devices,
so this section's designs are a parallel demonstration of an available
lever, not a re-fit of that table.

## 6.6 Spatially resolved, metal-dependent photocarrier collection

Section 6.4's model applies a single literature midpoint, EQE_bare =
0.15%, regardless of contact metal. But the "collection bottleneck"
described in Section 6.2 -- only photocarriers generated within a
built-in-field region near a contact survive graphene's ~1 ps carrier
lifetime to be collected -- is exactly the same contact-induced-doping
physics already modeled quantitatively for Chapter 4's contact
resistance work (`graphene_contact_doping_model.py`, Section 4.5). This
section closes the gap flagged since this chapter was first drafted
(2026-08-24): integrating that spatial doping-profile machinery into the
photodetector side (full literature review and design derivation in
`notes/2026-09-07-spatial-photocarrier-collection-model.md`).

**Model.** `graphene_photodetector_collection_model.py` superposes the
contact-doping-induced field (differentiating the same saturating
profile used for `graphene_contact_doping_model.py`'s density profile)
on the bare-device model's uniform bias field E_bias = V_bias/L_channel
= 5.0e5 V/m, then integrates the resulting position-dependent drift
velocity to get each photocarrier's transit time to the contact, and
from that an exponential-lifetime survival (collection) probability
p(x) = exp(-t(x)/tau_carrier). Averaging p(x) over the channel and
normalizing against a bias-field-only baseline (physically, a
work-function-matched contact -- already present in this thesis's metal
table as Cr, W = 4.50 eV vs. graphene's 4.5 eV) gives a metal-dependent
**collection enhancement factor**, plotted alongside the resulting
metal-resolved modeled EQE in
`photodetector_collection_efficiency_by_metal.png`.

**Mechanism support from the literature.** Xia et al. (*Nature
Nanotechnology* 4, 839 (2009)) demonstrated zero-external-bias graphene
photodetection using two *different*-work-function metal contacts (Ti,
Pd), which is only possible if contact-metal work function sets a real,
usable built-in field in the channel -- direct experimental support for
this section's core mechanism. Weiss & Duan (*NPG Asia Materials* 5, e74
(2013)) state the mechanism explicitly ("contacting graphene with metals
of variable work function... creates a potential offset... that
facilitate[s] the separation and transport of photocarriers") and make
an important scope-limiting point this model respects: identical
contacts produce "equal and opposing" built-in potentials that partially
cancel in the net terminal photocurrent. Mueller, Xia & Avouris (*Nature
Photonics* 4, 297 (2010)) independently estimate the built-in-field
collection region at ~100-200 nm, consistent in order of magnitude with
(though not identical to) the lambda_decay = 250 nm decay length reused
here from `graphene_contact_doping_model.py` (itself sourced from
Khomyakov et al. 2010's DFT-derived doping-decay length) -- a useful
independent cross-check, not a perfect match, and reported as such.

**Results.** Sweeping the seven contact metals already in this thesis's
`METAL_WORK_FUNCTIONS` table, the modeled collection-efficiency
enhancement ranges from 1.00x (Cr, the baseline) to 1.47x (Pt, the
largest work-function mismatch, |dV| = 1.15 V), i.e. a modest,
physically plausible boost -- not a multiple-order-of-magnitude effect,
since the doping field only dominates within roughly one lambda_decay of
the contact and the 200 nm channel modeled here is comparable in length
to that decay length. Commonly used real contact metals fall in between:
Ti and Cu (both common adhesion/contact metals) give ~1.21-1.23x; Ni,
Au, and Pd (higher work function, used specifically for their
historically low contact resistance -- Section 4.5/4.7) give ~1.39-1.41x.
This offers an additional physical rationale, alongside the
injection-resistance/quantum-capacitance reasoning of Chapter 4, for why
the high-work-function metals already favored for low contact resistance
are also a reasonable choice for photodetector contacts specifically.

**Scope limitation (stated explicitly).** This model represents only the
single contact whose doping field reinforces the bias-driven carrier
sweep direction; the real device's second contact, where the doping
field would partially oppose that sweep (per Weiss & Duan's "equal and
opposing" result, generalized to the biased case), is not modeled. A
fully self-consistent two-contact treatment is left as an open item
below rather than approximated away silently.

> **Resolved, and it overturns this section's ranking (2026-09-17).**
> Section 6.7 now models both contacts self-consistently. The result is
> not a small correction to the numbers below: for a *symmetric*
> two-terminal device the metal ranking is essentially reversed, and Pt
> -- this section's best metal at 1.47x -- becomes the worst, overstated
> here by a factor of 5.5. **The enhancement factors reported in this
> section remain valid only for a single junction in isolation, which is
> not the device geometry this chapter is otherwise about.** They are
> retained rather than deleted because the single-junction calculation
> is the correct building block and the limiting case that validates the
> two-contact model, but any device-level conclusion should be taken
> from Section 6.7.

**What this result can and cannot claim.** The relative, metal-to-metal
ordering (larger work-function mismatch -> shorter transit time -> higher
survival probability -> higher effective EQE) follows directly from the
mechanism above and is this section's main result. The *absolute*
EQE_bare = 0.15% literature midpoint this factor scales is not
attributed to any specific contact metal in the material reviewed for
Section 6.2, so the resulting absolute EQE_model(metal) numbers should
be read as relative-to-that-midpoint predictions, not validated absolute
values for each named metal.

A further correction, added 2026-09-17: the relative metal-to-metal
ordering described in the paragraph above -- the part called "this
section's main result" -- is precisely the part Section 6.7 overturns for
a two-terminal device. The chain "larger work-function mismatch ->
shorter transit time -> higher survival probability" is sound for one
junction, but in a symmetric device the same larger mismatch also
produces a stronger *opposing* field at the other contact, and the
second effect wins.

## 6.7 Self-consistent two-contact collection, and the reversal of Section 6.6's metal ranking

*(Heading retained for the record. The "reversal" it names was
downgraded to crossover-dependent by Section 6.8 on 2026-09-18; this
section's Result 2, by contrast, was confirmed and strengthened there.)*

Section 6.6 closed one gap and opened another, which it named: it models
a single contact and assumes every photocarrier is collected there. The
geometry makes that assumption untenable rather than merely approximate.
The modeled channel is L = 200 nm while the contact-doping decay length
reused from Chapter 4 is lambda_decay = 250 nm, so neither contact's
field has decayed appreciably by mid-channel -- the second contact is not
a perturbation on the first, it is comparable to it everywhere.

**Model.** `graphene_photodetector_two_contact_model.py` places contact
A (metal A) at x = 0 and contact B (metal B) at x = L, and defines a
signed total field, positive meaning "sweeps carriers toward contact A":

    F(x) = E_bias + g_A(x) - g_B(x)

where g_A is byte-for-byte the doping-field profile of Section 6.6 and
g_B is its mirror image about mid-channel. The minus sign on g_B is the
entire physical content that Section 6.6 was missing. Because F can
change sign inside the channel, the model assumes no destination: a
carrier drifts along F and is collected at A (signed +1) only if F > 0
along its whole path to x = 0, or at B (signed -1) only if F < 0 along
its whole path to x = L. A carrier whose path crosses a sign change
drifts into a **stagnation point** (F = 0) rather than a contact and is
counted as uncollected -- a real feature of the two-contact geometry,
since with two opposing doping fields and a modest bias there is
generically an interior null. Survival to the contact uses the same
tau = 1 ps photocarrier lifetime as Section 6.6. The figure of merit is
the signed, channel-averaged net response N(A,B), directly comparable to
Section 6.6's collection efficiency.

**Two validations.** Both are computed and printed, not asserted. First,
switching off contact B's doping field (by giving it graphene's own work
function) must collapse the model onto Section 6.6's result, since then
F > 0 everywhere and every carrier reaches A: it does, for all seven
metals, to a worst relative deviation of 7.8e-08. Second, at zero bias
two *identical* contacts must give exactly zero net response, because
F(x) = g(x) - g(L-x) is antisymmetric about mid-channel -- this is Weiss
& Duan's "equal positive and negative flow with a net zero
photocurrent", and Shimomura et al.'s observation that in symmetric devices
"the polarities of the photovoltages at each graphene/electrode
interface ... are canceled out under macroscopic light irradiation". The
model reproduces it to |N| <= 6.6e-17, i.e. machine precision, as an
emergent consequence rather than an input. The antisymmetry of the
metal-pair matrix under contact swap is a third consistency check.

**Result 1: the reversal.** At the working bias V_bias = 0.1 V
(E_bias = 0.5 MV/m), comparing Section 6.6's single-contact collection
efficiency against the symmetric two-contact device's actual net
response:

| metal (both contacts) | W (eV) | single-contact (6.6) | symmetric two-contact (6.7) | ratio |
|---|---|---|---|---|
| Ti | 4.33 | 0.776 | 0.672 | 0.87 |
| Cr | 4.50 | 0.632 | 0.632 | 1.00 |
| Cu | 4.65 | 0.766 | 0.673 | 0.88 |
| Ni | 5.04 | 0.877 | 0.321 | 0.37 |
| Au | 5.10 | 0.885 | 0.295 | 0.33 |
| Pd | 5.12 | 0.888 | 0.287 | 0.32 |
| Pt | 5.65 | 0.929 | 0.169 | 0.18 |

> **Downgraded by Section 6.8 (2026-09-18).** The reversal below is
> **crossover-dependent and is no longer claimed as a result.** Section
> 6.8's signed, carrier-resolved model reproduces this section exactly
> under its own restrictions (holes only, |ΔW|, crossover at 4.5 eV --
> agreement to 0.000e+00 over 98 comparisons), so the arithmetic here is
> right. But placing the n/p crossover at its physical value of 5.4 eV
> instead of graphene's 4.5 eV puts Pt **back at the top** (0.556, best)
> rather than the bottom. Three sections have now given three different
> answers to this one question, and the honest reading is that the
> symmetric metal ranking is not an established result of this thesis.
> The numbers below are retained because the calculation is correct for
> the convention it states; the *ranking* is not. What survives is
> Section 6.8.5's structural statement: in a symmetric device contact
> doping is purely parasitic, so the best metal is the one that dopes
> graphene least.

Section 6.6 ranked the metals Pt > Pd > Au > Ni > Ti > Cu > Cr. For a
symmetric device the ranking is essentially reversed: Cu ~ Ti > Cr > Ni >
Au > Pd > Pt. The mechanism is clear once the sign is right: a strongly
doping contact sweeps carriers toward *itself*, so two of them facing
each other produce a large opposing-field region and a mid-channel
stagnation point that the 0.5 MV/m bias cannot overcome against Pt's
~4.5 MV/m near-contact field. **Strong contact doping helps an isolated
junction and hurts a symmetric two-terminal device.** Cr is the one
metal unaffected (ratio exactly 1.00), because its work function matches
graphene's to within 0.01 eV and so contributes essentially no field to
cancel, leaving the bias to sweep the channel unopposed.

> **Restated by Section 6.8.6 (2026-09-18).** The tension as framed
> below -- that the low-contact-resistance metals are the worst for
> photoresponse -- does not survive the signed model: at the physical
> crossover Pd sits near the top of the symmetric ranking, not the
> bottom. The durable tension is narrower and different in kind:
> Chapter 4 optimises a *single* junction, whereas the two-terminal
> photoresponse depends on the *difference between two* junctions and is
> to first order blind to either one's own quality. Chapter 7 should
> synthesise that version.

This puts two of this thesis's own device-design arguments in tension,
which is worth stating plainly rather than smoothing over. Chapter 4
favours high-work-function metals (Ni, Au, Pd) for low contact
resistance, and Section 6.6 appeared to agree for photocarrier
collection. Section 6.7 shows the opposite for the photodetector: the
metals that are best for contact resistance are the worst for symmetric
two-terminal photoresponse. That is a real trade-off in device design,
not a modelling artefact, and a device-oriented treatment should make
the choice explicit rather than inherit it from the contact-resistance
chapter.

**Result 2: asymmetric contacts, at zero bias.** Because the symmetric
diagonal is exactly zero at zero bias, every nonzero entry isolates the
work-function-asymmetry mechanism:

| A \ B | Ti | Cr | Cu | Ni | Au | Pd | Pt |
|---|---|---|---|---|---|---|---|
| Ti | 0.000 | 0.596 | 0.067 | -0.753 | -0.794 | -0.802 | -0.903 |
| Cr | -0.596 | 0.000 | -0.563 | -0.835 | -0.849 | -0.853 | -0.917 |
| Cu | -0.067 | 0.563 | 0.000 | -0.778 | -0.805 | -0.812 | -0.905 |
| Ni | 0.753 | 0.835 | 0.778 | 0.000 | -0.078 | -0.103 | -0.581 |
| Au | 0.794 | 0.849 | 0.805 | 0.078 | 0.000 | -0.025 | -0.504 |
| Pd | 0.802 | 0.853 | 0.812 | 0.103 | 0.025 | -0.000 | -0.480 |
| Pt | 0.903 | 0.917 | 0.905 | 0.581 | 0.504 | 0.480 | 0.000 |

This is the quantitative version of the Xia et al. (2009) Ti/Pd
zero-bias device that Section 6.6 already cites as experimental support
for the mechanism: a work-function-asymmetric contact pair produces a
net zero-bias photoresponse where a symmetric pair produces none.

One non-obvious feature, worth flagging because it distinguishes a
mechanism model from a monotonic fit in |W_A - W_B|: the largest
response is Cr/Pt (0.917), not Ti/Pt (0.903), despite Ti/Pt having the
larger work-function difference (1.32 eV vs 1.15 eV). Ti carries its own
doping field (0.17 eV of mismatch) that sweeps carriers back toward Ti
and partially opposes Pt's, whereas Cr contributes no opposing field at
all. Under this model, the best zero-bias pairing is one strongly doping
contact against a *work-function-matched* one.

**The simplification that could overturn Result 2.** Both doping fields
are taken as sweeping carriers toward their own contact, through
|W_metal - W_graphene| -- the same magnitude convention Section 6.6 uses,
and the convention that makes identical contacts cancel exactly. But it
erases the n/p distinction: Ti and Cu (W < W_graphene) dope graphene
n-type while Pt, Pd and Au (W > W_graphene) dope it p-type, and a signed
treatment tracking electrons and holes separately would let an n/p pair
such as Ti/Pt *add* rather than partially cancel for one carrier
species. That would plausibly make Ti/Pt the strongest pairing and
invert the Cr/Pt-vs-Ti/Pt conclusion above. **Confirmed in Section 6.8.4
(2026-09-18): Ti/Pt is the strongest pair, |N| = 1.832, roughly double
the 0.917 quoted here, and the result is insensitive to where the
crossover is placed. The "Result 1 does not depend on it" clause that
follows, however, turned out to be wrong -- see the callout at the head
of Result 1 above.** Result 2's ordering should
therefore be read as provisional; Result 1's reversal does not depend on
it, since it concerns identical contacts, where the two conventions
agree. Full discussion in
`notes/2026-09-17-two-contact-self-consistent-collection.md`, Section 4.

Figure: `photodetector_two_contact_net_response.png` -- (a) the two
opposing fields and the Pt/Pt stagnation point at x = 100 nm, (b) the
single-contact vs symmetric two-contact comparison showing the reversal,
(c) zero-bias net response against W_A - W_B.

## 6.8 Signed, carrier-resolved contacts: the n/p crossover, and what survives of Section 6.7

Section 6.7 named its own principal simplification and predicted what
relaxing it would do. This section relaxes it. The prediction about
Result 2 is confirmed and strengthened; the prediction that Result 1
would be unaffected is **wrong**, and Section 6.7's Result 1 has to be
downgraded as a consequence.

### 6.8.1 Two things were wrong with `|W_metal − W_graphene|`, not one

The magnitude convention hid a sign, as Section 6.7 said. It also hid a
**crossover**, which Section 6.7 did not anticipate.

The sign of the contact doping is not `sign(W_metal − W_graphene)`.
Graphene's work function is 4.5 eV, but the n/p crossover for a metal
*on* graphene sits near **5.4 eV**, because the short-range
metal–graphene chemical interaction adds a ~0.9 eV interface dipole on
top of vacuum-level alignment [Giovannetti *et al.*, *Phys. Rev. Lett.*
**101**, 026803 (2008)]. Chapter 4's `METAL_WORK_FUNCTIONS` table
already carried this knowledge in its annotations — it calls Cu an
"n-type dopant" despite 4.65 > 4.5, and Pt "clearly above the ~5.4 eV
crossover" — but `|W − 4.5|` cannot express it.

Under the physical crossover, **six of the seven metals in this thesis's
table n-dope graphene; only Pt p-dopes it**:

| metal | W (eV) | ΔW vs 4.5 eV | type | ΔW vs 5.4 eV | type |
|-------|--------|--------------|------|--------------|------|
| Ti | 4.33 | −0.17 | n | −1.07 | n |
| Cr | 4.50 | +0.00 | n | −0.90 | n |
| Cu | 4.65 | +0.15 | **p** | −0.75 | n |
| Ni | 5.04 | +0.54 | **p** | −0.36 | n |
| Au | 5.10 | +0.60 | **p** | −0.30 | n |
| Pd | 5.12 | +0.62 | **p** | −0.28 | n |
| Pt | 5.65 | +1.15 | p | +0.25 | p |

Because the crossover changes which metals are which *type*, it is
carried as an explicit parameter of the model rather than chosen once,
and every result below is quoted under both values.

### 6.8.2 The model

With ΔW = W_metal − w_cross signed, the Dirac point in the channel is
rigidly shifted near each contact using the same saturating profile and
the same λ = 250 nm as Sections 6.6–6.7:

```
E_D(x) = ΔW_A/(1 + x/λ) + ΔW_B/(1 + (L − x)/λ)
```

A rigid band shift is an electron potential energy, so the field is
E(x) = (1/e) dE_D/dx, and **both carriers drift in that single field, in
opposite directions**: v_h = +μE, v_e = −μE. Each species is transported
by the same stagnation-aware logic as Section 6.7. The figure of merit
is the net charge delivered to contact A per absorbed photon,

```
N = (1/L) ∫ [ (+1)·s_h(x)·p_h(x) + (−1)·s_e(x)·p_e(x) ] dx
```

with s = +1 for collection at A, −1 at B, 0 for a carrier stranded at an
interior field null. **N now ranges over [−2, +2] rather than [−1, +1]**,
because both carriers can be collected — a possibility the single-species
magnitude convention could not represent at all.

### 6.8.3 Validation against three exactly known values

Section 6.7 recorded the lesson that an exactly-known validation beats a
plausible-range one, after an exact zero exposed an `inf − inf` bug that
a range check had passed. All three checks here are exact.

| # | Check | Worst deviation |
|---|-------|-----------------|
| 1 | Holes only, \|ΔW\|, w_cross = 4.5 eV must reproduce Section 6.7 across all 49 ordered pairs at both biases (98 comparisons) | **0.000×10⁰** (bitwise) |
| 2 | Identical contacts at zero bias ⇒ N = 0 exactly | 6.6×10⁻¹⁷ |
| 3 | Charge conjugation: N(−ΔW) = −N(ΔW) at zero bias | **0.000×10⁰** |

Check 1 is a genuine *reduction*, not an approximate agreement: under
those three restrictions the two models are algebraically the same
expression, E(x) = −F(x) of Section 6.7, so "a hole moves toward A" is
exactly "F > 0". Check 3 becomes statable only once the model is signed,
and would catch a carrier-mixing error that check 2 passes.

### 6.8.4 Result 2 confirmed, and roughly doubled

Section 6.7 predicted that a signed treatment "would plausibly make
Ti/Pt the strongest pairing and invert the Cr/Pt-vs-Ti/Pt conclusion".
It does:

| | best zero-bias pair | \|N\| |
|---|---|---|
| Section 6.7, magnitude convention | Cr/Pt | 0.917 |
| this section, w_cross = 5.4 eV | **Ti/Pt** | **1.832** |
| this section, w_cross = 4.5 eV | **Ti/Pt** | **1.830** |

The ordering inverts and the magnitude roughly doubles. The result is
**insensitive to the crossover** (1.832 vs 1.830) because Ti and Pt lie
on opposite sides of *both* candidate crossovers, which makes it the
most secure quantitative statement in this chapter.

The mechanism is visible in a stagnation audit at zero bias:

| pair | N | field nulls | holes stranded | electrons stranded | h→A | h→B | e→A | e→B |
|------|------|---|------|------|------|------|------|------|
| Pt/Pt (p/p) | +0.000 | 2 | 0.00 | **1.00** | 0.357 | 0.357 | 0 | 0 |
| Pd/Pd (n/n) | +0.000 | 2 | **1.00** | 0.00 | 0 | 0 | 0.368 | 0.368 |
| Ti/Pd (n/n) | −1.599 | **0** | 0.00 | 0.00 | 0 | 0.720 | 0.879 | 0 |
| Ti/Pt (n/p) | −1.832 | **0** | 0.00 | 0.00 | 0 | 0.907 | 0.925 | 0 |
| Cr/Pt (n/p) | −1.810 | **0** | 0.00 | 0.00 | 0 | 0.896 | 0.914 | 0 |

For an unequal pair the interior null **disappears entirely** and both
species are collected, at opposite ends. In a same-type symmetric pair,
one whole carrier species is stranded at the null while the other splits
evenly and cancels — so N = 0 there for two independent reasons at once,
which is a stronger statement of Weiss & Duan's symmetric-cancellation
result than Section 6.7 could make.

Ti/Pd reaches 1.599 with *both* metals n-type under the physical
crossover. What matters is therefore |ΔW_A − ΔW_B|, not the sign pair —
the practically important form of the statement, since Pt is the only
p-type metal available. This is also the design that Mueller, Xia &
Avouris (*Nature Photonics* **4**, 297 (2010)) actually built:
interdigitated fingers of Pd/Au against Ti/Au, chosen because identical
electrodes give a symmetric built-in field and zero total photocurrent,
yielding 6.1 mA W⁻¹ at 1.55 μm — a 15-fold improvement. The second-best
pair in the table is the pair the experiment used.

### 6.8.5 Result 1 contradicted: the symmetric ranking is not established

Section 6.7 asserted that its Result 1 "does not depend on" the signed
treatment, since identical contacts are where the two conventions agree.
That reasoning was incomplete: the conventions agree on the *sign
structure* for identical contacts, but the signed model also collects
the second carrier species, and under the physical crossover it places
the metals differently. The symmetric-device response at V_bias = 0.1 V:

| metal | Section 6.7 (magnitude) | signed, w_cross = 5.4 eV | signed, w_cross = 4.5 eV |
|-------|-------------------------|---------------------------|---------------------------|
| Ti | 0.672 | 0.180 | 0.987 |
| Cr | 0.632 | 0.210 | **1.264** |
| Cu | 0.673 | 0.245 | 1.101 |
| Ni | 0.321 | 0.436 | 0.321 |
| Au | 0.295 | 0.495 | 0.295 |
| Pd | 0.287 | 0.518 | 0.287 |
| Pt | **0.169 (worst)** | **0.556 (best)** | **0.169 (worst)** |

The same question has now received three answers in this thesis:
Section 6.6 (single contact) made Pt **best**; Section 6.7 (two contact,
magnitude) made Pt **worst**; this section makes Pt **best again** under
the physical crossover and **worst** under the vacuum one.

The correct conclusion is not that Pt wins. It is that **the symmetric
two-terminal metal ranking is not an established result of this thesis**:
it flips with a modelling choice — where the n/p crossover sits — that
Sections 6.6 and 6.7 never had to make explicit, because a magnitude
convention conceals it. Section 6.7's Result 1 should be read as
crossover-dependent, not as a settled reversal of Section 6.6. Its
numbers are retained in place, annotated, because the calculation is
correct for the convention it states.

What *does* survive, and is worth more than the ranking: under the
physical crossover the symmetric response is **monotone in |ΔW|**, from
Ti (|ΔW| = 1.07 eV, worst) to Pt (|ΔW| = 0.25 eV, best). The reason is
structural. In a symmetric device the contact doping field is *identical
at both ends*, so it contributes nothing to the net current and does
nothing but create a null that strands carriers. **For a symmetric
two-terminal device, contact doping is purely parasitic, and the best
metal is whichever one perturbs graphene least.** That statement is
crossover-independent in form — it is also exactly why Cr, with ΔW = 0,
wins the 4.5 eV column — and it is the form in which this result should
be carried into Chapter 7.

### 6.8.6 Design rules, and their effect on the Chapter 4 tension

Two rules, pointing opposite ways:

* **Symmetric device:** choose the metal that dopes graphene least
  (closest to the crossover).
* **Asymmetric device:** maximise |ΔW_A − ΔW_B|.

The asymmetric device wins outright — |N| ≈ 1.83 against at best 0.56 —
so for a photodetector the first rule is largely a statement about what
to avoid.

This **softens, and partly dissolves, the Chapter 4 vs Section 6.7
tension** that Section 6.7 flagged and that Chapter 1 lists as Chapter
7's synthesis task. Section 6.7 claimed the low-contact-resistance
metals (Ni, Au, Pd) are the *worst* for photoresponse. Under the signed
model that claim is crossover-dependent and, at the physical crossover,
reversed — Pd sits near the top, not the bottom. The durable tension is
narrower and different in kind: Chapter 4 optimises a *single* junction,
while the photodetector's figure of merit is a *difference between two*
junctions, so a contact metal chosen for its own low resistance is being
chosen on a criterion that the two-terminal photoresponse is, to first
order, blind to. Chapter 7 should synthesise *that*, not the metal-list
contradiction as Section 6.7 stated it.

### 6.8.7 What this section does not claim

The ΔW → doping-profile relation is taken as linear with a single λ for
every metal. Giovannetti *et al.* find it only roughly linear, and not
linear at all for the chemisorbed metals (Ti, Ni, Pd). Only the sign
structure and the crossover were changed here; the profile magnitude
convention is inherited unchanged from Chapter 4 Section 4.5.
Illumination is uniform, transport is drift-only with μ_e = μ_h, and
there is no photogain, photo-thermoelectric or bolometric contribution.
A genuine unreconciled discrepancy: the measured potential step at a
graphene/electrode interface is ≈ 0.12 eV with doping extending
0.2–0.3 μm into the channel [Mueller *et al.*, arXiv:0902.1479], well
below the 0.25–1.07 eV offsets this table assumes, because the
measurement is of the residual step in a gated device rather than the
flat-band charge transfer. Reconciling the two is left open.

Figure: `photodetector_signed_carrier_response.png` — (a) the field
profiles, showing that an n/p pair has no interior null where a
same-type pair has two; (b) zero-bias response against a Pt
counter-electrode under both conventions; (c) N against ΔW_A·ΔW_B,
separating reinforcing from cancelling pairs.

> **Annotation added 2026-09-21 (Section 6.12).** The single
> `W_CROSS_CHEM = 5.4 eV` used throughout this section and Sections
> 6.9–6.11 is now known to be the crossover *at the physisorbed separation*
> rather than a constant shared by all metals. Nothing in this section is
> retracted and no number below is changed, but every ΔW here carries an
> uncertainty of up to **17%** (and every same-sign pair up to **78%**) from
> that assumption alone — quantified, for the three metals it can be
> quantified for, in Section 6.12.3 and 6.12.5.

## 6.9 Non-uniform illumination: the collection kernel, and a ceiling on illumination engineering

Every two-contact model in Sections 6.6–6.8 declares uniform illumination
as an explicit simplification, and each one names the same experiment as
the reason it matters. This section removes the simplification, and the
removal turns out to cost almost nothing, because illumination never
enters the transport problem at all.

### 6.9.1 Illumination is a weight on a kernel the model already computes

Section 6.8's response is an integral over the channel:

    N = (1/L) ∫₀ᴸ [ q_h s_h(x) p_h(x) + q_e s_e(x) p_e(x) ] dx
      ≡ (1/L) ∫₀ᴸ k(x) dx

The bracket — the **collection kernel** k(x) — is the net charge delivered
to contact A per photon *absorbed at x*. It depends only on the two contact
metals, the bias, the crossover and τ. A generation weight g(x) therefore
requires no new transport physics whatever:

    N[g] = (1/L) ∫₀ᴸ g(x) k(x) dx ,    with (1/L) ∫₀ᴸ g(x) dx = 1

The normalisation fixes the total number of absorbed photons, which is the
only basis on which two illumination patterns can be honestly compared.

Three consequences follow at once, and all three are exactly checkable.

1. **g ≡ 1 must return Section 6.8's number bitwise.** Not to within a
   tolerance: under g ≡ 1 the expression is literally the same one.
2. **k(x) *is* the delta-spot photocurrent scan**, since N[δ(x−x₀)] = k(x₀).
   The kernel is not an intermediate quantity to be integrated away; it is
   the predicted position scan of a scanning-photocurrent measurement.
3. **max|k| is a hard ceiling.** For any normalised g, |N[g]| ≤ maxₓ|k(x)|,
   with equality only when the light is concentrated where |k| peaks. No
   mask, spot, grating or plasmonic pattern can beat it at fixed photon
   number. This is the first quantity in this thesis that bounds an entire
   design space rather than ranking points inside it.

For a **symmetric** pair at zero bias k is antisymmetric about mid-channel,
which gives the leaky shadow mask in closed form. With transmission T over
the left half, 1 over the right, and A ≡ (1/L)∫_{L/2}^{L} k dx,

    N(T) = 2A (1 − T)/(1 + T)   ⇒   **N(T)/N(0) = (1 − T)/(1 + T)**

### 6.9.2 Validation against five exactly known values

Continuing the practice established in Section 6.7 — an exact check beats a
plausible range — `graphene_photodetector_nonuniform_illumination_model.py`
is tested against five quantities known in advance:

| # | check | result |
|---|---|---|
| 1 | g ≡ 1 reproduces §6.8 over 196 (pair, crossover, bias) cases | **196/196 bitwise**, deviation 0.000e+00 |
| 2 | mirror-symmetric g on a symmetric pair at zero bias → 0 | \|N\| ≤ 1.3e−16 |
| 3 | mask left = −(mask right) | \|dev\| ≤ 2.2e−16 |
| 4 | leaky-mask closed form (1−T)/(1+T) | \|dev\| ≤ 2.2e−16 |
| 5 | \|N[g]\| ≤ max\|k\| over 343 (pair, pattern) cases | 0 violations |

Check 4 is the one that tests the physics rather than the plumbing: had the
equal-photon normalisation been dropped, the ratio would have come out as
(1 − T) instead, a 50% error at the experimentally relevant T = 0.5, and
none of checks 1–3 would have caught it.

### 6.9.3 The experiment, and a citation this thesis had wrong

The shadow-mask device is

> K. Shimomura, K. Imai, K. Nakagawa, A. Kawai, K. Hashimoto, T. Ideguchi
> and H. Maki, "Graphene photodetectors with asymmetric device structures
> on silicon chips", *Carbon Trends* **5**, 100100 (2021).

Sections 6.6–6.8 and the reference list previously attributed this to
"Suzuki *et al.*, *Carbon Trends* 5, 100115 (2021)" — wrong on both the
first author and the article number, introduced 2026-09-17 and copied
forward twice. Corrected here against the publisher record and the
Ideguchi group's publication list. Nothing in the physics changes; the
episode is recorded because a citation that is never re-fetched propagates.

Their device deposits 50 nm of nickel over one of the two
graphene/electrode interfaces, and their statement of the cancellation
modelled since Section 6.6 is explicit: "since the polarities of the
photovoltages at each graphene/electrode interface … are opposite, they
are canceled out under macroscopic light irradiation."

The number this thesis had not used is that **the mask leaks**: they report
the photovoltage under the mask as "about half of the opposite side", so
T ≈ 0.5. By the closed form that leaves **exactly one third** of a perfect
mask, not one half. The derivative of (1−T)/(1+T) at T = 0 is −2, so the
first few percent of leakage cost twice their face value.

### 6.9.4 Result: a mask rescues a symmetric device — and it is still the wrong strategy

Under uniform light a symmetric device gives identically zero at zero bias.
A mask over one interface is the only mechanism by which it can respond at
all — which matters more since Section 6.8 retracted the symmetric *metal
ranking*. Per absorbed photon, symmetric Ti/Ti with a perfect mask reaches
N = 0.841, against 1.832 for the best asymmetric pair (Ti/Pt) under uniform
light: 45.9%.

**That comparison is wrong, and the correct one reverses its emphasis.**
Responsivity is amps per incident watt. A mask that shades half the device
puts half the incident light into 50 nm of nickel, where the detector never
sees it. Per absorbed photon a mask only *redistributes*; per incident
photon it *discards*. The two accountings differ by exactly (1+T)/2:

    per absorbed photon : N(T)/N(0) = (1 − T)/(1 + T)
    per incident photon : N(T)/N(0) = (1 − T)

so a perfect mask gives exactly half as much per incident photon as the
per-absorbed-photon figure suggests. Corrected:

| device, per incident photon | vs. Ti/Pt under uniform light |
|---|---|
| symmetric Ti/Ti, perfect mask (T = 0) | **22.9%** |
| symmetric Ti/Ti, the mask actually built (T = 0.5) | **11.5%** |

**Asymmetric metallisation beats illumination engineering by roughly 4×,
and by about 8× against a mask anyone has fabricated.** Shimomura *et al.*'s
second design — a comb-shaped counter-electrode that *enlarges* one
interface rather than shading the other — avoids the incident-photon
penalty entirely and is, on this accounting, the more promising of their
two ideas.

### 6.9.5 Result: what illumination engineering can buy on an asymmetric pair, bounded

The obvious expectation is that masking cannot help a device that is
already asymmetric, because its kernel does not change sign. **That
expectation is false**, and it was written down before the model was run so
that it could be falsified: masking Ti/Pd's A side gains +0.7% and Ti/Pt's
B side +0.03% per absorbed photon. Sign is the wrong criterion — a
single-signed kernel is still not flat, and concentrating photons where
|k| is larger gains something.

The right statement is the bound, and it is exactly computable:

    (best achievable) / (uniform)  =  maxₓ|k(x)| / |N_uniform|

which is 1.0025 for Ti/Pt, 1.0159 for Ti/Pd, and under 1.02 for every
asymmetric pair in the table. The conclusion survives with a number
attached instead of a bad argument: **illumination engineering is worth
under 2% on an n/p pair, and negative per incident photon. Masks are for
symmetric devices only.**

> **RETRACTED 2026-09-20 — see Section 6.11.2.** The two numbers above are
> correct. The generalisation from them is not. Enumerating all 21
> asymmetric pairs under this same model gives **14 violations of the
> "under 1.02" claim, the worst being Au/Pd at 22.4**. The sample of two
> was drawn entirely from pairs that *straddle* the 5.4 eV crossover, and
> every straddling pair does satisfy the bound; same-sign pairs need not,
> because max|k|/|N_uniform| diverges as N_uniform → 0. The final sentence
> — "masks are for symmetric devices only" — is therefore **false**: a
> perfect mask gains 8.4× on Au/Pd and more than 2× on five pairs, per
> incident photon. What survives is the *practical* recommendation, on
> different grounds, in Section 6.11.2. The original numbers are left in
> place above rather than rewritten.

### 6.9.6 Result: the kernel as a predicted photocurrent scan

k(x) reproduces the qualitative signature the scanning-photocurrent
literature reports. For a symmetric pair the kernel runs from k(0) = +1.000
through exactly zero at mid-channel to k(L) = −1.000 — **opposite polarity
at the two interfaces**, which is Shimomura *et al.*'s "polarities … are
opposite" and the p–n–p structure Mueller *et al.* map with doping
extending 0.2–0.3 μm into the channel [arXiv:0902.1479]. For Ti/Pt the
kernel never changes sign, running between −1.837 and −1.830: a *flat*
scan, which is the experimental signature distinguishing an n/p pair from
a symmetric one.

### 6.9.7 What this section does not claim

**It is not a prediction of a measured SPCM trace.** Kasırga's review
[arXiv:2509.09390, 2025] stresses that measured scanning-photocurrent
signals are frequently photo-thermal rather than photovoltaic, and this
thesis models drift collection only — no photo-thermoelectric and no
bolometric term, a simplification standing since Section 6.4. The curves
here are the *photovoltaic contribution* to such a trace.

**Optical localisation is impossible at this device's geometry.** The
channel is L = 200 nm; Shimomura *et al.*'s focused spot is ~2 μm, ten
times the whole channel. A diffraction-limited spot illuminates both
contacts at once. The spot-size study is therefore a statement about
required channel length, not a proposal for this device: reaching 90% of
the delta-spot limit needs σ ≲ 0.05 L, so a ~1 μm spot would require a
channel of order 20 μm. The only non-uniform illumination realisable at
200 nm is a *lithographic* mask on the device — which is exactly what
Shimomura *et al.* built, and why they built it that way.

All of Section 6.8's simplifications are inherited unchanged: a linear
ΔW → doping-profile relation with a single λ, drift-only transport with
μ_e = μ_h, and no photogain. Two new ones are added: g(x) weights
absorption linearly (safe at graphene's 2.3%, not for a thick absorber),
and the mask is a pure transmission factor — no near-field scattering,
reflection off the nickel, or plasmonic response at the mask edge.

Figure: `photodetector_nonuniform_illumination.png` — (a) the collection
kernel as a delta-spot scan, showing the symmetric sign reversal and the
flat asymmetric kernel; (b) symmetric devices under uniform light, a
perfect mask and the measured mask; (c) the two accountings of mask
leakage, (1−T)/(1+T) per absorbed photon against (1−T) per incident photon.

## 6.11 The work-function → doping relation is a square root, and two things follow

Sections 6.6 through 6.9 all take the contact-induced doping profile to be
**linear in the work-function offset** dW = W_metal − 5.4 eV, with one decay
length λ = 250 nm for every metal. `total_field()` builds the field as
dW/λ/(1+x/λ)², so every result in this chapter is a statement about dW
rather than about the physical Fermi-level shift, *provided the two are
proportional*. This section removes that assumption. Literature and
derivation: `notes/2026-09-20-nonlinear-work-function-to-doping-relation.md`;
model: `graphene_contact_doping_nonlinear_model.py`.

### 6.11.1 The relation, and where it may be used

Charge leaving the metal enters graphene's *linear* density of states, so
n ∝ E_F², and the electrostatic step across the interface gap is therefore
quadratic in the Fermi-level shift. Equilibrium requires

    ΔW′ = φ + (α/2) φ²  ,        φ ≡ ΔE_F in eV,

which inverts to Khomyakov *et al.*'s Eq. 7 [PRB **79**, 195425 (2009)]:

    ΔE_F(ΔW′) = sgn(ΔW′) · ( √(1 + 2α|ΔW′|) − 1 ) / α ,

with ΔW′ = W_metal − W_graphene − Δ_c the offset *after* the short-range
chemical term. This is linear only as ΔW′ → 0 and asymptotically √ΔW′
beyond: **graphene resists being doped, increasingly, the harder one pushes.**
It also explains a constant this thesis has been hard-coding since
Section 6.8. The p/n crossover is not a property of graphene; it is
W₀(d) = W_graphene + Δ_c(d), which happens to be 4.5 + 0.9 = 5.4 eV *at the
physisorbed separation*.

The constant is derived rather than fitted,

    α = 2e³(d − d₀) / (ε₀ π ℏ² v_F²) = 2.39 eV⁻¹   for d − d₀ = 0.9 Å,

and is good to about 1.4× against the paper's Pt figure — the right size and
sign, not a quantitative calibration — so every result below is also reported
as a sweep over α from 0 to 5 eV⁻¹.

Five validations against exactly known values pass: ΔE_F(0) is bitwise zero;
the relation is bitwise odd; Eq. 7 inverts its own forward equation to
8.9 × 10⁻¹⁶ eV; **at α = 0 the new field function reproduces Section 6.8's
`total_field()` bitwise** on the asymmetric Ti/Pt pair, with the leading
residual confirmed cubic as the O(αΔW²) expansion requires; and the
symmetric-pair zero of Section 6.8.3 survives (1.3 × 10⁻¹⁶).

**Where it may not be used.** Khomyakov *et al.* find the chemisorbed metals
sitting at d_eq ≈ 2.05–2.3 Å, *below* d₀ = 2.4 Å, so Eq. 7's gap-capacitance
term is outside its own regime for them; this is Giovannetti *et al.*'s
stated conclusion that chemisorbed metals are not characterised by work
function alone. Of this thesis's seven metals, **Ti, Ni and Pd are
chemisorbed and Cr has no tabulated separation — four of seven are
excluded**, and the model refuses to extrapolate rather than producing a
number. This is not a technicality: Ti carries the largest |dW| in the table
(−1.07 eV) and is contact A of *both* of this chapter's headline pairs. The
Ti figures below are therefore reported as an indication of direction and
magnitude, not as corrected values.

### 6.11.2 Result: Ti/Pt survives; the Section 6.9 ceiling does not

A prediction was written into the note before the model was run —
"compression, not reversal", on the grounds that a monotone, sign-preserving
map cannot flip a p/n assignment. **It is half right, and the half that is
wrong is the informative half.**

| pair | straddles 5.4 eV? | N linear | N nonlinear | ratio | ceiling lin. | ceiling nonlin. |
|---|---|---|---|---|---|---|
| Ti/Pt | yes | −1.8323 | −1.7418 | 0.951 | 1.0025 | 1.0064 |
| Cr/Pt | yes | −1.8104 | −1.7207 | 0.950 | 1.0033 | 1.0077 |
| Cu/Pt | yes | −1.7855 | −1.6972 | 0.951 | 1.0043 | 1.0092 |
| Ti/Pd | no | −1.5988 | −0.6972 | **0.436** | 1.0159 | **1.4343** |
| Ti/Au | no | −1.5499 | −0.6628 | **0.428** | 1.0203 | **1.5088** |
| Ti/Ni | no | −0.8075 | −0.5705 | 0.707 | 1.2383 | 1.7527 |

For pairs straddling the crossover the prediction holds exactly as written:
**Ti/Pt compresses by 4.9%, from 1.832 to 1.742**, and the α-sweep keeps it
between 0.92 and 1.00 of the linear value across the whole range 0–5 eV⁻¹.
Section 6.8's headline is robust to this correction.

For same-sign pairs it fails badly — Ti/Pd loses 56% — for a reason the
prediction did not anticipate. A monotone map cannot reverse a *sign*, but
it can compress the *ratio* of two offsets, and when the two contacts dope
graphene the same way the response is a near-cancellation whose size is set
by exactly that ratio. Compression of a ratio is amplification of a
near-cancellation's fragility.

That is what exposed the Section 6.9 retraction. Ti/Pd's ceiling moving from
1.016 to 1.434 prompted enumerating all 21 asymmetric pairs **under the
original linear model**, where 14 violate the "under 1.02" claim and Au/Pd
reaches 22.4. The inequality |N[g]| ≤ max|k| that Section 6.9.2 validated is
a normalisation identity and remains true (0 violations in 343 cases); what
was over-generalised is its *tightness*, from a sample of two.

**And therefore so is the conclusion.** Scoring a perfect shadow mask per
*incident* photon — the accounting Section 6.9 itself insisted on:

| pair | N uniform | best masked (per incident) | gain | vs. unmasked Ti/Pt |
|---|---|---|---|---|
| Au/Pd | −0.0447 | 0.3765 | **8.42×** | 0.21× |
| Ti/Cr | −0.1370 | 0.4652 | 3.40× | 0.25× |
| Ni/Au | −0.1218 | 0.4089 | 3.36× | 0.22× |
| Cr/Cu | −0.1417 | 0.4595 | 3.24× | 0.25× |
| Ni/Pd | −0.1665 | 0.4140 | 2.49× | 0.23× |
| Ti/Pt | −1.8323 | 0.9165 | 0.50× | 0.50× |

> **Annotation added 2026-09-22 (Section 6.13.7).** The count *five* holds
> only at exactly `w_cross = 5.4 eV`. Re-scored across the scalar-crossover
> band `|δ| ≤ 0.2 eV` it runs **3 → 4 → 5 → 5 → 6**: at `δ = −0.2 eV` Ni/Au
> and Ni/Pd fall below 2×, and at `δ = +0.2 eV` Ti/Cu rises above it. The
> numbers in this table are **not** retracted — they are correct at `δ = 0` —
> but the *count* is, because it is a threshold evaluated at one point of a
> parameter both source papers give as "~5.4 eV". The qualitative claim below
> survives at every offset.

Five pairs gain more than 2×, so **"masks are for symmetric devices only" is
false**. The correct statement is that masks help whenever |N_uniform| is
small, which covers symmetric pairs (where it is zero) *and*
weakly-asymmetric same-sign pairs.

The *design* recommendation nevertheless survives, on different grounds than
those given on 2026-09-19. The best masked device in the table reaches 0.377
per incident photon; unmasked Ti/Pt reaches 1.832, a factor of 4.9 better. **A
large relative gain on a small number is still a small number**, and Section
6.9 conflated relative gain with absolute performance. Choose the metals
first; illuminate uniformly.

### 6.11.3 What this section does not claim

1. It does **not** supply corrected values for any Ti-, Ni-, Pd- or
   Cr-containing pair. Four of the seven metals fall outside Eq. 7's regime
   (Section 6.11.1), and the table above is run through the excluded metals
   only to show the *direction and size* of the correction.
2. α is derived, not calibrated: 1.4× against the one cross-check available.
   The α-sweep exists precisely so that no conclusion rests on its value,
   and none of the conclusions above changes over 0–5 eV⁻¹.
3. Δ_c is treated as a constant 0.9 eV, i.e. the crossover is still held at
   5.4 eV for all metals. Khomyakov *et al.* give Δ_c(d), so a metal at a
   different separation has a different crossover — a **per-metal w_cross**,
   which this chapter does not yet implement and which would move every
   dW in the table.
4. λ = 250 nm is still one value for all metals. The nonlinearity treated
   here is in the profile's *magnitude*; its spatial *extent* is untouched,
   and Khomyakov *et al.* 2010's screening nonlinearity is a separate
   question this chapter has not opened.
5. No per-metal ΔE_F from the PRB's Table I is used anywhere, because two
   extraction passes over that PDF disagreed with each other; see Section 3
   of the 2026-09-20 note.

Figure: `contact_doping_nonlinear_relation.png` — (a) ΔE_F(ΔW′) against the
linear model, with Ti, Pd and Pt marked; (b) the Ti/Pt collection kernel
under both models; (c) the α-sweep, showing compression without reversal.

## 6.12 The 5.4 eV crossover is a separation, not a constant: a per-metal p/n crossover

Sections 6.8–6.11 all place the p/n crossover at a single
`W_CROSS_CHEM = 5.4 eV` shared by every metal. That number is well sourced
[Giovannetti *et al.*, PRL **101**, 026803 (2008)] but it is a *derived*
quantity, and Khomyakov *et al.* [PRB **79**, 195425 (2009)] say what it is
derived from: the interface potential step carries a short-range **chemical**
term Δ_c(d) alongside the charge-transfer term, and charge neutrality puts
the crossover at

  w_cross(d) = W_G + Δ_c(d)     (6.12.1)

evaluated at the metal's own equilibrium separation. 5.4 eV is (6.12.1) at
the *physisorbed* separation d ≈ 3.3 Å, where Δ_c ≈ 0.9 eV. A metal that sits
elsewhere has a different crossover and therefore a different ΔW.

This is a different kind of correction from Section 6.11's. The square-root
relation **compressed** the existing offsets about a fixed zero. A per-metal
crossover **moves the zero**, separately for each metal, and a fixed shift in
the zero is felt *most* by the metals with the smallest |ΔW| — the opposite
of the usual intuition that small offsets are the safe ones.

### 6.12.1 The model, and why it has exactly one free parameter

Khomyakov *et al.*'s Eq. 4 parametrises Δ_c(d) = e^(−γd)(a₀ + a₁d + a₂d²).
Two extraction attempts on the ar5iv rendering returned the equation in
symbolic form with the four fitted constants absent — the second extraction
failure this thesis has recorded on this same paper. Rather than guess them,
Δ_c is written as a one-parameter family pinned to the one value every
extraction pass agrees on, Δ_c(3.3 Å) ≈ 0.9 eV:

  Δ_c(d ; ℓ) = 0.9 eV · exp(−(d − 3.3 Å)/ℓ)     (6.12.2)

with ℓ the short-range decay length, swept over 0.3–1.5 Å. Two exactness
properties make this defensible rather than merely convenient: Δ_c(3.3 Å) =
0.9 eV bitwise for every ℓ, and ℓ → ∞ returns 0.9 eV at every d, so the whole
per-metal apparatus collapses **bitwise** onto the 5.4 eV convention of
Sections 6.8–6.11. Equation (6.12.2) is a *simplification* of Eq. 4 — it drops
the polynomial prefactor — and Section 6.12.4 is where that bill comes due.

### 6.12.2 Validation against four exactly known values

Because `net_response()` takes one crossover for both contacts and forms ΔW
internally, it cannot express a per-metal crossover; a new entry point one
level lower was unavoidable. The first check is what licenses it.

1. **The ΔW-level entry point equals the Section 6.8 model** on **49/49
   ordered pairs, bitwise**, worst difference 0.000e+00. It is the same
   arithmetic with the subtraction moved outward, not a re-implementation.
2. **ℓ → ∞ reproduces the flat model bitwise**: Δ_c(d, ∞) = 0.9 eV bitwise at
   every tabulated separation, w_cross = W_G + 0.9 bitwise for every metal,
   and 49/49 pair responses bitwise identical. (4.5 + 0.9 also turns out to be
   the same double as the literal 5.4 — checked, not assumed.)
3. **The anchor survives the parametrisation**: Δ_c(3.3 Å, ℓ) = 0.9 eV
   bitwise for 61/61 swept values of ℓ.
4. **The exact symmetries of Section 6.8 survive**: a symmetric pair still
   gives N = 0 (worst 6.6 × 10⁻¹⁷ over 9 cases), and charge conjugation
   N(−ΔW) = −N(ΔW) holds to **exactly** 0.000e+00 over 27 cases. Stated at
   the ΔW level this is sharper than Section 6.8.3's version, which had to
   reflect W through a crossover and so could only test one crossover.

### 6.12.3 Result: three of seven metals, and a 17% shift on Cu

The model admits only the physisorbed metals — Cu, Au and Pt. Ti, Ni and Pd
are chemisorbed (Section 6.12.4); Cr has no tabulated separation and was not
guessed, the same stance Section 6.11 took.

| metal | W (eV) | d_eq (Å) | ΔW flat | ΔW range over ℓ | max shift |
|---|---|---|---|---|---|
| Cu | 4.65 | 3.26 | −0.7500 | −0.8784 … −0.7743 | **17.12 %** |
| Au | 5.10 | 3.31 | −0.3000 | −0.2940 … −0.2705 | 9.84 % |
| Pt | 5.65 | 3.30 | +0.2500 | +0.2500 … +0.2500 | 0.00 % |

No metal changes doping type anywhere in the range: Cu and Au stay n-doping,
Pt stays p-doping. The crossover moves, but never past a metal.

The size of the Cu shift is the result. **17% is more than three times the
4.9% compression Section 6.11 produced on the Ti/Pt headline**, and Section
6.11 treated ~5% as the reassuring outcome. The session pre-registered a
prediction of ≤ 10% — calibrated on exactly that precedent — and it was
**falsified**. Bisection puts the boundary at ℓ = 0.50 Å: the prediction holds
for every ℓ ≥ 0.50 Å and fails below, so it is wrong over roughly the
shortest sixth of the declared range. Declaring that range generously is what
made the failure visible; a narrower sweep would have returned agreement and
been worth less.

**Pt's exact 0.00% is an artefact and must not be read as robustness.**
d_eq(Pt) = 3.30 Å coincides with the anchor, so Pt is pinned by construction
for every ℓ. Every per-metal number in this section is a shift *relative to
Pt*, not an absolute one — a real limitation of a one-anchor parametrisation.

### 6.12.4 Result: the anchored exponential fails on the chemisorbed metals, for a reason unrelated to Section 6.11's

Extrapolating (6.12.2) inward to the chemisorbed separations gives crossover
work functions of **6.25 – 62.6 eV** (Ni 2.05 Å, Ti 2.10 Å, Pd 2.30 Å), every
one of them above 5.9 eV, the highest elemental work function there is. The
implied |ΔW| runs 1.13 – 57.5 eV against a largest *observed* Fermi-level
shift of 0.5 eV. The extrapolation therefore predicts that every chemisorbed
metal is maximally n-doping irrespective of its own work function, which is
false, and it does so across the **entire** admissible range of ℓ rather than
only at its short end.

This is a failure of (6.12.2), not of Khomyakov *et al.* Their polynomial
prefactor is precisely what allows Δ_c to turn over instead of running away at
short d, and Section 6.12.1 could not recover it. What makes the failure worth
reporting is that it reaches Section 6.11's conclusion **by an independent
route**: the gap-capacitance term fails on Ti, Ni and Pd because d_eq < d₀,
and the chemical term fails on the same three because Δ_c diverges there. Two
distinct parts of one framework failing on the same three metals for unrelated
reasons is a structural statement about which contacts graphene's physics can
be described by a work function at all — and Ti is in the failing set while
being contact A of both of this chapter's headline pairs.

### 6.12.5 Result: a 65-fold difference between same-sign and straddling pairs

All six ordered asymmetric pairs the model admits, enumerated:

| pair | kind | N flat | N range over ℓ | pair change | ΔW change | amplification |
|---|---|---|---|---|---|---|
| Cu/Au | same-sign | −0.6549 | −1.1653 … −0.6922 | **77.94 %** | 17.12 % | **4.55** |
| Au/Cu | same-sign | +0.6549 | +0.6922 … +1.1653 | **77.94 %** | 17.12 % | **4.55** |
| Cu/Pt | straddling | −1.7855 | −1.8072 … −1.7900 | 1.21 % | 17.12 % | 0.07 |
| Pt/Cu | straddling | +1.7855 | +1.7900 … +1.8072 | 1.21 % | 17.12 % | 0.07 |
| Au/Pt | straddling | −1.6402 | −1.6368 … −1.6228 | 1.06 % | 9.84 % | 0.11 |
| Pt/Au | straddling | +1.6402 | +1.6228 … +1.6368 | 1.06 % | 9.84 % | 0.11 |

A 17% shift in one contact's ΔW becomes a **78%** shift in the same-sign
pair's response and a **1.2%** shift in the straddling pair's — a factor of
**65 between the two kinds of pair, from identical inputs**. Section 6.11
inferred this ratio-amplification mechanism from a single same-sign pair under
a perturbation that *compressed* the offsets; finding it again, with the same
sign and comparable magnitude, under a perturbation that *moves their zero*,
is independent evidence that it is a property of near-cancelling pairs rather
than of either model.

The straddling pairs are correspondingly more robust than the per-metal number
suggests: Cu/Pt moves 1.2% while its own Cu contact moves 17%. A
near-cancellation amplifies; a reinforcing pair averages.

### 6.12.6 What this section does and does not change

* **It does not reach the Ti/Pt headline**, because Ti cannot be evaluated.
  Section 6.8.4's central number is now *untested* by two successive
  refinements rather than confirmed by them — a weaker position than Section
  6.11 left it in, not a stronger one.
* **The best evaluable asymmetric pair, Cu/Pt (−1.786), is robust to 1.2%.**
  The chapter's qualitative design rule — a straddling n/p pair is what works
  — survives cleanly.
* **Cu/Au should not be quoted to better than a factor of two** under any
  crossover convention.
* **W_CROSS_CHEM = 5.4 eV remains the right default for physisorbed metals**,
  now with a quantified error bar instead of an implicit claim of exactness:
  ±17% on ΔW for a metal 0.04 Å off the anchor, ±78% on a same-sign pair built
  from such metals. Sections 6.8–6.11 are **not retracted**; their numbers
  stand as committed, and this is the uncertainty to read them with.

### 6.12.7 What this section does not claim

Δ_c is (6.12.2), not Khomyakov *et al.*'s Eq. 4; its fitted constants could
not be obtained and are not guessed. ℓ is unmeasured and everything is
reported as a sweep. The anchor is a single datum at a single separation, so
all shifts are relative to whichever metal sits at it (here Pt, by
coincidence). Nothing here touches the chemisorbed metals, the 0.12 eV
measured-step discrepancy of Section 6.8.7, or the absence of a
photo-thermoelectric term. The per-metal crossover has **not** been propagated
back through Sections 6.9–6.11; those sections' numbers remain at the flat
5.4 eV convention and are to be read with the error bar above.

Figure: `per_metal_crossover.png` — (a) w_cross(d) for three decay lengths on
a log axis, with the seven metals placed and the chemisorbed region shaded as
out-of-regime; (b) the change in zero-bias net response against ℓ for all six
admitted pairs, separating the two same-sign curves from the four straddling
ones.

## 6.13 How far can the crossover move? A sensitivity band for Sections 6.8–6.12

Section 6.12 replaced the single 5.4 eV crossover with a per-metal one, and
in doing so raised a question it could not answer: **what does the rest of
the chapter look like if 5.4 eV is simply wrong by some amount?**

That question has to be asked separately, because the per-metal machinery
**cannot reach the pairs that matter**. It needs the metal–graphene
separation `d_eq`, which Khomyakov *et al.* tabulate for three of this
thesis's seven metals (Cu, Au, Pt). Ti, Ni and Pd are chemisorbed, where
Section 6.12.4 showed the anchored exponential diverges; Cr has no tabulated
separation at all. And every one of the five headline pairs of Section
6.11.2's retraction table — Au/Pd, Ti/Cr, Ni/Au, Cr/Cu, Ni/Pd — contains at
least one metal outside that set of three.

### 6.13.1 The perturbation, and why a scalar one is honest

Both papers the 5.4 eV comes from write it with a tilde:

> "the crossover from p-type to n-type doping occurs for a metal work
> function of **~5.4 eV**"
> — Giovannetti *et al.*, PRL **101**, 026803 (2008); repeated verbatim in
> Khomyakov *et al.*, PRB **79**, 195425 (2009).

Sections 6.8–6.12 have read that tilde as three significant figures. This
section puts an unknown scalar offset `δ` on the crossover, `w_cross = 5.4 +
δ`, and sweeps it. The offset has the *wrong shape* — it is not per-metal —
but the *right reach*: it needs no `d_eq`, so it applies to all seven metals.
The two treatments are complementary and neither supersedes the other.

It is also algebraically identical to a **common error in every entry of the
work-function table**, since only `dW = W − w_cross` enters. That is not a
hypothetical: Ni is entered at 5.04 eV from a range commonly quoted as
4.9–5.35 eV, and the 0.12 eV potential step Mueller *et al.* measure has been
an unreconciled discrepancy in this thesis since Section 6.8. A nominal band
of `|δ| ≤ 0.2 eV` is adopted, with a ±0.5 eV outer sweep; no result below
depends on that choice, because every result is reported as a function of δ.

### 6.13.2 Validation against five exactly known values

Four were planned; the fifth was found while trying to score a prediction and
turned out to be the section's main result.

| # | check | outcome |
|---|---|---|
| 1 | shift invariance: `N(W; 5.4+δ) == N(W−δ; 5.4)` | 215/231 bitwise, worst **1.78e-15** |
| 2 | `δ = 0` reproduces the Section 6.8 model | **bitwise**, 21/21 pairs |
| 3 | symmetric pair → `N == 0` at every δ | worst **1.3e-16**, 77 cases |
| 4 | charge conjugation `N(−dW) == −N(dW)` at every δ | **exactly 0.000e+00**, 539/539 |
| 5 | a straddling pair's `⟨\|dW\|⟩` is independent of `w_cross` | **exactly 0.000e+00**, 126 cases |

Check 1 is exact in real arithmetic but not in floating point, and the note
predicted that in advance — along with a bound of 1e-15, which the measured
1.78e-15 **misses by a factor of 1.8**. The bound was the wrong shape: the
right statement is a few ulp of `N` (1.78e-15 is 8 ulp of 0.593), not a fixed
absolute number. It is recorded because a pre-registration practice is worth
nothing if only the comfortable misses are recorded.

### 6.13.3 Result: the amplification ratio of Section 6.12 is not a property of a pair

Section 6.12.5's headline was a **65-fold** difference in how hard a
same-sign pair and a straddling pair amplify an input uncertainty, measured
on three pairs. Repeating that measurement here gives, for all six straddling
pairs, an amplification of **infinity** — not large, exactly undefined.

Check 5 is why. For a pair straddling the crossover,

    ⟨|dW|⟩ = [ (w_cross − W_A) + (W_B − w_cross) ] / 2 = (W_B − W_A) / 2

in which `w_cross` has cancelled identically. Moving the crossover changes
such a pair's response (`S` below is small but nonzero) while changing the
input measure by **exactly nothing**, so the ratio divides by zero.

**The amplification ratio is a property of the (pair, perturbation) couple,
not of the pair.** Section 6.12.5's number should not be quoted as a device
property, and this chapter came close to doing so.

### 6.13.4 Result: the dichotomy is real, at 30× rather than 65×

The statistic that survives is the logarithmic sensitivity
`S = d ln|N| / dδ`, which is well defined for every pair because it does not
divide by the input:

| group | range of `\|S\|` (eV⁻¹) | pairs |
|---|---|---|
| straddling | 0.0077 – 0.0218 | Ti/Pt, Cr/Pt, Cu/Pt, Ni/Pt, Au/Pt, Pd/Pt |
| same-sign | 0.6625 – 2.5628 | the other 15 |

`min(same-sign) / max(straddling) = **30.4**`, with no overlap anywhere in
the table. So the effect Section 6.12 saw on three pairs under a *per-metal*
perturbation reappears on all twenty-one under a *scalar* one — genuinely
independent evidence that it belongs to the pair's sign structure — but the
magnitude comes down from 65× to 30× when the sample is widened from three
pairs to twenty-one. That is the third time in this chapter that widening a
sample has shrunk a claim.

The most sensitive pair is **Au/Pd** (`|S| = 2.56 eV⁻¹`), whose contacts are
0.02 eV apart; the least is Pd/Pt (0.0077). This was predicted in advance and
held.

### 6.13.5 Result: no scalar offset can reverse a response sign

Write `dW_A = m − s` and `dW_B = m + s`, with `s = (W_B − W_A)/2` and
`m = (W_A + W_B)/2 − w_cross`. A scalar offset moves `m` and leaves `s`
**exactly** fixed, so the contact field separates as

    E(x) = m [ f(L−x) − f(x) ] + s [ f(x) + f(L−x) ],    f(x) = λ⁻¹(1+x/λ)⁻²

into an antisymmetric part carrying `m`, which contributes exactly zero to
`N` by the same argument as Validation 3, and a symmetric part carrying `s`,
which the perturbation cannot touch. **`sign(N)` is fixed by which contact
has the higher work function.** A scan over `|δ| ≤ 3 eV` (601 points × 21
pairs = 12621 evaluations) finds zero sign changes; the smallest `|N|`
anywhere is 5.1e-03, approached and never crossed. `m` only rescales `|N|`,
vanishing as `|m| → ∞` — and that limit *is* the near-cancellation of Section
6.11.2, now with an explicit mechanism.

This is a **stronger robustness statement than any in Sections 6.8–6.12**:
the *direction* of a two-terminal photoresponse, which is what a design rule
actually asserts, is immune to the entire "~" on 5.4 eV, and indeed to a
scalar error six times its size.

It is also the opposite of what was predicted. The pre-registered derivation
D5 held that a sign flip occurs whenever `w_cross` moves between a pair's two
work functions, and tabulated 21 flip windows accordingly. The premise is
wrong: `sign(dW_A · dW_B)` is not what sets `sign(N)`.

### 6.13.6 A methodological finding: the five exact checks could not see the bug

The first implementation of the sign-flip scan **reported 21 roots.** It
bisected between `δ = 0` and the midpoint of D5's window without ever
checking that a root was bracketed; with no sign change present, such a loop
walks its lower bound up to its upper bound and returns the **endpoint**.
Every printed root equalled `(W_A + W_B)/2 − 5.4`, every one was physically
plausible, and its headline — *"nearest sign flip: Pd/Pt at −0.0150 eV"* —
is exactly the kind of number that reaches a thesis unchallenged.

All five exact validations passed while it did so, because they validate the
*model* and the fault was in the *analysis wrapped around it*. This is a
third distinct failure class, and the chapter should keep all three separate:

| class | example | caught by |
|---|---|---|
| implementation error in the model | the `inf − inf` transit integral, 2026-09-17 | an exactly known zero |
| a claim quantified on an unrepresentative sample | the `< 1.02` ceiling, Section 6.9 → 6.11 | enumerating all 21 pairs |
| **error in the analysis layer** | this section's unbracketed bisection | **neither** — hand-checking one row |

### 6.13.7 Result: Section 6.11.2's "five pairs" is retracted as a count

Section 6.11.2's table is re-scored at each offset (best perfect mask, per
*incident* photon, divided by `|N_uniform|`):

| pair | δ=−0.20 | δ=−0.10 | δ=0 | δ=+0.10 | δ=+0.20 |
|---|---|---|---|---|---|
| Au/Pd | 2.985 | 5.741 | **8.421** | 11.055 | 13.663 |
| Ti/Cr | 2.735 | 3.066 | 3.396 | 3.727 | 4.056 |
| Ni/Au | 1.534 | **2.442** | 3.357 | 4.277 | 5.195 |
| Cr/Cu | 2.496 | 2.868 | 3.242 | 3.614 | 3.989 |
| Ni/Pd | 1.122 | 1.801 | **2.486** | 3.176 | 3.869 |
| Ti/Cu | 1.337 | 1.513 | 1.689 | 1.865 | **2.041** |

The count above 2× runs **3, 4, 5, 5, 6** across the band. The qualitative
replacement claim — masks help whenever `|N_uniform|` is small — holds at
every offset, and so does the design recommendation, since Au/Pd's best
masked response is 0.377 per incident photon against unmasked Ti/Pt's 1.832.
**What does not survive is the number five**, which is an artefact of
evaluating a threshold at one point of a parameter the literature gives with
a tilde. Section 6.11.2 is annotated in place rather than rewritten.

The Section 6.9.2 ceiling `|N[g]| ≤ max|k|` held at every offset: 0
violations in 210 checks.

### 6.13.8 Result: Ti/Pt is now the most-tested claim in the chapter

Over the whole nominal band Ti/Pt's `|N|` ranges 1.8241–1.8367 — a **0.45%**
worst-case change, against 4.9% under Section 6.11's nonlinearity and
1.2–1.8% under Section 6.12's per-metal crossover. Three independent
refinements have now failed to move it. Whatever the risk to Chapter 6's
headline is, it is not the value of `w_cross`.

### 6.13.9 What this section does and does not change

**Changes:** Section 6.12.5's 65× amplification is restated as 30× and
reclassified as perturbation-dependent rather than a device property; Section
6.11.2's "five pairs" becomes "three to six across the band"; Sections
6.8–6.12 gain a stated scalar-crossover band, which for straddling pairs is
**under 0.5%** and for same-sign pairs reaches **51%** at `|δ| = 0.2 eV`.

**Does not change:** every doping-type assignment, every response sign, the
Ti/Pt headline, the `max|k|` ceiling identity, and the design recommendation
to choose metals before illumination.

### 6.13.10 What this section does not claim

1. It does **not** determine the crossover. It measures sensitivity to an
   unknown offset; the offset remains unknown.
2. A scalar offset is **not** the per-metal structure of Section 6.12, and
   the robustness result of 6.13.5 does not transfer to it. A scalar offset
   cannot change `s = (W_B − W_A)/2`; a *differential* one can. **A per-metal
   crossover difference of 0.02 eV would flip Au/Pd's response sign, where a
   3 eV scalar offset cannot.** That is now the chapter's sharpest open
   question, and Section 6.12 can answer it for exactly three metals.
3. `λ = 250 nm` is still one value for all metals, and the
   photo-thermoelectric term is still absent.
4. The band is a sensitivity band, not an error bar on a measurement: it
   propagates one parameter's uncertainty through a model whose other
   approximations are unchanged.

Figure: `crossover_sensitivity.png` — (a) sensitivity of all 21 pairs on a
log axis, straddling in blue and same-sign in red; (b) `N(δ)` over ±0.5 eV
with the nominal band shaded, showing curves approaching but never crossing
zero; (c) the Section 6.11.2 ">2×" claim crossing its own threshold twice
inside ±0.2 eV.

## 6.14 The differential crossover: an exact sign-flip threshold, and a parity that makes three earlier results corollaries

Section 6.13 closed by naming its own blind spot. Its robustness theorem —
no scalar crossover offset, of any size, can reverse a two-terminal
photoresponse — is a statement about the **common mode** of the crossover
uncertainty and about nothing else. Section 6.12 does not predict one
crossover; it predicts one per metal. Two contacts therefore do not share a
`w_cross` at all, and the **difference** between their crossovers is exactly
the part Section 6.13's theorem cannot touch, because it is the only part
that moves `s`.

This section is that difference.

### 6.14.1 The two-variable decomposition

Give contact A the crossover `w_cross + c − τ/2` and contact B the crossover
`w_cross + c + τ/2`. In the offset variables

    dW_A = m − s,    dW_B = m + s,    s = (W_B − W_A)/2,
    m = (W_A + W_B)/2 − w_cross

this is exactly

    m → m − c,       s → s − τ/2.

`c` is Section 6.13's scalar offset renamed; `τ` is the **differential
crossover offset**. Together `(c, τ)` span the whole two-metal crossover
uncertainty, and they are not equally dangerous: `c` moves only `m`, while
`τ` is the only thing that moves `s`.

This is not a new model. It is the Section 6.8 signed-carrier field entered
at the `dW` level, as Section 6.12 already does, and Validations 1 and 2
below pin it bitwise to both predecessors.

### 6.14.2 The threshold, in closed form

    N = 0 exactly when s = 0, i.e. at   τ* = W_B − W_A                 (6.20)

independently of `c`, of `λ`, of `L`, and of every parameter in the transport
kernel. **The differential crossover offset required to reverse a pair's
photoresponse is exactly that pair's work-function gap.** No model evaluation
is needed to obtain it, and the entire scalar uncertainty band of Section
6.13 is exactly irrelevant to it.

(6.20) was derived and committed before the model existed (note
`2026-09-23-...`, commit `7175db7`) and is scored below as prediction P1.

### 6.14.3 Validation against six exactly known values

Five were planned. The sixth was found while trying to score a prediction —
the same way Section 6.13's fifth was — and it turned out to subsume three
earlier results.

| # | statement | result |
|---|---|---|
| 1 | `(c, τ) = (0, 0)` reproduces the Section 6.8 signed model | **bitwise**, 21 pairs |
| 2 | `τ = 0` reproduces Section 6.13's scalar-offset model | **bitwise**, 231/231 cases |
| 3 | `N = 0` at `τ = τ*`, at every `c` | worst `7.4 × 10⁻¹⁵`, 189 cases |
| 4 | charge conjugation under a *differential* offset | **exactly 0**, 315/315 |
| 5 | a bracketed bisection lands on (6.20) | worst `4.8 × 10⁻¹⁵` eV |
| 6 | parity factorisation (unplanned) | `2.3 × 10⁻¹⁵` |

Validation 5 carries the audit item Section 6.13.6 opened. The bug that
section paid for was an unbracketed bisection that walked its lower bound to
its upper bound and returned the endpoint as a root, 21 times, every one
plausible. Two independent detectors are now active on the same fault:
`bracketed_bisect()` raises on an unbracketed interval — it fired on **21/21**
deliberately root-free intervals, where last session's loop would have
returned a fake root on every one — and the root it does find is compared
against a closed form known in advance. The located threshold drifts by
`9.7 × 10⁻¹⁵` eV as `c` sweeps ±1 eV, which is Section 6.13's entire
uncertainty band expressed as machine noise.

### 6.14.4 Result: the response factorises by parity, and three earlier results collapse into it

Prediction P2 held that `|N|` would be visibly asymmetric about the flip,
since `m ≠ 0` there and the transport kernel is nonlinear. It came back
`0.00%` for all 21 pairs — the signature of an identity, not of a small
number. The identity is a **parity factorisation**:

    N(m, −s) = −N(m, s)      N is ODD  in the differential offset
    N(−m,  s) = +N(m, s)     N is EVEN in the common offset

both to `2.3 × 10⁻¹⁵` over a 13 × 13 grid of `(m, s)`. Three results this
chapter obtained separately are corollaries:

1. **Odd in `s` ⇒ `N(m, 0) = 0` exactly**, which *is* (6.20). The derivation
   in 6.14.2 went through the antisymmetry of the field; the parity is the
   stronger and simpler statement.
2. **Even in `m` ⇒ no common offset of any size can change `sign(N)`.** That
   is Section 6.13.5, which reached the same conclusion by scanning 12621
   points over ±3 eV. A scan is evidence over the interval scanned; a parity
   is not restricted to an interval.
3. **Odd in `s` ⇒ `|N|` depends on `s` only through `|s|`**, so `|N|` *must*
   be symmetric about the flip. P2 and P4 are therefore not near misses:
   they are **excluded**. The pre-registered note assumed the asymmetry
   without noticing it was assuming a symmetry away.

The device statement worth carrying forward: **at zero bias the response
magnitude is blind to the direction of a differential crossover error, while
its sign is entirely determined by it.** A fabrication run whose contact
chemistry is uncertain in an unknown direction has a predictable
signal *size* and an unpredictable signal *polarity* — which is the worse of
the two failure modes for a differential readout.

A note on where the parity stops being exact. Even-in-`m` is exact bitwise at
every bias, because the `−m` field is the literal spatial mirror of the `+m`
field and the grid `linspace(0, L, n)` is symmetric, so the mirror is exact
*on the grid*. Odd-in-`s` is exact at zero bias and departs at finite bias by
exactly `2/(n−1)`: measured at `n = 501, 1001, 2001, 4001, 8001` the residual
is `4.0×10⁻³, 2.0×10⁻³, 1.0×10⁻³, 5.0×10⁻⁴, 2.5×10⁻⁴`, i.e. residual `×
(n−1) = 2.0000` across a sixteen-fold range. That is one trapezoid cell — the
weight is 2 because the collection label jumps from `+1` to `−1` across it —
converging as `1/n`. **A residual that scales with the grid is the grid; a
residual that does not is physics.** This chapter has now needed that
distinction three times.

### 6.14.5 Result: the pairs most at risk are exactly the pairs the model cannot reach

Because (6.20) needs only `METAL_WORK_FUNCTIONS`, the margin table is exact
for all 21 pairs, whether or not Section 6.12 can say anything about them.

| pair | `τ*` (eV) | modelled `τ` range (eV) | margin | reachable? |
|---|---|---|---|---|
| Au/Pd | **0.020** | — | — | no (Pd chemisorbed) |
| Ni/Au | **0.060** | — | — | no (Ni chemisorbed) |
| Ni/Pd | **0.080** | — | — | no (both) |
| Cr/Cu | **0.150** | — | — | no (Cr has no `d_eq`) |
| Ti/Cr | 0.170 | — | — | no |
| Cu/Au | 0.450 | −0.158 … −0.030 | **2.85×** | yes |
| Au/Pt | 0.550 | +0.006 … +0.030 | 18.6× | yes |
| Cu/Pt | 1.000 | −0.128 … −0.024 | 7.79× | yes |

(Full 21-row table in `differential_crossover_output.txt`.)

Two things follow, and they point in opposite directions.

**The three pairs the model can reach keep their sign.** The worst margin is
Cu/Au at **2.85×** — the modelled differential offset would have to be nearly
three times its largest value across the whole `ℓ ∈ [0.3, 1.5] Å` sweep
before Cu/Au reverses.

That is survival, not comfort, and the right comparison makes the point
sharper than any ratio would. In the **common** channel there is no threshold
at all: by the even-in-`m` parity of 6.14.4, no offset of any magnitude
reverses a sign, so the margin is not large, it is undefined-because-infinite.
In the **differential** channel the margin is finite for every pair, and for
Cu/Au it is 2.85 — within a factor of three of the model's own largest
estimate, on the same physics and the same `ℓ` sweep. The two channels are
not quantitatively different; they are **qualitatively** different, and the
chapter's sign claims were previously defended only against the harmless one.
P3 predicted Cu/Au worst with a ratio in `[2, 4]`: **PASS**.

**The four smallest margins in the whole table belong to pairs Section 6.12
refuses.** Au/Pd, Ni/Au, Ni/Pd and Cr/Cu all sit below the largest
differential offset the model produces for the metals it *does* admit
(0.158 eV), and every one of them contains Ni, Pd or Cr — chemisorbed, or
without a tabulated separation. This is not a coincidence of this table; it
is structural. Chemisorbed metals sit ~1.2 Å below the physisorbed anchor,
which is *both* why Section 6.12.4's anchored exponential diverges on them
*and* why their crossovers should differ most from everyone else's. **The
model's reach and the thesis's risk are anti-correlated by construction.**

### 6.14.6 The conditional, stated as a conditional

**If** Ni and Pd have differential crossover offsets of the same order as the
physisorbed metals do (0.006–0.158 eV), **then** Au/Pd's response sign —
`τ* = 0.020` eV — is not merely uncertain but more likely wrong than right.

The antecedent is untested and Section 6.12 cannot test it. What can be said
is the direction of the error: 6.12.4 gives reason to expect the chemisorbed
differential offsets to be **larger**, not smaller. That strengthens the
conditional and is emphatically not licence to drop it. Au/Pd is stated here
as a **threshold that no available model can reach**, which is weaker than a
prediction of reversal and is the strongest thing the evidence supports.

Au/Pd is also, by Section 6.12.5, the most *sensitive* pair in the chapter,
and by Section 6.11.2 one of the five in the retraction table. Three
independent lines now converge on the same pair being the one this thesis
knows least about.

### 6.14.7 What this section does and does not change

**Changes.**
- Section 6.13.5's sign-robustness result is **strengthened in kind**: it was
  a 12621-point scan and is now a parity identity, valid without an interval.
- Section 6.13.10's closing remark — *"a per-metal crossover difference of
  0.02 eV would flip Au/Pd"* — is confirmed **exactly**, and generalised: the
  flipping offset is `W_B − W_A` for every pair.
- The chapter now has an explicit statement of which of its sign claims are
  safe (all of them against common offsets; three of twenty-one against
  differential ones at the modelled magnitudes) and which are unevaluable
  (the four smallest margins).

**Does not change.** No number in Sections 6.8–6.13 is retracted or altered.
Predictions P2 and P4 of this session's own note are falsified; they were
predictions about this model, not results of an earlier section.

### 6.14.8 What this section does not claim

1. It does **not** claim Au/Pd, Ni/Au or Ni/Pd are sign-reversed. It claims
   their thresholds are small and unreachable by the available model.
2. `τ_model` is **not** a measurement. It is Section 6.12's one-parameter
   anchored exponential with a swept decay length and no fitted prefactor;
   6.12.7's caveats apply unchanged.
3. The parity identities are established **numerically to machine precision**
   over the sampled grids, and analytically only in the reflection argument
   of 6.14.4. That is stronger than a plausible-range check and weaker than
   a proof for arbitrary transport kernels.
4. Zero bias is where odd-in-`s` was verified exactly. Under bias the grid
   study shows the parity surviving in the continuum, but every *device*
   statement above is a zero-bias statement.

Figure: `differential_crossover.png` — (a) `N(τ)` through the flip for the
three model-reachable pairs plus Au/Pd, with `τ*` marked; (b) all 21
thresholds ranked, greyed where Section 6.12 cannot reach the pair, against
the largest modelled differential offset; (c) the located threshold minus
`τ*` in units of `10⁻¹⁵` eV as `c` sweeps ±1 eV.

## 6.10 Further follow-on work

- ~~Model both contacts self-consistently~~ -- **done 2026-09-17,
  Section 6.7**, and it reversed Section 6.6's metal ranking for a
  symmetric device.
- Replace Section 6.7's |W_metal - W_graphene| magnitude convention with
  a signed, carrier-resolved treatment that distinguishes n-type
  (Ti, Cu) from p-type (Pt, Pd, Au) contacts and tracks electrons and
  holes separately. This is now the single most consequential open item
  in this chapter: it could invert Section 6.7's Result 2 ordering by
  making an n/p pair such as Ti/Pt add rather than partially cancel.
- ~~Model non-uniform illumination (a generation weight g(x))~~ --
  **done 2026-09-19, Section 6.9**, against Shimomura et al.'s shadow-mask
  device. It required no new transport physics (illumination is a weight on
  a fixed kernel) and it produced the chapter's first *bound*, max|k|,
  rather than another ranking.
- ~~Replace the linear dW -> doping profile relation~~ -- **done
  2026-09-20, Section 6.11**, and it both confirmed Section 6.8's Ti/Pt
  headline (4.9% compression) and retracted Section 6.9's ceiling claim and
  its "masks are for symmetric devices only" conclusion.
- **Implement a per-metal crossover w_cross = W_graphene + Delta_c(d)**, now
  the top open item in this chapter. Section 6.11.3 item 3: holding the
  crossover at 5.4 eV for all seven metals is the last unexamined constant
  in the chain, and unlike the magnitude nonlinearity it would move every
  dW in the table, not just compress them.
- **Find a description of Ti, Ni and Pd that does not go through their work
  function.** Section 6.11.1 excludes four of seven metals from Eq. 7, and
  Ti is contact A of both headline pairs -- so the chapter's central numbers
  now rest on the one assumption the literature most explicitly disowns.
- Reconcile Chapter 4's contact-metal recommendation with Section 6.7's:
  the two chapters now point in opposite directions for the same metals,
  and Chapter 7 (discussion/outlook) should resolve this explicitly
  rather than leaving the reader to notice it.
- Integrate Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery now used in both Section
  6.6 and Chapter 4 (`graphene_contact_doping_model.py`,
  `graphene_edge_contact_model.py`) -- a real device's plasmonic hot-spot
  and contact depletion region sit at related but distinct locations,
  which Section 6.5's spatially-uniform EQE multiplier does not capture.
- Model the photo-bolometric mechanism (dominant in the telecom bowtie/
  waveguide device reviewed in Section 6.5) as its own responsivity
  model -- dR/dT of graphene's resistance vs. absorbed power under bias
  -- rather than treating photogating as the only gain mechanism this
  thesis represents.
- Find or derive a measured (rather than assumed) plasmon resonance
  quality factor/FWHM for either Section 6.5 design, to replace the
  assumed Q = 6-8 values with literature-anchored ones.
- Revisit Park, Ahn et al.'s scanning-photocurrent-microscopy sign-map
  result (PubMed 19326919), which WebFetch could not retrieve this
  session (HTTP 429), as a more direct cross-check of Section 6.6's
  two-contact sign-reversal reasoning than the two sources it currently
  relies on (Xia et al. 2009, Weiss & Duan 2013).

## References

See `notes/2026-08-24-photodetector-responsivity.md` for the
full literature review and citation list supporting Sections 6.1-6.4,
`notes/2026-09-06-plasmonic-enhancement-graphene-photodetectors.md` for
Section 6.5's plasmonic-enhancement literature review and citations, and
`notes/2026-09-07-spatial-photocarrier-collection-model.md` for Section
6.6's spatially resolved collection-model literature review and design
derivation, and
`notes/2026-09-18-signed-carrier-resolved-contact-fields.md` for Section
6.8's signed carrier-resolved treatment, the 5.4 eV n/p crossover
[Giovannetti et al., PRL 101, 026803 (2008)] and the Mueller/Xia/Avouris
Pd-Ti asymmetric-finger experiment, and
`notes/2026-09-17-two-contact-self-consistent-collection.md` for Section
6.7's two-contact literature basis (Weiss & Duan 2013; Shimomura et al.,
*Carbon Trends* 5, 100100 (2021)), its two validations, and the
simplification that could overturn its Result 2, and
`notes/2026-09-19-non-uniform-illumination-and-the-collection-kernel.md`
for Section 6.9's generation-weight treatment, the Shimomura *et al.*
citation correction and the mask-leakage number, the SPCM references
[Kasırga, arXiv:2509.09390 (2025); Mueller *et al.*, *Phys. Rev. B* **79**,
245430 (2009), arXiv:0902.1479], and a corrigendum recording two claims
that session made and then falsified against its own model.
