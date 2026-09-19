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

## 6.9 Further follow-on work

- ~~Model both contacts self-consistently~~ -- **done 2026-09-17,
  Section 6.7**, and it reversed Section 6.6's metal ranking for a
  symmetric device.
- Replace Section 6.7's |W_metal - W_graphene| magnitude convention with
  a signed, carrier-resolved treatment that distinguishes n-type
  (Ti, Cu) from p-type (Pt, Pd, Au) contacts and tracks electrons and
  holes separately. This is now the single most consequential open item
  in this chapter: it could invert Section 6.7's Result 2 ordering by
  making an n/p pair such as Ti/Pt add rather than partially cancel.
- Model non-uniform illumination (a generation weight g(x)) rather than
  Section 6.7's uniform assumption -- Shimomura et al.'s shadow-mask device,
  which masks one of the two graphene/electrode interfaces, is precisely
  a non-uniform-generation experiment and is the natural validation
  target for it.
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
simplification that could overturn its Result 2.
