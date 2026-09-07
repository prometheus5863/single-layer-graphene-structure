# Spatially resolved photocarrier collection near a doped contact

**Date:** 2026-09-07
**Purpose:** literature review and design notes supporting a new
metal-dependent, spatially resolved collection-efficiency model for
`graphene_photodetector_model.py`, closing the Section 6.6 follow-on
item flagged since the chapter was first drafted (2026-08-24): *"Replace
the lumped EQE_bare = 0.15% parameter with a spatially resolved
diffusion-length collection model once the Chapter 4 spatially-resolved
contact-doping-profile follow-on... is implemented"* — that follow-on
(`graphene_contact_doping_model.py`) has existed since 2026-08-26, so
this session integrates it into the photodetector side as originally
planned.

## 1. Why the lumped EQE_bare parameter is unsatisfying

`graphene_photodetector_model.py` (Section 6.4) uses a single literature
midpoint, EQE_bare = 0.15%, applied identically regardless of contact
metal. But the physical mechanism behind the "collection bottleneck"
described in that section's own docstring — only photocarriers generated
within a built-in-field region near a contact are separated and
collected before graphene's ~1 ps carrier lifetime lets them recombine —
is explicitly metal-work-function-dependent: the built-in field comes
from the same contact-induced doping/Fermi-level-shift physics already
modeled quantitatively in `graphene_contact_doping_model.py` for Chapter
4's contact-resistance work, but that spatial profile has never been fed
into the photodetector side. This session closes that gap.

## 2. Literature

**Xia, Mueller, Lin, Valdes-Garcia & Avouris, "Ultrafast graphene
photodetector," *Nature Nanotechnology* 4, 839-843 (2009),
https://doi.org/10.1038/nnano.2009.292** (fetched successfully). This is
the original demonstration that a graphene MGM device with two
*different*-work-function metal contacts (Ti and Pd) generates photocurrent
at zero external bias, with photoresponse intensity-modulation
measurements showing no bandwidth degradation up to 40 GHz and an
estimated intrinsic bandwidth that "may exceed 500 GHz." The
zero-bias operation is only possible because the two contacts' built-in
fields are unequal (asymmetric metals) — this is the direct experimental
demonstration that contact-metal work function sets a real, usable
built-in field in the channel, motivating this session's model.

**Weiss & Duan, "Building potential for graphene photodetectors,"
*NPG Asia Materials* 5, e74 (2013), https://doi.org/10.1038/am.2013.64**
(fetched successfully). States the mechanism directly: *"Contacting
graphene with metals of variable work function changes its Fermi level
and creates a potential offset near graphene-metal contact that
facilitate the separation and transport of photocarriers."* Also makes
an important, honesty-relevant point this model must respect: *"symmetric
metal-graphene-metal structures create equal opposing potentials,
resulting in zero net photocurrent"* at zero bias — i.e., in a
two-identical-metal device the two contacts' built-in fields point in
mirror-image directions and their net contribution to the *terminal*
photocurrent partially cancels unless something breaks the symmetry
(unequal illumination position, asymmetric metals, or, as in this
thesis's existing bare-device model, an externally applied bias that
already picks out a preferred carrier-sweep direction). This paper does
not give a quantitative built-in-region length scale or a work-function-
resolved responsivity number, so it is used here for mechanism and
caveats, not calibration.

**Mueller, Xia & Avouris, "Graphene photodetectors for high-speed
optical communications," *Nature Photonics* 4, 297-301 (2010)** —
already cited in `notes/2026-08-24-photodetector-responsivity.md`
(Section 6.2) for two numbers reused directly here: the ~1 ps graphene
photocarrier lifetime, and a built-in-field collection-region estimate of
**~100-200 nm**. This 100-200 nm figure is a useful independent
cross-check against `graphene_contact_doping_model.py`'s
`lambda_decay = 250 nm` (sourced there from Khomyakov et al. 2010's
DFT-derived doping-decay length) — same order of magnitude, same
physical picture (contact-induced doping extending a few hundred nm into
the channel), different literature source. This session reuses
`lambda_decay` as-is rather than introducing a second, competing length
scale, and notes the ~100-200 nm vs. 250 nm difference explicitly instead
of quietly picking whichever number is more convenient.

**Access failure this session:** Park, Ahn et al., "Imaging of
photocurrent generation and collection in single-layer graphene" (the
scanning-photocurrent-microscopy paper that directly maps photocurrent
*sign reversal* between a graphene device's two contacts — exactly the
phenomenon Section 3 below reasons about qualitatively) was located via
WebSearch (PubMed ID 19326919) but the WebFetch attempt on
`pubmed.ncbi.nlm.nih.gov/19326919/` returned HTTP 429 (rate-limited) and
was not retried. The sign-reversal reasoning below is therefore built
from the two sources above (Xia 2009's asymmetric-contact zero-bias
device, Weiss & Duan 2013's explicit symmetric-contact cancellation
statement) plus standard drift-transport reasoning, not from a direct
citation of the SPCM sign-map result itself. Revisiting that PubMed
abstract (or the underlying Nano Letters paper) is a good candidate for
a future session if a more quantitative cross-check is wanted.

## 3. Model design

**Two contacts, one asymmetry already present.** The existing bare-device
model (`graphene_photodetector_model.py`) already includes a nonzero
external bias V_bias = 0.1 V across L_channel = 200 nm (reused from
`graphene_fet_model.py`), giving a uniform drift field E_bias = V_bias /
L_channel = 5.0e5 V/m that sweeps electrons and holes in opposite,
fixed directions everywhere in the channel. Superposed on this uniform
field, each contact's own doping-induced field points a definite
direction determined by that contact's metal (Section 2's mechanism).
For two *identical* contacts, symmetry means one contact's local doping
field reinforces the bias-driven sweep of its arriving carrier species
while the other's opposes it by the same amount (Weiss & Duan's
"equal and opposing" statement, generalized from the zero-bias case to
the biased case) — this is the compact-model explanation for the
well-known experimental result (Xia 2009 and, per the access-failure
note above, the SPCM literature) that graphene photocurrent maps show
opposite sign near the two contacts of one device.

**Scope decision.** Modeling both contacts self-consistently (one
reinforcing, one opposing, then solving for the net terminal
photocurrent) is the physically complete picture, but doing it honestly
needs a carrier-type-resolved, two-boundary drift-diffusion treatment
that goes beyond this thesis's compact-model scope established elsewhere
(e.g. `graphene_fet_model.py` and `graphene_interconnect_model.py` are
both explicitly *not* TCAD-grade solvers). This session instead models
the single **reinforcing** contact only — i.e., the collection-efficiency
enhancement available at whichever contact happens to have its doping
field aligned with the bias-driven sweep direction for a given carrier
type — and states plainly that the opposing contact's partially
cancelling contribution is not modeled. This is flagged as an explicit
open item in the updated Section 6.6 rather than silently assumed away,
following this repo's established practice (e.g. the Section 4.7
negative-R_extra finding, the Section 4.8 small-hole-diameter gap).

**Field and collection-probability model.**
Reusing `graphene_contact_doping_model.py`'s `lambda_decay` and each
metal's work-function offset dV = -(W_metal - W_graphene):

    E_doping(x) = |dV| / lambda_decay / (1 + x/lambda_decay)^2   (V/m)

(the field implied by differentiating the same saturating profile,
f(x) = 1/(1 + x/lambda_decay), already used for the *density* profile in
`graphene_contact_doping_model.doping_profile()` — reusing one profile
shape for both quantities is a modeling simplification, stated here
rather than hidden: a fully self-consistent treatment would derive both
n(x) and E(x) from one Poisson-like relation rather than positing the
same shape for each). This is superposed on the uniform bias field:

    E_total(x) = E_bias + E_doping(x)

Local drift velocity v(x) = mu * E_total(x) (mu = 0.4 m^2/(V.s), reused
from `graphene_fet_model.py`), and the transit time for a carrier
generated at position x to reach the contact at x=0 is the cumulative
integral

    t(x) = integral_0^x  dx' / v(x')

(evaluated numerically). The probability the carrier survives to be
collected, given graphene's ~1 ps photocarrier lifetime tau_carrier
(Mueller, Xia & Avouris 2010), is the standard exponential-lifetime
survival probability

    p(x) = exp(-t(x) / tau_carrier)

**Baseline and enhancement factor.** A "bias-field-only" baseline
p_baseline(x) is computed the same way with E_doping set to zero
(equivalent to a Cr contact, whose work function 4.50 eV is within
0.01 eV of graphene's own 4.5 eV, i.e. an already-present zero-doping
control case in `METAL_WORK_FUNCTIONS`, not a separately invented one).
Averaging p(x) over the channel (x in [0, L_channel]) gives a mean
collection efficiency eta_collect for each metal, and the ratio
eta_collect(metal) / eta_collect(baseline) is the modeled **collection
enhancement factor**. Multiplying `EQE_bare` by this factor gives a
metal-resolved modeled EQE.

**What this result can and cannot claim.** The relative, metal-to-metal
trend (larger |W_metal - W_graphene| -> shorter collection time -> higher
survival probability -> higher effective EQE) follows directly from the
mechanism in Section 2 and is the main result of this section. The
*absolute* EQE_bare = 0.15% literature midpoint this factor is applied to
is not attributed to any specific contact metal in the material reviewed
in `notes/2026-08-24-photodetector-responsivity.md`, so the resulting
absolute EQE_model(metal) numbers should be read as "what this compact
model predicts relative to whatever metal underlies that literature
midpoint," not as validated absolute predictions for each named metal.

## 4. Sanity check

By construction, `graphene_photodetector_model.TAU_TRANSIT` (computed
from the same L_channel/mu/V_bias parameters used here for E_bias) comes
out to 1.0 ps — numerically coincident with the assumed tau_carrier here,
since both this thesis's existing bare-device model and the independent
Mueller/Xia/Avouris ~1 ps lifetime describe the same "fast, transit-
time-limited" bare-device regime. This gives a bias-only baseline
collection efficiency of eta_baseline = 1 - exp(-1) ~= 0.632 (order-unity,
as expected for a device whose transit time already matches the carrier
lifetime), rising to ~0.93 for the largest-|dV| metal modeled here (Pt,
1.15 V), i.e. up to a ~1.47x enhancement — a modest, physically plausible
factor (not a multiple-order-of-magnitude effect, since the doping field
only dominates within roughly one `lambda_decay` of the contact and the
channel is comparable in length to that decay length here).

## 5. Web search availability

WebSearch and WebFetch were both available and used this session. Two of
four fetch attempts succeeded (Xia et al. 2009 abstract page, Weiss &
Duan 2013); one PMC article (`pmc.ncbi.nlm.nih.gov/articles/PMC7506932/`)
returned a Google reCAPTCHA interstitial instead of content (consistent
with this repo's prior PMC access failures, e.g. 2026-08-31); one PubMed
fetch (19326919) returned HTTP 429 rate-limiting, noted in Section 2
above rather than retried.
