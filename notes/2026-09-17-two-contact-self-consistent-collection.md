# 2026-09-17 -- Self-consistent two-contact photocarrier collection, and
# why it reverses the previous session's contact-metal recommendation

## 0. What this session closed

`AUTOMATION_LOG.md` has carried this as the first entry of "not yet
covered" since 2026-09-07, flagged there as "the natural next step":

> Self-consistent two-contact collection model (Section 6.6's single
> reinforcing-contact scope limitation)

`graphene_photodetector_collection_model.py` (2026-09-07) models ONE
contact. It superposes that contact's doping field on the uniform bias
field, integrates the drift transit time to the contact at x=0, and
assumes every photocarrier is collected there. Its own docstring and
Section 6.6 of the thesis draft both state the limitation explicitly
rather than hiding it: the second contact's "equal and opposing"
contribution is not modeled.

This session implemented the two-contact model
(`graphene_photodetector_two_contact_model.py`) and the result is not a
small correction. **It reverses the single-contact model's ranking of
contact metals for a symmetric photodetector.** Details in Section 3.

## 1. Literature basis

Web search was available; 3 of 4 fetch attempts succeeded.

**Source 1 -- Weiss & Duan, "Building potential for graphene
photodetectors", NPG Asia Materials 5, e74 (2013).**
https://www.nature.com/articles/am201364
Already cited in this repo since 2026-09-07 for its "equal and
opposing" claim; re-fetched this session for the precise wording.
Confirmed: "symmetric metal-graphene-metal devices generate an equal
positive and negative flow with a net zero photocurrent" -- identical
contacts create symmetric Schottky barriers whose opposing fields
cancel. And the escape route: "Using metals with asymmetric band
structures breaks this equilibrium." The article gives no numerical
built-in-potential or photovoltage values, so it anchors the model's
structure and sign convention but cannot calibrate its magnitude.

**Source 2 -- Shimomura, Imai, Nakagawa, Kawai, Hashimoto, Ideguchi &
Maki, "Graphene photodetectors with asymmetric device structures on
silicon chips", Carbon Trends 5, 100100 (2021).**
*(This note originally attributed the paper to "Suzuki et al., Carbon
Trends 5, 100115"; both the first author and the article number were
wrong. Corrected 2026-09-19 against the publisher page and the Ideguchi
group publication list. The physics quoted below is unchanged.)*
https://www.sciencedirect.com/science/article/pii/S2667056921000778
The experimental counterpart, and the more useful one, because its
entire device strategy exists to defeat the cancellation. In symmetric
two-electrode devices "the polarities of the photovoltages at each
graphene/electrode interface ... are canceled out under macroscopic
light irradiation". Two asymmetry mechanisms are built and measured:
(i) a 50 nm Ni shadow mask over one of the two graphene/Ti-electrode
interfaces (with a 120 nm Al2O3 insulating layer), so only one
interface is illuminated; (ii) unequal contact areas, pairing a
rectangular electrode with a comb-shaped one. Reported specific
detectivity ~1.7e5 cm.Hz^(1/2)/W at 690 nm and ~5.2e4 cm.Hz^(1/2)/W at
1310-1530 nm, under zero bias (photovoltage measurement). The authors
are candid that "it is very difficult to accurately determine the
quantitative performance", and report no photocurrent values -- so
again, structure yes, magnitude calibration no.

Worth noting what these two sources jointly establish and what they do
not. They establish that the cancellation is real, first-order, and the
central design obstacle for zero-bias M-G-M photodetectors. Neither
provides a measured net-response-vs-work-function-difference curve, so
the model below is **not calibrated against experiment** -- it is a
mechanism model whose structure is literature-anchored and whose
internal consistency is checked (Section 2). That distinction is stated
here so a later session does not mistake Section 3's numbers for
validated predictions.

**Fetch failure, recorded rather than worked around:** the NEGF-DFT
study on asymmetric metal contacts to bilayer-graphene quantum dots
(PubMed 34942635) returned a Google reCAPTCHA challenge page instead of
the abstract. This is the same PubMed/PMC reCAPTCHA block logged on
2026-09-07 and 2026-09-05, now on a third separate occasion -- it should
be treated as a standing access limitation of this environment, not
retried hopefully each session. It remains the most promising candidate
for an independent first-principles cross-check of the sign and
magnitude of asymmetric-contact photoresponse.

## 2. The model, and two validations it had to pass

Contact A (metal A) at x=0, contact B (metal B) at x=L. A positive total
field F(x) drifts a photocarrier toward contact A:

    F(x) = E_bias + g_A(x) - g_B(x)
    g_A(x) = |W_A - W_gr| / lam / (1 + x/lam)^2         (sweeps toward A)
    g_B(x) = |W_B - W_gr| / lam / (1 + (L-x)/lam)^2     (sweeps toward B)

g_A is byte-for-byte the `doping_field_magnitude()` of the 2026-09-07
module -- same saturating profile, same `lambda_decay` = 250 nm, same
`METAL_WORK_FUNCTIONS` table. g_B is its mirror about x = L/2. **The
minus sign on g_B is the entire physical content that was missing.**

Note the geometry: L_channel = 200 nm but lambda_decay = 250 nm, so
neither contact's field has decayed appreciably by mid-channel. The two
fields overlap everywhere. This is not a regime where the second contact
can be treated as a perturbation.

Because F can change sign inside the channel, the model does not assume
a destination. A carrier at x0 drifts along F and is collected at A
(signed +1) only if F > 0 on all of [0, x0], or at B (signed -1) only if
F < 0 on all of [x0, L]. If F changes sign along the path the carrier
drifts into a **stagnation point** (F = 0) rather than a contact and is
counted as uncollected. Survival uses the same tau = 1 ps photocarrier
lifetime (Mueller, Xia & Avouris, Nature Photonics 4, 297 (2010)) as
before. The figure of merit is the signed channel average

    N(A,B) = (1/L) integral_0^L s(x) exp(-t(x)/tau) dx,  s in {+1,-1,0}

**Validation 1 -- reduces to the model it supersedes.** Give contact B
graphene's own work function, so g_B == 0 exactly. Then F = E_bias + g_A
> 0 everywhere, every carrier reaches A, and the signed average must
collapse onto `mean_collection_efficiency()` from the single-contact
module. It does, for all seven metals, to a worst relative deviation of
**7.8e-08** (trapezoidal-grid resolution, not physics). This is run as
part of `__main__`, printing both models' numbers side by side, rather
than asserted in a comment.

**Validation 2 -- reproduces the literature's central claim.** At zero
bias, two identical contacts make F(x) = g(x) - g(L-x) antisymmetric
about x = L/2, so the net response must be exactly zero. Measured:
worst |N| over all seven symmetric pairs = **6.6e-17**, i.e. machine
precision. This is the Weiss & Duan / Shimomura cancellation, emerging from
the model rather than being put in by hand.

### 2.1 Validation 2 caught a real bug

Worth recording, because the bug was invisible to Validation 1 and to
every physical-plausibility eyeball check.

The first implementation computed one forward cumulative integral
`cum[i] = integral_0^{x_i} dx/|v|` and obtained the toward-B transit
time by subtraction, `cum[-1] - cum[i]`. That is wrong precisely when it
matters: 1/|v| diverges at a stagnation point, so a cumulative sum that
has crossed a null is infinite from there onward, and the subtraction
gives inf - inf for every point beyond the null. Effect: every carrier
on the far side of a stagnation point was silently marked uncollectable.

Validation 1 passed anyway (with g_B = 0 there is no null, so the
subtraction was never exercised). Validation 2 returned N = +0.46 for
Pt/Pt at zero bias instead of 0 -- a number that is not obviously
absurd, and which could easily have been rationalised as "some residual
asymmetry from the grid" by anyone who wanted the model to work. It was
the *exact* zero demanded by the antisymmetry that made the failure
unambiguous.

Fix: integrate each direction only over the interval the carrier
actually traverses (a separate forward and reverse cumulative sum),
which by construction contains no null for a collected carrier.

The general lesson, and the reason this is in the notes rather than just
the commit message: a validation against a *known exact* value (zero, by
symmetry) is worth far more than a validation against a plausible range.
The single-contact agreement check was the one that felt like the real
test; the symmetry check was the one that found the bug.

## 3. Results, and the reversal

### 3.1 Symmetric device at the working bias (V_bias = 0.1 V)

| metal (both contacts) | W (eV) | single-contact (2026-09-07) | symmetric two-contact | ratio |
|---|---|---|---|---|
| Ti | 4.33 | 0.776 | 0.672 | 0.87 |
| Cr | 4.50 | 0.632 | 0.632 | 1.00 |
| Cu | 4.65 | 0.766 | 0.673 | 0.88 |
| Ni | 5.04 | 0.877 | 0.321 | 0.37 |
| Au | 5.10 | 0.885 | 0.295 | 0.33 |
| Pd | 5.12 | 0.888 | 0.287 | 0.32 |
| Pt | 5.65 | 0.929 | 0.169 | 0.18 |

**This contradicts the 2026-09-07 conclusion and it must be stated as a
contradiction, not a refinement.** The single-contact model ranked the
metals Pt > Pd > Au > Ni > Ti > Cu > Cr, i.e. the larger the
work-function mismatch the better the collection (Pt best at 1.47x
enhancement). The two-contact model ranks a *symmetric* device almost
exactly the other way: Cu ~ Ti > Cr > Ni > Au > Pd > Pt, with Pt now
**worst**, and overstated by a factor of 5.5 by the old model.

The mechanism is straightforward once the sign is right. A strongly
doping contact produces a strong field that sweeps carriers toward
itself. Two of them, facing each other, produce a large opposing-field
region and a stagnation point near mid-channel that the modest bias
field (0.5 MV/m, against Pt's ~4.5 MV/m at the contact) cannot
overcome. Strong contact doping *helps* an isolated junction and *hurts*
a symmetric two-terminal device. Cr is the one metal unaffected (ratio
exactly 1.00) because its work function matches graphene's to within
0.01 eV, so it has essentially no doping field to cancel, and the bias
field sweeps the channel unopposed.

This matters for the thesis's device-design argument, which until now
pointed the same way in Chapter 4 (high-work-function metals for low
contact resistance) and Chapter 6 (high-work-function metals for
collection). Those are now in tension: the metals that are best for
contact resistance are the worst for symmetric-device photoresponse.
That tension is real device physics, not a modelling artefact, and it
is exactly the kind of trade-off a device-oriented thesis should be
making explicit.

### 3.2 Zero-bias net response for asymmetric pairs

At zero bias the symmetric diagonal is exactly zero, so every nonzero
entry is purely the work-function asymmetry -- the Weiss & Duan
mechanism, isolated:

| A \ B | Ti | Cr | Cu | Ni | Au | Pd | Pt |
|---|---|---|---|---|---|---|---|
| Ti | 0.000 | 0.596 | 0.067 | -0.753 | -0.794 | -0.802 | -0.903 |
| Cr | -0.596 | 0.000 | -0.563 | -0.835 | -0.849 | -0.853 | -0.917 |
| Cu | -0.067 | 0.563 | 0.000 | -0.778 | -0.805 | -0.812 | -0.905 |
| Ni | 0.753 | 0.835 | 0.778 | 0.000 | -0.078 | -0.103 | -0.581 |
| Au | 0.794 | 0.849 | 0.805 | 0.078 | 0.000 | -0.025 | -0.504 |
| Pd | 0.802 | 0.853 | 0.812 | 0.103 | 0.025 | -0.000 | -0.480 |
| Pt | 0.903 | 0.917 | 0.905 | 0.581 | 0.504 | 0.480 | 0.000 |

The matrix is antisymmetric, as it must be (swapping the contacts
reverses the current direction), which is a third internal consistency
check.

One non-obvious feature worth flagging, because it is the kind of detail
that distinguishes a mechanism model from a monotonic fit: the largest
|N| is **Cr/Pt at 0.917**, not Ti/Pt at 0.903 -- even though Ti/Pt has
the larger work-function difference (1.32 eV vs 1.15 eV). The reason is
that "asymmetry" here is not simply |W_A - W_B|: Ti has its own doping
field (|W_Ti - W_gr| = 0.17 eV) that sweeps carriers back toward Ti,
partially opposing Pt's, whereas Cr (0.00 eV) contributes no opposing
field at all. The best zero-bias pairing is therefore one strongly
doping contact against a *work-function-matched* one, not against an
oppositely-doping one -- at least under this model's magnitude
convention. Which is the immediate caveat: see Section 4.

## 4. The one simplification that could change the conclusion

Both doping fields are taken as sweeping carriers toward their own
contact, via |W_metal - W_gr|. That is the same magnitude convention the
2026-09-07 single-contact module already uses, and it is what makes
identical contacts cancel exactly. But it erases the n-type/p-type
distinction: Ti and Cu (W < W_gr) dope graphene n-type, while Pt, Pd and
Au (W > W_gr) dope it p-type, and a fully signed treatment would track
electrons and holes separately.

Under a signed treatment, an n/p pair such as Ti/Pt would ADD rather
than partially cancel for one carrier species, which would likely make
Ti/Pt the strongest pairing rather than the second-strongest, inverting
Section 3.2's Cr/Pt-vs-Ti/Pt conclusion. This is a genuine physical
extension, not a numerical refinement, and Section 3.2's ordering should
be treated as provisional until it is done. Listed as the top open item.

Other stated simplifications, unchanged from the single-contact module:
uniform illumination (the Shimomura shadow-mask device is precisely a
non-uniform-generation experiment, and a g(x) generation weight is the
natural way to model it); drift only, no diffusion; one effective
mobility; no photogain, photo-thermoelectric or bolometric contribution.

## 5. Files

- `graphene_photodetector_two_contact_model.py` -- the model, both
  validations, the result tables, the figure.
- `photodetector_two_contact_net_response.png` -- (a) the two opposing
  fields and the Pt/Pt stagnation point at x = 100 nm, (b) the
  single-contact vs symmetric two-contact bar comparison showing the
  reversal, (c) zero-bias net response vs W_A - W_B.
