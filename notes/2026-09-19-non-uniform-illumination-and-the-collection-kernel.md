# Non-uniform illumination, the collection kernel, and what a shadow mask can actually buy

**Date:** 2026-09-19
**Repo:** single-layer-graphene-structure
**Status of the open item:** this closes the item that has sat at the top of
"not yet covered" in `AUTOMATION_LOG.md` since 2026-09-18 — *"Non-uniform
illumination, a generation weight g(x)"*. Every two-contact model in this repo
(2026-09-07, 2026-09-17, 2026-09-18) states uniform illumination as an explicit
simplification, and each one names the shadow-mask experiment as the reason it
matters.

---

## 0. A citation correction before anything else

This repo has cited the shadow-mask experiment in six places as

> Suzuki *et al.*, *Carbon Trends* **5**, 100115 (2021)

(`graphene_photodetector_two_contact_model.py` lines 34 and 100 and its printed
output, `thesis_draft/06-graphene-photodetectors.md` §6.6/§6.9/References,
`notes/2026-09-17-two-contact-self-consistent-collection.md`).

**Both the first author and the article number are wrong.** The paper is

> Kenta Shimomura, Kaname Imai, Kenta Nakagawa, Akira Kawai, Kazuki Hashimoto,
> Takuro Ideguchi and Hideyuki Maki, **"Graphene photodetectors with asymmetric
> device structures on silicon chips"**, *Carbon Trends* **5**, 100100 (2021).
> <https://www.sciencedirect.com/science/article/pii/S2667056921000778>

Confirmed independently against the publisher page and against the Ideguchi
group's own publication list
(<https://www.ideguchi.ipst.s.u-tokyo.ac.jp/en/publications/>), which give the
same seven-author list, volume 5 and article number 100100. There is no Suzuki
in the author list, and 100115 is a different article. The error appears to have
been introduced on 2026-09-17 and copied forward twice without being rechecked.

This is worth more than the correction itself: the substantive physics
attributed to that citation *is* in Shimomura *et al.*, so nothing in
Sections 6.6–6.8 changes. But it is a reminder that a citation which is never
re-fetched propagates. All six occurrences are corrected in this session's
commits.

---

## 1. What the experiment actually does, and the number nobody in this repo used

Shimomura *et al.* build metal–graphene–metal detectors on silicon and make them
asymmetric in **geometry** rather than in metal. Two designs:

1. **A shadow mask** — 50 nm of nickel deposited over one of the two
   graphene/electrode interfaces, so that one interface is shaded and the other
   is not.
2. **A comb-shaped counter-electrode**, giving one interface much more contact
   perimeter than the other.

Their statement of the cancellation this repo has been modelling since
2026-09-07 is explicit:

> "since the polarities of the photovoltages at each graphene/electrode
> interface of two electrodes on both ends of the graphene are opposite, they
> are canceled out under macroscopic light irradiation."

That is the same zero that `validate_symmetric_cancellation()` checks, stated
experimentally. Design 1 breaks it by removing photons from one side rather than
by changing either contact.

**The number this repo has never used:** the mask is *leaky*. Shimomura *et al.*
report that

> "the photovoltage under the shadow mask was about half of the opposite side"

at 690 nm. So the physically realised mask has a transmission of roughly
**T ≈ 0.5**, not 0. That converts directly into a prediction below, and it turns
out to cost far more than the factor of two one would guess.

Device parameters for context: ~15-layer exfoliated graphene, two-terminal
resistances 1020 Ω (mask device) and 707.8 Ω (comb device), measured at 690 nm,
1310/1530 nm and 4.6 µm, with a ~2 µm focused spot for the position-resolved
measurements.

---

## 2. Scanning photocurrent microscopy is the same measurement, done with light
   instead of metal

The other way to illuminate non-uniformly is to focus the beam and move it —
scanning photocurrent microscopy (SPCM), "a method that allows the spatial
mapping of the photoresponse by raster scanning a focused laser beam over the
sample" (T. S. Kasırga, *Scanning photocurrent microscopy and its application to
one- and two-dimensional materials*, arXiv:2509.09390, 2025,
<https://arxiv.org/abs/2509.09390>).

For graphene transistors specifically, the reference already in this repo is

> T. Mueller, F. Xia, M. Freitag, J. Tsang and Ph. Avouris, "The role of contacts
> in graphene transistors: A scanning photocurrent study", *Phys. Rev. B* **79**,
> 245430 (2009). arXiv:0902.1479, <https://arxiv.org/abs/0902.1479>

whose abstract gives the number this repo already uses — the contact-induced
doping "extends 0.2–0.3 µm into the graphene channel" and a p–n–p structure
forms along the device. A p–n–p structure is precisely a device whose local
photoresponse **changes sign** between the two ends, which is the structure
predicted below.

**An honest caveat on the review, recorded rather than worked around.**
The Kasırga review's main thrust is that SPCM signals are frequently
*photothermal* rather than photovoltaic, and it "discuss[es] the shortcomings of
SPCM in determining the mechanisms leading to the photoresponse". This repo
models drift collection only (no photo-thermoelectric and no bolometric term;
see the standing simplification 4 in every photodetector module). So the
position-resolved curves computed here are the **photovoltaic contribution to**
an SPCM trace, not a prediction of a measured SPCM trace. Stating the difference
is the point; the repo has no PTE machinery with which to close it.

---

## 3. The structural fact that makes this cheap to model

Write the 2026-09-18 signed, carrier-resolved response as an integral over the
channel. Nothing about the illumination enters the transport at all:

    N = (1/L) ∫₀ᴸ [ q_h s_h(x) p_h(x) + q_e s_e(x) p_e(x) ] dx
      ≡ (1/L) ∫₀ᴸ k(x) dx

The bracket — call it the **collection kernel** k(x) — depends only on the two
contact metals, the bias, the crossover and τ. It is the net charge delivered to
contact A **per photon absorbed at x**. Non-uniform illumination therefore does
not require a new transport model at all; it requires only a generation weight:

    N[g] = (1/L) ∫₀ᴸ g(x) k(x) dx ,      with (1/L) ∫₀ᴸ g(x) dx = 1

normalised so that every illumination pattern is compared at **equal total
absorbed photon number**, which is the only comparison that means anything for
a detector.

Three consequences follow immediately, and all three are checkable exactly:

1. **g ≡ 1 must return the 2026-09-18 number bitwise.** Not approximately —
   the expression is literally the same one.
2. **k(x) is the delta-spot SPCM trace.** N[δ(x − x₀)] = k(x₀). The kernel is
   not an intermediate quantity; it is the predicted position scan.
3. **max|k| is a hard ceiling on illumination engineering.** For any normalised
   g, |N[g]| ≤ max_x |k(x)|, with equality only for illumination concentrated
   where |k| peaks. No mask, spot, grating or plasmonic pattern can exceed it
   at fixed photon number. This is the first quantity in this repo that bounds
   an entire design space rather than ranking points inside it.

And for a **symmetric** device at zero bias, k is antisymmetric about the
midpoint, which yields a closed form for the leaky mask. With mask transmission
T on the left half and 1 on the right, and A ≡ (1/L)∫_{L/2}^{L} k dx,

    N(T) = 2A (1 − T) / (1 + T)     ⇒     N(T)/N(0) = (1 − T)/(1 + T)

so the experimentally reported **T ≈ 0.5 leaves exactly one third** of what a
perfect mask would give, not one half. That is an exactly known value the code
can be tested against, and it is the strongest argument in this note for caring
about mask quality: halving the leakage from 0.5 to 0.25 buys a factor 1.8,
while going from a perfect mask to a merely good one (T = 0.1) already costs 18%.

---

## 4. What to expect, and what would falsify it

- A **symmetric** pair, which gives identically zero under uniform illumination
  at zero bias, must give a nonzero response under any mirror-asymmetric g. This
  is the only mechanism left by which a symmetric device can respond at all,
  and it matters more since 2026-09-18, when the symmetric *metal ranking* was
  retracted.
- An **n/p pair** (Ti/Pt) should gain little or nothing from masking: its kernel
  does not change sign, so masking half the channel only discards photons. If
  the model says masking helps Ti/Pt, the model is wrong.
- The kernel for a symmetric pair should show the experimentally reported
  **opposite polarities at the two interfaces** (Shimomura *et al.*'s "polarities
  ... are opposite"; Mueller *et al.*'s p–n–p).

**Geometry caveat, stated up front.** This repo's channel is L = 200 nm
(`graphene_photodetector_model.py`, inherited from `graphene_fet_model.py`).
Shimomura *et al.*'s focused spot is ~2 µm — **ten times the entire channel**.
Optical localisation of the kind SPCM performs is therefore *impossible* at this
repo's device geometry; a diffraction-limited spot illuminates the whole channel
and both contacts at once. The only realisable non-uniform illumination here is
a **lithographic** mask sitting on the device, which is exactly what Shimomura
*et al.* built. The spot-size study in the accompanying module is therefore run
as a statement about what channel length would be needed, not as a proposal for
this one.

---

## 5. Sources

- Kenta Shimomura, Kaname Imai, Kenta Nakagawa, Akira Kawai, Kazuki Hashimoto,
  Takuro Ideguchi, Hideyuki Maki, "Graphene photodetectors with asymmetric
  device structures on silicon chips", *Carbon Trends* **5**, 100100 (2021).
  <https://www.sciencedirect.com/science/article/pii/S2667056921000778>
  — shadow mask (50 nm Ni over one interface), comb electrode, the
  "canceled out under macroscopic light irradiation" statement, and the
  "about half of the opposite side" mask leakage.
- Ideguchi group publication list, used to confirm the author order and article
  number. <https://www.ideguchi.ipst.s.u-tokyo.ac.jp/en/publications/>
- T. S. Kasırga, "Scanning photocurrent microscopy and its application to one-
  and two-dimensional materials", arXiv:2509.09390 (2025).
  <https://arxiv.org/abs/2509.09390> — SPCM definition; the photothermal caveat.
- T. Mueller, F. Xia, M. Freitag, J. Tsang, Ph. Avouris, "The role of contacts in
  graphene transistors: A scanning photocurrent study", *Phys. Rev. B* **79**,
  245430 (2009), arXiv:0902.1479. <https://arxiv.org/abs/0902.1479>
  — contact doping extending 0.2–0.3 µm, p–n–p structure.
- T. Mueller, F. Xia, Ph. Avouris, "Graphene photodetectors for high-speed
  optical communications", *Nature Photonics* **4**, 297 (2010)
  — the Pd/Au vs Ti/Au asymmetric-metal device; already cited in this repo.

**Fetch record (honest):** 4 WebSearch queries and 4 WebFetch calls, all through
allowed domains. `arxiv.org/pdf/2509.09390` returned no machine-readable text
(the abstract page did), so the review is cited from its abstract only and no
number is attributed to its body. PubMed/PMC was **not** attempted: four
consecutive reCAPTCHA blocks (2026-09-05, 09-07, 09-17, noted 09-18) make it a
standing environment limitation, so the Park/Ahn SPCM sign map remains
unreachable — but the sign-map physics is covered here by Mueller *et al.* 2009
instead, which is open-access on arXiv.
