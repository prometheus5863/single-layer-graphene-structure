# Research notes: per-metal literature contact-resistance values and Rc recalibration

*Date: 2026-08-31. Closes the follow-on item flagged repeatedly in
`AUTOMATION_LOG.md` since 2026-08-21 and restated in the Chapter 1 status
table as of 2026-08-29: "use [the] doping model to recalibrate per-metal
Rc." WebSearch and WebFetch were both available and used for this
session's literature search.*

## 1. Goal

`graphene_fet_model.py`'s `Rc_total` is a single generic contact
resistance (300 Ω·µm) drawn from a literature *range* (~110-500 Ω·µm),
not tied to a specific contact metal. `graphene_contact_doping_model.py`
(added 2026-08-26) computes a metal-work-function-dependent *extra*
resistance contribution, `R_extra`, from the extended in-plane doping
junction beyond the contact edge -- but it was never checked against
real per-metal measured Rc values. This session sources individually
cited, metal-specific measured Rc values from the literature and uses
them to recalibrate/sanity-check the doping-junction model.

## 2. Sourced literature values

All four values below are for **top (surface) contacts** (not
edge-contacts, which are a distinct, generally lower-Rc geometry --
see Section 5) at **room temperature**, extracted via the standard
transfer length method (TLM) unless noted otherwise.

| Metal | Rc (Ω·µm) | Condition | Source |
|---|---|---|---|
| Cu | 184 | Two-terminal device, post-350°C anneal, unpatterned contact, gate-biased (on-state) | Smith, Franklin, Farmer & Dimitrakopoulos, "Reducing Contact Resistance in Graphene Devices through Contact Area Patterning," *ACS Nano* 7(4), 3661-3667 (2013), https://franklin.pratt.duke.edu/files/u9/Papers/Smith_ACSnano_2013.pdf |
| Ni | 470 | "Two-in-one" fabrication process (metal deposited before photolithographic patterning, minimizing resist-residue contamination), gate-biased (on-state) | Khosravi Rad, Mehrfar, Sadeghi Neisiani, Khaje & Eslami Majd, "Effect of fabrication process on contact resistance and channel in graphene field effect transistors," *Scientific Reports* 14, 9190 (2024), https://www.nature.com/articles/s41598-024-58360-9 |
| Au | 519 | Unpatterned surface contact, back-gate biased away from the Dirac point (on-state) | Passi, Gahoi, Marin, Cusati, Fortunelli, Iannaccone, Fiori & Lemme, "Ultra Low Specific Contact Resistivity in Metal-Graphene Junctions via Atomic Orbital Engineering," arXiv:1807.04772, https://arxiv.org/pdf/1807.04772 |
| Pd | 584 | Top-gated FET, unpatterned ("normal") contact, gate-biased (on-state) | Smith et al. 2013 (as above) |

Two additional data points, not used quantitatively (different
condition/geometry than the four above, kept for context):

- Smith et al. (2013) also report **patterned** (hole-etched) contacts
  reduce Rc further: Cu 125 Ω·µm (32% reduction vs. 184), Pd 457 Ω·µm
  (22% reduction vs. 584) -- consistent with edge-contact geometries
  generally outperforming pure top contacts (Section 5 below).
- Passi et al. (arXiv:1807.04772) also report Au **at the Dirac point**
  (1372 Ω·µm, unpatterned) and Au **edge contacts** (456 Ω·µm at Dirac
  point, 45 Ω·µm gate-biased via 200 nm etched holes) -- the edge-contact
  gate-biased number (45 Ω·µm) is dramatically lower than the top-contact
  number (519 Ω·µm) for the *same metal*, direct experimental evidence
  that contact geometry, not just metal choice, is a first-order lever
  on Rc (a point noted qualitatively in `notes/2026-08-21-*.md` but not
  previously backed by a same-metal, same-paper, geometry-only
  comparison).
- The earlier-cited Pd figure of ~110 ± 20 Ω·µm (2026-08-21 notes,
  PubMed 21297624) is at **6 K**, not room temperature, and is reported
  to *increase* with temperature by an anomalous mechanism -- not
  directly comparable to the room-temperature Smith et al. 584 Ω·µm
  figure above; both are retained in the thesis as genuinely different
  measurement regimes rather than reconciled into one number.

## 3. Search access notes

- `pmc.ncbi.nlm.nih.gov` (the "Electrical properties of graphene-metal
  contacts" PMC review, PMC5506027) returned a reCAPTCHA interstitial
  page, not article content -- blocked, not used.
- `researchgate.net` returned HTTP 429 (rate-limited) on two follow-up
  fetch attempts: a Ti-specific process-variability paper (Pallecchi et
  al./RSC, ResearchGate ID 321487931, also available at
  arXiv:1712.00331 but that abstract page did not surface a specific
  Ω·µm figure either) and a gold-graphene temperature-dependence paper
  (Gahoi et al., ResearchGate ID 320662706). Both remain unsourced this
  session.
- As a direct result, **Ti and Cr are not included in the quantitative
  Rc recalibration below** (Section 4) -- no individually-cited,
  single-value Ω·µm figure was obtained for either this session. Both
  remain in `METAL_WORK_FUNCTIONS` and the (Rc-measurement-independent)
  doping-profile analysis from 2026-08-26.
- **Ni's work function** (5.04 eV, added to `METAL_WORK_FUNCTIONS` this
  session) uses the commonly cited Michaelson (1977) polycrystalline
  tabulation; reported values in the device literature range roughly
  4.9-5.35 eV depending on facet and surface preparation, the same
  spread caveat already applied to the other entries in
  `METAL_WORK_FUNCTIONS`.

## 4. Recalibration result (see `graphene_contact_doping_model.py`,
`recalibrate_metal_rc()` / `print_rc_recalibration()` /
`plot_rc_recalibration()`, and `rc_recalibration.png`)

The original plan was a clean additive decomposition,
`Rc_measured = R_extra(doping junction, this model) + R_transmission
(interface tunneling + mode-limited injection, not modeled here)`,
which would let the doping-junction model "explain" part of each
metal's measured Rc and leave a residual attributable to the interface
transmission mechanism already discussed qualitatively in
`notes/2026-08-21-*.md`, Section 2.1.

**That decomposition does not hold up as posed.** Computing `R_extra`
at the same on-state bulk channel density (`n_bulk = 2.0e16 m^-2`) used
throughout `graphene_contact_doping_model.py`:

| Metal | Rc, measured | R_extra, computed | Naive R_transmission |
|---|---|---|---|
| Cu | 184.0 | 971.5 | **-787.5** |
| Ni | 470.0 | 830.5 | **-360.5** |
| Au | 519.0 | 609.2 | **-90.2** |
| Pd | 584.0 | 533.1 | +51.0 (~9% of total) |

For three of four metals, the computed doping-junction contribution
*alone exceeds the entire measured literature Rc*, forcing a negative
"implied transmission resistance" -- unphysical, since a resistance
cannot be negative. Only Pd gives a plausible small positive residual.
This is reported here as a genuine finding, not smoothed over: **the
doping-junction model, taken at face value with its existing
`lambda_decay = 250 nm` parameter and the on-state `n_bulk` used
elsewhere in this thesis, is not simply additive with TLM-measured
lumped Rc.**

Two non-exclusive candidate explanations, neither resolved this
session:

1. **TLM double-counting.** All four literature values were extracted
   via the transfer length method, which fits total two-terminal
   resistance vs. contact-to-contact spacing and extrapolates back to
   zero spacing to isolate a "contact resistance." If the near-contact
   doping gradient this model computes is spatially compact enough
   relative to the TLM spacing range used in a given experiment, some
   or all of it can already be absorbed into the TLM-fitted contact
   term rather than showing up as a spacing-dependent (and therefore
   TLM-separable) channel effect. In that case `R_extra` and
   `Rc_measured` are not independent, additive series terms at all --
   they at least partially double-count the same physical region. None
   of the four sourced papers reported their extracted transfer length
   `L_T` in the excerpts obtained this session, so this could not be
   checked quantitatively.
2. **`lambda_decay` overshoot.** The 250 nm decay length (Khomyakov et
   al. 2010's doped-graphene asymptotic regime, already used
   unmodified since 2026-08-26) may simply be too long for these
   specific literature devices' real channel/contact geometry (e.g. if
   their channel length is itself shorter than a few `lambda_decay`,
   the "far-channel bulk" reference density this model integrates
   against is never actually reached in the real device).

Both are flagged as open items for a future session (see
`AUTOMATION_LOG.md`) rather than resolved by adjusting `lambda_decay`
or the decomposition ad hoc to force agreement -- doing so without a
principled reason would just be curve-fitting the parameter to the
four data points collected today.

## 5. Secondary finding: edge contacts vs. top contacts (same metal)

Independent of the decomposition issue above, the Au top-contact vs.
edge-contact comparison within the same paper (Passi et al.,
arXiv:1807.04772; Section 2) is a clean, single-source, single-metal,
geometry-only comparison: **519 Ω·µm (top, gate-biased) vs. 45 Ω·µm
(edge, 200 nm holes, gate-biased)**, an ~11x reduction from geometry
alone. This is a much larger lever than the metal-to-metal spread
observed at fixed (top-contact) geometry in Section 2's table
(184-584 Ω·µm, ~3x). Not modeled quantitatively in
`graphene_contact_doping_model.py` (which has no notion of contact
geometry, only work function and bulk doping), but worth flagging in
`thesis_draft/04-graphene-fet-device-physics.md` as context for why
real fabrication processes increasingly favor edge or quasi-edge
(hole-patterned) contacts over simple top contacts.

## References

1. Smith, J. T., Franklin, A. D., Farmer, D. B., Dimitrakopoulos, C. D.
   "Reducing Contact Resistance in Graphene Devices through Contact Area
   Patterning." *ACS Nano* 7(4), 3661-3667 (2013).
   https://franklin.pratt.duke.edu/files/u9/Papers/Smith_ACSnano_2013.pdf
2. Khosravi Rad, B., Mehrfar, A. H., Sadeghi Neisiani, Z., Khaje, M.,
   Eslami Majd, A. "Effect of fabrication process on contact resistance
   and channel in graphene field effect transistors." *Scientific
   Reports* 14, 9190 (2024). https://www.nature.com/articles/s41598-024-58360-9
3. Passi, V., Gahoi, A., Marin, E. G., Cusati, T., Fortunelli, A.,
   Iannaccone, G., Fiori, G., Lemme, M. C. "Ultra Low Specific Contact
   Resistivity in Metal-Graphene Junctions via Atomic Orbital
   Engineering." arXiv:1807.04772. https://arxiv.org/pdf/1807.04772
4. Pallecchi, E. et al. "Titanium Contacts to Graphene: Process-Induced
   Variability in Electronic and Thermal Transport." arXiv:1712.00331
   (abstract page reviewed; full text not obtained this session --
   ResearchGate mirror rate-limited). https://arxiv.org/abs/1712.00331
5. `notes/2026-08-21-contact-resistance-and-quantum-capacitance.md`,
   `notes/2026-08-26-contact-induced-doping-profile.md` (this repo) --
   prior lumped-Rc literature range and doping-junction model
   derivation reused here.
