# Automation Log

This file tracks daily automated progress on the graphene thesis repo, so
each run can pick up where the last one left off without repeating work.

---

## 2026-08-21

**Status:** First automation run (no prior log existed).

**Repo state at start:** README + `graphene_band_structure.py`,
`graphene_transport_properties.py`, `optical_absorption`-related plots and
scripts covering intrinsic material properties only (band structure,
DOS, quantum Hall effect, universal optical absorption). No device-level
content (contacts, gate capacitance, transistor I-V, RF metrics) existed yet.

**Work done:**

1. **Research notes** (`notes/2026-08-21-contact-resistance-and-quantum-capacitance.md`):
   Literature review of metal-graphene contact resistance (physical origin,
   representative values from recent papers: Pd ~110 Ohm.um at 6K, Ni ~470
   Ohm.um from a 2024 fabrication-process study) and graphene quantum
   capacitance (functional form, role in vertical-scaling limits, and why
   it requires self-consistent I-V modeling). Full citations included.
   Web search (WebSearch tool) was available and used for this research.

2. **Code + results** (`graphene_fet_model.py`,
   `quantum_capacitance.png`, `gfet_transfer_characteristics.png`):
   New device-physics module implementing a thermally-broadened quantum
   capacitance model, a self-consistent channel carrier density combining
   C_ox and C_q, and a long-channel Id-Vg transfer characteristic model for
   a 200 nm GFET including a literature-calibrated series contact
   resistance term. Verified to run standalone and produce physically
   sensible output (V-shaped ambipolar transfer curve, correct suppression
   of on-current with increasing contact resistance).

**Not yet covered (candidates for future runs):**
- RF small-signal model / f_T, f_max figures of merit
- Graphene interconnect resistivity vs. linewidth (finite-size/edge
  scattering effects)
- Photodetector responsivity
- Metal contact work-function-dependent doping profile (currently only a
  lumped series resistance, not a spatially resolved p-n junction model)
- thesis_draft/ section expansion (no thesis_draft/ content written yet --
  next run could draft an introduction or a "device applications" literature
  review section using the notes gathered so far)
- Explicit contact-resistance-vs-channel-length crossover plot (planned in
  notes section 4, not yet implemented)

**Commits this run:** 2 (notes commit, code+results commit).

---

## 2026-08-22

**Status:** Second automation run.

**Repo state at start:** README + intrinsic-property scripts, plus
2026-08-21's device-physics addition (contact resistance / quantum
capacitance notes, `graphene_fet_model.py` with quantum-capacitance-
limited GFET transfer characteristics). No RF figures of merit and no
`thesis_draft/` content existed yet.

**Work done:**

1. **Research notes**
   (`notes/2026-08-22-rf-figures-of-merit-fT-fmax.md`): Literature review
   of graphene FET RF cutoff frequency (f_T) and maximum oscillation
   frequency (f_max) — representative values (intrinsic f_T up to
   ~300-427 GHz in record devices, f_max ~200 GHz for a 60 nm T-gate
   device, extrinsic f_T/f_max ~34/37 GHz for foundry-style 0.5-2 um
   gates), the physical reason f_max lags f_T by roughly an order of
   magnitude for graphene (lack of current saturation), and how contact
   resistance (from the 2026-08-21 notes) and quantum capacitance
   directly set the small-signal model parameters. Full citations
   included. WebSearch was available and used for this research.

2. **Code + results** (`rf_small_signal_model.py`,
   `rf_figures_of_merit.png`): New module implementing a hybrid-pi
   small-signal RF model on top of the existing GFET DC model — g_m and
   g_ds via numerical differentiation of `transfer_characteristic()`,
   C_gs reusing the existing quantum-capacitance-limited gate
   capacitance, C_gd as a fraction of C_gs, and a distributed gate-
   resistance estimate. Computes and plots f_T(V_g) and a simplified
   f_max(V_g), plus a second panel showing extrinsic f_T degradation
   across the literature contact-resistance range (0-500 Ω·µm). Verified
   to run standalone; peak f_T ~20 GHz for the 200 nm long-channel device
   is consistent with the literature range for non-exotic gate lengths.
   Added an explicit code-level caveat (printed at runtime) for the
   bias points where the simplified f_max estimate comes out above f_T,
   since that contradicts the literature trend — traced to this compact
   model's idealized (low) parasitic gate resistance rather than hidden.

3. **Writing** (`thesis_draft/01-introduction.md`): First
   `thesis_draft/` content (folder did not exist before today). Drafted
   Chapter 1 introduction motivating the material-physics vs. device-
   applications split, summarizing the three device-physics effects
   covered so far (contact resistance, quantum capacitance, absence of
   current saturation and its RF consequences), and a chapter-by-chapter
   status table kept in sync with the codebase.

**Not yet covered (candidates for future runs):**
- Contact-resistance-vs-channel-length crossover plot (still planned
  from 2026-08-21 notes, not yet implemented)
- Graphene interconnect resistivity vs. linewidth (finite-size/edge
  scattering effects) — Chapter 5 material
- Photodetector responsivity — Chapter 6 material
- Metal contact work-function-dependent doping profile (still only a
  lumped series resistance, not a spatially resolved p-n junction model)
- A more careful f_max model with realistic pad/interconnect parasitics,
  to resolve the f_max > f_T artifact flagged in today's code caveat
- Further thesis_draft/ chapters (2-3 could be drafted from the existing
  band-structure/transport code and results; 4 partially supported by
  today's RF work)

**Commits this run:** 3 (RF research notes, RF small-signal model code +
plot, thesis_draft/ Chapter 1 introduction).

---

## 2026-08-23

**Status:** Third automation run.

**Repo state at start:** README + intrinsic-property scripts, plus
2026-08-21/22 device-physics additions (contact resistance / quantum
capacitance notes and model, RF small-signal model, thesis_draft/
Chapter 1 introduction). No interconnect-resistivity content and no
Chapter 5 draft existed yet.

**Work done:**

1. **Research notes**
   (`notes/2026-08-23-interconnect-resistivity-vs-linewidth.md`):
   Literature review of graphene nanoribbon (GNR) interconnect
   resistivity vs. linewidth — physical origin of edge/line-edge-
   roughness scattering, the specularity-parameter edge-scattering model
   (adapted from Fuchs-Sondheimer thin-film theory), representative
   experimental data (Murali et al., arXiv:0906.0924: 15-25 uOhm.cm for
   18-52 nm GNRs, best individual GNR ~3x the phonon-scattering-limited
   intrinsic value, theoretical Naeemi & Meindl crossover projections),
   and a 2024 graphene-all-around-cobalt interconnect result (27%
   resistance reduction). Full citations included. WebSearch was
   available and used for this research.

2. **Code + results** (`graphene_interconnect_model.py`,
   `interconnect_resistivity_vs_linewidth.png`): New module implementing
   a Matthiessen's-rule combination of bulk, specularity-parameter edge-
   scattering, and residual-impurity mean free path terms to model
   graphene resistivity vs. linewidth, calibrated against the *best-case*
   Murali et al. result (kept deliberately distinct from the *average*
   measured cluster, which is overlaid as an empirical reference band
   rather than used as the fit target -- see the module docstring and
   notes for the reasoning). Includes a simplified Fuchs-Sondheimer-style
   Cu comparison model. Verified to run standalone and produce physically
   sensible, non-degenerate behavior across the specularity sweep (an
   earlier draft of the calibration had a bug where the impurity term was
   re-solved per specularity value, silently absorbing all p-dependence
   and making the curves nearly indistinguishable -- caught by inspecting
   the numerical output, not just the plot, and fixed before committing).
   Results: for realistic diffuse edges (p=0.15) the model crosses below
   the Cu projection at W ~ 130 nm; for idealized near-specular edges
   (p=0.9) graphene stays below Cu across the full 2-500 nm scan range.

3. **Writing** (`thesis_draft/05-graphene-interconnects.md`): First
   Chapter 5 content (previously "Not started" in the Chapter 1 status
   table, now updated to "In progress"). Drafted the interconnect-
   resistivity motivation, the specularity-based physical model, the
   realistic-vs-idealized crossover results from today's code, and an
   explicit methodological link back to Chapter 4's contact-resistance
   specularity/mode-counting discussion (both trace to the same
   underlying "atomically thin conductor is sensitive to boundary
   quality" physical picture) as a note for the eventual discussion
   chapter (Chapter 7).

**Not yet covered (candidates for future runs):**
- Contact-resistance-vs-channel-length crossover plot (still planned
  from 2026-08-21 notes, not yet implemented)
- Photodetector responsivity — Chapter 6 material
- Metal contact work-function-dependent doping profile (still only a
  lumped series resistance, not a spatially resolved p-n junction model)
- Graphene-all-around-metal (liner/cap) interconnect model, distinct from
  today's pure-GNR-wire model (flagged as a specific follow-on in
  today's Chapter 5 draft, Section 5.4)
- Copper comparison model does not yet include the liner/barrier-
  thickness effect (noted as a conservative simplification in today's
  work, Chapter 5 Section 5.2)
- A more careful f_max model with realistic pad/interconnect parasitics
  (still open from 2026-08-22, to resolve the f_max > f_T artifact
  flagged in that day's code caveat)
- Further thesis_draft/ chapters: 2-3 could be drafted from the existing
  band-structure/transport code and results; 6-7 await their respective
  research/code contributions

**Commits this run:** 3 (interconnect resistivity research notes,
interconnect resistivity model code + plot, thesis_draft/ Chapter 5 +
Chapter 1 status table update).

---

## 2026-08-24

**Status:** Fourth automation run.

**Repo state at start:** README + intrinsic-property scripts, plus
2026-08-21/22/23 device-physics additions (contact resistance / quantum
capacitance model, RF small-signal model, interconnect resistivity model,
thesis_draft/ Chapters 1 and 5). The contact-resistance-vs-channel-length
crossover plot had been flagged as planned since 2026-08-21 but not yet
implemented; no photodetector content (research, code, or thesis_draft/
Chapter 6) existed yet.

**Work done:**

1. **Code + results** (`contact_resistance_crossover.py`,
   `contact_resistance_crossover.png`): New module closing out the
   contact-resistance-vs-channel-length crossover analysis planned in the
   2026-08-21 notes (Section 4). Reuses the existing quantum-capacitance-
   limited carrier density and sheet conductivity functions from
   `graphene_fet_model.py` unmodified, sweeping channel length L instead
   of gate voltage. Computes total device resistance (channel + literature
   contact resistance) vs. L on a log scale, plus the contact-resistance
   fraction of total resistance vs. L, and reports the analytic crossover
   length L_x (closed form since sheet conductivity is L-independent) for
   each literature Rc value. Verified to run standalone and produce
   physically sensible output: L_x = 115 nm for the best-case literature
   contact (Pd, ~110 Ohm.um), 314 nm for a mid-range contact (300 Ohm.um),
   and 524 nm for the worst-case (500 Ohm.um) -- confirming that even the
   best literature graphene contacts are not negligible at the 200 nm
   channel length used elsewhere in this thesis's GFET model.

2. **Research notes**
   (`notes/2026-08-24-graphene-photodetector-responsivity.md`): Literature
   review of graphene photodetector responsivity -- the intrinsic
   absorption bottleneck (ties directly to the ~2.3% universal absorption
   already derived in `graphene_transport_properties.py`), the collection
   bottleneck (short photocarrier lifetime, ~100-200 nm built-in-field
   region), the photogating gain/response-speed tradeoff (representative
   values from ~10^3 A/W at ~400 ns to ~10^10 A/W at second-scale response
   times), and heterostructure strategies (graphene/Si Schottky ~510
   mA/W, plasmonic-enhanced graphene/Si ~1.9 A/W, graphene/TMD photogating
   up to ~4.4x10^6 A/W, a 2025 alternating-channel geometry reporting
   1.7x10^7 mA/W with 3-4 us response, and a heterostructure-engineered
   160 Gb/s zero-bias design). Full citations included. WebSearch was
   available and used for this research.

3. **Writing** (`thesis_draft/06-graphene-photodetectors.md`): First
   Chapter 6 content (previously "Not started" in the Chapter 1 status
   table, now updated to "In progress"). Drafted the absorption-limit
   motivation (explicitly connecting back to the Chapters 2-3 optical-
   absorption result), the two-bottleneck picture, the gain/speed
   tradeoff, and a specific follow-on plan (a quantitative responsivity/
   gain model reusing the existing optical-conductivity code). Also
   updated the Chapter 1 status table (Chapters 4 and 6 rows) to reflect
   today's work.

**Not yet covered (candidates for future runs):**
- Quantitative responsivity/gain model for graphene photodetectors
  (Chapter 6, Section 6.4 -- reusing `calculate_optical_conductivity`
  with a parametrized gain factor)
- Metal contact work-function-dependent doping profile (still only a
  lumped series resistance, not a spatially resolved p-n junction model)
- Graphene-all-around-metal (liner/cap) interconnect model, distinct from
  the existing pure-GNR-wire model (Chapter 5, Section 5.4)
- Copper comparison model does not yet include the liner/barrier-
  thickness effect (Chapter 5, Section 5.2)
- A more careful f_max model with realistic pad/interconnect parasitics
  (still open from 2026-08-22, to resolve the f_max > f_T artifact
  flagged in that day's code caveat)
- Further thesis_draft/ chapters: 2-3 could be drafted from the existing
  band-structure/transport code and results; Chapter 7 (discussion/
  outlook) awaits enough device-application chapters to synthesize

**Commits this run:** 3 (contact-resistance-vs-channel-length crossover
code + plot, photodetector responsivity research notes, thesis_draft/
Chapter 6 + Chapter 1 status table update).

---

## 2026-08-25

**Status:** Fifth automation run.

**Repo state at start:** README + intrinsic-property scripts, plus
2026-08-21 through 2026-08-24 device-physics additions (contact
resistance / quantum capacitance model, RF small-signal model,
interconnect resistivity model, contact-resistance-vs-channel-length
crossover analysis, thesis_draft/ Chapters 1, 5, and 6). Chapter 6's
"Planned follow-on work" (Section 6.4, as of 2026-08-24) explicitly
called for a quantitative responsivity/gain model; this was the clearest
next candidate since it closes a specific, already-scoped gap rather
than opening a new topic.

**Work done:**

1. **Code + results** (`graphene_photodetector_model.py`,
   `photodetector_responsivity_gain_tradeoff.png`): New module
   implementing the responsivity-vs-response-time / photoconductive-gain
   tradeoff described qualitatively in Chapter 6, Section 6.3. Combines
   the literature bare-device EQE (~0.15%) with a classic photoconductor
   gain G = tau_trap/tau_transit, where tau_transit is derived from the
   existing Chapter 4 GFET channel parameters (L=200nm, mu=0.4 m^2/Vs,
   Vds=0.1V, giving tau_transit=1.0 ps) rather than a new free parameter.
   Verified to run standalone; checked against the three literature
   anchor points from the 2026-08-24 notes (interfacial photogating @
   400ns/1e3 A/W, 2025 alternating-channel @ 3.5us/1.7e4 A/W, extended
   photogating @ 1s/1e10 A/W): the model reproduces the correct tradeoff
   slope across six decades of tau_trap, underpredicting absolute
   responsivity by a roughly constant 3.8x-15.0x factor attributed to
   device-specific trap/geometry effects outside the model's single
   lumped transit time. Also derives a model gain-bandwidth invariant
   (~159 GHz) and shows all three literature devices exceed it, with the
   2025 geometry-engineered device behaving like the trap-based points
   rather than showing a qualitatively different signature -- consistent
   with it decoupling gain from speed rather than using a different gain
   mechanism.

2. **Writing** (`thesis_draft/06-graphene-photodetectors.md`,
   `thesis_draft/01-introduction.md`): Added Section 6.4 (model
   description, results table, gain-bandwidth-invariant discussion, and
   an explicit caveat on the model's known absolute-magnitude offset),
   renumbered the former follow-on-work list to Section 6.5 with
   narrower remaining items, fixed a broken cross-reference to the
   2026-08-24 research notes filename, and updated the Chapter 1 status
   table's Chapter 6 row.

**Not yet covered (candidates for future runs):**
- Plasmonic-absorption-enhancement factor for the photodetector model
  (Chapter 6, Section 6.5), to close part of the 3.8x-15.0x
  literature-offset identified today
- Metal contact work-function-dependent doping profile (still only a
  lumped series resistance, not a spatially resolved p-n junction model)
  -- also flagged today as a prerequisite for a spatially resolved
  photodetector collection-efficiency model
- Graphene-all-around-metal (liner/cap) interconnect model, distinct from
  the existing pure-GNR-wire model (Chapter 5, Section 5.4)
- Copper comparison model does not yet include the liner/barrier-
  thickness effect (Chapter 5, Section 5.2)
- A more careful f_max model with realistic pad/interconnect parasitics
  (still open from 2026-08-22, to resolve the f_max > f_T artifact
  flagged in that day's code caveat)
- Further thesis_draft/ chapters: 2-3 could be drafted from the existing
  band-structure/transport code and results; Chapter 7 (discussion/
  outlook) awaits enough device-application chapters to synthesize

**Note:** WebSearch was not used this run -- today's work extended an
existing literature-anchored model with new derived quantities and
cross-checks rather than requiring new literature search, so no live
search was needed.

**Commits this run:** 2 (photodetector responsivity/gain-bandwidth model
code + plot, thesis_draft/ Chapter 6 Section 6.4 + Chapter 1 status
table update).

---

## 2026-08-26

**Status:** Sixth automation run.

**Repo state at start:** README + intrinsic-property scripts, plus
2026-08-21 through 2026-08-25 device-physics additions (contact resistance
/ quantum capacitance model, RF small-signal model, interconnect
resistivity model, contact-resistance-vs-channel-length crossover
analysis, quantitative photodetector responsivity/gain-bandwidth model,
thesis_draft/ Chapters 1, 5, and 6). "Metal contact work-function-dependent
doping profile (currently only a lumped series resistance, not a
spatially resolved p-n junction model)" had been flagged as not yet
covered in every run's log since 2026-08-21; this was the clearest
remaining gap in the Chapter 4 device-physics material, and Chapter 4
itself had no `thesis_draft/` file yet (only referenced from the Chapter 1
status table).

**Work done:**

1. **Research notes**
   (`notes/2026-08-26-contact-induced-doping-profile.md`): Literature
   review of metal-work-function-dependent doping of graphene at contacts
   -- the charge-transfer mechanism and its dependence on the metal/
   graphene work-function difference (Giovannetti et al., *Phys. Rev.
   Lett.* 101, 026803 (2008)), and the spatial extent and screening
   behavior of the induced doping (Khomyakov et al., *Phys. Rev. B* 82,
   115437 (2010), arXiv:0911.2027: x^-1 potential decay for doped
   graphene, doping extending hundreds of nm from the contact edge, n/p
   crossover work function ≈5.4 eV vs. graphene's own ~4.5 eV). Full
   citations included. WebSearch was available and used for this
   research.

2. **Code + results** (`graphene_contact_doping_model.py`,
   `contact_doping_profile.png`): New module computing a self-consistent
   contact-edge carrier density (reusing `quantum_capacitance()` from
   `graphene_fet_model.py` with an effective interface capacitance in
   place of the back-gate oxide capacitance) for six representative
   contact metals (Ti, Cr, Cu, Pd, Au, Pt), a saturating spatial decay
   profile matching the literature x^-1 asymptotic, and the resulting
   extra junction sheet resistance (width-normalized, comparable to the
   lumped Rc literature range) integrated across the junction region.
   Verified to run standalone; caught and fixed a divide-by-zero at the
   charge-neutrality crossing by applying the same disorder-puddle
   regularization already used in `graphene_fet_model.carrier_density()`.
   Results show Ti (weakest doping) contributes the smallest extra
   resistance (~55 Ω·µm); Cu/Au/Pd (p-n junctions against an n-type bulk
   channel) contribute several hundred Ω·µm; Pt's very strong p-type
   doping produces a net *negative* extra-resistance figure in this
   model, which is flagged explicitly at runtime as a real feature of the
   classical drift-conductance picture (an overdoped access region
   conducts well) rather than a claim that strongly-doping contacts are
   easier to make low-resistance in reality -- the model does not capture
   depletion/injection physics at the crossing itself.

3. **Writing** (`thesis_draft/04-graphene-fet-device-physics.md`,
   `thesis_draft/01-introduction.md`): First Chapter 4 content (the
   chapter previously existed only as a status-table row referencing the
   codebase, with no drafted text). Synthesizes the existing quantum-
   capacitance, GFET transfer-characteristic, contact-resistance, and
   contact-resistance-vs-channel-length-crossover material into a single
   narrative, then adds a new Section 4.5 covering today's contact-doping
   model and Section 4.6 covering the RF figures of merit (including the
   still-open f_max artifact from 2026-08-22). Updated the Chapter 1
   status table's Chapter 4 row (now "draft written") and Chapter 6 row
   (noting the new contact-doping model as an available, not-yet-
   integrated building block for the photodetector chapter's flagged
   "spatially resolved collection model" prerequisite).

**Not yet covered (candidates for future runs):**
- Using the new contact-doping model to *recalibrate* per-metal Rc values
  (currently diagnostic only; Section 4.7 open-items table)
- Plasmonic-absorption-enhancement factor for the photodetector model
  (Chapter 6, Section 6.5)
- Integrating the contact-doping spatial-profile machinery into the
  photodetector collection-efficiency model (Chapter 6, Section 6.5's
  other flagged prerequisite -- now more directly reachable given
  today's Section 4.5 work)
- Graphene-all-around-metal (liner/cap) interconnect model, distinct from
  the existing pure-GNR-wire model (Chapter 5, Section 5.4)
- Copper comparison model does not yet include the liner/barrier-
  thickness effect (Chapter 5, Section 5.2)
- A more careful f_max model with realistic pad/interconnect parasitics
  (still open from 2026-08-22, to resolve the f_max > f_T artifact
  flagged in that day's code caveat, and now also referenced from
  Chapter 4 Section 4.6)
- Further thesis_draft/ chapters: 2-3 could be drafted from the existing
  band-structure/transport code and results; Chapter 7 (discussion/
  outlook) awaits enough device-application chapters to synthesize

**Commits this run:** 4 (contact-doping research notes, contact-doping
model code + plot, thesis_draft/ Chapter 4 + Chapter 1 status table
update, this log entry).

---

## 2026-08-26 (second session)

**Status:** Seventh automation run. Note on dating: the previous "## 2026-08-26"
entry above was actually written at 2026-08-25 23:18 UTC / 2026-08-26 04:48
IST (the automation host's local clock is IST, which is ahead of UTC and
crossed midnight first) — this session started at 2026-08-26 17:38 UTC,
genuinely the next scheduled run, ~14 hours later. Disambiguated with
"(second session)" in the header rather than reusing an identical date
heading.

**Repo state at start:** Phase-equivalent state per the previous entry:
Chapters 1, 4, 5, 6 of `thesis_draft/` exist; `rf_small_signal_model.py`
(added 2026-08-22) had a long-standing flagged caveat that its simplified
f_max estimate sometimes exceeds f_T, treated in the code and in
Chapter 4 as an unresolved artifact requiring a fix. This was the clearest
concrete open item carried across multiple prior "not yet covered" lists
(2026-08-24, 2026-08-25, 2026-08-26 first session all listed it), and,
unlike some of the other open items (plasmonic enhancement, liner/cap
interconnect model), it corresponds to an actual modeling gap rather than
a new feature — so it was prioritized this session.

**Work done:**

1. **Research notes**
   (`notes/2026-08-26-fmax-parasitics-and-fT-fmax-ratio.md`): WebSearch
   and WebFetch were both available and used this session. Re-examined
   the 2026-08-22 note's claim that f_max is "typically about an order of
   magnitude lower than f_T" for graphene FETs against more specific
   literature: Feijoo et al., *Sci. Rep.* 6, 35717 (2016) — already cited
   in the 2026-08-22 notes but not read closely enough at the time —
   reports f_max/f_T ratios of 1.3-1.4 (f_max *exceeding* f_T) at every
   gate length measured, for a device with a deliberately engineered low
   gate resistance (~15 Ω). Also pulled the standard source/drain access-
   resistance term (missing from this repo's f_max formula since it was
   first written) from Wang/Hsu/Wu, arXiv:1112.4831, along with that
   paper's note that GSG pad capacitance is a materially larger effect on
   conductive-substrate (i.e. back-gated, like this repo's device)
   devices than on insulating-substrate ones. Two fetch attempts failed
   (IEEE Xplore PDF: HTTP 418; a ResearchGate page: rate-limited) and are
   recorded as such rather than silently skipped.

2. **Code fix** (`rf_small_signal_model.py`, `rf_figures_of_merit.png`):
   Added the missing source/drain access resistance term
   (R_s = Rc_total/2, reusing the existing literature-calibrated contact
   resistance) to the f_max denominator: `gds*(Rg+Rs) + ...` instead of
   `gds*Rg + ...`. Added an explicit intrinsic-vs-extrinsic (GSG pad
   capacitance) comparison mode. Discovered while implementing this that
   a naive pad-capacitance comparison at the repo's normalized 1 µm
   device width overstates the pad effect by >100x relative to the
   literature raw/de-embedded ratio (~0.6-0.7x), because a real
   RF-probed device is never 1 µm wide — fixed by adding an optional
   width/finger-count override (`compute_fT_fmax(W=..., N_fingers=...)`)
   that temporarily rescales `graphene_fet_model`'s module-level `W` and
   `Rc_total` (following the same save/restore idiom the 2026-08-22 code
   already used for its Rc sweep), plus the standard N²-reduction
   multi-finger gate-resistance formula, and evaluates the extrinsic
   comparison at a literature-representative 40 µm/8-finger device
   instead. At that scale, extrinsic degrades f_T and f_max to ~15-20% of
   their intrinsic values — same direction, larger magnitude than
   Feijoo et al.'s measured ratio, and reported as such rather than
   tuned to match. Also replaced the old blanket "f_max > f_T is
   unphysical" caveat text in `summary_numbers()` with the corrected
   framing from the research notes. Verified to run standalone; the
   updated plot now has three panels (intrinsic f_T/f_max, the existing
   Rc-vs-fT sweep, and the new intrinsic-vs-extrinsic comparison).

3. **Writing** (`thesis_draft/04-graphene-fet-device-physics.md`,
   `thesis_draft/01-introduction.md`): Rewrote Section 4.6 to describe
   the corrected model and the literature re-examination, updated the
   Section 4.7 status table row (RF figures of merit: "In progress" with
   a known artifact -> "Complete"), corrected the Chapter 1 Section 1.2
   motivation paragraph's blanket "f_max lags f_T by an order of
   magnitude" claim to the more accurate "design/parasitics-dependent"
   framing, and updated the Chapter 1 status table's Chapter 4 row.

**Not yet covered (candidates for future runs):**
- Using the contact-doping model (Section 4.5) to recalibrate per-metal
  Rc values (still open, Section 4.7)
- Plasmonic-absorption-enhancement factor for the photodetector model
  (Chapter 6, Section 6.5)
- Integrating the contact-doping spatial-profile machinery into the
  photodetector collection-efficiency model (Chapter 6, Section 6.5)
- Graphene-all-around-metal (liner/cap) interconnect model (Chapter 5)
- Copper comparison model liner/barrier-thickness effect (Chapter 5,
  Section 5.2)
- Confirming the "High fMAX/fT ratio in multi-finger embedded T-shaped
  gate graphene transistors" source referenced by title only this
  session (fetch was rate-limited) — would strengthen the multi-finger
  R_g discussion with a second, independent data point
- Further thesis_draft/ chapters: 2-3 could be drafted from the existing
  band-structure/transport code and results; Chapter 7 (discussion/
  outlook) awaits enough device-application chapters to synthesize

**Commits this run:** 3 (fmax-parasitics research notes, code fix +
regenerated plot, thesis_draft/ Chapter 4 Section 4.6 + Chapter 1
updates; this log entry commit makes 4).

---

## 2026-08-28

**Status:** Eighth automation run. No run is recorded for 2026-08-27 in
this log or in `git log` — the previous entry (2026-08-26, second
session) is the most recent prior activity, so this run picked up
directly from its "not yet covered" list rather than assuming an
intervening day's work exists.

**Repo state at start:** Chapters 1, 4, 5, 6 of `thesis_draft/` exist.
Chapter 5 (`thesis_draft/05-graphene-interconnects.md`) had, since
2026-08-23, explicitly flagged its copper comparison baseline as
surface-scattering-only and therefore "conservative in graphene's
favor," omitting the diffusion-barrier/adhesion-liner-thickness effect
— listed as a Section 5.4 follow-on item and repeated in every
subsequent "not yet covered" list through 2026-08-26. Unlike some other
open items (the graphene-all-around-metal hybrid architecture, the
plasmonic-enhancement factor), this one is a direct correction to an
already-implemented, already-drafted model rather than a new feature, so
it was prioritized this session.

**Work done:**

1. **Research notes**
   (`notes/2026-08-28-copper-liner-barrier-thickness-effect.md`):
   WebSearch and WebFetch were both available and used. Found
   Domenichini et al., "Selecting alternative metals for advanced
   interconnects," arXiv:2406.09106 (2024), which states the combined
   Cu barrier+liner thickness cannot scale below ~2-3 nm without losing
   function, describes in words exactly the effect Chapter 5 had
   flagged as missing (liner/barrier occupying an increasing volume
   fraction of shrinking cross-section), and gives an IRDS-based
   linewidth roadmap (23 nm in 2024/25 down to 12 nm by 2035). Also
   found "Mechanisms of Scaling Effect for Emerging Nanoscale
   Interconnect Materials," *Nanomaterials* 12(10), 1760 (2022), which
   gives concrete per-conductor liner thicknesses (Cu/TaN-Co: 3 nm; Ru:
   0.3 nm; Co: 1 nm; W: linerless) and states that below ~20 nm, Cu's
   resistance advantage over alternative conductors is significantly
   weakened by this effect. Neither source publishes a closed-form
   effective-resistivity formula (both use full resistance/Monte Carlo
   simulation instead), so the note derives an original simplified
   area-dilution model from their reported physical picture and
   thickness figures, stated as such rather than presented as a
   literature-derived equation.

2. **Code** (`graphene_interconnect_model.py`,
   `interconnect_resistivity_vs_linewidth.png`): Added
   `cu_resistivity_with_liner()`, treating the barrier/liner as
   non-conducting and consuming 3 nm from both in-plane dimensions of
   the drawn wire (literature value from the Nanomaterials review),
   giving an effective Cu core of `W_eff = W - 2*t_liner` whose
   resistivity is both computed via the existing surface-scattering
   model at the narrower width *and* diluted by the drawn/conducting
   area ratio squared. Below `W = 6 nm` (the liner-consumption floor)
   the function correctly returns NaN rather than a finite value, since
   no Cu core can physically exist there. Kept as an explicit second
   curve alongside (not replacing) the original surface-scattering-only
   Cu model, and extended `summary_numbers()` to report crossover widths
   against both. Verified to run standalone; regenerated the plot with
   both Cu curves overlaid. Result: at realistic diffuse-edge quality
   (p = 0.15), graphene now stays *below* the liner-aware Cu model
   across the entire valid modeled range (no crossover) — a materially
   different conclusion from the original ~130 nm crossover found
   against the surface-scattering-only baseline, though both remain
   simplified analytic models rather than full transport simulations.

3. **Writing** (`thesis_draft/05-graphene-interconnects.md`,
   `thesis_draft/01-introduction.md`): Rewrote Section 5.2 to describe
   both copper baselines and Section 5.3 to report the new crossover
   result with representative resistivity values at roadmap-relevant
   linewidths (12-23 nm), explicitly framing this as a qualitatively
   different conclusion from the earlier draft rather than a minor
   refinement, while being careful not to overclaim beyond what a
   simplified analytic model on both the graphene and copper sides can
   support. Marked the Section 5.4 liner-effect follow-on item done and
   added two narrower open items it surfaces (parallel liner conduction;
   a thinner Ru/Co-liner Cu scenario). Updated the Chapter 1 status
   table's Chapter 5 row.

**Not yet covered (candidates for future runs):**
- Parallel-conduction refinement to the liner/barrier model (currently
  treats the liner as strictly non-conducting; Section 5.4)
- Thinner-liner (Ru- or Co-liner-enabled Cu, ~1 nm or less) scenario,
  directly explorable via the new `t_liner_nm` parameter but not yet
  run/discussed (Section 5.4)
- Graphene-all-around-metal (liner/cap on a conventional Co/Cu core)
  interconnect model, distinct from the pure-GNR-wire model and from
  today's Cu-liner-effect correction (Chapter 5, Section 5.4)
- Using the contact-doping model (Chapter 4, Section 4.5) to recalibrate
  per-metal Rc values (still open, Section 4.7)
- Plasmonic-absorption-enhancement factor for the photodetector model
  (Chapter 6, Section 6.5)
- Integrating the contact-doping spatial-profile machinery into the
  photodetector collection-efficiency model (Chapter 6, Section 6.5)
- Confirming the "High fMAX/fT ratio in multi-finger embedded T-shaped
  gate graphene transistors" source referenced by title only in the
  2026-08-26 session (fetch was rate-limited at the time)
- Further thesis_draft/ chapters: 2-3 could be drafted from the existing
  band-structure/transport code and results; Chapter 7 (discussion/
  outlook) awaits enough device-application chapters to synthesize

**Commits this run:** 3 (liner-thickness research notes, code +
regenerated plot, thesis_draft/ Chapter 5 + Chapter 1 updates; this log
entry commit makes 4).
