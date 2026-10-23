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

---

## 2026-08-29

**Status:** Ninth automation run.

**Repo state at start:** Chapters 1, 4, 5, 6 of `thesis_draft/` exist.
`notes/2026-08-28-copper-liner-barrier-thickness-effect.md` (Section 6)
had flagged two specific, narrowly-scoped open items from that session's
liner/barrier work: (1) the liner was modeled as strictly non-conducting,
a simplification both source papers use but do not claim is exact; (2) a
thinner-liner (Ru- or Co-liner-enabled Cu) scenario was directly
explorable via the existing `t_liner_nm` parameter but not yet run or
discussed. Both were prioritized this session since they are direct,
already-scoped corrections/extensions to an existing model rather than
new topics, consistent with this repo's established practice of closing
flagged gaps before opening new ones.

**Work done:**

1. **Code + results** (`graphene_interconnect_model.py`,
   `interconnect_resistivity_vs_linewidth.png`): Added
   `cu_resistivity_with_liner_parallel()`, an explicit core-plus-liner
   parallel-conduction model that generalizes the 2026-08-28
   non-conducting-liner formula. Verified against two limiting-case
   self-checks run automatically at script start (not just asserted in a
   docstring): rho_liner -> infinity reproduces the non-conducting model
   to 12 significant figures; rho_liner == rho_core collapses to the
   plain core resistivity independent of liner thickness. Used the new
   model to evaluate the previously-flagged thinner-liner scenario for
   two materials -- Ru (0.3 nm) and Co (1.0 nm), thicknesses from the
   same Nanomaterials 12(10), 1760 (2022) review already cited on
   2026-08-28 -- with representative liner resistivities (TaN 400, Ru 30,
   Co 20 uOhm.cm). WebSearch/WebFetch were both available and used to try
   to source specific thin-film liner resistivity values; multiple
   promising sources (a PMC article, several ScienceDirect abstracts, a
   re-fetch of the already-cited Domenichini et al. arXiv HTML) failed
   (reCAPTCHA redirect, robots.txt blocks, or missing the specific
   numeric section) and are recorded as such in today's notes rather than
   silently worked around -- search did reliably confirm the qualitative
   TaN >> Ru, Co resistivity ordering across independent sources (imec's
   public Ru/Co liner materials, the Domenichini framing text that did
   load), so the values used are flagged explicitly as order-of-magnitude
   placeholders pending a future session's sourcing, not as literature
   figures. Regenerated the plot with all three new curves (TaN parallel,
   Ru parallel, Co parallel) alongside the existing four. Results: the
   TaN parallel-conduction correction is nearly negligible (TaN is still
   ~100x more resistive than the Cu core at relevant widths, so the
   2026-08-28 "no crossover" conclusion is unchanged); the thin-Ru case
   tracks close to the original liner-free baseline (crossover ~123 nm,
   vs. 130 nm liner-free); the thin-Co case is qualitatively different --
   because current can still flow through the conductive liner as the Cu
   core vanishes near the W = 2*t_liner floor, resistance stays low there
   instead of diverging, pushing the realistic-edge graphene crossover
   down to ~3 nm. This last result directly reproduces, from the model
   itself, the qualitative industry rationale (found in this session's
   search) for pursuing thin, low-resistance Ru/Co liners in the first
   place.

2. **Research notes**
   (`notes/2026-08-29-liner-parallel-conduction-and-thin-liner-scenarios.md`):
   Documents the parallel-conduction derivation, the search-access
   failures and what was/wasn't established from them, the representative
   liner-resistivity table and its caveats, and the full results
   discussion summarized above.

3. **Writing** (`thesis_draft/05-graphene-interconnects.md`,
   `thesis_draft/01-introduction.md`): Added a new discussion paragraph to
   Section 5.3 covering the parallel-conduction and thin-liner results,
   marked both Section 5.4 follow-on items done, and added two narrower
   open items they surface (sourcing the specific liner-resistivity
   values; spatial non-uniformity within the liner frame). Updated the
   Chapter 1 status table's Chapter 5 row to reflect that the
   realistic-edge crossover conclusion is now liner-choice-dependent
   (ranging from "no crossover" to "~3 nm") rather than a single number.

**Not yet covered (candidates for future runs):**
- Sourcing the specific liner-material resistivity values (TaN, Ru, Co)
  used today against a quantitative reference, rather than the
  order-of-magnitude placeholders used (Chapter 5, Section 5.4)
- Graphene-all-around-metal (liner/cap on a conventional Co/Cu core)
  interconnect model, distinct from the pure-GNR-wire model and from all
  of this session's and 2026-08-28's Cu-liner corrections (Chapter 5,
  Section 5.4)
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

**Commits this run:** 2 (parallel-conduction model code + notes +
regenerated plot, thesis_draft/ Chapter 5 + Chapter 1 updates; this log
entry commit makes 3).

## 2026-08-31

**Status:** Tenth automation run (first run after a two-day gap; no
2026-08-30 session occurred).

**Repo state at start:** Chapters 1, 4, 5, 6 of `thesis_draft/` exist.
Chapter 1's status table (Chapter 4 row) and `graphene_contact_doping_
model.py`'s module docstring both flagged the same standing open item
since 2026-08-21: "use the doping model to recalibrate per-metal Rc."
This was the only concretely-scoped, not-yet-started item in the
Chapter 1 status table's "still to add" lists (the interconnect and
photodetector rows' open items are all either done or require new
model architecture, not a recalibration of an existing one), so it was
prioritized this session per the repo's established practice of closing
already-flagged gaps before opening new topics.

**Work done:**

1. **Research notes**
   (`notes/2026-08-31-per-metal-contact-resistance-literature-and-
   recalibration.md`): WebSearch/WebFetch were available and used.
   Sourced four individually-cited, per-metal, room-temperature,
   top-contact, on-state measured Rc values: Cu 184 Ω·µm and Pd 584
   Ω·µm (Smith et al., ACS Nano 7(4), 3661 (2013)), Ni 470 Ω·µm
   (Khosravi Rad et al., Sci. Rep. 14, 9190 (2024) -- this is the same
   470 Ω·µm figure already cited generically in the 2026-08-21 notes,
   now confirmed as specifically a Ni value), and Au 519 Ω·µm (Passi et
   al., arXiv:1807.04772). Two fetch attempts failed and are recorded
   as such rather than silently worked around: pmc.ncbi.nlm.nih.gov
   returned a reCAPTCHA interstitial (blocked a PMC review article);
   researchgate.net returned HTTP 429 (rate-limited) on both a
   Ti-specific paper and an Au temperature-dependence paper -- as a
   result, Ti and Cr were NOT recalibrated this session (no clean
   individually-sourced figure obtained). Also documents a same-paper,
   same-metal, geometry-only Au top-vs-edge-contact comparison from
   Passi et al. (519 vs. 45 Ω·µm, an ~11x reduction) as a secondary
   finding.

2. **Code + results**
   (`graphene_contact_doping_model.py`, `rc_recalibration.png`): Added
   `METAL_LITERATURE_RC` (the four sourced values + citations/
   conditions), added a Ni work function (5.04 eV) to
   `METAL_WORK_FUNCTIONS`, and added `recalibrate_metal_rc()` /
   `print_rc_recalibration()` / `plot_rc_recalibration()`. The
   originally planned additive decomposition (Rc,measured = R_extra,
   computed doping-junction term + R_transmission, an unmodeled
   interface term) **does not hold**: for 3 of 4 metals (Cu, Ni, Au)
   the computed R_extra alone exceeds the entire measured literature
   Rc, forcing an unphysical negative "naive R_transmission." Only Pd
   gives a plausible small positive residual (~9% of its measured
   total). This is reported and plotted as a genuine negative/
   inconclusive result (grouped-bar comparison, not a stacked
   decomposition that would visually overstate validity) rather than
   adjusted (e.g. by retuning `lambda_decay`) to force agreement with
   only four data points -- candidate explanations (TLM
   double-counting of the same near-contact region; `lambda_decay`
   mismatch with the literature devices' real geometry) are recorded
   as open items, not resolved.

3. **Writing** (`thesis_draft/04-graphene-fet-device-physics.md`,
   `thesis_draft/01-introduction.md`): Added new Section 4.7
   ("Per-metal Rc recalibration: a genuine negative-residual result"),
   renumbered the prior 4.7 summary table to 4.8 and updated its
   recalibration row to reflect the attempted-but-inconclusive outcome
   rather than marking it simply "done." Also added the Au edge-vs-top
   contact secondary finding. Updated the Chapter 1 status table's
   Chapter 4 row accordingly.

**Not yet covered (candidates for future runs):**
- Isolating the root cause of the Section 4.7 negative-residual result
  (TLM double-counting vs. `lambda_decay` mismatch) -- needs either the
  literature devices' extracted transfer length L_T or a direct
  simulation of the TLM extraction procedure on top of the doping
  profile
- Ti and Cr per-metal Rc recalibration (blocked this session by
  ResearchGate rate-limiting; retry a future session, possibly via a
  different source than the two ResearchGate mirrors that failed)
- Contact-geometry dependence (top vs. edge contact) is not represented
  at all in `graphene_contact_doping_model.py`, despite Section 4.7's
  own finding that it is a larger lever (~11x) than metal choice
  (~3x) at fixed geometry
- Plasmonic-absorption-enhancement factor for the photodetector model
  (Chapter 6, Section 6.5)
- Integrating the contact-doping spatial-profile machinery into the
  photodetector collection-efficiency model (Chapter 6, Section 6.5)
- Confirming the "High fMAX/fT ratio in multi-finger embedded T-shaped
  gate graphene transistors" source referenced by title only in the
  2026-08-26 session (fetch was rate-limited at the time)
- Graphene-all-around-metal (liner/cap on a conventional Co/Cu core)
  interconnect model (Chapter 5, Section 5.4); sourcing specific
  liner-material resistivity values (Chapter 5, Section 5.4)
- Further thesis_draft/ chapters: 2-3 could be drafted from the
  existing band-structure/transport code and results; Chapter 7
  (discussion/outlook) awaits enough device-application chapters to
  synthesize

**Web search availability:** WebSearch/WebFetch were both available and
used this session (see Section 3 of today's notes for specific access
failures encountered and worked around by omission, not substitution).

**Commits this run:** 3 (per-metal literature Rc research notes,
recalibration code + plot + genuine negative-residual finding,
thesis_draft Chapter 4 Section 4.7 + Chapter 1 status table update;
this AUTOMATION_LOG.md entry commit makes 4).

---

## 2026-09-05

**Status:** Eleventh automation run (first run after a five-day gap; no
sessions occurred 2026-09-01 through 2026-09-04).

**Repo state at start:** Chapters 1, 4, 5, 6 of `thesis_draft/` exist.
Every "not yet covered" list since the 2026-08-31 run flagged the same
concrete gap: "Contact-geometry dependence (top vs. edge contact) is not
represented at all in `graphene_contact_doping_model.py`, despite Section
4.7's own finding that it is a larger lever (~11x) than metal choice (~3x)
at fixed geometry." This was prioritized this session since it is a
direct, already-scoped closure of a standing gap (consistent with this
repo's established practice) rather than a new topic, and the ~11x figure
already lived in this repo's own Section 4.7 text without a model behind
it.

**Work done:**

1. **Research notes**
   (`notes/2026-09-05-edge-vs-top-contact-geometry.md`): WebSearch and
   WebFetch were both available and used. Reviewed Wang et al., *Science*
   342, 614 (2013) (foundational true-1D edge-contact result, ~100 Ohm.um,
   direct fetch 403'd, corroborated via a secondary source rather than
   invented) and re-examined Passi et al., arXiv:1807.04772 (already
   partially cited in this repo for its Au top-contact value) for its
   edge/hole-pattern TLM table (5 hole diameters, 50-1000 nm, at both the
   Dirac point and on-state) and its DFT-computed Fermi-level shift
   comparison (0.35 eV at a graphene edge vs. 0.14 eV at the flat surface,
   for Au) -- the one quantitative, metal-specific, first-principles
   number found that converts directly into a carrier-density input for
   this repo's existing contact-doping machinery. A Wiley fetch of a
   second independent edge-contact paper (Lee et al. 2022) failed (403)
   and is recorded as not pursued further this session, not silently
   substituted.

2. **Code + results**
   (`graphene_contact_doping_model.py`, `graphene_edge_contact_model.py`,
   `edge_vs_top_contact.png`): Factored `junction_extra_resistance()`'s
   doping-profile integration out into a new
   `junction_extra_resistance_from_ncontact()` (no behavior change --
   verified by re-running the existing module standalone and confirming
   identical output to before the refactor) so a new module could reuse it
   with an n_contact computed a different way. `graphene_edge_contact_model.py`
   converts the two DFT Fermi shifts into carrier densities via graphene's
   linear dispersion (n ~ E_F^2), computes pure edge-mode vs. pure top-mode
   extra junction resistance for Au (379 vs. 982 Ohm.um -- edge lower, as
   expected, but only ~2.6x, not the full ~11x measured on Passi et al.'s
   actual patterned device), and adds a geometric hole-array model relating
   hole diameter/areal fill fraction to an "edge-influenced area fraction"
   via the existing `lambda_decay` (250 nm) parameter -- reusing an
   already-established length scale rather than introducing a new free
   parameter. Verified to run standalone. The patterned-contact sweep
   reproduces the large-hole-diameter (area-loss-dominated) branch of
   Passi et al.'s real non-monotonic data, but explicitly does NOT
   reproduce their small-diameter (50-100 nm) upturn -- flagged in both the
   code's `summary_numbers()` output and the notes as an open modeling
   gap (most plausibly a fabrication/lithography effect outside this
   compact model's scope) rather than hidden or fitted away. Passi et
   al.'s actual hole areal fill fraction is not reported in the fetched
   material, so the comparison is explicitly qualitative (swept over a
   few assumed fill fractions), not a point-by-point fit.

3. **Writing** (`thesis_draft/04-graphene-fet-device-physics.md`,
   `thesis_draft/01-introduction.md`): Added new Section 4.8 ("Contact
   geometry: edge vs. top contacts, quantitatively") plus a
   Section 4.8.1 on the patterned-contact geometry model, renumbered the
   prior 4.8 ("Summary and open items") to 4.9 and updated its status
   table with a new row for this session's work, corrected a
   now-outdated sentence in Section 4.7 that said contact geometry "is
   not modeled quantitatively here" to point forward to the new Section
   4.8, and updated the Chapter 1 status table's Chapter 4 row.

**Not yet covered (candidates for future runs):**
- Sourcing a second independent edge-contact dataset (Lee et al. 2022,
  Wiley 403'd this session) to cross-check the Au DFT-Fermi-shift-derived
  edge/top ratio against a source other than Passi et al.
- Extending the edge-mode carrier-density estimate to metals other than
  Au (no DFT edge-vs-surface Fermi-shift data found for other metals this
  session; would require either new literature or an untested assumption
  that Au's ~2.5x shift ratio generalizes)
- Isolating the small-hole-diameter (50-100 nm) upturn in Passi et al.'s
  data that this session's geometric model does not reproduce
- Isolating the root cause of the Section 4.7 negative-residual result
  (TLM double-counting vs. `lambda_decay` mismatch) -- still open from
  2026-08-31
- Ti and Cr per-metal Rc recalibration (still blocked by prior sessions'
  ResearchGate rate-limiting; not reattempted this session)
- Plasmonic-absorption-enhancement factor for the photodetector model
  (Chapter 6, Section 6.5)
- Integrating the contact-doping spatial-profile machinery into the
  photodetector collection-efficiency model (Chapter 6, Section 6.5)
- Graphene-all-around-metal (liner/cap on a conventional Co/Cu core)
  interconnect model (Chapter 5, Section 5.4); sourcing specific
  liner-material resistivity values (Chapter 5, Section 5.4)
- Further thesis_draft/ chapters: 2-3 could be drafted from the existing
  band-structure/transport code and results; Chapter 7 (discussion/
  outlook) awaits enough device-application chapters to synthesize

**Web search availability:** WebSearch/WebFetch were both available and
used this session (see notes file Section 5 for the two access failures
encountered: science.org 403, Wiley 403, and a Semantic Scholar fetch that
returned an empty page body).

**Commits this run:** 3 (edge-vs-top contact geometry research notes,
model code + plot, thesis_draft Chapter 4 Section 4.8 + Chapter 1 status
table update; this AUTOMATION_LOG.md entry commit makes 4).

## 2026-09-06

**Status:** Thirteenth automation run.

**Repo state at start:** Chapters 1, 4, 5, 6 of `thesis_draft/` exist.
Every "not yet covered" list since the photodetector model was first
implemented (2026-08-24/2026-08-25) flagged the same standing item:
"Plasmonic-absorption-enhancement factor for the photodetector model
(Chapter 6, Section 6.5)," and `graphene_photodetector_model.py`'s own
module docstring names it as the natural next step for closing the
absolute-magnitude offset between its simple transit-time gain-bandwidth
model and the three literature anchor points. This session closed that
item.

**Work done:**

1. **Research notes**
   (`notes/2026-09-06-plasmonic-enhancement-graphene-photodetectors.md`):
   WebSearch and WebFetch were both available and used; all fetches
   succeeded (no access failures to record this session). Reviewed
   Echtermeyer et al., Nature Communications 2, 458 (2011) (Ti/Au
   finger-grating plasmonic enhancement of graphene photovoltage, ~5x
   near-field amplitude / ~25x intensity enhancement, resonances at 514
   and 633 nm), Fang et al., Applied Physics Letters 105, 241114 (2014)
   (Au nano-antenna, 580 nm LSPR), and an arrayed bowtie-on-waveguide
   telecom-band graphene photodetector (arXiv:1808.10823, 8.5x simulated
   single-element absorption enhancement, 100 Gbit/s PAM-2/PAM-4
   reception). Explicitly distinguished which numbers are directly usable
   as a near-field/absorption enhancement factor (Echtermeyer's near-field
   simulation, the bowtie array's same-device absorption comparison) from
   one that is not (Fang et al.'s "four orders of magnitude" figure, a
   device-to-device responsivity comparison bundling in unrelated contact/
   bias/collection differences) -- following this repo's established
   practice of not force-fitting a not-apples-to-apples literature number
   into a compact model.

2. **Code + results**
   (`graphene_plasmonic_photodetector_model.py`,
   `plasmonic_photodetector_enhancement.png`): New module reusing
   `graphene_photodetector_model.py`'s `EQE_BARE`, `TAU_TRANSIT`,
   `responsivity_bare()`, and `photoconductive_gain()` unmodified. Models
   plasmonic enhancement as a Lorentzian intensity-enhancement factor
   applied to EQE as a function of wavelength, calibrated against the two
   directly-usable literature designs (F_max=25x at 514/633 nm;
   F_max=8.5x at 1550 nm). Resonance quality factor Q is an explicitly
   assumed (not measured or fitted) representative value for lossy Au/Ti
   nanostructures, since neither source paper reported FWHM in the
   material reviewed -- stated in the code's docstring and printed
   alongside the resulting FWHM in `summary_numbers()`'s output rather
   than hidden inside the Lorentzian shape. Verified to run standalone:
   on-resonance responsivity rises from 0.62-1.88 mA/W (bare) to
   15.6-19.2 mA/W (plasmonic) across the three resonances plotted.

3. **Writing** (`thesis_draft/06-graphene-photodetectors.md`,
   `thesis_draft/01-introduction.md`): Added new Section 6.5 ("Plasmonic
   near-field absorption enhancement"), folding in and completing the
   former Section 6.5 ("Further follow-on work"), renumbered to 6.6, with
   its first bullet (the plasmonics item) removed since it is now done
   and three new bullets added (spatial integration with the contact-
   doping machinery, a photo-bolometric mechanism model, and sourcing a
   measured Q). Updated the References section and the Chapter 1 status
   table's Chapter 6 row accordingly.

**Not yet covered (candidates for future runs):**
- No measured plasmon resonance FWHM/Q was found for either fitted
  design this session; the assumed Q=6-8 values remain open to
  replacement with literature-anchored ones if a future session finds
  supplementary simulation data reporting linewidth
- Integrating Section 6.5's spatially-uniform EQE-multiplier picture with
  the spatially-resolved contact-doping machinery
  (`graphene_contact_doping_model.py`, `graphene_edge_contact_model.py`)
  -- still open, now stated in both Chapter 4 and Chapter 6
- The photo-bolometric mechanism (dominant in the telecom bowtie/
  waveguide device) is not modeled anywhere in this repo -- a legitimate
  separate follow-on from the photogating-only picture currently modeled
- Isolating the root cause of the Section 4.7 negative-residual result
  (TLM double-counting vs. `lambda_decay` mismatch) -- still open from
  2026-08-31
- Ti and Cr per-metal Rc recalibration (still blocked by prior sessions'
  ResearchGate rate-limiting)
- Sourcing a second independent edge-contact dataset (Lee et al. 2022,
  Wiley 403'd 2026-09-05) -- still open
- Graphene-all-around-metal (liner/cap) interconnect model (Chapter 5,
  Section 5.4); sourcing specific liner-material resistivity values
- Further thesis_draft/ chapters: 2-3 could be drafted from the existing
  band-structure/transport code and results; Chapter 7 (discussion/
  outlook) awaits enough device-application chapters to synthesize

**Web search availability:** WebSearch/WebFetch were both available and
used this session; all fetches succeeded.

**Commits this run:** 3 (plasmonic-enhancement research notes, model code
+ plot, thesis_draft Chapter 6 Section 6.5 + Chapter 1 status table
update; this AUTOMATION_LOG.md entry commit makes 4).

## 2026-09-07

**Status:** Fourteenth automation run.

**Repo state at start:** Chapters 1, 4, 5, 6 of `thesis_draft/` exist.
Chapter 6's Section 6.6 (follow-on work) had listed, since the chapter
was first drafted on 2026-08-24, the same standing item repeated through
the 2026-09-06 plasmonic-enhancement session: "Replace the lumped
EQE_bare = 0.15% parameter with a spatially resolved diffusion-length
collection model once the Chapter 4 spatially-resolved contact-doping-
profile follow-on... is implemented" -- that Chapter 4 follow-on
(`graphene_contact_doping_model.py`) has existed since 2026-08-26 but had
never been integrated into the photodetector side. This session closed
that item.

**Work done:**

1. **Research notes**
   (`notes/2026-09-07-spatial-photocarrier-collection-model.md`):
   WebSearch and WebFetch were both available and used. Two of four
   fetch attempts succeeded: Xia, Mueller, Lin, Valdes-Garcia & Avouris,
   "Ultrafast graphene photodetector," *Nature Nanotechnology* 4, 839
   (2009) (zero-bias photodetection via asymmetric-work-function Ti/Pd
   contacts -- direct experimental evidence that contact metal sets a
   usable built-in field), and Weiss & Duan, "Building potential for
   graphene photodetectors," *NPG Asia Materials* 5, e74 (2013) (states
   the work-function-mismatch/potential-offset/photocarrier-separation
   mechanism explicitly, and the important caveat that identical
   contacts give "equal and opposing" potentials that partially cancel
   at zero bias). Two fetch attempts failed and are documented rather
   than silently worked around: a PMC review article returned a Google
   reCAPTCHA interstitial (consistent with this repo's prior PMC access
   failures, e.g. 2026-08-31), and a PubMed abstract (19326919, the
   scanning-photocurrent-microscopy sign-map paper) returned HTTP 429
   rate-limiting and was not retried. Also cross-checks the existing
   `lambda_decay = 250 nm` (Khomyakov et al. 2010, already used in
   `graphene_contact_doping_model.py`) against an independently-sourced
   ~100-200 nm built-in-field collection-region estimate already cited
   in this repo's own `notes/2026-08-24-photodetector-responsivity.md`
   (Mueller, Xia & Avouris 2010) -- same order of magnitude, different
   source, reported as a useful but imperfect cross-check rather than a
   match.

2. **Code + results**
   (`graphene_photodetector_collection_model.py`,
   `photodetector_collection_efficiency_by_metal.png`): New module
   superposing the contact-doping-induced field (reusing
   `graphene_contact_doping_model.py`'s `lambda_decay` and
   `METAL_WORK_FUNCTIONS`, and `graphene_fet_model.py`'s mobility) on the
   existing bare-device model's uniform bias field
   (`graphene_photodetector_model.py`'s `V_bias`/`L_channel`), integrating
   the resulting position-dependent drift velocity to get each
   photocarrier's transit time to the contact, then an exponential-
   lifetime collection (survival) probability against graphene's ~1 ps
   photocarrier lifetime (Mueller, Xia & Avouris 2010). Verified to run
   standalone. Result: collection-efficiency enhancement (relative to a
   work-function-matched, ~zero-doping baseline -- physically Cr, already
   present in `METAL_WORK_FUNCTIONS`) ranges from 1.00x (Cr) to 1.47x
   (Pt, largest work-function mismatch); common real contact metals Ti
   and Cu give ~1.21-1.23x, and the higher-work-function metals already
   favored in Chapter 4 for low contact resistance (Ni, Au, Pd) give
   ~1.39-1.41x. As a sanity check, the model's bias-field-only baseline
   transit time coincides almost exactly (1.0 ps) with the independently
   assumed carrier lifetime, giving a physically sensible order-unity
   baseline collection efficiency (1 - 1/e ~= 0.632) before any
   doping-field enhancement is added. Explicitly scoped as a
   single-contact (reinforcing-field-only) model -- the real two-terminal
   device's second, partially-cancelling contact (per Weiss & Duan's
   "equal and opposing" result) is not modeled, and this is stated as an
   open item rather than approximated away silently.

3. **Writing** (`thesis_draft/06-graphene-photodetectors.md`,
   `thesis_draft/01-introduction.md`): Added new Section 6.6 (model,
   literature mechanism support, results, explicit scope limitation),
   renumbered the prior Section 6.6 ("Further follow-on work") to 6.7
   with its first bullet resolved and a new bullet added (self-consistent
   two-contact modeling), added a bullet noting the unretrieved PubMed
   19326919 abstract as a future cross-check candidate, updated the
   References section, and updated the Chapter 1 status table's Chapter
   6 row.

**Not yet covered (candidates for future runs):**
- Self-consistent two-contact collection model (Section 6.6's single
  reinforcing-contact scope limitation -- the natural next step for this
  session's new model)
- Revisit Park, Ahn et al.'s scanning-photocurrent-microscopy sign-map
  paper (PubMed 19326919, WebFetch 429'd this session) for a more direct
  literature cross-check of the two-contact sign-reversal reasoning
- Integrating Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery now used in both Section
  6.6 and Chapter 4 -- still open
- The photo-bolometric mechanism (dominant in the telecom bowtie/
  waveguide device, Section 6.5) is not modeled anywhere in this repo
- Isolating the root cause of the Section 4.7 negative-residual result
  (TLM double-counting vs. `lambda_decay` mismatch) -- still open from
  2026-08-31
- Ti and Cr per-metal Rc recalibration (still blocked by prior sessions'
  ResearchGate rate-limiting; not reattempted this session)
- Sourcing a second independent edge-contact dataset (Lee et al. 2022,
  Wiley 403'd 2026-09-05) -- still open
- Isolating the small-hole-diameter (50-100 nm) upturn in Passi et al.'s
  patterned-contact data (Section 4.8.1) -- still open
- Graphene-all-around-metal (liner/cap) interconnect model (Chapter 5,
  Section 5.4); sourcing specific liner-material resistivity values
- Further thesis_draft/ chapters: 2-3 could be drafted from the existing
  band-structure/transport code and results; Chapter 7 (discussion/
  outlook) awaits enough device-application chapters to synthesize

**Web search availability:** WebSearch/WebFetch were both available this
session; 2 of 4 fetch attempts succeeded (see Section 2 of today's notes
file for the two failures: a PMC reCAPTCHA block, a PubMed 429).

**Commits this run:** 3 (spatial-collection-model research notes, model
code + plot, thesis_draft Chapter 6 Section 6.6 + Chapter 1 status table
update; this AUTOMATION_LOG.md entry commit makes 4).

---

## 2026-09-17 — Self-consistent two-contact collection: the result that reverses 2026-09-07

**Repo-history note:** the previous commit here is dated 2026-09-07, so
2026-09-08..09-16 has no sessions recorded. That was broken automation,
not a pause — see "Automation health" at the end of this entry.

**Work picked from the standing "not yet covered" list:** its first
entry, flagged on 2026-09-07 as "the natural next step" — the
self-consistent two-contact collection model.

**Work done:**

1. **Research** (`notes/2026-09-17-two-contact-self-consistent-
   collection.md`): literature basis for the two-contact problem. Weiss
   & Duan (*NPG Asia Materials* 5, e74 (2013)) re-fetched for the exact
   wording — "symmetric metal-graphene-metal devices generate an equal
   positive and negative flow with a net zero photocurrent", escaped by
   "using metals with asymmetric band structures". Suzuki et al.
   (*Carbon Trends* 5, 100115 (2021)) is the experimental counterpart
   and the more useful source, because its entire device strategy exists
   to defeat the cancellation: a 50 nm Ni shadow mask over one of two
   graphene/Ti interfaces, or unequal contact areas via a comb-shaped
   electrode; ~1.7e5 cm.Hz^(1/2)/W detectivity at 690 nm, zero bias.
   The notes state explicitly that **neither source gives a measured
   net-response-vs-work-function curve**, so this model is
   mechanism-anchored but *not* calibrated against experiment — recorded
   so a later session does not mistake the numbers below for validated
   predictions.

2. **Code** (`graphene_photodetector_two_contact_model.py`): signed
   total field F(x) = E_bias + g_A(x) - g_B(x), with g_A the existing
   single-contact doping-field profile (same `lambda_decay` = 250 nm,
   same `METAL_WORK_FUNCTIONS`) and g_B its mirror about mid-channel.
   The minus sign on g_B is the whole physical content that was
   missing. No assumed destination: a carrier drifts along F and is
   collected at A (+1) or B (-1) only if F keeps its sign along the
   whole path, and otherwise reaches a **stagnation point** and is
   counted uncollected. Note the geometry that makes this first-order
   rather than a correction: L = 200 nm against lambda_decay = 250 nm,
   so neither contact's field has decayed by mid-channel.

   **Validated against the model it supersedes, as required, rather
   than asserted:** switching off contact B's doping field must
   collapse this model onto 2026-09-07's `mean_collection_efficiency()`.
   It does, for all seven metals, worst relative deviation **7.8e-08**.
   Second validation, against physics rather than against the old code:
   at zero bias two identical contacts must give exactly zero net
   response by antisymmetry — measured worst |N| = **6.6e-17**. The
   antisymmetry of the pair matrix under contact swap is a third check.
   All three are computed and printed by `__main__`.

3. **Results** (`photodetector_two_contact_net_response.png`): three
   panels — the two opposing fields with the Pt/Pt stagnation point
   visible at x = 100 nm, the single-contact vs symmetric two-contact
   bar comparison, and zero-bias net response against W_A - W_B.

   **This contradicts the 2026-09-07 conclusion, and is reported as a
   contradiction rather than a refinement.** The single-contact model
   ranked metals Pt > Pd > Au > Ni > Ti > Cu > Cr, with Pt best at
   1.47x enhancement. For a *symmetric* two-terminal device at the same
   bias the ranking is essentially reversed — Cu ~ Ti > Cr > Ni > Au >
   Pd > Pt — with Pt now **worst** and overstated **5.5x** by the old
   model (0.929 -> 0.169). Mechanism: a strongly doping contact sweeps
   carriers toward itself, so two facing each other create a
   mid-channel stagnation point the 0.5 MV/m bias cannot overcome
   against Pt's ~4.5 MV/m near-contact field. Strong contact doping
   helps an isolated junction and hurts a symmetric device. Cr alone is
   unaffected (ratio exactly 1.00), its work function matching
   graphene's to 0.01 eV.

   Zero-bias asymmetric pairs give the Weiss & Duan mechanism in
   isolation (diagonal exactly zero): largest |N| is Cr/Pt at 0.917 —
   notably *not* Ti/Pt at 0.903, despite Ti/Pt's larger work-function
   difference, because Ti's own 0.17 eV mismatch sweeps carriers back
   toward Ti and partially opposes Pt's. Under this model the best
   zero-bias pairing is a strongly doping contact against a
   work-function-*matched* one.

4. **A bug the validation caught, worth recording** (notes Section 2.1):
   the first implementation obtained the toward-B transit time by
   subtracting one forward cumulative integral. Since 1/|v| diverges at
   a stagnation point, that gives inf-inf beyond any null and silently
   marked every carrier past a null uncollectable. Validation 1 passed
   anyway (with g_B = 0 there is no null). Validation 2 returned
   N = +0.46 for Pt/Pt at zero bias instead of 0 — plausible enough to
   rationalise as grid noise by anyone who wanted the model to work. It
   was the *exactly*-known zero demanded by antisymmetry that made the
   failure unambiguous. General lesson: a validation against an exactly
   known value beats one against a plausible range.

5. **Writing** (`thesis_draft/06-graphene-photodetectors.md`,
   `thesis_draft/01-introduction.md`): new Section 6.7 (model, both
   validations, both results, the simplification that could overturn
   Result 2), with the prior Section 6.7 renumbered to 6.8. Section 6.6
   was **not** quietly rewritten: its scope-limitation paragraph now
   carries an explicit callout that the limitation is resolved and the
   resolution overturns its ranking, and that its enhancement factors
   hold only for a single junction in isolation; its "what this can and
   cannot claim" paragraph gains a matching correction. The old numbers
   are retained because the single-junction calculation is the correct
   building block and the limiting case that validates the new model.
   Chapter 1's status table updated, including the Chapter 7 row, which
   now has a concrete synthesis task (see below).

**A real tension this created, for Chapter 7:** Chapter 4 favours
high-work-function metals (Ni, Au, Pd) for low contact resistance;
Section 6.7 finds those same metals are the worst for symmetric
two-terminal photoresponse. The two chapters now point in opposite
directions for the same design choice. That is genuine device physics,
not a modelling artefact, and Chapter 7 should resolve it explicitly
rather than let a reader notice it.

**Not yet covered (candidates for future runs):**
- **Signed, carrier-resolved two-contact treatment** — replace Section
  6.7's |W_metal - W_graphene| magnitude convention with one that
  distinguishes n-type (Ti, Cu) from p-type (Pt, Pd, Au) contacts and
  tracks electrons and holes separately. Now the most consequential
  open item in Chapter 6: it could invert Result 2's ordering by letting
  an n/p pair such as Ti/Pt add rather than partially cancel. Result 1's
  reversal does not depend on it (identical contacts, where the two
  conventions agree)
- Non-uniform illumination, a generation weight g(x) — Suzuki et al.'s
  shadow-mask device masks one of the two interfaces and is precisely
  a non-uniform-generation experiment, so it is the natural validation
  target
- **Chapter 7 (discussion/outlook)** is now genuinely actionable rather
  than waiting on more device chapters: it has a concrete Chapter 4 vs
  Section 6.7 contact-metal contradiction to synthesize
- Chapters 2-3 remain undrafted despite their band-structure/transport/
  optical computational results being complete — the largest remaining
  block of pure writing in the repo
- Integrating Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery — still open
- The photo-bolometric mechanism is still not modeled anywhere
- Isolating the root cause of the Section 4.7 negative-residual result
  (TLM double-counting vs. `lambda_decay` mismatch) — open since
  2026-08-31
- Ti and Cr per-metal Rc recalibration (blocked by ResearchGate
  rate-limiting; not reattempted)
- Second independent edge-contact dataset (Lee et al. 2022, Wiley 403'd
  2026-09-05) — still open
- Small-hole-diameter (50-100 nm) upturn in Passi et al.'s
  patterned-contact data (Section 4.8.1) — still open
- Graphene-all-around-metal (liner/cap) interconnect model (Chapter 5,
  Section 5.4) and specific liner-material resistivity values
- Park, Ahn et al.'s scanning-photocurrent-microscopy sign map
  (PubMed 19326919) — see the access note below before retrying

**Web search availability:** WebSearch available and used; 3 of 4
fetches succeeded. The failure: PubMed 34942635 (NEGF-DFT study of
asymmetric metal contacts on bilayer-graphene quantum dots) returned a
Google reCAPTCHA page. That is the **third** PubMed/PMC reCAPTCHA block
(also 2026-09-05, 2026-09-07), so it should now be treated as a standing
limitation of this environment rather than retried hopefully each
session — which also applies to the PubMed 19326919 item above. It
remains the best candidate for an independent first-principles
cross-check of asymmetric-contact photoresponse if another route to it
can be found.

**Automation health (needs Harsh's attention, reported separately):**
This session ran interactively, and diagnosed the 09-08..09-16 gap. The
scheduled task's device shell has no GitHub push credentials (`git push`
fails with `could not read Username for 'https://github.com'`; no
credential helper, no `gh`, and SSH cannot resolve github.com through
the HTTPS-only proxy), and the session VM is rebuilt per run so nothing
persists. Separately, `git` cannot run inside the connected folder at
all, because it must unlink its own `.git/*.lock` files and deletion in
connected folders is denied. Today's commits were made in the session's
own scratch clone and delivered as git bundles for Harsh to push
manually. Until credentials are resolved, unattended runs cannot push.
Also note `scipy` was absent from the device VM and had to be pip
installed for `graphene_contact_doping_model.py` to import.

**Commits this run:** 3 (two-contact literature/bug notes; the model +
figure; thesis_draft Chapter 6 Section 6.7 + Chapter 1 status table).
This AUTOMATION_LOG.md entry commit makes 4.

---

## 2026-09-18 — Signed, carrier-resolved contacts: one result confirmed and doubled, one retracted

**Status:** Full session. Closes the item that has sat at the top of
"not yet covered" since 2026-09-17 — the signed, carrier-resolved
two-contact treatment, described there as "the most consequential open
item in Chapter 6".

**Work done:**

1. **Research notes**
   (`notes/2026-09-18-signed-carrier-resolved-contact-fields.md`):
   literature basis for the signed treatment. The headline find is one
   the 2026-09-17 entry did **not** anticipate: the magnitude convention
   `|W_metal - W_graphene|` was hiding a *crossover* as well as a sign.
   The n/p crossover for a metal **on** graphene sits near **5.4 eV**,
   not at graphene's own 4.5 eV, because the short-range chemical
   interaction adds a ~0.9 eV interface dipole on top of vacuum-level
   alignment — Giovannetti, Khomyakov, Brocks, Karpan, van den Brink &
   Kelly, *Phys. Rev. Lett.* **101**, 026803 (2008),
   <https://arxiv.org/abs/0802.2267>. Under the physical crossover
   **six of the seven metals in `METAL_WORK_FUNCTIONS` n-dope graphene
   and only Pt p-dopes it**. The table's own inline comments already
   said so (Cu annotated "n-type dopant" despite 4.65 > 4.5; Pt
   "clearly above the ~5.4 eV crossover"); `|W - 4.5|` simply could not
   express it. Experimental anchor: Mueller, Xia & Avouris,
   *Nature Photonics* **4**, 297 (2010) — already cited in this repo for
   tau = 1 ps — built the Pd/Au vs Ti/Au interdigitated device precisely
   because identical electrodes give zero total photocurrent, reaching
   6.1 mA/W at 1.55 um, 15x. Supporting band-bending picture: Mueller
   *et al.*, arXiv:0902.1479 (0.12 eV potential step, doping extending
   0.2–0.3 um, p–n junction forming near the interface).

2. **Model** (`graphene_photodetector_signed_carrier_model.py`,
   `photodetector_signed_carrier_response.png`): signed offset
   dW = W_metal − w_cross drives a rigid Dirac-point shift near each
   contact; E(x) = (1/e) dE_D/dx; holes and electrons drift in that one
   field in opposite directions and are transported independently by the
   2026-09-17 stagnation-aware logic. N = net charge to contact A per
   photon, now over **[−2, +2]** because both carriers can be collected.
   The crossover is an explicit parameter and every result is quoted
   under both 4.5 eV and 5.4 eV rather than one being chosen silently.

3. **Three validations, all against exactly known values** (continuing
   the 2026-09-17 lesson that an exact check beats a plausible range):
   - reduction to the 2026-09-17 model (holes only, |dW|, crossover
     4.5 eV) across all 49 ordered pairs at both biases, 98
     comparisons — **0.000e+00, bitwise**. This is a true reduction, not
     an approximate agreement: under those restrictions the two modules
     are algebraically the same expression, E = −F_old.
   - symmetric pair at zero bias → |N| ≤ **6.6e−17**.
   - charge conjugation N(−dW) = −N(dW) at zero bias — **0.000e+00**.
     This one is only *statable* once the model is signed; it catches
     carrier-mixing errors that the symmetric check would pass.

4. **RESULT 2 — confirmed and roughly doubled.** The 2026-09-17 entry
   predicted this signed treatment "could invert Result 2's ordering by
   letting an n/p pair such as Ti/Pt add rather than partially cancel".
   It does. Best zero-bias pair goes from **Cr/Pt, |N| = 0.917** to
   **Ti/Pt, |N| = 1.832**. It is **insensitive to the crossover** —
   1.832 at 5.4 eV vs 1.830 at 4.5 eV — because Ti and Pt straddle both
   candidates, which makes it the most secure number in Chapter 6.
   Mechanism, from the stagnation audit: an unequal pair has **no
   interior field null at all** and both species are collected at
   opposite ends, while a same-type symmetric pair strands one entire
   species at the null and splits the other evenly. Ti/Pd reaches 1.599
   with both metals n-type, so the operative quantity is
   |dW_A − dW_B|, not the sign pair — which is the practically useful
   form, Pt being the only p-type metal in the table. Ti/Pd is also the
   pair Mueller *et al.* actually built.

5. **RESULT 1 — CONTRADICTS 2026-09-17, and retracts it.** The
   2026-09-17 entry asserted "Result 1's reversal does not depend on it
   (identical contacts, where the two conventions agree)". **That was
   wrong.** The conventions agree on the sign structure for identical
   contacts, but the signed model also collects the second carrier
   species, and the crossover choice reorders the metals. Symmetric
   device at V_bias = 0.1 V: Pt is **best (0.556)** under the 5.4 eV
   crossover and **worst (0.169)** under 4.5 eV. Three sections have now
   given three answers — 2026-09-07 Pt best, 2026-09-17 Pt worst,
   2026-09-18 either. The conclusion recorded is **not** "Pt is best
   after all" but that **the symmetric two-terminal metal ranking is not
   an established result of this repo**, because it flips with a
   modelling choice the earlier models never had to make explicit.
   Section 6.7's numbers are retained and annotated in place, since the
   arithmetic is correct for the convention it states.

   What survives, and is worth more than the ranking: under the physical
   crossover the symmetric response is **monotone in |dW|**. In a
   symmetric device the contact doping field is identical at both ends,
   contributes nothing to the net current, and does nothing but create a
   null that strands carriers. **For a symmetric two-terminal device
   contact doping is purely parasitic; the best metal is the one that
   perturbs graphene least.** Crossover-independent in form, and also why
   Cr (dW = 0 there) wins the 4.5 eV column.

6. **Writing** (`thesis_draft/06-graphene-photodetectors.md`,
   `thesis_draft/01-introduction.md`): new Section 6.8 with the full
   treatment (old 6.8 renumbered to 6.9). Section 6.7 annotated in
   place, not rewritten: a callout at Result 1 downgrading its reversal,
   a confirmation at Result 2, a note under the section heading, and no
   numbers deleted.

**Effect on Chapter 7's synthesis task (restated, not deleted).** The
2026-09-17 contradiction handed to Chapter 7 — Chapter 4's
low-contact-resistance metals (Ni, Au, Pd) being the worst for
photoresponse — **does not survive**: at the physical crossover Pd sits
near the top, and the symmetric ranking is no longer claimed at all. The
durable tension is narrower and more interesting: **Chapter 4 optimises
a single junction, while the two-terminal photoresponse depends on the
difference between two junctions and is to first order blind to either
one's own quality** (|N| = 1.832 for Ti/Pt vs at best 0.556 for any
symmetric pair). A second thread for the same chapter: the symmetric
design rule (dope graphene least) and the asymmetric one (maximise
|dW_A − dW_B|) point in opposite directions, and Chapter 5's
interconnect argument shares Chapter 4's single-junction framing.

**Methodological note worth carrying forward.** Two consecutive sessions
have now had a headline "reversal" overturned by the next session's
model. Both times the culprit was a convention that was stated as a
simplification but whose *consequences* were not enumerated — the
magnitude convention concealed a crossover nobody was looking for. The
practice that worked both times was the same: keep the superseded
numbers, state the contradiction in the commit message, the log and the
text, and prefer a structural statement (here: "doping is parasitic in a
symmetric device") over a ranking, because rankings are what flip.

**Not yet covered (candidates for future runs):**
- **Non-uniform illumination, a generation weight g(x)** — now the top
  open item in Chapter 6. Suzuki *et al.*'s shadow-mask device masks one
  of the two interfaces and is precisely a non-uniform-generation
  experiment, so it is the natural validation target. It is also the
  one remaining way a *symmetric* device can give nonzero response,
  which matters now that the symmetric metal ranking has been retracted
- **A non-linear dW → doping-profile relation.** Both the 2026-09-17 and
  2026-09-18 models take the profile magnitude as linear in dW with a
  single lambda for every metal. Giovannetti *et al.* find it only
  roughly linear and **not linear at all for the chemisorbed metals
  (Ti, Ni, Pd)** — which are exactly the metals carrying the largest
  |dW| under the physical crossover, so this is where the error is
  concentrated. Would test whether Ti/Pt's 1.832 is robust
- **Reconciling the 0.12 eV measured potential step** (Mueller *et al.*,
  arXiv:0902.1479) with the 0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS`
  assumes. The measurement is of the residual step in a gated device
  rather than the flat-band charge transfer, but the factor of 2–9 is
  unexplained in the repo
- **Chapter 7 (discussion/outlook)** — still actionable, with the
  restated single-junction-vs-two-junction synthesis above
- Chapters 2–3 remain undrafted despite their band-structure/transport/
  optical computational results being complete — still the largest
  remaining block of pure writing in the repo
- Whether Chapter 4's contact-resistance results should also be re-run
  at the 5.4 eV crossover. Section 4.5's doping profile uses the same
  `|W - W_graphene|` magnitude convention, and Chapter 4 has never been
  audited for it. **This is a new open item created by today's run and
  may affect published Chapter 4 numbers**
- Integrating Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery — still open
- The photo-bolometric mechanism is still not modeled anywhere
- Isolating the root cause of the Section 4.7 negative-residual result
  (TLM double-counting vs. `lambda_decay` mismatch) — open since
  2026-08-31
- Ti and Cr per-metal Rc recalibration (blocked by ResearchGate
  rate-limiting; not reattempted)
- Second independent edge-contact dataset (Lee *et al.* 2022, Wiley
  403'd 2026-09-05) — still open
- Small-hole-diameter (50–100 nm) upturn in Passi *et al.*'s
  patterned-contact data (Section 4.8.1) — still open
- Graphene-all-around-metal (liner/cap) interconnect model (Chapter 5,
  Section 5.4) and specific liner-material resistivity values
- Park, Ahn *et al.*'s scanning-photocurrent-microscopy sign map
  (PubMed 19326919) — PubMed remains reCAPTCHA-blocked, see below

**Web search availability:** WebSearch available and used. 3 of 3
fetches succeeded. One honest caveat: the `arxiv.org/abs/0802.2267`
fetch summary returned the n/p labels **inverted** (it claimed metals
above 5.4 eV n-dope graphene). The correct physics is the opposite — a
high-work-function metal withdraws electrons and p-dopes graphene — and
it was corrected against first principles before use. Had it been taken
at face value it would have flipped every sign in the new module. The
crossover *value* (5.4 eV) was the number actually used and is
independently confirmed by the repo's own table annotations. PubMed/PMC
was **not attempted**: three consecutive reCAPTCHA blocks (2026-09-05,
09-07, 09-17) established it as a standing environment limitation.

**Automation health:** Device reachable, `scipy` again absent from the
device VM and pip-installed (this is now the second occurrence; it is a
per-run cost, since the VM is rebuilt each session). Work done in the
session's own scratch clone outside the connected folder, per the
2026-09-17 finding that git cannot run inside it.

**Commits this run:** 3 (signed-carrier notes; the model + figure;
Chapter 6 Section 6.8 + 6.7 annotations + Chapter 1 status table). This
AUTOMATION_LOG.md entry makes 4.

---

## 2026-09-19 — Non-uniform illumination: a bound instead of a ranking, and two of this session's own claims falsified

**Status:** Full session. Closes the item that has sat at the top of
"not yet covered" since 2026-09-18 — the generation weight g(x) — which
has also been Simplification 2 of every two-contact module since
2026-09-07.

**Work done:**

1. **Research notes**
   (`notes/2026-09-19-non-uniform-illumination-and-the-collection-kernel.md`):
   literature basis, plus a **citation correction**. This repo had cited
   the shadow-mask experiment in six places as *"Suzuki et al., Carbon
   Trends 5, 100115 (2021)"*. **Both the first author and the article
   number were wrong.** The paper is Shimomura, Imai, Nakagawa, Kawai,
   Hashimoto, Ideguchi & Maki, *Carbon Trends* **5**, 100100 (2021),
   <https://www.sciencedirect.com/science/article/pii/S2667056921000778>,
   confirmed independently against the publisher page and the Ideguchi
   group's own publication list. The error entered on 2026-09-17 and was
   copied forward twice without being rechecked. The physics attributed
   to it is correct and nothing downstream changes. All six occurrences
   are fixed; the historical entries in *this file* are deliberately
   left alone, since they are a dated record.

   The substantive new number from the paper: **the mask leaks.** They
   report the photovoltage under the 50 nm Ni mask as "about half of the
   opposite side", i.e. T ≈ 0.5. This repo had used the paper for three
   sessions without ever using that number.

2. **Model** (`graphene_photodetector_nonuniform_illumination_model.py`,
   `photodetector_nonuniform_illumination.png`). The structural finding
   is that **illumination never enters transport at all**. Section 6.8's
   response is already an integral of a *collection kernel*
   k(x) = charge to contact A per photon absorbed at x, so

       N[g] = (1/L) ∫ g(x) k(x) dx ,   (1/L) ∫ g dx = 1

   with the normalisation fixing total absorbed photons. Three
   consequences, all exact: g ≡ 1 must return the 2026-09-18 number
   *bitwise*; **k(x) IS the delta-spot photocurrent scan**; and
   **|N[g]| ≤ max|k| for every pattern whatsoever** — the first quantity
   in this repo that bounds an entire design space instead of ranking
   points inside it.

   For a symmetric pair at zero bias, antisymmetry of k gives the leaky
   mask in closed form, **N(T)/N(0) = (1 − T)/(1 + T)**, so the measured
   T ≈ 0.5 leaves **one third**, not one half.

3. **Five validations, all against exactly known values:**
   - uniform g ≡ 1 reduces to the 2026-09-18 signed model over 196
     (pair, crossover, bias) combinations — **196/196 bitwise**,
     deviation 0.000e+00
   - mirror-symmetric g (uniform, centred spot, both-edges) on a
     symmetric pair at zero bias → |N| ≤ **1.3e−16**
   - mask-left = −(mask-right) → |dev| ≤ **2.2e−16**
   - the (1 − T)/(1 + T) closed form → |dev| ≤ **2.2e−16**. This is the
     one that tests physics rather than plumbing: dropping the
     equal-photon normalisation would give (1 − T), a 50% error at the
     experimentally relevant T = 0.5, and none of the other checks would
     have caught it
   - |N[g]| ≤ max|k| over 343 (pair, pattern) cases → **0 violations**

4. **RESULT — a mask rescues a symmetric device, and it is still the
   wrong strategy.** Per absorbed photon, symmetric Ti/Ti with a perfect
   mask reaches 0.841 against 1.832 for Ti/Pt under uniform light: 45.9%.
   **That comparison was this session's first headline and it was
   wrong.** It compares a masked device and an unmasked device *per
   absorbed photon*, which is the right comparison of collection
   mechanisms and the wrong comparison of detectors: responsivity is amps
   per incident watt, and a mask puts half the incident light into 50 nm
   of nickel. The two accountings differ by exactly (1 + T)/2 —
   (1 − T)/(1 + T) per absorbed photon versus **(1 − T) per incident
   photon** — so a perfect mask gives exactly half as much as the first
   figure suggests. Corrected: **22.9%** for a perfect mask and
   **11.5%** for the mask actually built. **Asymmetric metallisation
   beats illumination engineering by roughly 4×, and ~8× against a real
   mask.** Opposite emphasis from the draft number, and the one to quote.
   Corollary worth carrying: Shimomura *et al.*'s *other* design — a
   comb-shaped counter-electrode that enlarges one interface rather than
   shading the other — has no incident-photon penalty and is the more
   promising of their two ideas on this accounting.

5. **RESULT — a pre-registered prediction of this session's own note,
   falsified.** The note predicted, before the model was run, that
   masking "should gain little or nothing" on an already-asymmetric pair
   because its kernel does not change sign. It gains **+0.7% for Ti/Pd**
   and +0.03% for Ti/Pt. Sign was the wrong criterion: a single-signed
   kernel is still not flat, and at equal absorbed photons a mask
   *redistributes* rather than discards. The replacement is the computed
   bound **max|k| / |N_uniform| = 1.0025 (Ti/Pt), 1.0159 (Ti/Pd), < 1.02
   for every asymmetric pair**. The intended conclusion survives with a
   number attached instead of a bad argument: masks are for symmetric
   devices only.

6. **RESULT — the kernel reproduces the reported scan shape.** Symmetric
   pair: k runs +1.000 → 0 → −1.000, i.e. **opposite polarity at the two
   interfaces**, which is Shimomura *et al.*'s "polarities … are
   opposite" and Mueller *et al.*'s p–n–p. Ti/Pt: flat, between −1.837
   and −1.830, never changing sign — the signature distinguishing an n/p
   pair from a symmetric one. Note this covers the long-open
   "Park/Ahn SPCM sign map" item's *physics* from an open-access source,
   even though PubMed itself remains unreachable.

7. **Writing** (`thesis_draft/06-graphene-photodetectors.md`,
   `thesis_draft/01-introduction.md`): new Section 6.9 with the full
   treatment (old 6.9 → 6.10), stating both corrected claims against
   itself rather than presenting only the surviving version. Chapter 1's
   status table updated, and Chapter 7 handed a third thread.

**Honest caveats recorded, not worked around.**
(a) This is the **photovoltaic contribution** to a scanning-photocurrent
trace, not a prediction of a measured one: Kasırga's 2025 review
(arXiv:2509.09390) stresses that measured SPCM signals are frequently
photo-thermal, and this repo has no PTE or bolometric machinery.
(b) **Optical localisation is impossible at this geometry.** L = 200 nm;
Shimomura *et al.*'s spot is ~2 µm — ten times the whole channel. The
spot-size study is therefore a statement about required channel length
(90% of the delta limit needs σ ≲ 0.05 L, so a ~1 µm spot wants a ~20 µm
channel), not a proposal for this device. The only realisable non-uniform
illumination at 200 nm is a lithographic mask, which is what Shimomura
*et al.* built.

**Methodological note.** Three consecutive sessions have now had a
headline overturned — twice by the *next* session, and today for the
first time **within the same session**, by writing the prediction down
before running the model and by checking the accounting convention
against what responsivity actually means. Writing the falsifiable
version into the note first is cheap and worked; it is worth keeping.
The second failure (per-absorbed vs per-incident photon) was not a
physics error at all but a **units-of-comparison** error, and it would
have survived every one of the five exact validations, because all five
test the model against itself. Exact validations do not protect against
comparing the wrong two quantities.

**Not yet covered (candidates for future runs):**
- **A non-linear dW → doping-profile relation** — now the top open item
  in Chapter 6. Both the 2026-09-17 and 2026-09-18 models, and today's
  by inheritance, take the profile magnitude linear in dW with one λ for
  every metal. Giovannetti *et al.* find it not linear at all for the
  chemisorbed metals (Ti, Ni, Pd) — exactly the metals carrying the
  largest |dW| under the physical crossover. Would test whether Ti/Pt's
  1.832, and today's max|k| ceiling with it, are robust
- **Whether Chapter 4's contact-resistance results should be re-run at
  the 5.4 eV crossover.** Section 4.5's doping profile still uses the
  `|W − W_graphene|` magnitude convention and has never been audited for
  it. Open since 2026-09-18 and **may affect published Chapter 4 numbers**
- **Chapter 7 (discussion/outlook)** — now with three threads, the
  newest and sharpest being: Chapter 6 has a ceiling (max|k|), Chapters 4
  and 5 argue from unbounded single-junction optimisation; does an
  analogous ceiling exist for contact resistance and for interconnect
  resistivity?
- **Chapters 2–3 remain undrafted** despite their computational results
  being complete — still the largest remaining block of pure writing
- Reconciling the 0.12 eV measured potential step (Mueller *et al.*,
  arXiv:0902.1479) with the 0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS`
  assumes — open since 2026-09-18
- **A photo-thermoelectric term**, newly motivated: the Kasırga review
  makes the absence of one the main obstacle to comparing any of this
  chapter's position-resolved predictions with a measurement
- **Shimomura et al.'s comb-electrode design** (unequal contact
  *perimeter* rather than unequal metal or unequal illumination) — a new
  open item created today, and on the incident-photon accounting the
  more promising of their two geometries
- Integrating Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery. Today's max|k| ceiling
  now makes this *quantitative*: plasmonic patterning is an illumination
  pattern, so it is bounded by max|k| too
- Isolating the root cause of the Section 4.7 negative-residual result
  (TLM double-counting vs. `lambda_decay` mismatch) — open since
  2026-08-31
- Ti and Cr per-metal Rc recalibration (blocked by ResearchGate
  rate-limiting; not reattempted)
- Second independent edge-contact dataset (Lee *et al.* 2022, Wiley
  403'd 2026-09-05) — still open
- Small-hole-diameter (50–100 nm) upturn in Passi *et al.*'s
  patterned-contact data (Section 4.8.1) — still open
- Graphene-all-around-metal (liner/cap) interconnect model (Chapter 5,
  Section 5.4) and specific liner-material resistivity values

**Web search availability:** WebSearch available and used — 4 queries,
4 fetches. 3 of 4 fetches returned usable content; `arxiv.org/pdf/2509.09390`
returned no machine-readable text (its abstract page did), so the Kasırga
review is cited from its abstract only and no number is attributed to its
body. PubMed/PMC **not attempted**: four consecutive reCAPTCHA blocks
(2026-09-05, 09-07, 09-17, 09-18) make it a standing environment
limitation. The Park/Ahn SPCM sign map remains unreachable there, but its
physics is now covered from Mueller *et al.* 2009 on arXiv instead.

**Automation health:** Device reachable, folder connected. `scipy` again
absent from the device VM and pip-installed — **third consecutive
occurrence**, and it is a per-run cost because the VM is rebuilt each
session; `requirements.txt` already lists it, so this is an environment
fact rather than a repo defect. Work done in the session's own scratch
clone outside the connected folder, per the 2026-09-17 finding that git
cannot run inside it.

**Commits this run:** 4 (notes + citation correction; the model + figure;
the repo-wide citation fix; Chapter 6 Section 6.9; Chapter 1 status
table). This AUTOMATION_LOG.md entry makes 6.

---

## 2026-09-20 — The linear dW → doping assumption removed, and it took two 2026-09-19 claims with it

**Open item closed:** "A non-linear dW → doping-profile relation", the item
flagged top of the Chapter 6 list on 2026-09-19.

1. **Research** (`notes/2026-09-20-nonlinear-work-function-to-doping-relation.md`).
   Khomyakov, Giovannetti, Rusu, Brocks, van den Brink & Kelly, *Phys. Rev.
   B* **79**, 195425 (2009), [arXiv:0902.1203](https://arxiv.org/abs/0902.1203)
   — the companion paper to the PRL this repo has been citing for the 5.4 eV
   crossover. Its Eq. 7 gives the Fermi-level shift as
   `dE_F = sgn(dW')(sqrt(1 + 2*alpha*|dW'|) − 1)/alpha`, **linear only as
   dW' → 0**: charge enters graphene's linear DOS, so `n ~ E_F^2` and the
   interface potential step is quadratic in the shift. The paper also shows
   the 5.4 eV crossover this repo hard-codes is `W_0(d) = W_G + D_c(d)`
   *evaluated at the physisorbed separation*, not a constant of nature.

2. **Code** (`graphene_contact_doping_nonlinear_model.py`, 5 exact
   validations, all passing). `alpha = 2*e^3*(d−d0)/(eps0*pi*hbar^2*v_F^2)
   = 2.393 eV^-1` is **derived, not fitted** — deliberately, because the
   Table I extraction failed (below). Validations: `fermi_shift(0)` bitwise
   zero; bitwise odd symmetry; Eq. 1 ∘ Eq. 3 = identity to 8.9e-16;
   **alpha = 0 reproducing `signed_carrier_model.total_field()` BITWISE** on
   the asymmetric Ti/Pt pair with the residual confirmed cubic; and the
   symmetric-pair zero preserved at 1.3e-16.

3. **RESULT — Section 6.8's headline survives.** Ti/Pt goes **1.832 →
   1.742, a 4.9% compression**, and the alpha-sweep holds the ratio between
   0.92 and 1.00 over alpha = 0–5 eV^-1. Cr/Pt and Cu/Pt likewise ~5%. This
   was the point of the exercise and it is the reassuring half.

4. **RESULT — the session's own pre-registered prediction, half falsified.**
   The note predicted "compression, not reversal" everywhere, from
   monotonicity. Correct for pairs *straddling* the crossover; wrong for
   *same-sign* pairs, where **Ti/Pd loses 56%** (1.599 → 0.697) and Ti/Au
   57%. Monotonicity governs *signs*; a same-sign pair is a
   near-cancellation whose surviving magnitude is set by the *ratio* of the
   two offsets, and a concave map compresses ratios. **Compressing a ratio
   amplifies a cancellation's fragility** — the distinction the prediction
   missed.

5. **RESULT — RETRACTION of a 2026-09-19 bound, under that session's own
   linear model.** Ti/Pd's ceiling `max|k|/|N_uniform|` moving 1.016 → 1.434
   was implausible enough to prompt enumerating all 21 asymmetric pairs.
   2026-09-19 claimed "**< 1.02 for every asymmetric pair**"; **14 of 21
   violate it, worst Au/Pd at 22.4x**, with the linear model. Every
   straddling pair does satisfy it and every violator is same-sign, because
   the ratio diverges as `N_uniform → 0`. The *inequality* `|N[g]| ≤ max|k|`
   that session validated is a normalisation identity and remains true
   (0 violations, 343 cases); what was over-generalised, from a sample of
   two, is its **tightness**.

6. **RESULT — and therefore the 2026-09-19 conclusion is false too.**
   "Masks are for symmetric devices only": a perfect shadow mask gains
   **8.42x on Au/Pd**, 3.40x on Ti/Cr, and >2x on five pairs, scored per
   *incident* photon (2026-09-19's own accounting). The practical
   recommendation survives on **different grounds**, now stated explicitly:
   best masked device 0.377 per incident photon vs unmasked Ti/Pt's 1.832,
   so **a large relative gain on a small number is still a small number**.
   2026-09-19 conflated relative gain with absolute performance — the same
   class of error as its own per-absorbed vs per-incident correction, one
   level up.

7. **RESULT — four of seven metals are outside the model's regime, and the
   worst-placed one is Ti.** Khomyakov *et al.* put the chemisorbed metals
   at `d_eq` ≈ 2.05–2.3 Å, **below** `d0` = 2.4 Å, so Eq. 7's
   gap-capacitance term does not apply to **Ti, Ni or Pd**; **Cr** has no
   tabulated separation and was not guessed. `alpha_for_metal()` returns
   `None` with a reason rather than extrapolating. This matters because
   **Ti carries the largest |dW| in the table (−1.07 eV) and is contact A of
   both Chapter 6 headline pairs** — the chapter's central numbers rest on
   the one assumption the literature most explicitly disowns.

8. **Writing** (`thesis_draft/06-graphene-photodetectors.md`,
   `thesis_draft/01-introduction.md`). New Section 6.11 (three subsections,
   both result tables, and a five-item "what this does not claim").
   Section 6.9.5's numbers are **left in place** with a retraction block
   above them, per the repo's annotate-don't-delete rule. Chapter 1's status
   table updated and Chapter 7 handed a fourth thread.

**Honest reporting of a data-extraction failure.** Two WebFetch passes over
the *same* PRB PDF returned **mutually contradictory Table I contents** —
different signs, different metal-to-value assignment, different work
functions. At least one is a PDF-to-text artefact. **No per-metal `dE_F`
from that table is used anywhere.** Only quantities that agreed across both
passes are used (`d0` = 2.4 Å, `D_c(3.3 Å)` ≈ 0.9 eV, `W_0` = 5.4 eV, the
physisorbed/chemisorbed separations), and `alpha` is derived from first
principles instead — which is precisely why the derivation was worth doing
rather than a convenience. Cross-check against the paper's Pt figure:
**1.4x, i.e. right size and sign, not quantitative**, so every conclusion is
also reported as an alpha-sweep.

**Methodological note.** All five exact validations passed while the
retracted claim stood, because every one of them tests the model against
itself. **Exact validation protects against implementation error; it does
not protect against a claim quantified on an unrepresentative sample.** That
is now the second failure of this kind in four days (2026-09-19's was
comparing the wrong two quantities). The practice that caught today's was
neither a validation nor a pre-registered prediction but **enumerating the
whole space instead of tabulating two representative cases** — cheap here
(21 pairs), and worth making the default whenever a claim says "every".

**Not yet covered (candidates for future runs):**
- **A per-metal crossover `w_cross = W_G + D_c(d)`** — created today and now
  the **top** open item. `W_CROSS_CHEM = 5.4 eV` is held for all seven
  metals, but Khomyakov *et al.* give it as a function of separation. Unlike
  today's magnitude nonlinearity, this would **move every `dW` in the
  table** rather than compressing them, on the same evidence
- **A description of Ti, Ni and Pd that does not go through work function**
  — created today. Four of seven metals are outside Eq. 7's regime and Ti is
  contact A of both headline pairs, so this is now Chapter 6's central
  structural weakness rather than a refinement
- **Re-check whether other "for every ..." claims in the repo rest on small
  samples** — created today, directly from the retraction. Chapter 4's
  per-metal Rc recalibration and Chapter 5's liner scenarios both state
  general rules from tabulated subsets and have never been enumerated
- **Whether Chapter 4's contact-resistance results should be re-run at the
  5.4 eV crossover.** Section 4.5's doping profile still uses the
  `|W − W_graphene|` magnitude convention and has never been audited for it.
  Open since 2026-09-18 and **may affect published Chapter 4 numbers**
- **Chapter 7 (discussion/outlook)** — now with four threads, the newest
  being methodological: Chapter 6's claims have repeatedly failed on widened
  samples rather than on physics, and Chapters 4 and 5 have never been
  widened
- **Chapters 2–3 remain undrafted** despite their computational results
  being complete — still the largest remaining block of pure writing
- Reconciling the 0.12 eV measured potential step (Mueller *et al.*,
  arXiv:0902.1479) with the 0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS`
  assumes — open since 2026-09-18, and today's Δ_c finding is relevant to it
- **A photo-thermoelectric term** — the Kasırga review makes its absence the
  main obstacle to comparing any position-resolved prediction with a
  measurement
- **Shimomura et al.'s comb-electrode design** (unequal contact *perimeter*)
  — open since 2026-09-19, and today's retraction makes it more interesting,
  not less: it is another way to make `N_uniform` non-zero
- Integrating Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery
- Isolating the root cause of the Section 4.7 negative-residual result —
  open since 2026-08-31
- Ti and Cr per-metal Rc recalibration (blocked by ResearchGate rate-limiting)
- Second independent edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)
- Small-hole-diameter (50–100 nm) upturn in Passi *et al.*'s data (§4.8.1)
- Graphene-all-around-metal (liner/cap) interconnect model (Chapter 5, §5.4)

**Web search availability:** WebSearch available and used — 2 queries,
3 fetches. 2 of 3 fetches returned usable content; `arxiv.org/abs/0902.1590`
was fetched on a mis-guessed identifier and returned an unrelated paper
(recorded rather than silently discarded), and the IFW Dresden PRB mirror
returned usable equations but **contradicted itself between two passes on
Table I**. PubMed/PMC **not attempted**: five consecutive reCAPTCHA blocks
(2026-09-05, 09-07, 09-17, 09-18) make it a standing environment limitation,
and both papers are open access on arXiv anyway.

**Automation health:** Device reachable, folder connected. `scipy` again
absent from the device VM and pip-installed — **fourth consecutive
occurrence**; `requirements.txt` already lists it, so this is an environment
fact (the VM is rebuilt each session), not a repo defect. Work done in the
session's own scratch clone outside the connected folder, per the 2026-09-17
finding that git cannot run inside it.

**Commits this run:** 5 (the note; the model + figure; Chapter 6 Section
6.11 with the 6.9.5 retraction in place; Chapter 1 + Chapter 7; the note's
outcome section). This AUTOMATION_LOG.md entry makes 6.

## 2026-09-21 — The single 5.4 eV crossover removed, and a pre-registered prediction falsified at 17%

**Open item closed:** "A per-metal crossover `w_cross = W_G + D_c(d)`" — the
**top** item on 2026-09-20's list.

**Procedural change made this session.** The four predictions were written
into `notes/2026-09-21-...md` Sections 1–4 and **committed to git before the
model existed** (commit `588bb8d`, model in `a856bee`). The history is the
evidence that they were not written to fit the answer. 2026-09-20 observed
that this repo's claims have twice failed on widened samples rather than on
physics; pre-registration is the cheap half of the fix and enumeration is the
other, and both were applied today.

1. **Research** (`notes/2026-09-21-per-metal-crossover-from-the-chemical-interface-term.md`).
   Khomyakov *et al.*, *Phys. Rev. B* **79**, 195425 (2009) write the
   interface step as a charge-transfer term plus a **short-range chemical**
   term `Δ_c(d) = e^(−γd)(a₀ + a₁d + a₂d²)` (their Eq. 4). Neutrality puts
   the crossover at `w_cross(d) = W_G + Δ_c(d)`, so the **5.4 eV this repo
   has hard-coded since 2026-09-18 is that expression evaluated at the
   physisorbed separation** — not a constant shared by seven metals.

2. **Code** (`graphene_per_metal_crossover_model.py`, 4 exact validations).
   Eq. 4's fitted constants could **not** be extracted — two ar5iv fetches
   with differently-worded prompts both returned the equation symbolically
   with the numbers absent, the **second** extraction failure on this same
   paper (2026-09-20's Table I was contradictory rather than missing). They
   are not guessed. `Δ_c` is instead a one-parameter family pinned to the one
   value every pass agrees on, `Δ_c(3.3 Å) = 0.9 eV`, with the decay length
   `ℓ` swept over 0.3–1.5 Å.

   Validations: the new ΔW-level entry point equals the 2026-09-18 signed
   model on **49/49 ordered pairs bitwise**; `ℓ → ∞` collapses onto the flat
   convention with `Δ_c(d,∞) == 0.9` and `w_cross == W_G + 0.9` bitwise for
   every metal and **49/49 pair responses bitwise**; the anchor is bitwise
   for **61/61** values of `ℓ`; and the symmetric-pair zero (6.6e-17) and
   charge conjugation (**exactly** 0.000e+00, 27 cases) both survive.

3. **RESULT — the pre-registered ≤10% prediction is FALSIFIED at 17.12%.**
   Cu's `dW` goes −0.750 → −0.878 at the short end of the `ℓ` range;
   Au 9.84%; Pt 0.00%. That is **more than three times** 2026-09-20's 4.9%
   compression, which is the precedent the 10% was calibrated on — and it was
   the wrong precedent. **Moving the zero of `dW` is a larger perturbation
   than compressing `dW` about a fixed zero**, because the shift does not
   scale with `|dW|`: a metal *near* the crossover feels a fixed shift most,
   the opposite of the usual intuition. Bisection puts the boundary at
   **ℓ = 0.50 Å** — the prediction holds above it and fails below, i.e. over
   roughly the shortest sixth of the range. Declaring that range generously
   is what made the failure visible.

4. **RESULT — Pt's exact 0.00% is an artefact, reported as one.**
   `d_eq(Pt) = 3.30 Å` coincides with the anchor, so Pt is pinned by
   construction for every `ℓ`. **Every per-metal number this session produced
   is a shift relative to Pt, not an absolute one.**

5. **RESULT — a 65-fold difference between same-sign and straddling pairs.**
   All six admitted ordered pairs enumerated: Cu/Au moves **77.94%** while
   Cu/Pt moves **1.21%** and Au/Pt **1.06%**, from a 17.12% input shift —
   amplification 4.55 vs 0.07–0.11. 2026-09-20 inferred this
   ratio-amplification from **one** same-sign pair under a perturbation that
   *compressed* offsets; seeing it again, same sign and comparable size,
   under a perturbation that *moves their zero* is independent evidence that
   it belongs to near-cancelling pairs rather than to either model. **P3
   held; this is the session's most transferable result.**

6. **RESULT — the anchored exponential diverges on the chemisorbed metals,
   for a reason unrelated to 2026-09-20's.** Extrapolated inward, (N1) gives
   `w_cross` of **6.25–62.6 eV** for Ni/Ti/Pd — every value above 5.9 eV, the
   highest elemental work function anywhere — implying `|dW|` of 1.13–57.5 eV
   against a largest *observed* `dE_F` of 0.5 eV, **across the whole `ℓ`
   range**, not just its short end. This is a failure of the simplification
   (Eq. 4's polynomial prefactor, dropped because its coefficients were
   unavailable, is exactly what lets `Δ_c` turn over), and it is worth
   committing because it reaches 2026-09-20's conclusion **by an independent
   route**: the gap-capacitance term refuses Ti/Ni/Pd for `d_eq < d₀`, the
   chemical term refuses the same three for divergence. **Two distinct parts
   of one framework failing on the same three metals for unrelated reasons**
   is a structural claim about which contacts can be described by a work
   function at all — and Ti is contact A of both Chapter 6 headline pairs.

7. **Net effect on Chapter 6.** The Cu/Pt straddling rule is robust to 1.2%
   and the chapter's qualitative design conclusion survives cleanly. The
   Ti/Pt headline is **untested** by two successive refinements rather than
   confirmed by them — a weaker position than 2026-09-20 left it in.
   Sections 6.8–6.11 are **not retracted**; they now carry an explicit
   **±17% (per-metal) / ±78% (same-sign pair)** error bar, annotated in place
   at 6.8.7.

8. **Writing** (`thesis_draft/06-graphene-photodetectors.md`,
   `thesis_draft/01-introduction.md`). New Section 6.12 (seven subsections,
   three result tables, a does-and-does-not-change list and a
   does-not-claim list); Section 6.8.7 annotated in place; Chapter 1 status
   table updated; Chapter 7 handed a **fifth** thread.

**Methodological note, continuing 2026-09-20's.** That session concluded that
exact validation protects against implementation error but not against a
claim quantified on an unrepresentative sample. Today adds the complement:
**all four validations passed and the falsified prediction was still
falsified**, because the prediction was about the *size* of a physical
effect, which no self-consistency check can reach. What caught it was
declaring a generous range and a number in advance. The Chapter 7 framing is
now: a claim's reliability tracks **how wide a sample it was checked against
and whether its error bar was stated in advance**, not how carefully its
physics was derived.

**Not yet covered (candidates for future runs):**
- **Propagate the per-metal crossover back through Sections 6.9–6.11** —
  created today and now the **top** open item. Those sections' numbers are
  still at the flat 5.4 eV convention, and today's Cu/Au result says a
  same-sign pair under the flat convention can be wrong by 78%. The 2026-09-20
  retraction table (21 asymmetric pairs) is the obvious first target, since
  its worst violators were same-sign pairs
- **A description of Ti, Ni and Pd that does not go through work function** —
  now supported by **two independent failures** rather than one, and still
  Chapter 6's central structural weakness
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — created
  today. One anchor makes every result relative to whichever metal sits at
  it; a second would make the shifts absolute and would also constrain `ℓ`,
  which is currently only swept
- **Re-check whether other "for every …" claims in the repo rest on small
  samples** — open since 2026-09-20; Chapter 4's per-metal Rc recalibration
  and Chapter 5's liner scenarios still state general rules from subsets
- **Ask whether Chapter 4's `Rc` and Chapter 5's resistivity are
  near-cancellations** — created today, and the sharpest form the synthesis
  question has taken: today measured a 65x amplification of input uncertainty
  for a near-cancelling quantity versus a reinforcing one, and neither of
  those chapters has been asked which kind it is
- **Whether Chapter 4's contact-resistance results should be re-run at the
  5.4 eV crossover** — open since 2026-09-18, and today makes it sharper
  still: Section 4.5 uses the `|W − W_graphene|` magnitude convention, which
  is neither the flat crossover nor a per-metal one
- **Chapter 7 (discussion/outlook)** — now with five threads and, as of
  today, a positive methodological claim to make rather than only a
  self-critical one
- **Chapters 2–3 remain undrafted** despite their computational results being
  complete — still the largest remaining block of pure writing
- Reconciling the 0.12 eV measured potential step (Mueller *et al.*,
  arXiv:0902.1479) with the 0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS`
  assumes — open since 2026-09-18
- **A photo-thermoelectric term** — the Kasırga review makes its absence the
  main obstacle to comparing any position-resolved prediction with a
  measurement
- **Shimomura et al.'s comb-electrode design** (unequal contact *perimeter*)
  — open since 2026-09-19
- Integrating Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery
- Isolating the root cause of the Section 4.7 negative-residual result —
  open since 2026-08-31
- Ti and Cr per-metal Rc recalibration (blocked by ResearchGate rate-limiting)
- Second independent edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

## 2026-09-22 — Chapters 4 and 5 conditioned: the worst-cancelling step in the thesis is in the least-examined chapter

**Open item closed:** "Ask whether Chapter 4's `Rc` and Chapter 5's
resistivity are near-cancellations" — created 2026-09-21 and described there
as *the sharpest form the synthesis question has taken*. One branch of
"Isolating the root cause of the Section 4.7 negative-residual result"
(open since **2026-08-31**, the repo's oldest open item) is also closed.

**Procedure.** Six predictions written into
`notes/2026-09-22-condition-number-audit-of-chapters-4-and-5.md` and
**committed before the model existed** (note `968a9c0`, model `24c5fab`),
continuing 2026-09-21's practice. New this session: each prediction is
labelled `[blind]` or `[hand-derivable]`, and for the hand-derivable ones the
hand estimate is written into the note, so that a machine result disagreeing
with it exposes a bug rather than being silently accepted. **Four predictions
held; two failed.**

1. **Diagnostic** (`graphene_sensitivity_audit.py`). For `Q = Σ_i T_i`, the
   signed logarithmic sensitivity `S_i = T_i/Q` and the condition number
   `κ = max|S_i|`. Classification: `κ ≤ 1` reinforcing, `1 < κ < 3` mildly
   ill-conditioned, `κ ≥ 3` near-cancelling (Chapter 6's same-sign pairs sit
   at 4.55).

   **Four exact validations.** (V1) Euler's sum rule `Σ S_i = 1` held to
   **1.8e-15** over **15 decompositions** across three chapters — an
   identity a merely-approximate derivative cannot satisfy. (V2) known
   power-law exponents recovered to **1.7e-11**, including two that are
   **exactly** `+1` and `0.0`. (V3) degenerate limits: `Q = A − A` gives
   `1/κ = 0.000e+00` exactly, `Q = A + A` gives `κ = 0.5` exactly; scale
   invariance is **3/5 bitwise, not 5/5**, worst relative deviation 6.2e-16,
   **reported as such** rather than rounded to a pass. (V4, *not designed
   in* — it fell out of Family C) at `W = W_calibration` the model
   reproduces the calibration datum **bitwise** and **three** sensitivities
   are **exactly 0.0**.

2. **RESULT — the worst-conditioned step in this thesis is Chapter 5's
   calibration, not anything in Chapter 6.** `λ_impurity` is solved from
   `3.0000 − 1.0000 − 1.8478 = 0.1522`: `κ = 19.71`, against Chapter 4's
   11.46 and Chapter 6's 4.55. **A 5% error in the single datum
   `ρ = 3.6 µΩ·cm at 22 nm` moves the extracted impurity mean free path by
   99%.** It was invisible in the chapter's own formula, which is a
   manifestly reinforcing Matthiessen sum (`κ ≤ 0.66` at every width — P2's
   first half, confirmed) and propagates into `ρ(W)` damped but not removed
   (`κ` 0.88 → **1.551** over 18–52 nm; P2's second half predicted "above 1,
   below 3" and held). **A quantity can be well-conditioned in its own
   arguments and ill-conditioned in the quantities those arguments were
   derived from** — auditing the visible formula alone returns a clean bill
   of health.

3. **RESULT — the calibration point pins the model, and `S(ρ_bulk)` changes
   sign through it.** At 22 nm, `S(ρ_bulk) = S(p) = S(λ_bulk) = 0.0`
   exactly: the prediction there carries **no information** from the physical
   inputs, only from the calibration datum. `S(ρ_bulk)` runs `+0.120` (18 nm)
   → `0.000` (22 nm) → `−0.551` (52 nm). **The further from 22 nm a Chapter 5
   prediction is made, the more of its content comes from one measurement** —
   the opposite of how a calibrated model is usually read.

4. **RESULT — `κ` and sign-robustness run in OPPOSITE order in Section 4.7,
   and the negative-residual conclusion is *strengthened*.** The useful
   question is not how an error is amplified but how wrong `R_extra` would
   have to be to flip the residual's **sign**:

   | metal | residual | `κ` | `R_extra` must move by | sign |
   |---|---|---|---|---|
   | Pd | **+50.9** | **11.46** | **−9.6 %** | **fragile** |
   | Au | −90.2 | 6.76 | +14.8 % | intermediate |
   | Ni | −360.5 | 2.30 | +43.4 % | intermediate |
   | Cu | −787.5 | 1.23 | **+81.1 %** | **robust** |

   Pd — the only positive residual, and the one Section 4.7 called "a small,
   plausible positive residual", i.e. the single row consistent with the
   additive decomposition surviving — is the **weakest** row in the table.
   Cu, the largest apparent violation, is the **strongest**. So (i) the
   additive decomposition `R_c = R_extra + R_transmission` fails **robustly**
   rather than marginally, and (ii) **amplified input error is eliminated as
   the cause of the negatives**: had they been an amplification artefact, the
   largest negatives would sit at the largest `κ`; they sit at the smallest.
   Section 4.7's own two hypotheses (TLM double-counting, `λ_decay`) stand;
   a third that had never been named is removed. **First narrowing of that
   item since 2026-08-31.**

5. **RESULT — the Cu liner model's `W_eff = W − 2t` is an unremarked
   subtraction**, `κ` 1.130 → **1.500** as `W` falls to 18 nm, with
   `|S_t(ρ_eff)| = 1.167` there. Chapter 5's own two sources for `t` differ
   by ~17% (3 nm vs a "2–3 nm floor"), which converts to a **~20% bar on
   `ρ_eff` at 18 nm** — previously unstated, and in exactly the width range
   where Chapter 5's graphene-vs-Cu crossover argument is made. This is the
   one input in either chapter whose uncertainty is literature-available
   rather than hypothesised.

6. **Cross-chapter check (P5, blind, PASS).** The same general machinery
   reproduces 2026-09-21's Chapter 6 amplifications — 4.55 for the same-sign
   pair Cu/Au, 0.07–0.11 for the straddling pairs — to within **1.9%**, so
   the diagnostic is measuring the quantity that session measured and not a
   differently-normalised cousin of it.

7. **The two failures, reported in full.** **P4 FAILED**: "at least three of
   four metals near-cancelling" came in at 2/4 (Ni 2.30 and Cu 1.23 are only
   *mildly* ill-conditioned). The hand arithmetic behind it was correct for
   the two metals it was done for and was generalised from those two —
   **the third consecutive session in which a claim generalised from a
   partial enumeration failed on the full one**, and the first in which the
   failure was in a *prediction about* the model rather than in the model.
   The 2026-09-20 enumeration rule was applied to the model this session but
   not to the prediction. **P6 FAILED**: no published Chapter 4 or 5 number
   carries a >100% bar under a 5% single-input perturbation; the worst is
   Pd's 57%. The reason is the opposite of what P6 assumed — `κ` is large
   exactly where the *output* is small, so large `κ` here threatens **signs,
   not orders of magnitude**, and the sign-flip margin of (4.23) is the right
   instrument. P6 asked the wrong question and the audit had to supply the
   right one.

8. **Writing.** Chapter 4 Section 4.10 (five subsections, three tables,
   Section 4.7 annotated in place with **no number changed**); Chapter 5
   Section 5.5 (five subsections, four tables); Chapter 1 status table for
   Chapters 4, 5 and 7; **Chapter 7 handed a sixth thread.**

**Methodological note, continuing 2026-09-20's and 2026-09-21's.** Those two
sessions concluded that a claim's reliability tracks how wide a sample it was
checked against and whether its error bar was stated in advance. Today adds a
third term and a correction. The third term: **which of a quantity's terms
cancel**, which is a property of arithmetic that no amount of care about
physics detects. The correction is sharper and reframes threads four and five
of Chapter 7. The chapters that have publicly overturned their own headlines
are Chapter 6's; the worst-conditioned arithmetic in the thesis is Chapter
5's, `κ = 19.71` against Chapter 6's 4.55, and it had gone four weeks
unnoticed. **The difference between the chapters is not rigour and not
fragility — it is how much scrutiny each has received.** Chapter 7's claim
should therefore be about unequal examination, which is a claim about method,
rather than about Chapter 6's unreliability, which is a claim about graphene
and is not supported.

**Not yet covered (candidates for future runs):**
- **Condition the rest of Chapter 4 and Chapter 5** — created today and the
  **top** open item, because today's audit covered only three families.
  Section 4.8's edge-vs-top geometry model contains a **ratio of two computed
  resistances** that has never been conditioned; Section 4.6's `f_max` is a
  square root of a difference; Section 5.3's parallel-conduction refinement
  and the `p_cu = 0.6` Fuchs–Sondheimer baseline are both unaudited. On
  today's evidence the unexamined places are where the large `κ` live
- **A second calibration point for Chapter 5** — created today and now
  concrete: `κ = 19.71` is a direct consequence of having exactly one, and
  Section 5.5.3's three exact zeros show what a one-point fit costs. A second
  datum at a different width would over-determine `λ_impurity` and turn the
  audit's error bar into a residual
- **Propagate the per-metal crossover back through Sections 6.9–6.11** —
  top item on 2026-09-21's list, **not** attempted today and now understood
  to be partly blocked: the per-metal model admits only Cu, Au and Pt, so of
  the 21 asymmetric pairs in the 2026-09-20 retraction table only 3 can be
  recomputed. The blocker is the same admissibility problem as the next item
- **A description of Ti, Ni and Pd that does not go through work function** —
  supported by two independent failures, and today's finding makes it more
  pressing rather than less: it is what gates the item above
- **Apply the sign-flip margin (4.23) to Chapter 6's own results** — created
  today. Chapter 6 has error bars in `κ` but has never asked which of its
  conclusions survive a sign flip, and Section 4.7 showed the two orderings
  can be exactly reversed
- **Re-check whether other "for every …" claims in the repo rest on small
  samples** — open since 2026-09-20 and **now with a third instance** (P4);
  the rule keeps being applied to models and not to the claims made about them
- **Whether Chapter 4's contact-resistance results should be re-run at the
  5.4 eV crossover** — open since 2026-09-18; Section 4.5 still uses the
  `|W − W_graphene|` magnitude convention
- **Chapters 2–3 remain undrafted** despite their computational results being
  complete — still the largest remaining block of pure writing, and now the
  only chapters with no conditioning statement at all
- **Chapter 7** — six threads, and as of today a positive methodological
  claim (unequal scrutiny) with quantitative support from three chapters
- Reconciling the 0.12 eV measured potential step (Mueller *et al.*,
  arXiv:0902.1479) with the 0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS`
  assumes — open since 2026-09-18
- **A photo-thermoelectric term** — the Kasırga review makes its absence the
  main obstacle to comparing a position-resolved prediction with a measurement
- **Shimomura et al.'s comb-electrode design** (unequal contact *perimeter*)
  — open since 2026-09-19
- **A second anchor for `Δ_c`** at any separation other than 3.3 Å — open
  since 2026-09-21
- Integrating Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery
- Ti and Cr per-metal Rc recalibration (blocked by ResearchGate rate-limiting)
- Second independent edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)
- Small-hole-diameter (50–100 nm) upturn in Passi *et al.*'s data (§4.8.1)
- Graphene-all-around-metal (liner/cap) interconnect model (Chapter 5, §5.4)

**Web search availability:** WebSearch/WebFetch **not used and not needed** —
the session's question was entirely internal to the existing model and its
already-cited inputs. Recorded explicitly rather than left ambiguous. PubMed/
PMC not attempted (five consecutive reCAPTCHA blocks make it a standing
environment limitation).

**Automation health:** Device reachable, folder connected, 10:30 UTC firing;
neither repo had a 2026-09-22 entry, so a full session was run. `scipy` again
absent from the device VM and pip-installed — **fifth consecutive
occurrence**; `requirements.txt` already lists it, so this is an environment
fact (the VM is rebuilt each session), not a repo defect. Work done in the
session's own scratch clone outside the connected folder, per the 2026-09-17
finding that git cannot run inside it.

**Commits this run:** 6 (the pre-registered note; the audit model + machine
output + figure; the note's outcome section; Chapter 4; Chapter 5; Chapter 1
+ Chapter 7). This AUTOMATION_LOG.md entry makes 7.

## 2026-09-22 (second session of the day) — An exact identity makes yesterday's headline statistic undefined, and an unbracketed bisection prints 21 plausible fake roots

**Two automated sessions fired on 2026-09-22 and both did full work.** The
entry above (Chapters 4 and 5 conditioned, `κ = 19.7`) is the first; this is
the second. The redundancy check that normally prevents this compares
`AUTOMATION_LOG.md` and today's commits at the *start* of a run, and both
runs passed it because the first had not yet pushed when the second cloned
(clone 10:49 UTC, first run's push 10:54 UTC). **They lost a race, not the
check.** The two sessions turned out to be complementary rather than
duplicated — the first closed 2026-09-21's *near-cancellation* item on
Chapters 4 and 5, this one closed its *top* item on Chapter 6, and they
share no file except this log and the Chapter 1 status table — but that is
luck, and the scheduling should be fixed rather than relied on.

**Open item closed:** "Propagate the per-metal crossover back through Sections
6.9–6.11" — the **top** item on 2026-09-21's list. Closed by establishing
first that it **cannot be done as asked**, then doing the propagation that
can.

**Why it cannot be done as asked.** Section 6.12's per-metal crossover needs
the metal–graphene separation `d_eq`, tabulated for three of this repo's
seven metals (Cu, Au, Pt). Ti/Ni/Pd are chemisorbed, where 6.12.4 showed the
anchored exponential diverges; Cr has no tabulated `d_eq`. **All five
headline pairs of the 2026-09-20 retraction table — Au/Pd, Ti/Cr, Ni/Au,
Cr/Cu, Ni/Pd — contain a metal outside that set.** The tool that prompted the
question cannot answer it for a single one of them.

1. **Research** (`notes/2026-09-22-how-far-can-the-crossover-move.md`).
   Fetched both source papers. Giovannetti *et al.*, PRL **101**, 026803
   (2008) and Khomyakov *et al.*, PRB **79**, 195425 (2009) *both* write
   *"a metal work function of **~5.4 eV**"* — a tilde, no error bar. This
   repo has read it as three significant figures since 2026-09-18. The PRB
   abstract also gives the chemisorption shift as *"reduces considerably"*
   with no number: the **third** distinct extraction failure on this pair of
   papers, after 2026-09-20's contradictory Table I and 2026-09-21's
   twice-missing Eq. 4 coefficients. Recorded, not worked around.

   The perturbation adopted is a scalar offset `δ` on `w_cross` — wrong shape
   (not per-metal), right reach (needs no `d_eq`), and algebraically
   identical to a common error in every work-function entry, which is live
   given Ni's own 4.9–5.35 eV comment and the unreconciled 0.12 eV Mueller
   step. Four predictions and one derivation pre-registered and **committed
   before the model existed** (note `235b420`, model `a939852`).

2. **RESULT — P1 falsified by an exact identity, and this is the session's
   main result.** For a pair straddling the crossover,
   `⟨|dW|⟩ = ((w_cross − W_A) + (W_B − w_cross))/2 = (W_B − W_A)/2`, in which
   `w_cross` cancels identically. Verified to `0.000e+00` over 126 cases and
   promoted to **Validation 5** — it was not planned; it was found while
   trying to score P1. So a scalar offset changes a straddling pair's
   response while changing the input measure by **exactly nothing**, and
   2026-09-21's amplification ratio does not merely blow up for those pairs,
   it is **undefined**. **That ratio is a property of the (pair,
   perturbation) couple, not of the pair**, and yesterday's log came close to
   quoting "65×" as a device property.

3. **RESULT — the dichotomy is real, at 30× not 65×.** Restated in
   `S = d ln|N|/dδ`, which does not divide by the input: straddling pairs
   `|S| ∈ [0.0077, 0.0218]` eV⁻¹, same-sign pairs `[0.6625, 2.5628]`, **no
   overlap**, ratio **30.4**. Genuinely independent confirmation — three
   pairs under a per-metal perturbation, now twenty-one under a scalar one —
   and the third time in this chapter that widening a sample has shrunk a
   claim. P2 (Au/Pd most sensitive) held.

4. **RESULT — no scalar offset can reverse a response sign, and this is the
   chapter's first proof-shaped robustness result.** With `dW_A = m − s`,
   `dW_B = m + s`, a scalar offset moves `m` and leaves `s = (W_B − W_A)/2`
   **exactly** fixed; `m` multiplies the antisymmetric part of `E(x)`, which
   contributes exactly zero to `N`, and `s` multiplies the symmetric part,
   which sets `sign(N)`. Scan over ±3 eV, 601 points × 21 pairs = **12621
   evaluations, zero sign changes**; smallest `|N|` anywhere 5.1e-03,
   approached and never crossed. The *direction* of a two-terminal
   photoresponse — which is what a design rule asserts — is immune to the
   whole tilde and to a scalar error six times its size.

5. **RESULT — D5 falsified, and the bug that nearly hid it. A THIRD FAILURE
   CLASS.** D5 held that a sign flip occurs when `w_cross` moves between a
   pair's work functions. The premise is wrong. Worse, the **first
   implementation reported 21 roots**: it bisected between `δ = 0` and the
   midpoint of D5's window **without checking that a root was bracketed**,
   and with no sign change such a loop walks its lower bound to its upper
   bound and returns the **endpoint**. Every root equalled
   `(W_A + W_B)/2 − 5.4`; every one looked physical; the headline was
   *"nearest sign flip: Pd/Pt at −0.0150 eV"*. **All five exact validations
   passed while it printed that**, because they validate the model and the
   fault was in the analysis layer. Caught by hand-checking one row against
   what the field does at `dW_A = −dW_B`. A bracket assertion is now the
   first line of that function, and the bug is documented in its docstring
   rather than quietly fixed.

6. **RESULT — Section 6.11.2's "five pairs" retracted as a count.** Re-scored
   across the band, the number of pairs gaining >2× from a perfect mask runs
   **3 → 4 → 5 → 5 → 6** over `δ ∈ [−0.2, +0.2]` eV. The table's numbers are
   correct at `δ = 0` and are **not** retracted; the *count* is, being a
   threshold evaluated at one point of a tilde'd parameter. Annotated in
   place at 6.11.2. The qualitative claim and the design recommendation hold
   at every offset. Section 6.9.2's ceiling `|N[g]| ≤ max|k|`: 0 violations
   in 210 checks.

7. **RESULT — Ti/Pt moves 0.45%** over the whole nominal band, against 4.9%
   (2026-09-20) and 1.2–1.8% (2026-09-21). Three successive refinements have
   failed to move it; `w_cross` is not where its risk lives. P3 held.

8. **A pre-registered number missed, and recorded.** Shift invariance was
   predicted to agree to ≤1e-15; measured worst case **1.776e-15**, with
   215/231 combinations bitwise. Wrong by 1.8×, and wrong in shape: the right
   bound is a few ulp of `N` (1.78e-15 is 8 ulp of 0.593), not an absolute
   constant. Logged because the practice is worth nothing if only the
   comfortable misses are.

9. **Writing** (`thesis_draft/06-graphene-photodetectors.md`,
   `thesis_draft/01-introduction.md`). New Section 6.13, ten subsections,
   four tables, a does/does-not-change list and a does-not-claim list;
   Section 6.11.2 annotated in place; Chapter 1 status table updated; Chapter
   7 handed a **sixth** thread.

**Methodological note, continuing 2026-09-20's and 2026-09-21's.** 09-20:
exact validation protects against implementation error, not against an
unrepresentative sample. 09-21: pre-registration reaches what exact
validation cannot, namely a claim about the *size* of an effect. Today
closes the loop uncomfortably: **pre-registration does not reach an error in
the analysis layer either.** Three failure classes are now on record with a
distinct detector each —

| class | example | caught by |
|---|---|---|
| implementation error in the model | `inf − inf` transit integral (09-17) | an exactly known zero |
| claim quantified on a small sample | the `<1.02` ceiling (6.9 → 6.11) | enumerating all 21 pairs |
| error in the analysis layer | today's unbracketed bisection | **neither** — hand-checking one row |

— and the useful conclusion is that verification effort should be allocated
**by failure class**, not uniformly. Chapters 4 and 5 have had class (i)
attention only.

**Not yet covered (candidates for future runs):**
- **The DIFFERENTIAL part of the crossover uncertainty** — created today and
  the new **top** item, and the sharpest form the propagation question has
  taken. Today proves a *scalar* offset cannot change `s = (W_B − W_A)/2` and
  therefore cannot change `sign(N)`. A *per-metal* offset can. **A crossover
  difference of 0.02 eV would flip Au/Pd, where a 3 eV scalar offset cannot.**
  Section 6.12 can supply that difference for exactly three metals, which is
  enough for Cu/Au, Cu/Pt and Au/Pt — one same-sign pair and two straddling
  ones, i.e. exactly the 2026-09-21 sample, now with a sign test attached
- **A second anchor for `Δ_c`**, at any separation other than 3.3 Å — open
  since 2026-09-21, and now the prerequisite for the item above
- **A description of Ti, Ni and Pd that does not go through work function** —
  two independent failures on record; still Chapter 6's central structural
  weakness, and today it blocked the propagation that was asked for
- **Audit the repo's other analysis-layer code for the same bug class** —
  created today. Every root-finder, optimiser and threshold-crossing search
  in this repo should be asked whether it verifies its own bracket.
  `alpha_from_separation`, the `ell` bisection of 2026-09-21 and
  `contact_resistance_crossover.py` are the obvious first three, and `graphene_sensitivity_audit.py` from the first session
  of today is now a fourth
- **Re-check whether other "for every …" claims rest on small samples** —
  open since 2026-09-20; Chapter 4's per-metal `Rc` recalibration and Chapter
  5's liner scenarios still state general rules from subsets
- ~~Ask whether Chapter 4's `Rc` and Chapter 5's resistivity are
  near-cancellations~~ — **closed by the first session of today**, which
  measured `κ = 19.7` on the Chapter 5 calibration solve. Noted here because
  this session independently reached the same statistic from the other end:
  6.13.3 shows the amplification *ratio* is undefined for straddling pairs,
  so `S = d ln Q / d ln x` is the right quantity — and that is exactly the
  `S_i` the conditioning audit used. The two sessions agree on the
  diagnostic without having coordinated on it, which is worth more than
  either result alone
- **Whether Chapter 4's contact-resistance results should be re-run at the
  5.4 eV crossover** — open since 2026-09-18; Section 4.5 still uses the
  `|W − W_graphene|` magnitude convention
- **Chapter 7 (discussion/outlook)** — six threads, and as of today a
  taxonomy rather than a list
- **Chapters 2–3 remain undrafted** despite complete computational results —
  still the largest remaining block of pure writing
- Reconciling Mueller *et al.*'s 0.12 eV step (arXiv:0902.1479) with the
  0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS` assumes — open since
  2026-09-18, and today makes it quantitative: a common error there is
  exactly `δ`, hence bounded by the `S` table
- A photo-thermoelectric term (Kasırga review)
- Shimomura *et al.*'s comb-electrode design — open since 2026-09-19
- Integrating Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery
- Isolating the root cause of the Section 4.7 negative-residual result —
  open since 2026-08-31; **one branch closed by the first session of today**
- Ti and Cr per-metal `Rc` recalibration (ResearchGate rate-limiting)
- Second independent edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

---

## 2026-09-23

**Status:** Automated session. Live web search was **not used**; the session
was analytical and the inputs were already in the repo. The top item of
2026-09-22's "Not yet covered" list, taken verbatim.

**The question.** 2026-09-22 proved that no *scalar* crossover offset can
reverse a two-terminal photoresponse: a common offset moves
`m = (dW_A + dW_B)/2` and leaves `s = (dW_B − dW_A)/2` exactly fixed, and `s`
sets the sign. That theorem covers the **common mode** of the crossover
uncertainty and nothing else. Section 6.12 predicts a crossover *per metal*,
so two contacts do not share one, and their **difference** is the only part
that moves `s`. 09-22's log named this the top open item and attached a
number to it: *"a crossover difference of 0.02 eV would flip Au/Pd, where a
3 eV scalar offset cannot."*

**Pre-registration**, fourth consecutive session, and this time with a change
of practice: each prediction is labelled **[D]** where the note's own
reasoning already constrains it (so only a *failure* is informative) or
**[P]** where it is genuinely open. 09-22's four predictions were not
labelled this way and two of them were closer to [D] than they looked.
Committed at `7175db7`, before `graphene_differential_crossover_model.py`
existed.

**Work done:**

1. **The threshold, in closed form, and it is the session's cleanest
   result.** Parameterise the two contact crossovers as `w_cross + c ∓ τ/2`,
   so that exactly `m → m − c` and `s → s − τ/2`. Then

       N = 0 exactly when s = 0, i.e. at   τ* = W_B − W_A

   independently of `c`, of `λ`, of `L` and of every transport parameter.
   **The differential crossover offset needed to reverse a pair's
   photoresponse is exactly that pair's work-function gap** — no model
   evaluation required, and Section 6.13's whole scalar uncertainty band
   exactly irrelevant to it. P1 (bisection lands on it to ≤1e-12 eV, and
   does not drift as `c` sweeps ±1 eV): **PASS** at 4.8e-15 and 9.7e-15 eV.

2. **RESULT — a parity factorisation, unplanned, that makes three earlier
   results corollaries of one identity.** P2 predicted `|N|` would be
   visibly asymmetric about the flip. It came back **0.00% on all 21
   pairs** — the signature of an identity, not of a small number. The
   identity:

       N(m, −s) = −N(m, s)      odd in the differential offset
       N(−m,  s) = +N(m, s)     even in the common offset

   both to 2.3e-15 over a 13×13 `(m, s)` grid. Consequences: odd-in-`s`
   **is** the τ* derivation above, in stronger form; even-in-`m` **is**
   2026-09-22's sign-invariance result, which was obtained by scanning
   12621 points over ±3 eV — a scan is evidence over the interval scanned,
   a parity is not restricted to one; and odd-in-`s` forces `|N|` to depend
   on `s` only through `|s|`, so **P2 and P4 are not near misses, they are
   excluded**. The note assumed an asymmetry without noticing it was
   assuming a symmetry away. Found the same way 09-22's Validation 5 was
   found: while trying to score a prediction, not while planning.

3. **Where the parity stops being exact, and why that is the useful part.**
   Even-in-`m` is exact *bitwise* at every bias, because the `−m` field is
   the literal spatial mirror of the `+m` field and `linspace(0, L, n)` is
   symmetric, so the mirror is exact on the grid. Odd-in-`s` is exact at
   zero bias and departs at finite bias by exactly `2/(n−1)` — measured at
   `n = 501, 1001, 2001, 4001, 8001`, residual × (n−1) = **2.0000**
   throughout, a sixteen-fold range. One trapezoid cell, weight 2 because
   the collection label jumps `+1 → −1` across it, converging as `1/n`. **A
   residual that scales with the grid is the grid; a residual that does not
   is physics.** Third time this chapter has needed that distinction, and
   the first time it was decided by a convergence study rather than by
   argument.

4. **RESULT — the pairs most at risk are exactly the pairs the model cannot
   reach, and this is structural rather than accidental.** Because τ* needs
   only `METAL_WORK_FUNCTIONS`, the margin table is exact for all 21 pairs.
   The four smallest margins are Au/Pd **0.020 eV**, Ni/Au 0.060, Ni/Pd
   0.080, Cr/Cu 0.150 — and every one contains Ni, Pd or Cr, i.e. a metal
   Section 6.12 *refuses* (chemisorbed, or no tabulated `d_eq`).
   Chemisorbed metals sit ~1.2 Å below the physisorbed anchor, which is
   both why 6.12.4's anchored exponential diverges on them and why their
   crossovers should differ most from everyone else's. **The model's reach
   and the thesis's risk are anti-correlated by construction.**

5. **RESULT — the three reachable pairs keep their sign, and the two
   channels are qualitatively unlike.** Cu/Au 2.85×, Cu/Pt 7.79×, Au/Pt
   18.6× against the model's largest differential offset over
   `ℓ ∈ [0.3, 1.5] Å`. P3 (Cu/Au worst, ratio in [2, 4]): **PASS** at 2.85.
   P5 (no pair below the smallest modelled offset): **PASS**. The honest
   comparison is not a ratio: the common channel has **no** threshold at
   any magnitude, so quoting "N× more dangerous" against an infinite margin
   is meaningless — an error written into 6.14.5 earlier in this same
   session and corrected in `21245f0` rather than silently dropped.

6. **The 2026-09-22 audit item discharged, with two detectors instead of
   one.** Last session's unbracketed bisection walked its lower bound to its
   upper bound and returned the endpoint as a root, 21 times, every one
   plausible, while all five of that module's exact validations passed.
   This session's root-finder raises on an unbracketed interval — it fired
   on **21/21** deliberately root-free intervals — and the root it finds is
   cross-checked against a closed form known *in advance*. A root-finder
   whose answer is known before it runs is the best possible test of that
   fault class, and it is the reason this item was done here rather than as
   a separate sweep.

7. **Writing** (`thesis_draft/06-graphene-photodetectors.md`,
   `thesis_draft/01-introduction.md`). New Section 6.14, eight subsections,
   two tables; Chapter 1 status table updated; Chapter 7 handed an
   **eighth** thread.

**Methodological note, continuing the series.** 09-20: exact validation does
not protect against an unrepresentative sample. 09-21: pre-registration
reaches what exact validation cannot. 09-22: pre-registration does not reach
the analysis layer either, and three failure classes went on record. Today is
the first entry on the other side of the ledger — not a way claims *break*
but a way they *consolidate*. One parity absorbed a 12621-point scan, a
pre-registered derivation and two live predictions. The distinction it
forces: **a claim checked on a wide sample and a claim derived from a
symmetry are not the same kind of knowledge**, and this thesis has been
accumulating the first while quietly assuming it was acquiring the second.
Chapter 6 now holds exactly one symmetry-derived result against six
survey-derived ones.

**Predictions scored:** P1 **PASS**, P3 **PASS**, P5 **PASS**, P2
**FALSIFIED**, P4 **FALSIFIED** — both by exclusion rather than by margin.
Two of five failing is the same rate as 09-22's two of six; the difference
is that today's failures were caused by a fact about the model that neither
the note nor any earlier section had noticed.

**Not yet covered (candidates for future runs):**
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21 and now the **top** item, promoted by today's result.
  Every number in the margin table's "reachable" column is one anchored
  exponential with a swept decay length; a second anchor is the only thing
  that would turn `τ_model` from an illustration into an estimate, and it is
  what stands between Result 4 and a real error bar on the sign of a
  photoresponse
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record now, and today made it sharper than
  ever: the metals the framework cannot describe are *precisely* the metals
  whose pairs have the smallest sign margins. Still Chapter 6's central
  structural weakness
- **Ask which of Chapters 4 and 5's design rules could be restated as
  parities or bounds rather than rankings over a tabulated set** — created
  today, and the concrete form of the eighth Chapter 7 thread. Chapter 4's
  `Rc` ranking and Chapter 5's liner comparison are both rankings; a parity
  or a ceiling would survive a widened sample where a ranking has twice not
- **Audit the repo's remaining analysis-layer code for the bracket bug
  class** — `alpha_from_separation`, the `ell` bisection of 2026-09-21,
  `contact_resistance_crossover.py` and `graphene_sensitivity_audit.py`.
  Today's new root-finder is guarded; the four older ones are not known to be
- **Whether the parity survives a photo-thermoelectric term** — created
  today. Odd-in-`s` is a property of *this* collection kernel. A PTE term is
  even in the temperature gradient and would enter with a different parity,
  which would make it the first thing in the chapter that could break the
  one clean symmetry it has
- **Re-check whether other "for every …" claims rest on small samples** —
  open since 2026-09-20; Chapter 4's per-metal `Rc` recalibration and
  Chapter 5's liner scenarios still state general rules from subsets
- **Whether Chapter 4's contact-resistance results should be re-run at the
  5.4 eV crossover** — open since 2026-09-18; Section 4.5 still uses the
  `|W − W_graphene|` magnitude convention
- **Chapter 7 (discussion/outlook)** — eight threads, a taxonomy of how
  claims fail and now one of how they consolidate. Genuinely ready to draft
- **Chapters 2–3 remain undrafted** despite complete computational results —
  still the largest remaining block of pure writing
- Reconciling Mueller *et al.*'s 0.12 eV step (arXiv:0902.1479) with the
  0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS` assumes — open since
  2026-09-18; bounded by the `S` table of 6.13 in the common channel and
  **not** bounded in the differential one, which is a new gap as of today
- A photo-thermoelectric term (Kasırga review) — see the parity item above
- Shimomura *et al.*'s comb-electrode design — open since 2026-09-19
- Integrating Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery
- Isolating the root cause of the Section 4.7 negative-residual result —
  open since 2026-08-31
- Ti and Cr per-metal `Rc` recalibration (ResearchGate rate-limiting)
- Second independent edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

---

## 2026-09-24 — The bracket-bug audit: three of four suspects were not root-finders, the one that was never produced a wrong number, and the fix did

**Status:** Automated session. Live web search **not used** — the session was
an audit of the repo's own analysis layer and every input was already here.
Taken from 2026-09-23's "Not yet covered" list: *"Audit the repo's remaining
analysis-layer code for the bracket bug class."* Not the top item (that
remains the second `Δ_c` anchor); chosen over it because an unguarded
root-finder can silently corrupt any number the repo later computes,
including the ones a second anchor would produce.

**The question.** 2026-09-22 found an unguarded bisection in
`graphene_crossover_sensitivity_model.py` that walked one bound onto the other
when no sign change existed and reported that bound as a root — 21 times, all
21 physically plausible, while all five of that module's exact validations
passed and printed underneath. 2026-09-23 answered it with a
`bracketed_bisect` that raises, and left an item asking whether anything else
in the repo had the same hole, naming four suspects.

1. **RESULT — the audit list was wrong about three of its four entries, and
   it took a test rather than a re-reading to find that out.**
   `crossover_length` is a closed form (`L_x = R_c W σ_sheet`, one
   evaluation); `alpha_from_separation` is affine in `d`; Family B's
   "calibration solve" is a single division. None can iterate, so none can
   have a bracket. Validation 1 establishes this by identity rather than by
   inspection, and each identity is false for anything iterative:
   `crossover_length(2R_c) == 2·crossover_length(R_c)` **bitwise**, α's
   midpoint equals the mean of its endpoints to **1 ulp**, and
   `λ_imp == λ_bulk / Σ(terms)` **bitwise**. The two misclassified entries do
   carry hazards — α goes negative for `d < d₀` (already guarded by
   `alpha_for_metal`), and the Family-B denominator can vanish — and they are
   a *different* class. **An audit list assembled by recalling which
   functions "solve for something" collects functions that merely compute
   something.**

2. **RESULT — the one real instance, and no published number was ever
   wrong.** `p2_threshold()` in `graphene_per_metal_crossover_model.py` ran
   200 halvings with no bracket check from 2026-09-21. At the shipped
   defaults the interval **is** bracketed (`f(lo) = +1.371`,
   `f(hi) = −0.090`), so Chapter 6's ℓ = 0.4997 Å is a genuine root with
   residual `−9.6e-16`. The fault was **latent**: `tol`, `lo` and `hi` are
   keyword arguments, and the first caller to move any of them out of the
   bracketed region gets a plausible number and no warning.

3. **Naming an unguarded solver is not showing it can produce a wrong
   number**, so Validation 3 pre-registered the two wrong numbers and then
   measured them: `tol = 2.0` (above `worst(lo) = 1.4706`) returns 0.0500 Å,
   `tol = 0.005` (below `worst(hi) = 0.009638`) returns 5.00 Å. Both are
   reasonable decay lengths. Neither is a root.

4. **A CORRECTION to the 2026-09-22 wording, and it matters more than it
   looks.** That note recorded that an unguarded bisection "returns the
   ENDPOINT". It does in one direction and **does not** in the other: when
   the loop takes `lo = mid` it lands on `hi` exactly, but when it takes
   `hi = mid` it **stalls one ulp above `lo`** — once `hi` is the next double
   after `lo`, `0.5·(lo + hi)` rounds back up to `hi`. Both directions were
   pre-registered; the ulp count was the prediction and it measured 1.0. So
   `result == lo or result == hi` looks like a cheap test for this fault and
   is **false half the time while the answer is meaningless anyway**. The
   earlier wording would have reassured exactly the caller who went looking.

5. **RESULT — the fix introduced a second fault, of a class no guard can
   see, and this is the session's real finding.** The first version failed its
   own Validation 5. The cause was not the guard but the guard's **default
   tolerance**: `bracketed_bisect(..., tol=1e-14)` compares an **absolute**
   interval width, chosen in `graphene_differential_crossover_model` for τ, a
   variable of order 1 eV, where it is ~machine precision (48 halvings).
   `p2_threshold` bisects a decay length in **metres**, order 5e-11, where the
   same 1e-14 stops after **16** halvings at a relative precision of **2e-4**
   — landing `2.5e11` ulps from the root with residual `−3.4e-6`, wrong in the
   fifth significant figure, and **nothing raised**. The call now passes
   `tol=0.0, rtol=eps`; `rtol` defaults to 0.0 so every pre-existing caller is
   bitwise unaffected. **A guard moved to a new call site carries its defaults
   with it, and a default tolerance is a claim about the scale of the caller's
   variable.** Reviewing the guard for correctness would never have found
   this — the guard was correct. It was caught only because the unguarded loop
   had been returning the right answer for three days and could be used as an
   oracle. **Deleting the old body before validating against it would have
   destroyed the only evidence that the new one was worse**, which is the
   operational form of this repo's standing rule about not quietly rewriting
   superseded numbers.

6. **RESULT — a load-bearing assumption converted from asserted to proved.**
   `p2_threshold`'s docstring claimed since 2026-09-21 that the |ΔW| shift is
   monotone decreasing in ℓ. A checked bracket buys nothing without it: a
   bisection on a non-monotone function can bracket a root and return the
   wrong one of several. It is provable term by term — `d_m − d_anchor` has a
   fixed sign per metal, so each `|exp(−(d_m − d_anchor)/ℓ) − 1|` approaches
   zero monotonically and a pointwise max of monotone non-increasing functions
   is monotone non-increasing. The sharper test is an exact consequence rather
   than the proof: **Pt's `d_eq` equals `D_ANCHOR` exactly (3.30 Å), so Pt's
   term is `0.0` bitwise at all 801 grid points and the max is really over Cu
   and Au alone.** Continuing 2026-09-23's distinction, this is the repo's
   **second** conversion of a sampled claim into a derived one.

7. **Housekeeping, verified rather than assumed.** `bracketed_bisect` now
   lives in `bracketed_root.py`; `graphene_differential_crossover_model`
   imports it and its local copy is replaced by a comment recording why the
   function exists. All six validations and all four results of that module
   were run before and after the move: printed output **bitwise identical,
   107 lines**. (Its two `FAIL` lines are 2026-09-23's P2 and P4, falsified
   then and unchanged now.)

8. **A scoring near-miss worth recording, because it is the DV repo's
   verdict-vs-checking class appearing here.** The first pass at re-running
   the differential module scored each validation by the truthiness of its
   return value. `validate_zero_perturbation` returns its worst deviation —
   `0.0` on success — so the *passing* case scored as FAIL. A harness that
   scores a numeric return as a boolean inverts exactly the validations that
   succeed perfectly. It was caught in the same call by comparing against the
   original file, which failed identically.

**Methodological note, continuing the series.** 09-20: exact validation does
not protect against an unrepresentative sample. 09-21: pre-registration
reaches what exact validation cannot. 09-22: pre-registration does not reach
the analysis layer. 09-23: the other side of the ledger — how claims
consolidate. Today adds a fourth failure class and it is the first one
*created by* a fix: **a correct guard, correctly applied, silently degraded a
correct number because its default tolerance encoded the scale of its
original caller.** Every numeric default in this repo — tolerances, step
sizes, grid counts — is such a claim, made at the call site where each
function was born. `log_sensitivity`'s `rel_step=1e-5` is relative and
therefore safe; 09-23's grid counts were settled by a convergence study.
Nothing else has been looked at.

**Predictions scored:** V3-A **PASS** (one ulp above `lo`, not `lo` — the ulp
count was the prediction), V3-B **PASS** (`hi` exactly). V5 was written as a
bitwise claim and **honestly weakened to 1 ulp** when it could not be one:
the old loop returns `hi` and the shared one returns the final midpoint, which
differ by one ulp for any function, so a bitwise claim there would have been a
claim about which of two adjacent doubles a loop names. 7/7 validations pass.

**Not yet covered (candidates for future runs):**
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21 and still the **top** item, untouched today. Every number
  in 09-23's margin table's "reachable" column is one anchored exponential
  with a swept decay length; a second anchor is the only thing that would turn
  `τ_model` from an illustration into an estimate. Today sharpened the stake:
  the ℓ = 0.4997 Å boundary is now a *verified* root of a *one-anchor* model,
  so its numerical trustworthiness has outrun its physical content
- **Audit every other numeric DEFAULT in the repo for the scale assumption
  it encodes** — created today, and the direct successor to the item just
  closed. `rel_step`, grid counts, `max_iter`, the `h=H_DIFF` finite
  difference in `graphene_crossover_sensitivity_model`, `tol` everywhere. The
  bracket class is closed; this class has one measured member and no detector
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record, and 09-23 showed the metals the
  framework cannot describe are precisely the metals whose pairs have the
  smallest sign margins. Still Chapter 6's central structural weakness
- **Ask which of Chapters 4 and 5's design rules could be restated as
  parities or bounds rather than rankings over a tabulated set** — created
  09-23, the concrete form of the eighth Chapter 7 thread, and today added a
  second symmetry-derived result to argue from
- **Whether the parity survives a photo-thermoelectric term** — created
  09-23. Odd-in-`s` is a property of *this* collection kernel; a PTE term
  enters with a different parity
- **Re-check whether other "for every …" claims rest on small samples** —
  open since 2026-09-20; Chapter 4's per-metal `Rc` recalibration and
  Chapter 5's liner scenarios still state general rules from subsets
- **Whether Chapter 4's contact-resistance results should be re-run at the
  5.4 eV crossover** — open since 2026-09-18; Section 4.5 still uses the
  `|W − W_graphene|` magnitude convention
- **Chapter 7 (discussion/outlook)** — now **eight threads plus today's
  fourth failure class**, and genuinely ready to draft. It is the largest
  piece of writing whose material is entirely in hand
- **Chapters 2–3 remain undrafted** despite complete computational results —
  still the largest remaining block of pure writing
- **`__pycache__` is tracked in this repo**, so every script run dirties the
  working tree and every automated session has to work around it. Repo
  hygiene, no bearing on any result, deliberately not folded into an audit
  commit — created today
- Reconciling Mueller *et al.*'s 0.12 eV step (arXiv:0902.1479) with the
  0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS` assumes — open since
  2026-09-18; bounded in the common channel by 6.13's `S` table and **not**
  bounded in the differential one
- A photo-thermoelectric term (Kasırga review) — see the parity item above
- Shimomura *et al.*'s comb-electrode design — open since 2026-09-19
- Integrating Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery
- Isolating the root cause of the Section 4.7 negative-residual result —
  open since 2026-08-31
- Ti and Cr per-metal `Rc` recalibration (ResearchGate rate-limiting)
- Second independent edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

**Automation health:** Device reachable, folder connected; neither repo had a
2026-09-24 entry, so a full session was run. **New constraint measured
today, and it changed how the repo is obtained:** a full `git clone` of this
repo over the session VM's proxy would not complete — four attempts, two with
`early EOF` / `invalid index-pack output` and two silently stalled, at one
point measuring **238 B/s** from codeload while the other repo cloned in 8.7 s.
The working recipe is a **partial + shallow + sparse clone**, which never asks
for the plot blobs at all:

```
git init -b main && git remote add origin <url>
git config core.sparseCheckout true
git config remote.origin.promisor true
git config remote.origin.partialclonefilter blob:none
printf '/*\n!*.png\n!*.jpg\n!*.pdf\n!*.gif\n!*.npz\n' > .git/info/sparse-checkout
git fetch --depth 1 --filter=blob:none origin main     # 5.6 s
git update-ref refs/heads/main FETCH_HEAD && git checkout main   # 20 s
```

26 seconds instead of never. Two consequences for a future run: the checkout
reports "69% of tracked files present" and that is expected, not damage; and
a regenerated plot must be committed deliberately, since PNGs are
`SKIP_WORKTREE`. `scipy` is still absent from the device VM and still needs
installing (`pip install scipy`) — it failed once at 13.6/37.7 MB and
succeeded on retry, because pip resumes partial downloads. The bridge dropped
twice mid-session; both times the in-flight command had **not** executed, and
`git log` in each clone was the reliable check.

**Commits this run:** 5 (the shared guard module with the differential
module's migration; the `p2_threshold` fix with the tolerance-scale fix; the
audit with its recorded output; the study note; the Chapter 6.12.3
annotation). This AUTOMATION_LOG.md entry makes 6. The guard extraction and
the fix are separate commits because the second one carries a finding of its
own, and a fix that introduced a fault should be visible as that in the
history.

---

## 2026-09-25

**The item closed:** *"Audit every other numeric DEFAULT in the repo for the
scale assumption it encodes — the bracket class is closed; this class has one
measured member and **no detector**."* Created 2026-09-24 as the direct
successor to the item that session closed. The detector exists now, and it
convicted a default that has been feeding a published table for four days.

Code: `graphene_default_scale_audit.py` (929 lines), output
`default_scale_audit_output.txt`, figure `default_scale_audit.png`. Follow-on
code: `sensitivity_converged()` in `graphene_crossover_sensitivity_model.py`,
output `sensitivity_converged_output.txt`. Study note:
`notes/2026-09-25-numeric-default-scale-audit.md`. Thesis annotation: new
Section 6.13.11, plus an annotation block inside 6.13.4.

1. **THE DETECTOR, and why "no detector" was the operative half of that item.**
   The 09-24 fault was found only because the unguarded loop it replaced had
   been returning the right answer for three days and could serve as an
   **oracle**. Most defaults in this repo have no older correct implementation
   standing beside them, so a technique needing one detects nothing. The
   detector built today is a **unit-covariance probe**: restate a procedure's
   problem with its variable in different units and require the answer to
   transform covariantly, `P_λ(f_λ) == λ·P(f)`. A default that is a pure ratio
   satisfies this identically; one carrying the dimensions of the variable does
   not. It needs no oracle, no reference implementation and no correct answer —
   it tests a procedure against its own units — and it **reproduces the 09-24
   fault with nothing borrowed from it.**

2. **RESULT — the census, and the 09-24 item conflated two classes.** 19 numeric
   defaults enumerated and sorted: **5 TOLERANCE-like** (the probe reaches
   them), **4 COUNT-like** (a dimensionless count passes the probe trivially;
   only a convergence study reaches what it claims), **10 PHYSICAL** (`T=300 K`,
   `hopping=2.8 eV`, `sigma=20 nm` — statements about graphene, not about
   numerics). The PHYSICAL rows are enumerated anyway so the audit's boundary is
   auditable rather than asserted. Sharpest instance: `p2_threshold(tol=0.10)`
   *looks* like a numerical tolerance and is not — 10% is P2's definition, and
   tightening it changes the question. The two unstudied counts were studied:
   `n_points=400` is converged to 6.6e-06 against a 25 600-node reference (D6
   held), and `max_iter=200` has a 3.3× margin over the 61 halvings the widest
   interval here uses (D5 held).

3. **RESULT — the probe's five verdicts all matched pre-registration.** Two
   defaults are scale-bound: `bracketed_bisect`'s `tol=1e-14` (the known member,
   4.9e-02 at λ=1e-3) and **`H_DIFF`** (2.5e+01). The `rtol` call, the
   pure-halving call and `log_sensitivity`'s `rel_step` are scale-free at
   **exactly 0.0** deviation.

4. **A property of the probe, learned from a failed prediction.** D1 predicted
   the deviation band `1e-6 … 1e-2`; measured 4.9e-02 at λ=1e-3 and 4.8e-05 at
   λ=1e+3. The band was written from the λ>1 direction alone and **the probe is
   asymmetric in λ**: shrinking the variable makes an absolute tolerance coarse
   relative to the interval (catastrophic), while growing it makes the tolerance
   fine, so the residual deviation is then dominated by the error of the λ=1
   **baseline** rather than of the rescaled run. A deviation from this probe is a
   reliable yes/no and **not a calibrated error estimate.** Same asymmetry
   `bracketed_root.py` already records for the unguarded loop, in a new guise.

5. **RESULT — the screen is not the verdict, and separating them is what kept a
   correct number correct.** Failing covariance means a default *encodes* a
   scale, not that it is wrong. Stage (a), the screen, is cheap and needs
   nothing; stage (b), a plateau or convergence measurement at the shipped
   scale, is the only stage that can convict. Reporting (a) as the verdict would
   condemn correct numbers, which is the **mirror image** of the 09-24 fault and
   the more expensive mistake, because it moves published tables. Of the two
   flagged, `tol=1e-14` was acquitted (its one metre-scale caller passes `rtol`)
   and `H_DIFF` was convicted.

6. **RESULT — `H_DIFF = 1e-3 eV` sits 3.5 DECADES above the top of its own
   plateau.** The plateau of `d ln|N|/dδ` was measured for the first time today
   and runs ~`1e-10` to `3e-7 eV`. **D4 was written as a null result precisely
   so it could fail, and it failed.** The damage is uneven and the unevenness is
   the story: all six straddling pairs agree to better than 1 part in 1e5; four
   same-sign pairs move under 0.4%; the remaining **eleven move 5.9% to 44.4%**
   (Ti/Ni 44.38, Cr/Ni 38.90, Cr/Pd 31.99, Cr/Au 23.83, Cu/Au 21.65, Ti/Cu
   19.35, Ti/Cr 13.89, Cu/Ni 10.16, Cu/Pd 8.65, Cr/Cu 7.85, Ni/Au 5.90) — and
   **not one-signed**: Cr/Ni and Cu/Au were *under*-estimated, so no uniform
   correction factor would have rescued the table. **Section 6.13.4's dichotomy,
   that section's only load-bearing result, SURVIVES** at 26.7× rather than
   30.4×, with no overlap. That is the outcome to lead with.

7. **RESULT — prediction P2 is FALSIFIED, and in exactly the way its own
   rationale said it might be.** P2 (pre-registered 2026-09-22, recorded there
   as **HELD**, verdict row now annotated in place): the largest `|S|` at
   `δ = +0.1 eV` is Au/Pd, because amplification is near-cancellation and Au/Pd
   is the nearest-cancelling pair (0.02 eV apart). Converged: **Ti/Cu 2.828
   (0.32 eV apart), Cr/Au 2.518 (0.60 eV), Ni/Pd 2.291 (0.08 eV), Au/Pd fourth
   at 2.016** — identical at every step inside the plateau, over four decades.
   The two closest pairs place third and fourth; a pair **16× further apart
   wins**. Near-cancellation is **not monotone in contact separation at finite
   δ**. The 09-22 note had written the escape hatch into the prediction —
   *"near-cancellation could be non-monotone in separation once the kernel's
   spatial structure enters, in which case some other pair wins"* — and some
   other pair wins. That is what converts this from a surprise into a scored
   falsification, and it is the strongest evidence yet for the pre-registration
   habit adopted 09-21. **The mechanism generalises past this table:** Au/Pd's
   `|S|` drifts **2%** across five decades of step size while Ti/Cu's collapses
   **69%**, so Au/Pd led at `h=1e-3` not by being most sensitive but by being
   the pair that step suited. **A ranking taken outside the plateau ranks its
   entries by how well the step suits them.**

8. **RESULT — the fifth failure class, and the first in the series where the FIX
   for the previous class was itself the fault.** The obvious fix was already in
   the repo: `log_sensitivity` has returned `|S(h) − S(2h)|` with every value
   since 09-22. The first `sensitivity_converged()` adopted that pattern and it
   is **blind**: at `h=1e-3`, Ti/Ni's true error is **44.38%** and its
   step-doubling estimate is **2.33e-06** — under-reporting by **190 268×** (also
   Cu/Ni 10.16% vs 5.48e-05, Cr/Pd 31.99% vs 3.96e-04). The mechanism is exact:
   step-doubling measures `dE/d(log h)` of the error curve `E(h) = S(h) − S(0)`
   and reads **zero at any stationary point** of it, and `S(h)` for Ti/Ni is flat
   to four digits from `1e-3` to `1e-2`. **A shipped default has no reason to
   avoid a stationary point of its own error curve, and landing on one makes
   every local self-consistency test agree with it.** The criterion is therefore
   **anchored, not self-consistent** — accept a step only if the derivative is
   unchanged at `h/10` *and* `h/100`. Measured: **11 of 11** moved pairs caught
   at `H_DIFF`, **zero escapes**, Ti/Ni included; all 21 pass at
   `H_DIFF_CONVERGED = 1e-8`, worst residual 9.6e-05; the two populations are
   separated by **2.82 decades** and `rtol = 1e-3` is placed inside that measured
   gap rather than chosen for roundness — an odd thing to do on the day the repo
   learned a threshold is a default and a default is a claim.

9. **Two design points recorded in the code because they are not trivia.**
   `H_DIFF_CONVERGED` is `1e-8` and not `1e-7` because at `1e-7` two pairs fail a
   `1e-6` requirement — Cr/Ni 3.98e-05 and Cu/Au 5.07e-05, against 1e-9…1e-7 for
   the other nineteen — and both drop exactly 100× when `h` drops 10×, so it is
   ordinary second-order truncation and those two pairs simply have a third
   derivative ~100× the rest. **No single global step is right for every pair,
   and the only reason that is visible is that the function returns its
   estimates instead of asserting them.** And the residual 9.6e-05 at the
   accepted step belongs to Pd/Pt, whose `|S|` is 0.0077: two decades down at
   `1e-10` it reaches the **cancellation floor**, not any truncation.

10. **This audit's own validation committed the audited error.** Validation 4
    asserted a central difference of an exactly linear function recovers the
    slope to better than `1e-14` relative, and **measured 4.8e-11 and failed.**
    The arithmetic was right and the assertion wrong: cancellation leaves
    `eps·|c₀|` in the numerator, so the bound is `eps·|c₀|/(2h|c₁|)` — 7.7e-15 at
    `h=1e-3`, 7.7e-10 at `h=1e-8`. **The threshold `1e-14` was an absolute
    constant silently encoding `h ≈ 1e-3`: the audited fault class, one level up,
    inside the auditor.** Caught by the test failing, not by review, and recorded
    rather than quietly corrected.

11. **One validation retained BECAUSE IT CANNOT DISCRIMINATE.** Bisecting an odd
    function on a symmetric interval returns bitwise `0.0` after one evaluation
    at any tolerance, any width, with or without `rtol`. Exact, and it passes for
    the good default and the bad one alike — it proves the harness and nothing
    about the default. Recording a non-discriminating test **as**
    non-discriminating is cheaper than rediscovering that it was; this repo has
    twice mistaken one for the other (09-22's analysis-layer fault, which every
    model-level validation passed underneath; 09-24's truthiness-scored
    validation).

12. **Nothing was quietly rewritten.** `H_DIFF` stays at `1e-3` and
    `sensitivity()` is bitwise unchanged — the whole module's printed output is
    **identical before and after**, 100 lines. Section 6.13.4 keeps its table and
    gains an annotation; the 09-22 note keeps its `HELD` and gains a retraction
    beside it. The converged accessor is additive.

**Methodological note, continuing the series.** 09-20: exact validation does not
protect against an unrepresentative sample. 09-21: pre-registration reaches what
exact validation cannot. 09-22: pre-registration does not reach the analysis
layer. 09-23: the other side of the ledger — how claims consolidate. 09-24: a
correct guard silently degraded a correct number because its default encoded the
scale of its original caller. **09-25: a procedure asked whether it has
converged can answer yes and be wrong by 44%. Self-consistency near a stationary
point of the error curve is not evidence; only an anchored comparison against a
step known independently to be in the plateau distinguishes them.** And an
ordering worth stating, because it was not obvious in advance: the screen came
first and cost nothing, the margin measurement second and cost seconds, and the
published-claim check last — and it is the only one that produced a correction.
The screen alone flags `H_DIFF` without knowing whether to care; the margin study
alone shows a bad step without knowing which claim it reaches. Neither alone is
an audit.

**Predictions scored:** D2, D3, D5, D6 **PASS**. D1 **FAIL** on its magnitude
band and not its class call (see item 4). D4 **FAIL**, written as a null result
so that it could (see item 6). **5/5 validations pass**, after Validation 4 was
rebuilt on a derived bound following its first-form failure.

**Not yet covered (candidates for future runs):**
- **`output_conductance(dVds=1e-3)` at `Vds = 0.05 V`** — created today and the
  **top numerical item**. The identical construction to `H_DIFF`, a step 2% of
  its variable, found by the census rather than the probe. Deliberately not
  measured today because `g_ds` feeds `f_max`, a published number, and moving it
  needs its own before/after comparison. The instrument now exists and takes
  minutes per default
- **`log_sensitivity`'s own convergence estimate IS step-doubling** — created
  today and sharper than it looks: the Chapter 4 `R_transmission` decomposition
  and Chapter 5 liner scenarios rest on the pattern item 8 shows to be blind by
  five orders of magnitude, and have never been checked against an anchored
  criterion
- **Whether any Chapter 4 or 5 RANKING is a step artefact of the same kind** —
  created today. The 6.13 ranking was, and "a ranking orders its entries by how
  well a shared numerical choice suits each" is not specific to derivatives or
  to Chapter 6
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open since
  2026-09-21 and still the top *physics* item, untouched today. Every number in
  09-23's margin table's "reachable" column is one anchored exponential with a
  swept decay length; a second anchor is the only thing that would turn
  `τ_model` from an illustration into an estimate
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; still Chapter 6's central structural
  weakness. Today's result touches it: Ti/Cu and Ti/Ni are among the pairs whose
  sensitivities moved most, and Ti is one of the three metals the framework
  cannot describe
- **Ask which of Chapters 4 and 5's design rules could be restated as parities
  or bounds rather than rankings over a tabulated set** — created 09-23, and
  today is a strong new argument *for* it: a parity or a bound is not a ranking,
  and today showed a ranking can be an artefact of a shared numerical default
- **Whether the parity survives a photo-thermoelectric term** — created 09-23
- **Re-check whether other "for every …" claims rest on small samples** — open
  since 2026-09-20
- **Whether Chapter 4's contact-resistance results should be re-run at the
  5.4 eV crossover** — open since 2026-09-18
- **Chapter 7 (discussion/outlook)** — now **nine threads plus today's fifth
  failure class**, and still the largest piece of writing whose material is
  entirely in hand. Planned for today and displaced by the P2 falsification,
  which was worth the displacement
- **Chapters 2–3 remain undrafted** despite complete computational results
- **`__pycache__` is tracked in this repo** — open since 09-24 and it cost real
  time today: a `git stash` used for the before/after regression check refused
  to pop because the stale `.pyc` files had been rewritten by the intervening
  run. Recoverable (`git checkout -- __pycache__` then pop) but it is now a
  hazard and not only an annoyance
- Reconciling Mueller *et al.*'s 0.12 eV step (arXiv:0902.1479) with the
  0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS` assumes — open since 2026-09-18
- A photo-thermoelectric term (Kasırga review); Shimomura *et al.*'s
  comb-electrode design; integrating 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery; isolating the root cause of the
  Section 4.7 negative residual (open since 2026-08-31); Ti and Cr per-metal
  `Rc` recalibration (ResearchGate rate-limiting); a second independent
  edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

**Automation health:** Device reachable, folder connected; neither repo had a
2026-09-25 entry, so a full session was run. **The 09-24 clone constraint did
not reproduce:** a plain `git clone` of this repo completed normally and fast,
so the partial/shallow/sparse recipe recorded yesterday was not needed. That
makes yesterday's 238 B/s measurement look like a transient proxy condition
rather than a standing property of this repo — the recipe is kept in yesterday's
entry as a fallback, not a requirement. `git config user.name/user.email` was
again absent in the fresh clone and was set **in its own call**, per 09-24's
finding. `scipy` was again absent from the device VM and `pip install scipy`
succeeded first time (1.15.3). No live web search was used today: the session's
work was entirely internal — an audit of this repo's own numbers — and needed
none.

**Commits this run:** 4 (the audit with its output and figure; the
plateau-anchored accessor with its recorded output; the thesis 6.13 annotation
and P2 retraction; the study note with the 09-22 note's in-place annotation).
This AUTOMATION_LOG.md entry makes 5. The audit and the accessor are separate
commits because the accessor's first version was wrong in a way that is a result
of its own, and a fix whose obvious form was blind should be visible as that in
the history — the same reason 09-24 split its guard extraction from its fix.

---

## 2026-09-26 — The top numerical item closes as a NULL result, and the closure is only meaningful because the mechanism was shown to reach 1.32% elsewhere

**The item closed:** *"`output_conductance(dVds=1e-3)` at `Vds = 0.05 V` —
created today and the **top numerical item**. The identical construction to
`H_DIFF`, a step 2% of its variable... Deliberately not measured today because
`g_ds` feeds `f_max`, a published number, and moving it needs its own
before/after comparison."* Created 2026-09-25. It is measured now, and the
answer is that the default is innocent — together with a second default nobody
had noticed was entangled with it, a real latent bug in the surrounding
function, and a systematic error in the instrument that created the item.

Code: `graphene_gds_quadrature_audit.py` (~700 lines), output
`gds_quadrature_audit_output.txt`, figure `gds_quadrature_audit.png`. Fix:
`rf_small_signal_model.output_conductance`. Study note:
`notes/2026-09-26-gds-step-quadrature-audit.md`. Thesis: new Section 4.6.1,
plus an annotation inside 4.6 and a 4.9 status-table update. In-place
annotations: `graphene_default_scale_audit.py` RESULT 6 docstring and
`notes/2026-09-25-…` §10.

1. **THE STRUCTURAL DIFFERENCE FROM `H_DIFF`, which is why this was a session
   and not a footnote.** `H_DIFF` differentiates a closed form. `dVds`
   differentiates `transfer_characteristic()`, which is **itself a 50-point
   quadrature** over the drain-bias drop — and the quadrature grid *is set by
   `Vds`* (`linspace(0, Vds, n_segments)`). Differencing in `Vds` differences
   the quadrature error too. **Two defaults are entangled and not
   independent**, which nothing in the 09-22…09-25 series had encountered, and
   which makes a failure mode available that 09-25's fix does not reach: the
   anchored step criterion tests convergence **in the step**, so a step study
   of a discretised function converges to the derivative of *the
   discretisation you fixed* and is blind to its bias **at any tolerance, by
   construction.**

2. **RESULT — both defaults innocent, four to five decades clear.** Step vs.
   anchored-plateau step: **1.25e-6**. `n_segments = 50` vs. `n → ∞`
   (Richardson in 1/(n−1), n up to 6400): **8.4e-7**. Propagated to peak
   `f_max`: **1.1e-7** (0.00001%). `f_T` contains no `g_ds` and does not move
   at all. The quadrature error is confirmed **first order** in 1/n —
   difference ratios under n-doubling 2.031, 2.015, 2.008, 2.004, 2.002, 2.001
   — because `np.mean` over an endpoint-inclusive `linspace` is not the
   trapezoid rule. **Q2 FAILED on its magnitude band** (predicted ≥ 0.1%) and
   not on its class call. The item is **CLOSED, not carried** — the first
   clean exoneration in the series that began 09-22.

3. **WHY THE NULL RESULT IS NOT MERELY REASSURING, which is the actual content
   of the session.** "We refined the step and nothing moved" is *compatible
   with the criterion being blind*, so on its own it closes the item for the
   wrong reason. Two measurements fix that. **(a)** The `n=50` bias is
   **invariant under step refinement** on the real model: +3.089e-07 at
   `dVds=1e-2`, +3.039e-07 at 1e-3, +3.038e-07 at 1e-5, +3.038e-07 at 1e-7 —
   drifting **2.3e-4 of itself across four decades**. Refining the step
   removes none of it, exactly as item 1 predicts. **(b)** In an
   exactly-solvable case the same mechanism is **1.32%**: for
   `R(V_ch) = a + b·V_ch²` the discrete mean of `V_ch²` over an
   endpoint-inclusive grid is `V²(2n−1)/(6(n−1))` against a continuum `V²/3`,
   so the bias is exactly `V²/(6(n−1))`; measured `c_n` matched the closed
   form to 3.3e-16, the central difference converged to the **discretised**
   closed form at second order in `h` (order test 3.9995), and its residual
   against the **continuum** form converged to 1.3220563663e-02 where the
   exact gap is 1.3220563698e-02 — **to the gap, not to zero.** So the
   criterion can be satisfied to any tolerance while the answer is 1.32%
   wrong. **The exoneration is a property of graphene's nearly-linear
   `R(V_ch)` across a 50 mV drop, not a licence.** Validation C pins that
   down oracle-free: the endpoint-inclusive mean is exact for a constant *and*
   exact for a linear integrand (1.6e-16 over n ∈ {2,3,50,501}), failing only
   from curvature onward.

4. **RESULT — a real latent bug, convicted by an exact algebraic factor.**
   `max(Vds − dVds, 1e-4)` moved the **interval** without changing the
   **divisor**, so when it fired the function returned the true secant slope
   times exactly `(Vds+dVds−1e-4)/(2·dVds)`. Validation E reproduced that to
   **rel err 0.0**: measured 0.9158333333333334, exact 0.9158333333333334.
   Against the anchored reference at `V_g = 2.0 V`: **−8.46%**, silently.
   **Q4 FAILED on its magnitude band** (predicted > 10%) and not its class
   call. Reachable for `Vds ≤ dVds + 1e-4` (1.1 mV at the default step) **and
   from any upward step-size sweep — i.e. from exactly what a step-size audit
   does.** The bug was one careless sweep away from contaminating an audit of
   itself. It never fired at either operating point in use, so no published
   number was ever affected.

5. **The fix shrinks the step instead of widening the interval, and that is
   the substantive choice.** Widening was the smaller edit, but an asymmetric
   secant over `[1e-4, Vds+dVds]` is a second-order estimate of `dId/dVds` at
   that interval's **midpoint** — the caller who asked for `g_ds` at `Vds`
   would have received `g_ds` somewhere else, correctly computed. Shrinking
   keeps it centred: **+0.0017%** instead of −8.46% at the same point. A
   shrink now **warns** (09-24: a correct guard that silently degrades a
   correct number is worse than one that fails loudly) and `Vds ≤ 1e-4`
   raises. The non-firing path is unchanged expression by expression and
   verified **bitwise — not to a tolerance** — over
   `Vds ∈ {0.05, 0.1, 0.2} × dVds ∈ {1e-4, 1e-3, 1e-2}` on the 400-point
   sweep, with warnings promoted to errors so an unexpected shrink would have
   failed the check.

6. **UNPREDICTED RESULT — 09-25's census read a SIGNATURE, not a call site,
   and this is a systematic error across the whole list it produced.** This
   module's peak `f_T` came out **10.140 GHz** against Chapter 4's published
   **≈20 GHz**. Neither is wrong: `output_conductance`'s *signature* default
   is `Vds = 0.05 V`, while `plot_fT_fmax()` — the path that produces
   `rf_figures_of_merit.png` and the chapter's numbers — passes
   **`Vds = 0.1 V`**, giving 20.279 / 18.731 GHz. §4.6 had never named which.
   The consequence: the census scored `dVds` as "a step **2%** of its
   variable"; at the call site it is **1%**, low **by exactly the factor
   between the signature default and the call site** — 2× here and in general
   unbounded. **A default-scale census must read call sites, not signatures.**
   Recorded as unpredicted rather than folded into a prediction after the
   fact. Both the census code and the 09-25 note are annotated in place; the
   2e-2 stays.

7. **A number worth recording because it is counter-intuitive.** The `g_ds`
   term carries **0.9999** of the `f_max` denominator at peak `f_max`, so the
   dilution should have been the square root's 0.5. It measured **0.229**.
   The remainder is a property of reporting the peak value of a near-flat
   maximum, not of the physics: **a denominator weight of 0.9999 does not
   imply a sensitivity weight of 0.9999 in the reported figure of merit.**
   Q5 passed, for a reason partly other than the one predicted.

8. **A SECOND ITEM CLOSED, cheaply.** Q6: the 400-point `V_g` grid behind
   `g_m = np.gradient(Id, V_g)`, and hence behind the published peak `f_T`,
   moves it **0.001%** under 16× refinement (10.14017 → 10.14026 GHz). A
   third entangled default, created and closed in the same session.

9. **THE AUDITOR COMMITTED THE AUDITED ERROR, FOR THE SECOND TIME IN TWO
   DAYS.** Validation **[B]** required its residual below an absolute `1e-12`
   and measured **6.08e-10** — which **is** the derived cancellation floor
   `eps·|Id|/(2h|g|)` at `h=1e-9`, so the constant silently encoded
   `h ≳ 1e-7`. That is 09-25 item 10's fault class, in a fresh module, **one
   day after it was documented, by the same process that documented it.**
   Validation **[A]** likewise asserted a round "100× larger" and measured
   67×. Both rebuilt on derived quantities — an oracle-free order test and the
   exact algebraic gap — and both failures left in the source rather than
   quietly corrected. **The conclusion is not "be more careful":** knowing the
   class demonstrably does not prevent reproducing it, and what caught it both
   times was **a test that reports a number instead of asserting a verdict.**

10. **One validation retained BECAUSE IT CANNOT DISCRIMINATE** (continuing
    09-25 item 11). **[B]** is exact for every `(n, h)` tried, for a constant
    `R(V_ch)`, and passes for the shipped default and an absurd one alike: it
    proves the harness and nothing about the default. Labelled as such in the
    source. **[D]** states its cancellation bound as a derivation and reports
    measured/bound ratios (0.62, 1.06) rather than clearing a constant.

11. **Nothing was quietly rewritten.** `dVds` and `n_segments` are unchanged.
    §4.6's "≈20 GHz" sentence is kept verbatim with its operating point added
    beside it. The 09-25 census keeps its 2e-2 and gains an annotation. The
    only behavioural change in the repo is the guard fix, on a path that never
    ran.

**Methodological note, continuing the series.** 09-20: exact validation does
not protect against an unrepresentative sample. 09-21: pre-registration
reaches what exact validation cannot. 09-22: pre-registration does not reach
the analysis layer. 09-23: how claims consolidate. 09-24: a correct guard
silently degraded a correct number because its default encoded the scale of
its original caller. 09-25: a procedure asked whether it has converged can
answer yes and be wrong by 44%. **09-26: the anchored comparison is anchored
in ONE variable. Refining a step converges to the derivative of whatever other
discretisation was held fixed, and when two defaults are entangled —
and they are entangled whenever one sets the other's grid — convergence in one
is not evidence about the pair.** Second, and earned the hard way today: **a
null result is worth a session only if it is a null result about a mechanism
you have shown can be large.** Both formulations close the item; only one
closes it for a reason.

**Predictions scored:** Q1, Q3, Q5, Q6 **PASS**. Q2 **FAIL** on its magnitude
band (≥0.1% predicted, 0.0001% measured), class call correct. Q4 **FAIL** on
its magnitude band (>10% predicted, 8.42% measured), class call correct.
**5/5 validations pass**, after [A] and [B] each failed in their first form
and were rebuilt on derived bounds.

**Not yet covered (candidates for future runs):**
- **`log_sensitivity`'s own convergence estimate IS step-doubling** — created
  09-25 and now the **top numerical item** by default, since today closed the
  one above it. Chapter 4's `R_transmission` decomposition and Chapter 5's
  liner scenarios rest on the pattern 09-25 showed blind by five orders of
  magnitude, and have never been checked against an anchored criterion
- **Every remaining `h`-like step and tolerance in `rf_small_signal_model.py`
  and the photodetector modules** — open since 09-25, and today sharpened the
  method rather than the list: **score them at call sites.** Item 6 makes the
  census's signature-reading a systematic error across that entire list
- **Whether any Chapter 4 or 5 RANKING is a step artefact of the same kind** —
  created 09-23/09-25, untouched today, and today adds nothing for or against
  it
- **`n_segments = 50` at a bias with more curvature** — created today. Today
  exonerated it for `g_ds` at two bias points, but item 3 showed the
  exoneration is curvature-dependent, and `n_segments` sets **every** `Id` in
  Chapter 4. The follow-up is whether any other Chapter 4 quantity
  differentiates through that quadrature where `R(V_ch)` curves more
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21 and still the top *physics* item, **untouched for six
  consecutive sessions**. Every "reachable" number in 09-23's margin table is
  one anchored exponential with a swept decay length
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural
  weakness
- **Ask which of Chapters 4 and 5's design rules could be restated as parities
  or bounds rather than rankings over a tabulated set** — created 09-23
- **Whether the parity survives a photo-thermoelectric term** — created 09-23
- **Re-check whether other "for every …" claims rest on small samples** — open
  since 2026-09-20
- **Whether Chapter 4's contact-resistance results should be re-run at the
  5.4 eV crossover** — open since 2026-09-18
- **Chapter 7 (discussion/outlook)** — now **ten threads**, material entirely
  in hand, and displaced for a seventh session. It is the largest piece of
  writing whose material is complete, and the numerical thread has now
  produced two consecutive sessions that displaced it
- **Chapters 2–3 remain undrafted** despite complete computational results
- **`__pycache__` is tracked in this repo** — open since 09-24; it dirtied the
  working tree on every run today and cost real time on 09-25
- Reconciling Mueller *et al.*'s 0.12 eV step (arXiv:0902.1479) with the
  0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS` assumes — open since 2026-09-18
- A photo-thermoelectric term (Kasırga review); Shimomura *et al.*'s
  comb-electrode design; integrating 6.5's plasmonic near-field picture with
  the spatially-resolved contact-doping machinery; isolating the root cause of
  the Section 4.7 negative residual (open since 2026-08-31); Ti and Cr
  per-metal `Rc` recalibration (ResearchGate rate-limiting); a second
  independent edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

**Automation health:** Device reachable at the 04:34 UTC firing, folder
connected; neither repo had a 2026-09-26 entry and neither had commits since
midnight, so a full session was run. A plain `git clone` of both repos
completed normally and fast (the 09-24 partial/shallow recipe was again not
needed). `git config user.name/user.email` was again absent in the fresh
clones and set in its own call per 09-24's finding. `scipy` was again absent
from the device VM; `pip install scipy` succeeded first time (1.15.3), though
nothing in today's work needed it. No live web search was used: the session
was an internal audit of this repo's own numbers and needed none.
`__pycache__` again showed as modified in `git status` on every run (see the
open item).

**Commits this run:** 4 (the audit with its output and figure; the guard fix
with its bitwise regression check; the thesis 4.6 annotation and new 4.6.1;
the study note with the two in-place annotations). This AUTOMATION_LOG.md
entry makes 5. The audit and the fix are separate commits because the audit
is a null result and the fix is not — a session whose headline is "nothing
moved" should not bury a −8.46% bug inside that commit.

---

## 2026-09-27 — The convergence estimate reports machine precision where the answer is 14% wrong, and Chapter 7 is finally drafted

**Status:** Automated session. Live web search **not used** — the numerical half
was an internal audit of this repository's own estimator and the writing half
was a synthesis of results already in hand; neither needed it. Device reachable
at the **04:34 UTC** firing, folder connected; neither repo had a 2026-09-27
entry and neither had commits since midnight, so a full session was run.

**Pre-registration committed first** (commit `e2eed30`, before the audit module
existed): `notes/2026-09-27-log-sensitivity-preregistration.md`, seven questions
Q1–Q7 and five validations V1–V5, continuing the practice established 09-21.

**The item.** `graphene_sensitivity_audit.log_sensitivity(f, p, rel_step=1e-5)`
is the only finite-difference estimator in this thesis's conditioning work. It
reports `conv = |S(h) − S(2h)|` with every value — step-doubling, the pattern
09-25 caught answering *yes* while the answer was 44% wrong — and had been the
**top open numerical item since 09-25**.

1. **A correction to the session's own pre-registration, before any result.**
   Q1 and Q4 named "Chapter 4/5 call sites". **There are no Chapter 4 call
   sites.** Section 4.10 is entirely term-level and closed form
   (`S_i = T_i/Q`, no finite differences, no tolerance), so this estimator's
   whole exposure is Chapter 5. The pre-registration asserted a call site it had
   not checked for — a smaller version of 09-26 item 6, where a census read
   signatures instead of call sites. Recorded here rather than silently
   narrowed.

2. **RESULT — the estimator was never the central difference its docstring
   claimed, and this is 09-26 item 5's fault class on a DEFAULT path.** It
   sampled `p(1 ± h)`. In `L = ln p` those points sit at `L + ln(1+h)` and
   `L + ln(1−h)`, and `ln(1+h) ≠ −ln(1−h)`, so it was a **secant over an
   asymmetric log interval**, returning `dln|f|/dln p` at the interval midpoint
   `L + ½ln(1−h²) = L − h²/2 + O(h⁴)`. Its leading error carries
   `−(h²/2)g″` — proportional to `dS/dln p` — a term a centred difference's
   error does not involve at all. 09-26 found the identical class on a guard
   path that had **never fired**; this is the path behind **every** published
   sensitivity in Section 5.5.

3. **Proved by closed form, not asserted.** For `ln|f| = A(ln p)²` the old
   estimator's error is `A·ln(1−h²)` **exactly, at every order**, while a
   centred difference on the same function is **exactly zero** — a centred
   difference of a quadratic has no error. One function separates the two
   estimators by a closed form at every step, and the measured agreement sits at
   the derived cancellation bound throughout (worst measured/bound **2.08**).
   Magnitude at the shipped step: **9.47e-10** relative, worst over 23 call
   sites, so Q1's class call and magnitude band both **PASS** and no number
   moves.

4. **RESULT — THE HEADLINE, and it falsifies this session's own Q4.** Q4
   predicted **no** zero of `conv` at any real call site in `h ∈ [1e-8, 1e-1]`,
   reasoning that these are smooth rationals where `g″` dominates `g‴`. **Four
   zeros lie on the truncation branch**, located by bracketed bisection (09-24's
   rule) against a derived rounding floor, and **all four are in Family B — the
   λ_impurity calibration, `κ = 19.71`, the step Section 5.5 calls the
   worst-conditioned in this thesis**:

   | call site | `h` at `conv = 0` | `conv` there | error in `S` | relative |
   |---|---|---|---|---|
   | `λ_imp/ρ_calibration` | 2.928e-2 | 4.97e-14 | 2.76 | **14.0%** |
   | `λ_imp/ρ_bulk` | 3.084e-2 | 1.42e-14 | 2.62 | **13.3%** |
   | `λ_imp/W_calibration` | 4.391e-2 | 1.24e-14 | 1.83 | **15.1%** |
   | `λ_imp/λ_bulk` | 4.752e-2 | 1.95e-14 | 1.69 | **12.9%** |

   Understatement factor up to **1.8e14**. **At `h ≈ 3%` the only instrument the
   code offers reports convergence to machine precision while the sensitivity is
   wrong by 14%.** That is 09-25's finding reproduced on this estimator, at 14%
   rather than 44%. Q4's reasoning failed for the reason that made Family B the
   right place to look: `λ_imp = λ_bulk/residual` with a near-zero residual puts
   large higher derivatives into `ln|f|`, so `c₂` and `c₄` acquire opposite
   signs. **The near-cancellation that makes `κ = 19.71` is the same
   near-cancellation that creates the blind spot** — one mechanism, two
   symptoms, and Section 5.5 had measured only the first.

5. **What is NOT claimed, stated before anything else is drawn from item 4.**
   Nothing published moves. Every `S` in Sections 4.10 and 5.5 is stable to
   **1.3e-8** relative under 16× refinement, against the corrected estimator,
   and against an independent Richardson reference (**Q5 PASS**). The nearest
   zero is **≈3000× above** the shipped `h = 1e-5`. The shipped step is safe by
   a wide margin. What is unsafe is the **procedure**: a sweep over the natural
   range passes straight through all four points — and **the repo ran exactly
   such a sweep on 09-25 and 09-26**.

6. **The exactly-solvable half, which is what licenses item 4.** A blind spot
   was *constructed* rather than found: for `ln|f| = A L² + B L⁴` the secant is
   exact in closed form, and with `A = −2B + δ` at `L = 1`, `conv` vanishes near
   `h*² = 3δ/(20B)` with a relative error there of **exactly `(3/8)h*²`,
   independent of `A`, `B` and `δ`**. Measured ratio to that prediction
   **1.0668** — **Q3 FAILS on its magnitude band** (the `h⁶` term contributes
   6.7% at `h* = 0.0179`), class call correct. The construction also gives the
   honest *limit* on the mechanism: the ratio of understatement is unbounded
   while the absolute error hidden is of ordinary `O(h²)` size. Family B exceeds
   that limit only because its `h*` is 3% rather than 1.8%. **Q2 confirmed to
   four digits:** `conv/|error|` = 3.0000, 3.0001, 3.0006, 3.0054.

7. **UNPREDICTED RESULT — `conv` is minimised where the answer is worst, and
   this is more general than item 4.** Median ratio of the `h` that minimises
   the error to the `h` that minimises `conv`: **877**. Minimising `conv` gives a
   **worse** answer than the shipped step at **18 of 23** call sites, worst
   penalty **3.5e4×**. Far below the cancellation knee `conv` differences two
   noise samples and is small for that reason. Item 4's blind spots are four
   isolated points; **this is the whole low-`h` half of the range.**

8. **UNPREDICTED RESULT — for two of three families `conv` measures rounding,
   not convergence.** `conv` at the shipped step divided by the derived floor
   `4ε|ln Q|/h`: Family B median **431**, Family C median **0.4**, Family D
   median **0.5**. So Section 5.5's "worst convergence estimate over the table:
   2.7e-10" and "worst convergence estimate: 1.9e-10" lines are measurements of
   double precision and **would print roughly the same number for a model with
   any amount of curvature**. The reassurance scales with `|ln Q|` and `h`, not
   with the quality of the answer.

9. **Q7: eleven `log_sensitivity` call sites, ZERO of which compare the returned
   `conv` against anything.** Every one prints it. Class call right, **scope
   wrong** — one *unrelated* module does compare its own convergence quantity to
   a tolerance, so the repository knows how and simply did not here. A
   convergence estimate with no criterion attached is a number, not a check, and
   items 4, 7 and 8 are what that number was hiding.

10. **The fix, and what it costs — reported rather than absorbed.**
    `symmetric=True` (now the default) samples `p·exp(±h)`, exactly symmetric in
    `L`. The old path is kept as `symmetric=False`, not deleted, per 09-24.
    Before/after diff of `graphene_sensitivity_audit.py`: **every value printed
    to four decimals is identical**, sole exception one `0.0000` becoming
    `-0.0000`. But **Validation 4(c)'s "three independent exact zeros … the
    strongest form of check this repo uses" becomes 2/3 BITWISE**:
    `S(ρ_bulk)` at the calibration width is now `−1.11e-11` instead of exactly
    `0.0`. That value **is** the derived floor `4ε|ln ρ|/(2h) = 5.7e-11`, so all
    three zeros survive as mathematics and one has lost its *bitwise* exactness —
    which turns out to have depended on the **old** estimator's evaluation
    points happening to cancel. **An exact-zero check written as `v == 0.0`
    cannot distinguish structure from floating-point luck.** Both counts are now
    printed, the criterion is the derived floor, the annotation is in the source,
    and the structural fact is untouched.

11. **`validate_power_laws` annotated in place as NON-DISCRIMINATING, with its
    own output as the proof.** All four sub-checks are power laws or closed
    forms; a power law makes `ln|f|` linear in `ln p`, so every secant is exact
    and every step equally good — it passes identically for the old estimator,
    the corrected one, and an absurd `rel_step = 0.25`. **Switching to the
    provably more accurate estimator moved its "worst |error|" from 1.66e-11 to
    3.79e-11.** A validation whose number gets *worse* when the estimator gets
    *better* is measuring rounding. Retained deliberately and labelled, because
    a labelled blind check is more useful than a deleted one — continuing 09-25
    item 11 and 09-26 item 10, and now the third consecutive session with such a
    check identified.

12. **THE AUDITOR COMMITTED THE AUDITED ERROR THREE TIMES, IN ONE SESSION.**
    Three of five validations failed in their first form and **all three failed
    the same way: a round absolute tolerance — 1e-13, 1e-14, 1e-13 — that
    silently encoded a step size.** V2 required 1e-13 and measured 5.10e-08,
    which *is* the derived floor at `h = 1e-3`. V5 required 1e-13 and measured
    its own floor. V3 required 1e-14 and was wrong in **shape**: `ln|1/f| =
    −ln|f|` holds in exact arithmetic, but `1/f` is a separately rounded number,
    so `log(1/x)` is not `−log(x)` bitwise. Two boundary definitions for
    "truncation branch" were also tried and discarded before a derived one
    worked. **This is the class documented 09-25 item 10 and reproduced 09-26
    item 9 — now three consecutive sessions, and today three times inside the
    session auditing it.** 09-26 concluded that knowing the class does not
    prevent reproducing it. Today supports something stronger: **an absolute
    tolerance is the default way a numerical assertion gets written,
    exhortation does not fix it, and the only defence found so far is
    structural — report the measured value beside a derived bound and let the
    ratio be the verdict.** All three first forms are left in the source.

13. **WRITING — Chapter 7 is drafted, after seven consecutive sessions of
    displacement.** `thesis_draft/07-discussion-and-outlook.md`, 704 lines,
    eleven sections, discharging all ten synthesis threads. No new physics; it
    is the synthesis the other chapters were generating material for. Its four
    arguments: (i) the single-junction/two-junction split, with Section 6.14's
    parity making the blindness exact and the device-facing polarity statement
    carried verbatim; (ii) bounds and parities versus rankings over tabulated
    sets, with Section 6.9's retracted ceiling as the cautionary case that a
    claim's *grammar* is not evidence about its population; (iii) five failure
    classes with five **non-overlapping** detectors and a table of which ones
    Chapters 4–5 actually have; (iv) unequal scrutiny — `κ = 19.71` in Chapter 5
    against 4.55 in Chapter 6, and the only convergence blind spots in the
    repository are at that same Chapter 5 calibration, so the chapter that has
    failed most publicly is the one that has been *looked at* most. Chapter 4's
    `f_T` is the worked example: "≈20 GHz" and 10.140 GHz are both right, at
    `V_ds = 0.1 V` and `0.05 V`, and the chapter had never named which for a
    month. Chapter 1's status row rewritten to lead with the draft, with the
    ten-thread history kept below it. **Two threads added on the drafting day
    itself:** an untested prediction that any *differential* interconnect figure
    of merit is dominated by edge-scattering variance rather than mean
    resistivity, and the caution that **a bound derived from a near-cancelling
    quantity inherits the near-cancellation** — so "restate it as a bound" is
    not a way around conditioning but a way of making it explicit.

**Methodological note, continuing the series.** 09-20: exact validation does not
protect against an unrepresentative sample. 09-21: pre-registration reaches what
exact validation cannot. 09-22: pre-registration does not reach the analysis
layer. 09-23: how claims consolidate. 09-24: a default tolerance is a claim
about the scale of the caller's variable. 09-25: a procedure asked whether it has
converged can answer yes and be wrong by 44%. 09-26: the anchored comparison is
anchored in one variable; entangled defaults break it.
**09-27, and it is the sharpest form yet: an instrument can be systematically
smallest where the answer is worst. `conv` is minimised a median 877× below the
right step, and at four real points it reports 1e-14 where the answer is 14%
wrong. A convergence estimate is not a convergence criterion, and the two are
distinguished only by attaching a derived bound to it — which nobody in this
repository had done at any of eleven call sites in five weeks.** Second, and
earned three times over today: **the fault class of asserting a round absolute
tolerance is not a lapse of care but a default of authorship, and the fix has to
be structural.**

**Predictions scored:** Q1 (class and magnitude), Q2, Q5, Q7 (class) **PASS**.
Q3 **FAIL** on magnitude band (1.067 vs 1.000), class correct. Q4 **FAIL** on
its class call *and* as written. Q6 **FAIL** (Family B worse by 4.35×, not
>10×). Q7 **FAIL** as written (over-scoped). **Unpredicted:** items 7, 8, 10 and
the item-1 correction. **5/5 validations pass** after three of the five failed in
their first form and were rebuilt on derived bounds.

**Not yet covered (candidates for future runs):**
- **Whether `graphene_default_scale_audit`'s own `rel_step` sweeps are affected
  by today's item 7** — created today and the natural top numerical item, since
  that module swept this estimator's step on 09-25 and today showed that
  minimising `conv` walks 877× away from the right step. The question is whether
  any of 09-25's conclusions rest on a `conv`-guided step choice
- **Whether any OTHER near-cancellation in this repo would show the same
  `c₂`/`c₄` sign flip if it were differentiated rather than evaluated** —
  created today. Chapter 4 §4.7's residual and Chapter 6's same-sign pairs are
  the candidates; today's mechanism says near-cancellation and blind spot are
  the same mechanism, so `κ` is a *predictor* of where to look
- **Migrating every remaining validation in the repo to measured-value-beside-
  derived-bound form** — created today by item 12, and the only defence found
  against a fault class now at three consecutive sessions. This is the largest
  piece of purely mechanical numerical work outstanding
- **Whether Chapter 4 or 5 contains a RANKING that is a step artefact** —
  created 09-23/09-25, untouched today
- **`n_segments = 50` at a bias with more curvature** — created 09-26,
  untouched today
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open since
  2026-09-21, still the top *physics* item, and now **untouched for seven
  consecutive sessions**. Chapter 7 §7.9 lists it first for exactly that reason
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural weakness,
  and Chapter 7 §7.9 item 2
- **Chapters 2–3 remain undrafted** despite complete computational results —
  **now the largest remaining block of pure writing by a wide margin**, since
  Chapter 7 is drafted. Their computational results have been complete longest of
  anything in the repo
- **Whether Chapter 7's §7.3 differential-interconnect prediction holds** —
  created today, stated as a prediction in the chapter rather than a result
- **Asking which of Chapters 4 and 5's design rules could be restated as
  parities or bounds** — created 09-23; Chapter 7 §7.4.3 now names three
  concrete candidates, so the item is sharpened rather than closed
- **Whether the parity survives a photo-thermoelectric term** — created 09-23,
  and Chapter 7 §7.9 raises its stakes: if it does not, §7.2.3's polarity
  statement loses its exactness
- **Re-check whether other "for every …" claims rest on small samples** — open
  since 2026-09-20
- **Whether Chapter 4's contact-resistance results should be re-run at the
  5.4 eV crossover** — open since 2026-09-18
- **`__pycache__` is tracked in this repo** — open since 09-24; it dirtied the
  working tree again on every run today
- Reconciling Mueller *et al.*'s 0.12 eV step (arXiv:0902.1479) with the
  0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS` assumes — open since 2026-09-18;
  Chapter 7 §7.8.2 now states it as bounding everything in Chapter 6
- A photo-thermoelectric term (Kasırga review); Shimomura *et al.*'s
  comb-electrode design; integrating 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery; isolating the root cause of the
  Section 4.7 negative residual (open since 2026-08-31); Ti and Cr per-metal
  `Rc` recalibration (ResearchGate rate-limiting); a second independent
  edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

**Automation health:** Device reachable at the 04:34 UTC firing. A plain
`git clone` of both repos completed normally and fast (the 09-24 partial/shallow
recipe again not needed). `git config user.name/user.email` again absent in the
fresh clones and set in its own call per 09-24. `scipy` again absent from the
device VM at import time; `pip install scipy` succeeded immediately and **was
required today** — `graphene_sensitivity_audit` imports it transitively through
`graphene_contact_doping_model`. `__pycache__` again showed as modified on every
run. One self-inflicted syntax error (a multi-line f-string inside a patch
heredoc) and one heredoc-quoting collision (a `'''` sequence in prose inside a
`'''`-delimited Python string) cost two round-trips; both are authoring hazards
of patching Python through heredocs, not repo problems.

**Commits this run:** 5 (the pre-registration; the audit with its output and
figure; the estimator fix with its before/after regression check and the
Validation 4 annotation; the study note; Chapter 7 plus the Chapter 1 status
row). This AUTOMATION_LOG.md entry makes 6. The audit and the fix are separate
commits for the same reason as 09-26: the audit's Q5 result is that nothing
moved, and a session whose numerical headline is "nothing published moves"
should not bury a docstring-falsifying estimator bug inside that commit.

---

## 2026-09-28 — The band structure of graphene in this repository had no Dirac point, and a prediction's stated reason outlived the code it referred to

**Status:** Automated session. Live web search **not used** — both halves were
internal audits of this repository's own code and the writing half was a
synthesis of results in hand; neither needed it. Device reachable at the
**05:28 UTC** firing (a schedule of 04:33 delivered late), folder connected;
neither repo had a 2026-09-28 entry and neither had commits since midnight, so
a full session was run.

**Pre-registration committed first** (commit `6b5dc78`, before the audit module
existed): `notes/2026-09-28-covariance-probe-preregistration.md`, five
questions Q1–Q5 and five validations V1–V5, continuing the practice established
09-21. The second half of the session was **not** pre-registered, and why is
item 6.

### Half one — auditing the auditor

The logged top numerical item was whether `graphene_default_scale_audit`'s own
`rel_step` sweeps inherit 09-27 item 7. Reading the module to scope that
question turned up something sharper, so the session answered both.

1. **Q4 PASS — the logged top item closes as a NULL, and the null is
   load-bearing.** No conclusion in `graphene_default_scale_audit` rests on a
   `conv`-guided step choice. Its single `log_sensitivity` call site binds the
   second return value to `_`; the `H_DIFF` and RESULT 7 sweeps go through
   `sensitivity()`, which uses an **absolute** step and never computes `conv`.
   The census is by AST plus six regex patterns, and it is given a **positive
   control** — the same census run on a snippet that *does* minimise `conv`,
   which it detects. Without that control the null would be
   indistinguishable from a census that cannot see anything. This is
   mutation-testing's move applied to a static census. **None of 09-25's
   conclusions inherit 09-27 item 7. Item closed.**

2. **Q1 PASS — a prediction's stated reason outlived the code it referred
   to.** Prediction D2 of `graphene_default_scale_audit` reads: *"
   `log_sensitivity(rel_step=1e-5)` passes covariance EXACTLY (bitwise zero
   deviation), **because `rel_step` multiplies the parameter**."* On 09-27 the
   default path was changed so that `rel_step` **exponentiates** the parameter
   (`p·exp(±h)`). D2 still returns bitwise `0.0`, and returns it identically on
   both paths. **A test whose verdict is unchanged by the removal of its
   stated cause was never testing that cause.** D2's number stands; D2's
   explanation is falsified.

3. **UNPREDICTED, and it subsumes Q2 and Q3 — D2's check is a test of
   `(λa)/λ == a` and of nothing else.** When that round trip is exact, the
   scaled and unscaled branches are **literally the same floating-point
   expression**, so the deviation is bitwise zero for *any* function, estimator
   or step. Verified as an **exact set equality**, no tolerance anywhere: over
   λ = 1e1…1e15 the decades with non-zero deviation are exactly the decades
   with an inexact round trip — **{1e7} for the corrected estimator, {1e4} for
   the old one**, for both a power law and a log-quadratic. Different decades,
   one cause, and nothing mathematical distinguishes those decades. Three
   scored items fall out of this: **Q2 FAILS on magnitude** (non-zero at 7% of
   decades, predicted ≥ 20%; class correct); **Q3 FAILS outright** — a curved
   test function does *not* break the exactness, because the check never
   reaches the estimator; and **V5 FAILS as first written**, one level up,
   because `covariance_deviation` divides by λ too.

4. **The probe's docstring is wrong in both of its claims, measured.**
   *"Identically zero for a scale-free procedure"* — no, zero to **one ulp**:
   an exactly-proportional solver gives 1.29e-16 against a derived bound of
   `eps` = 2.22e-16, and V1's bitwise zero survives only because λ = 1 makes
   the round trip trivial. *"Grows with the mismatch"* — **only downward**.
   **Q5 FAILS: class call correct, direction and magnitude backwards.** Q5
   predicted the downward run would be the lenient one; it is the **severe**
   one, by **1.86e4×** (2.47e1 at λ = 1e-3 against 1.327e-3 at λ = 1e3). And
   the structure is not one curve: **upward it saturates** (λ = 1e3 and 1e6
   agree to 5e-5, so the probe has a **ceiling** and its magnitude cannot
   express severity) while **downward it diverges**, steepening faster than any
   power of λ — reported as a breakdown, not fitted, since two or three points
   are not a law. The convictions in RESULT 2 are **correct**; their
   magnitudes are not a severity scale.

5. **A true claim withdrawn by the wrong oracle, then recovered — and the
   failed oracle kept.** The saturated ceiling was identified with the shipped
   step's own truncation error. A Richardson reference from `(h, h/2)` implied
   1.53e-2, missed the saturated 1.327e-3 by **11.6×**, and the module printed
   *"the identification DOES NOT HOLD, and the claim is withdrawn."* Direct
   refinement to `h = 1e-6` gives 1.3288e-3 — the same number to **0.14%**.
   Richardson is the invalid one, and the refinement sequence proves it: the
   error **grows** from `h = 1e-3` (1.329e-3) to `h = 3e-4` (4.705e-3), so
   there is no smooth `h²` expansion to extrapolate. **Mirror image of 09-25
   and 09-27:** there an instrument reported success where the answer was
   wrong; here one reported failure where the claim was right. Both first
   forms stay in the source.

   Four annotations went into `graphene_default_scale_audit.py`, verified
   comment-only (full stdout byte-identical before and after), with every
   original wording retained above its annotation.

### Half two — the band structure, found by writing rather than by auditing

6. **THE HEADLINE, and it was not the planned work. The shipped band structure
   of graphene had NO DIRAC POINT.** `k_path_graphene()` wrote, three lines
   apart, `k_point = [4/(3√3), 0]` — **no π** — and `m_point = [π/(3√3), π/3]`
   — **with π** — while `graphene_hamiltonian()` builds its phases from
   nearest-neighbour vectors of unit length, for which the corner is at
   `4π/(3√3)`. The shipped K was **exactly a factor of π too small** (measured
   ratio 3.141592654), so `|φ(K)| = 2.5718` instead of 0:

   | quantity | as shipped | corrected | exact |
   |---|---|---|---|
   | gap at the point labelled K | **14.4019 eV** | 2.5e-15 eV | 0 |
   | minimum gap over the whole path | **11.2000 eV** | 2.5e-15 eV | 0 |
   | gap at the point labelled M | 11.2000 eV | 5.6000 eV | `2t` = 5.6 |

   Graphene's single defining electronic property was absent from
   `band_structure.png`, the repository's oldest figure, for **five weeks of
   daily automated sessions**. Superseded values kept as
   `HIGH_SYMMETRY_LEGACY`, broken spectrum still reproducible via
   `calculate_band_structure(legacy=True)`.

   **It was found by drafting Chapter 2**, not by any audit, validation,
   convergence study or pre-registration — which is why this half of the
   session has no pre-registration and should not have one. A chapter cannot
   contain the sentence *"the two bands touch at K"* without someone checking
   whether they do.

7. **Why five weeks of audits never touched it: three band/DOS figures exist
   here and NOT ONE was a calculation that could have failed.**
   `band_structure_simple.png` *draws* `E = ±v_F|k|`, the known answer.
   `calculate_density_of_states()` returns `|E|/(πt²)` analytically and **never
   calls the Hamiltonian at all** — demonstrated by **mutation**: replacing
   `graphene_hamiltonian` with the k-independent gapped matrix `diag(+7, −7)`,
   which has no Dirac point and no dispersion, leaves its output **bitwise
   identical**. The third computed, and was wrong. A DOS computed from the
   bands would have caught this on the first run. **The figure that could have
   falsified the calculation was instead the figure that certified it.** Both
   non-discriminating functions annotated in place and retained, per 09-25
   item 11 and 09-26 item 10.

8. **Five exact validations, two of which FAILED in their first form — both
   failures kept and explained.**
   **V2 (`v_F`) failed because of PHYSICS.** It assumed second-order
   convergence of the one-sided slope and failed by 2e4; measured errors were
   2.512e-3 / 2.501e-4 / 2.500e-5 at `dk` = 1e-2 / 1e-3 / 1e-4, clean **first**
   order. Derived along Γ–K: `|φ| = (3/2)q(1 − q/4 + O(q²))`, so the leading
   correction is linear in `q` with coefficient **exactly 1/4** — trigonal
   warping. Rewritten to check the derived constant: measured **0.250013**, and
   correcting by `(1 + q/4)` recovers `v_F = 9.062708e5 m/s` to **1.9e-9**.
   **V4 failed twice, and one failure was THE AUDITED FAULT COMMITTED BY THE
   AUDITOR, WITHIN THE HOUR.** The reciprocal vectors `(4π/3)(1,0)` and
   `(4π/3)(½,√3/2)` have the right magnitude and span the right area but are
   **not reciprocal lattice vectors** — `|φ|` is not periodic under them — so
   the sampled cell held **one** Dirac cone instead of two and halved the
   low-energy DOS. And the obvious check did not notice: *"the DOS integrates
   to 4 states per unit cell"* passed to **sixteen digits** with the wrong
   vectors in place, because that normalisation divides by the sample count and
   returns 4 per cell whatever the Hamiltonian is. **The same
   cannot-fail fault as the shipped DOS figure, inside the module written to
   audit it.** Kept and printed as **vacuous**, beside the two checks that
   replaced it: periodicity of `|φ|` under `b₁,b₂` (1.6e-15 vs a derived
   1.4e-14) and **a count of Dirac cones in the sampled cell** (2 required, 1
   found). Corrected, the computed DOS matches the exact slope
   `2q_e²/(πħ²v_F²) = 1.7891e18 eV⁻¹m⁻²` to **3.6%** (residual = the same `q/4`
   term, derived window bound 19.6%) and the van Hove peak sits at **2.771 eV**
   against `t` = 2.800. V1 (`|φ|` = 3, 1, 0 — three integers, three chances to
   fail), V3 (particle–hole symmetry **bitwise at 4000/4000** random k) and V5
   (the diagnosis reproduces the 14.401937320702 eV gap **bitwise** as
   `2t|φ(K_shipped)|`) passed first time.

9. **NEW PHYSICS, small but missing: the first quantitative validity window for
   the linear dispersion this thesis has ever had.** Carrier density from the
   full nearest-neighbour bands against `n = E_F²/(πħ²v_F²)`: ratio 0.9954 at
   0.1 eV, 1.0021 at 0.3, 1.0051 at 0.5, 1.0220 at 1.0, 1.0539 at 1.5. **1%
   over 0.1–0.5 eV, 5% to 1.0 eV.** Stated as a **window** not a limit, because
   the deviation is **non-monotone** — 1.0112 at 0.05 eV is worse than at 0.5 —
   and that is k-sampling (at 0.05 eV the cones cover 5.9e-5 of the cell, ~190
   samples on an 1800² grid), with the sample count printed so the artefact is
   legible. The **sign** is what the device chapters carry: the full bands hold
   *more* states than the cone above ~0.25 eV, so a linear-DOS model
   **understates** carrier density at high bias — conservative for drive
   current, against itself for quantum capacitance.

10. **NOTHING IN CHAPTERS 4–7 MOVES, and the reason is the finding.** Those
    chapters import `v_F`, the linear DOS and the absence of a gap, and import
    the k-path and the 2 × 2 `H(k)` **nowhere** (Chapter 2 §2.11 tabulates
    this). That is simultaneously why a 14 eV error was survivable for five
    weeks and why no device number changes now. `v_F` confirmed to 1e-6.

11. **WRITING — Chapter 2 is drafted**, `thesis_draft/02-electronic-structure.md`,
    459 lines, eleven sections, and the drafting **falsified its own status
    row**: "Computational results complete", in place since 2026-08-23. §2.8's
    quantum capacitance (1.43–14.33 µF/cm² over 0.05–0.5 eV) is the one place
    Chapter 2's physics is directly visible in Chapter 4's numbers, being
    comparable to a 1 nm-EOT `C_ox` rather than negligible beside it. §2.11 is
    a table of exactly what later chapters import, whose last row reads "the
    k-path and 2 × 2 `H(k)`: nothing". Chapter 1's status row rewritten to lead
    with the correction, ten-thread history kept below.

12. **`__pycache__` untracked and a `.gitignore` added — an item open since
    09-24, and it had teeth.** Five sessions of "shows as modified on every
    run" was not cosmetic: today it broke a `git stash` round trip mid-session,
    because the caches were rewritten while the stash was applied. That is a
    failure mode that loses work rather than adding noise.

**Methodological note, continuing the series.** 09-20: exact validation does
not protect against an unrepresentative sample. 09-21: pre-registration reaches
what exact validation cannot. 09-22: it does not reach the analysis layer.
09-23: how claims consolidate. 09-24: a default tolerance is a claim about the
scale of the caller's variable. 09-25: a procedure asked whether it has
converged can answer yes and be 44% wrong. 09-26: an anchored comparison is
anchored in one variable. 09-27: an instrument can be systematically smallest
where the answer is worst.
**09-28 gives two, and the first is the most transferable thing five weeks have
produced. (i) PROSE IS A DETECTOR, and it is the one this repository did not
have. A 14 eV error in its foundational figure was found by writing the
chapter, after five weeks of daily audits, exact validations, convergence
studies and pre-registrations had all passed over it — because those instruments
check what the code says about itself, and a chapter has to state what the code
says about GRAPHENE. The sessions that displaced Chapter 7 seven consecutive
times, and Chapters 2–3 for five weeks, were not merely accumulating a writing
backlog; they were accumulating an unmeasured error. (ii) A CHECK THAT CANNOT
FAIL IS WORSE THAN NO CHECK, because it certifies whatever sits beside it — and
the defence is not care but a deliberate provocation: mutate the Hamiltonian and
require the DOS to change; break the reciprocal vectors and require the cone
count to drop; hand the census a snippet containing the pattern it reports
absent. Every non-discriminating check found today — the shipped DOS figure,
`_logsens_invariance`, the 4-states-per-cell normalisation written this morning —
was found by that move and by no other. Three of them, one written today, in a
repository that has been auditing itself daily for five weeks.**
Secondary, and now twice in one day: **an instrument reporting FAILURE misleads
exactly as one reporting success does** — a Richardson oracle withdrew a true
claim, and a convergence-rate assumption mistook trigonal warping for a bug.
Both were settled by deriving the expected behaviour in closed form and testing
that, not the rate.

**Predictions scored:** **Q1 PASS** (D2's reason falsified by its own continued
success). **Q4 PASS** (null, with a positive control). **Q2 FAIL** on magnitude
(7% vs ≥ 20%), class correct. **Q3 FAIL** outright — and the reason it cannot
hold (RESULT 6) is stronger than the prediction was; the pre-registration
reached for the newest available criticism (09-27's power law, one day old)
rather than the sharpest one, which is recorded against this session.
**Q5 FAIL** — class call correct, direction and magnitude backwards (1.86e4×
more severe downward, not less). **Unpredicted:** RESULT 6, the saturation/
divergence structure, the Richardson reversal, and all of half two.
**Validations:** half one 5/5 against derived bounds, with V5's as-first-written
bitwise form failing and kept; half two 5/5, with V2 and V4 failing in their
first forms and both kept.

**Not yet covered (candidates for future runs):**
- **Whether any DEVICE figure in this repository is drawn from a known answer
  rather than computed** — created today, and **the new top item in the
  repository**. Three of three band-structure figures were. Chapters 4–6's
  figures have never been asked. The detector is known and cheap: mutate the
  model the figure claims to come from and require the figure to change
- **Chapter 3 is still undrafted, and after today it is a RISK item rather than
  a writing-backlog item** — it carries the same "computational results
  complete" label Chapter 2's did this morning, and that label has just been
  shown to mean nothing. Drafting it is now the highest-value work available on
  both the writing and the correctness axis
- **Whether other `== 0.0` exactness checks in this repo are round-trip
  tautologies** — created today by RESULT 6; the mechanism generalises to any
  check comparing a quantity against a rescaled recomputation of itself
- **Whether RESULT 2's other convictions are read from the saturated branch** —
  created today. Only the `tol=1e-14` conviction has an independent anchor
  (09-24)
- **A probe that reports both directions of λ** — created today; its whole
  dynamic range is in the direction the module never runs, and today's table is
  most of the work
- **Angular trigonal warping** (the `cos 3θ` structure) — created today; §2.5
  derives the coefficient along Γ–K only, and the angular part is what matters
  for anisotropic transport
- **Finite-temperature carrier density `n(E_F, T)`** — created today; §2.7 is
  the `T = 0` integral and Chapter 6 runs at room temperature
- **Whether `t′` can be EXCLUDED quantitatively** as the source of Chapter 4's
  electron–hole asymmetry, rather than argued away on magnitude — created today
- **Whether any OTHER near-cancellation would show the `c₂`/`c₄` sign flip if
  differentiated rather than evaluated** — created 09-27, untouched
- **Migrating every remaining validation to measured-value-beside-derived-bound
  form** — created 09-27 item 12, and today added two more instances of the
  fault class caught by that defence
- **Whether Chapter 4 or 5 contains a RANKING that is a step artefact** —
  created 09-23/09-25, untouched
- **`n_segments = 50` at a bias with more curvature** — created 09-26, untouched
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open since
  2026-09-21, still the top *physics* item, now **untouched for eight
  consecutive sessions**
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural weakness
- **Whether Chapter 7's §7.3 differential-interconnect prediction holds** —
  created 09-27
- **Which of Chapters 4 and 5's design rules could be restated as parities or
  bounds** — created 09-23, three candidates named in Chapter 7 §7.4.3
- **Whether the parity survives a photo-thermoelectric term** — created 09-23
- **Re-check whether other "for every …" claims rest on small samples** — open
  since 2026-09-20
- **Whether Chapter 4's contact-resistance results should be re-run at the
  5.4 eV crossover** — open since 2026-09-18
- Reconciling Mueller *et al.*'s 0.12 eV step (arXiv:0902.1479) with the
  0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS` assumes — open since 2026-09-18
- A photo-thermoelectric term (Kasırga review); Shimomura *et al.*'s
  comb-electrode design; integrating 6.5's plasmonic near-field picture with
  the spatially-resolved contact-doping machinery; isolating the root cause of
  the Section 4.7 negative residual (open since 2026-08-31); Ti and Cr
  per-metal `Rc` recalibration (ResearchGate rate-limiting); a second
  independent edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

**Automation health:** Device reachable at the **05:28 UTC** firing (scheduled
04:33, delivered late — the first of the day's three, so the redundancy was not
needed). Plain `git clone` of both repos completed normally and fast.
`git config user.name/user.email` again absent in the fresh clones and set in
its own call per 09-24. `scipy` again absent from the device VM and
`pip install scipy` succeeded immediately; it is required transitively by
`graphene_covariance_probe_audit` through `graphene_default_scale_audit` →
`graphene_contact_doping_model`. **New and worth recording: the documented
`GIT_ASKPASS` recipe failed with "Permission denied" on `/tmp/.tok` and
`/tmp/askpass.sh`, and was moved to `$HOME/.sess/`** (still outside `mnt/`,
still invisible to the user, deleted with the session).
**CORRECTION, appended later the same session after the other repository's half
diagnosed it properly — the first sentence written here said `/tmp` is not
writable, and that is wrong.** `/tmp` in the device VM *is* writable. What is
true is sharper and more consequential: **`/tmp` PERSISTS ACROSS SESSIONS, and
files left there by a previous session are owned by a different uid with mode
600**, so with `/tmp` sticky they cannot be overwritten, read or removed.
Yesterday's `/tmp/.tok` and `/tmp/askpass.sh` are still present, which is the
actual cause of the "Permission denied", and the 09-25 `iverilog` install at
`/tmp/iverilog_install` is still present for the same reason. The stale
`/tmp/.tok` is 93 bytes dated 2026-09-27 04:58 — the same length as the current
token — so **a plaintext copy of the PAT has outlived its session, the recipe's
`shred -u` did not take effect, and this session cannot remove it**; reported to
Harsh by notification, with a recommendation to rotate the token and move the
recipe to `$HOME/.sess/` permanently. The wrong first diagnosis is left above
rather than edited out, per this repository's own annotate-don't-rewrite rule,
because "the directory is unwritable" and "yesterday's files are still there and
belong to someone else" call for different fixes and only the second is right. Everything else in the push recipe worked unchanged
and every commit was pushed as it was made, per 09-25. **`__pycache__` cost a
round trip today and is now untracked** (item 12). Two self-inflicted
`print`-formatting errors (a `%`-format argument attached to the wrong `print`
of a multi-line block, twice) cost two round trips; the same authoring hazard
class as 09-27's heredoc collisions.

**Commits this run:** 8 (the pre-registration; the covariance-probe audit with
its output and figure; the four in-place annotations; the `__pycache__`
untracking; that half's study note; the band-structure fix with its audit,
output and regenerated figures; the linear-cone validity window; Chapter 2 plus
the Chapter 1 status row; that half's study note). This AUTOMATION_LOG.md entry
makes 9. The band-structure fix is separate from the validity-window result
because the first is a **correction to a published figure** and the second is a
**new measurement**, and a correction that arrives buried inside a new result is
one nobody reads — the same reasoning as 09-26 and 09-27.

---

## 2026-09-29 (automated session)

**Two halves, and the second replicated the first day's finding on a different
chapter.** Half one answered the repository's top open item; half two drafted
Chapter 3 and, in drafting it, falsified its status row for the second
consecutive chapter on the second consecutive day.

### Half one -- the top open item, answered

1. **`graphene_figure_provenance_audit.py`, `figure_provenance_audit_output.txt`,
   `figure_provenance_audit.png` -- NO device figure in Chapters 4-6 is drawn
   from a known answer.** The item created 09-28 and made the repository's top
   item is closed. Six device figures mutated against the models they claim to
   come from: `quantum_capacitance` responds to `v_F` and `T`;
   `gfet_transfer_characteristics` to `C_ox`, mobility and the puddle floor;
   `rf_figures_of_merit` to channel length, mobility and finger count;
   `edge_vs_top_contact` to the Au Fermi shift by 519x. **Four null controls
   held at exactly zero** (mobility and contact resistance mutated 10x against
   `C_q`, which is electrostatic and cannot contain them), so the detector is
   not reporting that everything moves.
   **[ANNOTATED 2026-09-30 -- this inference does not follow, and the numbers
   above are kept unchanged.** A `MUST_NOT_CHANGE` check passes whenever the
   mutation fails to arrive, even when the output is 100% dependent on the
   mutated constant, and its zero is bitwise the zero of a correct null
   control -- proved in `graphene_mutation_arrival_probe.py` SECTION 1. These
   four zeros therefore carry no evidence about the detector's specificity;
   that evidence came entirely from the `MUST_CHANGE` checks that passed in the
   same run. Three of the six null controls in the table are additionally
   identity re-assignments, and one of them, `ALPHA_ABS`, was run against a
   name read nowhere in its module. See the 2026-09-30 entry, items 2, 5 and
   6.]** **Four closed-form validations, three
   bitwise:** `R_bare` exactly linear in EQE (ratio 2.0, deviation 0.000e+00
   over 200 decades); `G*BW` exactly constant (spread 3.8e-16); `C_q` exactly
   `1/v_F^2` (deviation 0.000e+00 over 401 gate voltages); `C_q` exactly even
   in `V` (0.000e+00). **The band-structure fault is genuinely absent from the
   device chapters**, which was not the prior after 09-28 and is worth stating
   plainly.

2. **What the audit found instead -- THE WRITE-ONLY KNOB, a fault class new to
   this log.** Three MUST_CHANGE checks failed bitwise, and all three were one
   fault: a module-level constant captured as a DEFAULT ARGUMENT when `def`
   executed. Rebinding the module attribute rebinds the *name* and leaves the
   captured default alone, so the constant at the top of the file reads as a
   parameter and behaves as a comment. **This matters here more than it would
   elsewhere because this repository audits BY MUTATION:** a mutation applied
   the obvious way -- `setattr(module, CONST, value)` -- is silently a no-op,
   and the audit reports "no sensitivity" when it has measured "no mutation".
   **Census: 49 occurrences across 14 files, 14 of them inside audit modules**,
   with a positive control on the census itself per the 09-28 rule.

3. **Two instances worse than inert.** (i) `plot_resistivity_vs_linewidth`
   builds its legend from the LIVE global and its curve from the FROZEN
   default, so with the constant at 4 nm it draws the 3 nm curve (4.8980
   uohm.cm at W = 20 nm) under a legend reading "Cu, 4nm TaN/Co liner" -- the
   figure is wrong **in writing, on the figure**, by 42.9% in the plotted
   quantity, undetectable from the figure alone. (ii)
   `graphene_band_structure_audit.T_HOP` is captured by the MEASURING functions
   (`bands`, `fermi_velocity_analytic`, `dos_from_bands`) and LIVE in that same
   module's expected-value expressions, so doubling it moves the ORACLE
   (`2*T_HOP`: 5.600 -> 11.200 eV) and freezes the MEASUREMENT (8.379013 eV,
   bitwise unchanged). The audit would report a disagreement manufactured by
   its own mutation machinery, and an investigator following it would be
   debugging a Hamiltonian that never changed. **Yesterday's Hamiltonian
   mutation escaped this only because it replaced a FUNCTION rather than a
   constant -- luck, not design.**

4. **Nine sites fixed, under a bitwise no-change requirement.** Converted to
   late-bound defaults in `graphene_interconnect_model` (3 functions),
   `graphene_photodetector_model` (4) and `graphene_edge_contact_model` (2).
   The hazard in such a fix is quietly moving a published number, so it was
   made subject to an exact requirement rather than a plausibility check: all
   three figures (86400, 8112 and 14592 bytes of harvested curve data) and all
   three `summary_numbers()` transcripts captured before and after and compared
   exactly. **None moved.** Section A now reports 24/24 with the same mutations
   that failed an hour earlier, and Section B asserts the two properties rather
   than narrating them. **37 occurrences deliberately left, 13 in audit
   modules:** fixing an audit module from inside the audit that found it moves
   the instrument and the measurement in one step.

### Half two -- Chapter 3, and the same falsification again

5. **`graphene_transport_optical_audit.py`,
   `transport_optical_audit_output.txt` -- FOUR defects in
   `graphene_transport_properties.py`.** Chapter 3 carried "Computational
   results complete" from 2026-08-23, the same label and the same day as
   Chapter 2's. Every defect was caught by comparison against an exactly known
   value and none by a plausible range.
   - **T1.** Universal optical conductivity coded `pi*e^2/(2*hbar)` while its
     own docstring states `pi*e^2/(2*h)`. A factor 2*pi: 3.823530e-04 S against
     6.085337e-05 S, implying a single-layer absorption of **14.4044%** instead
     of 2.2925%. *(The 14.40 here and the 14.40 eV spurious gap of 09-28 are a
     numerical COINCIDENCE -- that was a reciprocal-lattice convention, this is
     `hbar` for `h`. Recorded so no future reader chases it.)*
   - **T2.** Conductance quantum coded `2e^2/hbar` (4.868270e-04 S against
     7.748092e-05 S, a value fixed by metrology). Minimum conductivity coded
     `4e^2/hbar` and commented "a hallmark of graphene"; the hallmark is
     `4e^2/(pi h)` and the shipped value was **19.74x** it -- wrong by 2*pi AND
     missing the ballistic `1/pi`.
   - **T3, the sharpest.** The Pauli-blocking switch compared `hbar_omega` in
     JOULES against a bare `0.2` intended as eV, so the branch required a photon
     energy above 1.25e18 eV. **Swept from 0.1 nm to 100 um the function
     returned exactly ONE distinct value** -- a constant `0.5*sigma_0`. The
     feature its comment describes had never operated, in either direction,
     since it was written. This is the 09-28 cannot-fail check in its purest
     form: not a check that always passes but **a branch that is never
     reached**, so the code appeared to model the physics while modelling
     nothing. Section 3.6's table is the first time this repository has actually
     computed it. Positive control included: written correctly the same sweep
     returns two values.
   - **T4.** Acoustic-phonon mobility used `T^-3/2`, the THREE-DIMENSIONAL
     deformation-potential exponent. Graphene's `rho_LA` is linear in T and
     density-independent, so `mu ~ 1/T` [Hwang & Das Sarma, PRB 77, 115449
     (2008)]. Measured exponent of the shipped curve 1.5000; understates the
     500 K mobility by 1.291x -- small enough to look plausible, which is why it
     survived.

6. **All four corrected in place**, superseded expressions and their numbers
   kept in comments beside the corrections. `G_min_experimental` (~`4e^2/h`) is
   now RETURNED ALONGSIDE the theoretical value rather than the two being
   conflated, because the gap between them is a real unresolved feature of
   graphene transport. The audit now reads values FROM the module instead of
   re-typing them, so its seven checks are regression guards rather than a
   transcript of what the auditor believed the code said -- the 09-28
   cannot-fail rule applied to the auditor. 7/7.

7. **One check was itself wrong, and is recorded rather than quietly retuned.**
   The `sigma_0/(eps_0 c) == pi*alpha` anchor first ran at `tol=1e-12` and
   FAILED at 4.5e-12. That residual is the slack between CODATA's independently
   MEASURED `alpha` and the `alpha` implied by its own `e`, `h`, `eps_0` and
   `c`, well inside the ~1.6e-10 relative uncertainty CODATA quotes. The fix
   was to check the algebraic identity (3.0e-16) and name the CODATA slack
   separately -- **not to loosen a number until the check passed**, which is the
   move this repository has to be most careful about given how many of its
   tolerances are hand-set.

8. **`thesis_draft/03-transport-and-optical-properties.md` drafted, ten
   sections, 355 lines.** Built around what Chapters 4-7 actually take from it
   (Section 3.8 tabulates it). **Section 3.2:** graphene's room-temperature
   mobility is set by its ENVIRONMENT, not by graphene -- the intrinsic phonon
   limit is ~25x Chapter 4's 4000 cm2/Vs, making the substrate the single
   largest engineering lever in the thesis. **Section 3.3:** theory
   (`4e^2/(pi h)` = 20.27 kOhm/sq) and experiment (`~4e^2/h` = 6.45 kOhm/sq)
   separated rather than conflated; this is the thesis's central NEGATIVE
   result, and Chapter 4's single-digit on/off ratio and Chapter 7's "no
   digital logic" verdict both descend from it. **Section 3.4:** the first
   Landau gap at 1 T is 32.9 meV against `kT` = 25.9 meV at 300 K, which is why
   graphene's QHE is a room-temperature effect. **Section 3.5:** universal
   absorption anchored two ways to 3.0e-16 and read as a BUDGET -- 97.7% of the
   light is lost, and every Chapter 6 architecture is an attempt to buy it back.
   **Section 3.6, genuinely new physical content:** the gate-tunable Pauli edge
   `hbar*omega = 2 E_F`, computed here for the first time in this repository --
   6199 nm at `E_F` = 0.1 eV down to **1550 nm at `E_F` = 0.4 eV**, so the
   telecom C-band lands inside Chapter 4's existing gate range. That is the
   physical basis of graphene electro-absorption modulators, and T3 means the
   code that claimed to compute it never had.

9. **Blast radius checked rather than assumed, and it is the Chapter 2 pattern
   again.** None of the module's figures are committed and nothing imports it
   -- but `graphene_photodetector_model.py`'s header CITES
   `calculate_optical_conductivity` as the source of its 2.3% absorption while
   independently re-typing the correct `pi/137.036`. **Chapter 6 was protected
   from a 6.28x error in its first-line input only because it did not use the
   result it cites.** Unlike Chapter 2, though, the protection was weaker: three
   of Chapter 3's four defects were in quantities the device chapters DO use
   (Section 3.8's first five rows), and they held only by re-typing.

10. **Chapter 1's status row rewritten** to lead with the correction, per the
    09-28 precedent.

**Methodological note, continuing the series.** 09-25: a procedure asked whether
it has converged can answer yes and be 44% wrong. 09-26: an anchored comparison
is anchored in one variable. 09-27: an instrument can be systematically smallest
where the answer is worst. 09-28: prose is a detector, and a check that cannot
fail is worse than no check.
**09-29 gives two. (i) A MUTATION THAT DOES NOT ARRIVE IS INDISTINGUISHABLE, IN
THE OUTPUT, FROM A SYSTEM THAT DOES NOT RESPOND.** "No sensitivity" and "no
stimulus" produce the same number, and every mutation-based instrument in this
repository reports that number identically. The defence is not care: a mutation
harness must carry a **positive control on the mutation itself** -- a paired
assertion that some quantity the mutation must reach did in fact move -- before
it is entitled to interpret a zero anywhere else. The null controls in today's
audit were built to stop the detector over-reporting; the missing control was
the opposite one, and it is the one the three failures needed.
**(ii) PROSE IS A DETECTOR -- REPLICATED, on a different chapter, one day
later.** 09-28 proposed it from a single instance, which is exactly the
unrepresentative-sample failure this repository has recorded four times. It now
has two, and the mechanism is identical in both: **an audit checks what the code
says about itself; a chapter has to state what the code says about GRAPHENE.**
Writing `sigma_0 = pi e^2/(2h)` in a sentence and then reading the line that
computes it is a comparison no test here was making, and neither was writing
"the interband transition switches on above 2 E_F" and then reading a `np.where`
that compares joules to a bare 0.2. Two consecutive chapters, two consecutive
days, both labelled "computational results complete" since 2026-08-23, both
falsified by the act of being written.
**Corollary that follows from (ii) and should shape the next weeks:** the only
remaining undrafted chapter is gone. Chapters 1-7 are all drafted as of today,
so the detector that found both of these has no further unused fuel, and the
next equivalent instrument has to be built rather than written.

**Validations:** half one 24/24 (4 closed-form, 4 null controls, 14 mutation
responses, 2 regression guards), after 3 failures that were real and are fixed.
Half two 7/7, after 5 failures that were real and are fixed, plus one check of
the auditor's own that was wrong (the CODATA tolerance) and is recorded as
wrong.

**Not yet covered (candidates for future runs):**
- **The 37 remaining write-only knobs, 13 of them in audit modules** -- created
  today and **the new top item**. `graphene_band_structure_audit.T_HOP` is the
  worst: captured by the measuring functions and live in the expected-value
  expressions, so a mutation of it moves the oracle and freezes the
  measurement. Fix the audit modules FIRST and from outside, then re-run every
  mutation-based result in the repository, because any of them could have been
  measuring nothing
- **A positive control on the mutation itself, retrofitted to every
  mutation-based audit here** -- created today by the methodological note; this
  is the generalisable defence and it is currently in none of them
- **Whether Chapter 4's `n_puddle` actually reproduces Section 3.3's measured
  6.45 kOhm/sq floor** -- created today, a one-line calculation, and the only
  place Chapter 3's central negative result touches Chapter 4's numbers
- **The Section 3.6 Pauli edge is absent from Chapter 6's model** -- created
  today, and it is a THIRD leg of the Chapter 4 / Chapter 6 contact-metal
  contradiction Chapter 7 already holds. Unlike the other two it is not a
  modelling disagreement but a physical mechanism: biasing to high `E_F` for
  low contact resistance switches off the detector's own absorption
- **A finite-temperature optical conductivity** -- created today; Section 3.6's
  edge is the `T = 0` step function, smeared over several `kT` at 300 K, which
  is what sets a modulator's extinction ratio
- **The remote-polar-phonon cap on SiO2 is asserted, not computed** -- created
  today; the third row of Section 3.2's table has no computational backing here
- **Angular trigonal warping** (the `cos 3theta` structure) -- created 09-28,
  untouched
- **Finite-temperature carrier density `n(E_F, T)`** -- created 09-28, and
  Section 3.6 now inherits the same `T = 0` assumption, so it is wanted twice
- **Whether other `== 0.0` exactness checks here are round-trip tautologies** --
  created 09-28, untouched
- **Whether RESULT 2's other convictions are read from the saturated branch** --
  created 09-28, untouched
- **A probe that reports both directions of lambda** -- created 09-28, untouched
- **Whether `t'` can be EXCLUDED quantitatively** as the source of Chapter 4's
  electron-hole asymmetry -- created 09-28, untouched
- **Whether any OTHER near-cancellation would show the `c2`/`c4` sign flip if
  differentiated rather than evaluated** -- created 09-27, untouched
- **Migrating every remaining validation to measured-value-beside-derived-bound
  form** -- created 09-27, and today's CODATA-tolerance incident is another
  instance of the class it defends against
- **Whether Chapter 4 or 5 contains a RANKING that is a step artefact** --
  created 09-23/09-25, untouched
- **`n_segments = 50` at a bias with more curvature** -- created 09-26,
  untouched
- **A second anchor for `Delta_c`, at any separation other than 3.3 A** -- open
  since 2026-09-21, still the top *physics* item, now **untouched for nine
  consecutive sessions**
- **A description of Ti, Ni and Pd that does not go through work function** --
  three independent failures on record; Chapter 6's central structural weakness
- **Whether Chapter 7's Section 7.3 differential-interconnect prediction
  holds** -- created 09-27
- **Which of Chapters 4 and 5's design rules could be restated as parities or
  bounds** -- created 09-23
- **Whether the parity survives a photo-thermoelectric term** -- created 09-23
- Re-check whether other "for every ..." claims rest on small samples -- open
  since 2026-09-20; whether Chapter 4's contact-resistance results should be
  re-run at the 5.4 eV crossover -- open since 2026-09-18; reconciling Mueller
  *et al.*'s 0.12 eV step (arXiv:0902.1479) with the 0.25-1.07 eV offsets
  `METAL_WORK_FUNCTIONS` assumes -- open since 2026-09-18; a
  photo-thermoelectric term (Kasirga review); Shimomura *et al.*'s
  comb-electrode design; integrating 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery; isolating the root cause of the
  Section 4.7 negative residual (open since 2026-08-31); Ti and Cr per-metal
  `Rc` recalibration (ResearchGate rate-limiting); a second independent
  edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

**Automation health:** Device reachable at the **04:34 UTC** firing, the first
of the day's three, so the redundancy was again not needed. Step 0's
already-ran check was clean: neither repo had a 2026-09-29 AUTOMATION_LOG entry
or a commit since midnight. Plain `git clone` of both repos completed normally.
`git config user.name/user.email` again absent in the fresh clones and set
per 09-24. `scipy` again absent and `pip install scipy` succeeded immediately
(1.15.3 against numpy 2.2.6). **The 09-28 correction worked:** the `GIT_ASKPASS`
recipe was placed in `$HOME/.sess/` from the start and had no permission
trouble, confirming that entry's revised diagnosis (stale root-owned `/tmp`
files from previous sessions, not an unwritable `/tmp`). Live web search was
not used this session; both halves were audits and drafting against material
already in the repository, and the one external citation added (Hwang & Das
Sarma, PRB 77, 115449 (2008)) is from established knowledge rather than a fetch.
**One authoring hazard cost two round trips and is worth recording as a class:**
a Python patch script that wrote the file *before* calling `ast.parse` on it
corrupted an untracked file beyond `git checkout`'s reach, and the second
attempt then spliced into the corrupted result. The same shape as 09-27's
heredoc collisions and 09-28's `print`-format errors -- **validate before
writing, not after** -- and the specific aggravating factor is that an untracked
file has no restore path, so the ordering matters more for a new file than for
an edit to a tracked one. Every commit was pushed as it was made, per 09-25.

**Commits this run:** 6 (the provenance audit with its output and figure; the
nine-site late-binding fix; the Chapter 3 defect audit; the four in-place
corrections; Chapter 3 plus the Chapter 1 status row; half one's study note).
This AUTOMATION_LOG.md entry makes 7.

## 2026-09-30 — The control was in the wrong place: a null control cannot certify its own mutation, and 15 of the 20 remaining knobs were in T_HOP's class

**The top open item of 09-29 is closed, and closing it falsified the
methodological note that created it.** Yesterday's item asked for the 37
remaining write-only knobs to be fixed, audit modules first and from outside,
and for "a positive control on the mutation itself" to be retrofitted to every
mutation-based audit. The first half is done. The second half turned out to be
unbuildable as specified, and the reason is more useful than the retrofit would
have been.

### The instrument, and why it is not an audit

1. **`graphene_mutation_arrival_probe.py`, `mutation_arrival_probe_output.txt`
   — 29/29, of which 8 are positive controls on the probe itself.** Not an
   audit module: it imports none, and nothing imports it. That satisfies the
   09-29 constraint that fixing an audit from inside the audit that found it
   moves the instrument and the measurement in one step.

2. **THE 09-29 DEFENCE IS UNSOUND, and this is proved rather than argued.**
   Yesterday: "a mutation harness must carry a positive control on the mutation
   itself — a paired assertion that some quantity the mutation must reach did
   in fact move." Two independent failures.
   - **It cannot be applied to a null control at all.** A `MUST_NOT_CHANGE`
     check asserts the output does *not* move, so there is by construction no
     quantity in it in which to require movement. The rule therefore exempts
     exactly the checks whose whole evidential content is a zero.
   - **Even for a `MUST_CHANGE` check it confounds two causes** — "did not
     arrive" and "genuinely insensitive" give the same zero, which is the
     09-29 finding restated as its own remedy.
   `null_controls_carry_no_arrival_evidence()` runs four cases against one
   synthetic model. The fourth is the result: **a null control on an output
   that is 100% dependent on the mutated constant PASSES when the mutation is
   inert**, and its zero is bitwise the zero of a correct null control.
   **Consequence for yesterday's entry:** its claim that "four null controls
   held at exactly zero, so the detector is not reporting that everything
   moves" is not supported by those four zeros. Whatever specificity the 09-29
   audit demonstrated came entirely from the `MUST_CHANGE` checks passing in
   the same run; the null controls were riding on them. Yesterday's §1 is
   annotated in place rather than rewritten, and its numbers are kept.

3. **Where the control belongs: at the point of DELIVERY, not downstream.**
   Whether the binding the callee will read has changed is decidable by
   introspection, with no model evaluation — which is exactly why it works
   where the downstream form fails. It is available for null controls, and for
   quantities the model is genuinely insensitive to. Six classes, all six
   required by the probe's own positive control: `LIVE`, `FROZEN`, `MIXED`,
   `IMPORT_USED`, `DEAD`, `ABSENT`. Two of them (`IMPORT_USED`, `DEAD`) were
   added after this file's **own first run returned `UNUSED` for two names
   with entirely different diagnoses**, which would have had them fixed the
   same wrong way — `graphene_fet_model.t_ox` is consumed at import to derive
   `C_ox` and is correctly handled by the 09-29 audit (verified mechanically,
   §3b, rather than taken from that audit's comment), whereas `ALPHA_ABS` was
   simply dead.

### Four findings against the repository

4. **A ONE-CHARACTER TYPO IS INDISTINGUISHABLE FROM INSENSITIVITY, on the real
   instrument.** Misspelling the first 09-29 `MUST_CHANGE` patch key by one
   character gives output **bitwise identical to no mutation over 2004
   harvested points**: `setattr` creates the misspelled attribute and returns
   normally. A renamed or mistyped mutation target is a write-only knob **with
   no syntactic trace anywhere in the source**, which the 09-29 AST census
   could not have found by construction.
   **Second-order, and the sharper half:** the probe's verdict on that key
   degrades from `ABSENT` to `DEAD` once the mutation has been applied, because
   `importlib.reload` re-executes the source into the same module dict without
   clearing it. **Applying a misspelled mutation destroys the evidence that it
   was misspelled.** Arrival must be classified before the mutation, and that
   ordering is a requirement, not a style choice.

5. **`ALPHA_ABS` was a decorative name, and the check policing it was vacuous.**
   `graphene_photodetector_model.ALPHA_ABS` was read by no function **and
   nowhere in the module body**. So (i) the module docstring's claim to "reuse
   the existing 2.3% universal-absorption result" was false — the only route
   absorption could take is that name, and `EQE_BARE = 0.0015` is an
   independent literature literal. This is **09-28's "prose is a detector" for
   the third time and INVERTED**: 09-28 and 09-29 found code contradicting its
   prose; here the prose asserts a dependency the code does not have. Distinct
   from 09-29 item 9, which found this module re-typing `pi/137.036` — the
   stronger statement is that it used the value by neither route. And (ii) the
   09-29 null control on it, logged as "a real finding if it holds: the figure
   must not double-count absorption", **could not have failed**: with the name
   unread it would have passed had the figure double-counted absorption,
   counted it once, or not counted it at all. The evidence for that was already
   in yesterday's log, one item away, unconnected.
   **The fix produces a physical number that was implicit until now.**
   `internal_collection_efficiency()` reads both constants live and reports
   `EQE_BARE / pi*alpha` = **6.5430%**: roughly fifteen of every sixteen
   absorbed photons are assumed lost before collection. That is the
   quantitative content of the phrase "collection bottleneck" the existing
   comment uses without a number, and it underlies every bare-device
   responsivity in the thesis — a Chapter 6 modelling assumption, not
   bookkeeping. `EQE_BARE` deliberately stays a literature literal rather than
   becoming derived, because the physical direction runs the other way.
   **No number moved:** `(EQE_BARE/ALPHA_ABS)*ALPHA_ABS == EQE_BARE` bitwise,
   and every pre-existing `summary_numbers()` line is unchanged (the two new
   lines are appended, not inserted).

6. **Three of the 09-29 null controls are identity re-assignments.**
   `graphene_fet_model.T -> 300.0` when `T` is already 300.0;
   `n_puddle -> 5e15` when it already is; `AU_FERMI_SHIFT_SURFACE_EV -> 0.14`
   when it already is. Each passes with `setattr` replaced by `pass`, with the
   constant frozen, and with the model deleted. **The 09-28 cannot-fail class
   living inside the 09-29 null controls**, which were themselves offered as
   the defence against over-reporting.

7. **`MIXED` is the majority class, not the exception — the top item was
   mis-sized.** Yesterday treated the 37 occurrences as one fault with `T_HOP`
   singled out as "the worst". At the level of **names**, **15 of 20 were in
   `T_HOP`'s class**, including `A_CC` in the very same file, which went
   unmentioned. The oracle/measurement split also exists outside the audits:
   `D0_SEPARATION`, `H_DIFF`, `W_CROSS_CHEM`, `ELL_DEFAULT` (4 frozen readers
   against 1 live) and `N_POINTS` in three separate photodetector modules.

### The fix, measured rather than described

8. **All 14 audit-module sites converted, under a bitwise no-change
   requirement.** `graphene_band_structure_audit` (4 sites / 3 functions),
   `gds_quadrature_audit` (6 / 4), `rootfinder_audit` (3 / 1),
   `sensitivity_audit` (1 / 1). All four audits run to completion before and
   after; all four stdout transcripts (8232, 9463, 6388, 8421 bytes) compare
   **IDENTICAL**, and each still reports its own full pass count.

9. **`T_HOP` before and after, by running the pre-fix file out of git beside
   the post-fix one** — 09-29 described this hazard; it is now measured:
   ```
   before  T_HOP MIXED: oracle 2*T_HOP  5.600 -> 11.200      (moved)
                        bands()[0]     -7.225687 -> -7.225687  (FROZEN)
   after   T_HOP LIVE : oracle 2*T_HOP  5.600 -> 11.200      (moved)
                        bands()[0]     -7.225687 -> -14.451374 (moved)
   ```
   Doubling the hopping integral used to double every expected value in the
   module and leave the computed band energies bitwise unchanged.

10. **A standing guard, so the fix does not depend on the next session having
    read this entry.** The probe's §7 fails if any `*audit*.py` in the
    repository captures a module-level constant as a default again.

11. **`notes/2026-09-30-arrival-not-response-...md`**, 209 lines, with the full
    argument, the four-row proof table and the delivery-class table.

**Methodological note, continuing the series.** 09-25: a procedure asked
whether it has converged can answer yes and be 44% wrong. 09-26: an anchored
comparison is anchored in one variable. 09-27: an instrument can be
systematically smallest where the answer is worst. 09-28: prose is a detector,
and a check that cannot fail is worse than no check. 09-29: a mutation that
does not arrive is indistinguishable, in the output, from a system that does
not respond.
**09-30: A CONTROL HAS TO SIT WHERE THE FAILURE ENTERS, NOT WHERE IT SHOWS.**
09-29 identified the failure correctly and then placed its control downstream,
in the output — the one place where, by the very property it had just
established, the two candidate causes are already indistinguishable. A defence
built at the point of measurement can only compare zeros. Arrival is a property
of the *binding*, so it has to be established at the binding, before the run;
and since that check requires no model evaluation, it costs nothing to put
there. **Corollary, and the uncomfortable one for this repository: the check
most likely to be vacuous is the one whose passing you find reassuring.** A
check that passes attracts no scrutiny, and a vacuous check always passes —
which is why today's three worst findings (§4, §5, §6) are all *passing* checks
from yesterday, not failing ones. Every previous entry in this series found its
fault in something that gave a wrong answer; today's came from three things
that gave the right answer for no reason.

**Validations:** 29/29 (8 positive controls on the probe itself; 4-row null-
control theorem; 2 real-instrument demonstrations; 19-target delivery census;
1 import-consumed coverage check; 5 before/after `T_HOP` measurements; 1
standing guard), plus 4 audit transcripts compared bitwise identical and 1
bitwise round-trip check on the `EQE_BARE`/`ALPHA_ABS` decomposition. Two
failures were found by this file against itself and are recorded in §3 and §4
rather than quietly fixed.

**Not yet covered (candidates for future runs):**
- **The 10 remaining knobs in 6 MODEL modules, 7 of them MIXED** — the
  audit-module half is done, and this is what is left of the 09-29 item and
  **the new top item**. `W_CROSS_CHEM` (9 frozen readers against 2 live) and
  `ELL_DEFAULT` (4 against 1) are the worst. Unlike the audit modules these are
  not currently mutated by anything, so the hazard is prospective: the next
  mutation-based result built on any of them starts out measuring nothing
- **`require_delivery` is not actually CALLED by any audit yet** — created
  today. The probe proves the rule and the audits do not yet obey it; wiring it
  into the 09-29 provenance audit and the four fixed audits is the obvious next
  step, and it must be done without moving their transcripts
- **Are there OTHER dead names in this repository?** — created today by §5.
  `ALPHA_ABS` was found only because a mutation happened to target it. A sweep
  for module-level constants read nowhere is a ten-line AST pass and would
  catch the rest of the class, including ones no audit touches. Each dead name
  is a provenance claim with no code behind it
- **Do any OTHER module docstrings claim a dependency the code does not have?**
  — created today, and the inverted form of the 09-28 detector. §5 found one by
  accident; the general check is to compare each module's stated inputs against
  the names it actually reads
- **Whether the three identity re-assignments should be replaced by real
  mutations** — created today by §6. Replacing them changes a passing check
  into one that can fail, which is the point, but it must be done knowing which
  of the three the model is genuinely insensitive to
- **Whether Chapter 4's `n_puddle` actually reproduces Section 3.3's measured
  6.45 kOhm/sq floor** — created 09-29, untouched, a one-line calculation, and
  the only place Chapter 3's central negative result touches Chapter 4's numbers
- **The Section 3.6 Pauli edge is absent from Chapter 6's model** — created
  09-29, untouched; a THIRD leg of the Chapter 4 / Chapter 6 contact-metal
  contradiction, and the only one that is a physical mechanism rather than a
  modelling disagreement: biasing to high `E_F` for low contact resistance
  switches off the detector's own absorption
- **Where does the 6.5430% internal collection efficiency come from?** —
  created today by §5. It is now a named number with no derivation; a transit-
  time-versus-recombination-lifetime estimate would either support it or make
  `EQE_BARE` the quantity in tension with the rest of Chapter 6
- **A finite-temperature optical conductivity** — created 09-29, untouched;
  Section 3.6's edge is the `T = 0` step function
- **The remote-polar-phonon cap on SiO2 is asserted, not computed** — created
  09-29, untouched
- **Angular trigonal warping** (the `cos 3theta` structure) — created 09-28,
  untouched
- **Finite-temperature carrier density `n(E_F, T)`** — created 09-28, wanted
  twice
- **Whether other `== 0.0` exactness checks here are round-trip tautologies** —
  created 09-28, untouched, and today's `EQE_BARE`/`ALPHA_ABS` round-trip check
  is an instance of the class that must not be allowed to become one
- **Whether RESULT 2's other convictions are read from the saturated branch** —
  created 09-28, untouched
- **A probe that reports both directions of lambda** — created 09-28, untouched
- **Whether `t'` can be EXCLUDED quantitatively** as the source of Chapter 4's
  electron-hole asymmetry — created 09-28, untouched
- **Whether any OTHER near-cancellation would show the `c2`/`c4` sign flip if
  differentiated rather than evaluated** — created 09-27, untouched
- **Migrating every remaining validation to measured-value-beside-derived-bound
  form** — created 09-27, untouched
- **Whether Chapter 4 or 5 contains a RANKING that is a step artefact** —
  created 09-23/09-25, untouched
- **`n_segments = 50` at a bias with more curvature** — created 09-26, untouched
- **A second anchor for `Delta_c`, at any separation other than 3.3 A** — open
  since 2026-09-21, still the top *physics* item, now **untouched for ten
  consecutive sessions**
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural weakness
- **Whether Chapter 7's Section 7.3 differential-interconnect prediction
  holds** — created 09-27
- **Which of Chapters 4 and 5's design rules could be restated as parities or
  bounds** — created 09-23; whether the parity survives a photo-thermoelectric
  term — created 09-23
- Re-check whether other "for every ..." claims rest on small samples — open
  since 2026-09-20; whether Chapter 4's contact-resistance results should be
  re-run at the 5.4 eV crossover — open since 2026-09-18; reconciling Mueller
  *et al.*'s 0.12 eV step (arXiv:0902.1479) with the 0.25-1.07 eV offsets
  `METAL_WORK_FUNCTIONS` assumes — open since 2026-09-18; a photo-thermoelectric
  term (Kasirga review); Shimomura *et al.*'s comb-electrode design; integrating
  6.5's plasmonic near-field picture with the spatially-resolved contact-doping
  machinery; isolating the root cause of the Section 4.7 negative residual (open
  since 2026-08-31); Ti and Cr per-metal `Rc` recalibration (ResearchGate
  rate-limiting); a second independent edge-contact dataset (Lee *et al.* 2022,
  Wiley 403'd)

**Automation health:** Device reachable at the **04:34 UTC** firing, the first
of the day's three, so the redundancy was again not needed. Step 0's
already-ran check was clean: neither repo had a 2026-09-30 entry or a commit
since midnight. `git clone` of both repos completed normally; `user.name`/
`user.email` again absent in the fresh clones and set per 09-24. `scipy` again
absent and `pip install scipy` succeeded immediately (1.15.3 against numpy
2.2.6). The `GIT_ASKPASS` recipe was placed in `$HOME/.sess/` per 09-28 and had
no permission trouble; every commit was pushed as it was made, per 09-25. Live
web search was not used: the session was an audit of material already in the
repository and adds no external citation. **The 09-29 authoring rule was
followed and it paid off twice:** every patch script parsed the modified source
with `ast.parse` BEFORE writing it to disk, including the 14-site late-binding
transformation, which was applied in two passes (signature edits, re-parse,
then guard insertion against recomputed line numbers) precisely because the
insertion invalidates the line numbers the first pass used.

**Commits this run:** 4 (the arrival probe with its output; the `ALPHA_ABS`
fix; the 14-site late-binding fix with the standing guard and the before/after
`T_HOP` measurement; the study note). This AUTOMATION_LOG.md entry makes 5.

## 2026-10-01 — The 09-30 fix's own before/after check could only pass in the run that wrote it, and a cited table column had been dead for ten sessions

**Status:** Automated session. **The 04:30 UTC firing did not reach the
machine; this run is the 06:59 one.** Step 0's already-ran check was clean:
neither repository had a 2026-10-01 entry or a commit since midnight. Live web
search **not used** — the session is an audit of material already in the
repository plus one unread column of an already-cited table, and it adds no
external citation. `scipy` was present without installing (1.17.1 against
numpy 2.4.4), for the first time in these sessions. See **Automation health**
for the one real operational problem, which was the device's network.

**The 09-30 top item closes, with two corrections to how it stated itself.**
That item was "the 10 remaining knobs in 6 MODEL modules, 7 of them MIXED".
It is **seven** modules — `graphene_differential_crossover_model`'s `U_SMALL`
was in the census listing and not in the module count — and **23 sites**, not
10. Ten was the number of *names*: `W_CROSS_CHEM` alone is captured at 8 sites
in one module and `ELL_DEFAULT` at 4 in another, which is exactly why 09-30
called those two the worst, since the 9-frozen-against-2-live and
4-against-1 ratios it quoted were **site** counts. The fix is per site.

### 1. All 23 sites converted, under a bitwise no-change requirement

| module | transcript |
|---|---|
| `graphene_contact_doping_nonlinear_model.py` | 7376 bytes IDENTICAL |
| `graphene_crossover_sensitivity_model.py` | 5909 bytes IDENTICAL |
| `graphene_differential_crossover_model.py` | 7153 bytes IDENTICAL |
| `graphene_per_metal_crossover_model.py` | 3434 bytes IDENTICAL |
| `graphene_photodetector_collection_model.py` | 1206 bytes IDENTICAL |
| `graphene_photodetector_signed_carrier_model.py` | 3728 bytes IDENTICAL |
| `graphene_photodetector_two_contact_model.py` | 2728 bytes IDENTICAL |

The probe's repository-wide census (Section 5) now reports `totals: {}`, down
from 7 MIXED and 3 FROZEN. **The hazard was prospective, not live** — nothing
currently mutates these names, so no committed number was wrong. What they
would have done is start the *next* mutation-based result built on any of them
at measuring nothing, which is why the item outranked the physics. Section 7b
extends 09-30's audit-module guard over the model modules so they cannot come
back. Authored to the 09-29 rule as extended on 09-30: `ast.parse` **and**
`compile` before writing, in two passes, with guard insertion against
**recomputed** line numbers.

### 2. THE 09-30 FIX'S OWN BEFORE/AFTER CHECK HAS BEEN FAILING SINCE IT WAS COMMITTED

Re-running the probe after the model-module fix reported **two failures** in
Section 6 — the before/after measurement of `T_HOP`'s manufactured
disagreement, which is the centrepiece of 09-30's §9. It read its `before`
source with `git show HEAD~2:graphene_band_structure_audit.py`.

**`HEAD~2` names a POSITION in history, not a state, and positions move.** On
09-30 it resolved to the pre-fix file because the fix was still in the working
tree and two commits had landed since the 09-29 log. The moment the fix was
itself committed, `HEAD~2` came to name the **post-fix** file — so `before`
and `after` became the same source, the oracle/measurement split collapsed to
two equal numbers, and both checks have failed ever since.

**Verified rather than argued:** a pristine clone checked out at `c99c48d`,
the 09-30 state exactly as committed, runs the probe to **27 passed, 2
failed**. **The 09-30 entry's "Validations: 29/29" is true of a working tree
that was never committed and false of every clone of this repository since.**

Fixed by pinning the reference to **content**: blob
`21a7e2ad96dd07f79d56807ca1f42eaccb8ecc29`, the pre-fix file as it stood at
`8af4537~1`. A blob hash is content-addressed, so it cannot come to mean a
different file. **Section 8** is the standing guard — no source file here may
reference `HEAD`, `HEAD~n`, `HEAD^` or `ORIG_HEAD` — because fixing one site
does not stop the next being written the same way.

Plus a **positive control on the reference**, 09-30's rule one level out: the
fetched blob must contain the capture the fix removed and must differ from the
file on disk. **Its first form fired, on me** — it looked for `'T_HOP=T_HOP'`
when the signature is `def bands(k, t=T_HOP)`, the parameter not being named
after the constant. Recorded rather than quietly re-tuned, because a control
that has never fired is a control with no evidence that it can.

**The test that matters is whether the fix survives its own commit.** A
pristine clone of the new HEAD runs the probe to **32 passed, 0 failed** —
the test the 09-30 version failed.

### 3. `graphene_dead_name_sweep.py`: a citation is a claim, and nothing was checking it

09-30 asked whether there are **other** dead names, because `ALPHA_ABS` was
found only because a mutation happened to point at it. 113 module-level
constants classified: **112 LIVE, 1 BODY_ONLY, 2 CROSS_MODULE, 1 EXEMPT, 0
DEAD** after four dispositions, 6/6 checks.

A dead constant is worse than an unused variable because each carries a
**citation**, so writing the name claims a result depends on the number — and
a false claim in the direction that flatters the work, since a reader counts
the citation as a dependency honoured.

- **`PASSI_RC_DIRAC_OHM_UM`** — six measured numbers whose own comment says
  the table is "on-state (V_BG = −40 V) **and** Dirac-point contact
  resistance". Only the on-state column was ever read. This is §4 below.
- **`D_CHEM_PHYS` = 0.9 eV** is a **second definition** of `DC_ANCHOR` =
  0.9 eV — one literature value, two copies, one unread. The physics is
  unaffected (the term is applied through `w_cross_for_metal` → `delta_c` →
  `DC_ANCHOR`), which is why it sat unread for ten sessions and why a revised
  extraction would move one copy silently. The modules cannot import each
  other, so the agreement check lives in the sweep.
- **`DELIVERED`** in the probe, defined for symmetry and read by nothing,
  given a reader: the census's conclusion stated positively.
- **`FANG_ANTENNA_TEST_WAVELENGTH_NM`** is **deliberately** unread. So the
  exemption is a source-level `# dead-name-exempt: <reason>` marker with the
  reason **required** — a reasonless marker reports as DEAD. It lives in the
  source and not the checker, because a checker that reports a documented
  decision as a fault is one that gets ignored, and that failure is silent.

**Two findings against the sweep itself, both from running it.** (i) **Its
cross-module check was a text search and its own docstring exercised the
loophole:** `re.search` cannot tell a *use* from a *mention*, so the paragraph
explaining the Fang exemption made that constant look alive — and this
module's own `CLASS_LIVE` came back CROSS_MODULE purely because the probe
defines a constant with the same spelling. **It was DEAD in here and the
loophole hid it**: a dead-name detector whose own dead name is concealed by
its own mechanism, one run from being committed. Now an AST pass over `Name`
loads, `Attribute` names and `ImportFrom` aliases. (ii) **Its exemption scope
was three lines** and the only real exemption here needs a four-line reason,
so it reported a declared exemption as DEAD.

Section 3 of the sweep is 09-28's "prose is a detector" **inverted**, as the
09-30 list asked: constants a module's own docstring or comments NAME while no
function reads them. It reports 0, and is deliberately a report rather than an
assertion — prose may mention a constant it does not depend on; what it may
not do is cite one as an input.

### 4. The physics: Chapter 4's ~11x was compared against the wrong column

**Rc(Dirac)/Rc(on-state)** is the *gate-tunable multiple* of the contact
resistance — near 1 means the interface sets it and the carrier density
beneath does not; large means graphene's own density of states in and near the
contact dominates, which the gate controls.

| hole D (nm) | on-state | Dirac | ratio |
|---|---|---|---|
| 0 | 519 | 1372 | 2.64 |
| 50 | 212 | 620 | 2.92 |
| 100 | 352 | 732 | 2.08 |
| **200** | **45** | **456** | **10.13** |
| 500 | 410 | 1354 | 3.30 |
| 1000 | 560 | 1590 | 2.84 |

2.08–3.30 at five of six diameters and **10.13** at D = 200 nm — the one
diameter whose 11x Chapter 4 quotes, and **3.7x the mean of the other five**.

**So the headline reduction is a gate-bias statement, not a geometry
statement.** D = 0 → D = 200 nm is **11.53x** in the on state and **3.01x** at
the Dirac point, for the same two devices and the same etched geometry.

**And 3.01x is the comparison the model's own physics selects.** `R_extra` is
computed from the *contact-induced* carrier density alone and the gate does not
enter it; at the Dirac point the gate contributes no channel carriers, so the
measured access resistance is dominated by exactly the term the model
isolates. **Against that column the model's 2.59x (982.1 → 379.0 Ω·µm) is
within 14% of measured, where against the on-state column it was short by a
factor of 4.5.** Chapter 4 §4.8 held a gate-independent model against an
on-state measurement, which is a category error; **§4.8.2 annotates it in
place and withdraws no number**, and the §4.9 status row is annotated to
match. Two caveats kept rather than absorbed: the 2.59x is a *pure-mode* ratio
against a device that is a mixture at a fill fraction Passi et al. do not
report, and §4.7's absolute-magnitude discrepancy is untouched and still open.

**And §4.8.1's open item is sharpened rather than closed.** Whatever produces
the D = 200 nm optimum is **3.5x more of an effect in the on state than at the
Dirac point** (4.71x against 1.36x, taking D = 50 nm as reference), so it
scales with **carrier density**, not interface area. A model built entirely
from interface geometry — which §4.8.1's is, edge length per unit area at
fixed fill fraction — cannot produce a minimum whose depth depends on the
gate. **Its failure was the wrong class of model, not a missing refinement,
and the evidence was in the column the chapter cited and did not read.**

**Methodological note, continuing the series.** 09-25: a procedure asked
whether it has converged can answer yes and be 44% wrong. 09-26: an anchored
comparison is anchored in one variable. 09-27: an instrument can be
systematically smallest where the answer is worst. 09-28: prose is a detector,
and a check that cannot fail is worse than no check. 09-29: a mutation that
does not arrive is indistinguishable from a system that does not respond.
09-30: a control has to sit where the failure enters, not where it shows.
**10-01: AND IT HAS TO NAME WHAT IT COMPARES AGAINST IN A WAY THAT CANNOT COME
TO MEAN SOMETHING ELSE.** Three instruments failed today and **not one had a
wrong comparison**. The before/after measurement compared correctly, against a
revision expression that came to denote a different file. The cross-module
check compared correctly, against a name resolved by text so a mention counted
as a use. The exemption check compared correctly, against a marker whose scope
was three lines and whose reason was four. **Every one is a naming failure,
and a naming failure is invisible at the comparison, because the comparison is
doing exactly what it says.** 09-30 put the control where the failure enters;
*where* turns out to be half of it, and the other half is that the thing being
pointed at has to stay the thing being pointed at.
**The uncomfortable corollary:** the one that mattered most was in the
instrument 09-30 built *to enforce* this series' own rules, and it broke **by
being committed** — working in the run that wrote it and in no run afterwards,
which is the worst possible failure schedule, since it passes under scrutiny
and fails only once nobody is looking. **A check validated only in the session
that wrote it has not been validated**, and the cheapest test for it — run the
suite from a pristine clone of the commit — is one this repository had never
run. It ran twice today and found the fault immediately both times.

**Validations:** probe 33/33 (was 29 checks reporting 27/2 as committed),
including the new positive control on the reference, Section 7b over the model
modules and Section 8's relative-revision guard; sweep 6/6, including a
five-class positive control on a synthetic module, a check that a reasonless
exemption reports as DEAD, and the two-copy agreement check on the chemical
term; seven module transcripts compared bitwise identical; `graphene_edge_
contact_model.py`'s existing output bitwise identical with the new block
appended; and two pristine-clone runs, one at `c99c48d` to establish the
regression and one at the new HEAD to establish that the fix survives commit.
**Three faults of my own were found and are recorded rather than quietly
corrected:** the reference control's first form looked for the wrong text, the
sweep's cross-module check hid its own dead constant, and the sweep's
exemption scope was shorter than its only real exemption's reason.

**Not yet covered (candidates for future runs):**
- **Run every suite in this repository from a pristine clone of HEAD, not from
  the working tree** — created today and **the new top item**, because today's
  worst finding is that the difference between those two is where a check goes
  to die, and nothing here has ever checked it. The probe was validated this
  way twice today; the other audits never have been
- **`require_delivery` is still not CALLED by any audit** — created 09-30,
  untouched. The probe proves the rule and the audits do not obey it; wiring it
  into the 09-29 provenance audit and the four fixed audits is the obvious next
  step and must not move their transcripts
- **Do the four existing audits have references that are not content-pinned?** —
  created today by §2. Section 8 guards against relative *git* revisions; a
  reference to "the current value of X" or to a file by path is the same class
  of hazard and is not guarded
- **Is `D_CHEM_PHYS` the only duplicated literature value?** — created today by
  §3. One sweep found one duplicate pair by accident, which is exactly the
  position `ALPHA_ABS` put us in a session ago. A census of numerically equal
  module-level constants across modules is the general check
- **Whether the three identity re-assignments should be replaced by real
  mutations** — created 09-30, untouched
- **Whether Chapter 4's `n_puddle` actually reproduces Section 3.3's measured
  6.45 kΩ/sq floor** — created 09-29, untouched, a one-line calculation, and
  the only place Chapter 3's central negative result touches Chapter 4's numbers
- **The Section 3.6 Pauli edge is absent from Chapter 6's model** — created
  09-29, untouched; the third leg of the Chapter 4 / Chapter 6 contact-metal
  contradiction and the only one that is a physical mechanism
- **Where does the 6.5430% internal collection efficiency come from?** —
  created 09-30, untouched; a transit-time-versus-recombination-lifetime
  estimate would either support it or put `EQE_BARE` in tension with Chapter 6
- **Chapters 2 and 3 are drafted but §4.8.2's form of error is not audited for
  elsewhere** — created today: the fault was comparing a gate-independent model
  against a gate-biased measurement, and Chapters 5 and 6 contain several
  model-against-literature comparisons whose bias conditions are not stated
- **A finite-temperature optical conductivity** (09-29); **the remote-polar-
  phonon cap on SiO2 is asserted, not computed** (09-29); **angular trigonal
  warping** (09-28); **finite-temperature carrier density `n(E_F, T)`**
  (09-28, wanted twice); whether other `== 0.0` exactness checks are
  round-trip tautologies (09-28); whether RESULT 2's other convictions are
  read from the saturated branch (09-28); a probe reporting both directions of
  lambda (09-28); whether `t'` can be EXCLUDED quantitatively as the source of
  Chapter 4's electron-hole asymmetry (09-28); whether any other
  near-cancellation shows the `c2`/`c4` sign flip if differentiated (09-27);
  migrating every remaining validation to measured-value-beside-derived-bound
  form (09-27); whether Chapter 4 or 5 contains a RANKING that is a step
  artefact (09-23/09-25); `n_segments = 50` at a bias with more curvature
  (09-26)
- **A second anchor for `Delta_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, now **untouched for eleven
  consecutive sessions**. Today's §3 touched `DC_ANCHOR` without adding an
  anchor, which is worth saying plainly
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural weakness
- **Whether Chapter 7's Section 7.3 differential-interconnect prediction
  holds** (09-27); which of Chapters 4 and 5's design rules could be restated
  as parities or bounds (09-23); whether the parity survives a
  photo-thermoelectric term (09-23); re-check whether other "for every ..."
  claims rest on small samples (09-20); whether Chapter 4's contact-resistance
  results should be re-run at the 5.4 eV crossover (09-18); reconciling Mueller
  *et al.*'s 0.12 eV step (arXiv:0902.1479) with the 0.25–1.07 eV offsets
  `METAL_WORK_FUNCTIONS` assumes (09-18); a photo-thermoelectric term (Kasirga
  review); Shimomura *et al.*'s comb-electrode design; integrating 6.5's
  plasmonic near-field picture with the spatially-resolved contact-doping
  machinery; isolating the root cause of the Section 4.7 negative residual
  (08-31); Ti and Cr per-metal `Rc` recalibration (ResearchGate rate-limiting);
  a second independent edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

**Automation health.** Device reachable and folder connected, but **the real
problem today was the DEVICE'S NETWORK, not the toolchain.** `git clone` of
either repository from inside `device_bash` **could not complete**: raw
throughput to GitHub measured **~13 KB/s** (1.3 MB of tarball in 97 s), so a
full clone exceeded the 180 s per-call shell limit repeatedly, and `nohup`'d
background clones **do not survive between `device_bash` calls** — each call
is a fresh shell and the children were reaped. A blobless `--no-checkout`
clone finished in **3 seconds**, which locates the fault precisely: git
protocol negotiation is fine, bulk packfile transfer is throttled. The session
was restructured rather than abandoned, and the recipe is recorded because it
will be needed again:
1. Clone and do all work in the **cloud sandbox** (2 s for both repositories).
2. `git bundle create <f> <old_origin_main>..main`.
3. Ship the bundle to the device with `device_commit_files` into
   `C:\scheduled harsh\_transfer\`.
4. On the device, `git clone --depth 1 --filter=blob:none --no-checkout` (3 s),
   `git fetch <bundle> refs/heads/main:refs/remotes/incoming/main`, then
   `git push origin refs/remotes/incoming/main:refs/heads/main`.
**Pushing from a shallow, blobless, no-checkout clone WORKS**, and this was
tested on the other repository's first commit **before** the rest of the
session's work was done, precisely so a broken push path would be found early
rather than at the end. A push is a few KB, so 13 KB/s is no obstacle to it;
only the clone was. The `GIT_ASKPASS` recipe ran from `$HOME/.sess/` per 09-28
and the 09-29 correction (**`/tmp` may carry another session's files**), and
the temp token copy was shredded. **Deviation from the 09-25 push-as-you-go
rule, stated deliberately:** each push now costs a bundle, a file transfer and
a device round-trip, so commits were made locally and pushed in one batch per
repository, then verified against the GitHub API rather than against git's own
output. `user.name`/`user.email` were again absent in the fresh clone and set
per 09-24. `scipy` was already present (1.17.1), so no install was needed.

**Commits this run:** 5 (the 23-site late-binding fix; the pinned reference
with its positive control and the relative-revision guard; the dead-name sweep
with its four dispositions; Chapter 4 §4.8.2; the study note). This
AUTOMATION_LOG.md entry makes 6. The Design-Verification-Roadmap repository
took 5.

## 2026-10-02 — The 10-01 top item closes, and the day's worst finding was inside the one suite that agreed with its transcript perfectly

**Status:** Automated session, the 04:30 UTC firing. Step 0's already-ran check
was clean: neither repository had a 2026-10-02 entry or a commit since
midnight. Live web search **not used** — the session is an audit of material
already in the repository and adds no external citation. `scipy` present
without installing (1.18.1 against numpy 2.5.3). See **Automation health**: the
device's network was again unusable for cloning, as on 10-01, and the recorded
workaround was used unchanged.

**The 10-01 top item closes.** That item was "run every suite in this
repository from a pristine clone of HEAD, not from the working tree", created
because 10-01's worst finding was that the difference between those two is
where a check goes to die. `graphene_pristine_transcript_audit.py` does it for
every suite and asks one step more than 10-01 asked — not "does it still pass"
but **does the committed code still produce the committed transcript.** The
transcripts here are captured stdout written by a shell redirect, so a session
that changes a module is free to forget it, and **`git status` cannot see the
result** because the transcript is unmodified and it is the code that moved.

Opening state at `2d51dd8`: **4 STALE, 1 NUMERIC_DRIFT**, and one of the new
instrument's own controls failing on a premise of its own. Closing: **0 STALE**.

### 1. THE STRONGEST VERDICT WAS WRONG, AND THAT IS THE DAY'S FINDING

`graphene_band_structure_audit.py` classified **IDENTICAL** — code and
transcript agreeing byte for byte — and the transcript said, under a heading
reading `RESULT 1 -- the shipped band structure has NO DIRAC POINT`:

```
  minimum gap over the whole shipped path = 0.0000 eV
  gap at the index labelled K             = 0.0000 eV
  ... The shipped figure shows a 0.0 eV gap at its own K label.
```

**A 0.0 eV gap at K *is* the Dirac point.** The sentence offers the absence of
the bug as the evidence for the bug. `RESULT 1` reports the *superseded*
spectrum and was reading the *corrected* path, so **this repository's own bug
report has been quoting the fixed value as the fault for five weeks.** Now
prints superseded beside corrected: **14.4019 eV** and **2.487e-15 eV** at K,
both on the record, neither deleted.

The instrument built today could not have found this. It asks whether two
artefacts agree; these agreed perfectly, on something false.

### 2. V5 could only FAIL, and the oracle it needed was built in the same commit as the fix

| commit | date | V5 | tally |
|---|---|---|---|
| `64c6347` | 09-28 | **PASS** | 5/5 |
| `fef532d9` | 09-28 | **FAIL** | 4/5 — and in every clone since |

Same mechanism as §1. `legacy=True`, which preserves the superseded path
exactly, **was added in `64c6347`, the fix's own commit**, explicitly under the
09-24 rule that a deleted implementation destroys the only oracle for judging
its replacement — and the one validation whose purpose is to exercise that
oracle was never pointed at it. **The oracle was preserved and then not used.**

It survived because the suite prints FAIL on a summary line and exits 0. **A
standing FAIL in a tally is worse than a missing check**: it teaches a reader
that 4/5 is this suite's normal state, which is the exact condition under which
a real regression here would be invisible. 09-28's "a check that cannot fail is
worse than no check" has a mirror, and the mirror is louder and therefore
easier to stop hearing. V5 now reads the legacy path; new **V5b** asserts the
default path *has* a Dirac point. Before, V5 could only fail; now V5 fails if
the oracle is lost and V5b fails if the fix is reverted. **4/5 → 8/8.**

### 3. The third site: the 09-28 missing π was still live in the repository

`plot_fermi_surface()` carried its own inline copy of the wrong zone corner,
`k_K = 4.0/(3*np.sqrt(3))` = 0.76980, π times too small — so *"Fermi Surface of
Single-Layer Graphene (Around K Point)"* was centred where the gap is **14.40
eV**: a figure of the Dirac cone with no Dirac cone in it.

| centre | value | gap there |
|---|---|---|
| superseded (inline literal) | 0.7698003589195 | **14.401937 eV** |
| corrected | 2.4183991523123 | 2.487e-15 eV |

**Nothing could have caught it:** the function is called by no suite and
`fermi_surface.png` is not committed, so the only artefact that could have
disagreed with the code does not exist. **A literal is not fixed by fixing the
constant it duplicates** — the dead-name sweep guards unread *names*, and this
is the mirror, a read *value* that is not a name and so is invisible to every
name-based instrument here. Now reads `fermi_surface_grid_centre()`; **V6**
requires the gap at that centre to sit on the floor, **V6b** forbids an inline
copy outside `HIGH_SYMMETRY_LEGACY`. V6b mutation-tested: reinjecting the
literal takes the suite 8/8 → 7/8 and names the line.

### 4. A ranked table of 21 equal keys

The audit's first run reported `graphene_differential_crossover_model.py` STALE
on a *structural* difference: the committed transcript named `Ti/Pd` in a table
and a pristine run of the **same commit** named `Cr/Cu`, with every printed
number identical. The cause is printed four lines below the table by the
function's own prose — `asym` is zero for all 21 pairs to machine precision
(Validation 6: N is exactly odd in s) — so `sort(key=-asym)` sorted 21 equal
keys and `rows[:6]` took an arbitrary six, chosen by the last bits of a
quantity already proved zero. Measured spread **5.574e-15**. No number and no
verdict moves; P2 and P4 are still falsified by the identity. What was wrong is
that **a reader counts a ranked table as a claim that the top row differs from
the bottom one.** Tie-broken lexicographically and the degeneracy printed. This
answers the 09-23 "is there a RANKING that is a step artefact" item in its
**noise** form; the step form stays open.

### 5. CROSS_MODULE resolved on SPELLING, and this session's own module walked into it

The new module defined `SELF` and never read it, and the sweep did **not**
report it DEAD — it reported CROSS_MODULE, resolved against a **function-local
variable of the same spelling** in `graphene_log_sensitivity_step_audit.py`, a
module that does not import the new one. 10-01 recorded this exact fault in the
sweep and **fixed half of it**: the mention-vs-use half, `re.search` → AST. The
half left standing is that an AST reference in another module is a liveness
claim *only if that module can reach this one*. New **Section 2b** checks every
resolution against a real import, with a positive control: 3 rows, **1
unsupported**. The two supported are real and untouched (`per_metal` imports
`CHEMISORBED`; the sweep imports `nl.D_CHEM_PHYS`). `SELF` was given a *reader*
rather than a deletion — a runtime assertion that the module is not in `SUITES`,
since the exclusion had been carried by a comment and a comment is not
checkable. Verified by adding it to `SUITES` in a copy: the assertion fires.

### 6. Two transcripts that were already stale, neither visible to `git status`

- **`figure_provenance_audit_output.txt`** — stale since 10-01. Still reported
  **37 frozen default-capture sites across 11 files** (the sites 10-01
  converted; the census now reads 0 of 36) and still showed the `T_HOP`
  mutation **frozen** at `8.379013 -> 8.379013` where it now **arrives** at
  `-> 16.758026`. The transcript was a session out of date on exactly the fault
  10-01 fixed, which is the sharpest possible demonstration that regenerating
  transcripts by hand is partial by nature.
- **`log_sensitivity_step_audit_output.txt`** — stale since 09-28, when the
  covariance probe audit landed: **11 call sites → 18**, and
  discarded-convergence-flag comparisons **1 → 3**. A reader of the committed
  transcript would have concluded the repository was cleaner than it is.

### 7. One fault of my own, recorded rather than quietly re-tuned

Control **C2** bumped the last printed digit of a `%.6f` number and asserted
`NUMERIC_DRIFT`. For `7.500000 → 7.500001` that is a relative change of
**1.3e-7**, a hundred times *outside* `NUMERIC_RTOL = 1e-9`, so `STALE` was the
correct verdict and the control's premise was false. **"The last digit" is a
fact about a format string; the tolerance is a fact about the number.** The
classifier was right about the real data — the drift it exists to catch is
2.6e-16 — and only the synthetic control was wrong. C2 now perturbs by
`NUMERIC_RTOL/100`, and **C2b** asserts the perturbation is genuinely inside
the tolerance so the fault cannot return.

### 8. The FAIL that is kept deliberately

`graphene_rootfinder_audit.py` classifies NUMERIC_DRIFT at **2.586e-16**: the
root agrees to a few ULP but the last three printed digits differ by
environment (`4.99733219460098000e-11` committed against `...97870e-11` here;
residuals identical). So *"every transcript is exactly reproducible"* is **false
for this repository and should keep saying so.** These are reproducible
**results**, not reproducible **artefacts**, and byte comparison cannot certify
them across library versions. Lowering the bar to make the suite green would
delete the only statement of that limitation.

**Methodological note, continuing the series.** 09-25: a procedure asked
whether it has converged can answer yes and be 44% wrong. 09-26: an anchored
comparison is anchored in one variable. 09-27: an instrument can be
systematically smallest where the answer is worst. 09-28: prose is a detector,
and a check that cannot fail is worse than no check. 09-29: a mutation that
does not arrive is indistinguishable from a system that does not respond.
09-30: a control has to sit where the failure enters, not where it shows.
10-01: and it has to name what it compares against in a way that cannot come to
mean something else. **10-02: AND WHEN IT AGREES, THAT IS A FACT ABOUT TWO
ARTEFACTS, NOT ABOUT THE WORLD.** Today's instrument asks whether code and
transcript agree, and the day's worst finding sat inside the one suite that
agreed **perfectly** — `RESULT 1` was reproducible, pinned, byte-identical and
false. 10-01 established that a reference must name its target in a way that
cannot drift; 10-02 adds that **a reference that cannot drift can still point
at the wrong thing from the day it was written**, and agreement between two
artefacts is silent about which.
**The second thread, and the more useful one:** §2, §3 and §5 are all the same
shape — **a repair that was built and then not connected.** 09-28 preserved the
legacy path and left the validation pointed elsewhere. 09-28 fixed the
zone-corner constant and left a duplicate of its *value* three functions away.
10-01 diagnosed the spelling loophole and fixed one of its two halves. In every
case the session understood the fault correctly, wrote the right mechanism, and
**stopped one wiring step short** — and in every case what remained was
invisible *because the fix's own write-up read as complete*. The remedy is not
more care. It is that **the last step of a repair is a check that the repair is
reachable from where the fault was**, and none of these three had one.

**Validations:** pristine transcript audit 9/10 as committed (9 IDENTICAL, 0
STALE, 1 NUMERIC_DRIFT deliberately left failing per §8), with all six
controls firing correctly including C6's end-to-end injected defect; band
structure audit **8/8**, up from 5 checks reporting 4/5, with V5 repointed,
V5b, V6 and V6b added and V6b mutation-tested 8/8 → 7/8; dead-name sweep
**8/8**, up from 6/6, with Section 2b and its positive control, and
unsupported CROSS_MODULE resolutions 1 → 0; differential crossover reproducible
across two consecutive runs with every number and both verdicts unchanged;
mutation-arrival probe 33/33 unchanged with the new file present, so Section 8's
relative-revision guard stays green — the new module reads `.git/HEAD` as a file
and speaks only in 40-hex ids rather than naming the revision. **Two faults of
my own are recorded rather than quietly corrected:** control C2's first form was
mis-specified (§7), and the new module's own `SELF` was a dead name that the
sweep mis-classified, which is what exposed §5.

**Not yet covered (candidates for future runs):**
- **`require_delivery` is still not CALLED by any audit** — created 09-30,
  untouched for three sessions and **now the top item**. The probe proves the
  rule and the audits do not obey it. Today's second thread makes it sharper
  than it looked: this is itself a repair built and not connected, which is
  the exact failure shape §2, §3 and §5 all share
- **Audit every REMAINING check for the §1 fault: does it read the artefact its
  own prose names?** — created today, and the direct successor to the top item
  that just closed. V5 and RESULT 1 both read the corrected path while
  describing the superseded one, and both were found by eye, not by instrument.
  The general form is mechanisable: a function whose docstring says "shipped",
  "superseded", "legacy" or "pre-<date>" must not call the default path
- **Do the four existing audits have references that are not content-pinned?** —
  created 10-01, untouched. Section 8 guards relative *git* revisions; a
  reference to "the current value of X" or to a file by path is the same class
- **Is `D_CHEM_PHYS` the only duplicated literature value?** — created 10-01,
  untouched, and today's §3 raises it from a name question to a **value**
  question: a census of numerically equal module-level constants would not have
  found the Fermi-surface literal, because that one is not a constant at all.
  The general check is duplicated literal VALUES, not duplicated names
- **Is `plot_fermi_surface` the only uncommitted figure?** — created today. The
  reason §3 hid for five weeks is that no artefact existed to disagree with the
  code. Every plotting function whose output is not committed is in that
  position, and the count is not known
- **Whether Chapter 4's `n_puddle` actually reproduces Section 3.3's measured
  6.45 kΩ/sq floor** — created 09-29, untouched for four sessions, a one-line
  calculation, and the only place Chapter 3's central negative result touches
  Chapter 4's numbers
- **Whether the three identity re-assignments should be replaced by real
  mutations** — created 09-30, untouched
- **The Section 3.6 Pauli edge is absent from Chapter 6's model** — created
  09-29; the third leg of the Chapter 4 / Chapter 6 contact-metal contradiction
  and the only one that is a physical mechanism
- **Where does the 6.5430% internal collection efficiency come from?** —
  created 09-30, untouched
- **Chapters 2 and 3 are drafted but §4.8.2's form of error is not audited for
  elsewhere** — created 10-01; Chapters 5 and 6 contain several
  model-against-literature comparisons whose bias conditions are not stated
- **Whether Chapter 4 or 5 contains a RANKING that is a step artefact** —
  created 09-23/09-25, and today's §4 closes only its *noise* form. The step
  form is untouched, and §4 is evidence the class is real here
- **A finite-temperature optical conductivity** (09-29); **the remote-polar-
  phonon cap on SiO2 is asserted, not computed** (09-29); **angular trigonal
  warping** (09-28); **finite-temperature carrier density `n(E_F, T)`** (09-28,
  wanted twice); whether other `== 0.0` exactness checks are round-trip
  tautologies (09-28); whether RESULT 2's other convictions are read from the
  saturated branch (09-28); a probe reporting both directions of lambda
  (09-28); whether `t'` can be EXCLUDED quantitatively as the source of
  Chapter 4's electron-hole asymmetry (09-28); whether any other
  near-cancellation shows the `c2`/`c4` sign flip if differentiated (09-27);
  migrating every remaining validation to measured-value-beside-derived-bound
  form (09-27); `n_segments = 50` at a bias with more curvature (09-26)
- **A second anchor for `Delta_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, now **untouched for twelve
  consecutive sessions**. Today was again an audit session and did not touch
  it, which is worth saying plainly rather than letting the streak go unnamed
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural weakness
- **Whether Chapter 7's Section 7.3 differential-interconnect prediction
  holds** (09-27); which of Chapters 4 and 5's design rules could be restated
  as parities or bounds (09-23); whether the parity survives a
  photo-thermoelectric term (09-23); re-check whether other "for every ..."
  claims rest on small samples (09-20); whether Chapter 4's contact-resistance
  results should be re-run at the 5.4 eV crossover (09-18); reconciling Mueller
  *et al.*'s 0.12 eV step (arXiv:0902.1479) with the 0.25–1.07 eV offsets
  `METAL_WORK_FUNCTIONS` assumes (09-18); a photo-thermoelectric term (Kasirga
  review); Shimomura *et al.*'s comb-electrode design; integrating 6.5's
  plasmonic near-field picture with the spatially-resolved contact-doping
  machinery; isolating the root cause of the Section 4.7 negative residual
  (08-31); Ti and Cr per-metal `Rc` recalibration (ResearchGate rate-limiting);
  a second independent edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

**Automation health.** Device reachable and folder connected. **The device's
network was again the only real operational problem, and worse than on 10-01.**
Measured throughput to GitHub from inside `device_bash`: **9.3 KB/s** (1.34 MB
of tarball in 143 s), against 13 KB/s on 10-01, and a bare `api.github.com`
request took **11.5 s**. `git clone` could not complete inside the 180 s shell
limit at either full or `--depth 20`, and `nohup`'d background clones were
confirmed again **not to survive between `device_bash` calls** — each call is a
fresh shell and the children are reaped, which was re-tested rather than
assumed. The 10-01 recipe was used unchanged: clone and work in the cloud
sandbox, ship an incremental bundle to the device with `device_commit_files`,
then clone `--depth 1 --filter=blob:none --no-checkout` on the device (3 s) and
push from there. Push payload is a few hundred KB, so 9.3 KB/s is survivable
for the push even though it is not for the clone. The `GIT_ASKPASS` recipe ran
from the session VM's own temp space and the token copy was shredded; the token
was never written into `.git/config`, a remote URL, any repository file, or the
connected folder. **One deviation worth stating:** `pushed_at` from the GitHub
API was used as a cheap cross-check on Step 0's "has today's work already run"
question before the clones finished, which is a second, independent source for
that decision rather than a replacement for the AUTOMATION_LOG check.

---

## 2026-10-03

**Status:** Audit + correction session. The top open item of 09-30, 10-01 and
10-02 closes, and the day's real finding is a **device-physics number that was
wrong by 44 %**, not an audit-hygiene point. Live web search was not used: this
was an internal-consistency session on code and numbers already in the
repository, and nothing below rests on a literature fetch.

**Top open item closed:** *"`require_delivery` is still not CALLED by any audit
— created 09-30, untouched for three sessions. The probe proves the rule and
the audits do not obey it."*

### 1. The rule, obeyed. It passes, and that is the honest headline of Section 1

`graphene_cross_module_delivery_audit.py` (new, 25/25) runs `require_delivery`
over every mutation target `graphene_figure_provenance_audit.py` actually
applies — the only audit here that rebinds attributes and then interprets
numbers. **19 mutation applications, 15 distinct `(module, name)` targets: 18
`LIVE`, 1 `IMPORT_USED`** (`graphene_fet_model.t_ox`, whose derived child
`C_ox` is co-patched in the same patch dict, and the proviso is now checked
mechanically rather than trusted). **None of the figure-provenance verdicts
were manufactured by a mutation that failed to arrive.** Three sessions of the
rule going uncalled hid nothing in the audit it was written for.

### 2. And it could not have. A single-module instrument, a cross-module question

`classify(M, N)` asks whether the readers of `N` **inside `M`** read it at call
time. Three of the nineteen applications patch `graphene_fet_model.L`, `.mu`
and `.n_puddle` while harvesting a figure drawn by `rf_small_signal_model`. For
those, a `LIVE` verdict in the **defining** module is silent about a frozen
capture in the **consuming** one. One existed:

```python
def gate_resistance(L=gfet.L, W=gfet.W, N_fingers=1):   # until today
```

`classify(graphene_fet_model, 'W')` returns `LIVE` — **correctly**; every
reader of `W` inside that module is live. The frozen reader is a reader of a
*value copied out of the module before the audit existed*, not a reader of the
name. **No amount of care applying the 09-30 rule as written would have found
it.**

### 3. It was not latent. It defeated the model's own override, and Chapter 4's numbers

`compute_fT_fmax(..., W=W_RF)` exists to rescale the device from this repo's
1 µm per-width convention to a 40 µm / 8-finger geometry, and does it by
**rebinding `gfet.W`**. `R_g` kept the 1 µm it had captured, so the gate
resistance entering both `f_max` denominators was too small by **exactly
`W_RF/W = 40`** — 0.208 Ω where it should be 8.333 Ω.

| peak, 40 µm / 8 fingers, V_ds = 0.1 V | superseded | corrected | factor |
|---|---|---|---|
| f_max intrinsic | 18.786 GHz | **13.094 GHz** | 1.4347 |
| f_max extrinsic (+15 fF/pad) | 3.1885 GHz | **2.2196 GHz** | 1.4365 |
| f_max / f_T there | 0.928 | **0.647** | — |
| f_T, every bias, both geometries | — | — | **bitwise identical** |

`f_T` is untouched because `R_g` does not appear in `f_T = |g_m|/(2π C_gs)`.
The **normalized 1 µm path** — panels 1–2 of `rf_figures_of_merit.png` and the
≈20 GHz / 18.731 GHz figures §4.6 quotes — is **bitwise unchanged, verified
against the pinned pre-fix blob** rather than assumed. §4.6 is annotated in
place with the superseded numbers kept.

**The uncomfortable half:** the corrected model sits at f_max/f_T = 0.65, i.e.
**further** from Feijoo *et al.*'s de-embedded 1.3–1.4, not closer. §4.6 read
as though model and literature had been brought into rough agreement; they have
not been. The gap is at least attributable now: at 8 fingers the corrected
`R_g = 8.33 Ω` is already in that paper's engineered-low range, so the
shortfall is in `g_ds` and the `R_g·C_gd` feedback term, not the gate
resistance. **Not quantified — now an open item.**

### 4. Why it survived five weeks, and the third reason is the general one

**(a) The transcript printed the RIGHT `R_g` beside an `f_max` computed from
the wrong one.** `summary_numbers()` calls `gate_resistance(W=W_RF, ...)` with
`W` **explicit** and printed `R_g (8-finger) = 8.33 Ohm`.
`_compute_fT_fmax_core` called the same function with the default and got
0.208 Ω. The printed diagnostic and the number it was supposed to diagnose came
from two different values of the same quantity, two lines apart.

**(b) Every check in reach was a RATIO, and the defect was multiplicative in
both terms.** intrinsic/extrinsic peak f_max **5.8918 → 5.8994 (0.13 %)**;
extrinsic/intrinsic degradation **0.1697 → 0.1695 (0.1 %)**; both f_max values
individually **≈44 %**. The module's cited sanity check against Feijoo *et al.*
is a ratio. **A ratio cannot see a defect it shares.**

**(c) `MUST_CHANGE` is a test against ZERO.** The figure-provenance audit's
`N_FINGERS_RF x2` mutation acts through `R_g` alone (`R_g ∝ 1/N²`). With `R_g`
frozen 40× too small it barely mattered in the denominator:

```
  max rel change, as committed : 0.01033
  max rel change, corrected    : 0.28620     suppression 27.7x
```

**Both are PASSes.** The knob had lost **96 %** of its authority and still read
as alive. The frozen default corrupted *the audit's measurement of its own
sensitivity*, and the audit had no way to say so.

### 5. Three faults of my own, recorded rather than quietly re-tuned

- **E5's exactness premise was false.** It asserted
  `R_g(after) == R_g(before) * 40.0` exactly because 40 is representable; it
  failed at the last bit (`...332` vs `...334`). Representability of the
  *factor* says nothing about the *rounding sequence*: scaling `W` rounds in
  the numerator, `r0*40` rounds an already-rounded quotient. **An exactness
  claim must name the operation it is exact under.** E5 is now an identity
  about *delivery* (default path == explicit path, bitwise) and **E5b** keeps
  the 1-ULP measurement that exposed it.
- **Section 3's first assertion was simply false** — it claimed `f_T` was
  affected too, and failed on its own false premise. Now a **paired**
  `MUST_NOT_CHANGE`/`MUST_CHANGE`, making the defect's scope an assertion
  instead of prose.
- **`BLOB_PREFIX_RF` was a dead name in the new module**, and
  `graphene_dead_name_sweep.py` reported it `DEAD` on the first run — the same
  fault 10-02 hit with its own `SELF`, caught by the same sweep, in the session
  that created the name. It now has a real reader:
  `historical_claim_is_checkable()` reads the **pinned pre-fix blob**
  (`c0c620ac…`, a blob hash, not a revision expression) out of the object store
  and requires this module's claim *about history* to hold of it, reporting
  `SKIP` rather than `PASS` when the object store cannot answer. **A pinned
  reference with no reader is decoration.**

### 6. Census, and what "dormant" does not mean

**10 cross-module captured defaults this morning, 8 now.** The remaining eight
are `graphene_sensitivity_audit.py`'s `_calib_terms` / `_rho_terms` capturing
`graphene_interconnect_model` constants, and they are **dormant, not faults**:
that audit perturbs by explicit keyword and never rebinds `icm.*`. Calling them
faults would be the 10-01 error of letting a classifier's name drift from what
it measures. **Dormant means no call site rebinds the name *today*** — one
future `setattr` converts all eight at once, so they are enumerated rather than
dismissed. The census also reports that **two `setattr` calls use a computed
name**, so it cannot be complete by construction: said out loud, so that
"0 faults" does not come to mean "none possible".

### 7. The same-session wiring step, done on purpose

10-02's second thread was *a repair that was built and then not connected*, and
its remedy: the last step of a repair is a check that the repair is reachable
from where the fault was. So both new `(code, transcript)` pairs were
registered in `graphene_pristine_transcript_audit.py`'s `SUITES` **in the same
session that created them** —
`graphene_cross_module_delivery_audit.py → cross_module_delivery_audit_output.txt`
and `rf_small_signal_model.py → rf_small_signal_output.txt`. The second is not
an audit, but its transcript carries the `f_T`/`f_max` numbers Chapter 4 quotes
and this session changed them, which is exactly the condition the pristine
audit exists to detect. An unregistered transcript is a claim on disk that
nothing checks.

**Methodological note, continuing the series.** 09-28: prose is a detector, and
a check that cannot fail is worse than no check. 09-29: a mutation that does not
arrive is indistinguishable from a system that does not respond. 09-30: a
control has to sit where the failure enters. 10-01: and it has to name what it
compares against in a way that cannot drift. 10-02: and when it agrees, that is
a fact about two artefacts, not about the world. **10-03: AND A PASS/FAIL AT
ZERO IS SILENT ABOUT MAGNITUDE.** Every `MUST_CHANGE` here measures a response
size, prints it, and throws it away; the verdict retains only the sign, so a
knob can lose 96 % of its authority — or, in a limit this session did not test,
99.9 % — without one check changing colour. **The second strand is narrower and
sharper:** a single-module instrument was asked a cross-module question and
returned *the right answer to the wrong question*. The only reason that was
discoverable is that somebody asked what the instrument's **scope** was rather
than what its **verdict** said. 10-02 said agreement between two artefacts is
silent about the world; 10-03 adds that **agreement between an instrument and
its own specification is silent about whether the specification covers the
case**.

**Validations:** new cross-module delivery audit **25/25**, including six
tolerance-free exact checks (E1 the frozen/live `R_g` ratio is `W_RF/W`
exactly; E2 the no-override path is bitwise unchanged by the fix; E3 explicit
arguments bitwise unchanged; E4 an identity rebind moves `R_g` by exactly 0;
E5 the rebind now arrives, bitwise; E6 the frozen replica is inert under the
same rebind, which is what makes E5 a test of the fix rather than of
arithmetic), seven controls on the new AST detector, a mutation test that
converts the synthetic capture to a live read and requires C1 to stop firing,
and the exact `R_g(N=16) == R_g(N=8)/4` 1/N² law. Dead-name sweep **8/8** (7/8
intermediate, the FAIL real and mine). Figure-provenance audit re-run and
green, with its `N_FINGERS_RF` sensitivity up 27.7× for the reason in §4(c).
`rf_figures_of_merit.png` regenerated; the normalized-device panels verified
bitwise identical to the pre-fix blob.

**Not yet covered (candidates for future runs):**

- **Record the MAGNITUDE of every `MUST_CHANGE`, not just its sign** —
  created today and **the top item**, because §4(c) shows the current form
  passed a knob that was 96 % dead. The mechanisable version is a committed
  per-mutation sensitivity baseline, with a check that fires when a response
  shrinks by more than some factor, not merely when it reaches zero. Note the
  trap this session walked into: the obvious threshold is a magic number, so
  the baseline has to be *measured and committed*, which makes it a transcript
  and puts it under the 10-02 rule
- **Which other checks in this repository are RATIOS of two quantities the same
  defect would scale?** — created today, the direct successor to §4(b). Three
  were found by eye in one module; Chapters 4–6 compare model to literature by
  ratio in several places, and a ratio is blind to any factor common to both
  terms. The general census is mechanisable: every comparison whose two sides
  share a module-level constant in their dependency closure
- **Audit every REMAINING check for the 10-02 §1 fault: does it read the
  artefact its own prose names?** — created 10-02, untouched
- **Do the four existing audits have references that are not content-pinned?**
  — created 10-01, untouched. Today added one correctly-pinned reference
  (`BLOB_PREFIX_RF`) and no census
- **Is `D_CHEM_PHYS` the only duplicated literature VALUE?** — created 10-01,
  untouched; the check is duplicated literal values, not duplicated names
- **Is `plot_fermi_surface` the only uncommitted figure?** — created 10-02,
  untouched
- **Quantify the `g_ds` vs `R_g·C_gd` split in the remaining f_max shortfall**
  — created today, and the first *device-physics* item this series has produced
  in a while: with `R_g` corrected, this model is at f_max/f_T = 0.647 against
  Feijoo *et al.*'s 1.3–1.4, and the two candidate causes are separable by
  zeroing each term in turn
- **Whether Chapter 4's `n_puddle` actually reproduces Section 3.3's measured
  6.45 kΩ/sq floor** — created 09-29, **untouched for five sessions**, still a
  one-line calculation, and still the only place Chapter 3's central negative
  result touches Chapter 4's numbers
- **Whether the three identity re-assignments should be replaced by real
  mutations** — created 09-30, untouched. Today's E4 shows the identity form
  does have one legitimate use (an exact zero-by-symmetry control), which
  narrows the item rather than closing it
- **The Section 3.6 Pauli edge is absent from Chapter 6's model** — created
  09-29; the third leg of the Chapter 4 / Chapter 6 contact-metal contradiction
  and the only one that is a physical mechanism
- **Where does the 6.5430 % internal collection efficiency come from?** —
  created 09-30, untouched
- **Chapters 2 and 3 are drafted but §4.8.2's form of error is not audited for
  elsewhere** — created 10-01; Chapters 5 and 6 contain several
  model-against-literature comparisons whose bias conditions are not stated
- **Whether Chapter 4 or 5 contains a RANKING that is a step artefact** —
  created 09-23/09-25; the noise form closed 10-02, the step form is untouched
- **A finite-temperature optical conductivity** (09-29); **the
  remote-polar-phonon cap on SiO₂ is asserted, not computed** (09-29);
  **angular trigonal warping** (09-28); **finite-temperature carrier density
  `n(E_F, T)`** (09-28, wanted twice); whether other `== 0.0` exactness checks
  are round-trip tautologies (09-28) — **sharpened today**, since E5's first
  form was an exactness claim with a false premise, so the census should ask
  what operation each `== 0.0` is exact *under*; whether RESULT 2's other
  convictions are read from the saturated branch (09-28); a probe reporting
  both directions of lambda (09-28); whether `t'` can be EXCLUDED
  quantitatively as the source of Chapter 4's electron-hole asymmetry (09-28);
  whether any other near-cancellation shows the `c2`/`c4` sign flip if
  differentiated (09-27); migrating every remaining validation to
  measured-value-beside-derived-bound form (09-27); `n_segments = 50` at a bias
  with more curvature (09-26)
- **A second anchor for `Delta_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, now **untouched for thirteen
  consecutive sessions**. Today was again an audit session and did not touch
  it, which is worth saying plainly rather than letting the streak go unnamed
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural weakness
- **Whether Chapter 7's Section 7.3 differential-interconnect prediction
  holds** (09-27); which of Chapters 4 and 5's design rules could be restated
  as parities or bounds (09-23); whether the parity survives a
  photo-thermoelectric term (09-23); re-check whether other "for every …"
  claims rest on small samples (09-20); whether Chapter 4's contact-resistance
  results should be re-run at the 5.4 eV crossover (09-18); reconciling Mueller
  *et al.*'s 0.12 eV step (arXiv:0902.1479) with the 0.25–1.07 eV offsets
  `METAL_WORK_FUNCTIONS` assumes (09-18); a photo-thermoelectric term (Kasirga
  review); Shimomura *et al.*'s comb-electrode design; integrating 6.5's
  plasmonic near-field picture with the spatially-resolved contact-doping
  machinery; isolating the root cause of the Section 4.7 negative residual
  (08-31); Ti and Cr per-metal `Rc` recalibration (ResearchGate rate-limiting);
  a second independent edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

**Automation health.** Device reachable and folder connected at the 04:30
firing. **The device network problem that dominated 10-01 and 10-02 was absent
today:** both repositories cloned fully in under 15 s from inside
`device_bash`, against 10-02's measured 9.3 KB/s and a `git clone` that could
not finish inside the 180 s shell limit at all. The 10-01/10-02 workaround
(clone in the cloud sandbox, ship a bundle to the device, push from a
blobless partial clone) was therefore **not needed and not used**; work was
done entirely in the device VM's scratch space at `$HOME/work`, outside the
connected folder, because git still cannot create its lock files inside a
connected folder. `scipy` was absent from the device VM as expected and
`pip install`ed in 20 s. The `GIT_ASKPASS` recipe ran from the session VM's
own temp space and the token copy was shredded; the token was never written
into `.git/config`, a remote URL, any repository file, or the connected
folder. **One deviation worth stating:** Step 0's "has today's work already
run" check was answered from the freshly cloned `AUTOMATION_LOG.md` and
`git log --since=midnight` only, with no GitHub-API `pushed_at` cross-check,
because the clones finished fast enough that the cheap pre-check 10-02 added
had no latency to hide.

---

## 2026-10-04 — f_max was computed correctly, under the wrong name

**Status.** The device-physics item created 10-03 — *quantify the g_ds vs.
R_g·C_gd split in the remaining f_max shortfall* — is **closed**, and the
answer retires the thread that produced it. New module
`graphene_fmax_shortfall_decomposition.py` (**30/30, 6 of 6 mutants
killed**, five tolerance-free checks), new figure
`fmax_shortfall_decomposition.png`, new note
`notes/2026-10-04-fmax-was-a-resistance-ratio.md`, Chapter 4 **§4.6.2** and
Chapter 7 **§7.8.1a** drafted, two defects of my own fixed, and one latent
crash in a standing guard fixed that had been live since 10-01.

### 1. The split: one of the two named causes carries 0.22 %

The f_max denominator is a sum, so the split needs no perturbation theory.
40 µm / 8 fingers, V_ds = 0.1 V, peak-f_T bias:

| term | value | share |
|---|---|---|
| A = g_ds(R_g+R_s) | 6.033531e-01 | **99.7846 %** |
| — g_ds·R_g | 3.175543e-01 | 52.5182 % |
| — g_ds·R_s | 2.857988e-01 | 47.2664 % |
| B = 2π f_T C_gd R_g | 1.302534e-03 | **0.2154 %** |

Deleting C_gd **entirely** moves f_max/f_T from 0.643007 to 0.643701:
**+0.108 %** against a shortfall of a factor of 2.2. Its share reaches only
1.9908 % at V_ds = 1 V. Chapter 4's annotation named two residual causes and
**one of them is not a cause**; the sentence is struck through in place with
its numbers intact.

### 2. R_g was never the lever, and that is a verdict on five weeks of work

| counterfactual | f_max/f_T |
|---|---|
| as committed | 0.643007 |
| C_gd = 0 | 0.643701 |
| 16 fingers (R_g/4) | 0.827025 |
| 64 fingers (R_g/64) | 0.927229 |
| **R_g = 0 exactly** | **0.935276** |
| R_s = 0 exactly | 0.885467 |
| g_ds ÷ 4.8483 | 1.410000 |

A **perfect gate** lands 28.06 % below the bottom of Feijoo *et al.*'s
1.3–1.4 band, and 64 fingers is within 1 % of that ceiling, so the
multi-finger mechanism is exhausted at about sixteen. The 10-03 session
found a real 40× delivery bug in R_g and was right to. **The belief that
R_g was the lever did not survive measurement**, and no check in this
repository could have said so: every one asked whether R_g was computed
correctly, none asked whether R_g mattered.

### 3. Why, and the number the gap is now worth

`transfer_characteristic()` computes `Id = Vds/(R_channel + Rc_total)` — a
bias-dependent resistor, no velocity saturation, no pinch-off. So g_ds **is**
the channel conductance: measured `g_ds·R_total` ∈ [1.0026, 1.0487] across
V_ds = 0.05–1 V. Term A therefore collapses to a resistance ratio and
`f_max/f_T ≈ ½√(R_total/(R_g+R_s))` reproduces the exact value to **0.36 %**
(0.645352 vs 0.643007). Reaching 1.41 needs R_total/(R_g+R_s) ≥ 7.95; this
device has 1.67. **f_max is primarily a measurement of output resistance and
this model has none to measure.** The g_ds that would reach the band at the
model's own R_g, R_s, C_gd is 0.1965 mS/µm against 0.9527 mS/µm — the missing
physics is worth **×4.8483**, so "add saturation" is now a target with a
tolerance rather than a direction.

**Independent corroboration of the mechanism, found after the derivation.**
The Chalmers doctoral thesis on GFETs for high-frequency applications states
it outright: *"f_max is mostly limited by the high drain conductance g_ds due
to the lack of current saturation"*, with the reason graphene cannot pinch
off (no bandgap, so the carrier type inverts along the channel) and a second,
non-intrinsic cause in interface-state emission at high fields. Saturation
velocity there is 0.2–0.8 v_F theoretical, ≈1.4×10⁷ cm/s measured. **This
corroborates the mechanism, not the number** — and the distinction is kept,
because R5's 0.2 threshold is set from this model's own range, making it a
drift detector rather than a comparison against the world (its label now
says so).

### 4. A withdrawn inference about f_T, costing more than it looks

Peak f_T scales **×19.8 over a ×20 drain-bias range** (10.140 GHz at
V_ds = 0.05 V to 200.440 GHz at 1 V), because g_m ∝ V_ds in a resistor model.
So the ≈20 GHz that §4.6 called "consistent with the literature range for
non-exotic gate lengths" is a fact about the bias `plot_fT_fmax()` passes;
the same model matches the 200 GHz record devices at V_ds = 1 V. **No §4.6
number is retracted. What is retracted is reading their proximity to the
literature as corroboration.** A model with one knob worth ×20, compared
against a literature range spanning ×20, will agree somewhere — and this is
the 10-02 rule in its most expensive form yet, because the second artefact
was the published literature rather than another file here. Free diagnostic
from the same table: f_max/f_T **falls** with bias here (0.6442 → 0.6218)
where a saturating device's rises, so the sign is wrong, not just the
magnitude.

### 5. Two faults of my own, both caught by this repository's own guards

- **The module written about frozen defaults contained four.**
  `decompose()`'s first committed form captured `VG_SWEEP`, `VDS_REF`,
  `rf.W_RF` and `rf.N_FINGERS_RF` as function-parameter defaults; the last
  two are the **exact cross-module shape of the 10-03 bug**, committed by the
  session that had just read the whole note about it. Caught by
  `graphene_mutation_arrival_probe.py` §7b on the first run after the module
  landed — **the first time that census has had anything to report since
  10-01**. All four late-bound; the sentinel is the builtin `Ellipsis`, not a
  module-level `_DEFAULT = object()`, because a module-level sentinel is
  itself a captured module-level name and §7b reports it correctly by its own
  rule. The rule was left as written and the code moved (10-01: a
  classifier's name must not drift from what it measures). Every transcript
  number is **bitwise unchanged** by the fix, because nothing here rebinds
  `rf.W_RF` between import and call — *latent*, which is what 09-29 said
  about ten names, one of which became the 40× bug eleven sessions later.
- **The guard detected it and then crashed.** §7b's renderer unpacked **five**
  fields from rows the census has only ever built with **four**, so any
  non-empty census raised `IndexError` on its first element and the module
  exited 1 before Sections 8 onward. It had never fired, because the census
  has returned zero rows every run since 10-01: **a reporting path that only
  runs when the guard fires had never run.** The operator saw a stack trace
  where four named sites were the intended output. A crash is a detection but
  a strictly weaker one — no diagnosis, every later section lost, resting on
  an arithmetic accident downstream rather than on an assertion. Fixed as
  `_print_census()`, and **C7b now renders a synthetic three-row census
  (FROZEN, MIXED, IMPORT_FAILED) on every run**, so the path is covered by the
  suite rather than by the day the guard fires; without it the fix would be
  untested until the next captured default, which on this week's evidence is
  three days to never. Probe **33/33 → 34/34**, exit 1 → 0.

### 6. Validations

New module **30/30, 6/6 mutants killed**, with control A (unmutated battery
green) and control B (an identity rebind is *not* caught, so "all mutants
die" is a claim about the mutants and not about a battery that fails on
anything). Five tolerance-free checks, each naming the operation it is exact
under per 10-03: **E1** bitwise reconstruction of the model's own f_max,
exact under *summation order* — re-associating the identical algebra breaks
it, and **M1 confirms that** rather than leaving it asserted; **E2/E2b/E2c**
multiplication by zero; **E2d** every parasitic deleted gives a denominator
of exactly 0.0, then the model's own 1e-30 clip, reported as the guard
constant it is rather than printed as a counterfactual f_max/f_T of 5e14;
**E3/E3b/E3c** the 1/N² finger law at powers of two, with R_s exactly
invariant; **E5/E5b** the **overlapping-limit validation** — this module's
default-geometry path must agree with `rf.compute_fT_fmax` bitwise, which is
the superseded-model comparison 10-03 required. **E4 is a round trip and is
checked to 1 ULP rather than claimed exact**, because sqrt is involved.
**Magnitude-form controls**, per the top 10-03 item: C4 asserts a *band*
(deleting C_gd moves the answer by between 0 and 0.5 %, i.e. the knob is
alive and its authority is negligible), C5 asserts the R_g = 0 limit still
falls short, C6 requires the two responses not to be aliased (they differ
421.32×). Repository-wide: pristine transcript audit **10/10** with the new
transcript **IDENTICAL** in a pristine clone at 20.4 s; dead-name sweep 8/8;
figure provenance 26/26; cross-module delivery 25/25.

**Worth reading off the mutation report.** M4 — the 10-03 frozen-default bug
itself, 40× in R_g, reinstated — is killed by **exactly one** check, R7, and
R7 is the only check in the battery that retains a **magnitude**. Every
PASS/FAIL-at-zero check is blind to it, because shrinking R_g moves this
module's conclusion in the direction it already argues. **That is 10-03's top
open item reproduced from the other side, on the very bug that created it.**

### 7. What the top open item got, and did not get

The 10-03 top item asks for a *committed repository-wide per-mutation
sensitivity baseline*. It did **not** get one. What it got is one module
written in the form it asks for — every MUST_CHANGE prints its response size,
and the two checks whose point is a magnitude assert bands — plus the M4
evidence above that the form is the right one. **Said plainly so the item is
not quietly downgraded**: it remains the top methodological item, now with a
worked example and an argument, which is less than a mechanism.

**Methodological note, continuing the series.** 09-28: prose is a detector.
09-29: a mutation that does not arrive is indistinguishable from a system
that does not respond. 09-30: a control has to sit where the failure enters.
10-01: and name what it compares against in a way that cannot drift. 10-02:
when it agrees, that is a fact about two artefacts, not about the world.
10-03: a PASS/FAIL at zero is silent about magnitude. **10-04: AND A QUANTITY
CAN BE COMPUTED CORRECTLY UNDER THE WRONG NAME.** Every check here asked
whether f_max was computed correctly from g_ds, R_g, R_s and C_gd. It was —
bitwise, reproducibly, with the 40× bug found and fixed. None asked **what
the number is a measurement of**, and the answer was the device's
access-to-channel resistance ratio. There is no tolerance that detects this,
no mutation that kills it, and no transcript that goes stale: the arithmetic
was right the whole time. The only thing that found it was asking the
denominator which of its terms mattered — a question about **attribution**,
not about **correctness**. 10-03's corollary was that an instrument's
agreement with its own specification is silent about whether the
specification covers the case; **today's is that an instrument's correctness
is silent about whether it measures the quantity its name claims** — and, in
the form that is about effort rather than error, that the R_g thread was
*correct work aimed at a non-binding constraint*, which auditing could never
have told anyone, because auditing asks whether a number is right and never
whether it is worth having.

**Not yet covered (candidates for future runs):**

- **A saturation term in `transfer_characteristic()`** — created today and
  **the new top item**, and the first open item in this series that arrives
  with a measured acceptance criterion: it must divide g_ds by ≈4.85 at the
  peak-f_T bias AND reverse the sign of d(f_max/f_T)/dV_ds, which currently
  runs the wrong way (0.6442 → 0.6218 over V_ds = 0.05–1 V). The literature
  route is a saturation velocity v_sat in the 0.2–0.8 v_F range with the
  I_d = W·q·n·v_sat ceiling, or the drain-field-dependent mobility of the
  compact large-signal GFET models (arXiv:1605.08235). **Validate it against
  the resistor model in the low-V_ds overlapping limit**, where the two must
  agree, per the rule 10-03 set and E5 implements
- **Record the MAGNITUDE of every `MUST_CHANGE`, not just its sign** —
  created 10-03, **still the top methodological item**, advanced by one
  worked example and not mechanised (§7 above). Today's M4 result is new
  evidence for it: the only detector that caught a historical 40× bug is the
  only one in its battery that keeps a magnitude
- **Which other checks here are RATIOS of two quantities the same defect
  would scale?** — created 10-03, untouched. Today adds a second form of the
  same blindness worth folding in: not only ratios, but any check whose
  *name* asserts more than its *arithmetic* tests
- **Which other quantities in this repository are computed correctly under a
  name that claims more than they measure?** — created today, the direct
  successor to §Methodological note. The f_max case was found by asking a
  denominator which term dominated; the mechanisable version is a census of
  every reported figure of merit whose value is dominated (>95 %) by a single
  term, since such a quantity is a restatement of that term. Candidates by
  construction: Chapter 5's interconnect figures of merit, Chapter 6's
  responsivity–gain tradeoff
- **Audit every REMAINING check for the 10-02 §1 fault: does it read the
  artefact its own prose names?** — created 10-02, untouched
- **Do the four existing audits have references that are not content-pinned?**
  — created 10-01, untouched
- **Is `D_CHEM_PHYS` the only duplicated literature VALUE?** — created 10-01,
  untouched
- **Is `plot_fermi_surface` the only uncommitted figure?** — created 10-02,
  untouched
- **Whether Chapter 4's `n_puddle` actually reproduces Section 3.3's measured
  6.45 kΩ/sq floor** — created 09-29, **untouched for six sessions**, still a
  one-line calculation, and still the only place Chapter 3's central negative
  result touches Chapter 4's numbers. Today's work makes it *more* relevant,
  not less: `n_puddle` sets R_channel, and R_channel is now known to be half
  of what f_max measures
- **Whether the three identity re-assignments should be replaced by real
  mutations** — created 09-30; narrowed 10-03, and narrowed again today,
  since control B here is a deliberate identity rebind with a stated purpose
- **The Section 3.6 Pauli edge is absent from Chapter 6's model** — created
  09-29; the third leg of the Chapter 4 / Chapter 6 contact-metal
  contradiction and the only one that is a physical mechanism
- **Where does the 6.5430 % internal collection efficiency come from?** —
  created 09-30, untouched
- **Chapters 2 and 3 are drafted but §4.8.2's form of error is not audited
  for elsewhere** — created 10-01, untouched
- **Whether Chapter 4 or 5 contains a RANKING that is a step artefact** —
  created 09-23/09-25; the noise form closed 10-02, the step form untouched
- **A finite-temperature optical conductivity** (09-29); **the
  remote-polar-phonon cap on SiO₂ is asserted, not computed** (09-29) —
  **today makes this one sharper**, since a remote-polar-phonon-limited
  velocity is one of the two standard sources of the saturation term the new
  top item needs; **angular trigonal warping** (09-28); **finite-temperature
  carrier density `n(E_F, T)`** (09-28, wanted twice); whether other `== 0.0`
  exactness checks are round-trip tautologies (09-28, sharpened 10-03 to ask
  what operation each is exact *under* — today's five exactness checks each
  name theirs, so the form exists and the census does not); whether RESULT 2's
  other convictions are read from the saturated branch (09-28); a probe
  reporting both directions of lambda (09-28); whether `t'` can be EXCLUDED
  quantitatively as the source of Chapter 4's electron-hole asymmetry
  (09-28); whether any other near-cancellation shows the `c2`/`c4` sign flip
  if differentiated (09-27); migrating every remaining validation to
  measured-value-beside-derived-bound form (09-27); `n_segments = 50` at a
  bias with more curvature (09-26) — **now load-bearing**, because g_ds is a
  V_ds difference of that quadrature and g_ds is 99.78 % of the f_max
  denominator
- **A second anchor for `Delta_c`, at any separation other than 3.3 Å** —
  open since 2026-09-21, still the top *physics* item, now **untouched for
  fourteen consecutive sessions**. Today was a device-physics session rather
  than an audit session and still did not touch it, which is the more
  pointed version of the streak
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural
  weakness
- **Whether Chapter 7's Section 7.3 differential-interconnect prediction
  holds** (09-27); which of Chapters 4 and 5's design rules could be restated
  as parities or bounds (09-23); whether the parity survives a
  photo-thermoelectric term (09-23); re-check whether other "for every …"
  claims rest on small samples (09-20); whether Chapter 4's
  contact-resistance results should be re-run at the 5.4 eV crossover
  (09-18); reconciling Mueller *et al.*'s 0.12 eV step (arXiv:0902.1479) with
  the 0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS` assumes (09-18); a
  photo-thermoelectric term (Kasirga review); Shimomura *et al.*'s
  comb-electrode design; integrating 6.5's plasmonic near-field picture with
  the spatially-resolved contact-doping machinery; isolating the root cause of
  the Section 4.7 negative residual (08-31); Ti and Cr per-metal `Rc`
  recalibration (ResearchGate rate-limiting); a second independent
  edge-contact dataset (Lee *et al.* 2022, Wiley 403'd)

**Automation health.** Device reachable and folder connected at the first
firing checked today; Step 0 found neither repository carrying a 2026-10-04
entry and no commits since midnight, so this was a full session rather than
10-03's one-repo branch. Both repositories cloned fully in well under 15 s,
so the network problem of 10-01/10-02 remains absent for a second
consecutive session — two clean sessions still do not establish it is gone,
and the bundle-via-`_to_push` recipe stays documented. `scipy` was absent
from the device VM as expected and `pip install`ed in ~20 s. Work was done in
the device VM's scratch space at `$HOME/work`, outside the connected folder,
because git still cannot create its lock files inside one. **Live web search
was available and used**, from inside the device VM, and the Chalmers
full-text PDF fetched successfully; no fetch failures to record, and PubMed/PMC
were not attempted, per the standing note that they return reCAPTCHA here.
The `GIT_ASKPASS` recipe ran from `$HOME/.sess/`, the token was never written
into `.git/config`, a remote URL, any repository file or the connected
folder, and the temp copy is shredded at the end of the session.
`user.name`/`user.email` were again absent in the fresh clone and set per
09-24. **Push-as-you-go was followed**: the first five commits were pushed and
verified against the GitHub API before the thesis and log work began. **One
deviation worth stating:** `log_sensitivity_step_audit.png` was re-rendered as
a side effect of re-running that audit for its transcript, differing in bytes
without differing in content, and was restored to its committed bytes rather
than committed — the committed figure stays the pinned artefact instead of
being churned by an unrelated run.

**Commits this run:** ~~8~~ **9** (the module with its figure and transcript;
the late-binding fix with the SUITES registration; the §7b renderer fix with
C7b; the four censuses the module's existence moved; the pristine-audit
transcript; Chapters 4 and 7; **the notes file**; the R5 label correction with
its own transcript; the cross-module transcript). This AUTOMATION_LOG.md entry
makes **10**. **Corrected in place after verifying against the GitHub API
rather than against this list:** the first form said 8 and 9, having
enumerated eight items and omitted the notes commit from its own list while
counting the rest correctly. It is a trivial error and it is corrected rather
than overwritten, because the whole content of the 10-01 entry is that a
number quoted from memory and a number read from the artefact are different
things, and the API answer was ten.

## 2026-10-05 — The saturation term arrived, works, and is not the binding constraint

**Status.** The 2026-10-04 top item — *a saturation term in
`transfer_characteristic()`* — is **opened and reported against all three of
its acceptance criteria: B met, C met, A missed by 4.3×.** New module
`graphene_velocity_saturation_model.py` (**26/26**, eight tolerance-free
checks), new mutation harness (**7 of 7** after a first run of 5 of 7 that
found two real holes), new figure `velocity_saturation_model.png`, new note,
Chapter 4 **§4.6.3** and Chapter 7 **§7.8.1b** drafted. Two long-standing
items close, one of them with a negative result; two predictions of mine were
disproved by my own decomposition and are recorded rather than edited out.

### 1. The item was the best this repository has written, because it had a number in it

2026-10-04 did not ask for velocity saturation. It asked for a term that
**(A)** divides g_ds by ≈4.85 at the peak-f_T bias, **(B)** reverses the sign
of d(f_max/f_T)/dV_ds, and **(C)** validates against the resistor model in the
low-V_ds overlapping limit. That is a criterion that can be *missed*, and it
was — which is the whole of today's result. Stated plainly, because the
temptation to report this as a success is real: **the term works, the sign is
fixed, and the item is not closed.**

### 2. The model, and why v_sat was not allowed to be a free parameter

Caughey–Thomas soft saturation at β = 1 (Feijoo *et al.*'s γ = 1 form). Current
continuity makes the channel-length integral separable, so I_d stays
closed-form even with v_sat varying along the channel:

    I_d = μWe·Q / (L + μ·S),   Q = ∫n dV,   S = ∫dV/v_sat

closed on the contacts by V_ds,ch = V_ds − I_d·R_c. S = 0 recovers
drift-diffusion **along the same code path**, so criterion C compares a model
against its own limit rather than against a second implementation. Added
*alongside* `transfer_characteristic()`, which is checked bitwise unchanged
(G4): five modules and every committed transcript read it, and rewiring is a
2026-10-02 rule change that should not ride along in someone else's commit.

Criterion A names a number and v_sat is the only knob that moves it, so
**fitting v_sat to the criterion would have been 2026-10-02's fault in its
purest available form.** v_sat comes instead from optical-phonon emission,
(2/π)Ω/√(πn), with ħΩ = 0.10 eV taken as published. At this device's density
that gives 6.63e7 cm/s — *above* Dorgan *et al.*'s measured 1–3e7 cm/s band on
SiO₂, i.e. the generous end. β = 1 is generous too: Dorgan's own best fit is
β = 2, which saturates *less* hard (24.6–26.0 % higher velocity at the largest
field here). The shortfall below is bounded from both sides.

### 3. Criterion B — met, and it is the clean win

| V_ds | resistor | saturated |
|---|---|---|
| 0.05 | 0.644461 | 0.663187 |
| 0.10 | 0.642966 | 0.683186 |
| 0.20 | 0.640511 | 0.728061 |
| 0.50 | 0.633440 | 0.911727 |
| 1.00 | 0.621775 | **1.384374** |

Slope −0.0227 → **+0.7212**, a reversal 31.8× the magnitude of what it
replaces. **The V_ds = 1 V row lands in Feijoo's 1.3–1.4 band and is
explicitly NOT offered as corroboration**, for the reason 2026-10-04 recorded:
a ×20 knob against a ×20 literature band agrees somewhere. The claim is the
sign.

### 4. Criterion A — missed at 23.4 %, and the miss locates the constraint

| | resistor | saturated | factor |
|---|---|---|---|
| g_ds at the resistor's peak-f_T bias | 3.811e-2 S | 3.355e-2 S | **1.1361** |
| g_ds at the saturated model's peak | 3.767e-2 S | 3.354e-2 S | 1.1230 |
| required | | | **4.8483** |

The term is not weak — μS/L = 0.3036 would divide the *channel* conductance by
1.699 — it is **diluted**. At the 40 µm geometry R_c,total is 15.0 Ω of a
≈26 Ω device, so about half of g_ds is a contact resistance no channel
mechanism can reach. **The contacts, not the channel, now cap g_ds.**

Inverted, so the gap is falsifiable rather than rhetorical: criterion A needs
v_sat = **6.50e6 cm/s, below the measured band**, and a sweep over the four
available phonon energies (0.059 / 0.100 / 0.149 / 0.196 eV) reaches **1.2491
at best**. No phonon energy in the physical range satisfies it. That table also
discharges the 2026-09-29 item *the remote-polar-phonon cap on SiO₂ is
asserted, not computed*, which 2026-10-04 predicted this work would sharpen:
Eq. (6) is linear in Ω, so the cap's authority over g_ds is now a measured
1.07–1.25.

**Three levers, three non-binding results.** R_g exhausted (R_g = 0 exactly
still falls 28 % short, 10-04). C_gd negligible (deleted entirely: 0.11 %,
10-04). v_sat diluted by the contacts (today). Each was proposed in turn as the
missing piece and **all three point at the access resistance** — §4.5, §4.7,
§4.8, the part of this thesis with the most content. Chapter 7 §7.8.1b records
what that does to the RF verdict: §7.8.1a's reason was a fact about the
*instrument* (it measured a resistance ratio); the reason now available is a
fact about the *device* (the figures of merit are contact-limited). **That is a
claim this thesis is entitled to make, because the mechanism was built and
measured rather than assumed.**

### 5. Criterion C — met, and it disproved two of my own predictions

Four artefacts made to converge as V_ds → 0, separating three mechanisms whose
different orders in V_ds are what make them separable at all:

| mechanism | at V_ds = 0.1 V | order |
|---|---|---|
| quadrature rule (50-sample mean vs. trapezoid) | −4.3e-05 % | 2.04 |
| **profile domain** | **+0.253555 %** | **1.00** |
| Jensen gap ⟨1/n⟩ vs. 1/⟨n⟩ | +3.8e-04 % | 2.00 |

**(i) I predicted `np.mean` over `np.linspace(0, Vds, 50)` carried an O(1/N)
sample-mean bias. It does not** — a mean of equally spaced samples including
both endpoints is exact for a linear integrand and second order for a smooth
one. `n_segments = 50`, open since 2026-09-26 and called **load-bearing** on
2026-10-04 because g_ds is a V_ds difference of exactly that quadrature, is
hereby measured at 4.3e-05 % and **closes with a negative result after nine
days.**

**(ii) I predicted the residual was the Jensen gap. It is 674× not.**
`transfer_characteristic()` profiles the channel over the **full V_ds**, but
over half of V_ds is dropped across the contacts and never appears across the
channel, so it evaluates n(V_ch) over a range about twice too wide. First order
in V_ds, and the real content of criterion C.

**No committed number is withdrawn.** The total is one-signed, 0.2539 % at
V_ds = 0.1 V, and makes I_d *larger* — so every shortfall against literature in
Chapter 4 is if anything understated. A single aggregate residual would have
shown one number and had it attributed to whichever cause was in mind; the
attribution is only available because the three orders differ.

### 6. The mutation harness found two real holes, and M6 is the headline

First run **5 of 7**. Both survivors were genuine.

- **M6 survived.** It multiplies v_sat by π/2 and does nothing else: no
  structure, no exception, no sign flip, 57 % on the one knob that decides
  whether criterion A was measured or fitted. **Every check in the battery as
  first committed was a sign, a zero or a ratio — all invariant under a
  constant rescale of v_sat.** The standing top methodological item landing on
  the most load-bearing quantity in the item being closed.
- **M7 survived.** It clips the reported physical ceiling to 0.5 instead of
  reporting it; X6 asserted only `max < 1.0`, which a guard that lies about its
  own report satisfies exactly as well as one that does not. Same family as
  10-04's census renderer that detected and then crashed.

Three checks added — **X7** (a magnitude pin on v_sat against an anchor
computed outside the module at 30 significant digits, so it is not Eq. (6)
agreeing with itself), **X6b** (the reported ceiling must equal an independent
recomputation), **X6c** (it must respond across the ladder). Battery **23/23 →
26/26**, harness **7 of 7**, control B still survives. **X7's own first form
used a 6-digit hand value of Ω and FAILED at rel. err 1.9e-04** — it caught my
arithmetic before it caught any mutant, and that is recorded in the module.

**Read M6's killer column.** M6 dies to **exactly one** check, X7, and X7 is
the only check in the battery that retains a magnitude. **That is the identical
shape 2026-10-04 reported for M4/R7, now reproduced independently on a
different module, a different quantity and a different defect.** Two
independent instances is no longer an anecdote.

### 7. Validations

**26/26.** Eight tolerance-free checks, each naming the operation it is exact
under per 10-03: **X1** multiplication by zero; **X2** the solver against the
closed form, reported in ULP (47) rather than claimed exact; **X3** with L = 0
the ceiling I_d = Wenv_sat holds to 1 ULP independent of μ and V_ds, exact
under *cancellation* — and X3's first form drove μ → large at finite L and
compared in ULP, which was the **wrong instrument** for a quantity that
converges as L·v_sat/(μV_ds) and failed on its own tolerance rather than on the
model; **X3b** the finite-L law 1/(1 + L·v_sat/(μV_ds)) to 1 ULP; **X4**
doubling v_sat doubles I_d **bitwise**; **X5** constant n gives Q = n·V_ds,ch at
0 ULP; **X6b** the reported ceiling against an independent recomputation at
rel. err 0.00e+00; **X7** the v_sat magnitude pin at rel. err 5.3e-13.
**X6** reports the physical ceiling of Eq. (3) (max u/v_sat = 0.721) rather
than clipping it. **Magnitude-form controls** per the standing top item: G1
asserts a band, and **its own first lower bound (1.5) was wrong and is recorded
as such in the transcript** rather than quietly relaxed; G3 halving v_sat and
halving μ are not aliased (ratio 0.8112); G6 the two models' f_max/f_T differ
in monotonicity, so saturation is not a relabelled mobility. Repository-wide:
dead-name sweep 8/8, figure provenance audit exit 0.

**Two crashes of my own, both in reporting paths.** `fT_fmax_saturated()`
raised `IndexError` inside `np.gradient` on a single-bias array, twice — the
10-04 shape (a crash is a detection, but a strictly weaker one) reproduced by
the session that had just read the entry about it. Fixed by `gds_at_bias()`,
which needs no gradient.

**Methodological note, continuing the series.** 09-28: prose is a detector.
09-29: a mutation that does not arrive is indistinguishable from a system that
does not respond. 09-30: a control has to sit where the failure enters. 10-01:
and name what it compares against in a way that cannot drift. 10-02: when it
agrees, that is a fact about two artefacts, not the world. 10-03: a PASS/FAIL
at zero is silent about magnitude. 10-04: and a quantity can be computed
correctly under the wrong name. **10-05: AND A CRITERION WITH A NUMBER IN IT
CAN BE MISSED IN A WAY THAT LOCATES THE REAL CONSTRAINT, WHICH A CRITERION
WITHOUT ONE CANNOT.** Had 10-04 written "add velocity saturation", today would
have added it, watched the sign reverse, noticed the V_ds = 1 V row meets
Feijoo's band, and closed the item — every step defensible, the conclusion
wrong, because the quantity that would have gone unmeasured is the 1.1361, and
the 1.1361 is what says the contacts are binding. The corollary is about where
such numbers come from: 4.8483 was not available until f_max had been
*decomposed* rather than computed. **The practice that makes items falsifiable
is the same practice that found 10-04's result: attribution, not correctness.
An item that can only be *done* is an item written before the attribution
was.**

**Not yet covered (candidates for future runs):**

- **f_max/f_T at R_c = 0 EXACTLY, in the saturated model** — created today and
  **the new top item**, named in Chapter 7 §7.8.1b as the single counterfactual
  that would overturn today's annotation. Exactly parallel to 10-04's R_g = 0
  and R_s = 0 rows. If a perfect contact also falls short of 1.3, four levers
  are exhausted and §7.8.1a's structural verdict is final; if it reaches the
  band, the RF case becomes a contact-engineering problem with a quantified
  target. **Either outcome is more useful than what this thesis currently
  carries, and it is a one-counterfactual calculation.**
- **MECHANISE the magnitude item** — created 10-03, worked example 10-04,
  **second independent instance today**, still not mechanised. Two modules, two
  quantities, two defects, same result: the only detector that catches a
  pure-magnitude change is the only check that keeps a magnitude. The
  mechanisable form is a repository-wide census of every MUST_CHANGE recording
  its response size. **Today is the last session at which "one more worked
  example" is a defensible answer.**
- **Rewire `transfer_characteristic()` to the saturated form, or decide not
  to** — created today. Now priced: five modules read it, their transcripts
  move, and §5's profile-domain correction moves them a further one-signed
  0.25 %. **Recommended yes, as its own session**, because it is a 10-02 rule
  change and nothing else should share the commit.
- **Why does the peak-f_T bias MOVE BRANCH between the two models?** — created
  today. −1.1128 V (hole) to +2.8571 V (electron). Both reported, neither
  explained. Either a real electron–hole asymmetry — which Chapter 4 has
  independent reason to expect — or an artefact of where g_m peaks. The `t'`
  item (09-28) is the natural place to look.
- **A β = 2 implementation** — created today; destroys the separability of
  Eq. (4) and needs quadrature on x rather than V. Bounded today at 24.6–26.0 %
  and shown to move the answer the *wrong* way, so completeness, not a lever.
- **Whether Chapter 4's `n_puddle` reproduces Section 3.3's 6.45 kΩ/sq floor**
  — created 09-29, **untouched for seven sessions**, and today raises it again:
  n_puddle sets n at the Dirac point, n sets v_sat through Eq. (6), and the
  resistor model's peak-f_T bias sits near the Dirac point.
- **A second anchor for `Delta_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, **untouched for fifteen
  consecutive sessions**. Today was a device-physics session and again did not
  touch it.
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural
  weakness, and §7.8.1b makes it a limitation on the RF case too.
- **Which other quantities here are computed correctly under a name that claims
  more than they measure?** (10-04); **which checks are RATIOS of two quantities
  the same defect would scale?** (10-03); **audit every remaining check for the
  10-02 §1 fault** (10-02); **do the four existing audits have references that
  are not content-pinned?** (10-01); **is `D_CHEM_PHYS` the only duplicated
  literature VALUE?** (10-01); **is `plot_fermi_surface` the only uncommitted
  figure?** (10-02); **the Section 3.6 Pauli edge is absent from Chapter 6's
  model** (09-29); **where does the 6.5430 % internal collection efficiency come
  from?** (09-30); **whether the three identity re-assignments should be real
  mutations** (09-30, narrowed again today since control B is a deliberate
  identity rebind with a stated purpose); **Chapters 2 and 3 are drafted but
  §4.8.2's form of error is not audited for elsewhere** (10-01); whether
  Chapter 4 or 5 contains a RANKING that is a step artefact (09-23/09-25, step
  form untouched); a finite-temperature optical conductivity (09-29); angular
  trigonal warping (09-28); finite-temperature `n(E_F, T)` (09-28); whether
  other `== 0.0` checks are round-trip tautologies (09-28); whether RESULT 2's
  other convictions are read from the saturated branch (09-28); a probe
  reporting both directions of lambda (09-28); whether `t'` can be EXCLUDED
  quantitatively as the source of Chapter 4's electron–hole asymmetry (09-28);
  whether any other near-cancellation shows the `c2`/`c4` sign flip (09-27);
  migrating every remaining validation to measured-value-beside-derived-bound
  form (09-27); whether §7.3's differential-interconnect prediction holds
  (09-27); which of Chapters 4 and 5's design rules could be restated as
  parities or bounds (09-23); whether the parity survives a photo-thermoelectric
  term (09-23); re-check whether other "for every …" claims rest on small
  samples (09-20); whether Chapter 4's contact-resistance results should be
  re-run at the 5.4 eV crossover (09-18); reconciling Mueller *et al.*'s 0.12 eV
  step with the 0.25–1.07 eV offsets `METAL_WORK_FUNCTIONS` assumes (09-18); a
  photo-thermoelectric term (Kasirga review); Shimomura *et al.*'s comb-electrode
  design; integrating 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery; isolating the root cause of the
  §4.7 negative residual (08-31); Ti and Cr per-metal `Rc` recalibration
  (ResearchGate rate-limiting); a second independent edge-contact dataset
  (Lee *et al.* 2022, Wiley 403'd)

**Automation health.** Device reachable and folder connected at the 04:30
firing; Step 0 found neither repository carrying a 2026-10-05 entry and no
commits since midnight, so this was a full session. **The network problem of
10-01/10-02 is back, and worse than either of those days.** Sustained
throughput from the device VM measured **13.9–18.8 kB/s** against sub-second
API latency, so small requests are fine and bulk transfer is not: a full clone,
a `--depth 1` clone and a `--filter=blob:none --depth 1` clone all exceeded the
180 s `device_bash` ceiling, and background processes do not survive between
calls, so a clone cannot be split across them. Two clean sessions did not
establish the problem was gone, and the recipe that was kept documented for
exactly this case was needed today.

**The route that worked, recorded for the next occurrence.** The cloud
container clones both repositories in ~2 s, so the work was done there; but its
git proxy **refuses to push** to these repositories (`access denied by the git
proxy: … is not in this session's authorized repository set`, HTTP 403), which
is a hard constraint and not a credential problem — the token was never the
issue. The working route is therefore **neither shell alone**: build the
commits in the container, `git bundle create <repo>.bundle main` (11.9 MB for
graphene, full history), move it to the connected folder with
`device_commit_files` — which goes over the desktop bridge rather than the
proxy and completed **effectively instantly** — then on the device VM
`git clone <bundle>`, set the GitHub remote, and push. **The push itself took
5.5 s**, because it transfers only the new objects. The bundle leaves
`_to_push/` populated as a side effect, which is the documented fallback
artefact anyway. Note `git clone` from a bundle leaves HEAD unresolved and
lands on an empty `master`; `git checkout -B main refs/remotes/origin/main`
fixes it.

The `GIT_ASKPASS` recipe ran from `$HOME/.sess/` on the device VM; the token
was never written into `.git/config`, a remote URL, any repository file or the
connected folder, and both the device-VM copy and the transient container copy
are shredded at the end of the session. `user.name`/`user.email` were again
absent in the fresh clone and set per 09-24. **Push verified against the GitHub
API rather than git's own output**, per 09-26. **Push-as-you-go was followed**:
the first two commits were pushed and API-verified before the notes and thesis
work began, and `git fetch origin main` was run immediately before the push per
the 10-03 race item (`origin/main` unmoved, clean fast-forward). **Live web
search was available and used**, from the cloud container; both literature PDFs
fetched successfully and PubMed/PMC were not attempted per the standing note.
`scipy` was present in the container and no install was needed. One deviation
to state: `figure_provenance_audit.png` re-rendered as a side effect of running
that audit and was **restored to its committed bytes rather than committed**,
per the 10-04 precedent.

**Commits this run:** 4 (the model with its transcript and figure; the mutation
harness with its report and the three checks it forced; the session note;
Chapter 4 §4.6.3 with Chapter 7 §7.8.1b). This AUTOMATION_LOG.md entry makes
**5**.

---

## 2026-10-06 — The counterfactual answered both branches of its own dichotomy, because the question was a ratio

**Status.** The 2026-10-05 top item — *`f_max/f_T` at `R_c = 0` exactly, in
the saturated model* — is **closed, and its two-branch dichotomy turned out
not to be exhaustive: both branches are true at once.** New module
`graphene_perfect_contact_counterfactual.py` (**23/23**, five tolerance-free
checks), new mutation harness (**7 of 7** after a first run of 3 of 7 that
found three real holes, all the same fault), new figure, new note, Chapter 4
**§4.6.4** and Chapter 7 **§7.8.1c** drafted. One prediction of 2026-10-05
is confirmed to 0.45 % from the other side of the mechanism it explained.
One literature paper **cuts against the premise of the whole 10-04/10-05
thread** and is recorded in place rather than absorbed.

### 1. The answer, and why the item could not see it

The item was well-formed by every standard this series has accumulated — a
named quantity, a number, both branches of the decision stated in advance:

> If a perfect contact also falls short of 1.3, four levers are exhausted and
> §7.8.1a's structural verdict is final; if it reaches the band, the RF case
> becomes a contact-engineering problem with a quantified target.

Measured, at the 40 µm / 8-finger geometry and `V_ds = 0.1 V`, fixed bias:

| quantity | baseline | `R_c = 0` coherent | change |
|---|---|---|---|
| `f_max/f_T` | 0.683186 | 0.759363 | **+11.15 %** (58.41 % of 1.3) |
| `f_T` [Hz] | 2.07398e+10 | 7.46504e+10 | **×3.5994** |
| `f_max` [Hz] | 1.41692e+10 | 5.66868e+10 | **×4.0007** |

Both branches hold. The ratio falls short, so the **fourth** lever is
exhausted and §7.8.1a's verdict on `f_max/f_T` is final rather than
provisional. And the same counterfactual is worth **a factor of four in
`f_max` itself**, so the RF case is also a contact-engineering problem with a
quantified target. The item could not see the second because **`f_max/f_T`
divides out exactly the quantity the contact dominates**: the contact limits
`I_d`, `I_d` sets `g_m`, `g_m` sets `f_T`, and the same contact sets `R_s` in
the `f_max` denominator. Removing it moves numerator and denominator together.

Residual after the perfect contact: `g_ds` must still come down **×2.9955**
(R1, with R1b recomputing it by bisection sharing no algebra with the closed
form), and 10-05 bounded what `v_sat` can supply at ×1.2491 over the whole
physically available phonon range. Term B is 1.11 % of the denominator, so it
is still term A that has to move.

### 2. The coherence trap, which inverts the sign of the answer

`R_c` enters this model in **two** places: the self-consistent solve for the
intrinsic channel drop, and `R_s = R_c/2` in the `f_max` denominator, reached
through another module that reads `gfet.Rc_total` as a **global**. Passing
`Rc_total=0.0` through `**kw` — the obvious route, and the one
`fT_fmax_saturated` would propagate — zeroes only the first:

| counterfactual | `f_max/f_T` | change |
|---|---|---|
| `R_c = 0` current path only | 0.552349 | **−19.15 %** |
| `R_c = 0` coherent | 0.759363 | **+11.15 %** |

**The half-removed perfect contact answers the item with the opposite sign,
and the half-removed answer is the plausible one** — a 19 % degradation is
exactly what "the contacts are binding" predicts if read without care, and
nothing about the number looks wrong. C3 asserts the sign inversion, so the
trap sits in the transcript rather than in a comment. This is the 10-03
two-names fault (a frozen cross-module default) in a new location: a
counterfactual that arrives in one of the two places the quantity enters.

### 3. The two models disagree on the SIGN of the contact lever

The same coherent counterfactual in the superseded resistor model:
**0.642901 → 0.581707, −9.52 %**, against the saturated model's +11.15 %
(C4). Not a magnitude discrepancy — a sign disagreement, asserted and left
unreconciled.

It is 10-04's resistance-ratio identity acquiring a falsifiable consequence.
`f_max/f_T ≈ (1/2)√(R_total/(R_g+R_s))` carries `R_total` in the
**numerator**, so a model obeying it is *rewarded* for a worse contact.
10-04 read that off the algebra; today it has a sign test. **A figure of merit
that improves when the device gets worse is not measuring the device.**

And Section 6b found the saturated model is not exempt, only better-behaved
over a range: **`f_max/f_T` is NON-MONOTONIC in `R_c`** (Q3), bottoming at
0.680047 at 470 Ω·µm and rising to 0.698615 at 4000 Ω·µm, while `f_max`
falls **37×** across those same two rows. The ratio and the figure of merit it
summarises move in **opposite directions** over part of this design space.

### 4. What 10-05 predicted, measured from the other side

10-05 explained its ×1.1361 against a required ×4.8483 by dilution:
`µS/L = 0.3036` "alone would divide the CHANNEL conductance by 1.699", and
the contacts absorb the rest. That is an untested quantitative prediction —
at `R_c = 0` the factor must be 1.699. **Measured: 1.70661 (100.45 %)** at
the item's own bias, 1.69621 (99.84 %) at the other, from a finite difference
of a solved terminal current rather than from the channel integral that
produced the 1.699 (D1). The prediction was recorded before the measurement
existed and the two computations share no machinery, which is a stronger form
of evidence than the agreement-between-two-artefacts 10-02 warned about.

### 5. The reachability sweep, and the objection it closes

Section 6b sweeps per-contact resistances that have been **measured**:

| `R_c` [Ω·µm] | `f_max` [GHz] | `f_max/f_T` | share of the `R_c = 0` gain |
|---|---|---|---|
| 0 | 56.687 | 0.759363 | 100.0 % |
| 65 | 39.526 | 0.722202 | 59.6 % |
| 165 | 24.308 | 0.695514 | 23.8 % |
| **300 (this model)** | **14.169** | **0.683186** | 0.0 % |
| 470 | 8.268 | 0.680047 | −13.9 % |
| 4000 | 0.222 | 0.698615 | −32.8 % |

Q1: the lowest reported contact resistance captures **59.6 %** of the perfect
contact's `f_max` gain, so most of the ×4 is reachable with contacts that
exist. Q2 is the load-bearing row: **this model's 300 Ω·µm is 1.8× the
165 Ω·µm Feijoo *et al.* report for the very devices whose 1.3–1.4 band
Chapter 4 is measured against.** The whole `f_max` thread could therefore have
been nothing but a worse contact. Given Feijoo's own contact, this model
reaches 0.695514, not 1.3. That objection is now closed by measurement.

### 6. The literature cuts against the premise of the 10-04/10-05 thread

**Feijoo, Pasadas, Bonmann *et al.*, "Does carrier velocity saturation help to
enhance f_max in graphene field-effect transistors?", *Nanoscale Advances* 2
(2020), DOI 10.1039/c9na00733d** — a paper that asks today's question in its
title — answers it in the negative: "the largest `f_max` are located at biases
close to the onset of bipolar conduction and **far from the saturated velocity
regime**", at ~45 % of `v_sat`, with the **diffusion** contribution to the
current comparable to drift, concluding that "our results do not support that
operating in the regime of velocity saturation results in the highest
`f_max`". They also measure **self-heating alone degrading `f_max` from 65 to
40 GHz**.

Three consequences, all annotated in place rather than absorbed:

1. **10-05's criterion A is narrowed, not retracted.** Its arithmetic stands
   (×4.8483 required, ×1.1361 delivered). What it licenses does not: it was
   constructed from this model's own algebra and is a statement about this
   model, not about real GFETs.
2. **This model has no diffusion current.** Eq. (4) is drift-only. If
   diffusion is ~40 % of the current near the drain at the peak-`f_max` bias
   in a measured device, a drift-only model cannot be asked about that bias at
   all, and every `V_ds` ladder of 10-05 and today is a drift-only slice.
   **This is the RF thread's new binding question, ahead of saturation.**
3. **This model has no self-heating**, a mechanism worth 65 → 40 GHz. Today's
   ×4.0007 is a ratio between two runs of one model and is largely insulated;
   the absolute 56.687 GHz is not.

**Fetch failures, recorded honestly.** `pubs.aip.org` **403'd** on *AIP
Advances* 11, 045220 (not retried). PubMed/PMC **not attempted** per the
standing reCAPTCHA note. `arxiv.org/pdf/1110.0978` fetched but is **unusable
for this question** — it states it includes "not a model for contact
resistance, but only a contact resistance parameter" — recorded so no future
session re-fetches it for the same purpose. The 65 Ω·µm row of Section 6b is
a **secondary citation** (Liu *et al.* 2019 via Khosravi Rad *et al.*,
*Sci. Rep.* 14, 9190 (2024)) and is labelled as such. Publishers are now the
fourth distinct 403 source for this repository, and **every literature number
here that came from a 403'd source came from a secondary citation** — the
10-01 fault in its most ordinary form.

### 7. The mutation harness, and three holes that were one fault

Final **7 of 7**, control A green, control B survived. First run **3 of 7**.
A counterfactual is a claim about what was switched off, so the mutants attack
**arrival and coherence** rather than arithmetic. The three survivors were one
fault in three places — a quantity or a claim nothing in the battery read:

- **M5** (one-sided `g_m` stencil) survived because `f_max/f_T` is nearly
  **blind** to `f_T` here: term B is 1.11 % of the denominator, so a 1e-5
  shift in `f_T` moves the ratio by ~5e-8. *That is the same fact as §1's.*
  Pinning `f_T` against the committed transcript would **not** have worked —
  the committed `f_T` is an `np.gradient` on a 0.015 V gate grid while this
  module differences at 1e-4, a different **estimator**, differing 1.4e-5 on
  the unmutated module. That gap is the committed `f_T`'s `O(h²)` gate-grid
  discretisation and **had never been quantified**. Repairs: **X2b**, a
  Richardson ratio on stencil order (3.9999 against 4 for central, 2 for
  one-sided), and **X2c**, which reports the estimator gap instead of
  asserting it away.
- **M6** (term B dropped from the residual formula) survived because R1
  asserted only finiteness. Repair: **R1b**, bisection on `f_max/f_T`.
- **M3** (`RC_ZERO` → 1e-12) survived **two** runs because X1a was written
  against a **literal** `0.0` while the counterfactual used `RC_ZERO` — an
  exactness check pointed at a value the module under test did not use. The
  09-30 arrival fault **inside an exactness check**, a new location for it.

Each of M3–M6 is now killed by **exactly one** check, and in each case it is
the only check that reads the quantity the mutant changed. Battery 17/17 →
20/20 → 23/23. **Fifth consecutive harness to improve the suite it was
pointed at rather than bless it.**

### 8. Methodological note, continuing the series

09-28: prose is a detector. 09-29: a mutation that does not arrive is
indistinguishable from a system that does not respond. 09-30: a control has to
sit where the failure enters. 10-01: and name what it compares against in a
way that cannot drift. 10-02: when it agrees, that is a fact about two
artefacts, not the world. 10-03: a PASS/FAIL at zero is silent about
magnitude. 10-04: and a quantity can be computed correctly under the wrong
name. 10-05: and a criterion with a number in it can be missed in a way that
locates the real constraint. **10-06: AND A CRITERION CAN BE DECIDED
CORRECTLY ON THE QUANTITY IT NAMES AND STILL MISS THE FINDING, IF THAT
QUANTITY IS A RATIO AND THE MECHANISM SCALES BOTH ITS ARGUMENTS.**

The 10-05 item did everything right. It was answered on exactly the quantity
it named, in advance, with both branches stated. And the ×4.0007 — the only
part of the result an industrial reader would act on — appears in neither
branch. **The practical rule is checkable, not advisory:** for every criterion
written on a ratio, ask what the mechanism does to numerator and denominator
*separately*, and if it does the same thing to both, the criterion is the
wrong instrument whatever number it carries. Q3 is the proof by construction —
a region of this very design space where satisfying a ratio criterion means
losing on the figure of merit.

**Not yet covered (candidates for future runs):**

- **A DIFFUSION TERM IN EQ. (4)** — created today and **the new top item**,
  and it displaces velocity saturation as the RF thread's binding question.
  Feijoo *et al.* 2020 measure diffusion as comparable to drift at the
  peak-`f_max` bias, in devices at ~45 % of `v_sat`. Eq. (4) is drift-only, so
  **every `V_ds` ladder in this repository's RF thread is a drift-only slice**
  and the peak-`f_max` bias is a bias this model cannot be asked about. This
  is the one item that would let the RF thread say something about a device
  rather than about a model. Chapter 7 §7.9 item 7.
- **CENSUS EVERY CRITERION AND CHECK IN THIS REPOSITORY THAT IS A RATIO** —
  created today, and it is the 10-03 magnitude item's sibling rather than a
  new thread. 10-03 asked which checks are ratios of two quantities the same
  defect would scale; today shows the same fault one level up, in the
  *criteria* rather than the checks, and gives it a worked example with a
  measured non-monotonicity behind it. The mechanisable form: for every
  criterion written on a ratio, record what the mechanism under test does to
  numerator and denominator separately.
- **MECHANISE the magnitude item** — created 10-03, worked examples on 10-04,
  10-05 and **twice more today** (M5 and M6 were both pure-magnitude holes).
  **Five modules, six defects, same result.** The mechanisable form is a
  repository-wide census of every MUST_CHANGE recording its response size.
  10-05 said "today is the last session at which one more worked example is a
  defensible answer"; today produced two more and is therefore **overdue, not
  pending**.
- **A self-heating term** — created today. Feijoo *et al.* 2020 measure
  65 → 40 GHz in `f_max` from self-heating alone, which bounds every absolute
  frequency in Chapter 4. Ratios between model runs are much less exposed, so
  this is a limitation to state rather than a lever. Chapter 7 §7.9 item 8.
- **Rewire `transfer_characteristic()` to the saturated form, or decide not
  to** — created 10-05, priced there (five modules read it, §5's
  profile-domain correction moves them a further signed 0.25 %), **and today
  makes the decision harder rather than easier**: the saturated form is now
  known to be drift-only against a literature result that says diffusion
  matters at the bias of interest, so rewiring would propagate a model whose
  known defect is named. Recommended: **defer until the diffusion item is
  settled**, and that is a change of recommendation from 10-05's "yes".
- **Why does the peak-f_T bias MOVE BRANCH between the two models?** — created
  10-05, −1.1128 V (hole) to +2.8571 V (electron). Today adds a datum: at
  `R_c = 0` there is **no interior peak at all** in either model, so the peak
  is a contact artefact in both and the branch question may be a question
  about where `R_c` caps `I_d` rather than about electron–hole asymmetry.
  **Narrowed, not closed.**
- **A second anchor for `Delta_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, **untouched for sixteen
  consecutive sessions**. Today was a device-physics session and again did not
  touch it.
- **Whether Chapter 4's `n_puddle` reproduces Section 3.3's 6.45 kΩ/sq floor**
  — created 09-29, **untouched for eight sessions**.
- **A description of Ti, Ni and Pd that does not go through work function** —
  three independent failures on record; Chapter 6's central structural
  weakness, and §7.8.1c keeps it a limitation on the RF case too.
- **A β = 2 implementation** (10-05, bounded at 24.6–26.0 % and moving the
  answer the wrong way, so completeness not a lever); **which other quantities
  here are computed correctly under a name that claims more than they
  measure?** (10-04); **which checks are RATIOS of two quantities the same
  defect would scale?** (10-03, and see the census item above); **audit every
  remaining check for the 10-02 §1 fault** (10-02); **do the four existing
  audits have references that are not content-pinned?** (10-01); **is
  `D_CHEM_PHYS` the only duplicated literature VALUE?** (10-01); **is
  `plot_fermi_surface` the only uncommitted figure?** (10-02); **the Section
  3.6 Pauli edge is absent from Chapter 6's model** (09-29); **where does the
  6.5430 % internal collection efficiency come from?** (09-30); **whether the
  three identity re-assignments should be real mutations** (09-30); **Chapters
  2 and 3 are drafted but §4.8.2's form of error is not audited for
  elsewhere** (10-01); whether Chapter 4 or 5 contains a RANKING that is a
  step artefact (09-23/09-25); a finite-temperature optical conductivity
  (09-29); angular trigonal warping (09-28); finite-temperature `n(E_F, T)`
  (09-28); whether other `== 0.0` checks are round-trip tautologies (09-28);
  whether RESULT 2's other convictions are read from the saturated branch
  (09-28); a probe reporting both directions of lambda (09-28); whether `t'`
  can be EXCLUDED quantitatively as the source of Chapter 4's electron–hole
  asymmetry (09-28); whether any other near-cancellation shows the `c2`/`c4`
  sign flip (09-27); migrating every remaining validation to
  measured-value-beside-derived-bound form (09-27); whether §7.3's
  differential-interconnect prediction holds (09-27); which of Chapters 4 and
  5's design rules could be restated as parities or bounds (09-23); whether
  the parity survives a photo-thermoelectric term (09-23); re-check whether
  other "for every …" claims rest on small samples (09-20); whether Chapter
  4's contact-resistance results should be re-run at the 5.4 eV crossover
  (09-18); reconciling Mueller *et al.*'s 0.12 eV step with the 0.25–1.07 eV
  offsets `METAL_WORK_FUNCTIONS` assumes (09-18); a photo-thermoelectric term
  (Kasirga review); Shimomura *et al.*'s comb-electrode design; integrating
  6.5's plasmonic near-field picture with the spatially-resolved
  contact-doping machinery; isolating the root cause of the §4.7 negative
  residual (08-31); Ti and Cr per-metal `Rc` recalibration (ResearchGate
  rate-limiting); a second independent edge-contact dataset (Lee *et al.*
  2022, Wiley 403'd); **the primary source for the 65 Ω·µm contact resistance**
  (Liu *et al.* 2019, cited only secondarily today)

**Automation health.** Device reachable and folder connected at the 04:30
firing (05:16 UTC); Step 0 found neither repository carrying a 2026-10-06
entry and no commits since midnight, so this was a full session.

**The 10-05 network problem is GONE, and the measurement is recorded so the
pattern can be read.** The device VM cloned both repositories in a single
`device_bash` call and `pip install scipy` completed in **7.4 s**, against
10-05's sustained 13.9–18.8 kB/s where no clone of any depth fit the 180 s
ceiling. So the problem is intermittent across days (10-01, 10-02, 10-05 bad;
10-03, 10-04, 10-06 fine) and the documented bundle route remains the thing to
keep rather than a workaround to retire.

**A NEW device-VM constraint, and it is not the network.** The device VM runs
**Python 3.10.12**, and this repository's committed modules use PEP 701
nested-quote f-strings, which are a **3.12+** feature:
`graphene_velocity_saturation_model.py` raises `SyntaxError` there. So the
graphene battery cannot be run on the device VM at all, whatever the network
does, and the cloud container (Python 3.13.16, numpy 2.5.3, scipy 1.18.1) is
now the only place this repository executes. That is a **standing** constraint,
unlike the network, and it means the container-plus-bundle route is the only
route for this repository rather than the fast route. `scipy` was present in
the container and no install was needed there.

The `GIT_ASKPASS` recipe ran from the device VM; the token was never written
into `.git/config`, a remote URL, any repository file or the connected folder,
and every transient copy is shredded at the end of the session.
`user.name`/`user.email` were again absent in the fresh clones and set per
09-24. **Push verified against the GitHub API rather than git's own output**,
per 09-26, and `git fetch origin main` was run immediately before the push per
the 10-03 race item. **Live web search was available and used**, and two of the
three literature fetches succeeded. **One deviation from 10-05's practice,
stated:** push-as-you-go was **not** followed commit-by-commit this session.
With the container as the only execution environment and the device VM as the
only push point, each push costs a bundle transfer, so the five graphene
commits were pushed as one batch at the end of the repository's work rather
than in two. The incremental bundle was **197 kB**, against 10-05's 11.9 MB
full-history bundle, because `git bundle create <old_origin_main>..main`
carries only the new objects and the device VM already held a full clone.
**That is a better recipe than 10-05's and is the one to reuse.**

**Commits this run:** 5 (the counterfactual with its transcript and figure;
the mutation harness with its report and the four checks it forced; Section 6b
with the reachability sweep and the non-monotonicity; the session note;
Chapter 4 §4.6.4 with Chapter 7 §7.8.1c). This AUTOMATION_LOG.md entry makes
**6**.

---

## 2026-10-07 — The diffusion term was built exactly as specified, and its premise did not survive the build

**Status.** The 2026-10-06 top item — *a diffusion term in Eq. (4)* — is
**closed, and it closed by moving its own premise.** New module
`graphene_diffusion_current_model.py` (**30/30**, nine tolerance-free
checks), new mutation harness (**7 of 7**, two controls surviving), new
figure, new note, Chapter 4 **§4.6.5** and Chapter 7 **§7.8.1d** drafted,
§7.9 item 7 **rewritten rather than ticked off**, three new §7.9 items. One
claim this repository had committed is **withdrawn**. The device VM's
standing "cannot run graphene at all" constraint is **retired**, and
replaced by a narrower and measured one.

### 1. The answer, in the order the item asked

| the item said | the measurement says |
|---|---|
| Eq. (4) is drift-only | **it is undetermined** — §3 |
| diffusion is comparable to drift | **±0.50 % of `I_d` here** — §2 |
| the peak-`f_max` bias cannot be asked about | **true, for a different reason** — §4 |

### 2. The derivation, and the thing it turned out to be

Zebrev's generalised Einstein relation (`arXiv:1102.2348`, Eq. 20) with
`eps_D = E_F/2` for graphene's linear dispersion gives `D = mu*E_F/(2e)`, and
with `V_F = E_F/e = A_F*sqrt(n)` the diffusion term collapses **exactly** to
`mu*n*dV_F/dx`. So the drift-diffusion current is `mu*n*d(V - V_F)/dx` —
transport down the gradient of the quasi-Fermi potential — and the `beta = 1`
separability that makes Eq. (4) closed-form **survives**. `lambda = 0`
recovers Eq. (4) **bitwise**, 0.0 ULP at twelve (`V_g`, `V_ds`) points (X1):
the overlapping-limit rule of 2026-10-03 in its strongest available form.

**`int n dV_F` has an exact antiderivative**, `(A_F/3)[n^{3/2}]`, so the
diffusion contribution to the channel integral is a **boundary term** — it
reads only the source and drain densities and not the profile between them.
That is worth more than the number it produces. 2026-10-05's criterion C
found this model's profile-domain handling carries a one-signed error **674×**
the Jensen gap it was looking for; the diffusion term **cannot inherit it**.
**X6a asserts this structurally rather than arguing it**: a deliberately
lopsided interior grid with the same endpoints leaves the term **bitwise**
unchanged, and **X6c** is the contrast — the quadrature route to the same
number moves 0.4 % on the same regridding. A structural claim with a
structural check is a better artefact than a structural claim with a
numerical one, and this is the first time in this series the distinction has
been available.

The size, at the two fixed biases 10-06 reported at, `V_ds = 0.1 V`:

| bias | `I_d` (λ=0) | `I_d` (λ=1) | share | `Qd/Q` |
|---|---|---|---|---|
| saturated-peak, `+2.857143 V` | 9.007791e-05 | 9.052996e-05 | **+0.5018 %** | +1.1621 % |
| resistor-peak, `−1.112782 V` | 8.910103e-05 | 8.865209e-05 | **−0.5038 %** | −1.1546 % |

**The correction is SIGNED and its sign flips between branches** (G1b). It
follows `d|n|/dV`, and `carrier_density()` returns a magnitude, so on the hole
branch `|n|` rises toward the drain and the term subtracts. A model that has
discarded the carrier's sign cannot say whether that is right for holes;
**new §7.9 item 9**, and it bears on the still-open "why does the peak-`f_T`
bias move branch" question. The share of `I_d` is half of `Qd/Q` because
`S_d > 0` — Eq. (12) opposes Eq. (11) — which the mutation harness found, not
the module (§5).

**Zebrev's closed form as an independent check.** His Eq. 62 gives
`kappa = C_ox/(C_Q + C_it)`: capacitor algebra against a boundary term, no
shared machinery. Interior mean ratio **0.6927** (0.011699 against 0.016677 at
the source end). Agreement to 31 % for two routes this different is
agreement, and the residual — this model's series factor and puddle floor —
is reported rather than tuned away.

### 3. The premise, and the claim that is withdrawn

Pasadas and Jiménez, *IEEE TED* 2016 (`arXiv:1605.08235`) — **the same group
as the Feijoo *et al.* paper that created the item** — write the GFET channel
current as `I = -W*Q_tot*v` with `v = mu*F`, `F = -dV/dx`, and state
explicitly that *"V(x) is the quasi-Fermi level along the graphene channel"*.
A drift expression whose driving potential is the quasi-Fermi level **is
already the complete drift-diffusion current.** Under that reading
`lambda = 1` is a **double count**, not missing physics.

`gfet.carrier_density(V_g, V_ch)` contains evidence for both readings and
settles neither: `V_ch` enters an electrostatic charge relation
`C_ox*(V_g - V_ch - V_dirac)`, and is then corrected by a quantum-capacitance
series factor — which is exactly what one applies when `V_ch` is the
quasi-Fermi level and the graphene drop `E_F/e` is separate. **The two
readings differ by exactly `lambda*V_F`**, the term itself. So ±0.50 % is
simultaneously the diffusion current under one reading and **the size of an
ambiguity that was already in every committed `I_d` in this repository.**

**Withdrawn:** the claim, made in §4.6.4 and §7.9 item 7 and repeated in the
10-06 log, that this repository's RF ladders are *drift-only slices*. They are
slices of a model whose driving potential was never named, bounded at ±0.50 %.
The numbers stand; the description of what they are did not. §7.9 item 7 is
**rewritten** to the actionable form — *name the potential* — rather than
ticked off, because building the term did not answer it and cannot.

`lambda = 1` is deliberately **not** wired into any other module. That is the
10-05 discipline (a rewiring is a priced decision, not a side effect) applied
to a term whose correctness is undecided rather than to one whose cost is
known.

### 4. A prediction written down first, refuted, with its mechanism measured

Zebrev's `kappa = C_ox/C_Q` **diverges** at charge neutrality and Feijoo
*et al.* put the peak-`f_max` bias there. Both predict the share **peaks at
the Dirac point**, and that was written into the module before the sweep ran.

Measured: it **collapses** there. `−0.006925 %` at the nearest grid point to
the Dirac point against `0.5053 %` at the maximum — a factor **73** — with the
maximum two volts away at `V_g = −1.272 V`.

G4b measures the mechanism instead of arguing it:
`n = sqrt(n_eff^2 + n_puddle^2)` makes `d|n|/dV = (n_eff/n)*dn_eff/dV`
**exactly 0.0** at `n_eff = 0`, while `C_Q` stays finite at `n_puddle` —
`0.000000e+00` against `−2.290589e+15 m^-2/V` unfloored. **The puddle floor
flattens precisely the region where the capacitor ratio diverges.** The
diffusion-dominated regime is not small in this model, it is **absent**, and
that is a fact about the regularisation rather than about graphene. 10-06 was
right that the peak-`f_max` bias cannot be asked about; the reason is the
floor, and **no amount of care about the transport term recovers it.** New
§7.9 item 11 asks which *other* near-Dirac results rest on a derivative the
floor has set to zero.

### 5. The harness, and the survivor that was not a hole

**7 of 7** after a first run of 6 of 7; control A green, controls B and C
surviving. Battery **26/26 → 30/30**. The mutants attack the **structure** of
the claim: **M6** swaps the closed form for the quadrature route to the same
number — a <0.1 % numerical change that destroys the only claim in the module
title — and X6a is the only check that reads it.

Two checks forced:

- **X9a/b/c.** λ enters in **two** places and M2 arrives at only the first. It
  was killed only by **G5c**, a check written for something else; *an
  incidental detector is not a detector.* X9 gave Eq. (12) its own check with
  a magnitude, and in doing so **found the `S_d > 0` fact that explains §2's
  otherwise unexplained factor of two.** A repair that explains an existing
  observation is a stronger outcome than one that only closes a hole.
- **G5d.** The stencil is pinned against the **committed** 0.683186
  (reproduced to 1.010e-07), not only against `vsat` on the same stencil. G5a
  compares λ=0 against `vsat` *on the same stencil*, so a stencil change moves
  both together and **G5a is blind to it by construction** — 10-06's M5 fault
  exactly.

**AND THE PART WORTH THE MOST.** The first run's survivor was the stencil
*spacing*, `6/399 → 6/400`. The reflex — five harnesses deep — is to file a
survivor as a hole and add a check. Measured before filing: the change moves
`f_max/f_T` by **−6.560e-11** relative and `f_T` by 5.275e-08, against a value
committed to **six** figures. A central difference is `O(h^2)` accurate
precisely *so that* its answer does not depend on `h`, so the mutant is
**semantically inert at the precision the claim is made at**, and a check
sensitive to it would be a check on an artefact. It is **control C** now —
it must survive — carrying its measured size in its own description, and M7
was replaced by a mutation that is not inert (the stencil *centre* becoming
its one-sided *edge*, killed by G5a and G5d).

That is **the verification repository's 2026-10-01 rule arriving here for the
first time**: *a constant that can be mutated with no observable effect is not
a suite weakness.* The two repositories' methodological series have run
separately since 08-21, and this is the first rule to cross between them. It
crossed because the same reflex produced the same wrong filing in both, which
is the best available evidence that these rules are about verification rather
than about either subject matter.

Sixth consecutive harness to improve the suite it was pointed at rather than
bless it.

### 6. A second finding, annotated in place

`gfet.quantum_capacitance()` evaluates `C_Q` at the **gate overdrive** in the
slot the dispersion wants `E_F/e`. At the RF bias: **2.057 V against
0.0977 V, a ratio of 21.05.** So this repository's `C_q` overstates `C_Q` and
would understate `kappa` by that factor. **No committed number is withdrawn**:
`C_q` is used downstream only through `C_q/(C_q+C_ox)`, where `C_q >> C_ox`
makes the factor ≈ 1 either way. Chapter 4 §4.6.5 carries the annotation
beside §4.3's numbers rather than recomputing them, per the 10-02 rule.

### 7. The RF consequence

| quantity | λ = 0 | λ = 1 | change |
|---|---|---|---|
| `f_T` | 2.073963e+10 | 2.085305e+10 | +0.5469 % |
| `f_max` | 1.416902e+10 | 1.421070e+10 | +0.2941 % |
| `f_max/f_T` | 0.683186 | 0.681468 | **−0.2514 %** |
| `g_ds` | 3.374491e-02 | 3.391519e-02 | +0.5046 % |

`f_max` rises and the ratio **falls**. 10-06's finding in a **third**
mechanism, and this time the two move in *opposite directions* rather than
merely by different amounts. §7.8.1a's verdict survives: 52.42 % of 1.3 either
way, so the term that was supposed to threaten the structural conclusion does
not.

### 8. Methodological note, continuing the series

09-28: prose is a detector. 09-29: a mutation that does not arrive is
indistinguishable from a system that does not respond. 09-30: a control has to
sit where the failure enters. 10-01: and name what it compares against in a
way that cannot drift. 10-02: when it agrees, that is a fact about two
artefacts. 10-03: a PASS/FAIL at zero is silent about magnitude. 10-04: and a
quantity can be computed correctly under the wrong name. 10-05: and a
criterion with a number in it can be missed in a way that locates the real
constraint. 10-06: and a criterion can be decided correctly on the quantity it
names and still miss the finding, if that quantity is a ratio and the
mechanism scales both its arguments. **10-07: AND AN ITEM CAN BE WELL-POSED,
CORRECTLY MOTIVATED FROM THE LITERATURE AND BUILT EXACTLY AS WRITTEN, WHILE
ITS PREMISE — THAT THE MODEL LACKS THE TERM — IS NOT A FACT ABOUT THE MODEL
BUT AN ARTEFACT OF A VARIABLE WHOSE MEANING WAS NEVER FIXED.**

The checkable rule: **before adding a term to a transport equation, ask what
the potential in it is.** A drift expression in the quasi-Fermi potential and
a drift-diffusion expression in the electrostatic potential are the *same
equation*, and no amount of care about the term distinguishes them. The
mechanisable form is new §7.9 item 10 — a census of every potential-like
variable in this repository recording which potential it is. `V_ch` was the
first one asked and the answer was "both".

### 9. A near-miss in this session's own method, recorded because it is the genre

An attempt to verify that `graphene_velocity_saturation_model.py` reproduces
its committed transcript on the **device VM** was cut off by the device
shell's 180 s ceiling at 175 s (exit 124). `diff -q` against the committed
transcript then reported **no differences** — because the module writes its
transcript only at the end, so the file still held the committed copy the
timed-out run had never overwritten. **A diff against a file the run did not
write is not a reproduction check**, and it was one command away from being
logged as a PASS. This is 09-29's shape (a mutation that does not arrive looks
like a system that does not respond) in the **verification of an environment**
rather than of a model, and the defence is the same one: *the control has to
sit where the failure enters* — delete the artefact first, then check that it
was recreated.

### 10. Environment: one standing constraint retired, a narrower one measured

10-06 recorded as **standing** that this repository cannot execute on the
device VM at all: it runs Python 3.10.12 and
`graphene_velocity_saturation_model.py` raised `SyntaxError` because two table
headers containing an apostrophe were written as PEP 701 nested same-quote
f-strings (3.12+). The consequence recorded was two interpreters, two
execution environments and a four-file bridge crossing that delivered one
stale copy.

**Retired today.** The two headers are bound to names and interpolated;
`ast.parse` under 3.10 now accepts **every** `.py` file in the repository, and
the whole of today's work — module, full-resolution transcript, figure and a
nine-run mutation harness — ran **on the device VM**, in 14.5 s and 66 s
respectively. Nothing crossed the bridge, so 10-06's stale-copy failure mode
did not exist this session.

**Verified output-neutral rather than asserted**, which is the right control
for a string-formatting change: the patched module was run at full resolution
on **Python 3.13.16** in the cloud container and
`velocity_saturation_output.txt` came back **byte-identical** to the committed
transcript. Cross-version byte-identity on 3.10 is **not** claimed — see §9
for why the attempt does not support it.

**The replacement constraint, measured.** `graphene_velocity_saturation_model.py`
takes **2 m 11 s** in the cloud container and **more than 2 m 55 s** on the
device VM, against a 180 s `device_bash` ceiling, and background processes do
**not** survive a `device_bash` call (`nohup` and `setsid` were both tried and
both were killed). So that one module's full transcript cannot be regenerated
on the device VM in a single call. That is narrower than a `SyntaxError` — it
affects one module rather than the repository — and the right response is to
regenerate that transcript in the container when it is ever needed, not to
move the repository back there.

**Automation health.** Device reachable and folder connected at the **04:30**
firing (04:57 UTC); Step 0 found neither repository carrying a 2026-10-07
entry and no commits since midnight, so this was a full session. `pip install
scipy` took **3 s** (10-06: 7.4 s; 10-05: unusable), so the network was
healthy for the fourth session in five. **Live web search was available and
used**; two of three literature fetches succeeded, and `link.springer.com`
returned a **429 with an explicit do-not-retry** — a new refusal class, and
publishers remain the fourth distinct refusal source. PubMed/PMC not
attempted per the standing note, though one plausibly relevant PMC result was
visible in search and is left for a session with another route to it.

**Not yet covered (candidates for future runs):**

- **NAME THE POTENTIAL IN EQ. (4)** — created today and **the new top item**,
  and it replaces the diffusion item rather than following it. `V_ch` is
  undetermined between the electrostatic and quasi-Fermi readings, which
  differ by exactly the term built today, so every `I_d` in this repository
  carries a ±0.50 % ambiguity that no further term can remove. The actionable
  form: declare which potential `V_ch` is, consistently across
  `carrier_density()`, the series factor and Eq. (4), and re-derive. Chapter 7
  §7.9 item 7, rewritten.
- **A CENSUS OF EVERY POTENTIAL-LIKE VARIABLE IN THIS REPOSITORY** — created
  today, the mechanisable form of the item above, and §7.9 item 10. Nothing
  has asked this of Chapter 5's interconnect channel potential, Chapter 6's
  two-contact photodetector channel potential, or the local Dirac-point shift
  in §4.4. One ask produced one "both"; the base rate is unknown.
- **A SIGNED-CARRIER FET MODEL** — created today, §7.9 item 9. Today's
  correction flips sign between branches because it follows `d|n|/dV` and
  `carrier_density()` returns a magnitude, so **the hole-branch sign of every
  branch-asymmetric Chapter 4 result is unverified.** This repository already
  has a signed-carrier treatment for photodetectors and nothing equivalent for
  the FET. This is also the strongest new datum on the open
  peak-`f_T`-moves-branch question.
- **WHETHER `n_puddle` REMOVES OTHER REGIMES, NOT ONLY THE
  DIFFUSION-DOMINATED ONE** — created today, §7.9 item 11. The floor sets a
  *derivative* exactly to zero, which deleted the near-Dirac regime two
  independent literature sources point at. Which other near-Dirac results rest
  on that derivative is not known. Sibling of the §4.9 item on whether
  `n_puddle` reproduces §3.3's 6.45 kΩ/sq floor — created 09-29, **untouched
  for nine sessions**, and today is the first session to give it a reason that
  is not bookkeeping.
- **MECHANISE the magnitude item** — created 10-03; worked examples on 10-04,
  10-05, twice on 10-06 and **again today** (M4 was a pure-magnitude hole
  killed by X5 alone, and G2's factor-of-3 band was blind to a 50 % change).
  **Six modules, seven defects, same result.** 10-05 said one more worked
  example was the last defensible answer; 10-06 called it overdue; today makes
  it **two sessions overdue**. The mechanisable form is a repository-wide
  census of every MUST_CHANGE recording its response size.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 10-06,
  **widened today**: `f_max` rose while `f_max/f_T` fell, so the ratio and the
  figure of merit moved in *opposite* directions for the third distinct
  mechanism. The 10-06 instance was a magnitude mismatch; this one is a sign
  mismatch, which is strictly worse and strictly easier to detect.
- **A self-heating term** — created 10-06, unchanged. Feijoo *et al.* measure
  65 → 40 GHz from self-heating alone, which bounds every absolute frequency
  in Chapter 4 while leaving ratios between model runs much less exposed.
  **Demoted below the naming item today**: a limitation to state rather than a
  lever.
- **Rewire `transfer_characteristic()` to the saturated form, or decide not
  to** — created 10-05, deferred 10-06 pending the diffusion question.
  **Still deferred, and for a firmer reason**: the question the deferral waited
  on is now known not to be answerable by adding a term. It waits on the
  naming item instead.
- **Why does the peak-`f_T` bias MOVE BRANCH between the two models?** —
  created 10-05, narrowed 10-06. Today adds a second datum: the diffusion
  correction itself flips sign between the branches, so there are now two
  independent branch asymmetries and one of them is known to be an artefact of
  carrying `n` as a magnitude. **Narrowed again, not closed**, and now
  coupled to item 9.
- **A second anchor for `Delta_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, **untouched for seventeen
  consecutive sessions**. Today was a device-physics session and again did not
  touch it.
- **Whether Chapter 4's `n_puddle` reproduces Section 3.3's 6.45 kΩ/sq floor**
  — created 09-29, **untouched for nine sessions**, and now the sibling of a
  live item rather than an isolated one.

**Commits this run:** 5 (the portability fix retiring the device-VM
constraint; the diffusion module with its transcript and figure; the mutation
harness with its report and the two checks it forced; the session note;
Chapter 4 §4.6.5 with Chapter 7 §7.8.1d and the §7.9 rewrite). This
AUTOMATION_LOG.md entry makes **6**.

---

## 2026-10-08 — The potential is named, and naming it exposed a question 1.4 times larger

Two items closed, both created or rewritten on 2026-10-07: §7.9 item 10 (the
census) and §7.9 item 7 (the declaration). Doing the census first was the right
order, and it is the census that produced the surprise.

### 1. The census: one defect in 86 names, and three predictions all wrong

`graphene_potential_census_audit.py` collects every potential-like name in the
repository — **86 across 25 modules** — and requires each to carry an explicit
adjudication in a registry inside the file. A name present in the repository
and absent from the registry is reported as a fault in its own right, so the
census **cannot go stale silently**. It proved that on the very next file
added: committing `graphene_potential_declaration.py` made the completeness
check FAIL until that module's nine names were adjudicated.

The name pattern is deliberately too wide — it catches `mu_acoustic`, a
mobility — because a pattern narrow enough to have no false positives has
already decided which names are potentials, which is the thing being audited.
`MOBILITY_NOT_A_POTENTIAL` is one of the readings, and 12 names carry it.

**The 10-07 entry named three sibling defects it expected. All three were
checked and NONE is one, for three different reasons:**

| 10-07's candidate | Found |
|---|---|
| Chapter 5's interconnect channel potential | **does not exist** — `graphene_interconnect_model.py` has no potential-like name at all |
| Chapter 6's two-contact channel potential | imported and printed, **enters no computation** in its own module |
| §4.4's local Dirac-point shift | appears **only in a docstring** derivation the evaluated closed form never uses |

The first is kept as a check **expected to fail**, and failing on every run,
because a pass there would mean a potential variable had appeared in Chapter 5.
Measured base rate: **1 defect in 86 names.** Seven names remain
`UNADJUDICATED` — open work, not a pass — and none is channel-swept, so none
can carry the defect.

### 2. Thirteen false positives, and the gate they forced

The first run reported **thirteen** dual-role names, including `V_g`, `V_hi`,
`mu_t` and `V_SAT_ANCHOR`. Every one was a false positive and instructively so.
`V_g` really does appear in a charge relation and a current relation in the same
function, and that is **not** a defect: a terminal bias is the same physical
quantity in both, so it has no second reading.

What made `V_ch` ambiguous was never that it appeared twice. It is that `V_ch`
is **channel-swept** — swept from 0 to the drain bias *inside one current
evaluation* — which is the only situation in which the electrostatic and
quasi-Fermi readings differ. The verdict now requires both roles **and**
channel-swept, detected by following `profile = linspace(0, Vds, n)` through to
`for V_ch in profile`; the benign case is reported as `DUAL_ROLE_TERMINAL` and
ten names sit there. **Mutants M4a/M4b exist so the gate is not decoration:**
M4a strips the declaration and requires `V_ch` to fall back to `UNDETERMINED`
(proving the *declaration* cleared the defect, not a classifier change); M4b
strips the declaration *and* the drain-bias name and requires `V_ch` to leave
the defect class. Both pass.

Two first-version faults worth recording as this repository's own classes:
role markers fired on **prose** (a comment about `sheet_conductivity`, and each
function's own `def` line matching its own name — comments, docstrings, strings
and signatures are now blanked first); and three registry rows had **no code
behind them** because the scanner saw only `ast.Name` and missed every
`gfet.V_dirac`, i.e. cross-module use through an `Attribute`. The dead-row
check earned its keep on its first run.

### 3. The declaration, and the literature that does not make it

> **`V_ch` is the quasi-Fermi (electrochemical) potential of the channel
> carriers, in volts, measured from the source.**

Declared in `graphene_fet_model.py` where `V_ch` is defined, as a
machine-readable in-source marker, following `graphene_dead_name_sweep.py`'s
convention. Three grounds: Pasadas and Jiménez (*IEEE TED* 2016,
arXiv:1605.08235) define exactly that variable in the compact model Eq. (4)
descends from; classical long-channel theory uses the same variable and heads
the current equation it feeds "both drift and diffusion"; and
`carrier_density()` already applies the quantum-capacitance series factor,
which is the correction for exactly that reading.

**And a literature result worth having.** Lu *et al.*'s review of GFET compact
models (arXiv:1703.09759) was fetched specifically to see whether the
identification is standard. **It is not made there**: §2.5 calls the imref
splitting "a.k.a. channel voltage" but never connects it to the `V(x)` of the
compact-model sections; Eq. (12) carries a diffusion term that Eq. (15) writes
with the same expression as drift, unexplained; and the §5 compact model uses a
drift form while only *assuming* drift-diffusion. **The ambiguity 10-07 found
is inherited from the literature rather than invented here** — which does not
excuse it and does raise the value of declaring it.

### 4. 10-07's term is a double count, measured three ways

Under the declaration, `λ = 0` is already the complete drift-diffusion current,
so `λ = 1` adds a term already present. At `V_g` = 2.0 V, `V_ds` = 0.1 V:

| Level | λ=0 | λ=1 | difference |
|---|---|---|---|
| intrinsic | 1.432116e-04 A | 1.442091e-04 A | **+0.6965 %** |
| terminal | 8.190217e-05 A | 8.225263e-05 A | **+0.4279 %** |
| 10-07, at its own saturated-peak bias | — | — | +0.5018 % |

The terminal figure **must** be the smaller of the first two because the series
contact resistance is negative feedback; that inequality is checked, and its
sign is fixed by the circuit rather than by a fit. Nothing is rewritten:
§4.6.5's numbers stand, and `λ` is retained as the instrument that *measures*
the difference between readings.

### 5. The question underneath, and it is larger

The declaration promises that the charge relation is the one `V_ch` belongs to:

    V_g − V_dirac − V_ch = e·n/C_ox + E_F(n)/e        (4.28, exact, closed form)
    n = (C_ox·dV/e)·C_q/(C_q + C_ox)                  (4.29, what the repo computes)

Running Eq. (4) on (4.28) with nothing else changed moves `I_d` at the RF bias
by **−0.6023 %** — **1.4 times** the term just retired. **Two sessions went
into a +0.43 % transport term while a −0.60 % question about the same chapter's
charge relation sat underneath it, reachable by replacing one function.** That
is §7.7's unequal scrutiny in a third instance, and the cause is nameable: the
attention was captured by a literature paper's framing rather than by an audit
of this repository's own equations.

**A number this session declined to quote.** The raw (4.28)-vs-(4.29) charge
disagreement reaches **+68 % at `dV` = 0.01 V** — in a regime `n_puddle` =
5e15 /m² swamps, where no committed result lives. In what `carrier_density()`
actually returns the worst is **+2.22 %**, and in `I_d` it is **−0.60 %**. The
raw figure would have overstated the declaration's cost by a **factor of 31**.
10-06's rule, applied to this session's own result.

### 6. Validation against exactly known values

- (4.28)'s closed-form root satisfies (4.28) to **3.8e-16** relative — an
  algebraic identity, not a plausible range.
- `dV = 0` gives `n = 0` **exactly** in both forms before the puddle floor.
- `n(+dV) == n(−dV)` **bitwise** in both forms (both carry `n` as a magnitude).
- The large-drive limit is measured and monotone for both, not asserted.
- `graphene_fet_model.py`'s output is verified **bitwise unchanged** by this
  session's edit over a 41-point sweep. The edit is comments plus one docstring
  annotation; 10-02 established byte-identity is not correctness, and it is
  exactly the right instrument for a comment-only change.

### 7. A new defect found on the way, measured and ranked last on purpose

`quantum_capacitance()` returns `inf`/`NaN` above `dV` = **18.3493 V** (found
by bisection), because `log(2(1+cosh η))` overflows once `η = dV/(kT/e)`
exceeds ~710; `kT/e` = 25.85 mV here. The exact large-`η` limit is
`|η| + log 2`. The maximum drive swept anywhere in this thesis is ~2.7 V, so
**no committed number is affected** — it is a *latent* defect, is §7.9 item 13,
and is deliberately ranked last. Not patched today because that module owns
committed transcripts across three chapters.

### 8. Two census rows adjudicated, one against the census's own guess

The census left `V_F` `UNADJUDICATED` on the guess that it was a Fermi
*velocity*. In `graphene_diffusion_current_model.py` it is not —
`fermi_voltage(n)` returns `E_F/e` in **volts**, and it is the very quantity
that converts between the two readings. In
`graphene_contact_doping_nonlinear_model.py` the same name **is** a velocity,
1.0e6 m/s. One name, two quantities, two modules: 10-01's fault in a third
place. The third row stays `UNADJUDICATED` because its module was not read
today, which is cheaper than guessing again.

### 9. Methodological note, continuing the series

09-28: prose is a detector. 09-29: a mutation that does not arrive is
indistinguishable from a system that does not respond. 09-30: a control has to
sit where the failure enters. 10-01: and name what it compares against in a way
that cannot drift. 10-02: when it agrees, that is a fact about two artefacts.
10-03: a PASS/FAIL at zero is silent about magnitude. 10-04: a quantity can be
computed correctly under the wrong name. 10-05: a criterion with a number in it
can be missed in a way that locates the real constraint. 10-06: and decided
correctly on the quantity it names while missing the finding, if that quantity
is a ratio. 10-07: an item can be well-posed, correctly motivated and built
exactly as written while its premise is an artefact of an unnamed variable.
**10-08: AND CLOSING SUCH AN ITEM CAN SHOW THE AMBIGUITY WAS THE SMALLER HALF.
NAMING THE POTENTIAL SETTLED ±0.50 % AND EXPOSED −0.60 % THAT THE NAMING
QUESTION HAD BEEN STANDING IN FRONT OF.**

The checkable rule: **when a declaration resolves an ambiguity, measure what
the declaration itself costs in the quantity you quote.** A declaration
promises the rest of the model belongs to the declared variable, and that
promise is computable — here, one monkey-patched function and one re-run.
Mechanisable form: new §7.9 item 14. Items 2, 3 and 6 all turn on a named
quantity and none has been costed this way.

Second rule, from §1: **a predicted class is not a measured class.** Three
siblings were predicted on 10-07 from one instance; the census found zero, with
three different reasons. Predicting where a defect class lives is cheap; the
census that measures it cost one session and is now standing.

### 10. Environment and automation health

Device reachable and folder connected at the **04:30** firing (04:50 UTC);
Step 0 found neither repository carrying a 2026-10-08 entry and no commits
since midnight, so this was a full session. `pip install scipy` succeeded (the
device VM still lacks it; 10-07: 3 s), so the network was healthy for the fifth
session in six. The whole session ran **on the device VM** — the 10-06
`SyntaxError` constraint stayed retired, and nothing crossed the bridge.
**Live web search was available and used**: three fetches attempted, two
returned usable content (arXiv:1703.09759, the long-channel MOSFET notes) and
one returned **table-of-contents text only** — Pasadas's 2017 dissertation
(arXiv:1709.01324), whose §§3.4.3–3.4.4 are the sections that would have
answered the question directly. Recorded as a *partial* fetch rather than a
failure: the extraction succeeded and the content was not in it. PubMed/PMC not
attempted, per the standing note.

**Not yet covered (candidates for future runs):**

- **RE-DERIVE CHAPTER 4 ON THE EXACT CHARGE RELATION (4.28)** — created today,
  §7.9 item 12, and **the new top item**. Worth −0.6023 % in `I_d` at the RF
  bias, 1.4× the term two sessions went into. **Blocked by the signed-carrier
  item**: both charge forms carry `n` as a magnitude, so re-deriving first
  would rebuild the branch-sign flaw in a new equation.
- **A SIGNED-CARRIER FET MODEL** — created 10-07, §7.9 item 9, **promoted today
  from sibling to prerequisite** of the top item. This is the first session in
  which it blocks something rather than sitting beside it, and it remains
  untouched.
- **MEASURE WHAT EACH DECLARATION COSTS** — created today, §7.9 item 14, the
  mechanisable form of today's second half. One function patched, one re-run.
  Items 2, 3 and 6 all turn on a named quantity; none costed.
- **MECHANISE the magnitude item** — created 10-03; worked examples on 10-04,
  10-05, twice on 10-06, on 10-07 and **again today** (the +68 % raw charge
  disagreement was a pure-magnitude trap that a sign or PASS/FAIL check would
  have waved through, and the factor between the quotable and unquotable
  version was 31). **Seven modules, eight defects, same result. Three sessions
  overdue.** The mechanisable form is a repository-wide census of every
  MUST_CHANGE recording its response size — and today's census is the proof
  such a sweep is buildable in one session, which removes the last excuse.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 10-06, widened
  10-07. Untouched today. Today's census is a worked template for it: the
  registry-plus-completeness-check pattern transfers directly.
- **WHETHER `n_puddle` REMOVES OTHER REGIMES** — created 10-07, §7.9 item 11.
  **Strengthened today**: the floor is now known to dominate the near-Dirac
  charge of *both* (4.28) and (4.29) — 0.001 % disagreement at `dV` = 0.01 V
  once floored, against +68 % raw — so whatever it silenced, it silenced in the
  replacement model too. Sibling of the §4.9 item on whether `n_puddle`
  reproduces §3.3's 6.45 kΩ/sq floor, created 09-29 and **untouched for ten
  sessions**.
- **THE `cosh` OVERFLOW** — created today, §7.9 item 13, **ranked last on
  purpose**: latent, no committed number affected, and a numerics change in
  `graphene_fet_model.py` needs its own output-neutrality pass.
- **The seven `UNADJUDICATED` census rows** — created today. `V_SAT_ANCHOR`,
  `mu_resp`, two `phi` and one `V_F` among them. None is channel-swept so none
  can carry the 10-07 defect; what is undecided is whether they are potentials
  at all. Cheap, and the census will keep reporting them until they are decided.
- **A self-heating term** — created 10-06, demoted 10-07, **unchanged and now
  third in the RF queue** behind items 12 and 9. Still a limitation to state
  rather than a lever.
- **Rewire `transfer_characteristic()` to the saturated form, or decide not
  to** — created 10-05, deferred 10-06 and 10-07. **The deferral's reason is
  now gone**: it waited on the naming question, which is closed. It waits on
  item 12 instead, which is a stronger reason and a worse position.
- **Why does the peak-`f_T` bias MOVE BRANCH between the two models?** —
  created 10-05, narrowed 10-06 and 10-07. Untouched today; still coupled to
  item 9, which is now blocking.
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, **untouched for eighteen
  consecutive sessions**. Today was a naming-and-audit session and again did
  not touch it. Eighteen sessions is no longer a backlog entry; it is a
  statement about what this automation preferentially does, and §7.7 is the
  chapter that should be made to answer for it.

### 11. Addendum: the device dropped mid-session, and the work was nearly lost

Written after the fact, because the §10 health note above was committed before
the interruption and would otherwise read as a clean session.

**The device dropped off the bridge at roughly 05:25 UTC**, immediately after
this entry was committed, with **all five of this repository's commits and the
Design-Verification-Roadmap's first commit sitting local and unpushed**. Three
retries failed. The session sent one notification saying so and stopped, on the
correct assumption that the scratch clone was lost and a later firing would
redo the work from GitHub. **The connection came back**, `$HOME/work` had
survived intact in the session VM, and the first action on reconnect was to
push everything already committed and verify it against the GitHub API before
resuming.

Two things worth keeping:

- **The three-firing redundancy would have worked.** Nothing had been pushed,
  so both repositories on GitHub were untouched — no half-pushed state — and a
  later firing's Step 0 would correctly have found no 2026-10-08 entry and run
  a full session. The design held under its first real test.
- **PUSH EARLIER, and this is a standing correction to how these sessions are
  run.** Five commits accumulated here before anything was pushed, which left
  roughly forty minutes of work in a place only a live bridge could reach. The
  remaining work was pushed in batches instead. A commit that exists only in
  the session VM is not preserved work, and the cost of pushing more often is
  one command.

**Commits this run:** 5 (the census with its transcript; the declaration, the
double-count measurement and the two adjudicated registry rows; the session
note; Chapter 4 §4.6.6 with Chapter 7 §7.8.1e and the §7.9 rewrite; this
addendum). This AUTOMATION_LOG.md entry makes **6**.

---

## 2026-10-09 — The sign was computed and thrown away, and the floor was never a floor

**Session type:** full session. Step 0 at the **04:30** firing (04:35 UTC) found
the device reachable, the "scheduled harsh" folder connected, **neither**
repository carrying a 2026-10-09 entry and no commits since midnight.

**Item worked:** §7.9 item 9, *a signed-carrier FET model* — created 10-07,
promoted 10-08 from sibling to **prerequisite** of the top item 12, and
untouched until today. Closed. Closing it also closed the §4.9 floor item and
partly answered item 11.

### 1. Why item 9 had to precede item 12, more sharply than 10-08 put it

10-08's reason — re-deriving first "would rebuild the branch-sign flaw in a new
equation" — is correct and understates the obstruction. **Both terms on the
right of (4.28) are odd in `n`**: `e·n/C_ox` trivially, and
`E_F(n)/e = A_F·sign(n)·√|n|` because a hole gas's Fermi level lies below the
Dirac point. A relation with two odd sides has an **odd root**, and an odd root
cannot be produced by a solver whose first line is `dV = np.abs(dV)` — which is
the first line of `graphene_potential_declaration.n_quasi_fermi()`. That solver
returns the right *magnitude*; nothing in the repository could return **which
branch the channel is on at a given point**.

`carrier_density()` is where it is lost. It **computes** the sign — `n_eff` is
odd in `dV` because `C_q` is even — and discards it in
`√(n_eff² + n_puddle²)`. Harmless for conductivity, which takes `abs(n)`;
fatal for (4.28).

`n_net_exact()` restores the branch. Validated against exactly known values:

- the root is odd **BITWISE** over 801 points (residual exactly 0);
- substituting it back into (4.28) leaves **7.690e-16** relative over
  `dV ∈ [−2, 2] V`.

**§7.9 item 12 is UNBLOCKED, and for the first time in three sessions the top
item has nothing in front of it.**

### 2. The decomposition, and the identity it exposes

    n_net = n_e − n_h  (electrostatics, signed, odd in dV)
    S = √(n_net² + n_puddle²),   n_e = (S + n_net)/2,   n_h = (S − n_net)/2

| identity | residual | what it means |
|---|---|---|
| `n_e + n_h` == `carrier_density()` | 3.231e-16 rel | the floored quantity is the **TOTAL** density |
| `n_e·n_h` == `n_puddle²/4` | 3.603e-16 rel | the August quadrature floor **is** a mass-action law |
| `n_e(+dV)` == `n_h(−dV)` | **bitwise 0** in `dV` | ambipolar symmetry, machine-checked |

**Contract honoured: no number in Chapter 4 moved.** Two-carrier `σ` equals the
shipped `σ` to 3 ulp; the shipped sweep is re-run after the monkey-patching
with `n_puddle` and `Rc_total` verified restored.

**Literature position recorded honestly rather than favourably.** Live web
search was available; two queries and one fetch, all successful, no
PubMed/PMC attempted per the standing note. [Wiedmann *et al.*, *Coexistence of
electron and hole transport in graphene*, **Phys. Rev. B 84, 115314 (2011)**,
arXiv:1107.3929](https://arxiv.org/pdf/1107.3929) splits net-from-gate and
total-from-disorder exactly as this module does, and argues from the `1/R_H`
divergence at neutrality that `n + p` must stay finite — **independent,
measurement-side support for identifying `n_puddle` as the total**. It imposes
**no** electron–hole product; `n` and `p` are outputs of their Hall fit. So the
mass-action law is a **stronger assumption than the literature makes**,
inherited from a smoothness choice, and is logged as a limitation rather than a
result. Their sample B gives a total puddle density of 8.4e14 /m², i.e. the
shipped `n_puddle` is **5.95×** an externally measured value — noted, not acted
on, since §4.6.7's comparison is internal to this thesis and is the stronger
finding.

### 3. Three exactness checks failed as first written, and none was a code defect

This is the session's methodological content.

**(a) The ambipolar symmetry failed "bitwise", by 2.0 /m².** The symmetry is not
broken — it is bitwise in `dV`. `(V_dirac − x) − V_dirac` is **not**
`−((V_dirac + x) − V_dirac)` in binary floating point: worst case 2.220e-16 V at
`x` = −2 V, where `0.8 − 2.0 − 0.8` returns `−1.9999999999999998`. The `C_ox/e`
lever is 2.395e15 /(V·m²). **Every sweep in this repository is a `linspace` over
`V_g` followed by a subtraction of `V_dirac`**, so any bitwise symmetry claim
about any of them must name its parameterisation; `n_net_of_dV()` exists so
such claims can be made in the variable the model is a function of. The residual
is **accounted for, not thresholded**: feeding the reconstructed `dV` back into
the `dV`-level model reproduces the `V_g`-route species bitwise, which proves
all 2.0 /m² is the argument's and none is the model's. A "≤ 2 ulp" tolerance
would also have passed, and would have been a number tuned until green.

**(b) "`n_e − n_h` recovers `n_net`" failed at 1.772e-14 and passes at
1.5e-16** — the same 1-ulp-of-`S` residual, **119× apart**, because `n_net → 0`
at neutrality while `S` stays pinned at the floor. **The denominator decided
PASS from FAIL.** That is 10-06's ratio-criterion item occurring in a check
written *by the audit*, not found in the audited code — its first instance of
that kind. Design rule added to the module: never reconstruct the net density
from the species pair.

**(c) The `σ` identity failed "bitwise" at 3 ulp** and was **not** repaired by
widening a ulp budget. 09-29's α anchor set the policy: state the identity and
name the slack. The criterion is now 1e-12 relative, justified as six orders
below the four significant figures Chapter 4 quotes and eleven above float
noise — a statement about the chapter, not about the test. 10-02 said byte
identity is not correctness; today adds the converse.

### 4. §4.9's floor item, open since 2026-09-29 and untouched for ELEVEN sessions

| | value | sheet resistance |
|---|---|---|
| Chapter 4: `σ_min = e·μ·n_puddle` | 3.204353e-04 S/sq | **3.1208 kΩ/sq** |
| Chapter 3 §3.3, measured `4q_e²/h` | 1.549618e-04 S/sq | **6.4532 kΩ/sq** |

**Chapter 4 is 2.0678× MORE CONDUCTIVE at neutrality than Chapter 3** — its
channel at the Dirac point is less resistive than any measured graphene sheet.
Reconciling `n_puddle` = 2.418e11 cm⁻² against the shipped 5e11 cm⁻². §3.3's
6.45 kΩ/sq is recomputed here from CODATA (to 0.0032 kΩ/sq) rather than
re-typed, so the comparison cannot drift.

**It is decidable only because item 9 is closed.** `n_puddle` floors a
*magnitude*, and a magnitude cannot say whether it floors the net density, the
total, or each species — three readings differing by factors of 1, 1 and 2 in
`σ_min`, **the same order as the answer**. Item 9 settles it as the total, and
the per-species reading moves the gap the **wrong way** (4.1357×), so the
convention cannot be the explanation.

**Cost, intrinsic and terminal** (`V_ds` = 0.1 V, on-state `V_g` = 3.5 V):

| | shipped | Ch.3-matched | on/off multiplier |
|---|---|---|---|
| terminal (`R_c` = 600 Ω) | 1.2408 | 1.7820 | **1.4362×** |
| intrinsic (`R_c` = 0) | 1.6146 | 2.8024 | **1.7357×** |

`I_off` falls 35.19 % at the terminals, 51.58 % intrinsically; the contacts
attenuate the error by a further **1.2085×**. **10-08's negative-feedback
structure on a completely different parameter**, which upgrades it from a
property of the charge relation to a property of the device topology. §4.4's
numbers stand with the limitation attached. §7.2 undisturbed: 1.78 is no more a
switch than 1.24 — and §7.8.1's reason is *strengthened* by the floor having
been understated.

**Not quietly rewritten:** §3.3's table is untouched and §3.3.1 corrects only
the word "encodes" — `n_puddle` encodes the same *picture* as §3.3's floor, not
the same *number*.

### 5. Item 11 answered, and its premise was wrong

Item 11 asked which other **near-Dirac** regimes the floor had removed. The
floor is not a near-Dirac term:

| `V_g` (V) | \|n_net\| (m⁻²) | `n_puddle`/\|n_net\| | floor adds |
|---|---|---|---|
| 0.90 | 2.357e+14 | 21.22 | +2023.9 % |
| 2.00 | 2.870e+15 | 1.742 | +100.9 % |
| **3.50** | **6.462e+15** | **0.774** | **+26.4 %** |

At the largest overdrive swept anywhere in this thesis the *residual* puddle
density is **77 % of the gate-induced density**; its contribution drops below
1 % only at `dV` = **14.73 V**, 5.46× anything swept. It does not remove
regimes, **it adds a parallel conducting density to all of them**, and
`carrier_density()`'s docstring calling it a Dirac-point regularizer describes a
regime the thesis never enters. The original `d|n|/dV` = 0 question stays open.

**Found by a check written to PASS.** §5 asserted `I_on` is untouched by an
off-state parameter, on the docstring's authority. It failed at **−6.92 %**.
09-28: prose is a detector. **Today: so is a sanity check you expected to pass.**

### 6. Item 13 produced a wrong number on its first use outside its own module

10-08 ranked the `cosh` overflow (`dV` > 18.3493 V) last as *latent*, "no
committed number affected". Above the onset `C_q` is `inf`, so
`C_q/(C_q + C_ox)` is `inf/inf = nan`; `nan > 0` is **False**; and the first
draft of §4.6.7's 1 %-crossover bisection scored the `NaN` branch as "the floor
dominates", walked away from the root and **returned its own 200 V bracket
ceiling as the answer**. Correct value 14.73 V. Caught before commit because the
figure was 74× the swept range. **Latent and harmless are not the same
property.** Ranking unchanged (the output-neutrality reason stands); mitigation
added — callers outside the module must treat non-finite `C_q` as a hard stop.

### 7. Methodological note, continuing the series

09-28: prose is a detector. 09-29: a mutation that does not arrive is
indistinguishable from a system that does not respond. 09-30: a control has to
sit where the failure enters. 10-01: and name what it compares against in a way
that cannot drift. 10-02: when it agrees, that is a fact about two artefacts.
10-03: a PASS/FAIL at zero is silent about magnitude. 10-04: a quantity can be
computed correctly under the wrong name. 10-05: a criterion with a number in it
can be missed in a way that locates the real constraint. 10-06: and decided
correctly on the quantity it names while missing the finding, if that quantity
is a ratio. 10-07: an item can be well-posed, correctly motivated and built
exactly as written while its premise is an artefact of an unnamed variable.
10-08: and closing such an item can show the ambiguity was the smaller half.
**10-09: AN EXACTNESS CHECK THAT FAILS IS A MEASUREMENT OF WHATEVER IT WAS
WRITTEN ON TOP OF.** Three failed today; one measured the repository's sweep
parameterisation, one measured its own denominator, one measured a docstring.
None measured what it was aimed at, and none was a defect in the new code.

The checkable rule: **when an exactness check fails, find out what it measured
before deciding what to change.** Three repairs were available at the moment of
failure — loosen the tolerance, fix the model, or find out — and only the third
produced the findings; the first would have silently accepted all three.
Mechanisable form: **new §7.9 item 15**, sweep every tolerance in this
repository and flag any set *after* a failure rather than derived from the
precision being claimed.

Second rule, from §4: **a one-line comparison can be blocked for eleven
sessions by an undeclared signature in one of its two inputs.** The floor item
was never hard; it was undecidable, because the three readings of its input
differ by the same factor as its answer. 10-01 found the fault was in the
naming. Today's variant: the fault was in the **signature**.

### 8. Environment and automation health

Device reachable and folder connected at the **04:30** firing. `pip install
scipy` succeeded (3 s; the device VM still lacks it) — network healthy for the
sixth session in seven. The whole session ran **on the device VM**; nothing
crossed the bridge. **Live web search available and used**: two searches, one
fetch, all successful. **10-08's push-earlier correction was followed**: every
commit was pushed and verified against the GitHub API before the next one was
written, so at no point did more than one commit exist only in the session VM.

**Not yet covered (candidates for future runs):**

- **RE-DERIVE CHAPTER 4 ON THE EXACT CHARGE RELATION (4.28)** — §7.9 item 12,
  **now the top item with NOTHING IN FRONT OF IT** for the first time in three
  sessions. Worth −0.6023 % in `I_d` at the RF bias. Its prerequisite closed
  today, and it now *blocks* the §4.6.5 rewiring item rather than merely
  outranking it. **The next session should do this unless something larger
  appears.**
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — the floor question that today's work
  did *not* answer. Today established the floor is global rather than
  near-Dirac; the derivative question is separate and untouched. Now the
  cheapest open physics item.
- **AUDIT EVERY TOLERANCE, FLAG ANY SET AFTER A FAILURE** — §7.9 item 15,
  created today. The natural instrument for item 6's magnitude census *and*
  10-06's ratio-criterion item, since a tolerance set after a failure and a
  criterion written as a ratio are the two ways a check stops being a check
  without anyone editing it.
- **MECHANISE the magnitude item** — created 10-03; worked examples on 10-04,
  10-05, twice on 10-06, 10-07, 10-08 and **twice again today** ((b) above, and
  the 2.0678 factor where a PASS/FAIL would have been silent). **Eight modules,
  nine defects, same result. Four sessions overdue**, and item 15 now subsumes
  part of it, which is either the way it finally gets built or the way it gets
  quietly dropped.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 10-06, widened
  10-07, **untouched for the third session**, and today it produced an instance
  *inside the audit itself* rather than in audited code. That is the strongest
  argument yet for building it.
- **MEASURE WHAT EACH DECLARATION COSTS** — §7.9 item 14, created 10-08,
  untouched. Items 2, 3 and 6 all turn on a named quantity; none costed.
- **THE `cosh` OVERFLOW** — §7.9 item 13. Ranking unchanged (last), description
  amended: it is a trap for every new caller, and it caught one today.
- **WHETHER `n_puddle` SHOULD ACTUALLY CHANGE** — created today, and
  deliberately *not* an item yet. Adopting 2.418e15 /m² would move §4.4, §4.7
  and §4.10 and require its own output-neutrality pass; the question of whether
  two floors that are not required to agree *should* be reconciled is a physics
  judgement this automation should not make silently.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created today. The two-carrier form
  `σ = e(μ_e n_e + μ_h n_h)` is built and deliberately run at `μ_e = μ_h`, so
  it moves no number; contact-doped graphene has a real mobility asymmetry and
  this is now a one-parameter experiment away.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** — split out
  of item 9 today. It is now *expressible* and still unverified; it was never
  the same question as the blocker.
- **Why does the peak-`f_T` bias MOVE BRANCH between the two models?** — also
  split out of item 9 today, and no longer coupled to a blocking item.
- **The seven `UNADJUDICATED` census rows** — created 10-08, untouched. Cheap,
  and the census will keep reporting them.
- **A self-heating term** — created 10-06, **third in the RF queue** behind
  items 12 and the §4.6.5 rewiring. Still a limitation to state.
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, **untouched for NINETEEN
  consecutive sessions**. Today was an audit-and-signature session and again did
  not touch it. §7.7 gained a fourth instance today and this is the standing
  one.

**Commits this run:** 4 (the signed-carrier module with its transcript and
figure; the session notes; Chapter 4 §4.6.7 with Chapter 3 §3.3.1; Chapter 7
§7.9 with §7.8.1f). This AUTOMATION_LOG.md entry makes **5**.

## 2026-10-10 — The series factor was a derivative, and the number it moved was two numbers cancelling

**Session type:** full session. Step 0 at the **first firing of the day**
(06:19 UTC) found the device reachable, the "scheduled harsh" folder
connected, **neither** repository carrying a 2026-10-10 entry and no commits
since midnight.

**Item worked:** §7.9 **item 12**, *re-derive Chapter 4 on the exact charge
relation (4.28)* — created 10-08, blocked by item 9, unblocked 10-09, the top
item with nothing in front of it. **Closed.** Closing it unblocked the §4.6.5
rewiring item and opened one new item, already costed.

### 1. The item asked for a number. The finding is that the chapter had the pair of equations mis-described.

Chapter 4 has called (4.29) "the linearisation of (4.28)" in three places
since 10-08 — §4.6.6, the §4.3 annotation, and §7.9 item 12 itself. **It is
checkable in one line and it is false.** A linearisation shares its function's
leading order. Measured log-log exponents `d ln n/d ln dV`:

| `dV` (V) | (4.28) | (4.29) |
|---|---|---|
| 10⁻¹⁰ | **1.99999854** | **1.00000000** |
| 10⁻⁴ | 1.943766 | 1.000000 |
| 10⁻² | 1.274493 | 1.002201 |

In (4.28), `a·n + b·√n = dV`, so `n → (dV/b)²` as `dV → 0` — **quadratic**. In
(4.29), `C_q(0)` is finite — **linear**.

**What (4.29) actually is.** Differentiating (4.28) gives
`d(dV)/dn = e/C_ox + e/C_Q(n)` with `C_Q(n) = e²dn/dE_F = 2e√n/A_F`, hence

    dn/d(dV) = (C_ox/e)·C_Q(n)/(C_Q(n) + C_ox)                      (4.31)

and (4.31) **is** (4.29) with `dn/d(dV)` where (4.29) writes `n/dV`. The
series factor is the charge relation's **exact derivative**, used as though it
were algebraic: a **one-point rectangle rule** for its own integral. A
truncated expansion improves as the argument shrinks; a rectangle rule
degrades where the slope varies fastest, which is near the Dirac point — which
is why 10-08 measured +68 % there.

Checked, not asserted: `e·dn/d(dV)` equals the series combination of `C_ox`
with `C_Q(n)` to **1 ulp**; integrating (4.31) back to the closed-form root
converges at observed order **3.994** with the `h⁴` extrapolation landing
**1.5e-13** relative.

### 2. Two independent errors, separated and sized for the first time

| `dV` (V) | total error of (4.29) | quadrature alone | **argument alone** | sum |
|---|---|---|---|---|
| 0.01 | +68.2193 % | +32.1722 % | **+68.0927 %** | +100.2649 % |
| 1.20 (RF bias) | **+5.2064 %** | **+2.6396 %** | **+4.7416 %** | +7.3812 % |
| 2.70 | +3.4729 % | +1.7523 % | **+3.2193 %** | +4.9716 % |

The second error is `quantum_capacitance()` evaluating `C_q` at the gate
overdrive instead of at `E_F/e` — noticed in §4.6.5's 10-07 annotation and
never costed. **The argument error is the larger single contributor.** Neither
is the explanation alone. And they are **not additive**: the sum overstates by
1.42× at the RF bias and 1.47× at `dV` = 0.01 V, because they compose
multiplicatively in the series factor. A decomposition reported as a sum would
have erred in the flattering direction.

The puddle floor attenuates all of it by **4.3×** (+5.21 % raw → +1.22 % in
what `carrier_density()` returns). That is a property of `n_puddle`, which
10-09 measured to be 2.0678× off Chapter 3 — so it is not reassurance.

### 3. The re-derivation itself

| quantity | (4.29) | (4.28) | change |
|---|---|---|---|
| `n` at 3.5 V (m⁻²) | 8.170466e15 | 8.000036e15 | −2.0859 % |
| `R_ch` at the Dirac point (Ω) | 624.1509 | 624.1509 | **0.0000 %** |
| `L_x` Pd/300/500 (nm) | 115.20/314.17/523.62 | 112.79/307.62/512.70 | −2.0859 % |
| terminal on/off | 1.240756 | 1.230617 | −0.8171 % |
| intrinsic on/off | 1.614570 | 1.581261 | −2.0630 % |
| `I_d` saturated, RF (A) | 8.190217e-05 | 8.140885e-05 | **−0.6023 %** |

The last line **reproduces §4.6.6's published −0.6023 % to four decimal places
from an independently written patch path**. By 10-02's rule that is a fact
about two artefacts, not a validation; what it buys is the assurance that the
(4.29) column is the committed model and not a re-run with drift.

Every qualitative claim of §§4.3–4.4 survives: the 200 nm crossover verdict,
the on/off ratio still between 1 and 2, the resistance ordering. 0.0000 % at
the Dirac point because the floor supplies the whole density there.

**On/off moves less at the terminals (−0.8171 %) than intrinsically
(−2.0630 %) — the THIRD parameter in three sessions** to show that structure,
after the diffusion term (10-08) and `n_puddle` (10-09). Recorded as a
property of the **device topology**, not of any one parameter.

### 4. The item was mis-scoped, and the number it moved was a cancellation

Under (4.31) the small-signal gate capacitance **is** the charge relation's
derivative, so re-deriving the charge forces `C_gs` to move too:

| | (4.29) | charge only | charge + `C_Q` |
|---|---|---|---|
| peak-`f_T` bias (V) | −1.1115 | −1.1992 | −1.1867 |
| `f_T` (GHz) | 20.2787 | 19.8734 | 20.2666 |
| change in `f_T` | — | **−1.9986 %** | **−0.0594 %** |

**The −0.06 % is a cancellation of two ~2 % terms, not a small effect.**
Reporting the charge alone overstates by 34×; reporting only the total calls a
pair of 2 % corrections negligible. §7.4's bound-versus-ranking distinction
gains a third member — *a number that is small because two larger numbers
oppose* — and **nothing in this repository would detect one, because every
check ever written here asks about one quantity at a time.** That is a gap in
the method and it is logged as such.

**The peak's LOCATION moves 6.76 % against 0.0594 % in height, a ratio of
114.** Both peaks stay on the hole branch, which answers §4.6.7's peak-`f_T`
branch-crossing question **NO** for the (4.28)-vs-(4.29) pair and leaves it
open for the pair it was asked about. §7.8.1a untouched for the right reason:
`f_max/f_T` moves **+0.0228 %**, 0.643007 → 0.643154, still below the 1.3–1.4
band a *perfect* gate cannot reach. **Five levers now tried** against the RF
verdict and the binding constraint is still that `transfer_characteristic()`
is a resistor with no output resistance.

### 5. Five checks failed as first written, and today's addition to the series

**(a)** The quadratic-exponent check, swept only to `dV` = 10⁻⁴ V, returned
1.9438 against ±0.01. It measured the **crossover drive**
`dV_x = 2A_F²C_ox/e` = **6.519 mV**: 10⁻⁴ V is 65× below it, where the
correction is O(√(dV/dV_x)) ≈ 12 %.

**(b)** The integral identity failed at 1.1e-07, order 3.49, on a uniform
grid — and measured **the same 6.519 mV**. The integrand turns over inside a
layer of that width at the lower limit and 4001 uniform nodes put 22 in it. A
graded grid gives order 3.994 with 12× fewer points. **Two exactness checks,
different subsections, different purposes, both failing by measuring the same
number** — the strongest evidence yet that 10-09's rule is a rule.

**(c)** "Every number moves by under 2 %" failed at 2.0859 %. The 2 % was a
round number chosen before the run, so it was **deleted**, not widened:
widening to 2.5 % is precisely the tuned tolerance item 15 exists to find. The
magnitudes are now printed with no threshold and the check tests the
qualitative claims. Separately, a guessed 10⁻¹⁰ bound on (b) was replaced by
**the observed convergence order and the absence of a plateau** — a criterion
the method supplies rather than the author.

**(d), (e)** Two checks were wrong in their **prose** and failed on it: one
asserted the capacitance half was the *larger* half (it is 0.970×), one
asserted `f_T` was "two orders below the 100–300 GHz RF requirement" (about
one order, and §7.8.1a does not rest on a frequency threshold at all).
09-28: prose is a detector. Today it detected the author's own checks.

**Today's addition.** 10-09 named three repairs available at the moment an
exactness check fails — loosen, fix the model, find out what it measured.
**There is a fourth, and (c) needed it: delete the criterion and check the
claim instead.** (c) was not a measurement wanting a better threshold; it was
a threshold with no claim underneath it. The checkable form: **a criterion
whose failure would not change any sentence in the thesis is not a check**,
and item 15's sweep should flag those as well as tuned tolerances.

### 6. A delivery fault found on the way

`contact_resistance_crossover.py` opens with `from graphene_fet_model import
carrier_density`, so it holds its own reference bound at import. Rebinding
`graphene_fet_model.carrier_density` does not reach it, and a re-derivation
patching only the owning module would have left §4.4 on (4.29) **while
reporting it as re-derived**. 10-03's frozen-default fault in monkey-patch
form. An **arrival control** (09-30's rule) now precedes the §4.4 numbers.

### 7. Literature

**Live web search available and used: one search, one fetch, both
successful.** No PubMed/PMC attempted, per the standing note.
[arXiv:2209.00388](https://arxiv.org/pdf/2209.00388), the
Pasadas/Jiménez-lineage compact-model review, supports **both** diagnoses as
its own conventions: its `C_q` is `∂Q_net/∂V_c` evaluated at the chemical-
potential shift and never at the overdrive, and its series-capacitance
expression appears as `dV/dV_c = 1 + C_q(V_c)/(C_t + C_b)`, explicitly a
differential relation used for charge partitioning rather than to solve the
static gate-to-channel relation — which it treats as implicit and hands to the
simulator. It also states the Fermi–Dirac form "is not convenient for a
compact model" and uses `C_q ≈ k·c₁·√(1 + (V_c/c₁)²)`, whose **`c₁ → 0` limit
is exactly (4.30)**. So **(4.28) is the T = 0 limit of a standard literature
form**, and it has a closed-form root *because* it drops the thermal rounding
the literature keeps. The convenience and the limitation are the same choice,
and that is now Chapter 4's stated limitation rather than a claim of exactness.

### 8. Environment and automation health

Device reachable and folder connected at the first firing. `pip install scipy`
succeeded (the device VM still lacks it) — network healthy for the seventh
session in eight. The whole session ran **on the device VM**; nothing crossed
the bridge. **Every commit was pushed and verified against the GitHub API
before the next was written**, per 10-08's correction, so at no point did more
than one commit exist only in the session VM.

**Not yet covered (candidates for future runs):**

- **FIX `quantum_capacitance()`'s ARGUMENT** — §7.9 **item 16**, created
  today, and **the first item on that list to arrive already costed**:
  +4.7416 % of the raw charge at the RF bias, larger than the quadrature error
  and larger than anything else open on Chapter 4. **Not latent** — it is
  active in every committed number in Chapters 4, 5 and 6 — and it deviates
  from the literature's own convention rather than being a choice this thesis
  is entitled to. Ranked low only because the module owns committed
  transcripts and the change needs an output-neutrality pass. **It would be
  the first change in this thesis to MOVE committed numbers rather than
  annotate them**, and the floor's 4.3× attenuation is why that is survivable.
  The next session should do this unless something larger appears.
- **A CHECK THAT CAN SEE A CANCELLATION** — created today, out of §4.6.8's
  `f_T` result. Every check in this repository asks about one quantity; a
  number that is small because two larger numbers oppose is invisible to all
  of them. The mechanisable form is a *two-knob* sweep: for each committed
  number, move each of its inputs separately and report whether the combined
  move is smaller than either. Sibling of item 14 (what a declaration costs),
  and the first methodological item here that is about **interaction** rather
  than about one quantity being named or bounded wrongly.
- **THE §4.6.5 REWIRING ITEM**, unblocked today by item 12's closure and with
  nothing in front of it now.
- **AUDIT EVERY TOLERANCE** — §7.9 item 15, created 10-09, **now five
  sessions overdue** and today carrying two worked repairs (delete the bound
  and check the claim; let the method set the criterion) plus a widened remit:
  also flag any criterion whose failure would change no sentence in the thesis.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — untouched for the second session;
  still the cheapest open physics item.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 10-06,
  **untouched for the fourth session**. Today produced no new instance, which
  is the first session since its creation that has not.
- **MECHANISE the magnitude item** — created 10-03, worked examples on seven
  subsequent sessions, **five sessions overdue**. Items 15 and the new
  cancellation item now both subsume parts of it.
- **MEASURE WHAT EACH DECLARATION COSTS** — §7.9 item 14, created 10-08,
  untouched for the second session.
- **THE `cosh` OVERFLOW** — §7.9 item 13, ranking unchanged (last).
- **WHETHER `n_puddle` SHOULD ACTUALLY CHANGE** — created 10-09, still
  deliberately not an item. Item 16 would move committed numbers first, and
  doing both in one session would make neither attributable.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 10-09, untouched. One
  parameter away.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** — still
  expressible and unverified, untouched for the second session.
- **The seven `UNADJUDICATED` census rows** — created 10-08, untouched.
- **A self-heating term** — created 10-06, now second in the RF queue behind
  item 16.
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, **untouched for TWENTY
  consecutive sessions**. Today was a numerics-and-wording session and again
  did not touch it. §7.7 gained a fifth instance today and this is the
  standing one.

**Commits this run:** 4 (the re-derivation module with its transcript and
figure; Chapter 4 §4.6.8 with the in-place annotations to §4.6.6, §4.3 and
`graphene_fet_model.py`; Chapter 7 items 12/15/16 and §7.8.1g; the session
notes). This AUTOMATION_LOG.md entry makes **5**.

---

## 2026-10-11 — The argument was a state, the item ranked last was the blocker, and the floor and the temperature are the same effect twice

**Target:** §7.9 **item 16**, "fixing `quantum_capacitance()`'s argument" —
created 2026-10-10 already costed at +4.7416 % of the raw charge at the RF
bias, the largest cost then open on Chapter 4, and the first item that would
*move* committed numbers rather than annotate them. Yesterday's note said to
do it unless something larger appeared. Nothing larger appeared; what appeared
was inside the item.

**Modules:** `graphene_quantum_capacitance_argument.py` (20 checks, 20 passed,
six against exactly known values), `graphene_quantum_capacitance_overflow.py`
(5/5). **Thesis:** §4.6.9 (a)–(h), §4.6.10, §7.9 items 13/16/17/18, §7.10.
**Figure:** `quantum_capacitance_argument.png`. **Notes:**
`notes/2026-10-11-the-argument-was-a-state-and-the-item-ranked-last-was-the-blocker.md`.

### 1. The item understated the work, not the size

`quantum_capacitance(V)` takes a potential; the dispersion requires it at
`E_F/e = A_F·√n`, and `n` is what `C_q` is being used to compute. **The
correct argument is a state, not a voltage**, so the fix is the
self-consistent solve (4.32), not an edit to a call site. That is why the
largest open cost on Chapter 4 sat for a session as a one-line description.

`E_F/e = 0.0625 V` at the RF bias, **5.2082 %** of the 1.2 V overdrive. The
two candidate arguments put `C_q` a factor **17.9** apart and the series
factor only **2.2966 %** apart, because both are far above `C_ox`. The error
is small because **the device is not quantum-capacitance-limited at this
bias**, not because the argument is nearly right.

### 2. The obvious fix deletes the temperature, and yesterday's number priced it

`C_Q(n) = 2e√n/A_F` is a T = 0 expression. `quantum_capacitance()` is
thermally rounded for a reason its own docstring states. Substituting the
dispersion gives `C_q(0) = 0`, a series factor of zero, and an induced charge
of **exactly zero** at the Dirac point — the floor would supply 100 % of the
density there instead of 99.9 %. Sharper: with the dispersion, `n ≡ 0` solves
`dn/d(dV) = (C_ox/e)·f(C_Q(n))` through the origin, so the initial-value
problem is **degenerate**. With the thermal form it is not. **Temperature is
what makes the charge relation well-posed at neutrality.**

2026-10-10 sized item 16 with that substitution, so **+4.7416 % was the cost
of a fix that is wrong in a second way.** A new failure mode for this
repository's collection: the previous ones were checks that measured the wrong
thing; this was a *cost estimate that priced the wrong repair*.

The form that keeps the temperature (`B_T`: thermal `C_q` at `E_F/e`, solved
self-consistently) reproduces `C_q(0) = PREFACTOR·(kT/e)·ln 4` to **0 ulp**
and reduces to the dispersion as `T → 0`.

### 3. Every number in yesterday's decomposition was inflated

| | vs T = 0 (4.28) | vs 300 K reference | overstated |
|---|---|---|---|
| argument error | +4.7416 % | **+2.6408 %** | 1.80× |
| quadrature error | +2.6396 % | **+0.7603 %** | 3.47× |
| both wrong | +5.2064 % | **+3.0963 %** | 1.68× |

(4.28) is a T = 0 relation and the model runs at 300 K, so every percentage
contained a temperature difference as well as the knob it named. The reference
is the exact integral of the model's own differential relation, solved as an
IVP and **validated to −1.5e-14 against (4.28)'s closed-form root**. The
**ordering survives** — the argument is still the larger error, so item 16 was
ranked correctly — and the magnitudes do not. Annotated in place in four
locations (§4.6.8's table, `graphene_fet_model.py`'s marker, §7.9 items 12 and
16).

**The trap.** Today's 300 K *argument* error (+2.6408 %) agrees with
yesterday's T = 0 *quadrature* error (+2.6396 %) to **0.0012 %**. Different
quantities, different references; read in sequence it looks like a
reproduction. 2026-10-02's rule applied for the first time to a **numerical**
coincidence rather than to a re-run.

### 4. The item ranked last was gating the item ranked largest — item 13 CLOSED

The exact-limit check on `B_T` is `T → 0`, and it could not be run. `η = eV/kT`
grows as `T` falls, `cosh(η)` overflows at **η = 710.4758**, and the series
factor returns `inf/inf = nan`. That is item 13 — added 2026-10-08, judged
latent, **ranked last** because η = 710 needs 18.35 V of overdrive at 300 K.
At T = 1 K and the ordinary RF bias, η = 719.

> **Latency had been judged against the inputs the model takes, not against
> the inputs its validations take.**

Now a standing caution on any latency claim and part of item 15's remit. Fixed
with the identity `ln(2(1+cosh η)) = |η| + 2·log1p(exp(−|η|))` from
`2(1+cosh x) = 4cosh²(x/2)`, installed in a new `_thermal_rounding()`. With it
the limit runs: `B_T` reproduces the dispersion fixed point to **7e-16** at
T = 10⁻³ K, against `nan` before.

### 5. The neutrality pass measured the wrong thing three times

1. **Literal constants.** `scipy`'s `ħ = h/2π = 1.05457181764615647e-34`; the
   literal `1.054571817e-34` differs at the tenth digit, ~1.2e-9 once squared
   in the prefactor. The pass reported **3983 of 4003** values differing by up
   to **1.1e7 ulp**. Sixth instance of 2026-10-09's rule, and the first in
   which the substituted quantity was a *physical constant* rather than a model
   parameter.
2. **Association.** `prefactor·(kT/e)·log(…)` against the module's
   `(kT/e)·log(…)` then `prefactor·…`. Floating multiplication is not
   associative; the regrouping alone moved **118 of 400** values by up to 2 ulp.
3. **A silent no-op.** The patch meant to fix (2) matched nothing and reported
   success, so draft 2 ran again and was nearly read as confirmed.

> **A function that "restores the old spelling for comparison" is a
> re-implementation, and every detail of it not copied exactly — the
> constants, the grouping, the order of operations — becomes part of the
> difference it reports.**

The repair is not a better transcription but making the changed thing
**substitutable**: the pass now swaps one function and leaves the rest of the
call path as shipped code. That form reproduces the direct file-against-file
result exactly — **4 of 400, 1 ulp** — where all three transcriptions did not.

Stated as a measurement, not as neutrality: `C_q` moves by 1 ulp at 4 of 400
biases within 0.093 V of the Dirac point; **no committed quantity moves by a
bit** over 3603 values.

### 6. Two self-corrections inside the session

**(a) A criterion I chose, which failed.** "≤ 1 ulp" held on 11 hand-picked
values and failed at **2 ulp** on 401. Not widened — 2026-10-10's fourth
repair applies, and for an *identity* a ulp count cannot answer the question
that matters. `mpmath` at 60 digits adjudicates, and **the answer is not the
expected one**: the installed spelling is **strictly worse at 81 of 401
samples** and better at only 45, with a slightly larger mean error (5.162e-17
vs 4.643e-17). All it has is the smaller worst case (1.527e-16 vs 1.793e-16).
**Accuracy is not a reason to install it**; being finite where the other is
not is the only reason, which is the only reason item 13 ever gave. Recorded
in §4.6.10 because today's first commit message states the ≤ 1 ulp figure and
a pushed message cannot be amended.

**(b) An arrival control that failed on a healthy run.** One counter,
incremented only in the capacitance patch, gating a section that exercises the
*charge* patch and never asks for `C_gs`. Zero capacitance calls was the right
answer and the control called it a fault. 2026-10-10's rule in its other
failure mode: **a control that counts a different arrival from the one it
gates reports a fault that does not exist.** One counter per patched name now,
each checked where that name is used.

### 7. What moves

| | (4.29), `C_q(dV)` | `C_q(E_F/e)` | change |
|---|---|---|---|
| raw charge at `dV` = 1.2 V (m⁻²) | 2.86980420e15 | 2.80477739e15 | **−2.2659 %** |
| `n` at 3.5 V (m⁻²) | 8.170466e15 | 8.087032e15 | **−1.0212 %** |
| `R_ch` at the Dirac point (Ω) | 624.150907 | 624.150907 | **0.0000 %** |
| `L_x` Pd / 300 / 500 (nm) | 314.17 / 314.17 / 523.62 | 310.96 / 310.96 / 518.27 | −1.0212 % |
| terminal on/off | 1.240756 | 1.235801 | −0.3993 % |
| intrinsic on/off | 1.614570 | 1.598207 | −1.0135 % |
| `I_d` at the RF bias (A) | 8.190217e-05 | 8.169499e-05 | **−0.2530 %** |

The floor attenuates the raw move by **0.451**. Every qualitative claim of
§§4.3–4.4 survives. **The correction moves numbers and no verdicts.**

`f_T`: charge alone **−1.0911 %**, the consistent `C_gs` **+1.8976 %**, total
**+0.8065 %** — opposite signs, the **second independent instance** of the
structure 2026-10-10 recorded on one instance as a property of the device
*topology*. That claim now has what it needed. `f_max/f_T` moves **+0.0009 %**:
§7.8.1a untouched by a **sixth** lever, binding constraint still that
`transfer_characteristic()` is a resistor with no output resistance.

**Item 16 does not close.** The default is still (4.29). Flipping it moves
every committed number in Chapters 4, 5 and 6 and needs its own session.

### 8. The two-knob check's first use fails in the direction it cannot see

| reference | argument only | quadrature only | sum | both | both/sum |
|---|---|---|---|---|---|
| 300 K, consistent | −2.2659 % | −0.4418 % | −2.7077 % | −3.0033 % | **1.109** |
| T = 0 (2026-10-10) | −2.4398 % | −0.4418 % | −2.8816 % | −4.9488 % | **1.717** |

The knobs **reinforce**, so a check testing whether the combined move is
*smaller* than either would not have flagged it. It must report the **ratio**
and flag any departure from 1. And **the interaction measure inherits the
reference error of whatever it is differenced against**: 0.717 against T = 0
versus 0.109 at 300 K, an inflation of **6.6×**. A two-knob check must be run
at the model's own operating point. Amended, not closed → **item 17**.

### 9. Literature, and a new item from it

**Live web search available and used: one search, one fetch, both
successful.** No PubMed/PMC attempted, per the standing note.

A **second, independent** lineage —
[the Southampton bilayer-GFET model](https://eprints.soton.ac.uk/366757/1/GFetIEEEtran.pdf)
— confirms §1 from a different source than 2026-10-10's
[Pasadas/Jiménez review](https://arxiv.org/pdf/2209.00388): `C_q` is a
function of the **surface potential**, never of the overdrive, and the layer
potentials form a **coupled system solved together**. Item 16's claim is now
held by two groups rather than inferred from one.

But it regularises `C_q` at neutrality with **disorder, not temperature**:
`C_q,min = e²√n₀/(√π·ħ·v_F)`, with temperature never entering its `C_q`. And
at this thesis's `n_puddle` the two are **the same size**: 9.710793e-03 against
8.437399e-03 F m⁻², ratio **1.1509**. They agree exactly at `n₀ = 3.77e15 m⁻²`
(0.755×) or at **T = 345.28 K** (1.151×), so this is not an accident of one
parameter — **SiO₂-supported graphene at room temperature sits where disorder
and temperature regularise `C_q` equally.**

**This thesis applies both**, and the literature form shows they are two
spellings of one effect at this `n₀`. Using both is a **double count worth
about a factor 2 in `C_q(0)`**, never costed → **item 18**, the first item on
that list to arrive from a literature source rather than from an internal
audit. It would not hold for suspended graphene or at 4 K, which is why
neither convention is obviously wrong and why the overlap is easy to miss.

### 10. Environment and automation health

Device reachable and folder connected at the first firing. `pip install scipy`
succeeded again (the device VM still lacks it); `mpmath` was already present
and served as the 60-digit adjudicator. The whole session ran **on the device
VM**; nothing crossed the bridge. **Every commit was pushed and verified
against the GitHub API before the next was written.**

**Not yet covered (candidates for future runs):**

- **COST THE DOUBLE REGULARISATION** — §7.9 **item 18**, created today, and
  the **first item ever to arrive from a literature source** rather than from
  an internal audit. The thermal rounding and the puddle floor are the same
  physical effect to within 15 % at the operating point, and the thesis applies
  both. Cheapest route in: compare `C_q(0)` and the Dirac-point density under
  each regularisation alone against both. Ranked first because it is the only
  open item that could change a *physical* claim rather than a numeric one, and
  because it reframes the standing `n_puddle` question.
- **FLIP THE DEFAULT TO `B_T`** — item 16, open. Costed on the right form now
  (−2.2659 % raw, −1.0212 % after the floor, no verdict moved) and unblocked.
  Needs a session of its own for transcript regeneration across Chapters 4, 5
  and 6. Second only because item 18 might change what `B_T` should even be.
- **A THIRD KNOB: THE ROUNDING FORM** — created today out of §4.6.9's 2×2,
  which holds the rounding fixed. Fermi–Dirac integrals have never been
  compared against `ln(2(1+cosh η))` at all. Item 18 is the same question from
  the other side.
- **BUILD ITEM 17 AS A SWEEP** — one worked instance and two requirements
  (report the ratio; run at the model's operating point), still no sweep over
  committed numbers.
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, now **six sessions
  overdue**, today carrying a sixth worked instance: a ≤ 1 ulp bound chosen
  before the run, which failed, and whose repair was to let 60-digit
  arithmetic set the criterion.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — untouched for the third session;
  still the cheapest open physics item.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10 with nothing in
  front of it, untouched for a second session.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  **untouched for the fifth session**, and today produced a new instance:
  `f_max/f_T` used again as the structural verdict's instrument.
- **MECHANISE the magnitude item** — created 2026-10-03, **six sessions
  overdue**; items 15 and 17 subsume parts of it.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, untouched for the third
  session.
- **WHETHER `n_puddle` SHOULD ACTUALLY CHANGE** — still deliberately not an
  item, and item 18 is now the reason to look at it.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** — still
  expressible and unverified, untouched for the third session.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item, **untouched for TWENTY-ONE
  consecutive sessions**. Today was again a numerics-and-method session. §7.7
  gained a sixth instance today and this is the standing one.

**Commits this run:** 4 (the item-16 module with its transcript and figure;
item 13's installation with its neutrality pass; §§4.6.9–4.6.10 with the
in-place annotations and §7.9's revisions; the literature section, item 18 and
the session notes). This AUTOMATION_LOG.md entry makes **5**.

---

## 2026-10-12 — item 18 closed by being wrong twice, and the "factor 2" was unreachable in principle

**Session type:** one item, costed and closed. Three commits plus this entry.
Device reachable and folder connected at the first firing; `pip install scipy`
succeeded again (the device VM still lacks it); the whole session ran on the
device VM and nothing crossed the bridge. Every commit was pushed and verified
against the GitHub API before the next was written.

**Item taken:** §7.9 **item 18** — "costing the double regularisation of `C_q`
at neutrality". It was ranked **first** for this session because it was the
only open item that could change a *physical* claim rather than a numeric one,
and because it was the **first item on that list ever to arrive from a
literature source** rather than from an internal audit.

It closes. It closes by being wrong twice over, with a third finding that puts
an upper bound on the entire class of error it described.

### 1. In the committed relation the two regularisers never meet

The premise was "this thesis applies both". For (4.29) that is **false**, and
the reason is structural rather than numerical.

`C_q`'s argument is the gate overdrive `dV`, and `dV` is never floored, so
`n_puddle` cannot reach `C_q`. The converse is sharper: `n_eff =
(C_ox·dV/e)·f(C_q)` has `dV` as a **factor** and `f ∈ (0,1)`, so at `dV = 0`
the product is zero whatever `f` returns. Measured on the shipped function with
four `C_q` conventions substituted — the thermal rounding, the T = 0
dispersion, the Southampton floor, and this thesis's own `C_Q(n_puddle)` —
`n_eff(0)` is **`0.0` bitwise in all four**, and `carrier_density(V_dirac)`
returns **exactly** `n_puddle`.

Both halves are checked against **exactly-known values**, `0.0` and
`n_puddle`, not against tolerances — which is what 2026-09-17 established and
what a plausible-range check would have passed under either premise.

> **Applying one regulariser to each of two different quantities is a reuse of
> one parameter for two jobs — a modelling choice to defend. It is not a
> double count, which is an arithmetic error.**

The generalisable form, new today: **a cross-model comparison establishes that
two conventions exist; it does not establish that one model contains both.**
Item 18 was created by comparing two *literature forms* without first asking
which quantity each term in **this** model reaches.

### 2. Where they do meet, the cost is +2.5262 % and not a factor 2

One spelling does put both regularisers on one quantity: form B_T (item 16)
with the floor inside `C_q`'s argument, `C_q(E_F(√(n² + n₀²))/e)`. It is the
most double-counting spelling available and arguably the consistent one.

| `C_q(0)` in form B_T | thermal OFF | thermal ON |
|---|---|---|
| floor outside `C_q` | 0.000000000e+00 | 8.437398908e-03 |
| floor inside `C_q` | 1.942158579e-02 | **1.991221401e-02** |

**both/max = 1.025262**, both/sum = 0.714750; thermal on top of the floor
**+2.5262 %**, floor on top of thermal +135.9994 %.

> **The "factor 2" was the ratio of the two candidate values. That ratio
> prices *choosing wrongly between* the conventions; it is not the cost of
> *applying both*.**

First instance of the 2026-10-06 ratio item — *census every criterion that is
a ratio* — **in which the ratio was not a criterion but an estimate of a
magnitude**. Five sessions open as a question about verdicts; it was a
question about magnitudes too.

And the structural reason, stated as a rule: **two regularisers entering
through a monotone function of a quadrature sum compose max-like, not
additively — the larger swallows the smaller.**

Downstream, at the model's own biases (item 17's second requirement):

| `ΔV` | neither | thermal | disorder | both | both/max |
|---|---|---|---|---|---|
| 0 V | 5.000000e+15 | 5.000000e+15 | 5.000000e+15 | 5.000000e+15 | 1.000000 |
| 1.2 V | 5.730515e+15 | 5.732955e+15 | 5.741244e+15 | 5.741709e+15 | 1.000081 |
| 2.7 V | 8.085749e+15 | 8.087032e+15 | 8.095579e+15 | 8.096222e+15 | 1.000079 |

At the Dirac point **all four agree bitwise** — the quadrature floor on the
*result* overwrites every `C_q` convention, so **the quantity item 18 was
created to protect is the one quantity in the model that cannot move.**
`σ_min = 3.204353268e-04` S/sq in all four, identical to the last bit: no
committed conductivity, resistance or on/off number in Chapters 4, 5 or 6 is
exposed at neutrality. At the RF bias the whole double count is **+0.1527 %**
in `n`, against −1.0212 % for item 16 itself and −0.60 % for the exact charge
relation. **Real but small, and smaller at the operating point than at
neutrality — the opposite of how it was ranked.**

Second session running in which **a cost estimate priced something other than
the repair it named** (2026-10-11's lesson, one level up).

### 3. The floor's own spelling carries an exact factor 2

The Southampton floor `C_q,min = e²√n₀/(√π ħ v_F)` is **exactly half** this
thesis's own DOS-exact `C_Q(n₀) = 2e²√n₀/(√π ħ v_F) = 2e√n₀/A_F` — **0.000 ulp
from ½** with one consistent `ħ` — because graphene's `D(E) = 2|E|/(π ħ² v_F²)`
carries the `g = 4` spin–valley degeneracy and the literature form does not.

So yesterday's ratio **1.1509** is `kT` against the halved spelling; on this
thesis's own convention it is **2.3018**, and "disorder and temperature
regularise `C_q` equally here" weakens from *the same size* to *the same
order*. Yesterday's figure is **annotated in place and not rewritten**, in
`graphene_quantum_capacitance_argument.py` §9 and in §4.6.9(h).

### 4. Literature: the source never defines `n₀`, and its own composition law bounds the error

**Live search and fetch available and used: one fetch, successful.** No
PubMed/PMC attempted, per the standing note.
[`GFetIEEEtran.pdf`](https://eprints.soton.ac.uk/366757/1/GFetIEEEtran.pdf)
re-read for the one residue arithmetic could not reach: what `n₀` means there.

**It does not say.** At zero surface potential it states only that "the channel
has been reported to have a charge density, `n₀`", cites two references, and
never calls it total, per-species, or a disorder parameter. **No degeneracy
factor appears anywhere and it gives no DOS at all** — it goes straight from
`n` to `C_q`.

> **The comparison is blocked not by disagreement but by the source having
> left a symbol undefined and delegated it to a citation.** There is no stated
> convention to convert; the residue stands and both figures stay in the record.

What the reading *did* settle, in the measurement's favour. Their Eq. (3),

```
C_q(φ_s) = (C_q,min² + 2h²)/√(C_q,min² + h²),      h = C_q,var/2
```

returns `C_q,min` **bitwise** at `φ_s = 0`, tends to `C_q,var` to twelve
digits when that dominates, and is therefore a **smooth max and not a sum**.
Its excess over pure max is maximal exactly where the two terms are **equal**:

```
sup C_q/max(C_q,min, C_q,var) = 1.5/√1.25 = 3/√5 = 1.341640786
```

measured at `C_var/C_min = 1.000000` over 200 001 points against that closed
form.

**(a) The literature applies `n₀` twice itself** — Eq. (2) floors `C_q`,
`n = √(n₀² + n*²)` floors the density, no temperature anywhere. So this is a
comparison between two **allocations** of regularisers to quantities, not
between a model that double-counts and one that does not. This thesis spends
`kT` on `C_q` and `n₀` on `n`; the source spends `n₀` on both.

**(b) And yesterday's intuition was right about *where* the worst case is and
wrong about *how big*.** A smooth-max composition is furthest from pure max
exactly when its terms are equal — precisely the condition §9 identified at
this `n₀` and this `T`. That reasoning was sound. But the excess there is
**+34.164 %**, and this thesis sits at **+2.5262 %** because its terms are
2.30:1 rather than 1:1.

> **A quadrature-composed regulariser pair cannot exceed the larger of its two
> terms by more than 34.2 %. The "factor 2" was unreachable in principle, not
> merely unmeasured.**

An upper bound on the whole class of error the item named, derived from the
*source's* composition law rather than this thesis's — and quoted as context
for the item's estimate, not as a bound on this thesis's own nested spelling,
whose +2.5262 % is measured.

Also: the source's fitted **`n₀ = 1.2e16 m⁻²` is 2.4×** this thesis's
`n_puddle = 5e15`. A second literature anchor on the standing `n_puddle`
question, pointing the same way as §4.6.7's: this floor may be **low**.

### 5. Item 17 amended again — and this time the amendment is to the null

Item 17's prescribed form reports combined/(sum of separate) and flags any
departure from 1. Today's pair gives **0.7148** against that null, a departure
of 28.5 %, which it would flag as a large interaction. **It is not one.**
Against `max(alone)` the same measurement is **1.0253**, a departure of 2.53 %.

> **A two-knob check must state which composition law it is testing against,
> because the same numbers read as a 28.5 % interaction or a 2.5 % one
> depending on a null that was never chosen explicitly.**

Three requirements now, still unbuilt as a sweep: report the ratio, run it at
the model's own operating point, name the null.

### 6. And the module's first draft failed its own substitutability check

Recorded because the ordering is new: the check was written **before** the
measurement, not during a neutrality pass afterwards — the first time that
ordering has held here.

Draft 1 built its switches on `graphene_diffusion_current_model`'s literal
`E_CHARGE`/`HBAR` while `quantum_capacitance()` uses `scipy`'s, and reproduced
the shipped `carrier_density(2.0)` only to **2375 ulp**. **Eighth** instance of
the 2026-10-09 rule; **second** with a physical constant.

The repair is 2026-10-11's repair and not a better transcription: each switch
**substitutes one shipped name** — `gfet._thermal_rounding`, `gfet.n_puddle` —
and leaves the rest of the call path as shipped code, so form A agrees with
the shipped function at **0 ulp** by construction rather than by care. One
arrival counter per patched name, each checked at a call that reads that name.

One caveat kept in the output rather than smoothed over: the `T → 0`
validation reports **exactly** zero departure, which is *weaker* than `1e-16`
would be, because at `η ~ 10⁷` the term `2·log1p(exp(−η))` underflows to
exactly `0.0` and the two expressions reduce to the same floating-point
operations. The informative version of that check is §4.6.10's 60-digit
adjudication at moderate `η`.

### 7. Environment and automation health

22/22 checks in the new module; `graphene_quantum_capacitance_argument.py`
re-run after annotation, 20/20 and the check set unchanged against its
committed transcript. Nothing deleted anywhere: yesterday's §9 numbers, item
18's premise text and item 17's two earlier forms are all annotated in place.

**Not yet covered (candidates for future runs):**

- **FLIP THE DEFAULT TO `B_T`** — item 16, open, and now **unblocked on both
  sides**. Costed on the right form (−2.2659 % raw, −1.0212 % after the floor,
  no verdict moved), its blocking defect gone since 2026-10-11, and the item
  that was ranked in front of it closed today. It needs a session of its own
  for transcript regeneration across Chapters 4, 5 and 6. **Ranked first:
  nothing is now in front of it.**
- **WHAT `n_puddle` ACTUALLY MEANS** — promoted today from "deliberately not
  an item" to the top *physics* question on this list after the top one.
  Reached from three directions now: §4.6.7's floor being 2.07× Chapter 3's,
  item 18's literature use of the same symbol for DOS broadening, and today's
  fitted `1.2e16 m⁻²` at 2.4× this thesis's value. The cheapest route in is
  the two references the Southampton paper delegates `n₀` to, which were not
  followed today.
- **WHETHER THE FLOOR BELONGS INSIDE `C_q`'s ARGUMENT** — created today. Form
  B_T with `floor_in_cq=True` was built to *measure* the double count, and
  whether it is the right physics was never asked. It is the same unexamined
  question as the one above, from the model's side: the DOS should be evaluated
  at the density the channel has, and whether the puddle density is part of
  that density is exactly what "what `n_puddle` means" decides. Costed at
  +0.1527 % in `n`, so it is cheap to flip and cheap to leave.
- **A THIRD KNOB: THE ROUNDING FORM** — created 2026-10-11, untouched.
  Fermi–Dirac integrals have never been compared against `ln(2(1+cosh η))` at
  all. Today's smooth-max finding is a reason to expect the comparison to be
  small: the rounding form is itself a smooth max between `kT ln4` and `|V|`.
- **BUILD ITEM 17 AS A SWEEP** — now carrying **three** requirements and two
  worked instances, still no sweep over committed numbers. The third
  requirement (name the null) arrived today and is the one that changes what
  the sweep must *report*, not just where it must run.
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, now **seven
  sessions overdue**. Today adds a seventh worked instance from the other
  direction: a validation whose departure was *exactly* zero, which looked
  stronger than a 1e-16 and was weaker.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  **untouched for a sixth session**, and today produced its first instance
  outside its own stated scope: a ratio read as a magnitude rather than used
  as a criterion. The census as specified would not have caught it.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — untouched for the fourth session;
  still the cheapest open physics item.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, untouched for a
  third session.
- **MECHANISE the magnitude item** — created 2026-10-03, seven sessions
  overdue; items 15 and 17 subsume parts of it.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, untouched for a fourth
  session.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** — still
  expressible and unverified, untouched for a fourth session.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item in the thesis, and now
  **untouched for TWENTY-TWO consecutive sessions**. Today was again a
  numerics-and-method session, with one literature fetch that served the
  numerics rather than the physics. §7.7's unequal-scrutiny finding gained a
  seventh instance today and this is the standing one.

**Commits this run:** 3 (the module with its transcript and figure; §4.6.11
with items 17 and 18 annotated and yesterday's §9 annotated in place; the
literature reading with the smooth-max bound, §4.6.11(f) and the notes file).
This AUTOMATION_LOG.md entry makes **4**.

---

## 2026-10-13 — The harness could not reach the second module, so the item ranked first was costed with half of Chapter 4 switched out

**Target:** §7.9 **item 16**, ranked first since 2026-10-12, carrying the note
that flipping the default "moves every committed number in Chapters 4, 5 and
6 and needs its own session". This was that session. **The default is not
flipped.** The blast radius was measured first, and the sentence is wrong in
both directions: **smaller across the thesis, larger inside Chapter 4**, and
inside Chapter 4 it reaches a verdict.

`graphene_contact_argument_blast_radius.py`, **24/24**, with
`contact_argument_blast_radius_output.txt` and its figure. Thesis §4.6.12,
plus annotations in place to §4.7 and §4.10.

### 1. The hole: a second module, holding its own import-time reference

`qca.on_consistent_argument()` rebinds `gfet.carrier_density`,
`gfet.quantum_capacitance` and `crx.carrier_density` — the third added
2026-10-10 after exactly this class of fault. It does not rebind
`graphene_contact_doping_model.quantum_capacitance`, bound at import time on
that module's line 50 and called by `contact_edge_carrier_density()` on line
105 as `quantum_capacitance(dV, T=T)` with `dV = -(W_metal - W_graphene)` —
**the item 16 defect verbatim, in a second module, against a work-function
potential instead of a gate overdrive**.

Measured, not read from the source: under the patch that binding is reached
**11 times** and **exactly 0** of those calls arrive at the corrected
function. So **−2.2659 % raw, −1.0212 % in `n`, −0.2530 % in `I_d`** were all
measured with §§4.5, 4.7 and 4.8 **switched out of the circuit**.

> **Nothing errored and no number moved, which is the whole difficulty: a cost
> table built through this harness reports a model that was ABSENT as a model
> that is INSENSITIVE.**

The 2026-10-10 delivery fault was repaired for `crx.carrier_density` and the
sweep that found it was never re-run across the other modules. **A repair to
one binding is not a repair to the class.**

### 2. And Chapters 5 and 6 were never in the circuit at all

Same counters, on **nine committed summary routines** —
`graphene_interconnect_model.summary_numbers()` and
`self_check_parallel_limit()`, the collection model's `summary_numbers()` and
`modeled_eqe_by_metal()`, the signed-carrier model's three tables, the
two-contact model's two. **0 charge-path calls between them.** Chapter 5
imports nothing from the charge model; Chapter 6 imports `mu`,
`METAL_WORK_FUNCTIONS`, `W_GRAPHENE`, `lambda_decay` — parameters, not the
charge path.

> **Item 16 is a Chapter 4 operation that reaches two models, not a
> three-chapter one that reaches one.** The scope sentence was written from
> the chapter list and never from a call graph.

### 3. The correction at the contact is 4.7× to 12.2× the channel's

`C_interface = ε₀/3.5 Å = 2.529768e-02 F m⁻²` is **65.93×** `C_ox`, and `C_q`
is comparable to it rather than small against it: the series factor runs
**0.555–0.806** at the contact against **0.999** in the channel at 3.5 V. The
channel's near-unity factor is what attenuated item 16 by 0.451. And the
potential is 0.15–1.15 V of work-function difference rather than a few hundred
mV of overdrive — `E_F/e` lands at 0.315 V where `ΔV` is 0.620 V for Pd.

| | Ti | Cu | Ni | Au | Pd | Pt |
|---|---|---|---|---|---|---|
| Δn_contact | −6.2578 % | −4.7696 % | −12.4233 % | −12.4878 % | **−12.4964 %** | −11.8738 % |

against **−1.0212 %** in the channel. (Cr has `ΔW_F = 0` exactly and both
forms return `-0.0`.)

### 4. And it changes a verdict — the one §4.10 had already flagged as fragile

§4.7 reports `R_transmission = Rc,measured − R_extra`, three negative and
**only Pd positive** at **+50.949 Ω·µm** — the only evidence the additive
decomposition might survive. §4.10's condition-number audit found that row
**flips sign on a 9.6 % error in `R_extra`**, Cu's on 81.1 %.

| | Cu | Ni | Au | Pd |
|---|---|---|---|---|
| ΔR_extra | −1.2583 % | +20.8352 % | +39.2281 % | **+47.6925 %** |
| R_tr, form A | −787.461 | −360.462 | −90.159 | **+50.949** |
| R_tr, form B_T | −775.236 | −533.490 | −329.120 | **−203.276** |

**Five times the margin.** Under form B_T the additive decomposition fails on
**4 of 4** metals.

> **§4.10's 9.6 % was a margin against INPUT error. It is consumed by a
> MODELLING error that was already in the committed code. "Pd is the fragile
> one" was right about which row and wrong about what would break it.**

Unchanged, and checked rather than assumed: the metal ordering by residual
(Cu < Ni < Au < Pd under both forms) and the n/p type assignment for all
seven metals. And the direction is the physical one — `n_contact` down means
a less strongly doped, more resistive junction, so `R_extra` up; Cu is the
exception at −1.26 % because its contact density sits *below* the 2.0e16 m⁻²
bulk reference, so reducing it moves the profile toward the bulk.

§4.7's conclusion is **strengthened**: a four-of-four failure fits its stated
TLM double-count explanation better than three-of-four did.

### 5. Three instrument defects found on the way, all kept in place

**(a) This module's own first draft counted the wrong arrival.** It counted
calls to the *name* `cd.quantum_capacitance` and asserted zero under the
patch. They are 11 whether or not anything was substituted into that name.

> **A control that counts calls to a NAME cannot see whether the right
> FUNCTION is bound to it.** The counter has to live inside the replacement.

**Ninth** instance of the 2026-10-10 rule, **second** in its "counts a
different arrival from the one it gates" direction.

**(b) Item 13's overflow is closed in the shipped function and open in the
harness.** `qca.cq_thermal` is still the committed `ln(2(1+cosh η))` spelling
and still overflows at `η = 710.4758`. A fixed point built on it takes the
*electrostatic* density as its first iterate, and that iterate grows with the
series capacitance — so the solver returns `nan` for any `C_s` above
≈ 100 F m⁻². `C_interface` is 2.53e-02 and the worst first-iterate `η` here is
**19.2**, so no number is affected. Item 13's lesson for a **third** time:
*a defect unreachable from the device's operating range is reachable from a
check on the device's model.* Form B_T is therefore built on the **shipped**
`gfet.quantum_capacitance`, which also removes a **1.2254e-09** constants gap
(prefactor ratio `qca/gfet = 1.000000001225439`) from the A→B_T comparison —
so the two forms differ in the argument and in nothing else by construction,
not by tolerance. **Ninth** instance of the 2026-10-09 substitutability rule,
**third** with a physical constant.

**(c) A tolerance chosen before the run failed.** The large-`C_s` limit check
demanded independence to `1e-9` and measured `5.254e-04`. The limit is
approached at **first order** in `C_q/C_s`, so no fixed tolerance states what
it claims. Replaced by the **rate**: the departure must fall by ten per decade
of `C_s`, measured **9.9933** and **9.9993**. **Item 15's eighth worked
instance**, repaired the 2026-10-11 way — let the structure set the criterion.

### 6. Validations against exactly known values

The zero null at `ΔW_F = 0` (exactly `-0.0`, both forms); bitwise
electron/hole antisymmetry `n(W_G+δ) = −n(W_G−δ)` at 0 ulp, both forms; an
**identity-patch neutrality pass** in which every counted binding is wrapped
with a wrapper that returns the shipped value and every number must be
bit-identical — it is, so the counters in §§1–2 are not the perturbation;
reproduction of 2026-10-11's **−1.0212 %** in `n` at 3.5 V to four decimals
(−1.021162 %) through this session's own solver; and the solver proved to BE
`qca._fixed_point` by giving 0 ulp against `qca.n_arg_EF_thermal` when handed
qca's own `C_q`. The fixed-point residual is **0.000e+00** at all seven work
functions.

Every patch substitutes **one shipped name** and leaves the rest of the call
path as shipped code: form A at the contact is the shipped
`contact_edge_carrier_density()` with no patch at all, form B_T is the same
function with `cd.quantum_capacitance` rebound, and the solver is
`qca._fixed_point` run with `gfet.C_ox` substituted by `C_interface` — because
`qca.series_factor` reads that name at call time, so the contact solve and the
channel solve are literally the same function.

### 7. Environment and automation health

Device reachable and folder connected at the first firing. `pip install scipy`
succeeded again (the device VM still lacks it). No web search was needed or
attempted: today's work was entirely internal. The whole session ran **on the
device VM**; nothing crossed the bridge. **Every commit was pushed and
verified against the GitHub API before the next was written.** Nothing was
deleted: §4.7's table, its "only Pd" sentence and §4.10's 9.6 % are all
annotated in place.

**Not yet covered (candidates for future runs):**

- **FLIP THE DEFAULT TO `B_T`** — item 16, still open, still **ranked first**,
  and now **correctly scoped for the first time**. Today removed the reason it
  kept being deferred (nobody knew what it touched) and added two
  preconditions nobody knew about: (i) `contact_edge_carrier_density()` must
  flip in the SAME COMMIT as `carrier_density()`, or Chapter 4 holds a
  corrected channel against an uncorrected contact; (ii) §4.7's prose and
  §4.10's annotation need **rewriting, not re-running** — "only Pd gives a
  small, plausible positive residual" becomes false. It is a Chapter 4
  operation only: Chapters 5 and 6 are measured at 0 arrivals.
- **SWEEP EVERY IMPORT-TIME BINDING OF EVERY PATCHED NAME** — created today,
  and the generalisation of the only defect that has now appeared in three
  modules. `crx.carrier_density` was repaired 2026-10-10 as a one-off;
  `cd.quantum_capacitance` was found today by hand. The repo has four
  module-level bindings of the two charge-path names and nobody has ever
  enumerated the bindings of `sheet_conductivity`, `n_puddle`, `C_ox`, `mu`,
  `lambda_decay` or `W_GRAPHENE`, any of which a future audit will want to
  substitute. Cheapest route: an AST pass over every module, emitting the
  binding table, with a test that every name a context manager rebinds has
  ALL of its bindings rebound. **Ranked second**: it is what makes item 16's
  flip verifiable rather than hopeful.
- **THE CONTACT MODEL'S OWN APPROXIMATIONS** — created today out of
  §4.6.12(f). A 12.5 % argument correction was found in a model that assumes a
  vacuum-like 3.5 Å gap with no interface dipole layer and one
  `lambda_decay = 250 nm` for every metal. Both are plausibly larger than
  12.5 % and neither has ever been costed. §4.7 already names the second as a
  candidate explanation for its negative residuals and left it open.
- **WHAT `n_puddle` ACTUALLY MEANS** — carried from 2026-10-12, untouched.
  Three directions converge on it (§4.6.7's floor at 2.07× Chapter 3's, item
  18's literature use of the symbol for DOS broadening, and the fitted
  1.2e16 m⁻² at 2.4× this thesis's value). Cheapest route in is the two
  references the Southampton paper delegates `n₀` to, still not followed.
- **WHETHER THE FLOOR BELONGS INSIDE `C_q`'s ARGUMENT** — carried from
  2026-10-12, costed at +0.1527 % in `n`, cheap either way.
- **A THIRD KNOB: THE ROUNDING FORM** — created 2026-10-11, untouched for a
  second session.
- **BUILD ITEM 17 AS A SWEEP** — three requirements, three worked instances
  (today's §4.7 table is the third: a two-form comparison run at the model's
  own operating point with the null stated), still no sweep.
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, now **eight
  sessions overdue**, today carrying its eighth worked instance from today's
  own check (e).
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  untouched for a **seventh** session.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — untouched for the fifth session.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, untouched for a
  fourth session.
- **MECHANISE the magnitude item** — created 2026-10-03, eight sessions
  overdue.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, untouched for a fifth
  session.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** — untouched
  for a fifth session.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.
- **A second anchor for `Δ_c`, at any separation other than 3.3 Å** — open
  since 2026-09-21, still the top *physics* item in the thesis, and now
  **untouched for TWENTY-THREE consecutive sessions**. Today was a
  scope-and-plumbing session with no literature fetch at all, which is the
  §7.7 unequal-scrutiny finding in its purest form to date: a day spent
  entirely on what the code does and none on whether the physics is anchored.

**Commits this run:** 2 (the blast-radius module with its transcript and
figure; §4.6.12 with §§4.7 and 4.10 annotated in place). This
AUTOMATION_LOG.md entry makes **3**.

---

## 2026-10-14 — Three hand-found bindings were a rate; the class is a runtime object, and only 1 of 22 misses can fire

**Target:** §7.9's **second-ranked** item, "sweep every import-time binding of
every patched name", created 2026-10-13 on the stated ground that it is "what
makes item 16's flip verifiable rather than hopeful". It is now largely
closed, and it corrected a claim this session itself made four commits earlier.

`graphene_import_binding_sweep.py` (**6/6**, with
`import_binding_sweep_output.txt` and its figure) and
`graphene_binding_class.py`. Thesis §4.6.13 (new), §7.9 item 19 (new), and an
annotation to item 16.

### 1. The exposure a hand-maintained list leaves behind

An AST pass over all 53 modules: **149** module-level `from`-imports of
repository names across **73** distinct `(owner, name)` pairs, of which **25**
sites are left stale by the repository's **14** hand-listed rebinding sites.
Three had been found by hand (2026-10-03, 2026-10-10, 2026-10-13). **Twenty-two
had not.**

Of the ten distinct stale bindings, **seven are measured USED** — the holder's
own numbers move when only its local copy is perturbed:

| stale binding | holder's own numbers move by | |
|---|---|---|
| `cd.quantum_capacitance` | **9.8781** (988 %) | USED |
| `cd.n_puddle` | 35.764 % | USED |
| `crx.carrier_density` | 27.007 % | USED |
| `crx.W` | 27.007 % | USED |
| `pcm.mu` / `psc.mu` / `ptc.mu` | 9.1743 / 5.3067 / 3.9039 % | USED |
| `cd.mu`, `cd.C_ox`, `cd.W` | exactly 0 | never reach a number |

> **The usedness measurement needs no manager to be entered, which is what
> makes it applicable to the ten inline `try/finally` rebinding sites that
> cannot be entered at all.**

### 2. The repair: a computed class, not a longer list

A list is exactly what was wrong all three times, so the mechanism had to stop
being one. `from M import n` copies a **reference** at import time, so the set
to rebind is the set of names that at runtime hold the same **object** —
discoverable by identity, needing no import graph and no list.
`graphene_binding_class.rebind()` computes it at entry and is installed at all
three charge-path sites: `on_consistent_argument`, `on_exact`, and
`graphene_potential_declaration.main`, which had no second line at all.

**Behaviour unchanged, verified not asserted:** all three committed transcripts
regenerate **byte-identically**, and each was confirmed byte-reproducible
*before* the change, so the diff is a control and not a hope.

### 3. The class is necessary and NOT sufficient — and 2026-10-13 said this wrongly

`cd.quantum_capacitance` is deliberately **not** rebound. It is reached at a
contact where the series capacitance is `C_interface`, **65.93×** `C_ox`, and
the series factor runs **0.555–0.806** against **0.999** in the channel. The
object the owner holds after item 16's flip solves for `n` against `C_ox`, so
copying it to the contact installs a **channel solver at the contact** — a
different wrong answer, not a repair.

> **2026-10-13's precondition (i) — "`contact_edge_carrier_density()` must flip
> in the SAME COMMIT as `carrier_density()`" — was written as though one object
> served both. It does not. Item 16 needs TWO replacements, not one rebinding
> done in one commit.**

The exclusion lives once, in code, with its reason attached, and `rebind()`
**refuses an exclusion without a stated reason** — an unstated exclusion is how
all three 2026-10 faults were introduced.

### 4. And the static count over-reports, which retires this session's own claim

A miss needs the stale holder to be **imported in the patching process**.
Measured per site in a **fresh interpreter**, so the audit's own imports cannot
supply the answer: of the 22 remaining static misses the holder is **PRESENT
for 1** and **ABSENT for 21**. `import graphene_velocity_saturation_model`
brings in three repository modules and none of the four `mu` holders is among
them. `graphene_potential_declaration` does not import
`contact_resistance_crossover` at all.

> **So this morning's commit calling that site "the 2026-10-10 fault
> unrepaired in a third module" was WRONG ABOUT THIS REPOSITORY while right
> about the shape. There was no stale binding to leave behind and no wrong
> number.** The 25 is retained in place as an **upper bound**; the honest
> standing is three independent axes — **PRESENT, USED, REACHED** — and nothing
> currently in the repository clears all three. The sweep found a **landmine
> count, not a defect count.** The one survivor,
> `blast_radius.series_capacitance` leaving `cd.C_ox` stale, is also measured
> UNUSED, so it is cosmetic as well.

### 5. Four instrument faults, all four kept in place

**(a) A repair that blinds its own audit is not a repair** — §7.5.4's "created
by a fix" class in its purest form yet. Retrofitting `rebind()` made the three
sites it repaired **invisible** to the sweep, whose only notion of a rebinding
site was `alias.name = ...`, and validation (2) — which exists to keep
2026-10-13's finding on the record — **failed while nothing was wrong**. The
sweep now recognises both forms; the validation is rewritten to assert the same
fact in the form the repair left it in, and its original no-silent-miss clause
is kept rather than replaced.

**(b) The planted-binding control failed on the instrument.** The module
scanner skipped every file beginning with an underscore, so the sweep reported
**zero** bindings in a file it had itself just created with one. **Fourth**
occurrence here of a control failing on the instrument rather than on the thing
measured.

**(c) The identity control found a blind spot in the AST rule, and the rule was
extended rather than the control dropped.** **Twelve** second-order *renaming*
aliases (`dr.C_OX is gfet.C_ox`, `dr.N0_SHIPPED is gfet.n_puddle`,
`psc.W_CROSS_VACUUM is cd.W_GRAPHENE`) are invisible to any rule reading import
statements. What identity yields is not a directed edge but the **equivalence
class** of names sharing one object — and for "what must a patch rebind", the
class is the safe criterion and the import graph the unsafe one.

**(d) The arrival counter carries a positive control, so LATENT is a
measurement.** Both charge-path managers reach `cd` **zero** times under their
own workloads — agreeing with 2026-10-13 — but a broken counter reads zero too.
The same counter reads **exactly 7** when `cd` is exercised over its seven
metals.

### 6. Validations against exactly known values

Six, all passing: 0 misses for the 2026-10-10 pair, and both managers now cover
it by construction; 2026-10-13's finding surviving as a reasoned exclusion with
zero silent misses; a planted binding found **1 / 1 / 2** (AST rule, transitive
rule, identity — the import and the rename); identity-patch neutrality over
**98** wrapped sites at **0 ulp**; the usedness instrument's own null, a
perturbation by a factor of exactly 1.0, at **exactly 0**; and the arrival
counter's positive control at **7 of 7**.

### 7. Environment and automation health

Device reachable and folder connected at the **first** firing; the two later
firings will find both repositories already carrying today's entry and stop.
`pip install scipy` succeeded again (the device VM still lacks it). **No web
search was attempted** — today's work was entirely internal, which is the
§7.7 unequal-scrutiny finding for a second consecutive session. The whole
session ran **on the device VM**; nothing crossed the bridge. Every commit was
pushed and **verified against the GitHub API** before the next was written.
Nothing was deleted: the morning's 25 and its 27.0 % are both annotated in
place.

**Not yet covered (candidates for future runs):**

- **FLIP THE DEFAULT TO `B_T`** — item 16, still open, still **ranked first**,
  and its plumbing preconditions are now **discharged**. What is left is
  exactly two things, both named: (i) `graphene_contact_doping_model` needs its
  **own contact-correct replacement**, because today established that the
  contact cannot take the channel's object — the flip is two replacements, not
  one rebinding; (ii) §4.7's prose and §4.10's annotation need **rewriting**,
  since "only Pd gives a small, plausible positive residual" becomes false.
  Neither is a measurement; both are work.
- **CONVERT THE ELEVEN SURVIVING HAND-LISTED REBINDING SITES** — created today,
  **ranked second**, and deliberately not done: `mu`, `W`, `n_puddle` and
  `C_ox` sites in `graphene_velocity_saturation_model`,
  `graphene_diffusion_current_model`, `rf_small_signal_model`,
  `graphene_fmax_shortfall_decomposition`,
  `graphene_perfect_contact_counterfactual`,
  `graphene_fet_signed_carrier_model` and `graphene_double_regularisation`.
  Each is **latent by import today**, so converting one changes nothing until
  it is not — which is precisely the condition under which it will be
  forgotten. Cheapest route: convert them all mechanically and require the
  transcripts to be byte-identical, exactly as today's three were.
- **THE PRESENT/USED/REACHED TABLE AS A TEST, NOT A REPORT** — created today.
  The sweep is a module someone has to run. A test that fails when a *new*
  hand-listed rebinding site appears, or when a converted site regresses, is
  what would make today's repair durable. Related to item 17's unbuilt sweep.
- **A SECOND ANCHOR FOR `Δ_c`, AT ANY SEPARATION OTHER THAN 3.3 Å** — open
  since 2026-09-21, still the top *physics* item in the thesis, and now
  **untouched for TWENTY-FOUR consecutive sessions**. Two consecutive sessions
  have now run entirely inside the code with no literature fetch at all. This
  is §7.7's unequal-scrutiny finding not as a diagnosis but as a habit.
- **THE CONTACT MODEL'S OWN APPROXIMATIONS** — carried from 2026-10-13: a
  vacuum-like 3.5 Å gap with no interface dipole layer, and one
  `lambda_decay = 250 nm` for every metal, both plausibly larger than the
  12.5 % correction already found and neither ever costed.
- **WHAT `n_puddle` ACTUALLY MEANS** — carried from 2026-10-12. Three
  directions converge on it; today adds a fourth datum, that `cd.n_puddle` is
  measured **USED at 35.8 %** of the contact model's own numbers, so the
  parameter's meaning is load-bearing in two chapters and not one. Cheapest
  route in is still the two references the Southampton paper delegates `n₀` to.
- **WHETHER THE FLOOR BELONGS INSIDE `C_q`'s ARGUMENT** — carried, costed at
  +0.1527 % in `n`, cheap either way.
- **A THIRD KNOB: THE ROUNDING FORM** — created 2026-10-11, untouched for a
  third session.
- **BUILD ITEM 17 AS A SWEEP** — three requirements, three worked instances,
  still no sweep.
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, now **nine sessions
  overdue**.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  untouched for an **eighth** session.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — untouched for the sixth session.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, untouched for a
  fifth session.
- **MECHANISE the magnitude item** — created 2026-10-03, nine sessions overdue.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, untouched for a sixth
  session.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** — untouched
  for a sixth session.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.

**Commits this run:** 3 (the sweep; the binding-class repair with the
correction to the sweep's own headline; the thesis). This AUTOMATION_LOG.md
entry makes **4**.

---

## 2026-10-15 — The second anchor was an inequality all along, and identity cannot see an operating point

**Target:** the two items the 2026-10-14 entry ranked first and named as the
standing habit problem. The **top physics item** (a second anchor for `Δ_c`,
open since 2026-09-21, **twenty-four consecutive sessions untouched**, with two
consecutive sessions having run with no literature fetch at all) and **item
16's precondition (i)** (the contact's own contact-correct replacement). Both
moved; neither closed; and each turned out to be a different thing from what
its own description said.

`graphene_second_anchor_bound.py` (**6/6**), `graphene_contact_consistent_argument.py`
(**9/9**), a new `bc.rebind_by_operating_point()`, thesis §4.6.14 (new), and
annotations to §7.9 items 1, 16 and 19.

### 1. The second anchor exists, it is at 5.0 Å, and it is stated in words

The item asked for a second `(d, Δ_c)` pair. **There is none, and there was
never going to be one.** Khomyakov et al., *Phys. Rev. B* **79**, 195425
(2009), [arXiv:0902.1203](https://arxiv.org/abs/0902.1203) — the same paper
(N1)'s single anchor comes from — introduce the short-range term

```
Δ_c(d) = e^(−γd)(a₀ + a₁d + a₂d²)                               (their Eq. 4)
```

in Section III.3 **symbolically and never evaluate it numerically in the
text**: no `γ`, no `aᵢ`, no `Δ_c` tabulated at any separation. Section III.3
disposes of it with "|Δ_c(d)| becomes negligible for large d anyway", which is
why the paper does not need the numbers and the reader cannot have them.

> **Twenty-four sessions carried an item whose cheapest route was a search that
> cannot succeed. Recording that is the part of today's finding with the
> longest shelf life.**

What the paper *does* contain is a **second separation**. Its charge-transfer
calculation is run at **`d = 5.0 Å`** as well as at `d_eq ≈ 3.3 Å` (Fig. 4's
two point sets), and of the far one Section III.2 says the crossover
"decreases to its expected value, `W_M ∼ W_G`, for large separations". Through
`W_0(d) = W_G + Δ_c(d)`, that sentence **is** `Δ_c(5.0 Å) ≤ ε` for an `ε` the
paper never quantifies.

### 2. Which bounds the swept parameter, and excludes two thirds of the range

```
ℓ ≤ ℓ_max(ε) = (5.0 Å − 3.3 Å) / ln(0.9 eV / ε)
```

Reported as a **family** over `ε`, because the paper quantifies none:

| `ε` / eV | `ℓ_max` / Å | sweep excluded | `ELL_DEFAULT = 0.6 Å` survives? |
|---|---|---|---|
| 0.20 | 1.1303 | 30.8 % | yes |
| 0.10 | 0.7737 | 60.5 % | yes |
| **0.08** | **0.7024** | **66.5 %** | **yes** ← the paper's own scale |
| 0.05 | 0.5882 | 76.0 % | **NO** |
| 0.02 | 0.4466 | 87.8 % | **NO** |

`ε = 0.08 eV` is **borrowed from a sentence about a different quantity** (the
non-rigid-shift correction the paper then neglects) — it is the only number
the paper attaches to "small" nearby, and the table exists so that nothing
rests on the borrowing.

> **ANNOTATED 2026-10-18.** The row's marker "← the paper's own scale" is
> **wrong and is kept in place**, as is every number in the table. The paper
> quantifies **no** `ε` at 5.0 Å, so no row here is the paper's: the family is
> the result. And the caveat immediately above understates the borrowing — the
> `∼0.08 eV` is the deviation at `d ≈ 3.3 Å`, the separation about which the
> same sentence says Eq. (1) *holds*, and it measures the failure of the
> rigid-shift picture that `Δ_c` itself is. Wrong separation, wrong kind of
> quantity. Found while checking why the *transcript* of
> `graphene_second_anchor_bound.py` carried the attribution and none of this
> paragraph's qualification; see the 2026-10-18 entry. Read forward, the top of the committed sweep puts the
5.0 Å crossover **0.2898 eV above `W_G`**, which is not "≈ `W_G`" under any
reading of the source.

**The committed `ELL_DEFAULT = 0.6 Å` survives — by luck.** It was set as "the
midpoint of the swept range", a choice with no physics in it; had the midpoint
been 0.9 Å it would now be excluded. It is checked against literature for the
first time rather than merely stated.

**And no doping verdict moves.** Measured across the six metals with a
tabulated `d_eq`, the sign of `W_M − W_0` is identical at the bound and at the
old top of the range. **So this is a tightening of the error bar on `W_0`, not
a correction:** no committed number in Chapter 4 or Chapter 6 moves, Section
6.12's margin table is unchanged, and nothing was rewritten on its strength.

The item **narrows** rather than closes. Still open: a **two-sided**
constraint — nothing here bounds `ℓ` from below and the lower half of the
range survives untouched.

### 3. The contact gets a replacement, and identity turns out to be the wrong criterion for *to what*

Item 16's precondition (i) is **discharged**. The contact-correct numbers were
already measured (2026-10-13, −4.7696 % to −12.4964 % in `n`, against
−1.0212 % in the channel), and having them was not having an object, for three
reasons of which only the third was a surprise: the object belonged to the
audit; the harness rebound by hand, the mechanism 2026-10-14 retired; and —

> **`bc.rebind()`'s vocabulary is `{same replacement, left alone}`. It has no
> word for "rebind this member to a DIFFERENT object", which is exactly what a
> contact with an object of its own needs. IDENTITY ANSWERS WHICH NAMES MOVE
> TOGETHER AND DOES NOT ANSWER TO WHAT.**

`quantum_capacitance`'s class spans the channel (`C_ox`, series factor 0.999)
and the contact (`C_interface` = 65.93 × `C_ox`, series factor 0.555–0.806).
An operating point is not something identity can see, so the class must be
**partitioned** by one and the partition must come from the caller.
`bc.rebind_by_operating_point()` does that; `bc.rebind()` is **untouched** and
its three transcripts regenerate **byte-identically**. 2026-10-14's refusal
carries across in the form the same lesson takes here: there an exclusion
without a stated reason is refused, here an **unassigned class member** is,
because silently leaving a member at the shipped object is the identical fault
reached by omission from a map instead of from a list.

`on_consistent_argument_everywhere()` now installs four sites with an **empty**
left-alone set, asserted by the mechanism and backed by a negative control
that plants a missing member and requires the refusal to fire.

**2026-10-14's finding (ii) is thereby superseded in its mechanism and upheld
in its physics.** It said the contact must be *excluded* because the channel's
object would be wrong there. That is still true; what is no longer true is
that exclusion is the only thing the mechanism can do about it. The exclusion
in `bc.CONTACT_OPERATING_POINT` is retained in place, with its reason, rather
than deleted.

### 4. Validations against exactly known values

Fifteen across the two modules, all passing. The five that carry weight:

- **The new machinery fed the OLD argument reproduces the shipped
  `contact_edge_carrier_density()` bitwise over all seven metals (0 ulp).**
  This is what establishes that form A and form `B_T` at the contact differ in
  `C_q`'s **argument** and in nothing else. Without it every per-metal number
  could be the residue of some other difference between two implementations.
- **Today's owned object equals 2026-10-13's independent construction at 0 ulp
  over seven metals.** Two constructions of the same physics that disagree are
  one too many.
- **The quantum limit `n = 4(dV/A_F)²` approached at first order in `1/C_s`,
  error ratio exactly 100.000 when `C_s` is raised 100×.** At finite `C_s` the
  solve differs from the limit by `O(C_q/C_s)`, so a single threshold would
  have been a plausible-range check wearing a closed form's clothes; the exact
  statement is the **convergence order**.
- **A planted-`ℓ` round trip through the bound's inversion, 65 values, worst 2
  ulp**, and a **positive control on the exclusion counter** reading exactly
  0.0 at the top of the range and exactly 1.0 at the bottom — without which
  §2's "no verdict moves" and a broken counter would read identically.
- `ΔWF = 0` → `n = 0` **exactly** under both forms, and form `B_T` **odd in
  `dV` at 0 ulp**.

### 5. Two instrument faults, both kept in place

**(a) A check that read two knobs attributed the effect to the wrong one —
second consecutive session.** The oddness check was first written through the
shipped `contact_edge_carrier_density()` and failed at **3 ulp**. It then
failed at **2–3 ulp under form A as well**, which located the cause: the
shipped `dWF = work_function_metal - W_GRAPHENE` is **not antisymmetric in
double precision** — `(4.5+0.05)−4.5 = 0.050000000000000266` against
`(4.5−0.05)−4.5 = −0.04999999999999982`. The solver is exactly odd. A
pre-existing defect in the shipped contact model was briefly read as one
introduced today; it is cosmetic at 2–3 ulp against a 12 % effect, and it is
kept as its own validation row rather than dropped.

**(b) A closed-form check that was a tolerance check.** Draft (d) asserted
`n = 4(dV/A_F)²` to `1e-9` at one large `C_s` and failed at `2.3e-8` — because
the deviation is a real `O(C_q/C_s)` physical remainder, not numerical error.
The threshold was the defect, not the solver. Replaced by the convergence-order
form above.

### 6. Environment and automation health

Device reachable and folder connected at the **first** firing of the day
(08:17 UTC); the later firings will find both repositories carrying today's
entry and stop. `scipy` was already importable in the device VM this session
(the first session not to need `pip install scipy`). **Live web search WAS used
— the first literature fetch in three sessions**, and it is the reason §1
exists; §7.7's unequal-scrutiny habit is interrupted rather than cured.
`arXiv:0902.1203` was read through the ar5iv HTML rendering, which
**truncated three times** at the opening of Section III.3 — including at
explicit character offsets of 38 000 and 52 000 — so Table 1's remainder and
Fig. 4's point values were **never retrieved**, and §1's claim that the fitted
coefficients are absent is a claim about the **retrieved portion**. No
PubMed/PMC fetch was attempted (three logged reCAPTCHA occurrences). The whole
session ran **on the device VM** except the literature fetch; nothing else
crossed the bridge. Nothing was deleted.

**Not yet covered (candidates for future runs):**

- **FLIP THE DEFAULT TO `B_T`** — item 16, still open, still **ranked first**,
  and now blocked on **exactly one thing, which is writing**: precondition
  (ii), §4.7's prose and §4.10's annotation. "Only Pd gives a small,
  plausible positive residual" cannot survive a −12.4964 % correction to Pd's
  `n_contact`. Both preconditions for the *plumbing* are discharged and
  `on_consistent_argument_everywhere()` exists; what is left is rewriting two
  passages and then flipping. **Cheapest route: do precondition (ii) and the
  flip in the same session, since the prose is only rewritable once the new
  numbers are the committed ones.**
- **A LOWER BOUND ON `ℓ` — THE TWO-SIDED HALF OF THE SECOND ANCHOR** —
  created today, **ranked second**, and the top *physics* item now that the
  upper bound exists. Cheapest route named: **digitise Fig. 4's upper (5.0 Å)
  point set** in `arXiv:0902.1203`, which turns the verbal "≈ `W_G`" into a
  number with an error bar and bounds `ℓ` from both sides. The fetch could not
  retrieve the figure today; a PDF-page read of that figure is the obvious
  route and is one of the few genuine reasons to cross the bridge to the cloud
  sandbox.
- **A `Δ_c` FORM FOR THE CHEMISORBED METALS** — sharpened today. (N1) is
  outside its regime for Ti, Ni and Pd (`d_eq` = 2.05–2.30 Å < `d_0`)
  whatever `ℓ` is, and today's bound does not rescue them: `W_0(Ni) = 62.5 eV`
  at `ℓ = 0.3 Å` is (N1) diverging, not a prediction. These three metals are
  §7.9 item 2 and contact A of *both* Chapter 6 headline pairs, so the
  divergence is load-bearing.
- **CONVERT THE ELEVEN SURVIVING HAND-LISTED REBINDING SITES** — carried from
  2026-10-14, untouched. `mu`, `W`, `n_puddle`, `C_ox` sites in
  `graphene_velocity_saturation_model`, `graphene_diffusion_current_model`,
  `rf_small_signal_model`, `graphene_fmax_shortfall_decomposition`,
  `graphene_perfect_contact_counterfactual`,
  `graphene_fet_signed_carrier_model`, `graphene_double_regularisation`. Each
  is latent by import, so converting one changes nothing until it is not —
  precisely the condition under which it is forgotten. **Today adds a reason
  to do them with the new mechanism rather than the old: `mu` and `W` are read
  at both the channel and the contact, so at least two of the eleven are
  operating-point-split exactly as `quantum_capacitance` was, and converting
  them with `bc.rebind()` would install one object across two operating
  points — the fault §4.6.14 was written to prevent.**
- **THE PRESENT/USED/REACHED TABLE AS A TEST, NOT A REPORT** — carried from
  2026-10-14. Now wants a fourth column: whether a class member's operating
  point is assigned. Related to item 17's unbuilt sweep.
- **THE 2–3 ULP ANTISYMMETRY IN THE SHIPPED `dWF` SUBTRACTION** — created
  today, cosmetic, and deliberately ranked low. Worth a line in item 15's
  tolerance audit rather than a session.
- **THE CONTACT MODEL'S OWN APPROXIMATIONS** — carried from 2026-10-13: a
  vacuum-like 3.5 Å gap with no interface dipole layer, and one
  `lambda_decay = 250 nm` for every metal, both plausibly larger than the
  12.5 % correction already found and neither ever costed. **Today raises its
  rank: the contact's series factor is 0.555–0.806, so the contact is where
  the model is least forgiving and also where it is least examined.**
- **WHAT `n_puddle` ACTUALLY MEANS** — carried from 2026-10-12; four
  directions now converge on it. Cheapest route in is still the two
  references the Southampton paper delegates `n₀` to.
- **WHETHER THE FLOOR BELONGS INSIDE `C_q`'s ARGUMENT** — carried, costed at
  +0.1527 % in `n`, cheap either way.
- **A THIRD KNOB: THE ROUNDING FORM** — created 2026-10-11, untouched for a
  fourth session.
- **BUILD ITEM 17 AS A SWEEP** — three requirements, three worked instances,
  still no sweep.
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, now **ten sessions
  overdue**, and today produced a second instance of the exact fault it
  targets (§5(b): a threshold that was the defect).
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  untouched for a **ninth** session.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — untouched for the seventh session.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, untouched for a
  sixth session.
- **MECHANISE the magnitude item** — created 2026-10-03, ten sessions overdue.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, untouched for a seventh
  session.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** —
  untouched for a seventh session.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.

**Commits this run:** 3 (the second anchor as a bound; the contact's own
replacement and the operating-point partition; the thesis). This
AUTOMATION_LOG.md entry makes **4**.

## 2026-10-16 — Three censuses of the same fault disagreed for twelve days, and the pattern turned out to be two classes with opposite behaviour

**Target:** the 2026-10-15 entry ranked item 16's flip first and named its
cheapest route (precondition (ii) plus the flip, in one session). That route
was taken up and **abandoned within the first twenty minutes, for a reason
that became the session**: establishing the flip's blast radius means
regenerating committed transcripts and reading what moved, so the first step
was to run the instrument of record —
`graphene_pristine_transcript_audit.py`, built 2026-10-02 on the principle
that *a committed transcript is a claim about the code beside it*. It had not
been run since **2026-10-04**. Asked today it reported **four of its own
subjects STALE or ERROR**.

`graphene_frozen_default_semantics.py` (new, **24/24**), amendments to
`graphene_mutation_arrival_probe.py` and
`graphene_pristine_transcript_audit.py`, four transcripts regenerated,
thesis §4.6.15 (new), §7.9 item 20 (new) and a precondition (iii) on item 16.

### 1. Four instruments, one date, and a fault class that regrew behind them

Since 2026-09-28 this repository has tracked the **inert knob**: a
module-level constant captured as a function parameter default, bound once
when the `def` executes, so a rebinding — the one mechanism by which every
mutation audit and every patch manager here delivers a change — never
arrives. The audit then measures nothing and reports a pass. Three separate
instruments census it. Their committed transcripts read **0**, **0** (an
empty census behind a `[PASS]` on "no module-level constant is captured as a
default anywhere") and **8**. All three were last written **2026-10-04**.
Run today: **21**, **20**, **16**.

| instrument | committed (2026-10-04) | today |
|---|---|---|
| `figure_provenance_audit` Section C | 0 across 0 of 38 files | **21 across 5 of 56** |
| `mutation_arrival_probe` Section 7b | empty census, 34 passed / 0 failed | **10 site-groups, 8 MIXED**, 32 / 2 |
| `cross_module_delivery_audit` | 8 captures | **16** |
| `pristine_transcript_audit` | every subject IDENTICAL | **4 STALE or ERROR** |

> **A fault class driven to zero on 2026-10-04 regrew over twelve days in
> front of three censuses, none of which was being read — and the fourth
> instrument, whose only purpose is to notice that, was last run on the same
> day. This is §7.7's unequal-scrutiny habit one level out: the repository
> has spent three weeks learning that a guard which is not run is not a
> guard, and had not applied the rule to the guard that checks the guards.**

`git status` cannot see any of it. The transcripts are unmodified; it is the
code that moved out from under them — the exact fault the pristine audit was
built for, which is also why it is the instrument that found it.

### 2. The disagreement was the finding: 21 and 20 are counts of different things

Over one path set both censuses name the **same five files** and **four of
the five agree exactly**. They differ by one occurrence, at
`graphene_velocity_saturation_model.py:589` —
`_R_harmonic(V_g, V_hi, n_quad=N_QUAD)`, defined **inside another function's
body**. Measured in all three configurations against exactly known values,
with `K` rebound 1.0 → 99.0:

| where the `def` is | `def` executes | returns after rebind | class |
|---|---|---|---|
| module body | once, at import | exactly **1.0** | FROZEN |
| **class body** | once, at import | exactly **1.0** | **FROZEN** |
| function body | **every call** of the enclosing function | exactly **99.0** | **LIVE** |

> **THE SYNTACTIC PATTERN `def f(x=K)` IS TWO CLASSES WITH OPPOSITE BEHAVIOUR
> UNDER THE ONE OPERATION THIS REPOSITORY USES TO DELIVER EVERY CHANGE. A
> FUNCTION-NESTED CAPTURE IS NOT AN INERT KNOB. IT IS LIVE.**

So neither census is simply right. **figure_provenance's 21** counts the
syntax correctly and over-counts the *fault* by exactly the live captures —
and its own historical figures ("49 occurrences across 14 files, 14 of them
inside the audit modules … 37 remain, and the 13 in audit modules are the
next run's top item") are syntactic counts of a class that is not uniformly
the fault, never split. **mutation_arrival_probe's 20** counts the fault
correctly *for the wrong reason*: `classify()` confirmed each static site
with `getattr(module, fname)` — added correctly, to drop a function a later
redefinition or decorator had replaced — and a function-nested function is
not a module attribute, so the site was dropped **silently**: no row, no
class, no reason.

**And the accident does not generalise.** The same line also drops a
**class-body** capture, which is *frozen* — so the next capture of that shape
would be **undercounted**. That is the dangerous direction: the signature of
this fault class is that a mutation silently fails to arrive and the audit
reports a pass, and an instrument that silently omits a true site has the
identical signature one level up. The repository has no class-body capture
today; **that is luck, not a property of the mechanism.**

### 3. Reconciled to zero residue

```
21 syntactic  =  20 FROZEN (the genuine inert knobs)  +  1 LIVE
```

Four checks, because a reconciliation with a remainder means a third cause is
present and unaccounted for: the arbiter reproduces figure_provenance's total
exactly (R1), mutation_arrival_probe's confirmed total exactly (R2), accounts
for the whole gap with nothing left over (R3), and agrees **per `(file,
constant)`** and not only in total (R4). R4 was added deliberately — two
totals can agree while the rows behind them disagree in compensating
directions — and the mutation tests below confirm it is the *only* check that
catches a single dropped row.

Per file: velocity_saturation 9 = 8 + 1, diffusion_current 5 = 5 + 0,
exact_charge_rederivation 3 = 3 + 0, perfect_contact_counterfactual 2 = 2 +
0, second_anchor_bound 2 = 2 + 0 — **the last of these added by yesterday's
session.**

### 4. The repairs

**(a) `classify()` takes the verdict from the SCOPE instead of discarding
what it cannot resolve.** `_frozen_sites_static` now descends explicitly and
records `in_function_body` (fifth tuple field, appended not inserted, so both
existing readers are unaffected). Function-enclosed → `live_by_scope`;
otherwise → `frozen_unresolved`, which **still counts as frozen** because the
static evidence stands and `getattr` merely cannot reach it. The fix is a
verdict from the scope and **not a widening of the resolver** — resolvability
was never the question being asked. **COUNT-NEUTRAL:** the frozen total is 20
before and after, and the suite's output differs in exactly the two lines
that now carry `1 rescoped-live`. The `getattr` confirmation's original
purpose is kept, not removed (V6: a module attribute rebound after import
genuinely no longer carries the captured default).

**(b) The pristine audit could not tell a crashed suite from a correctly
failing one — and this session's own decision made that load-bearing.** Every
suite here ends `return 0 if failed == 0 else 1`, so a red guard is a nonzero
exit *by design*. `run_suite` discarded stdout on any nonzero exit and
reported ERROR **without comparing the transcript at all**. In the real
twelve-day window it did fire, but blind: `exit 1: ` with empty stderr, no
diagnosis, both transcripts unexamined. Harmless while every suite was green;
not harmless from today, when this repository began committing **red
transcripts in preference to stale green ones** — two suites would then have
had transcripts nothing checks, permanently. ERROR is now reserved for a run
producing no usable stdout. **New check:** a transcript's own `TOTAL: N
passed, M failed` line and its suite's exit code are two claims about one
run, and must agree.

**(c) Four transcripts regenerated, and two of them are now red.**
`mutation_arrival_probe` 32/2 and `dead_name_sweep` 6/2, against committed
34/0 and 8/0. **Committing red is the point, not a concession:** the
alternative was keeping a green transcript from 2026-10-04 over a guard that
has been failing since, which is the shape of fault this repository has spent
three weeks learning to refuse. Baseline after: **13 IDENTICAL, 0 STALE, 0
NUMERIC_DRIFT, 0 ERROR**, two suites `exit 1 (guard red, transcript still
checked)`.

### 5. A defect in this session's own first draft, kept on the record

`graphene_frozen_default_semantics.py`'s first draft classified by **nesting
depth** — depth 0 frozen, depth ≥ 1 live — and its own control planted a
method in a class body, called it depth 1, and asserted *in prose* that it
was therefore not an inert knob. A class body executes once at import, so
that is wrong (measured, V3). The repository contains no class-body capture,
**so the draft's count would have been right by luck — the same luck
`classify()` was running on**, which is what made it worth recording rather
than quietly fixing. Controls C3/C4/C5 now pin all three cases including a
class defined *inside* a function, where a depth rule and a class rule both
get it wrong. The draft was caught by reading its own control's prose, not by
a check: no check could fail on it, because the repository has no instance.

### 6. Mutation tests, run outside the repository

So the instrument and the measurement do not move in one step (the
2026-09-28 rule). All five mutants detected:

| mutant | detected by |
|---|---|
| the depth rule restored | C1, C2, C3, C4 |
| every capture called frozen | C1, C2, C5, R2, R3, R4 |
| one frozen row dropped | **only R2 and R4** — which is why R4 exists |
| `run_suite` reverted to discarding stdout | C7, and the two red suites revert to ERROR |
| 2026-10-04's green TOTAL restored over today's red guard | the new consistency check, **by name**: "transcript claims 0 failed, suite exits 1" |

### 7. Environment and automation health

Device reachable and folder connected at the **first** firing of the day
(03:4x UTC); the later firings will find both repositories carrying today's
entry and stop. `scipy` was **not** importable in the device VM and was
`pip install`ed (1.15.3) — the 2026-10-15 entry's note that it was already
present did not hold today. **No live web search was used and none was
needed**: this session touched no physics and made no claim about the
literature, which is the first time that has been true without also being a
failure to look. §7.7's unequal-scrutiny habit is therefore neither
interrupted nor cured today. The whole session ran **on the device VM**;
nothing crossed the bridge. Nothing was deleted.

**Not yet covered (candidates for future runs):**

- **CONVERT THE 20 INERT KNOBS, AND DISPOSITION THE FOUR DEAD NAMES** —
  created today as §7.9 item 20 and **ranked first**, ahead of item 16,
  which carries much the larger physics cost. The ordering is a
  **measurement constraint, not a change of priority**: the baseline is
  clean of staleness but now carries two deliberately red suites, and a flip
  run over it would mix its own effect with theirs exactly as it would have
  mixed with the drift. The work is mechanical (`def f(x=CONST)` →
  `None`-sentinel) and is now **specified correctly for the first time**:
  **20, not 21**, with the live capture named so it is not converted by
  mistake. Required of each conversion: the affected transcript regenerates
  **byte-identically**, since the change is behavioural only under
  rebinding. The four names are
  `graphene_potential_census_audit.CHANNEL_POTENTIAL_READINGS` (genuinely
  DEAD) and three **hidden** dead names — CROSS_MODULE resolutions resolved
  on *spelling* whose claimed readers do not import the module at all:
  `graphene_contact_consistent_argument.E_CHARGE` (**added yesterday**),
  `graphene_diffusion_current_model.HBAR`,
  `graphene_quantum_capacitance_argument.A_OX`, over 12 unsupported rows.
  The 2026-10-01 precedent applies: give a reader, or declare exempt with a
  stated reason.
- **FLIP THE DEFAULT TO `B_T`** — item 16, still open, now **ranked second**
  for the reason above and for no other. Preconditions (i) and (iii) are
  discharged; (ii) — §4.7's prose and §4.10's annotation — is untouched and
  is still writing, not plumbing. "Only Pd gives a small, plausible positive
  residual" cannot survive a −12.4964 % correction to Pd's `n_contact`.
  **Cheapest route unchanged: precondition (ii) and the flip in the same
  session, after item 20 clears the baseline.**
- **A LOWER BOUND ON `ℓ` — THE TWO-SIDED HALF OF THE SECOND ANCHOR** —
  carried from 2026-10-15, untouched, and still the top *physics* item.
  Cheapest route named there: **digitise Fig. 4's upper (5.0 Å) point set**
  in `arXiv:0902.1203`. The ar5iv HTML truncated three times yesterday, so a
  **PDF-page read of that figure** remains one of the few genuine reasons to
  cross the bridge to the cloud sandbox.
- **A `Δ_c` FORM FOR THE CHEMISORBED METALS** — carried, untouched. (N1) is
  outside its regime for Ti, Ni and Pd (`d_eq` = 2.05–2.30 Å < `d_0`)
  whatever `ℓ` is; `W_0(Ni) = 62.5 eV` at `ℓ = 0.3 Å` is (N1) diverging, not
  a prediction. §7.9 item 2, and contact A of *both* Chapter 6 headline
  pairs, so the divergence is load-bearing.
- **CONVERT THE ELEVEN SURVIVING HAND-LISTED REBINDING SITES** — carried
  from 2026-10-14, untouched for a second session. **Today sharpens the
  prerequisite rather than the task**: these sites are latent *by import*,
  and today showed that a capture's liveness also depends on the **scope of
  its `def`**, so a site's conversion has two independent ways of arriving
  or not. Any sweep over the eleven should report both, and item 20 should
  be done first so the two are not confounded.
- **THE PRESENT/USED/REACHED TABLE AS A TEST, NOT A REPORT** — carried from
  2026-10-14. Wanted a fourth column (is a class member's operating point
  assigned); today adds a fifth (is a capture frozen or live by scope).
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, now **eleven
  sessions overdue**. Today produced no new instance, which is the first
  session in three that did not.
- **THE 2–3 ULP ANTISYMMETRY IN THE SHIPPED `dWF` SUBTRACTION** — carried
  from 2026-10-15, cosmetic, deliberately ranked low; a line in item 15's
  audit rather than a session.
- **THE CONTACT MODEL'S OWN APPROXIMATIONS** — carried from 2026-10-13: a
  vacuum-like 3.5 Å gap with no interface dipole layer, and one
  `lambda_decay = 250 nm` for every metal, both plausibly larger than the
  12.5 % correction already found and neither ever costed.
- **WHAT `n_puddle` ACTUALLY MEANS** — carried from 2026-10-12; four
  directions converge on it. Cheapest route in is still the two references
  the Southampton paper delegates `n₀` to.
- **WHETHER THE FLOOR BELONGS INSIDE `C_q`'s ARGUMENT** — carried, costed at
  +0.1527 % in `n`, cheap either way.
- **A THIRD KNOB: THE ROUNDING FORM** — created 2026-10-11, untouched for a
  fifth session.
- **BUILD ITEM 17 AS A SWEEP** — three requirements, three worked instances,
  still no sweep.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  untouched for a **tenth** session.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — untouched for the eighth session.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, untouched for a
  seventh session.
- **MECHANISE the magnitude item** — created 2026-10-03, eleven sessions
  overdue.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, untouched for an eighth
  session.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** —
  untouched for an eighth session.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.

**Commits this run:** 6 (the semantics module and the reconciliation;
`classify()`'s scope verdict and the first red transcript; the three
remaining stale transcripts; the pristine audit's red/crashed conflation; its
regenerated transcript; the thesis). This AUTOMATION_LOG.md entry makes **7**.

---

## 2026-10-17 — The pattern was also two spellings, the census that could see the other one has been reporting it as "dormant" for sessions, and this session claimed its finding was invisible

**Status:** Item 20 cleared, and the half of the fault class it did not cover
found, converted, and then found to have been on the record all along.

**Entry point.** Yesterday ranked item 20 first for a measurement reason — the
baseline carried two deliberately red suites, so a flip run over it would mix
its own effect with theirs. That ordering held: today is the conversion, and
the baseline it leaves is clean.

### 1. Item 20, cleared: the 20 inert knobs converted

`def f(p=CONST)` → `p=None` plus a body read, at exactly the 20 FROZEN sites,
leaving alone the one LIVE capture (`main._R_harmonic`,
graphene_velocity_saturation_model line 589) that yesterday's module
separated out. Per file: diffusion_current 5, exact_charge_rederivation 3,
perfect_contact_counterfactual 2, second_anchor_bound 2,
velocity_saturation 8 of 9.

**An independent AST pass, written for this session, reproduced the arbiter's
split before any edit** — 21 total, 20 frozen, the LIVE one at line 589 — so
the site list the conversion ran on is not the arbiter's output taken on
trust. Yesterday's requirement was that each affected transcript regenerate
**byte-identically**, and all five do; the four fast ones were also confirmed
to reproduce their *committed* transcripts before the edit, so "identical"
means identical to the claim on disk.

### 2. The same fault in a second spelling, which no census this repository
reconciles was counting

All three reconciled censuses resolve a default by looking for an `ast.Name`
whose identifier is a module-level **`ast.Assign`** target in the same file.
`import x as y` is an `ast.Import`, so `y` is never in that set, and

    import other as o
    def f(x=o.K): ...

is outside all three — while being frozen by exactly the same mechanism.
Thirteen such sites, in three files: graphene_sensitivity_audit 8
(`icm.rho_bulk_uohm_cm`, `p_default`, …), graphene_diffusion_current_model 3
(`vsat.OMEGA_OP`, `rf.W_RF`, `rf.N_FINGERS_RF`), graphene_velocity_saturation
2 (`rf.W_RF`, `rf.N_FINGERS_RF`). **So the inert-knob total for 2026-10-16 was
20 + 13 = 33, not 20**, and the "census reads zero" of 2026-10-01 and
2026-10-04 — and of this morning's own commit — was a statement about one of
two spellings.

New suite `graphene_attribute_default_capture.py`, **committed RED (25/2)
with its transcript before the conversion**, so the before/after is two
transcripts rather than one claim that a repair worked; green at **34/0**
after. That was only possible because 2026-10-16 stopped `run_suite`
reporting a nonzero exit as ERROR with its transcript uncompared — the
amendment earning its keep one day later.

**The two worst sites are the two whose docstrings deny them.**
`graphene_velocity_saturation_model._at_literature_geometry` carried
`W=rf.W_RF, N_fingers=rf.N_FINGERS_RF` under *"Late binding throughout — **no
module-level name is captured as a default** (the 2026-10-03 frozen-default
bug, and the 2026-10-04 repeat of it)"*. **True as written and false in
substance**: no module-level *name* was captured, because what was captured
was an *attribute*. The claim was audited by spelling, passed, and the thing
it was auditing for was in the signature it annotates. Its twin in
graphene_diffusion_current_model says *"Late binding only — the 2026-10-03
frozen-default fault"* over the same two defaults. Both sentences are **kept
verbatim** inside the corrected docstrings with the date and the reason.
Section 4 measures the repaired function taking `W = 7e-06` while a
reproduction of the old signature keeps `4e-05`, exactly.

### 3. THE RETRACTION — and it is this session's own claim being retracted

The first two commits today state that **no census here has ever been able to
see the attribute spelling. That is false.**
`graphene_cross_module_delivery_audit`'s Section 6 has censused exactly this
pattern since it was written — **16 rows in its committed transcript** — its
`cross_module_captures()` docstring states the `classify()` frame problem the
new suite's Section 5 re-derives from scratch, and its own prose already
reads: *"Dormant is not 'safe': it means no call site rebinds the name TODAY.
Eight of these are graphene_sensitivity_audit.py, which perturbs by passing
keywords explicitly; a single future setattr converts all of them at once,
which is why they are enumerated rather than dismissed."* Three of Section
5's conclusions were reached independently and claimed as new. It was
falsified by running the suites this session's own conversion invalidated —
i.e. by the repository's own standing instruments, not by a review.

**What survives the retraction, which is not nothing:**

1. **The repair.** Those rows were "enumerated rather than dismissed" and then
   left in place across sessions. 13 are converted today. A census that
   enumerates a fault indefinitely without it ever being repaired is
   2026-10-02's fault — a repair built and not connected — in slow motion, and
   it is the one failure mode this repository has not yet named for itself.
2. **The two false docstrings.** cross_module_delivery_audit lists
   `_at_literature_geometry` twice as `dormant` and does not notice that both
   functions' docstrings explicitly deny being captures.
3. **The reconciliation nobody did.** 2026-10-16 named three instruments as
   censusing "the pattern", reconciled 21 against 20 *to the occurrence*, and
   recorded the third as agreeing at a **different number** — 8 committed, 16
   live — without reconciling it. 8 and 16 are not reconcilable against 21 and
   20 **because they count a different spelling**. The session that split the
   pattern into two *classes* (frozen/live by scope) did not notice it was
   also two *spellings*, and the 2026-10-04 "agreement" of 0, 0 and 8 was
   agreement read **across a partition**. New Section 3b asserts it, and T4
   reconciles the committed 16 exactly: **13 constants + 3 `2.0 / np.pi` rows,
   no remainder** — the `np.pi` three left alone and named, since nothing
   rebinds `np.pi`.
4. **The mechanism.** cross_module_delivery_audit censuses the spelling; it
   does not establish its semantics. Section 1 does, against exactly known
   values: module-body attribute capture FROZEN, class-body FROZEN,
   function-body LIVE, **body read LIVE** — so the fault is the default and
   not cross-module access, which is what makes the sentinel the right fix —
   and `setattr` and `mod.attr =` the same mutation. **A7 is the decisive
   one**: the shipped Name census is handed a file whose only content is this
   pattern and required to report 0. It does.

### 4. Three defects in this session's own work, all caught by checks

- **R1's first draft** counted the suite's own Section-4 scaffolding in the
  repository total and **FAILED**. R1 now excludes SELF.
- **S1's first draft** asserted zero captures in the suite file and **FAILED**
  on the same two rows.
- **Section 4's reproduction, written as a nested `def` in the suite file, was
  a real cross-module captured default in the tree.** Because Section 4
  rebinds `rf.W_RF`, cross_module_delivery_audit's census classified it **LIVE
  FAULT** and turned that suite **red, 25/0 → 23/1, on this module's own
  scaffolding.** It is a temp plant now, and S1b is strengthened from "no
  FROZEN capture" to "no capture of either class": **passing one's own
  self-census is not the same as not planting the pattern in the tree other
  audits read.** That is a new rule, and it generalises past this file.
- A first Section 5 used a **textual grep** for the rebinding route and
  reported "no" for all thirteen; the route is `setattr(mod, name, v)` over
  *computed* mod and name in seven modules here, so the grep answered a weaker
  question. Replaced by putting the question to the repository's own oracle,
  on a plant.

### 5. The suite that would have gone vacuous, and the baseline that fixes it

The first draft had, in a sharper form, the exact defect it was written
about: **once the conversion landed, R1/R2 and the whole of Section 5 passed
vacuously over an empty row set, and the SUMMARY's "20 + 13 = 33" recomputed
itself to "20 + 0 = 20"** — a transcript that erases its own finding the
moment the finding is acted on, which is precisely how the Name census came
to be read as a clean bill of health for five sessions. Three fixes:
`BASELINE_2026_10_17` records the 13 by (file, qualname, param, default
expression) with **no line numbers** (a baseline that goes stale on an
unrelated edit teaches the operator to ignore it); **R4** requires every one
to be gone — *a count can be satisfied by the wrong 13* — and **R5** is the
standing half, no attribute capture may appear that the baseline did not
name; Section 5 measures the oracle's frame on a **plant**, where it cannot
evaporate; and R1/R2 say in their own detail lines that they read vacuously
today.

### 6. Six transcripts regenerated, every one of which changes

| transcript | before | after |
|---|---|---|
| mutation_arrival_probe | 32 / 2, "10 site-groups remaining", 8 MIXED | **33 / 1**, 1 site-group, 0 MIXED |
| cross_module_delivery_audit | 25 / 0, census 16 | 25 / 0, census **3** |
| figure_provenance_audit | 21 occurrences in 5 files | **1 in 1** |
| frozen_default_semantics | 21 syntactic / 20 frozen | **1 / 0** |
| dead_name_sweep | 6 / 2 | 6 / 2, rows changed |
| log_sensitivity_step_audit | STALE at HEAD | regenerated |

`graphene_frozen_default_semantics.py`'s SUMMARY is **amended rather than left
to print a stale sentence**: it used to end *"NOT done here: converting the N
genuine inert knobs, ranked first in AUTOMATION_LOG for the next run"*, which
at N = 0 reads as an instruction to convert nothing, and its *"all three
censuses read zero"* line is now followed by **READ THIS NUMBER NARROWLY** and
the cross-reference to the other spelling. That sentence is the exact form of
words today's finding is about, sitting in the transcript of the module that
established the count.

**`graphene_pristine_transcript_audit` was run, twice.** The first run is what
found five STALE transcripts; after the regenerations it is **13 passed, 0
failed — 15 suites IDENTICAL, 0 STALE, 0 NUMERIC_DRIFT, 0 ERROR**, the first
fully clean baseline since 2026-10-04, with the new suite registered in
`SUITES` in the same commit that created it (the 2026-10-03 wiring rule) and
IDENTICAL on its first registered run.

### 7. Environment and automation health

Device reachable and folder connected at the **second** firing (10:3x UTC);
04:30 found it asleep. **Two environment facts worth carrying forward.**
First, `scipy` was again absent in the device VM and was `pip install`ed
(1.15.3) — second session in a row, so the 2026-10-15 note that it is present
should be treated as false. Second, **background processes do not survive a
`device_bash` call**: `setsid nohup … &` was tried three times and every
child was reaped at the call boundary, so the 180 s per-call cap is a **hard**
cap on any single command. `graphene_velocity_saturation_model.py` takes 5
minutes, so its pristine-vs-converted pair was run in the **cloud sandbox**,
with the converted sources there confirmed by **sha256** to be the same bytes
as the ones committed; its pristine run also reproduced the committed
transcript exactly, on a different scipy (1.18.1). Nothing else crossed the
bridge and nothing was deleted. **No live web search was used and none was
needed:** this session touched no physics and made no claim about the
literature — the second such session, and §7.7's unequal-scrutiny habit is
therefore untouched rather than cured. One artifact recorded: the clone was
`--depth 30`, and mutation_arrival_probe's T_HOP measurement reads a pinned
blob with `git cat-file`, which **failed with CalledProcessError(128) in a
shallow clone and silently lost four checks behind the failure**. The clone
was unshallowed before anything was committed.

**Not yet covered (candidates for future runs):**

- **A CENSUS THAT ENUMERATES A FAULT FOREVER IS NOT A GUARD — AUDIT EVERY
  "DORMANT", "ENUMERATED" AND "REPORTED, NOT ASSERTED" ROW IN THE
  REPOSITORY** — created today and **ranked first**, because it is the
  generalisation of today's retraction and the one failure mode this
  repository has not named for itself. cross_module_delivery_audit's 16 rows
  sat behind a `[PASS]` with the correct diagnosis written beside them, and
  nothing converted them until a session arrived at the same finding from the
  other direction and claimed it was new. The work: sweep every suite for a
  row that is *reported* rather than *checked*, and for each either turn it
  into a check with a baseline (today's R4/R5 shape) or record a stated reason
  it must stay a report. Cheapest entry: `grep` for `dormant`, `reported, not
  asserted`, `enumerated`, `NOTE:` across the suites — the
  cross_module_delivery_audit NOTE about 20 computed-name `setattr` calls
  ("the census cannot be complete by construction") is a second instance
  already visible.
- **A MODULE MUST NOT PLANT THE PATTERN IT CENSUSES IN THE TREE OTHER AUDITS
  READ** — created today from §4's third defect, and the rule is wider than
  the file that generated it. Every suite here plants controls; the check
  should be mechanical, not a habit. Candidate: a guard that runs each
  suite's *own* arbiter plus every *other* census over the suite files.
- **FLIP THE DEFAULT TO `B_T`** — item 16, now **ranked second**. Yesterday's
  reason for ranking it behind item 20 is discharged: item 20 is done and the
  baseline is clean of staleness *and* of deliberate red, so a flip run's
  effect is now separable from everything else. Preconditions (i) and (iii)
  were already discharged; **(ii) — §4.7's prose and §4.10's annotation — is
  still untouched and is still writing, not plumbing.** "Only Pd gives a
  small, plausible positive residual" cannot survive a −12.4964 % correction
  to Pd's `n_contact`. **Cheapest route unchanged: precondition (ii) and the
  flip in the same session**, and there is now nothing ahead of it on
  measurement grounds.
- **A LOWER BOUND ON `ℓ` — THE TWO-SIDED HALF OF THE SECOND ANCHOR** — carried
  from 2026-10-15, untouched for a third session, and still the top *physics*
  item. Cheapest route named there: digitise Fig. 4's upper (5.0 Å) point set
  in `arXiv:0902.1203`. The ar5iv HTML truncated three times, so a **PDF-page
  read of that figure** remains one of the few genuine reasons to cross the
  bridge to the cloud sandbox — and today establishes a second: any suite over
  ~170 s.
- **CHAPTERS 2 AND 3 ARE STILL UNDRAFTED** — their computational results have
  been complete for weeks. This is the largest remaining block of pure
  writing and has now gone three sessions without being touched, while four
  consecutive sessions have gone to instrument work. Worth ranking against
  item 16 deliberately rather than by default.
- **A `Δ_c` FORM FOR THE CHEMISORBED METALS** — carried, untouched. (N1) is
  outside its regime for Ti, Ni and Pd (`d_eq` = 2.05–2.30 Å < `d_0`) whatever
  `ℓ` is; `W_0(Ni) = 62.5 eV` at `ℓ = 0.3 Å` is (N1) diverging, not a
  prediction. §7.9 item 2, and contact A of *both* Chapter 6 headline pairs.
- **THE FOUR DEAD NAMES** — carried from 2026-10-16, untouched today, and
  `dead_name_sweep` is still red 6/2 on exactly them.
  `graphene_potential_census_audit.CHANNEL_POTENTIAL_READINGS` (genuinely
  DEAD) plus three **hidden** dead names whose CROSS_MODULE resolutions
  resolved on *spelling*: `graphene_contact_consistent_argument.E_CHARGE`,
  `graphene_diffusion_current_model.HBAR`,
  `graphene_quantum_capacitance_argument.A_OX`. The 2026-10-01 precedent
  applies: give a reader, or declare exempt with a stated reason. **Today
  sharpens why this one matters:** the hidden three are the same
  resolve-on-spelling loophole as the docstrings in §2, one layer down.
- **CONVERT THE ELEVEN SURVIVING HAND-LISTED REBINDING SITES** — carried from
  2026-10-14, untouched for a third session. Today adds a third way a site's
  conversion can arrive or not: by import, by the scope of its `def`, and now
  by **spelling**.
- **THE PRESENT/USED/REACHED TABLE AS A TEST, NOT A REPORT** — carried, and it
  is an instance of the item ranked first above.
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, **twelve sessions
  overdue**. Today produced no new instance; every check added today is an
  exact equality.
- **THE 2–3 ULP ANTISYMMETRY IN THE SHIPPED `dWF` SUBTRACTION** — carried,
  cosmetic, a line in item 15's audit rather than a session.
- **THE CONTACT MODEL'S OWN APPROXIMATIONS** — carried from 2026-10-13: a
  vacuum-like 3.5 Å gap with no interface dipole layer, and one
  `lambda_decay = 250 nm` for every metal, both plausibly larger than the
  12.5 % correction already found and neither ever costed.
- **WHAT `n_puddle` ACTUALLY MEANS** — carried from 2026-10-12; four
  directions converge on it.
- **WHETHER THE FLOOR BELONGS INSIDE `C_q`'s ARGUMENT** — carried, costed at
  +0.1527 % in `n`.
- **A THIRD KNOB: THE ROUNDING FORM** — created 2026-10-11, untouched for a
  sixth session.
- **BUILD ITEM 17 AS A SWEEP** — three requirements, three worked instances,
  still no sweep.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  untouched for an **eleventh** session.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — untouched for the ninth session.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, untouched for an
  eighth session.
- **MECHANISE the magnitude item** — created 2026-10-03, twelve sessions
  overdue.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, untouched for a ninth
  session.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** —
  untouched for a ninth session.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.

**Commits this run:** 5 (the 20-knob conversion; the new suite committed red;
the 13-site conversion and the two corrected docstrings; the retraction, the
suite restructure and six regenerated transcripts; the clean pristine
baseline). This AUTOMATION_LOG.md entry makes **6**.

## 2026-10-18 — The list that ranks the work had never been audited, two of its items denied the tree, and one of them came from outside the repository

**Status:** The artefact every session reads to choose its work was audited for
the first time. Two of its twenty-five items assert the opposite of what the
tree contains; one of those was not stale but **revived**, from a wording this
list had itself retired eleven entries earlier, and its probable source is the
scheduled brief rather than the repository.

**Entry point, and it was not the ranked-first item.** 2026-10-17 ranked the
"census that enumerates a fault forever" item first. Before taking it, the two
items it offered as alternatives were checked against the tree — the cheapest
possible first act, never previously performed here — and both failed.

### 1. A1: Chapters 2 and 3 have been drafted for nineteen days

The 2026-10-17 list carries **"CHAPTERS 2 AND 3 ARE STILL UNDRAFTED — their
computational results have been complete for weeks. This is the largest
remaining block of pure writing"**, closing *"Worth ranking against item 16
deliberately rather than by default"*. Chapter 2 was drafted **2026-09-28**
(`5bbd5f3`, eleven sections), Chapter 3 **2026-09-29** (`c1d5abc`, ten
sections); the tree holds 375 and 330 non-blank lines and 11 and 9 `##`
sections. So the one item put forward as a competitor for first rank was
describing writing committed three weeks earlier.

**And the list had already got it right**, which is what makes this a
different fault from staleness. Measured over every dated entry:

| wording | entries |
|---|---|
| "Chapters 2–3 remain / are still undrafted" | 2026-09-17 … 09-27 (12) **and 2026-10-17** |
| "Chapters 2 and 3 **are drafted**, but §4.8.2's form of error is not audited for elsewhere" | 2026-10-01 … 10-06 (6) |
| neither | 2026-10-07 … 10-16 (10) |

Every "undrafted" entry to 09-27 was **true when written**. The list tracked
the drafting correctly for six entries, carried nothing for ten, and the entry
after those ten returned the pre-09-28 wording. **Staleness is cured by
re-checking an item; revival is cured by making retirement checkable** — and
only the second fix would have caught this.

### 2. A2: item 16's precondition (ii) was recorded open the day after it was done

2026-10-17: *"(ii) — §4.7's prose and §4.10's annotation — is still untouched
and is still writing, not plumbing."* §7.9, 2026-10-14: *"The only precondition
still open is the one that was never numerical."* Commit **`7604c05`,
2026-10-13**, is titled *"the two sentences in Sections 4.7 and 4.10 that the
measurement makes false, annotated in place"*, and both annotations are in
Chapter 4 carrying −1.2583 % (Cu), +20.8352 % (Ni), +39.2281 % (Au),
**+47.6925 % (Pd)** and `+51.0 → −203.276`.

**The narrowing matters more than the correction.** Both annotations are
written in the subjunctive and defer explicitly to the flip ("pending the item
16 flip"). What remains is the writing only the flip makes possible — the point
at which that very phrase goes stale, the committed table becomes the
superseded column, and §4.10's margin-against-input-error becomes a statement
about a margin already consumed. So **precondition (ii) cannot be discharged
ahead of the flip**: a session believing it is missing can discharge it by
asserting a state of the code that does not exist, which is §7.5.4's
"created by a fix" class in prose. 2026-10-17's operational advice is right and
its reason is wrong. **Item 16's remaining content is the flip and nothing
else.**

### 3. The cause is outside the repository, and that is the part no audit here could find

Both items trace to one mechanism: **the standing list is not a closed
system.** The scheduled brief that drives these sessions carries a frozen
description of repository state dated **2026-09-17**, containing the A1
sentence almost verbatim — and, for the other repository, claims about Phase 4
that `Design-Verification-Roadmap` retired on 09-18 and 09-19. A session that
reads its own brief as current state will re-derive items the list has retired,
and will do so **with the authority of an instruction rather than the
tentativeness of a memory**. That asymmetry is the whole mechanism.

It is recorded as a `NOT_MECHANISABLE` row with its reason rather than
silently omitted, because the brief is not in the tree. The compensating rule
splits: **(a) procedural — a claim about the state of the tree arriving from
outside the tree is re-verified against the tree before it ranks anything;
(b) mechanical — the audit below.**

**The other repository already had the instrument.**
`Design-Verification-Roadmap` built `tools/open_item_staleness_audit.py` for
this failure on 2026-10-13, and on **2026-10-17 used it to retire an item for
the identical reason** — *"it names work finished before the item started being
carried"*. It caught its instance on the same day this repository revived its
instance. The difference between the two days was one module. **An instrument
built in one place does not protect the other**, which is §7.5.2 across
repositories instead of across modules, and it generalises to every tool
either repository has that the other lacks.

### 4. The audit, and why its retired rows are the load-bearing half

`graphene_standing_list_audit.py`, **committed RED (20/1) with its transcript**
and registered in `SUITES` in the same commit (2026-10-03 wiring rule). Each
baseline row pairs an item's assertion with a **predicate over the tree**, so a
verdict is gated on the mechanism it reports rather than on text similarity —
Design-Verification-Roadmap's 2026-10-16 rule.

The trap this had to avoid is §4.6.15's: a suite that erases its own finding
the moment the finding is acted on. Once A1 and A2 leave the list, "absent"
must not read as a pass. A **retired** row therefore asserts two things — that
the wording stays gone (`REVIVED` otherwise, a fault class, and the check that
would have caught 2026-10-17 on the day) and that the tree-side fact which
justified retirement still holds (`RETIREMENT-STALE` otherwise, so truncating
Chapter 2 later re-opens the item instead of leaving a silent retirement).
`ABSENT` is also a fault class, so a list the parser cannot read cannot pass as
clean. Both plant sections force their retirement state explicitly, so the
controls are invariant under retiring a row.

### 5. The second finding: the caveat was in the docstring and the claim was in the transcript

A literature pass on the top *physics* item — a **lower** bound on `ℓ` — did
not obtain the lower bound, and found something else. `ε = 0.08 eV`, the
`graphene_second_anchor_bound` bound's headline scale, has always been
borrowed, and **four artefacts say so**: the module's docstring ("borrowed from
a sentence about a DIFFERENT quantity"), the 2026-10-15 notes file, the
2026-10-15 entry of this log in bold, and — more weakly — thesis §7.9 item 2.
**The transcript carried none of it.** Zero occurrences of `borrow`,
`different quantity` or `non-rigid`, while printing the attribution twice in
stronger words than any prose artefact used: *"ℓ ≤ 0.7024 Å at the paper's own
0.08 eV scale"*, *"66.5 % of the committed swept range is excluded by the
source"*, and a row marked `<- paper`.

The transcript is what `pristine_transcript_audit` compares, what a reader of
results sees, and what the standing item has been quoting for three sessions —
in the unqualified form. **2026-10-17 found two docstrings that denied what
their code did; this is that defect inverted, and 2026-10-01's dead-name
formula inverts with it: not a citation with no code behind it, but a
qualification with no output behind it.**

**And the caveat was itself too weak**, which is the substantive half. Read
verbatim, the source sentence is *"while Eq. (1) holds for a relatively large
separation d = 5.0 Å …, there is a small deviation of ∼0.08 eV at the
equilibrium separation d ≈ 3.3 Å which can be traced to a perturbation [… that
cannot be described as a rigid shift]"*. So the number is **(i) attached to the
other separation** — the one about which that sentence says Eq. (1) *holds* —
and **(ii) a measure of the failure of the rigid-shift picture, while `Δ_c`
is a rigid shift.** Wrong separation and wrong kind of quantity, compounding
rather than cancelling. The consequence none of the four artefacts drew:
**there is no privileged row.** The family over `ε ∈ [0.02, 0.20] eV` is the
result; the SUMMARY's single headline number should never have existed.

**No number moved.** `ℓ_max(0.08 eV) = 0.7024 Å` and the 66.5 % exclusion
stand, and the transcript's numeric multiset is a strict **superset** of its
predecessor's, verified by diff — every differing line an addition. What
changed is a printed `PROVENANCE` section giving each sourced number its
section, verbatim quote, **the separation it applies at** and the quantity it
measures; `<- paper` → `<- BORROWED, see PROVENANCE`; and the two false
attributions replaced with the superseded wording **quoted in the transcript**
rather than deleted. Four new validations, one of which is the rule made
mechanical: *no `print()` may assert the source attribution unqualified, and
quoting it as superseded is allowed only within the same string literal.*
10 passed, 0 failed.

### 6. Three defects in this session's own work, all caught by checks, all one cause

- **`p_item16_ii_untouched`'s first draft FAILED its own plant, on this
  session's own subject matter.** It tested `phrase in file.read()` for a
  sentence the markdown **wraps**, and for numbers written with a Unicode
  minus. It therefore reported §4.10 "untouched" about the file carrying the
  amendment — it **confirmed the item it was written to contradict**. Fixed by
  normalising blockquote markers, dashes and whitespace before matching; the
  FAIL is kept in the commented record above the fix.
- **The provenance check's first draft FAILED** because it forbade a phrase
  anywhere in a `print()` argument and so could not tell *asserting* "the
  paper's own 0.08 eV scale" from *quoting it as retracted* — which
  `report_provenance()` does deliberately. A check that cannot make that
  distinction forces a choice between keeping the record and passing the
  check, i.e. it pushes toward **deleting** the retracted wording, the one
  move the annotation rule exists to refuse. Scoping the supersession marker
  to the individual string literal fixes it, and has the right shape: the
  retraction travels in the same breath as the retracted words.
- **`graphene_standing_list_audit.py` shipped a DEAD module-level constant.**
  `PLANT_LOG_TRUE` was written and never read — 2026-10-01's class exactly,
  introduced by the module whose entire subject is claims nothing checks.
  Caught by running `dead_name_sweep` before committing. Given a **reader**
  rather than deleted: as `PLANT_LOG_PARTIAL` it backs new control P3b, which
  requires A1 `CONTRADICTED` with A2 `ABSENT`, pinning row independence.

**All three are the same cause, and it is three for three today: this
repository's instruments fail first by answering a weaker question than the one
asked.** With §4.6.15's two docstrings and 2026-10-17's textual grep, that is
now five recorded instances of *resolve-on-spelling*, across censuses,
docstrings, greps and checks. It is the repository's dominant failure mode and
it has never been given an item of its own. It gets one below.

### 7. Environment and automation health

Device reachable and folder connected at the **first** firing of the day
(06:5x UTC). `scipy` was absent again and was `pip install`ed (1.15.3) — **third
consecutive session**, so the 2026-10-15 note that it is present remains false
and should be treated as a standing fact: it is never there.

**New environment fact, and it closes a route that three sessions have named as
the cheapest one.** `arxiv.org` is **refused by the egress proxy in both
shells** — the device VM returns `curl: (56) Received HTTP code 403 from proxy
after CONNECT`, the cloud sandbox `connect_rejected`. So the PDF of
`arXiv:0902.1203` cannot be fetched into either place, and *"digitise Fig. 4's
upper (5.0 Å) point set"* is blocked **by policy rather than by effort**. The
ar5iv HTML was reachable and truncated for the **fourth** time, ending
mid-derivation of `Δ_tr(d)` before **Eq. (7)** — the equation Fig. 4's fit
lines follow. Recovered verbatim nonetheless: Eq. (4)'s form
`Δ_c(d) = e^(−γd)(a₀ + a₁d + a₂d²)`, the confirmation that `γ, a₀, a₁, a₂` are
reported nowhere in the text, and the `ε` sentence quoted in §5. A copy of the
PDF placed in the connected folder would open the item immediately.

`pristine_transcript_audit` was run **before** committing the transcript
regenerations and is what identified the three STALE ones; after regeneration
the only non-`IDENTICAL` row is the new suite's own deliberate red. Nothing
was deleted. One stale-file trap recorded for future sessions: a redirect to
`/tmp/pta.txt` failed with `Permission denied` while `tail` of that path
succeeded and returned **a previous session's output**, naming commit
`3212ed3` — i.e. a file in the device VM's `/tmp` outlived its session and a
read of it looked like a successful run. Session scratch belongs under `$HOME`.

**Not yet covered (candidates for future runs):**

- **FLIP THE DEFAULT TO `B_T`** — item 16, **ranked first**, and today changes
  what it contains rather than where it ranks. Precondition (ii) is **not**
  outstanding writing (§2 above); preconditions (i) and (iii) were discharged
  on 2026-10-14 and 2026-10-17. **The flip and nothing else is left**, and it
  is a Chapter 4 operation only — `graphene_contact_argument_blast_radius`
  measured zero arrivals from Chapters 5 and 6, and two models inside
  Chapter 4. Two replacements, not one rebinding: the channel takes
  `bc.rebind()`, the contact needs `bc.rebind_by_operating_point()` because
  its series capacitance is `C_interface` (65.93× `C_ox`). Costs already
  measured and committed: `n` −1.0212 %, `I_d` −0.2530 %, `n_contact`
  −4.7696 % to −12.4964 %, §4.7 `R_extra` −1.2583 % to +47.6925 %, **Pd's
  residual flips sign** and the additive decomposition fails 4-of-4.
  **Honest scope warning, which is why it was not started today:** the flip
  moves ~15 registered transcripts and dozens of quoted numbers across a
  3100-line chapter, and a half-landed flip leaves the baseline mixed between
  forms A and B_T — the one state every instrument here is least able to
  diagnose. Cheapest safe route: land the code change and the transcript
  regenerations in one session with the thesis numbers annotated-but-not-yet-
  rewritten, then the prose in the next.
- **RESOLVE-ON-SPELLING IS THIS REPOSITORY'S DOMINANT FAILURE MODE AND HAS
  NEVER HAD AN ITEM** — created today, **ranked second**. Five recorded
  instances: §4.6.15's two docstrings (audited by spelling, passed, while
  being the thing audited for), 2026-10-17's textual grep for a `setattr`
  route over computed names, and today's three — a predicate defeated by line
  wrapping, a check that could not tell assertion from quotation, and (the
  2026-10-01 original) a `CROSS_MODULE` resolution resolved against a
  function-local of the same name. The pattern is always the same: a question
  about *meaning* answered by a test on *text*. The work: enumerate every
  check in the repository whose predicate is a string operation, and for each
  either put the question to an oracle (an AST pass, a plant, a runtime
  identity test) or record a stated reason the textual form is sufficient.
  Cheapest entry: `grep -n "in _src\|in txt\|\.find(\|re\.search" *.py` over
  the suites.
- **A CENSUS THAT ENUMERATES A FAULT FOREVER IS NOT A GUARD** — carried from
  2026-10-17, where it was ranked first, **untouched today and now ranked
  third**. Today is a second instance of its shape from the other direction:
  `second_anchor_bound`'s docstring *stated* the borrowing for three sessions
  and nothing ever acted on it, which is "enumerated rather than dismissed"
  applied to a caveat instead of to a row. The 2026-10-17 work statement
  stands verbatim, and today adds its own targets: every `HONESTY`,
  `NOTE:` and `what this does not claim` paragraph in a module docstring
  should be checked for whether its qualification reaches that module's
  output.
- **A MODULE MUST NOT PLANT THE PATTERN IT CENSUSES IN THE TREE OTHER AUDITS
  READ** — carried from 2026-10-17, untouched. Today supplies a second
  instance: `graphene_standing_list_audit.py`'s plant strings are
  standing-list-shaped text now sitting in the tree, inert only because the
  parser reads `AUTOMATION_LOG.md` and not `*.py`. That is inertness by
  coincidence of scope, not by construction.
- **A LOWER BOUND ON `ℓ`** — carried, **and its status changed today**. Still
  the top *physics* item and still unobtainable from the paper's text, whose
  5.0 Å statement is a tilde with no tolerance. It needs Fig. 4's upper point
  set, which needs the PDF, which needs `arxiv.org` — **refused by the egress
  proxy in both shells** (§7). The named cheapest route is therefore closed
  until the allowlist changes or a copy of the PDF is placed in the connected
  folder, and three sessions of carrying it as "cheapest route: digitise
  Fig. 4" were carrying a route nothing could take.
- **IS `ℓ` EVEN THE RIGHT PARAMETER TO BOUND?** — created today. Under the
  source's Eq. (4) the single-exponential `ℓ` of (N1) is a **chord** decay
  length over [3.3, 5.0] Å, while `w_cross_for_metal` evaluates `delta_c` at
  each metal's `d_eq` — 2.05–3.31 Å, mostly *below* the anchor and outside the
  chord. `graphene_per_metal_crossover_model`'s docstring already records that
  (N1) drops Eq. (4)'s polynomial prefactor and that "RESULT 2 is precisely
  where it breaks"; what is not recorded is that the bound and the per-metal
  evaluation use the parameter over **different intervals**. Not costed. Cheap
  to cost: the exact relation is `1/ℓ_chord = γ − ln(P(d₂)/P(d₁))/(d₂ − d₁)`.
- **CHAPTERS 2 AND 3 — RETIRED 2026-10-18, successor named.** The undrafted
  wording is retired as false (§1) and must not return; `A1` in
  `graphene_standing_list_audit.py` now fails on its revival. The 2026-10-01
  form is the live successor: **§4.8.2's form of error is not audited for
  elsewhere**, created 10-01, untouched since 10-06, and that is the item —
  not the drafting. Both chapters may of course be *extended*; at 375 and 330
  lines against Chapter 4's 3100 they are the thinnest substantive chapters,
  but that is a judgement about depth and is not what the retired item said.
- **A `Δ_c` FORM FOR THE CHEMISORBED METALS** — carried, untouched. (N1) is
  outside its regime for Ti, Ni and Pd (`d_eq` = 2.05–2.30 Å < `d_0`) whatever
  `ℓ` is; `W_0(Ni) = 62.5 eV` at `ℓ = 0.3 Å` is (N1) diverging, not a
  prediction. §7.9 item 2, and contact A of *both* Chapter 6 headline pairs.
  **Today's new item above is its upstream half**: if `ℓ` is a chord length,
  the chemisorbed metals are being evaluated with a parameter fitted over an
  interval they are not in, independently of the `d_0` regime question.
- **THE FOUR DEAD NAMES** — carried from 2026-10-16, untouched for a third
  session, `dead_name_sweep` still red 6/2 on exactly them.
  `graphene_potential_census_audit.CHANNEL_POTENTIAL_READINGS` (genuinely
  DEAD) plus three **hidden** dead names resolved on *spelling*:
  `graphene_contact_consistent_argument.E_CHARGE`,
  `graphene_diffusion_current_model.HBAR`,
  `graphene_quantum_capacitance_argument.A_OX`. Today adds a fourth near-miss
  caught before commit, and the hidden three are instances of the
  ranked-second item above.
- **CONVERT THE ELEVEN SURVIVING HAND-LISTED REBINDING SITES** — carried from
  2026-10-14, untouched for a fourth session.
- **THE PRESENT/USED/REACHED TABLE AS A TEST, NOT A REPORT** — carried, and an
  instance of the ranked-third item.
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, **thirteen sessions
  overdue**. Today added five checks and every one is an exact equality or a
  set membership, so no new instance.
- **THE 2–3 ULP ANTISYMMETRY IN THE SHIPPED `dWF` SUBTRACTION** — carried,
  cosmetic, a line in item 15's audit rather than a session.
- **THE CONTACT MODEL'S OWN APPROXIMATIONS** — carried from 2026-10-13: a
  vacuum-like 3.5 Å gap with no interface dipole layer, and one
  `lambda_decay = 250 nm` for every metal, both plausibly larger than the
  12.5 % correction already found and neither ever costed.
- **WHAT `n_puddle` ACTUALLY MEANS** — carried from 2026-10-12; four
  directions converge on it.
- **WHETHER THE FLOOR BELONGS INSIDE `C_q`'s ARGUMENT** — carried, costed at
  +0.1527 % in `n`.
- **A THIRD KNOB: THE ROUNDING FORM** — created 2026-10-11, untouched for a
  seventh session.
- **BUILD ITEM 17 AS A SWEEP** — three requirements, three worked instances,
  still no sweep.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  untouched for a **twelfth** session.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — untouched for the tenth session.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, untouched for a
  ninth session.
- **MECHANISE the magnitude item** — created 2026-10-03, thirteen sessions
  overdue.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, untouched for a tenth
  session.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** —
  untouched for a tenth session.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.

**Commits this run:** 5 (the audit committed red with its transcript and
registered; the thesis's §4.6.16 and §7.9 item 21; the second-anchor
provenance amendment and today's notes; the in-place annotation of the borrowed
`ε` in all three prose artefacts; the row-independence control and three
regenerated census transcripts). This AUTOMATION_LOG.md entry makes **6**, and
retiring A1/A2 in the audit's baseline — which turns the suite green now that
this entry is the latest — makes **7**.

---

## 2026-10-20 — gamma cancels, and the parameter is right at no metal at all

### 1. What was taken

Not the item ranked first. **FLIP THE DEFAULT TO `B_T`** carries an explicit
scope warning — ~15 transcripts and dozens of quoted numbers across a
3100-line chapter, with a half-landed flip leaving the baseline mixed between
forms A and `B_T`, the one state every instrument here is least able to
diagnose. Today's session began under a device bridge that had already proved
unreliable (§7), and starting a flip that must land whole in one session under
those conditions is how the mixed state gets created. Taken instead:
**IS `ℓ` EVEN THE RIGHT PARAMETER TO BOUND?**, created 2026-10-18, never
costed, with its entry point already named. The choice is recorded because
skipping a ranked-first item needs a reason and "it looked big" is not one.

### 2. Gamma cancels, which is what makes the item costable

(N1) is a single exponential; Khomyakov et al.'s Eq. (4) is
`e^(−γd) P(d)` with `P` quadratic. (N1) therefore has one decay length where
Eq. (4) has a different one at every separation, and `γ, a₀, a₁, a₂` are
reported nowhere in the paper's text — four extraction passes, with
`arxiv.org` refused by the egress proxy in both shells since 2026-10-18, so
Fig. 4 cannot be digitised either. The obvious objection is that the item
cannot be costed without the numbers.

It can. With `L = ln P`,

```
1/ℓ_chord(d₁,d₂) = γ − ⟨L'⟩        1/ℓ_local(d) = γ − L'(d)
1/ℓ_chord − 1/ℓ_local(d) = L'(d) − ⟨L'⟩
```

**`γ` cancels exactly**, and since `L' = P'/P` is scale-invariant `a₀` cancels
with it. **Two of the four unreported coefficients leave the question
entirely**, and the gap is computable from `P`'s two shape ratios, which can
be swept. The relative error needs no `γ` either: `ℓ_chord` **is** the
quantity (N1) calls `ℓ`, so `ℓ_local/ℓ = 1/(1 − ℓ·Δ)`.

### 3. The result is an exact separation, and it is stronger than the item

By the mean value theorem `ℓ_chord = ℓ_local(ξ)` at one interior separation.
Measured over 204 admissible quadratic shapes:

```
ξ ∈ [3.8482, 4.4952] Å      largest tabulated d_eq = 3.41 Å (Al)
                            margin +0.4382 Å
```

**`ξ` lies above every tabulated metal, for every admissible shape.** So the
swept parameter is the local decay length at **no metal in the table**, and
the four physisorbed metals that sit *inside* the chord interval are not
spared: the gap vanishes only at `ξ`, not across the interval.

**The item was wrong in both directions.** Its `d_eq` range "2.05–3.31 Å" is
wrong at the top — Al sits at 3.41 Å — and its "outside the chord" is wrong
for Al, Ag, Au and Pt, which sit inside `[3.3, 5.0]`. But the correct charge
is a **universal over all nine metals**, not the partial one it made. Both
corrections are printed by `report_intervals()` and pinned by Validation 6c.

The cost, with `Δ` `ℓ`-independent by construction:

| metal | d_eq (Å) | Δ range (1/Å) | ℓ_local/ℓ @ 0.6 Å | bad @ 1.5 Å |
|---|---|---|---|---|
| Al | 3.41 | −0.182 … 0.767 | 0.90–1.85 | 3 |
| Au | 3.31 | −0.202 … 0.797 | 0.89–1.92 | 5 |
| Cu | 3.26 | −0.212 … 0.811 | 0.89–1.95 | 5 |
| Pd | 2.30 | −0.389 … 0.988 | 0.81–2.45 | 14 |
| Ni / Co | 2.05 | −0.433 … 1.056 | 0.79–2.73 | 19 |

"bad" counts shapes whose `ℓ_local` at that metal is **infinite or negative** —
Eq. (4)'s step is not decaying there, which no single exponential represents
at any `ℓ`. **None at `ℓ = 0.3` or `0.6 Å`; 92 at `ℓ = 1.5 Å`**, the top of
the committed sweep, and they reach the physisorbed metals too. It **compounds
with** the `d₀` regime failure rather than replacing it: Ti, Ni, Co and Pd now
fail for two independent reasons.

**15 validations, all against exactly known values**: a pure exponential gives
a bitwise-zero gap over two different chord intervals; `γ` cancels to 1e-13
across four orders of magnitude, recomputed through the two length formulae
rather than through the cancelled form; the MVT root is bracketed with a zero
residual and `1/ℓ_chord` lies between the extremes of `1/ℓ_local`; the
`d₂ → d₁` limit is first order with monotone ratios reaching 2; and
`Eq. (4)|P=const` reproduces the existing `delta_c` **bitwise**. Transcript
byte-identical across two consecutive runs.

### 4. Annotated in place, not rewritten

Thesis §7.9 item 2 listed three independent failures for the chemisorbed
metals; the fourth is annotated beside them as **widening** the item, since
unlike the other three it indicts the physisorbed metals too. The per-metal
model's docstring already conceded that (N1) drops Eq. (4)'s prefactor and
that RESULT 2 is where that breaks; what it did not say — and now says beside
the existing concession — is that dropping the prefactor changes what the
surviving parameter **means**. No number in either artefact is revised, and
item 1's bound stands: it is correct as arithmetic about the chord, and what
today adds is that the bounded quantity is not the one the per-metal
evaluation uses.

### 5. Five defects in this session's own work, all one shape

- **Validation 5 claimed `BITWISE` and got `2.13e-14 eV`.** The algebra is
  exact; the arithmetic is not the same arithmetic — `delta_c()` divides in
  metres, this file multiplied by a reciprocal in Å. Bitwise was never
  available and **the claim was the error**. Split into `5a` (matched
  arithmetic, bitwise) and `5b`, which **measures** the units route at
  **12 ULP** rather than hiding it behind a tolerance.
- **Validation 4 demanded an asymptotic order at a finite step**, requiring
  every error ratio within 0.05 of 2 and getting 2.286 at `h = 0.4 Å` — a
  textbook first-order sequence. The assertion's **shape** changed to monotone
  convergence, not its tightness; loosening the tolerance would have passed
  anything.
- **Validation 3 silently skipped an inadmissible planted shape**, reporting
  on two cases while claiming three.
- **The cost table contradicted its own prose**: `diverges = 0` on all nine
  rows while the paragraph beneath said the divergences "cluster at the SMALL
  separations". Now swept over `ℓ`, where they actually appear.
- **Validation 6c's first draft expected a shape-dependent third verdict** for
  the four metals inside the chord interval, on the strength of a throwaway
  probe reporting `ξ` reaching 3.3000 Å. **The probe was wrong** — it counted
  an endpoint where `Δ` evaluates to exactly 0.0 as an interior root, so its
  minimum was the interval's own lower bound. Computed properly `ξ ≥ 3.8482 Å`
  and there is no third verdict. **The correction made the result stronger
  than the draft assumed**, which is the rarer direction.

**All five are the same shape and it is this repository's recorded dominant
failure mode: a check or a claim answering a weaker question than the one it
was written to ask.** With the sibling repository's two today, that is seven
instances in one day across both repositories, and the ranked-second item
*resolve-on-spelling* is the general form of it.

**And a sixth, which only looking could catch.** The figure was unreadable
twice. Admissibility only requires `P > 0` on `[2.0, 5.2] Å`, so a shape may
approach a zero of `P` just inside the right edge where `L' = P'/P` diverges;
one such shape drove the left panel's autoscale to `−100 1/Å` and flattened
every other curve onto the axis. Fixed, then the metal labels overprinted the
x-axis label with Ni clipped to read "li". **Neither defect was reachable from
the transcript, which said `15 passed, 0 failed` throughout.** Both views are
clipped now with the clipping stated on the axis; nothing is dropped from the
computation, only from the view. This is the first session to render a figure
and inspect it rather than trusting the validation count, and it found two
defects on the first try — worth making routine.

### 6. Environment and automation health

**The bridge dropped mid-session and the session survived, which is the
operational finding of the day and is shared with the sibling repository's
entry.** Three sibling-repository commits in, the device went offline while a
command was running and did not return across ~6 minutes of retries. Nothing
had been pushed; the work was in the session's scratch clone under `$HOME`,
invisible to the user and lost with the session. When the bridge returned the
clone was intact and the first action was `git push`, before any further work.
**Every commit after that point was pushed immediately rather than batched to
the end of the session**, and that is now the recommended pattern: a push per
commit costs nothing and bounds the loss at one commit instead of a session.

`scipy` was absent again and was `pip install`ed (1.15.3) — **fifth
consecutive session**. It is never there; the 2026-10-15 note that it is
present remains false. One new detail worth recording: backgrounding the
install with `&` in a `device_bash` call does **not** survive the call's
return, so three attempts earlier in this session produced an empty log and no
package. Installs must run in the foreground.

A figure was copied into the connected folder to be staged and inspected, and
`_figcheck_chord.png` could not be deleted afterwards — deletion inside
connected folders is blocked and requesting permission was denied by the
safety classifier on 2026-09-17. It has been moved to `_to_delete/` for Harsh
to remove. `arxiv.org` was not retried; it is refused by policy and the
2026-10-18 finding stands. No live web search was used and none was needed —
today's result is algebra plus the repository's own tables.

**Not yet covered (candidates for future runs):**

- **FLIP THE DEFAULT TO `B_T`** — item 16, **still ranked first**, and
  untouched for a second session. Preconditions (i)–(iii) are all discharged;
  the flip and nothing else is left, and it is a Chapter 4 operation only.
  Two replacements, not one rebinding: the channel takes `bc.rebind()`, the
  contact needs `bc.rebind_by_operating_point()`. Costs already measured and
  committed: `n` −1.0212 %, `I_d` −0.2530 %, `n_contact` −4.7696 % to
  −12.4964 %, §4.7 `R_extra` −1.2583 % to +47.6925 %, **Pd's residual flips
  sign**, additive decomposition fails 4-of-4. Cheapest safe route unchanged:
  land the code change and the transcript regenerations in one session with
  the thesis numbers annotated-but-not-yet-rewritten, then the prose next.
  **Added today: it needs a session with a stable bridge, because a
  half-landed flip is the one state the instruments cannot diagnose.**
- **RESOLVE-ON-SPELLING** — carried from 2026-10-18, **ranked second**, and
  today adds five more instances (§5) for a running total of ten across both
  repositories. The work is unchanged: enumerate every check whose predicate
  is a string operation and either put the question to an oracle or record a
  stated reason the textual form suffices. Cheapest entry:
  `grep -n "in _src\|in txt\|\.find(\|re\.search" *.py` over the suites.
  **The sibling repository worked an instance of it today and found the same
  pattern in a control written to repair that very pattern.**
- **RENDER AND INSPECT EVERY COMMITTED FIGURE** — created today, **ranked
  third**. §5's sixth defect: two unreadable figures, neither visible in a
  transcript reading `15 passed, 0 failed`, both found in one look. This
  repository has ~60 committed PNGs and no session on record has looked at
  one. The work: a pass that renders each figure-producing module and
  inspects the output, starting with the ones whose data can diverge
  (anything sweeping an admissible family). Cheap, and the first sample had a
  100 % hit rate.
- **IS THE ADMISSIBLE SHAPE BOX THE RIGHT ONE?** — created today. Today's
  sweep takes `p₁ ∈ [−0.60, 0.30] Å⁻¹`, `p₂ ∈ [−0.10, 0.12] Å⁻²` with 24 of
  228 refused for `P ≤ 0`, and **the box has no justification from the
  paper** — it is wide enough to contain plausible shapes and no argument says
  it is the right width. The `ξ ≥ 3.8482 Å` margin is the headline and it is a
  **minimum over this box**. Work: either derive a constraint on `P`'s shape
  from something the paper does state, or report the headline's sensitivity to
  the box, which is a one-line sweep of the box's corners.
- **A LOWER BOUND ON `ℓ`** — carried. Still blocked on Fig. 4, which is
  blocked on `arxiv.org`, refused by the egress proxy. A copy of the PDF in
  the connected folder opens it immediately. **Today changes its value:** a
  lower bound on a chord length constrains a quantity that, per §3, is the
  local decay length at no tabulated metal — so the item is worth less than
  its ranking suggested and should be read alongside the `Δ_c` form item
  below rather than ahead of it.
- **A `Δ_c` FORM FOR THE CHEMISORBED METALS** — carried, and **promoted in
  practice by today's result.** The two items are now one question seen from
  two ends: (N1) is outside its regime for Ti, Ni, Co and Pd because
  `d_eq < d₀`, *and* its parameter is a chord length no metal sits on. The
  honest fix for both is a form with the polynomial prefactor retained, which
  needs Eq. (4)'s coefficients, which needs the PDF.
- **A CENSUS THAT ENUMERATES A FAULT FOREVER IS NOT A GUARD** — carried from
  2026-10-17, untouched for a third session.
- **A MODULE MUST NOT PLANT THE PATTERN IT CENSUSES IN THE TREE OTHER AUDITS
  READ** — carried from 2026-10-17, untouched here, but **the sibling
  repository supplied the first worked instance of avoiding it today**, by
  planting control cases as text passed to a predicate rather than as files.
  That technique transfers directly.
- **CHAPTERS 2 AND 3** — the undrafted wording stays retired as false; `A1`
  in `graphene_standing_list_audit.py` fails on its revival. The live
  successor is **§4.8.2's form of error is not audited for elsewhere**,
  created 10-01, untouched since 10-06.
- **THE FOUR DEAD NAMES** — carried from 2026-10-16, untouched for a fourth
  session; `dead_name_sweep` still red 6/2 on exactly them.
- **CONVERT THE ELEVEN SURVIVING HAND-LISTED REBINDING SITES** — carried from
  2026-10-14, untouched for a fifth session.
- **THE PRESENT/USED/REACHED TABLE AS A TEST, NOT A REPORT** — carried.
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, **fourteen sessions
  overdue**, and today supplies its first real instance rather than another
  clean sheet: Validation 4's `0.05` was a tolerance standing in for a
  statement about convergence order, and replacing the tolerance with the
  order is exactly what this item asks for. One worked example now exists.
- **THE 2–3 ULP ANTISYMMETRY IN THE SHIPPED `dWF` SUBTRACTION** — carried,
  cosmetic, and today adds a sibling measurement: the metres/Å conversion
  route costs **12 ULP** in `delta_c`. Both belong in item 15's audit.
- **THE CONTACT MODEL'S OWN APPROXIMATIONS** — carried from 2026-10-13.
- **WHAT `n_puddle` ACTUALLY MEANS** — carried from 2026-10-12.
- **WHETHER THE FLOOR BELONGS INSIDE `C_q`'s ARGUMENT** — carried, costed at
  +0.1527 % in `n`.
- **A THIRD KNOB: THE ROUNDING FORM** — created 2026-10-11, untouched for an
  eighth session.
- **BUILD ITEM 17 AS A SWEEP** — three requirements, three worked instances,
  still no sweep.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  untouched for a **thirteenth** session.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — untouched for the eleventh session.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, untouched for a
  tenth session.
- **MECHANISE the magnitude item** — created 2026-10-03, fourteen sessions
  overdue.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, untouched for an
  eleventh session.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** —
  untouched for an eleventh session.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.

**Commits this run:** 3 (the chord/local module with its transcript and
figure; today's notes; the in-place annotations of thesis §7.9 item 2 and the
per-metal model's docstring). This AUTOMATION_LOG.md entry makes **4**. Each
was pushed immediately after commit rather than batched, for the reason in §6.

---

## 2026-10-21 — Item 16's two halves do not compose, and the composition failure is a third delivery-fault mode

**Ranked-first item worked:** item 16, "flip the default to form `B_T`",
untouched for two sessions. The 2026-10-20 entry recorded that "the flip and
nothing else is left". That sentence was wrong, and what it left out is worth
more than the flip.

### 1. The result

`qca.on_consistent_argument()` (channel) and
`blast.contact_on_consistent_argument()` (contact) are each committed, each
validated and each correct **alone**. Installing both at once does **not**
give form `B_T` everywhere. `n_contact` lands **4.86× (Ni) to 10.77× (Cu)**
past the effect item 16 exists to install — −51.3850 % to −61.9143 % against
the committed contact-only form `B_T`, where the item's own effect is
−4.7696 % to −12.4964 % — **and in the same direction**, so it is not a sign
error, not a `nan`, and not out of range for a quantity in units of 10¹⁶ m⁻².
A session that landed the flip and regenerated transcripts would have
committed it as the answer.

The cause is **one attribute read**. `blast.cq_shipped_of_n` builds the
contact solve on the SHIPPED `gfet.quantum_capacitance` — a choice that is
right, for two measured reasons its own docstring gives — and reads that name
by attribute lookup at call time. The channel half rebinds exactly that name.
So with both installed, the contact's self-consistency is solved against a
`C_q` that has itself already been solved, against the wrong series
capacitance.

### 2. The general form, and why neither existing guard sees it

- **2026-10-14, `bc.rebind()`**: the set a patch must reach is the name's
  runtime identity class. Guarantees *no binding of a rebound name is stale*.
- **2026-10-17, `bc.rebind_by_operating_point()`**: a class can span two
  operating points, so it must be partitioned, and an unassigned member
  raises. Guarantees *no class member is left at the shipped object*.

Both are statements about **which names move**. Neither says anything about a
**replacement that reads one of the moved names**. Such a replacement is
**neither stale nor unassigned** — it is correctly installed, computed by the
right mechanism, reported as installed by that mechanism's own report, and it
reads a different object than the one it was validated against. **That is a
third failure mode beside the two those mechanisms were built for, and it is
invisible to both.**

### 3. Two mechanisms, and they partially cancel

(i) the inner `C_q` is the channel's `B_T` solve; (ii) that inner solve
reaches `qca.series_factor`, which reads `gfet.C_ox` at call time, and the
contact solve has `C_ox` substituted by `C_interface` for its duration.
Isolated, (i) alone is −52.4054 % to −64.2107 % and **(ii) offsets it by
+1.0204 to +2.2965 percentage points for every metal**. (ii) is the more
obviously wrong-looking of the two — a module global mutated across a nested
call is visible on inspection — so **repairing the conspicuous part alone
would have been a regression**, and the only thing distinguishing that from
progress is having measured both parts separately.

### 4. The repair, and what did NOT happen

Bind the contact replacement to the shipped object **captured at manager
entry**. `graphene_default_flip.on_default_flip()` then installs the whole
flip through the two managers the 2026-10-20 note named, and is the **first
caller ever to give `rebind_by_operating_point` a non-`None` contact
assignment** — until 2026-10-17 the contact had no replacement of its own and
the exclusion was correct exactly as long as that was true.

**18 validations, 0 failures, all against exactly known values**: the identity
assignment is bitwise neutral on all nineteen committed numbers; the contact
half reproduces the committed contact-only `B_T` **bitwise** on all seven
metals; the channel half reproduces `on_consistent_argument()` **bitwise** on
all eight census entries; `n_contact(ΔW_F = 0)` is exactly `0.0`; the
antisymmetry is bitwise; an incomplete assignment raises and rebinds nothing;
an assignment with no stated reason is refused; interleaved channel/contact
calls under one flip leave both bitwise (the stated reentrancy hazard,
measured rather than assumed); and **V8 is a positive control on the repair**
— it asserts the unrepaired flip *does* differ, so if the two ever agree the
premise is known to have expired rather than silently satisfied. Transcript
byte-identical across two consecutive runs.

**Every cost item 16 has ever quoted is reproduced exactly** by the composed
flip: `n` −1.0212 %, `I_d` −0.2530 %, `n_contact` −4.7696 % to −12.4964 %,
§4.7 `R_extra` −1.2583 % to +47.6925 %, Pd's residual +50.9489 → −203.2764
Ω·µm, additive decomposition 4-of-4 negative. That reproduction is the
evidence the defect lives in the composition and in neither half.

**`graphene_fet_model.py` is unedited and the chapter default is still
(4.29).** The flip is now a single installable validated costed object, so
what remains is a decision about what the thesis claims rather than a build.

### 5. Annotated in place, not rewritten

§4.6.17 is new. §4.6.9(g) item 1 and Chapter 7 item 16 are **annotated beside
their existing text**, with no number revised — which is the right call here
in a stronger sense than usual, since the composed flip reproduces every one
of those numbers exactly.

### 6. Three defects in this session's own work

- **An unpaired ratio.** §1's first draft divided the smallest-magnitude
  unrepaired deviation by the **largest**-magnitude committed effect and the
  other way round, quoting "4.1× to 13.0×" from two pairs of numbers
  belonging to different metals. The per-metal answer is 4.86× to 10.77×.
  **Same shape as this repository's recorded dominant failure mode** — a
  claim answering a weaker question than the one it was written to ask.
- **Two literal `%%` leaks** in `print` calls carrying no `%` operator.
- **A FALSE FIGURE TITLE.** Panel (a) read *"overshoots the effect it installs
  by 4-5x"* — the figure of the superseded unpaired calculation. The prose was
  repaired and the title was not, **because a title is not in the transcript
  and nothing reads it.** It is now computed from the same `ratios` the table
  prints. Beside it, panel (b)'s legend sat on top of Cu's two tallest bars
  and its `R_tr` annotation.

All three were caught by **reading the transcript and looking at the figure**,
not by any check in the file. This is the **second consecutive session** in
which rendering a figure and inspecting it found defects on the first try, and
the first in which one of them was a **false claim** rather than an
unreadable one.

### 7. Environment and automation health

The bridge was stable throughout, unlike 2026-10-20. Every commit was pushed
immediately after commit regardless; that is now the standing pattern and it
cost nothing. `scipy` absent again and `pip install`ed (1.15.3) — **sixth
consecutive session**; the 2026-10-15 note that it is present remains false.
`numpy` 2.2.6, `matplotlib` 3.10.9 present. No live web search was used and
none was needed — today's result is this repository's own code measured
against itself. `arxiv.org` was not retried; the 2026-10-18 finding stands. A
figure was copied into the connected folder to be staged and inspected, then
moved to `_to_delete/` for Harsh to remove, since deletion inside connected
folders is blocked and requesting permission was denied by the safety
classifier on 2026-09-17.

**Not yet covered (candidates for future runs):**

- **DOES ANY REPLACEMENT READ A NAME ITS OWN FLIP REBINDS?** — created today,
  **ranked first**, and it is ranked first because today's defect was found by
  accident rather than by a guard. The two existing mechanisms answer *which
  names move*; nothing answers *what the replacements read*. The work is one
  question per replacement — which module-level names does this function read
  at call time, and does any flip that installs it rebind one of them — and
  the repository has only a handful of replacements, so it is **cheap**.
  Cheapest entry: the replacements are `qca._carrier_density_consistent`,
  `qca._quantum_capacitance_consistent`, `blast._cq_consistent_contact` and
  `flip.make_contact_cq`; `blast._cq_consistent_contact` is the known
  positive, and the question is whether it is the only one.
- **FLIP THE DEFAULT TO `B_T`** — item 16, now **ranked second** rather than
  first, and for the first time the demotion is not a deferral. Every
  precondition is discharged, composition is repaired, the flip is one
  installable object and all costs are reproduced. What is left is editing
  `graphene_fet_model.py`, which moves committed Chapter 4 numbers and is a
  decision about what the thesis claims. **It should not be landed in the
  same session as the item ranked first**, because a flip landed over an
  unaudited replacement set is exactly the state today showed is
  undiagnosable from a green transcript.
- **RENDER AND INSPECT EVERY COMMITTED FIGURE** — carried from 2026-10-20,
  **ranked third**, and **upgraded today**: the question is no longer only
  "is the figure legible" but "does the figure assert anything the text no
  longer says". Today's false title is the first instance of the second kind,
  and it is the more dangerous one, since an unreadable figure announces
  itself and a confidently wrong one does not. ~62 committed PNGs, two
  sessions of sampling, a 100 % hit rate on both. **Any figure carrying a
  number in its title or annotations should derive it from the computation
  rather than from a string literal**, which is one mechanical sweep.
- **RESOLVE-ON-SPELLING** — carried from 2026-10-18, and today adds three more
  instances (§6) for a running total of thirteen across both repositories.
  Unchanged work: enumerate every check whose predicate is a string operation
  and either put the question to an oracle or record why the textual form
  suffices. Cheapest entry: `grep -n "in _src\|in txt\|\.find(\|re\.search"
  *.py`.
- **IS THE ADMISSIBLE SHAPE BOX THE RIGHT ONE?** — carried from 2026-10-20,
  untouched. `p₁ ∈ [−0.60, 0.30] Å⁻¹`, `p₂ ∈ [−0.10, 0.12] Å⁻²` with no
  justification from the paper; the `ξ ≥ 3.8482 Å` headline is a minimum over
  that box. Work: derive a constraint on `P`'s shape, or report the
  headline's sensitivity to the box — a one-line sweep of its corners.
- **A LOWER BOUND ON `ℓ`** — carried, still blocked on Fig. 4, which is
  blocked on `arxiv.org`. A copy of the PDF in the connected folder unblocks
  it immediately. Worth less than its old ranking suggested, per 2026-10-20.
- **A `Δ_c` FORM FOR THE CHEMISORBED METALS** — carried, and read alongside
  the item above: both need Eq. (4)'s coefficients, which need the PDF.
- **A CENSUS THAT ENUMERATES A FAULT FOREVER IS NOT A GUARD** — carried from
  2026-10-17, untouched for a fourth session.
- **A MODULE MUST NOT PLANT THE PATTERN IT CENSUSES IN THE TREE OTHER AUDITS
  READ** — carried from 2026-10-17; the sibling repository's technique of
  planting controls as text passed to a predicate transfers directly.
- **CHAPTERS 2 AND 3** — the undrafted wording stays retired as false; `A1` in
  `graphene_standing_list_audit.py` fails on its revival. The live successor
  is **§4.8.2's form of error is not audited for elsewhere**, created 10-01,
  untouched since 10-06.
- **THE FOUR DEAD NAMES** — carried from 2026-10-16, untouched for a fifth
  session; `dead_name_sweep` still red 6/2 on exactly them.
- **CONVERT THE ELEVEN SURVIVING HAND-LISTED REBINDING SITES** — carried from
  2026-10-14, untouched for a sixth session, and **today raises its value**:
  a hand-listed site is one more place a replacement can read a rebound name
  unnoticed.
- **THE PRESENT/USED/REACHED TABLE AS A TEST, NOT A REPORT** — carried.
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, fifteen sessions
  overdue; 2026-10-20 supplied its first worked example.
- **THE 2–3 ULP ANTISYMMETRY IN THE SHIPPED `dWF` SUBTRACTION** — carried,
  cosmetic, and belongs in item 15's audit with 2026-10-20's 12-ULP
  metres/Å measurement.
- **THE CONTACT MODEL'S OWN APPROXIMATIONS** — carried from 2026-10-13.
- **WHAT `n_puddle` ACTUALLY MEANS** — carried from 2026-10-12.
- **WHETHER THE FLOOR BELONGS INSIDE `C_q`'s ARGUMENT** — carried, costed at
  +0.1527 % in `n`.
- **A THIRD KNOB: THE ROUNDING FORM** — created 2026-10-11, ninth session.
- **BUILD ITEM 17 AS A SWEEP** — three requirements, three worked instances,
  still no sweep.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  **fourteenth** session untouched, and **today is an instance of it**: §6's
  unpaired ratio is exactly the fault this item exists to find.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — twelfth session untouched.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, eleventh session.
- **MECHANISE the magnitude item** — created 2026-10-03, fifteen sessions
  overdue.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, twelfth session.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** — twelfth
  session untouched.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.

**Commits this run:** 3 (the flip module with its transcript and figure; the
thesis annotations and §4.6.17; today's notes). This AUTOMATION_LOG.md entry
makes **4**. Each was pushed immediately after commit.

---

## 2026-10-22 — The item's own wording clears the defect that created it, and a capture is one level deep

**Ranked-first item worked:** *does any replacement read a name its own flip
rebinds?* — created 2026-10-21 and ranked first because that session's defect
was found **by accident**. The item is now discharged, and two of this
session's three results are things the item did not ask about.

### 1. The item, answered as asked

Of the four replacements the item names, **exactly one** reads a name that is
rebound while it is installed: `blast._cq_consistent_contact`. So the
2026-10-21 diagnosis was **complete as well as correct** — the accident found
the whole of the exposure as the item framed it. That is now a measurement
rather than a hope, and it is the least interesting of the three results.

### 2. THE ITEM'S OWN WORDING ANSWERS *EMPTY* FOR THAT SAME REPLACEMENT

The item, and §4.6.17(c) it came from, both ask whether **the flip that
installs it** rebinds a name it reads. For `blast._cq_consistent_contact`:

```
installed by      : blast.contact_on_consistent_argument()
that flip rebinds : cd.quantum_capacitance
it reads          : gfet.quantum_capacitance
intersection      : EMPTY
```

**A guard written to the item's own wording would have passed the composed
flip and reported nothing.** The defect is **cross-flip**: the name read is
rebound by a *different* manager, concurrent only because item 16 composes
them. The relation is replacement × every flip installable **concurrently**,
and the composed column is the one that answers the question.

Same shape as this repository's recorded dominant failure mode — a check
answering a weaker question than the one it was written to ask — reached this
time through **the scope of a quantifier** rather than through a ratio. The
sentence in §4.6.17(c) is annotated beside and left unrevised, since it is the
record of what was asked.

### 3. Three readers, and the controls earned their place

**A** static bytecode closure over live lookups, with closure cells resolved
from the real function object so *captured* is told apart from *live*;
branch-insensitive. **B** identity sentinels at each rebind-set member
(object substitution). **C** module-type `__getattribute__`/`__setattr__`
interception, which sees non-callables and **writes**.

B and C are blind in opposite directions — a read made from inside a name's
own module is a globals-dict lookup that never reaches the module type — and
agree **exactly** on all five replacements. **Three of the four planted
controls are flagged by B and invisible to C.** The five real replacements do
not exercise that path, so without the controls C would have looked
interchangeable with B and the blind spot would have been invisible. Two
readers agreeing because neither was asked a question it could fail is not
corroboration.

Rebind sets are derived from each manager's own `setattr` calls (entering it
under reader C), not re-derived by calling `identity_class()` again — so a
manager that grows a name is reported with it and the audit is not edited.

### 4. NEW: a write/read relation, outside every existing question

All five replacements **read** `gfet.C_ox` at call time; **three write it**.
`gfet.C_ox` is in no identity class, so it is invisible to `bc.rebind()`, to
`bc.rebind_by_operating_point()`, **and** to this item's own question — all
three are about names a flip *rebinds*. What governs composition is therefore
not a set of names but a **bipartite relation between what replacements write
and what replacements read**; 2026-10-21's mechanism (ii) is one edge of it.

### 5. 2026-10-21's two mechanisms are ONE

Channel-replacement arrivals during the contact census, with `gfet.C_ox`
recorded at each arrival: **repaired flip 0**; **unrepaired flip 192, of which
181 with `C_ox` already substituted**. (ii) cannot occur without (i) — it is a
property of the context (i)'s arrivals land in. One cause, two consequences,
181/192 carrying both. That is *why* one change repairs both, which the
percentage table of 2026-10-21 could not have said and which was being read
as a coincidence. No percentage is revised; the reading is.

### 6. NEW: CAPTURE IS ONE LEVEL DEEP

A cell protects the **name** it holds, not the names the held **object** reads
at call time.

| captured object | its own live reads | safe? |
|---|---|---|
| `gfet.quantum_capacitance` | `{gfet._thermal_rounding}` | **yes** |
| `gfet.carrier_density` | `{gfet._thermal_rounding, gfet.quantum_capacitance}` | **no** |

The two captures are the *same construction*, indistinguishable at the call
site. **§4.6.17's repair is correct and its correctness is contingent**: had
the contact needed the shipped charge relation rather than the shipped
capacitance, the same repair would have been written, would have looked right,
and would have left the defect in place. A planted captured-leaky control
holding `gfet.carrier_density` is still flagged, which demonstrates the limit
rather than arguing it.

The family, in the form its members take:

- **10-14** a repair to one binding is not a repair to the class.
- **10-17** a class can span two operating points, so it must be partitioned.
- **10-21** a replacement that *reads* a rebound name is neither stale nor
  unassigned.
- **10-22** a capture that fixes one name does not fix the names the captured
  object reads. The safe criterion is the **transitive live-read set of the
  captured object**.

Second limit, from reader A: the live read is still *in* `make_contact_cq` on
the `live=True` branch, and that function is public. **A default is not a
removal.**

### 7. Validation

**30/30, all against exactly known values.** Both instruments bitwise neutral
on the seven committed contact densities, checked before any count was used;
the repaired flip reproduces the committed contact-only form `B_T` **bitwise**
on all seven metals and the unrepaired one does not (positive control,
51.3850 %–61.9143 %); Cr's exact `0.0` at ΔW_F = 0 is **excluded from the
ratio rather than quietly divided by**, after the first draft raised
`ZeroDivisionError` on it; planted positive, negative, recursive, captured-
clean and captured-leaky controls; and V17, which re-measures reader C with
reader B's sentinels installed so that neither reader is creating the read it
reports. Transcript byte-identical across two consecutive runs.

### 8. Three defects in this session's own work — all in the figure

1. Panel (a)'s legend sat on the tallest bar (y-limit set after the legend).
2. **A FALSE LEGEND.** Panel (b)'s swatch for the "inside the `C_interface`
   window" overlay was green; the overlay it labelled was drawn dark red.
   **Second consecutive session with a figure asserting something the data
   does not say.** 2026-10-21's prescription — *derive any number in a title
   or annotation from the computation* — was followed and held; what got
   through was a **colour**, which that prescription does not cover. Repaired
   by building the legend from the bar **artists** so the two cannot drift.
   Generalisation: **any property of a mark that the legend restates should be
   read from the mark.**
3. An exact zero drawn as a zero-height bar, indistinguishable from missing
   data, in the panel whose whole point is that the repair gives exactly 0.

**Third consecutive session in which rendering a figure and looking at it
found a defect on the first try.**

### 9. Environment and automation health

Bridge stable. Every commit pushed immediately after commit and verified
against the GitHub API rather than against git's own output. `scipy` absent
again and `pip install`ed (1.15.3) — **seventh consecutive session**; the
2026-10-15 note that it is present remains false and should be treated as
settled. `numpy` 2.2.6, `matplotlib` 3.10.9 present. No live web search used
and none needed — today's subject is this repository's own code measured
against itself. A figure was copied into the connected folder to be staged and
inspected, then moved to `_to_delete/`, since deletion inside connected
folders is blocked and requesting permission was denied by the safety
classifier on 2026-09-17.

**Not yet covered (candidates for future runs):**

- **GUARD THE WRITE/READ RELATION, DO NOT JUST REPORT IT** — created today,
  **ranked first**. §4 found that `gfet.C_ox` is written by three
  replacements and read by all five, and that it is in no identity class, so
  no existing mechanism can see it. Today's audit *reports* that relation; it
  does not refuse a composition that contains an edge. The work is to turn
  the bipartite relation into a check that fails: enumerate (writer, reader)
  pairs reachable under a given flip and require each to be either impossible
  or stated with a reason, in the vocabulary `bc.rebind_by_operating_point()`
  already uses for assignments. Cheapest entry: reader C already returns the
  write set; the missing piece is reachability, and §5's arrival probe is
  the mechanism.
- **AUDIT THE CAPTURES, NOT JUST THE READS** — created today, **ranked
  second**, and it is the direct consequence of §6. Every capture in the
  repository should be checked against the **transitive live-read set of the
  captured object**, which reader A computes today for free. There are few
  captures; `flip.make_contact_cq`'s `shipped_cq` is the only one audited.
  Cheapest entry: `grep -n "def make_\|= gfet\.\|= qca\." *.py` and run
  `static_reads` on each captured object.
- **FLIP THE DEFAULT TO `B_T`** — item 16, **ranked third**. Every
  precondition is discharged, composition is repaired, the audit that was
  not to share a session with the flip is **done as of today**, and all
  costs are reproduced. What is left is editing `graphene_fet_model.py`,
  which moves committed Chapter 4 numbers and is a decision about what the
  thesis claims rather than a build. It is below the two items above only
  because both bear directly on whether the composed flip is sound, and
  §6 showed its repair is contingent in a way nobody had stated.
- **RENDER AND INSPECT EVERY COMMITTED FIGURE** — carried, **ranked fourth**,
  and today extends it: the question is not only "is it legible" and "does it
  assert a number the text no longer says" but "**does any property of a mark
  that the legend restates disagree with the mark**". ~63 committed PNGs,
  three sessions of sampling, 100 % hit rate on all three. Two mechanical
  sweeps are available: derive every number in a title from the computation,
  and build every legend from the artists.
- **RESOLVE-ON-SPELLING** — carried from 2026-10-18; running total thirteen
  across both repositories, unchanged today. Cheapest entry:
  `grep -n "in _src\|in txt\|\.find(\|re\.search" *.py`.
- **IS THE ADMISSIBLE SHAPE BOX THE RIGHT ONE?** — carried from 2026-10-20,
  untouched for a second session. `p₁ ∈ [−0.60, 0.30] Å⁻¹`,
  `p₂ ∈ [−0.10, 0.12] Å⁻²` with no justification from the paper; the
  `ξ ≥ 3.8482 Å` headline is a minimum over that box. One-line corner sweep.
- **A LOWER BOUND ON `ℓ`** — carried, still blocked on Fig. 4, which is
  blocked on `arxiv.org`. A copy of the PDF in the connected folder unblocks
  it immediately.
- **A `Δ_c` FORM FOR THE CHEMISORBED METALS** — carried; needs the same PDF.
- **A CENSUS THAT ENUMERATES A FAULT FOREVER IS NOT A GUARD** — carried from
  2026-10-17, fifth session, and **today is the strongest instance yet**:
  §4's write/read relation is currently a census, which is exactly why the
  item ranked first above exists.
- **A MODULE MUST NOT PLANT THE PATTERN IT CENSUSES IN THE TREE OTHER AUDITS
  READ** — carried from 2026-10-17. Today's module plants four controls as
  *functions*, not as text, so it does not add to the exposure.
- **CHAPTERS 2 AND 3** — the undrafted wording stays retired as false; `A1`
  in `graphene_standing_list_audit.py` fails on its revival. Live successor:
  **§4.8.2's form of error is not audited for elsewhere**, created 10-01,
  untouched since 10-06.
- **THE FOUR DEAD NAMES** — carried from 2026-10-16, sixth session;
  `dead_name_sweep` still red 6/2 on exactly them.
- **CONVERT THE ELEVEN SURVIVING HAND-LISTED REBINDING SITES** — carried from
  2026-10-14, seventh session, and today raises its value again: a
  hand-listed site is a place a replacement can read a rebound name
  unnoticed, and today's audit explicitly does **not** cover them.
- **THE PRESENT/USED/REACHED TABLE AS A TEST, NOT A REPORT** — carried.
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, sixteen sessions
  overdue.
- **THE 2–3 ULP ANTISYMMETRY IN THE SHIPPED `dWF` SUBTRACTION** — carried,
  cosmetic, belongs inside item 15.
- **THE CONTACT MODEL'S OWN APPROXIMATIONS** — carried from 2026-10-13.
- **WHAT `n_puddle` ACTUALLY MEANS** — carried from 2026-10-12.
- **WHETHER THE FLOOR BELONGS INSIDE `C_q`'s ARGUMENT** — carried, costed at
  +0.1527 % in `n`.
- **A THIRD KNOB: THE ROUNDING FORM** — created 2026-10-11, tenth session.
- **BUILD ITEM 17 AS A SWEEP** — three requirements, three worked instances,
  still no sweep.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  **fifteenth** session untouched.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — thirteenth session untouched.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, twelfth session.
- **MECHANISE the magnitude item** — created 2026-10-03, sixteen sessions
  overdue.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, thirteenth session.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** —
  thirteenth session untouched.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched.
- **A self-heating term** — created 2026-10-06.

**Commits this run:** 3 (the audit module with its transcript and figure; the
thesis annotations and §4.6.18; today's notes). This AUTOMATION_LOG.md entry
makes **4**. Each was pushed immediately after commit and verified against the
GitHub API.

## 2026-10-23 — Two findings are two edges of one relation, and a refusal nothing reaches is not a refusal

**Ranked-first item worked:** *guard the write/read relation, do not just
report it* — created 2026-10-22 §4. It is discharged, and with it the
standing item carried since 2026-10-17, *a census that enumerates a fault
forever is not a guard*, whose strongest instance 2026-10-22 recorded as §4
itself. Module: `graphene_write_read_guard.py`, 21/21.

### 1. The vocabulary is deliberately not new

This repository already had the refusal twice. `bc.rebind()` refuses an
exclusion **without a stated reason**; `bc.rebind_by_operating_point()`
refuses an **unassigned class member**, raises before anything is rebound,
and makes the reason mandatory for *left alone* as well as for *replaced*.
The guard is that same refusal applied to an **edge**: every
`(writer, name, reader)` carries a verdict from
`{GOVERNED, IMPOSSIBLE, OCCURS}` and a mandatory reason, and an unadjudicated
edge raises with nothing run.

### 2. WHAT IS NEW: EVERY VERDICT IS MEASURED, AND IN BOTH DIRECTIONS

Completeness is static — `rebind_by_operating_point` can check it because it
knows the identity class. A verdict is a claim about **reachability**, which
no static set contains, so a guard that *accepted* the claim would be a
census with a stricter tone, i.e. the item not discharged. Each verdict is
therefore checked against an order probe:

| verdict | refused when |
|---|---|
| `GOVERNED` | the name is in **no** rebind set |
| `IMPOSSIBLE` | the probe **observes** the edge |
| `OCCURS` | the probe **never** observes the edge |

**The third row is the half a census cannot have.** A guard that only catches
a wrongly-*permissive* claim answers a weaker question than it was written to
ask — this repository's recorded dominant failure mode. An `OCCURS` reason
describing an edge that has stopped occurring is a sentence about nothing,
and it looks true indefinitely because nothing ever asks it the opposite
question.

### 3. THE TWO FINDINGS ARE ONE RELATION, AND ONE CHANGE REMOVES BOTH

A name is written either by a **replacement** (10-22's `gfet.C_ox`) or by the
**flip itself**, which is every member of the rebind set (10-21's rebound
names). Two findings, two sessions — **two sources of writes into one
relation**, and composing them makes the composed column computable rather
than argued. Measured, reads traced and each edge derived from the read:

```
unrepaired (live_contact=True), 3 observed
    THE FLIP          --gfet.quantum_capacitance--> make_contact_cq   [10-21]
    make_contact_cq   --gfet.C_ox--> qca._quantum_capacitance_consistent [10-22]
    make_contact_cq   --gfet.C_ox--> itself
repaired   (live_contact=False), 1 observed
    make_contact_cq   --gfet.C_ox--> itself
```

The survivor is the contact substituting `C_ox` for `cd.C_interface` and
reading it back **inside its own call** — same operating point, restored
before return — so it is the substitution, not a leak. **Both
cross-replacement edges go away with one change.** That is 10-22 §5's "one
cause, two consequences" as a measurement instead of a reading of a
coincidence, and it is why 181/192 was never two phenomena. No number from
10-21 or 10-22 is revised; what is revised is the relation's status.

### 4. THE ITEM, ANSWERED (V21)

The repaired composition's adjudication — **the same stated reasons,
unchanged** — is **accepted** for `live_contact=False` and **refused** for
`live_contact=True`, naming three disagreeing verdicts. The difference
between the two compositions is now a `raise`, not a percentage table a
reader has to interpret.

### 5. Validation, 21/21, anchored to exactly-known values

**V4–V6** reproduce the 10-22 arrival counts rather than citing them: **0**
repaired, **192** unrepaired, exactly **181** with `C_ox` already
substituted — an exact zero a symmetry requires plus two integers measured by
a previous session independently of this code. **V8** exact zero: no
cross-replacement edge survives the repair. **V9** the unrepaired composition
has exactly 2, and they are the two reported defects. **V7** every observed
edge lies inside the static relation, so the static set is a sound
over-approximation and not a different answer. **V2–V3** the probe is
deterministic — census and edge set bitwise identical across consecutive
runs — checked before any count was used.

**V10–V13 put the DERIVER under exactly-known traces**, because a deriver
that returns few edges and a repository that *has* few are indistinguishable
otherwise: write-then-read gives the one edge, read-**before**-write gives
none (ordered, not a set intersection), a nested reader is attributed to
itself, a write outside every replacement is attributed to `THE FLIP`.
**V14–V20** are the refusals, one planted control each.

### 6. Two defects in this session's own work

1. **THE FIRST PROBE DID NOT TRACE READS.** It closed an edge on every
   *arrival* of a reader following a write — *did this replacement run after
   that name moved*, which is weaker than *did it read that name*, and wrong
   in the permissive direction. It reported **twelve** observed edges under
   the unrepaired flip, four on names the reader does not read at all, and
   the observed set was **not even a subset** of the relation. A probe whose
   over-report must be corrected by intersecting it with a static set has not
   measured reachability. Same shape as the dominant failure mode, reached
   this time through **the choice of which events to instrument** — not a
   ratio (10-22), not a quantifier's scope (10-21), but an instrument never
   asked whether it could see the thing its verdict named.
2. **A PLANTED CONTROL FOUND A CRASH IN A REFUSAL PATH THAT HAD NEVER RUN.**
   V18 — a false `GOVERNED` — raised `TypeError` instead of
   `FalseEdgeVerdict`: a name is a `(module, attr)` **tuple** and the message
   interpolated it into a single `%s`. Nothing else reaches that branch, so
   without the control the guard would have shipped with one of its three
   refusals unable to report itself **while looking correct**, because the two
   refusals the live subject exercises both work.

The family, in the form its members take:

- **10-14** a repair to one binding is not a repair to the class.
- **10-17** a class can span two operating points, so it must be partitioned.
- **10-21** a replacement that *reads* a rebound name is neither stale nor
  unassigned.
- **10-22** a capture that fixes one name does not fix the names the captured
  object reads.
- **10-23** **a refusal that nothing reaches is not a refusal.** Every verdict
  in a vocabulary needs its own planted control, not only the verdicts the
  live subject happens to exercise.

### 7. No figure, stated rather than silent

The result is a 3-edge set and a 1-edge set; a table carries that exactly and
a plot would decorate it. The last three sessions each found a defect in a
committed figure, two of them figures asserting something the data did not
say. Adding a fourth opportunity for a result this small is not worth it.
**This does not discharge the standing render-and-inspect item**, which is
about the ~63 figures already committed.

### 8. Environment and automation health

`scipy` absent again and installed (1.15.3) — **eighth consecutive** session;
the 2026-10-15 note that it is present remains false and is settled. `numpy`
2.2.6, `matplotlib` 3.10.9 present. **Two new environment faults, both of
which blocked the documented working pattern and both worked around:** the
device VM's `/sessions` filesystem is **100 % full** (9.3 G of 9.8 G, ~100
stale session directories owned by `nobody`, unreadable and undeletable from
here), so `git clone` into `$HOME/work` fails with *No space left on device*
and `pip install` fails even with `--target`, because pip builds under
`TMPDIR`, which also points there. Both are fixed by using `/tmp`, a separate
filesystem (`/dev/sda1`, 4.4 G free): clone into `/tmp`, and
`TMPDIR=/tmp/piptmp pip install --no-cache-dir --target /tmp/pylibs scipy`
with `PYTHONPATH=/tmp/pylibs`. Also: **`pgrep -f` is unusable in this VM** —
the agent proxy puts a base64 copy of the running command into `HTTP_PROXY`,
so `pgrep -f <anything>` matches its own wrapper and always reports RUNNING.
The bridge dropped once mid-session and the session resumed with `/tmp`
intact. Every commit pushed immediately after commit and verified against the
GitHub API. No live web search used and none needed — today's subject is this
repository's own composition measured against itself.

**Not yet covered (candidates for future runs):**

- **AUDIT THE CAPTURES, NOT JUST THE READS** — carried from 2026-10-22, now
  **ranked first**, and today raises its value rather than leaving it where
  it was. §6 of 10-22 showed a capture protects the name it holds and not the
  names the held object reads at call time, and that `§4.6.17`'s repair is
  therefore *correct and contingent*. Today's guard covers **reads of written
  names**; a captured object's own live reads are reached by neither the guard
  nor `bc.*`, so this is now the largest remaining unguarded surface rather
  than a tidy-up. Cheapest entry unchanged:
  `grep -n "def make_\|= gfet\.\|= qca\." *.py` and `static_reads` on each
  captured object — and the natural shape is a fourth verdict in today's
  vocabulary rather than a separate module.
- **GIVE EVERY VERDICT IN EVERY VOCABULARY A PLANTED CONTROL** — created
  today, **ranked second**, and it is the generalisation of §6.2 rather than
  a repeat of it. V18 found a crash in the one refusal path the live subject
  never exercises. The question is not about this module: it is **which
  refusal branches anywhere in the repository have never been reached by a
  test**. `bc.UnassignedClassMember` has three raise sites and the
  not-a-pair and no-reason ones are plausibly unexercised; `bc.rebind()`'s
  bare-set refusal likewise. Cheapest entry:
  `grep -n "raise " *.py` and, for each, ask what reaches it.
- **FLIP THE DEFAULT TO `B_T`** — item 16, **ranked third**, unchanged in
  content and strengthened in precondition. Every precondition was already
  discharged; today adds that the composition is not merely repaired but
  **guarded**, so the flip can no longer reintroduce either cross-replacement
  edge silently. What is left is editing `graphene_fet_model.py`, which moves
  committed Chapter 4 numbers — still a decision about what the thesis
  claims, not a build.
- **RENDER AND INSPECT EVERY COMMITTED FIGURE** — carried, **ranked fourth**.
  ~63 committed PNGs, three sessions of sampling, 100 % hit rate on all
  three. Today produced no figure, so the backlog is unchanged and the two
  mechanical sweeps are still available: derive every number in a title from
  the computation, and build every legend from the artists.
- **ADJUDICATE THE RELATION FOR THE FLIPS THIS SESSION DID NOT COMPOSE** —
  created today. The guard was exercised on `flip.on_default_flip` at both
  settings. `qca.on_consistent_argument()` and
  `blast.contact_on_consistent_argument()` have rebind sets and appear in
  `rra.FLIPS`, but no adjudication was written for compositions headed by
  them. Stated so today is not read as having guarded every composition.
- **RESOLVE-ON-SPELLING** — carried from 2026-10-18; running total **fourteen**
  across both repositories, and today's instance is in the sibling repository
  rather than this one. Cheapest entry:
  `grep -n "in _src\|in txt\|\.find(\|re\.search" *.py`.
- **IS THE ADMISSIBLE SHAPE BOX THE RIGHT ONE?** — carried from 2026-10-20,
  **third** session untouched. `p₁ ∈ [−0.60, 0.30] Å⁻¹`,
  `p₂ ∈ [−0.10, 0.12] Å⁻²` with no justification from the paper; the
  `ξ ≥ 3.8482 Å` headline is a minimum over that box. One-line corner sweep.
- **A LOWER BOUND ON `ℓ`** — carried, still blocked on Fig. 4, which is
  blocked on `arxiv.org`. A copy of the PDF in the connected folder unblocks
  it immediately.
- **A `Δ_c` FORM FOR THE CHEMISORBED METALS** — carried; needs the same PDF.
- **A MODULE MUST NOT PLANT THE PATTERN IT CENSUSES IN THE TREE OTHER AUDITS
  READ** — carried from 2026-10-17. Today's module plants its controls as
  synthetic **traces**, not as functions or text, so it does not add to the
  exposure at all — the cleanest instance yet.
- **CHAPTERS 2 AND 3** — the undrafted wording stays retired as false; `A1`
  in `graphene_standing_list_audit.py` fails on its revival. Live successor:
  **§4.8.2's form of error is not audited for elsewhere**, created 10-01,
  untouched since 10-06.
- **THE FOUR DEAD NAMES** — carried from 2026-10-16, **seventh** session;
  `dead_name_sweep` still red 6/2 on exactly them.
- **CONVERT THE ELEVEN SURVIVING HAND-LISTED REBINDING SITES** — carried from
  2026-10-14, **eighth** session. Today raises it again: a hand-listed site is
  a place a replacement can read a written name unnoticed, and today's guard
  derives its relation from measured writes and reads, so it does **not**
  cover a site that is never entered.
- **THE PRESENT/USED/REACHED TABLE AS A TEST, NOT A REPORT** — carried, and
  today is a worked reference for exactly that conversion.
- **AUDIT EVERY TOLERANCE** — item 15, created 2026-10-09, **seventeen**
  sessions overdue.
- **THE 2–3 ULP ANTISYMMETRY IN THE SHIPPED `dWF` SUBTRACTION** — carried,
  cosmetic, belongs inside item 15.
- **THE CONTACT MODEL'S OWN APPROXIMATIONS** — carried from 2026-10-13.
- **WHAT `n_puddle` ACTUALLY MEANS** — carried from 2026-10-12, and today
  gives it a concrete hook: `gfet.n_puddle` is read by
  `qca._carrier_density_consistent` and is in no rebind set, so it sits in the
  same structural position `gfet.C_ox` did before §4 of 10-22 — read across a
  composition and governed by nothing. It has no writer today, which is the
  only reason it is not an edge.
- **WHETHER THE FLOOR BELONGS INSIDE `C_q`'s ARGUMENT** — carried, costed at
  +0.1527 % in `n`.
- **A THIRD KNOB: THE ROUNDING FORM** — created 2026-10-11, **eleventh**
  session.
- **BUILD ITEM 17 AS A SWEEP** — three requirements, three worked instances,
  still no sweep.
- **CENSUS EVERY CRITERION AND CHECK THAT IS A RATIO** — created 2026-10-06,
  **sixteenth** session untouched.
- **THE `d|n|/dV` = 0 HALF OF ITEM 11** — **fourteenth** session untouched.
- **THE §4.6.5 REWIRING ITEM** — unblocked since 2026-10-10, **thirteenth**
  session.
- **MECHANISE the magnitude item** — created 2026-10-03, seventeen sessions
  overdue.
- **MEASURE WHAT EACH DECLARATION COSTS** — item 14, **fourteenth** session.
- **ASYMMETRIC ELECTRON/HOLE MOBILITY** — created 2026-10-09, untouched.
- **The hole-branch sign of Chapter 4's branch-asymmetric results** —
  **fourteenth** session untouched.
- **The seven `UNADJUDICATED` census rows** — created 2026-10-08, untouched,
  and the name is now awkward: today's module uses *unadjudicated* for an
  edge with no verdict. Worth renaming one of the two.
- **A self-heating term** — created 2026-10-06.

**Commits this run:** 3 (the guard module with its transcript; today's notes;
the thesis §4.6.19 and the 4.6.18 annotation). This AUTOMATION_LOG.md entry
makes **4**. Each was pushed immediately after commit and verified against the
GitHub API.
