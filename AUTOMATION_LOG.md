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
