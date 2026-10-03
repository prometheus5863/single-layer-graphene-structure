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
