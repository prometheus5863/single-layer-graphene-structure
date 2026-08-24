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
