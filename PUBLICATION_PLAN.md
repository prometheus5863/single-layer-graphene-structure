# Publication plan

Written 2026-10-16, after the Phase 1 reorganisation. Section numbers refer to `thesis/` at that date. Every
claim is quoted with the number it rests on and the section that states it; where the code and the thesis
disagree, both are given. "Novelty" is my reading of the closest literature found in a targeted search
(2008–2026); it is a starting point for a proper literature review, not a substitute for one.

## Summary

| # | candidate | novelty | readiness | effort | verdict |
|---|---|---|---|---|---|
| D | failure taxonomy of an AI-assisted computational pipeline | moderate — taxonomies exist (2025–26); a longitudinal, git-pinned single-project record does not | data exist, unanalysed | 5–6 wk | **go** (arXiv + workshop) |
| A | signed two-contact photocarrier model, kernel as a predicted scan, per-metal crossover | low–moderate — asymmetric metallisation and opposite-polarity scans are known; closed-form parity/threshold are new | headline revised three times; no PTE; no data comparison | 8–10 wk | **conditional go** — only after a PTE check and one digitised dataset |
| B | metal-specific R_c from a doping profile; additive decomposition fails | low — two-component contact models benchmarked on several metals exist (2015) | verdict changed under a known defect; 0.12 eV unresolved | 10–12 wk | **no-go now** |
| C | liner-dependent graphene/Cu crossover | essentially known | model ready, inputs not sourced | 3 wk | **drop** as a paper; keep in thesis |
| E | the package as a software paper | n/a | fails JOSS scope and testing bars today | 8+ wk | **no-go** |

Ranking by novelty × readiness ÷ effort (each scored 1–5, effort in weeks): **D 1.6 > A 0.4 > C 0.3 > B 0.2 > E ≈0**.

---

## Candidate A — signed, carrier-resolved two-contact collection; the kernel as a scan map; per-metal p/n crossover

**Claim.** For a two-contact graphene photodetector, the zero-bias response is a fixed, signed collection
kernel `k(x)` weighted by the illumination, so a symmetric pair gives opposite-polarity interface response
(`k` from +1.000 to −1.000) and an n/p pair such as Ti/Pt gives a flat, single-sign scan (−1.837 to −1.830,
§6.9.6), with `|N| = 1.832` for Ti/Pt (§6.8.4), robust to 0.45 % over the crossover's nominal band (§6.13.8);
`N` is exactly odd in the differential offset and even in the common one (measured to 2.3e-15, §7.6), which
gives a closed-form sign-flip threshold `τ* = W_B − W_A` (§6.14.2).

**Closest work.**
1. Mueller, Xia, Freitag, Tsang, Avouris, *PRB* 79, 245430 (2009), arXiv:0902.1479 — scanning photocurrent
   maps of contact-induced doping: p–n–p structure, doping 0.2–0.3 µm into the channel, ≈0.12 eV interface
   step. The qualitative scan signature A predicts is their measurement.
2. Mueller, Xia, Avouris, *Nat. Photonics* 4, 297 (2010), arXiv:1009.4465 — **asymmetric metallisation
   (Pd/Au vs Ti/Au interdigitated fingers) to break the mirror symmetry of the built-in field.** The design
   rule "use asymmetric metals" is therefore established practice, 16 years old.
3. Xu, Gabor, Alden, van der Zande, McEuen, *Nano Lett.* 10, 562 (2010), arXiv:0907.3173 — the
   **photothermoelectric (PTE) effect dominates** the photocurrent at a graphene interface junction.
4. Song, Rudner, Marcus, Levitov, *Nano Lett.* 11, 4688 (2011), arXiv:1105.1142 — hot-carrier transport
   theory of the photocurrent; PTE and photovoltaic contributions have distinct sign patterns.
5. Echtermeyer et al., *Nano Lett.* 14, 3733 (2014), arXiv:1402.1266 — separates PTE and photoelectric
   contributions in metal–graphene–metal detectors.
6. Shimomura et al., *Carbon Trends* 5, 100100 (2021) — shadow-mask asymmetric detector; leakage T ≈ 0.5
   (used in §6.9.3).

**Assessment: incremental reframing with two new pieces.** The opposite-polarity scan and the
asymmetric-metal design are known (2, 1). What is not in the papers above, as far as this search found: the
explicit parity of `N` and the closed-form threshold `τ*` (§6.14), the bound that illumination enters only
as a weight on a fixed kernel so no mask beats metallisation on a reinforcing pair (§6.9.1, surviving per
§7.4.2), and the 26.7× same-sign/straddling sensitivity dichotomy (§6.13.11). The decisive weakness is that
**Ch 6 is photovoltaic only and never mentions PTE**, which (3–5) show dominates zero-bias photocurrent in many
devices. A referee who knows (3) will ask first whether the kernel predictions survive a PTE term; at present
the thesis cannot answer.

**What is missing.**
* A PTE term, or a demonstrated regime (biased operation; long channel vs cooling length) where it is
  negligible, with the Seebeck profile `S(x)` driven by the same doping profile.
* Digitise and fit: Mueller 2009 line profiles (photocurrent vs position for Pd contacts and their gate
  dependence); Mueller 2010 responsivity of the Pd/Ti device; Shimomura 2021 masked vs unmasked
  photovoltage. Fit only `w_cross` or the screening length; report residuals.
* Uncertainty bands that include the 0.12 eV vs 0.25–1.07 eV offset discrepancy (§7.8.2), which currently
  sits outside every error bar in Ch 6.
* Ti, Ni, Pd: the per-metal crossover model refuses all three chemisorbed metals (§6.12.4), and the four
  smallest reversal margins (Au/Pd 0.020 eV, Ni/Au, Ni/Pd, Cr/Cu) are exactly in that unreachable set (§7.2.4).
  A paper must either restrict itself to physisorbed metals or solve this.

**Self-audit history.** Heavily revised — this is the chapter §7.7 calls the most-scrutinised: §6.7 reversed
§6.6's metal ranking; §6.8.5 contradicted §6.7's symmetric ranking; §6.9's `< 1.02` ceiling was retracted when
14 of 21 pairs violated it (§6.11.2, §7.4.2); §6.11.2's "five pairs" count was retracted (§6.13.7); the
dichotomy ratio went 65× → 30.4× → 26.7× (§6.12.5, §6.13.4, §6.13.11); a pre-registered ≤ 10 % was falsified at
17.12 % (§6.12.3); a citation was wrong for two sessions (§6.9.3). **The headline number has moved three times;
by the brief's own rule that is not ready.** What has *not* moved: the sign structure, Ti/Pt's `|N|` (three
refinements, ≤ 4.9 %), and the parity identity.

**Venue.** If PTE and one dataset are added: *Journal of Applied Physics* or *Physical Review Applied*
(full paper; a compact-model result with experimental comparison fits both), or *IEEE Journal of Quantum
Electronics*. Not a letter: the result is a framework plus a bound, not a single surprising number. Without
the PTE check: arXiv only.

**Effort and recommendation.** 8–10 weeks (PTE term 3, digitising and fitting 3, writing 2–3). **Conditional
go**, second priority: start only after a two-week PTE feasibility check shows the kernel's sign pattern
survives.

---

## Candidate B — metal-specific contact resistance from a work-function doping profile; failure of the additive decomposition

**Claim.** With the junction resistance `R_extra` computed from a spatially resolved work-function doping
profile, `R_c,measured − R_extra` is negative for Cu, Ni and Au (−787.5, −360.5, −90.2 Ω·µm) and positive only
for Pd (+51.0) (§4.7); the negatives are robust (Cu would need an 81.1 % error in `R_extra`, Pd only 9.6 %,
§4.10.4); with the `C_q`-argument defect corrected, all four are negative (Pd −203.3, §4.6.12(d)).

**Closest work.**
1. Giovannetti et al., *PRL* 101, 026803 (2008) and Khomyakov et al., *PRB* 79, 195425 (2009) — metal-induced
   doping of graphene, physisorbed vs chemisorbed; the source of the model's ΔW and crossover.
2. Xia, Perebeinos, Lin, Wu, Avouris, *Nat. Nanotechnol.* 6, 179 (2011) — origins and limits of metal–graphene
   junction resistance: contact plus junction contributions, gate dependence.
3. Chaves, Jiménez, Sagade, Kim, Riikonen, Lipsanen, Neumaier, arXiv:1507.08412 (2015) — **a physics-based
   two-component (tunnelling + potential-step transmission) contact model benchmarked against Pd, Ni and Ti
   data.** This is the closest prior work and it reports agreement, not a negative residual.
4. Russo, Craciun, Yamamoto, Morpurgo, Tarucha, arXiv:0901.0485 (2009) — Ti/Au `R_c` ≈ 800 ± 200 Ω·µm,
   independent of gate voltage and layer number.
5. Song et al., *Nano Lett.* 12, 3887 (2012) — work function of graphene under a metal electrode and its role
   in `R_c`.

**Assessment: low novelty, and the negative result is probably a modelling artefact.** A negative residual
from subtracting a computed term from a TLM-extracted `R_c` is most simply explained by §4.7's own first
hypothesis — the TLM `R_c` already contains part of the near-contact junction — or by `λ_decay = 250 nm`
over-integrating. Neither is a finding about graphene. Reference 3 shows a two-component model *can* be made
consistent with multi-metal data; a paper reporting the opposite must explain the difference, which requires
the source papers' transfer lengths (§4.7: not available).

**What is missing.** (i) Make form B_T the default (item 16) and re-run §4.7/§4.10 on it; (ii) simulate the TLM
extraction on top of the doping profile, so `R_extra` and `R_c,measured` are defined on the same region;
(iii) resolve 0.12 eV vs 0.25–1.07 eV (§7.8.2); (iv) one multi-metal dataset from a single lab and process
(the four §4.7 values come from different papers and processes). Without (ii) the result cannot distinguish
physics from bookkeeping.

**Self-audit history.** The one row that supported the decomposition (Pd) changed sign under a defect already
in the code (§4.6.12(d)); a pre-registered "at least three near-cancelling" was falsified (§4.10.3); the
counter-hypothesis "amplified input error" was eliminated (§4.10.4). The conclusion has strengthened, but its
numbers are about to change again when item 16 lands.

**Venue.** After (i)–(iv): *IEEE Transactions on Electron Devices* (compact contact models are its core
audience). Before: none.

**Effort and recommendation.** 10–12 weeks, most of it (ii) and (iv). **No-go** for the next semester; keep it
as the thesis's Ch 4 contribution, stated as a negative result about the decomposition, not about graphene.

---

## Candidate C — liner-dependent graphene-vs-Cu interconnect crossover

**Claim.** At edge specularity p = 0.15 the crossover is ≈130 nm against a liner-free Cu model, ≈123 nm with a
0.3 nm Ru liner, ≈3 nm with a 1.0 nm Co liner, and absent against a 3 nm TaN/Co liner (§5.3); so "a paper quoting
one crossover has chosen a baseline" (§7.8.3).

**Closest work.** Naeemi & Meindl (GNR interconnect compact models, *IEEE EDL/TED* 2007–2009); Xu, Li, Banerjee,
*IEEE TED* 56, 1567 (2009) — GNRs match Cu only for near-specular edges, against an ITRS Cu baseline; Barmak et
al., arXiv:2004.00971 — Cu resistivity at narrow widths and Ru/Co alternatives; the industry's liner/barrier
scaling literature that Ch 5 itself cites as the reason thin liners are pursued (§5.3).

**Assessment: already known in substance.** That Cu's effective resistivity at narrow widths is dominated by
liner/barrier consumption, and that benchmarks therefore depend on the assumed Cu stack, is the standard
framing of the beyond-Cu literature. The thin-liner resistivities used are "representative, order-of-magnitude"
values with no quantitative source (§5.3), the calibration step has `κ = 19.71` (§5.5), and the liner model's
`W − 2t` subtraction carries a ~20 % bar at 18 nm (§5.5.4). Per the brief: **drop** as a paper.

**Self-audit history.** 130 nm → "no crossover" (2026-08-28) → 3 nm / 123 nm (2026-08-29). Also an internal
inconsistency (logged below): §7.8.3 lists the thin-Ru case under "no crossover", §5.3 gives ≈123 nm.

**Venue / effort.** None as a standalone paper. Keep §5.5.3's structural result (three sensitivities exactly
zero at the calibration width) in the thesis. 0 weeks.

---

## Candidate D — failure taxonomy of an AI-assisted computational research pipeline

**Claim.** Over 43 dated sessions of a daily AI-assisted modelling pipeline, every correction of record falls
into one of five classes — implementation error, unrepresentative sample, analysis-layer error, error created
by a fix, self-certifying procedure — each with a distinct, non-overlapping detector (§7.5.1–7.5.5, table
§7.5.7); the class that recurs despite being known is the absolute tolerance that silently encodes a step size
(three sessions running, three times in one session, §7.5.6); and the only defence found is structural:
report the measured value beside a derived bound (§7.5.6). Supporting evidence: a convergence estimate
reporting 4.97e-14 on a 14.0 % error (§7.5.5); 21 fake roots from an unbracketed bisection while all five exact
validations passed (§7.5.3); a pristine-clone audit that found committed transcripts failing since their
own commit (note 2026-10-02).

**Closest work.**
1. Luo, Kasirzadeh, Shah, arXiv:2509.08713 (2025; NeurIPS 2025 listing) — four pitfalls of AI-scientist
   systems (benchmark selection, data leakage, metric misuse, post-hoc selection) in controlled experiments;
   auditing traces and code raises detection from ~55 % to ~82 %.
2. *False Summit and Silent Drift: A Failure Taxonomy and Efficiency Analysis of LLM-Assisted Multiphysics
   Simulation in an Open-Source Framework*, arXiv:2606.21841 (2026) — title located; text not machine-readable
   in this search. **The closest competitor; read it first** — if its taxonomy matches §7.5, D's novelty
   shrinks to the longitudinal record.
3. Yang, *A Judge Agent Closes the Reliability Gap in AI-Generated Scientific Simulation*, arXiv:2603.25780
   (2026) — classical mathematical checks on AI-generated simulation code cut undetected failures from 42 % to
   1.5 %; remaining failures cluster where results cannot be certified.
4. *Your Simulation Runs but Solves the Wrong Physics: PDE-Grounded Intent Verification for LLM-Generated
   Multiphysics Simulation Code*, arXiv:2605.09360 (2026) — title located.
5. ReplicationBench, arXiv:2510.24591 (2025) — agents replicating astrophysics papers.

**Assessment: moderately new.** Taxonomies of AI-generated scientific code failures now exist (2–4), mostly
from benchmark tasks. What this repository has that they do not: a **longitudinal record of one real
project** (2026-08-21 → 2026-10-15), every failure dated, with pre-registered predictions committed to git
before the code ran (§7.5.2), mutation harnesses with positive and null controls, and audits of the auditors
(notes 2026-09-28 → 2026-10-15). Class (v), self-certification, and the "a control that never fired has no
evidence it can" rule (note 2026-09-30) are the most transferable findings. The weaknesses are n = 1, no
counterfactual (would a human-only pipeline fail the same way?), and the author as both subject and
analyst.

**What is missing.** (i) A coding of every correction of record in `AUTOMATION_LOG.md` and `notes/` into the
five classes by someone other than the pipeline (the author, blind to Ch 7's assignments, then compared);
counts per class, detection latency in sessions, and which detector caught each. (ii) The repository's own
current state as data: at `main` 77a2f3a, 12 of 35 committed transcripts do not reproduce and 3 standing audits fail
(docs/INVENTORY.md §3) — that the verification layer itself drifts is the most important result and must be
in the paper. (iii) A clear AI-use methods section (below). (iv) Position against (2) explicitly.

**Self-audit history.** This candidate *is* the audit history; its own claims have been revised within Ch 7
(§7.5.6 grew from one instance to five). Its risk is the opposite of A's: it is a narrative assembled by the
system it describes.

**Venue.** arXiv (cs.SE + physics.comp-ph) first; then a workshop on AI for science or agents for science at a
2027 ML conference, or *Computing in Science & Engineering* (IEEE/AIP) as a department article. A 2025 venue
built on AI-first authorship, Agents4Science, exists but listed no 2026 edition when checked; its rules
(AI authorship *required*) do not fit a human-authored case study.

**Effort and recommendation.** 5–6 weeks. **Go — first priority.** It needs no new physics, the data already
exist, and it is the one candidate whose closest competitors are months, not years, old.

---

## Candidate E — the package itself (JOSS-style)

**Assessment: no-go.** JOSS requires an obvious research application with evidence of research use, a
public development history, tests and CI, and excludes "a one-off tool for a single analysis"; its 2026
policy permits generative AI with a mandatory disclosure section. This repository has no tests directory, CI
only as of this PR, 12 of 35 transcripts not reproducing, a flat module layout, and one user. Revisit only if
another group uses the contact or photodetector models.

---

## Ranking and 12-week plan

| | novelty (1–5) | readiness (1–5) | effort (wk) | score |
|---|---|---|---|---|
| D | 3 | 3 | 5.5 | **1.6** |
| A | 2 | 2 | 9 | 0.4 |
| C | 1 | 3 | 3 (as a note) | 0.3 — dropped on novelty |
| B | 2 | 1 | 11 | 0.2 |
| E | 1 | 1 | 8+ | ≈0 |

Week 1 is the first week of the thesis-only semester (≈ January 2027 — confirm the date with the department).

| week | D (lead) | A (second) |
|---|---|---|
| 1 | Read arXiv:2606.21841 and 2603.25780 in full; write a one-page positioning note. Freeze the session set (2026-08-21 → 2026-10-15). | — |
| 2 | Extract every correction of record from the log and notes into a table (date, claim, class, detector, sessions to detection). | Read Xu 2010, Song 2011, Echtermeyer 2014 for the PTE formulation. |
| 3 | Blind re-coding by the author (without Ch 7 open); agreement statistic; resolve disagreements in writing. | PTE feasibility: add `S(x)` from the existing doping profile; does the symmetric/n-p scan sign pattern survive? **Go/no-go for A at end of week 4.** |
| 4 | Quantify the drift of the verification layer itself (INVENTORY §3) and the cost of each detector (runtime, lines). | PTE feasibility continued; decision. |
| 5 | Draft: introduction, pipeline description, AI-use methods section. | (if go) Digitise Mueller 2009 line profiles and their gate dependence. |
| 6 | Draft: results by class, the class (v) cases, limitations (n = 1, no counterfactual). Internal review by supervisor. | Digitise Mueller 2010 and Shimomura 2021. |
| 7 | Revise; post to arXiv; submit to the chosen workshop or CiSE. | Fit `w_cross`/ℓ per dataset; residuals; uncertainty band including the 0.12 eV issue. |
| 8 | Buffer / respond to workshop reviews. | Physisorbed-only results table (Cu, Au, Pt); state the Ti/Ni/Pd restriction explicitly. |
| 9 | — | Draft: model, kernel, parity and threshold. |
| 10 | — | Draft: comparison with data, PTE discussion. |
| 11 | — | Supervisor review; revise. |
| 12 | — | arXiv + submission (J. Appl. Phys. or Phys. Rev. Applied). |

Thesis writing runs alongside; the plan assumes ~50 % of each week on it.

## What to stop doing

**Stop adding §4.6.x audit sub-sections by default.** §4.6 now runs from 4.6.1 to 4.6.15 (4.6.15 was added on `main` the same day this plan was written) and Ch 4 is ~27k words;
most of 4.6.5–4.6.15 is the item-16/`C_q`-argument thread and the audits of the instruments used to cost it. It
records real defects found and fixed, but it reads as a lab notebook inside a chapter. Recommendations:

* **Thesis body (Ch 4) keeps:** §4.1–4.5 (model), §4.6 core RF result with §4.6.2's resistance-ratio finding
  stated once, §4.7 and §4.10 (the negative result and its conditioning), §4.8, and a one-page summary of the
  item-16 correction and its effect on §4.7 (4 of 4 negative).
* **Move to an appendix ("Verification record for Chapter 4"):** §4.6.1, §4.6.5–4.6.11 and §4.6.13–4.6.15 in
  their current form, unchanged, so the annotate-don't-rewrite rule still holds — the body points to them.
* **Move to the paper-D track (as data, not as thesis text):** §4.6.12–4.6.13's import-binding and blast-radius
  measurements, and §7.5–7.7's taxonomy.
* **Change the automation's default target.** The highest-value open items are not new audits: (1) make
  item 16 (form B_T) the default and regenerate Ch 4–6 numbers in one commit; (2) regenerate or annotate the 15
  non-reproducing transcripts and fix the three failing standing audits, so the verification layer is green at
  HEAD; (3) the PTE feasibility check for A. A new audit section should be written only when an audit changes
  a number a chapter reports.

## Thesis-scope check — confirm with your supervisor

The thesis is formally titled around **quantum transport in twisted bilayer graphene–superconductor
junctions**. This repository is **single-layer graphene device physics** (contacts, interconnects,
photodetectors) and contains no bilayer, twist, moiré, superconducting-pairing or Andreev-reflection code.
That is a mismatch an examiner will see on page one. Questions to settle:

1. Is the formal title changing to match this work, or is the TBG–SC topic still the registered thesis? What
   is the department's deadline and procedure for a title change?
2. If TBG–SC remains the topic: is this repository Part I (single-layer preliminaries) and is there a TBG–SC
   part still to be written next semester? Which chapters survive, and does Ch 7's "industrial application
   space" framing still fit?
3. If the title changes: does the supervisor's group and the M.Tech regulations accept a thesis whose central
   contribution is device modelling plus a verification methodology (candidate D's content)?
4. Is publishing D (a paper about the AI-assisted process of the thesis) compatible with the department's
   rules on AI use in thesis work, and does it need to be declared in the thesis itself?

## Authorship and AI-use disclosure

The repository's commits carry Claude as a co-author. No venue below allows an AI tool to be an author; all put
responsibility on the human authors and require disclosure of substantive use. What each asks (checked
2026-10-16; confirm on the venue's page at submission):

| venue | rule | where to state it |
|---|---|---|
| arXiv | AI tools cannot be authors; authors take full responsibility; report significant use of generative AI per field norms | in the paper (methods / acknowledgements) |
| APS (Phys. Rev. Applied, PRB) | AI cannot be an author; disclose substantive use — tool, version, how it assisted, how output was verified; code that materially affects results counts | Methods for research use; figure captions for figures; acknowledgements otherwise (per a 2026 summary of APS policy — read APS's own page) |
| AIP (J. Appl. Phys.) | not checked in this pass | check before choosing |
| IEEE (TED, JQE, CiSE) | disclose AI-generated content: the system, the sections, and the level of generation, for text, figures, images and code | acknowledgements |
| Elsevier (Carbon) | AI not an author; separate "Declaration of generative AI and AI-assisted technologies" at submission; AI used in research methods described in Methods with tool and version | dedicated declaration + Methods |
| ACS (Nano Lett.) | authors responsible; disclose tool and how used; for code/analysis confirm verification and link scripts | acknowledgements and/or Methods |
| JOSS | AI allowed with a mandatory disclosure section: tools and versions, where used, scope, and an affirmation that humans reviewed, modified and validated all AI output | AI usage disclosure section of the paper |

**What you must state, for any of these.** That the models, analysis code, audits and much of the thesis and
note text were drafted by an AI assistant (Claude, Anthropic; give the model versions from the commit
trailers) in daily automated sessions under your direction; which parts you wrote or substantially rewrote;
how outputs were verified (exact-value validations, mutation harnesses, pristine-clone transcript audit,
independent re-derivation of the headline numbers you rely on); and that you take full responsibility. For
candidate D the AI involvement *is* the subject, so it goes in the methods in full, not in a footnote.
Remove Claude from any author list; credit it in the disclosure.

## Defects found during reorg

Logged, not fixed: Phases 1–2 made no physics or numeric change.

1. **Committed transcripts that do not reproduce.** At `2a72139`: 15 of 34 (docs/INVENTORY.md §2–3) — STALE (7):
   `cross_module_delivery_audit`, `dead_name_sweep`, `mutation_arrival_probe`, `figure_provenance_audit`,
   `log_sensitivity_step_audit`, `potential_census`, `potential_declaration` (overflow onset 18.3493 V →
   40.0000 V after the item-13 fix); DRIFT (5): `default_scale_probe_audit`, `quantum_capacitance_overflow`,
   `quantum_capacitance_argument`, `fet_signed_carrier`, `rootfinder_audit`; records (2): `audit_output.txt`,
   `sensitivity_converged_output.txt`; by design (1): `pristine_transcript_audit`. The 2026-10-16 session on
   `main` regenerated three of the stale ones; **at `77a2f3a`: 12 of 35** (§8) — and it made a new one stale,
   `import_binding_sweep`, by adding a module the sweep counts.
2. **Three standing audits FAIL** (at both commits). `graphene_dead_name_sweep.py`: 1 DEAD constant
   (`CHANNEL_POTENTIAL_READINGS`, `graphene_potential_census_audit.py:310`) and 12 unsupported CROSS_MODULE
   resolutions. `graphene_mutation_arrival_probe.py`: 8 MIXED + 2 FROZEN captured defaults in modules added
   after 2026-09-30. `graphene_potential_census_audit.py`: potential-like names missing from
   `POTENTIAL_REGISTRY`. (`graphene_quantum_capacitance_overflow.py` also exits 1, on a FAIL it declares
   deliberately.)
3. **`graphene_pristine_transcript_audit.py` exits 0 while reporting FAIL checks** — 3 of 10 at `2a72139`,
   2 of 13 at `77a2f3a` (1 STALE, 1 DRIFT) — so no caller relying on the exit code notices. The 2026-10-16 log
   entry records "13 IDENTICAL, 0 STALE" for this audit; in the reference environment it is 11 IDENTICAL and
   1 STALE, so that session's run and this one were not in the same environment.
4. **`requirements.txt` could not run the code** (numpy 1.24.3 lacks `np.trapezoid` and has no wheel for the
   Python ≥ 3.12 the code needs). Fixed in Phase 1 because it is packaging, not physics.
5. **Two writers for one file.** `generate_plots.py` and `graphene_band_structure.py` both wrote
   `band_structure.png`, `density_of_states.png`, `honeycomb_lattice.png`; `generate_plots.py` and
   `graphene_transport_properties.py` both wrote `quantum_hall_effect.png`. The committed honeycomb and QHE
   figures — shown in the old README — were the legacy versions. Resolved by the move (legacy figures in
   `legacy/figures/`; canonical ones regenerated into `results/`).
6. **Figure provenance in Ch 6.** §6.1 attributes `optical_absorption.png` to `graphene_transport_properties.py`;
   that module writes `optical_properties.png`, and `optical_absorption.png` is drawn only by
   `legacy/generate_plots.py`.
7. **Ch 5 vs Ch 7 inconsistency.** §7.8.3 lists the thin-Ru-liner case under "no crossover at all"; §5.3 gives
   ≈123 nm for it.
8. **Unannotated downstream claims.** §7.8.2 still states "fails for three of four metals" without the
   §4.6.12(d) note that the corrected model gives four of four; §7.8.4's "≈4× per incident photon" rests on the
   §6.9 argument that §7.4.2 says was retracted (the recommendation survives on other numbers, which §7.8.4
   does not cite).
9. **Fossil figures** with no producer at HEAD: `band_structure_audit.png`, `default_scale_probe_audit.png`,
   `test_plot.png`.
10. **`graphene_fmax_shortfall_decomposition.py` line 21**: invalid escape sequence `\_` in a docstring —
    a `SyntaxWarning` now, scheduled to become an error in a future Python.
11. **Docstrings that over-claim outputs**, e.g. `graphene_double_regularisation.py` says it writes
    `double_regularisation_output.txt`; it prints, and the transcript was made with `tee`.
12. **Six canonical outputs were never committed** (`fermi_surface.png`, `linear_dispersion.png`,
    `mobility_vs_temperature.png`, `conductance_map.png`, `klein_tunneling.png`, `optical_properties.png`);
    they are now in `results/`.
13. **Committed figures stale against the code that claims to produce them** (pixels differ from a HEAD
    regeneration that is itself bit-reproducible): `log_sensitivity_step_audit.png` (25 % of pixels),
    `velocity_saturation_model.png` (15 %), `differential_crossover.png` (14 %),
    `perfect_contact_counterfactual.png` (13 %), `gds_quadrature_audit.png` (2 %); a sixth,
    `cross_module_delivery_audit.png`, was regenerated on `main` on 2026-10-16. The committed versions are kept;
    `scripts/reproduce_all.py` regenerates the current ones.
14. **Five legacy demo scripts cannot leave `src/` without changing three audit transcripts**, because those
    audits count every `*.py` in their directory (`docs/INVENTORY.md` §6). This is a design coupling, not a
    bug, but it means the audits' populations are defined by a directory listing rather than a declaration.

## Sources

Journal volume/page numbers given in the text without a link were written from memory and should be checked against the DOI before citing.


* [Mueller et al., The role of contacts in graphene transistors: a scanning photocurrent study (arXiv:0902.1479)](https://arxiv.org/abs/0902.1479)
* [Mueller, Xia, Avouris, Graphene photodetectors for high-speed optical communications (arXiv:1009.4465)](https://arxiv.org/abs/1009.4465)
* [Xu et al., Photo-thermoelectric effect at a graphene interface junction (arXiv:0907.3173)](https://arxiv.org/abs/0907.3173)
* [Song et al., Hot carrier transport and photocurrent response in graphene (arXiv:1105.1142)](https://arxiv.org/abs/1105.1142v1)
* [Echtermeyer et al., Photo-thermoelectric and photoelectric contributions… (arXiv:1402.1266)](https://arxiv.org/abs/1402.1266)
* [Chaves et al., A physics based model of gate tunable metal-graphene contact resistance… (arXiv:1507.08412)](https://arxiv.org/abs/1507.08412)
* [Xia et al., The origins and limits of metal–graphene junction resistance (Nat. Nanotechnol.)](https://www.doi.org/10.1038/NNANO.2011.6)
* [Russo et al., Contact resistance in graphene-based devices (arXiv:0901.0485)](https://arxiv.org/abs/0901.0485)
* [Xu, Li, Banerjee, Modeling, analysis, and design of graphene nano-ribbon interconnects (IEEE TED 2009)](https://web.ece.ucsb.edu/Faculty/Banerjee/pubs/TED_2009_GNR.pdf)
* [Barmak et al., Epitaxial metals for interconnects beyond Cu (arXiv:2004.00971)](https://arxiv.org/abs/2004.00971)
* [Luo, Kasirzadeh, Shah, The More You Automate, the Less You See (arXiv:2509.08713)](https://arxiv.org/abs/2509.08713)
* [False Summit and Silent Drift (arXiv:2606.21841)](https://arxiv.org/pdf/2606.21841)
* [Yang, A Judge Agent Closes the Reliability Gap… (arXiv:2603.25780)](https://arxiv.org/abs/2603.25780)
* [Your Simulation Runs but Solves the Wrong Physics (arXiv:2605.09360)](https://arxiv.org/pdf/2605.09360)
* [ReplicationBench (arXiv:2510.24591)](https://arxiv.org/pdf/2510.24591)
* [Agents4Science 2025](https://agents4science.stanford.edu)
* [JOSS policies](https://joss.readthedocs.io/en/latest/policies.html) · [JOSS submitting](https://joss.readthedocs.io/en/latest/submitting.html)
* [arXiv policy on generative AI tools](https://blog.arxiv.org/2023/01/31/arxiv-announces-new-policy-on-chatgpt-and-similar-tools)
* [IEEE author guidelines for AI-generated text](https://open.ieee.org/author-guidelines-for-artificial-intelligence-ai-generated-text/)
* [Elsevier generative AI policies for journals](https://elsevier.com/about/policies-and-standards/generative-ai-policies-for-journals)
* [Summary of APS / PRL AI policy (third-party)](https://manusights.com/blog/physical-review-letters-ai-policy) · [APS announcement, June 2026](https://www.aps.org/about/news/2026/06/releases-updated-ai-policy-journals)
* [Summary of ACS Publications AI policy (BioRender policy check)](https://policycheck.tools.biorender.com/publisher/acs-publications)
