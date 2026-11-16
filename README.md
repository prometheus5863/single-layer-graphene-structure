# Single-layer graphene device models: contacts, interconnects, photodetectors

Computational part of an M.Tech thesis (Engineering Physics, IIT (BHU) Varanasi, 2027).
The thesis text is in [`thesis/`](thesis/); every claim below cites the section it comes from.

## Abstract

This repository contains compact, analytical device models of single-layer graphene for three industrial
applications — metal contacts and the back-gated FET they limit (Ch 4), interconnects compared with scaled
copper (Ch 5), and two-contact photodetectors (Ch 6) — built on a nearest-neighbour tight-binding description
of the material (Ch 2–3). Alongside the models sits an explicit verification layer: validations against
exactly known limits, mutation harnesses that require a test battery to fail on an injected defect, and
audits of the audit code itself, whose findings are kept in the thesis as corrections of record (§7.5–7.7).
These are compact models checked against exact limits and against published measurements. They are not DFT
or NEGF calculations, and none of their numbers should be read as device predictions to better than the
stated error bars.

## Headline results

Each number is quoted as it stands in the cited section, including where a later section revised it: on
2026-10-16, except the first bullet, which was corrected on 2026-11-16. Section numbers did not change when the
thesis was restructured (2026-11). The audit sections they cite — §2.10, most of §4.6.x, §4.10, §6.13 and others —
now sit in [Appendix A](thesis/appendix-A-verification-record.md) under their original numbers, and §7.5–7.7 in
[Appendix B](thesis/appendix-B-failure-taxonomy.md).

* **The contact-resistance decomposition `R_c = R_transmission + R_extra` is sign-consistent, not refuted
  (corrected 2026-11-16).** With the contact doping measured from the 5.4 eV n/p crossover, Cu and Au — the
  physisorbed metals with a measured `R_c` — n-dope the contact, so `R_extra` is negative (Cu −445.04, Au −161.14
  Ω·µm) and the naive residuals are positive (Cu +629.04, Au +680.14 Ω·µm). Ni and Pd chemisorb and are outside
  the model (§4.7; `results/transcripts/ch4/contact_doping_model_output.txt`). The residuals are well-conditioned:
  Cu's sign would need a 141.3 % error in `R_extra` to flip (`results/transcripts/ch4/sensitivity_audit_output.txt`).
  The earlier headline here, negative residuals for 3 of 4 metals, was an artefact of measuring the doping from
  graphene's 4.5 eV, and it is withdrawn (Appendix A.1). The corrected result is consistency, not confirmation:
  for Cu the decomposition needs an interface term 3.4× the measured total, and Au's doping sign is not
  established (§4.7).
* **Whether graphene beats copper interconnects depends on the copper baseline, not only on graphene.** At
  realistic edge specularity (p = 0.15) the crossover is ≈130 nm against a liner-free Cu model, ≈123 nm with a
  0.3 nm Ru liner, ≈3 nm with a 1.0 nm Co liner, and absent against a 3 nm TaN/Co liner (§5.3). The
  `λ_impurity` calibration behind `ρ(W)` is the worst-conditioned step in the thesis, `κ = 19.71` (§5.5, §7.5.2).
* **The p/n crossover work function is metal-specific.** Replacing the shared 5.4 eV crossover with a
  per-metal one from the chemical interface term moves Cu's contact offset by up to 17.12 % — falsifying a
  pre-registered ≤ 10 % — and changes no metal's doping type (§6.12.3).
* **Same-sign contact pairs are an order of magnitude more sensitive than straddling pairs.** The ratio was
  65× on three pairs (§6.12.5), 30.4× on all 21 pairs (§6.13.4), and is 26.7× with the finite-difference step
  converged (§6.13.11). The separation itself, with no overlap between the two groups, has held at each step.
* **The signed collection kernel is a predicted scanning-photocurrent map.** For a symmetric pair it runs from
  +1.000 to −1.000, opposite polarity at the two contacts, which is the measured signature; for Ti/Pt it never
  changes sign (−1.837 to −1.830), a flat scan (§6.9.6). No scalar offset of the crossover can reverse a
  response sign (§6.13.5), and the sign-flip threshold of a pair is closed-form (§6.14.2).
* **This model supports no RF conclusion.** In the drift-only model `f_max/f_T ≈ ½·√(R_total/(R_g+R_s))` to
  0.36 %, so its `f_max` restates a resistance ratio (§4.6.2, §7.8.1a). A velocity-saturation term fixes the
  sign of the bias trend but supplies only 1.1361 of the ×4.8483 the comparison needs (§4.6.3, §7.8.1b), and
  the summary judgement for RF remains "unaddressed" (§7.8.5).
* **The band structure originally shipped in this repository had no Dirac point** — a minimum gap of
  11.2000 eV over Γ–K–M–Γ — and no figure could have caught it, because none was computed in a way that
  could fail (§2.10).

## Thesis map

| Ch | Topic | Status (Ch 1 §1.3) | Key modules (`src/`) | Key figures (`results/figures/`) |
|---|---|---|---|---|
| 1 | Introduction | drafted | — | — |
| 2 | Electronic structure | drafted 2026-09-28 | `graphene_band_structure`, `graphene_band_structure_audit` | `ch2/band_structure.png`, `ch2/density_of_states.png` |
| 3 | Transport and optics | drafted 2026-09-29 | `graphene_transport_properties`, `graphene_transport_optical_audit` | `ch3/quantum_hall_effect.png`, `ch3/optical_properties.png` |
| 4 | FET and contact device physics | in progress | `graphene_fet_model`, `graphene_contact_doping_model`, `graphene_edge_contact_model`, `rf_small_signal_model`, `graphene_velocity_saturation_model` | `ch4/rc_recalibration.png`, `ch4/contact_doping_profile.png`, `ch4/edge_vs_top_contact.png`, `ch4/rf_figures_of_merit.png` |
| 5 | Interconnects | in progress | `graphene_interconnect_model` | `ch5/interconnect_resistivity_vs_linewidth.png` |
| 6 | Photodetectors | in progress | `graphene_photodetector_*_model`, `graphene_per_metal_crossover_model`, `graphene_differential_crossover_model` | `ch6/photodetector_signed_carrier_response.png`, `ch6/photodetector_nonuniform_illumination.png`, `ch6/per_metal_crossover.png` |
| 7 | Discussion and outlook | drafted 2026-09-27, annotated since | — (synthesis) | `ch7/second_anchor_bound.png` |
| 8 | Towards TBG–superconductor junctions | outline only | — | — |
| A, B | Verification record; failure taxonomy | moved out of the chapter bodies 2026-11 | audits in `src/` | `audits/` |

Every module, the chapter it feeds, the outputs it writes and the repository modules it imports are in
[`docs/INVENTORY.md`](docs/INVENTORY.md).

## Verification approach

Models are validated against exactly known values — the Dirac point and the `2t` gap at M (§2.4), limits in
which each device model reduces to a closed form (§6.8.3, §6.9.2, §6.12.2), and bitwise identities where the
mathematics guarantees one (§7.6). Chapter 7 classifies the corrections of record into five failure
classes, each caught by a different kind of check (§7.5):

| class | example | detector |
|---|---|---|
| (i) implementation error | `inf − inf` in a transit-time integral | an exactly known value |
| (ii) unrepresentative sample | §6.9's `< 1.02` ceiling, violated by 14 of 21 pairs | full enumeration, pre-registration |
| (iii) analysis-layer error | an unbracketed bisection printing 21 fake roots | hand-checking one row |
| (iv) created by a fix | a bisection tolerance chosen for an eV-scale variable, reused on a 5e-11 m length | the superseded path kept as an oracle |
| (v) self-certifying procedure | a convergence estimate reporting 4.97e-14 on a 14 % error | a derived bound beside the measurement |

The audit, mutation and census modules live beside the models in `src/` (classed `audit` in the inventory).
`graphene_pristine_transcript_audit.py` re-runs registered suites in a clean clone of the committed state and
compares their output with the committed transcripts.

**Reproducibility as of this commit.** `python scripts/reproduce_all.py --check` regenerates every output and
compares it with `docs/baseline_hashes.txt`. The baseline was taken on Windows (Python 3.13.5, pinned
requirements) at `main` @ `e4fbc26`, and the 21 rows the 2026-11-16 crossover correction changed were then
re-hashed. In that environment every row is reproduced byte-for-byte, except the two audits whose stdout
contains timings by design (docs/INVENTORY.md §9.4). On Linux CI, PNG bytes differ because font rasterisation
differs between platforms. Several audits exit 1 on their own recorded FAIL verdicts, as they did on `main`.
`tests/test_contact_crossover.py` pins the contact model's exact zeros.

## Known limitations

* **The contact-doping offsets are not reconciled with measurement.** The measured graphene/electrode
  interface step is ≈0.12 eV with doping 0.2–0.3 µm into the channel; the metal table assumes 0.25–1.07 eV.
  Until that is resolved every ΔW-derived number in Ch 6 carries an unquantified systematic on top of its
  stated ±17 % (§7.8.2).
* **The photodetector model is photovoltaic only.** §6.6–6.14 model drift collection of photocarriers; the
  photothermoelectric effect is neither modelled nor discussed in Ch 6, and §6.5 notes the bolometric
  mechanism is not modelled either (see PUBLICATION_PLAN.md, candidate A, for why this matters).
* **Nearest-neighbour tight binding only** (§2.3, §2.9): no next-nearest-neighbour hopping, no trigonal
  warping, no substrate or strain.
* **No RF claim**, for the reasons above (§7.8.1a–g).
* **The `C_q`-argument correction (item 16) is costed but not yet the default** (§4.6.12–4.6.15), so the
  committed Ch 4 numbers are form A.
* **The interconnect model is analytic on both sides**: a phenomenological specularity for graphene and an
  isotropic, liner-subtracted cross-section for copper (§5.2, §5.5.4).

## Reproduce

Python ≥ 3.12 (the committed results were made with 3.13.16).

```bash
pip install -r requirements.txt
pip install -e .
python scripts/reproduce_all.py --check     # ~35 min; regenerates results/ and compares with docs/baseline_hashes.txt
```

A single module: `cd src && python graphene_interconnect_model.py`. Outputs go to
`results/figures/<chapter>/` and `results/transcripts/<chapter>/`. CI runs the same check
(`.github/workflows/reproduce.yml`). A passing check shows the code still produces the recorded outputs; it
does not show those outputs are correct (note 2026-10-02).

## Repository layout

```
src/                      every model, audit, helper and legacy module (flat; see docs/INVENTORY.md §6)
src/graphene/paths.py     results_path(chapter, name): where a module writes its outputs
scripts/reproduce_all.py  regenerate everything; --check compares with the baseline
results/figures/<ch>/     figures, by the chapter of the module that produced them
results/transcripts/<ch>/ *_output.txt transcripts
thesis/                   chapters 01–07
notes/                    dated research notes (unchanged)
AUTOMATION_LOG.md         dated session log
docs/                     INVENTORY.md, PATH_MAP.md, baseline_hashes.txt
legacy/figures/           figures drawn by the superseded demo scripts (the scripts stay in src/, INVENTORY §6)
PUBLICATION_PLAN.md       what in this repository could be published, and what it would take
```

## Background figures

The textbook figures that made up the previous version of this README (honeycomb lattice, linear dispersion,
half-integer quantum Hall plateaux, universal 2.3 % absorption) were drawn by `src/generate_plots.py` and
`src/simple_graphene_plots.py` (legacy), mostly from known closed forms rather than computed from the model (§2.10
explains why that mattered). They are kept in [`legacy/figures/`](legacy/figures/). The computed versions are
in `results/figures/ch2/` and `results/figures/ch3/`.

## Citation

See [`CITATION.cff`](CITATION.cff). Parts of the code and text were drafted with an AI assistant (Claude,
Anthropic) under the author's direction; the author reviewed and is responsible for all of it.

## License

MIT — see [`LICENSE`](LICENSE).
