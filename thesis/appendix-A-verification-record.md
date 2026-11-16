# Appendix A: Verification record

This appendix holds the thesis's verification record: the audit of each model, the audits of those audits, and every annotation that recorded a change. It was moved out of the chapter bodies on 2026-11-16 so that the body can be read in one sitting. Nothing here was rewritten. Each moved section keeps its original number, which is also its link anchor (§4.6.1 is `#s4-6-1`). Where a moved section changed a number reported in the body, the body carries one sentence and a link here.

<a id="corrections"></a>

## A.1 Corrections of record

Every correction, retraction and supersession recorded in the draft, in date
order. *Old value* is what the thesis or code reported before the correction;
*new value* is what it reports after. The full account of each is at the linked
section, which keeps the superseded text as it was written.

| date | claim | old value | new value | cause | detecting check | record |
|---|---|---|---|---|---|---|
| 2026-09-17 | §6.6 metal ranking for photocarrier collection (Pt best) | Pt best at 1.47× | for a symmetric device Pt is worst, overstated 5.5× | single-junction model applied to a two-contact device | self-consistent two-contact model (`graphene_photodetector_two_contact_model.py`) | [§6.6](#s6-6) |
| 2026-09-18 | §6.7 symmetric metal ranking (reversed) | Cu ~ Ti > … > Pt | not claimed; at the 5.4 eV crossover Pt returns to the top (0.556) | crossover placed at graphene's 4.5 eV | signed, carrier-resolved model (`graphene_photodetector_signed_carrier_model.py`) | [§6.7](#s6-7) |
| 2026-09-20 | §6.9.5 mask-gain ceiling "under 1.02" for every asymmetric pair | 2 pairs checked, both under 1.02 | 14 of 21 pairs violate it; worst Au/Pd at 22.4 | sample drawn only from straddling pairs | full enumeration of the 21 pairs (§6.11.2) | [§6.9 note](#s6-9-note) |
| 2026-09-21 | Cu ΔW moves ≤10% under a per-metal crossover (pre-registered) | ≤10% | 17.12% | prediction calibrated on §6.11's 4.9% | per-metal crossover model (`graphene_per_metal_crossover_model.py`) | §6.12.3 (body) |
| 2026-09-22 | §6.12.5 same-sign vs straddling amplification | 65× | 30× (well-defined statistic, all 21 pairs) | the ratio divides by an exactly zero input change for straddling pairs | exact validation found mid-session (§6.13) | [§6.13](#s6-13) |
| 2026-09-22 | §6.11.2 "five pairs gain more than 2×" | 5 | 3 → 4 → 5 → 5 → 6 across δ = −0.2…+0.2 eV | a threshold count evaluated at one point of an uncertain parameter | scalar-offset sensitivity band (`graphene_crossover_sensitivity_model.py`) | [§6.13.7](#s6-13-7) |
| 2026-09-25 | §6.13.4 same-sign sensitivity range and dichotomy | 0.6625 – 2.5628; 30.4× | 0.5817 – 2.5662; 26.7× | finite-difference step 3.5 decades above the derivative's plateau | numeric-default audit (`graphene_default_scale_audit.py`) | [§6.13.11](#s6-13-11) |
| 2026-09-28 | Ch 2 band structure: gap at the point labelled K | 14.4019 eV | 0 (2.5 × 10⁻¹⁵ eV) | K point coded without its factor π | exact-value audit (`graphene_band_structure_audit.py`) | [§2.10](#s2-10) |
| 2026-09-29 | Ch 3 universal optical conductivity (T1) | 3.8235 × 10⁻⁴ S (14.4044% absorption) | 6.0853 × 10⁻⁵ S (2.2925%) | `ħ` written for `h` | exact-value audit (`graphene_transport_optical_audit.py`) | [§3.9](#s3-9) |
| 2026-09-29 | Ch 3 conductance quantum and minimum conductivity (T2) | G₀ = 4.8683 × 10⁻⁴ S; G_min 19.74× the hallmark | G₀ = 7.7481 × 10⁻⁵ S; G_min = 4e²/(πh) | 2π twice | same audit | [§3.9](#s3-9) |
| 2026-09-29 | Ch 3 Pauli-blocking branch (T3) | constant 0.5 σ₀ at every wavelength | branch reachable; §3.6 computed for the first time | joules compared with a bare eV number | same audit | [§3.9](#s3-9) |
| 2026-09-29 | Ch 3 mobility temperature exponent (T4) | 1.5000 | 1 | three-dimensional exponent | same audit | [§3.9](#s3-9) |
| 2026-10-01 | §4.8 edge-contact reduction compared against the model's 2.59× | 11.53× (on-state column) | 3.01× (Dirac-point column, the one the model's physics selects) | the measured table's second column was never read | dead-name sweep (`graphene_dead_name_sweep.py`) | [§4.8.2](#s4-8-2) |
| 2026-10-03 | §4.6 literature-scale intrinsic peak f_max | 18.786 GHz | 13.094 GHz | gate resistance captured W as a definition-time default (40× too small) | cross-module delivery audit (`graphene_cross_module_delivery_audit.py`) | [§4.6 note b](#s4-6-note-b) |
| 2026-10-03 | §4.6 literature-scale extrinsic peak f_max | 3.1885 GHz | 2.2196 GHz | same | same | [§4.6 note b](#s4-6-note-b) |
| 2026-10-03 | §4.6 f_max/f_T at the literature-scale device | 0.928 | 0.647 | same | same | [§4.6 note b](#s4-6-note-b) |
| 2026-10-04 | §4.6 attribution of the f_max shortfall to "g_ds and R_g·C_gd" | two causes | R_g·C_gd carries 0.2154 %; the shortfall is g_ds | attribution written without measurement | exact decomposition (`graphene_fmax_shortfall_decomposition.py`) | [§4.6.2](#s4-6-2) |
| 2026-10-08 | §4.6.5 diffusion term as missing physics | +0.50 % missing | a double count, +0.4279 % at the terminals | the channel potential had two readings | declaration of V_ch (`graphene_potential_declaration.py`) | [§4.6.6](#s4-6-6) |
| 2026-10-10 | (4.29) described as the linearisation of (4.28) | "linearisation" | the exact differential relation used algebraically | the two forms do not share a leading order | re-derivation (`graphene_exact_charge_rederivation.py`) | [§4.6.8](#s4-6-8) |
| 2026-10-11 | §4.6.8 cost decomposition at the RF bias (argument / quadrature / total) | +4.7416 % / +2.6396 % / +5.2064 % | +2.6408 % / +0.7603 % / +3.0963 % | each figure differenced against a T = 0 relation | 2×2 at one temperature (`graphene_quantum_capacitance_argument.py`) | [§4.6.9](#s4-6-9) |
| 2026-10-13 | §4.7 Pd naive transmission residual | +51.0 Ω·µm | −203.276 Ω·µm under the C_q-argument correction | C_q evaluated at the work-function potential, outside the harness | blast-radius audit (`graphene_contact_argument_blast_radius.py`) | [§4.6.12](#s4-6-12), [§4.7 notes](#s4-7-note-a) |
| 2026-10-29 | §4.6.7 disorder floor e·μ·n_puddle called exact | exact | exact only at μ_e = μ_h; overestimates by μ̄/√(μ_eμ_h) off it | equal branch mobilities assumed | two-carrier minimisation (`graphene_branch_asymmetric_mobility.py`) | §4.11.3 (body) |
| 2026-11-03 | §4.11.5 crossover density p ≈ 3 × 10¹² cm⁻² read as a material constant | 3 × 10¹² cm⁻² | proportional to the charged-impurity density; high by at least 13 % | n_c moves with the gate stack (A = n_c·μ̄_L is the invariant) | bound from Zion *et al.*'s own numbers (`graphene_kappa_impurity_bound.py`) | [§4.11.5a](#s4-11-5a) |
| 2026-11-16 | §4.7 Table 4.7 naive transmission residuals; the additive decomposition "fails robustly" | Cu −787.5, Ni −360.5, Au −90.2, Pd +51.0 Ω·µm (Pd −203.276 under the C_q-argument correction) | Cu +629.04, Au +680.14 Ω·µm (physisorbed); Ni +695.99, Pd +719.78 outside the model; conclusion withdrawn | contact doping measured from graphene's 4.5 eV, not the 5.4 eV n/p crossover §4.5 states | review of §4.5 against the code; now pinned by `tests/test_contact_crossover.py` | [§4.7 superseded](#s4-7-superseded), [§7.8.2](#s7-8-2-superseded) |

## A.2 Chapter 1: Introduction

<a id="s1-3"></a>

*Moved from Chapter 1, §1.3: the chapter-status table as it stood on 2026-11-16, with its dated history.*

### 1.3 Thesis structure (current status)

| Chapter | Topic | Status |
|---|---|---|
| 2 | Electronic band structure of single-layer graphene (tight-binding model, Dirac dispersion, density of states) | **DRAFTED 2026-09-28** (`thesis/02-electronic-structure.md`, eleven sections), and the drafting immediately falsified this row's previous status, which read "Computational results complete (`graphene_band_structure.py`)" and had said so since 2026-08-23. They were not complete: **the shipped band structure had no Dirac point.** `k_path_graphene()` wrote the zone corner as `4/(3√3)` -- a factor of π too small for the bond-length convention `graphene_hamiltonian()` actually uses -- three lines above an M point written *with* π factors, so `|φ(K)| = 2.5718` instead of 0 and `band_structure.png` showed a **14.40 eV gap at its own K label**, minimum 11.20 eV over the whole path, for five weeks of daily runs. Corrected: gap at K = 2.5e-15 eV, gap at M = 5.6000 eV = `2t` exactly. The superseded path is retained as `HIGH_SYMMETRY_LEGACY` and reproducible via `calculate_band_structure(legacy=True)`, per the annotate-don't-delete rule. **Why five weeks of automated audits never touched it is the more useful finding, and Section 2.10 states it:** three band/DOS figures existed and not one was a calculation that *could* have failed -- `band_structure_simple.png` draws `E = ±v_F|k|`, the known answer; `calculate_density_of_states()` returns `|E|/(πt²)` analytically and never calls the Hamiltonian at all, demonstrated by mutation (replacing H with the k-independent gapped matrix `diag(+7,−7)` leaves its output **bitwise identical**); and the third computed and was wrong. A figure that hard-codes the expected answer certifies whatever sits beside it. Five exact validations now hold, two of them after their first form failed and with both failures kept in the source: `|φ|` = 3, 1, 0 at Γ, M, K; `v_F` = 3`t`a_cc/2ℏ = 9.0627e5 m/s to 1e-6, **after** the first form mistook trigonal warping for a convergence failure and was rewritten to check the exactly derived leading coefficient `q/4` (measured 0.250013); exact particle-hole symmetry bitwise at 4000/4000 random k; and a DOS computed from the bands reproducing the exact slope `2q_e²/(πℏ²v_F²)` to 3.6% and the van Hove peak at `E = t` to 1.0% -- whose own first normalisation check ("integrates to 4 states per cell") is **vacuous**, passed to 16 digits with the reciprocal vectors wrong, and is retained *labelled as vacuous* beside the two checks that replaced it. **Genuinely new physical content:** Section 2.7 bounds the linear-dispersion approximation every device chapter rests on for the first time -- 1% over 0.1-0.5 eV, 5% to 1.0 eV, stated as a window rather than a limit because the deviation is non-monotone at small `E_F` for a k-sampling reason the audit prints -- and gives its **sign**: the full bands hold *more* states than the cone above ~0.25 eV, so a linear-DOS device model **understates** carrier density at high bias, erring conservatively for drive current and against itself for quantum capacitance. **Nothing in Chapters 4-7 moves**, and Section 2.11 tabulates why: those chapters import `v_F`, the linear DOS and the absence of a gap, and import the k-path and the 2×2 `H(k)` **nowhere** -- which is simultaneously why the error survived and why no device number depends on it. |
| 3 | Transport and optical properties (mobility, minimum conductivity, quantum Hall effect, universal optical absorption, Pauli blocking) | **DRAFTED 2026-09-29** (`thesis/03-transport-and-optical-properties.md`, ten sections), and as with Chapter 2 one day earlier, the drafting falsified this row's previous status, which read "Computational results complete (`graphene_transport_properties.py`, optical absorption scripts)" and had said so since 2026-08-23. **Four defects, all in that module, all caught by comparison against an exactly known value and none by a plausible range** (`graphene_transport_optical_audit.py`, `transport_optical_audit_output.txt`, Section 3.9): (T1) the universal optical conductivity was coded `pi*e^2/(2*hbar)` while its own docstring states `pi*e^2/(2*h)` -- a factor 2*pi, implying a single-layer absorption of **14.4044%** instead of 2.2925%; (T2) the conductance quantum was coded `2e^2/hbar` (4.8683e-04 S against 7.7481e-05 S) and the Dirac-point minimum conductivity `4e^2/hbar`, commented "a hallmark of graphene" -- the hallmark is `4e^2/(pi h)`, and the shipped value was **19.74x** it, wrong by 2*pi *and* missing the ballistic `1/pi`; (T3) the Pauli-blocking switch compared `hbar*omega` in JOULES against a bare `0.2` intended as eV, so the branch required a photon energy above 1.25e18 eV and **was never taken at any frequency** -- swept from 0.1 nm to 100 um the function returned exactly ONE distinct value, a constant `0.5*sigma_0`, so the physics its comment describes had never been computed; (T4) the acoustic-phonon mobility used the three-dimensional `T^-3/2` deformation-potential exponent, where graphene's `rho_LA` is linear in T and density-independent so `mu ~ 1/T` [Hwang & Das Sarma, PRB 77, 115449 (2008)] -- measured exponent of the shipped curve 1.5000, understating the 500 K mobility by 1.291x. All four corrected in place with the superseded expressions and their numbers kept in comments; the audit now reads values FROM the module rather than re-typing them, so its seven checks are regression guards. **One check was itself wrong and is recorded rather than retuned:** the `sigma_0/(eps_0 c)` == `pi*alpha` anchor first ran at tol=1e-12 and FAILED at 4.5e-12, which is the slack between CODATA's *measured* alpha and the alpha implied by its own e, h, eps_0 and c, inside the 1.6e-10 uncertainty CODATA quotes; the fix was to check the algebraic identity (3.0e-16) and name the slack, not to loosen a number until it passed. **Blast radius checked rather than assumed:** none of this module's figures are committed and nothing imports it -- but `graphene_photodetector_model.py`'s header CITES `calculate_optical_conductivity` as the source of its 2.3% figure while independently re-typing the correct `pi/137.036`, so **Chapter 6 was protected from a 6.28x error in its first-line input only because it did not use the result it cites.** Same shape as Chapter 2: what was wrong was what nothing read. **Genuinely new physical content:** Section 3.6 computes the gate-tunable Pauli edge `hbar*omega = 2 E_F` for the first time in this repository -- 6199 nm at E_F = 0.1 eV down to **1550 nm at E_F = 0.4 eV**, i.e. the telecom C-band lands inside Chapter 4's existing gate range, which is the physical basis of graphene electro-absorption modulators -- and Section 3.3 separates the theoretical `4e^2/(pi h)` = 20.27 kOhm/sq floor from the measured `~4e^2/h` = 6.45 kOhm/sq one instead of conflating them, the latter being the direct origin of Chapter 4's single-digit on/off ratio and of Chapter 7's "no digital logic" verdict. **Two consecutive chapters, drafted on consecutive days, whose "results complete" label was falsified by the act of writing them** -- the 2026-09-28 finding that PROSE IS A DETECTOR, replicated. |
| 4 | Device physics: contact resistance, quantum capacitance, GFET transfer characteristics, RF figures of merit | In progress -- DC model, RF small-signal model (incl. access resistance + extrinsic pad-capacitance estimate), contact-resistance-vs-channel-length crossover analysis, spatially-resolved work-function-dependent contact doping model, and a DFT-Fermi-shift-derived edge-vs-top contact geometry model implemented (`graphene_fet_model.py`, `rf_small_signal_model.py`, `contact_resistance_crossover.py`, `graphene_contact_doping_model.py`, `graphene_edge_contact_model.py`); literature review of contact resistance, quantum capacitance, f_T/f_max (incl. the 2026-08-26 access-resistance/pad-capacitance revision), contact doping, and contact geometry (2026-09-05) complete; draft written (`thesis/04-graphene-fet-device-physics.md`); per-metal Rc recalibration against 4 literature-sourced values (Cu, Ni, Au, Pd) attempted 2026-08-31 -- additive decomposition found NOT to hold for 3 of 4 metals (Section 4.7), an open modeling question rather than a completed recalibration; Ti/Cr not yet recalibrated (no literature Rc sourced); edge-vs-top geometry model (Section 4.8) reproduces the qualitative direction and large-hole-diameter trend of Passi et al.'s patterned-contact data but not their small-diameter upturn or the full measured ~11x device-level reduction **Conditioned 2026-09-22 (Section 4.10):** the Section 4.7 recalibration table audited for near-cancellation. `kappa` and sign-robustness run in *opposite* order -- Pd's positive residual, the only row consistent with the additive decomposition surviving, flips sign on a 9.6% error in `R_extra`, while Cu's large negative needs 81%. The additive decomposition therefore fails **robustly**, and amplified input error is eliminated as the cause of the negatives (one branch of an item open since 2026-08-31, narrowed by measurement for the first time). |
| 5 | Graphene interconnects: resistivity vs. linewidth, comparison to scaled copper | In progress -- literature review and resistivity-vs-linewidth model complete (`notes/2026-08-23-interconnect-resistivity-vs-linewidth.md`, `graphene_interconnect_model.py`); copper baseline includes a literature-calibrated liner/barrier-thickness effect and a parallel-conduction refinement, plus thin-Ru/Co-liner scenarios (`notes/2026-08-28-copper-liner-barrier-thickness-effect.md`, `notes/2026-08-29-liner-parallel-conduction-and-thin-liner-scenarios.md`) -- the realistic-edge-quality crossover conclusion now ranges from "no crossover" (bare 3nm TaN/Co, or thin-Ru-liner) to "~3 nm" (thin-Co-liner), i.e. materially liner-choice-dependent rather than a single number; draft written and updated (`thesis/05-graphene-interconnects.md`); still to add: graphene-all-around-metal liner/cap model (a distinct architecture from the Cu-liner corrections above), sourcing the specific liner-material resistivity values against a quantitative reference **Conditioned 2026-09-22 (Section 5.5):** the resistivity formula is reinforcing at every width (`kappa` <= 0.66), but the `lambda_impurity` calibration *behind* it is a three-term subtraction with `kappa = 19.71` -- **the worst-conditioned step in this thesis**, larger than anything in Chapter 4 (11.46) or Chapter 6 (4.55). It propagates into `rho(W)` damped but not removed (`kappa` 0.88 -> 1.55 over 18-52 nm). At the calibration width three sensitivities are **exactly zero** and `S(rho_bulk)` changes sign through it. The liner model's `W_eff = W - 2t` carries a previously unstated ~20% bar on `rho_eff` at 18 nm. |
| 6 | Graphene photodetectors: responsivity, gate-tunable absorption | In progress -- literature review complete (`notes/2026-08-24-photodetector-responsivity.md`, `notes/2026-09-06-plasmonic-enhancement-graphene-photodetectors.md`, `notes/2026-09-07-spatial-photocarrier-collection-model.md`); quantitative responsivity/gain-bandwidth model implemented (`graphene_photodetector_model.py`, `photodetector_responsivity_gain_tradeoff.png`); literature-calibrated plasmonic near-field absorption-enhancement model added (`graphene_plasmonic_photodetector_model.py`, `plasmonic_photodetector_enhancement.png`: 25x/8.5x on-resonance enhancement for two literature designs); spatially resolved, metal-dependent photocarrier collection model added (`graphene_photodetector_collection_model.py`, `photodetector_collection_efficiency_by_metal.png`: 1.00x-1.47x collection-enhancement factor across this thesis's 7-metal table, single-reinforcing-contact scope), integrating the Chapter 4 Section 4.5 contact-doping-profile machinery into the photodetector model as planned since 2026-08-24; self-consistent TWO-CONTACT collection model added 2026-09-17 (`graphene_photodetector_two_contact_model.py`, `photodetector_two_contact_net_response.png`, `notes/2026-09-17-two-contact-self-consistent-collection.md`), which **reverses** the single-contact model's metal ranking for a symmetric device (Pt goes from best to worst, overstated 5.5x) and reproduces the literature's symmetric-cancellation result to machine precision; signed, carrier-resolved two-contact model added 2026-09-18 (`graphene_photodetector_signed_carrier_model.py`, `photodetector_signed_carrier_response.png`, `notes/2026-09-18-signed-carrier-resolved-contact-fields.md`), which places the n/p crossover at its physical 5.4 eV rather than graphene's 4.5 eV [Giovannetti et al., PRL 101, 026803 (2008)] and tracks electrons and holes separately: it **confirms and doubles** the asymmetric-pair result (best pair inverts from Cr/Pt at |N| = 0.917 to Ti/Pt at |N| = 1.832, insensitive to the crossover) and **contradicts** the 2026-09-17 symmetric ranking, which turns out to be crossover-dependent and is no longer claimed -- all three of its validations are against exactly known values, two passing bitwise; draft written (`thesis/06-graphene-photodetectors.md`, Sections 6.1-6.9); still to add: non-uniform illumination, bolometric-mechanism model, and a reconciliation of the 0.12 eV measured contact potential step with the 0.25-1.07 eV offsets the metal table assumes; **Section 6.9 added 2026-09-19** removes the uniform-illumination simplification standing since Section 6.6, showing illumination enters only as a weight g(x) on a fixed collection kernel k(x), and produces the chapter's first *bound* rather than a ranking: |N[g]| <= max|k| for every illumination pattern at fixed photon number. Its practical verdict is that asymmetric metallisation beats illumination engineering by ~4x per incident photon (~8x against the shadow mask Shimomura et al. actually built), so masks are worth using only on symmetric devices. **Section 6.11 added 2026-09-20** removes the linear dW -> doping assumption standing since Section 6.6, implementing Khomyakov et al.'s (PRB 79, 195425 (2009)) self-consistent square-root relation (`graphene_contact_doping_nonlinear_model.py`, `contact_doping_nonlinear_relation.png`, `notes/2026-09-20-nonlinear-work-function-to-doping-relation.md`). It **confirms** Section 6.8's Ti/Pt headline (1.832 -> 1.742, a 4.9% compression, stable over alpha = 0-5 eV^-1) but **retracts two claims of Section 6.9**: the 'max|k|/|N_uniform| < 1.02 for every asymmetric pair' bound (14 of 21 pairs violate it under 6.9's *own* linear model, worst Au/Pd at 22.4x) and the conclusion that masks are for symmetric devices only (a perfect mask gains 8.4x on Au/Pd per incident photon). The design recommendation survives on different grounds: 8.4x of a small number is still 4.9x worse than unmasked Ti/Pt. It also establishes that **four of this thesis's seven metals -- Ti, Ni, Pd (chemisorbed, d_eq < d0) and Cr (untabulated) -- lie outside the regime of the relation the chapter depends on**, Ti most consequentially, so a description of the chemisorbed metals that does not go through their work function is now the chapter's central open problem  **Section 6.12 added 2026-09-21** removes the last of the chapter's three hidden conventions -- the single 5.4 eV p/n crossover shared by all seven metals (`graphene_per_metal_crossover_model.py`, `per_metal_crossover.png`, `notes/2026-09-21-per-metal-crossover-from-the-chemical-interface-term.md`). Khomyakov et al. put the crossover at w_cross(d) = W_G + Delta_c(d), so 5.4 eV is that expression at the *physisorbed* separation, not a constant. Four exact validations (49/49 ordered pairs bitwise against the Section 6.8 model; 49/49 bitwise in the ell -> infinity limit; the anchor bitwise for 61/61 values of ell; charge conjugation exactly zero over 27 cases). Results: no metal changes doping type, but **Cu's dW moves 17.12%** -- more than three times Section 6.11's 4.9% compression, **falsifying this session's own pre-registered <=10% prediction**, which had been calibrated on that precedent; and the same-sign pair Cu/Au moves **77.9%** against 1.06-1.21% for every straddling pair, a **65-fold difference from identical inputs** that independently reconfirms Section 6.11's ratio-amplification mechanism under a perturbation of a different kind. The anchored exponential also **diverges** on the chemisorbed metals (w_cross of 6.25-62.6 eV, above the highest elemental work function anywhere), refusing Ti, Ni and Pd for a reason wholly unrelated to Section 6.11's d_eq < d0 -- two distinct parts of one framework failing on the same three metals. Net effect on the chapter: the Cu/Pt straddling design rule is robust to 1.2%, the Ti/Pt headline remains **untested** by two successive refinements rather than confirmed by them, and Sections 6.8-6.11 now carry an explicit +-17% (per-metal) / +-78% (same-sign pair) error bar rather than an implicit claim of exactness **Section 6.13 added 2026-09-22** (the second of two sessions that fired the same day; Sections 4.10/5.5 above are the first, and the two are independent) asks the complementary question Section 6.12 could not: what if 5.4 eV is simply wrong by some amount? The per-metal crossover cannot answer it, because it needs d_eq -- available for 3 of 7 metals -- and all five headline pairs of Section 6.11.2's retraction table contain a metal outside that set. A scalar offset delta on w_cross has the wrong shape but reaches all seven (`graphene_crossover_sensitivity_model.py`, `crossover_sensitivity.png`, `notes/2026-09-22-how-far-can-the-crossover-move.md`), and is the literal meaning of the '~' both source papers put on 5.4 eV. Five exact validations, one of them found mid-session and decisive: a straddling pair's <|dW|> equals (W_B - W_A)/2 EXACTLY, independent of w_cross, so **Section 6.12's 65x amplification ratio is not merely large for such pairs but UNDEFINED** -- it divides by an exactly zero input change. Restated in the well-defined statistic S = dln|N|/d(delta) the dichotomy survives on all 21 pairs at **30x**, down from 65x on three. The chapter's strongest robustness result also lands here: writing dW = m -/+ s, a scalar offset moves m and leaves s = (W_B - W_A)/2 exactly fixed, and s multiplies the symmetric part of E(x) that sets sign(N) -- so **no scalar offset whatever can reverse a response sign**, confirmed by 12621 evaluations over +-3 eV with zero sign changes. Ti/Pt moves 0.45% over the nominal band, its third successive survival. Retracted: Section 6.11.2's *count* of five >2x mask-gain pairs, which runs 3-4-5-5-6 across the band (the table's numbers stand at delta = 0; it is annotated in place). The session's methodological finding is a THIRD failure class: the first sign-flip implementation bisected without a bracket test and printed 21 plausible fake roots while all five exact validations passed, because they validate the model and the fault was in the analysis layer **Section 6.14 added 2026-09-23** closes the blind spot Section 6.13 named in its own closing remark. 6.13's robustness theorem covers only the COMMON MODE of the crossover uncertainty; Section 6.12 predicts a crossover PER METAL, so the two contacts of a pair do not share one, and their DIFFERENCE is the only part that can move s (`graphene_differential_crossover_model.py`, `differential_crossover.png`, `notes/2026-09-23-the-differential-crossover-and-the-sign-flip-threshold.md`). Result in closed form: **the differential crossover offset needed to reverse a pair's photoresponse is EXACTLY that pair's work-function gap, tau* = W_B - W_A**, independent of the common offset, of lambda, of L and of every transport parameter -- so Section 6.13's entire scalar uncertainty band is exactly irrelevant to the sign. Six exact validations (two bitwise against the 6.8 and 6.13 models; N = 0 at tau* to 7.4e-15 over 189 cases; charge conjugation exactly zero 315/315; a BRACKETED bisection landing on the closed form to 4.8e-15 eV and drifting 9.7e-15 eV across +-1 eV of common offset). The sixth was unplanned and subsumes three earlier results: **N is exactly ODD in the differential offset s and exactly EVEN in the common offset m** (2.3e-15), which yields tau* as a corollary, upgrades 6.13.5's 12621-point sign scan from evidence-over-an-interval to a parity valid without one, and **excludes** two of this session's own pre-registered predictions rather than merely missing them -- odd-in-s forces |N| to be symmetric about the flip. The two channels are qualitatively unlike: the common one has no sign threshold at any magnitude, the differential one has a finite threshold for every pair, worst margin Cu/Au at 2.85x the model's largest estimate. And the chapter's four smallest margins -- Au/Pd 0.020 eV, Ni/Au 0.060, Ni/Pd 0.080, Cr/Cu 0.150 -- belong to exactly the four pairs Section 6.12 REFUSES, because chemisorbed metals sit ~1.2 A below the anchor, which is both why 6.12.4's exponential diverges on them and why their crossovers should differ most: **the model's reach and the thesis's risk are anti-correlated by construction**. Au/Pd is stated as a threshold no available model can reach, not as a prediction of reversal. The 2026-09-22 analysis-layer audit item is also discharged here: the new root-finder refuses an unbracketed interval and fired on 21/21 deliberately root-free intervals where last session's loop returned a plausible fake root on every one |
| 7 | Discussion and outlook: graphene's realistic near-term application space in the semiconductor industry | **DRAFTED 2026-09-27** (`thesis/07-discussion-and-outlook.md`, eleven sections). All ten synthesis threads listed below are discharged; the record of how each was handed to this chapter is kept in place rather than deleted, per the repo's annotate-don't-rewrite rule. The draft's four arguments: (1) Chapters 4-5 optimise a SINGLE interface while Chapter 6's figure of merit is a DIFFERENCE between two, which is a structural split rather than a disagreement and yields the device-facing statement that a differential contact-chemistry error leaves signal magnitude predictable and POLARITY unpredictable; (2) a bound or parity is a different kind of knowledge from a ranking over a tabulated set, and Section 6.9's retracted '<1.02 ceiling' shows the grammatical FORM of a claim is not evidence about the population it was checked against; (3) five failure classes with five non-overlapping detectors, so verification effort must be allocated BY CLASS -- with the recurring absolute-tolerance class now at three consecutive sessions and, on 2026-09-27, three instances inside the session auditing it, supporting the structural conclusion that exhortation does not fix it and only reporting a measured value beside a DERIVED bound does; (4) the chapter that has publicly overturned its own headlines most often is not the one whose arithmetic is most fragile but the one that has been LOOKED AT most -- kappa = 19.71 in Chapter 5 against 11.46 in Chapter 4 and 4.55 in Chapter 6, and the only convergence blind spots found anywhere in the repository are at that same Chapter 5 calibration. Application verdicts: digital logic NO (no gap, not fixable without losing the mobility); RF/analogue and interconnect CONDITIONAL (the interconnect crossover ranges from 'none' to ~3 nm on liner choice alone); photodetector the strongest near-term case; contacts the enabling contribution. Two threads added on the drafting day itself: an untested prediction that any DIFFERENTIAL interconnect figure of merit is dominated by edge-scattering variance rather than mean resistivity, and the caution that a bound derived from a near-cancelling quantity inherits the near-cancellation. **Historical record of the ten threads follows.** Synthesis task **restated 2026-09-18**. As framed on 2026-09-17 (Chapter 4's low-contact-resistance metals Ni/Au/Pd being the worst for photoresponse) the contradiction does not survive Section 6.8's signed model: at the physical n/p crossover Pd sits near the top of the symmetric ranking, and the symmetric ranking itself is no longer claimed. The durable and more interesting tension is narrower: Chapter 4 optimises a **single** junction, whereas the two-terminal photoresponse depends on the **difference between two** junctions (Section 6.8.4: |N| = 1.832 for Ti/Pt vs at best 0.556 for any symmetric pair) and is to first order blind to either junction's own quality. A second, now sharper thread for this chapter: the symmetric-device design rule (pick the metal that dopes graphene least) and the asymmetric one (maximise |dW_A - dW_B|) point in opposite directions, and Chapter 5's interconnect argument shares Chapter 4's single-junction framing. **Third thread, added 2026-09-19:** Section 6.9 turns the design question into a bounded one -- max|k| caps what any illumination pattern can achieve -- while Chapters 4 and 5 still argue from unbounded single-junction optimisation. Whether an analogous ceiling exists for contact resistance and for interconnect resistivity is the sharpest form of the synthesis question this chapter has yet been handed. **Fourth thread, added 2026-09-20:** that third thread is now *itself* the example, because Section 6.11 retracted the ceiling it was built on. Three of Chapter 6's last four sessions have overturned the previous one's headline, always by widening a sample rather than by finding a physics error -- which is the methodological point Chapter 7 should make in its own voice: this thesis's chapters differ not in rigour but in how many cases each claim was checked against, and Chapters 4 and 5 have never had that widening done to them  **Fifth thread, added 2026-09-21:** Section 6.12 supplies the first clean counter-example to that pessimistic reading and sharpens it into something Chapter 7 can actually argue. The session pre-registered four predictions *before* writing the model and committed them to git first; three held and one was falsified at 17.12%, and the falsification was informative rather than embarrassing precisely because the prediction had been stated with a number attached. The methodological claim Chapter 7 should make is therefore not 'this thesis's claims keep failing' but the sharper and more useful one: **a claim's reliability tracks how wide a sample it was checked against and whether its error bar was stated in advance, not how carefully its physics was derived** -- Chapter 6 now has six sections' worth of evidence for that, and Chapters 4 and 5 state general design rules from tabulated subsets with no error bars at all. The concrete synthesis task: Section 6.12 shows a near-cancelling quantity amplifies an input uncertainty 65-fold relative to a reinforcing one. Chapter 4's contact resistance and Chapter 5's liner scenarios both contain differences of comparable quantities, and neither has been asked whether it is a near-cancellation **Sixth thread, added 2026-09-22 -- and the fifth thread's concrete synthesis task is now DONE.** The fifth thread asked whether Chapter 4's contact resistance and Chapter 5's liner scenarios are near-cancellations. They were audited (Sections 4.10, 5.5) and the answer is yes, with two findings Chapter 7 should carry. (i) **The worst-conditioned step in the thesis is not in Chapter 6 but in Chapter 5's calibration**, `kappa = 19.71` against Chapter 6's 4.55 -- so the chapter that has failed most often publicly is not the one whose arithmetic is most fragile; it is the one that has been *looked at* most. That reframes threads four and five from a story about Chapter 6's unreliability into a story about **unequal scrutiny**, which is a claim about method rather than about graphene. (ii) `kappa` alone is the wrong instrument: in Section 4.7 it ran *opposite* to the robustness of the conclusions drawn from it, and the useful quantity was the **sign-flip margin** (4.23). Chapter 7's methodological claim should accordingly be stated as: a result's reliability tracks how wide a sample it was checked against, whether its error bar was stated in advance, **and which of its terms cancel** -- with the honest rider that two of this session's own six pre-registered predictions failed, both by generalising from a partial enumeration, the same failure mode as 2026-09-20 and 2026-09-21. **Seventh thread, added 2026-09-22 (second session of that day, independent of the sixth):** Section 6.13 turns the methodological claim into a taxonomy Chapter 7 can structure itself around, because it produced a failure that pre-registration does not catch. Three classes are now on record with a distinct detector each: (i) an implementation error inside the model, caught by an exactly known value (2026-09-17's inf-inf transit integral); (ii) a claim quantified on an unrepresentative sample, caught only by enumerating the full population (Section 6.9's <1.02 ceiling, 14 of 21 pairs violating); and (iii) an error in the ANALYSIS layer wrapped around a validated model -- an unbracketed bisection that returned its own endpoint 21 times, caught by hand-checking one row against the physics while every exact validation passed. Pre-registration addresses (ii) and nothing else, so verification effort has to be allocated BY FAILURE CLASS rather than uniformly -- which sharpens the sixth thread's finding, since the kappa = 19.7 conditioning result says WHERE Chapters 4 and 5 are fragile and this says WHICH KIND of check would find it. Section 6.13 also supplies the thesis's first PROOF-shaped robustness result rather than a survey one -- sign(N) is immune to any scalar crossover offset because the perturbation cannot touch the symmetric part of the contact field -- which is the right model for what a design rule in Chapters 4 and 5 ought to look like **Eighth thread, added 2026-09-23:** Section 6.14 gives Chapter 7 the first example of the opposite of a retraction -- a result that made three earlier results *stronger* by absorbing them. One parity identity (N odd in the differential offset, even in the common one) implies 6.13.5's sign-invariance scan, implies 6.14's own pre-registered derivation, and excludes two predictions that had been stated in good faith. The methodological point is the counterpart of the seventh thread's taxonomy: the three failure classes say how claims break, and this says how they consolidate -- **a claim checked on a wide sample and a claim derived from a symmetry are not the same kind of knowledge, and this thesis has been accumulating the first while assuming it was acquiring the second**. Chapter 6 now has exactly one symmetry-derived result and six survey-derived ones. The sharper form of the synthesis question: which of Chapters 4 and 5's design rules could be restated as parities or bounds rather than as rankings over a tabulated set? Section 6.14 also supplies a device-facing statement Chapter 7 should carry verbatim, because it is the kind of thing an industrial reader acts on: at zero bias a differential contact-chemistry error of unknown direction leaves the signal MAGNITUDE predictable and the signal POLARITY unpredictable -- the worse of the two failure modes for a differential readout, and invisible to every single-junction figure of merit in Chapters 4 and 5 |

This draft will be expanded chapter-by-chapter as the corresponding
research notes and computational studies are completed; see
`AUTOMATION_LOG.md` for the dated record of what has been added and when.

## A.3 Chapter 2: Electronic structure

<a id="s2-10"></a>

*Moved from Chapter 2, §2.10. Original section number kept.*

### 2.10 Correction of record: this repository's band structure had no Dirac point

**What was wrong.** From the creation of `graphene_band_structure.py` until
2026-09-28, the module's high-symmetry points were written as

```
k_point = [4.0/(3·√3), 0]          ← no π
m_point = [π/(3·√3), π/3]          ← π
```

three lines apart in the same function, while `graphene_hamiltonian()` builds
its phases from nearest-neighbour vectors of unit length — the bond-length
convention of Section 2.2, in which the zone corner is at `4π/(3√3)`. The
shipped K was therefore **exactly a factor of π too small** (the measured
ratio `|K_correct|/|K_shipped|` is 3.141592654), and the shipped M was wrong in
both components.

**What it cost.** `|φ(K_shipped)| = 2.5718` instead of 0. The consequences,
measured:

| quantity | as shipped | corrected | exact value |
|---|---|---|---|
| gap at the point labelled K | **14.4019 eV** | 2.5 × 10⁻¹⁵ eV | 0 |
| minimum gap over the whole Γ–K–M–Γ path | **11.2000 eV** | 2.5 × 10⁻¹⁵ eV | 0 |
| gap at the point labelled M | 11.2000 eV | 5.6000 eV | `2t` = 5.6 eV |

Graphene's single defining electronic property — the gapless linear crossing —
was **absent from `band_structure.png`**, the oldest figure in the repository
and the one this chapter is built on, through five weeks of daily automated
sessions. The superseded values are retained in the module as
`HIGH_SYMMETRY_LEGACY` and the broken spectrum is still reproducible via
`calculate_band_structure(legacy=True)`, per the standing rule that deleting a
superseded implementation destroys the only oracle available for judging its
replacement.

**Why nothing caught it, which is the more useful half.** The repository
contained three band-structure or density-of-states figures and **not one of
them was a calculation that could have failed**:

- `band_structure.png` was computed, and was wrong;
- `band_structure_simple.png` *draws* `E = ±v_F|k|` — the known answer,
  plotted rather than derived;
- `density_of_states.png` came from a function that returns `|E|/(πt²)`
  analytically and **never calls the Hamiltonian at all**. This is demonstrated
  in the audit by mutation: replacing the Hamiltonian with the k-independent
  gapped matrix `diag(+7, −7)` — which has no Dirac point and no dispersion
  whatsoever — leaves the shipped DOS output **bitwise identical**.

A density of states computed from the bands would have caught this on the
first run, because it would have shown no states near zero energy and its van
Hove peaks in the wrong place. The figure that could have falsified the
calculation was instead the figure that certified it. That function is now
annotated in place and retained, labelled as non-discriminating, with a pointer
to `dos_from_bands()`.

**What does and does not move.** Nothing in Chapters 4 through 7 changes. Those
chapters import `v_F` and the linear density of states (Section 2.11) and never
the k-path, which is simultaneously the reason the error was survivable for
five weeks and the reason no device number depends on it. `v_F` itself is
confirmed to 1 × 10⁻⁶ relative by the corrected computation.

**The methodological point, stated narrowly.** The audit that found this also
committed the same fault inside itself, twice, within an hour: its first
normalisation check on the computed DOS — that the density of states integrates
to four states per unit cell — passed to sixteen digits *while the reciprocal
lattice vectors were wrong*, because the normalisation divides by the sample
count and returns four states per cell whatever the Hamiltonian is. It is kept
in the code and printed as **vacuous**, beside the two checks that replaced it
and do have power: periodicity of `|φ|` under the reciprocal vectors, and a
count of the Dirac cones inside the sampled cell (which must be two, and was
one). The lesson this chapter carries forward is therefore not "check your
conventions" but something more specific and more actionable: **a check that
cannot fail is worse than no check, because it certifies whatever sits beside
it**, and the way to tell the two apart is to make the check fail on purpose.

## A.4 Chapter 3: Transport and optical properties

<a id="s3-3-1"></a>

*Moved from Chapter 3, §3.3.1. Original section number kept.*

### 3.3.1 Annotation (2026-10-09): Chapter 4's floor is 2.07× this one, and the two had never been compared

This section states that the measured floor is `4q_e²/h` = 6.45 kΩ/□, and then
says the repository "encodes that picture as `n_puddle = 5 × 10¹⁵ m⁻²` in
`graphene_fet_model.py`". **Those two statements had never been checked against
each other**, in the eleven sessions since this chapter was drafted. They
disagree.

Chapter 4's disorder floor is `σ_min = e·μ·n_puddle` — a one-line calculation,
available since August, decidable only once §4.6.7 established that `n_puddle`
is the **total** carrier density `n_e + n_h` rather than a net or per-species
density. It gives 3.204353 × 10⁻⁴ S/sq = **3.1208 kΩ/□**, against the
1.549618 × 10⁻⁴ S/sq = 6.4532 kΩ/□ in the table above:

> **Chapter 4's channel at neutrality is 2.0678× more conductive than the floor
> this section calls measured.** The `n_puddle` that would reconcile them is
> 2.418 × 10¹¹ cm⁻², against the shipped 5 × 10¹¹ cm⁻².

Nothing in the table above is changed and no number in it is deleted. What is
corrected is the *word* "encodes": `n_puddle = 5 × 10¹⁵ m⁻²` is not an encoding
of this section's 6.45 kΩ/□ floor, it is an independent parameter that happens
to be of the same order and lands a factor of two on the conductive side. The
sentence should be read as describing the same *picture* — charge puddles
averaging over a non-uniform local density — and not the same *number*.

Two consequences recorded rather than acted on:

- **§4.4's on/off ratio is understated** by up to 1.44× at the terminals, if
  this section's measured floor is the right one for that device. §4.6.7
  measures it: 1.2408 → 1.7820 at `V_ds` = 0.1 V. It does not disturb this
  section's central negative result, nor §7.2's conclusion — 1.78 is no more a
  switch than 1.24.
- These two floors are **not the same physics** and are not required to agree.
  `4q_e²/h` is a quantum/ballistic minimum conductivity observed
  experimentally; `e·μ·n_puddle` is a diffusive disorder floor. The defect is
  not that they differ, it is that a thesis quoted both in chapters that feed
  each other without ever putting them side by side — and that the diffusive
  one came out **below** the measured one, which no graphene sheet does.

An external cross-check, for scale and not as an anchor:
[Wiedmann et al., Phys. Rev. B **84**, 115314 (2011), arXiv:1107.3929] extract
both species from the Hall coefficient and report `n = p ≈ 4.2 × 10¹⁴ m⁻²` at
neutrality for their sample B, i.e. a total puddle density of ≈ 8.4 × 10¹⁴ m⁻².
The shipped `n_puddle` is **5.95×** that, and even the Chapter-3-matched value
is 2.88× it. Different sample, unstated substrate quality, so this is context
rather than a calibration — §4.6.7's comparison is internal to this thesis and
is the stronger finding.

Source: `graphene_fet_signed_carrier_model.py` §4,
`fet_signed_carrier_output.txt`.

<a id="s3-3-2"></a>

*Moved from Chapter 3, §3.3.2. Original section number kept.*

### 3.3.2 Annotation (2026-10-29): the floor is a geometric mean, so the matched puddle density is an underestimate

§3.3.1 divides Chapter 3's measured floor by `e·μ` to obtain the puddle density
that would reproduce it, and that division assumes the two branches carry the
same mobility. They do not, in any real device (§1 of
`notes/2026-10-29-the-floor-is-a-geometric-mean.md` tabulates the two competing
mechanisms and the measurements behind each). Minimising the two-carrier Drude
conductivity over the signed net density — derived in §4.11 — gives

```
  σ_min = e · n_p · sqrt(μ_e μ_h)
```

**the geometric mean of the two branch mobilities, not the arithmetic mean**,
and the minimum does not sit at charge neutrality but displaced onto the
low-mobility branch at `n* = -n_p r/sqrt(1-r²)`, `r = (μ_e-μ_h)/(μ_e+μ_h)`.

Consequences for this section, both second order and both small:

1. The divisor in `n_puddle_matching_ch3()` should be `e·sqrt(μ_e μ_h)`, not
   `e·μ`, so the **matched puddle density is an underestimate** by
   `μ̄/sqrt(μ_eμ_h)` — **+0.88 %** at `μ_e/μ_h = 1.3`, +2.06 % at 1.5.
2. The 2.88× and 5.95× ratios against Wiedmann et al. above move by the same
   factor and **no conclusion in §3.3 or §3.3.1 changes**. The numbers above
   are left exactly as they were; this annotation does not revise them.

The reason the effect is this small is worth stating, because it is the same
reason the floor literature is so much cleaner than the doping literature:
`sqrt(1-r²) = 1 - r²/2 + O(r⁴)`, so a branch asymmetry reaches the **floor only
at second order**, while it moves the apparent Dirac point at **first** order —
by a third of a volt at `μ_e/μ_h = 1.3` (§4.11). §3.3's central negative result
— that graphene's off state is a 6.45 kΩ/□ resistor and therefore that graphene
cannot switch off — is robust to the asymmetry for exactly this reason.

Source: `graphene_branch_asymmetric_mobility.py` §3,
`branch_asymmetric_mobility_output.txt`.

<a id="s3-9"></a>

*Moved from Chapter 3, §3.9. Original section number kept.*

### 3.9 Correction of record: four defects in this chapter's computational results

`graphene_transport_properties.py` carried the label "computational results
complete" from 2026-08-23. Drafting this chapter and auditing each claim
against an exactly known value
(`graphene_transport_optical_audit.py`, run recorded in
`transport_optical_audit_output.txt`) found four defects. All are corrected in
place, with the superseded expressions and their numeric values kept in
comments beside the corrections.

**T1 — universal optical conductivity, a factor 2π.** The function's own
docstring states `σ₀ = πe²/(2h)`; the code computed `np.pi * e**2 / (2 * hbar)`.
Shipped value 3.8235 × 10⁻⁴ S against the correct 6.0853 × 10⁻⁵ S, implying a
single-layer absorption of **14.4044%** instead of 2.2925%. *(For the record, so
that no future reader chases it: the 14.40 here and the 14.40 eV spurious gap at
the K label found on 2026-09-28 are numerically a coincidence. That one was a
reciprocal-lattice convention; this one is `ħ` written for `h`. They share no
mechanism.)*

**T2 — the conductance quantum and the minimum conductivity, a factor 2π
twice.** `G₀` was coded `2e²/ħ` (4.8683 × 10⁻⁴ S against 7.7481 × 10⁻⁵ S), and
`G_min` was coded `4e²/ħ` and commented "a hallmark of graphene". The hallmark
is `4e²/(πh)`; the shipped value was **19.74×** it — wrong by 2π *and* missing
the `1/π` of the ballistic self-consistent result.

**T3 — a Pauli-blocking branch that could never be taken.** The code computed
`hbar_omega` in joules and compared it against `2 * mu` where `mu = 0.1` was a
bare number intended as eV. The branch therefore required a photon energy above
0.2 J = 1.25 × 10¹⁸ eV. Swept across the entire electromagnetic spectrum from
0.1 nm to 100 µm, the function returned **exactly one distinct value**: a
constant `0.5 σ₀`. The feature the comment describes — "interband transitions
become significant when photon energy > 2μ" — had never operated in either
direction since it was written. This is the 2026-09-28 "check that cannot fail"
in its purest form: not a check that always passes, but a *branch that is never
reached*, so the code appeared to model the physics of §3.6 while modelling
nothing. Section 3.6's table is the first time that physics has actually been
computed in this repository.

**T4 — a three-dimensional exponent in a two-dimensional material.** Discussed
in §3.2. Measured exponent of the shipped curve: 1.5000, where graphene requires
1.

**Blast radius, checked rather than assumed.** None of this module's figures are
committed to the repository, and no other module imports it. But
`graphene_photodetector_model.py`'s header **cites**
`graphene_transport_properties.py::calculate_optical_conductivity` as the source
of its 2.3% absorption figure, while independently re-typing the correct
`ALPHA_ABS = π/137.036`. Chapter 6 was protected from a 6.28× error in its
first-line input only by the fact that it did not use the result it cites.

**What this makes two of.** Two consecutive chapters, drafted on consecutive
days, whose "computational results complete" label was falsified by the act of
writing the chapter, after five weeks in which daily audits, exact validations,
convergence studies and pre-registrations all passed over the same code. The
2026-09-28 entry proposed that **prose is a detector** and that it is the one
this repository lacked. Chapter 3 is the replication. The mechanism is the same
in both cases and worth stating precisely: an audit checks what the code says
about itself, and a chapter has to state what the code says *about graphene* —
so writing `σ₀ = πe²/(2h)` in a sentence and then reading the line that computes
it is a comparison no test in this repository was making.

---

#### Follow-on items generated by this chapter

- **The §3.6 Pauli edge is not in Chapter 6's model.** A photodetector biased to
  high `E_F` for low contact resistance is switching off its own interband
  absorption. Chapter 4 favours high-work-function metals for exactly that
  reason, so this is a *third* leg of the contradiction Chapter 7 is already
  holding between Chapter 4's contact-metal ranking and §6.7's photoresponse
  ranking — and unlike those two it is not a modelling disagreement but a
  physical mechanism currently absent from the model.
- **`σ_min` enters Chapter 4 through `n_puddle`, not through `4e²/h`.** Whether
  the puddle density Chapter 4 assumes actually reproduces the measured
  6.45 kΩ/□ floor has never been checked, and it is a one-line calculation.
- **The remote-polar-phonon cap on SiO₂ is asserted here, not computed.** The
  ~40 000 cm²/V·s figure is from the literature; the module models only
  acoustic-phonon and impurity scattering, so the third row of §3.2's table has
  no computational backing in this repository.
- **A finite-temperature optical conductivity.** §3.6's Pauli edge is the
  `T = 0` step function; at 300 K it is smeared over several `k_BT`, which
  matters for a modulator's extinction ratio and is not modelled.
- **`n(E_F, T)` at finite temperature** — carried over unresolved from §2.7, and
  §3.6's table inherits the same `T = 0` assumption.

## A.5 Chapter 4: GFET device physics

<a id="s4-6-note-a"></a>

*Annotation moved from the body of Chapter 4, §4.6.*

> **Operating point, made explicit 2026-09-26.** The ≈20 GHz above is
> evaluated at **V_ds = 0.1 V** — the value `plot_fT_fmax()` passes, and
> therefore the value behind `rf_figures_of_merit.png`. Measured precisely:
> peak f_T = 20.279 GHz, peak f_max = 18.731 GHz. This sentence previously
> quoted the number without naming the bias, which invited exactly the
> confusion described in Section 4.6.1: `output_conductance`'s own signature
> default is V_ds = 0.05 V, at which the same model gives peak
> f_T = 10.140 GHz and peak f_max = 9.355 GHz. Neither number is wrong and
> neither supersedes the other; they are one model at two bias points, and
> the ≈20 GHz figure quoted throughout this chapter is the V_ds = 0.1 V
> one.

<a id="s4-6-note-b"></a>

*Correction moved from the body of Chapter 4, §4.6.*

> **Correction, 2026-10-03: the literature-scale f_max numbers are revised
> downward; the normalized-device numbers above are NOT.** The
> distributed gate resistance `gate_resistance()` took its channel length and
> gate width as *function-parameter defaults* read out of
> `graphene_fet_model` at definition time. `compute_fT_fmax`'s
> literature-scale option rescales the device to the 40 µm / 8-finger
> geometry by *rebinding* `graphene_fet_model.W`, and R_g kept the 1 µm it
> had captured — so the gate resistance entering both f_max denominators at
> that geometry was too small by exactly W_RF/W = 40 (0.208 Ω instead of
> 8.33 Ω). The three numbers that change, with the superseded values kept:
>
> | quantity, 40 µm / 8 fingers, V_ds = 0.1 V | superseded | corrected |
> |---|---|---|
> | intrinsic peak f_max | 18.786 GHz | **13.094 GHz** |
> | extrinsic peak f_max (+15 fF/pad) | 3.1885 GHz | **2.2196 GHz** |
> | f_max / f_T at that geometry | 0.928 | **0.647** |
>
> Peak f_T is **bitwise** unchanged at every bias and both geometries,
> because R_g does not enter f_T = g_m/(2π C_gs) at all; and the
> **normalized 1 µm device** — the path behind the ≈20 GHz and 18.731 GHz
> figures quoted above and behind panels 1 and 2 of
> `rf_figures_of_merit.png` — is bitwise unchanged too, verified against the
> pinned pre-correction source rather than assumed. Only the literature-scale
> comparison moves.
>
> **Two consequences for the argument of this section, in opposite
> directions.** First, the extrinsic *degradation* claim below (f_T and f_max
> falling to roughly 15–20 % of their intrinsic values) survives unchanged:
> measured, it moves from 0.1697 to 0.1695. So does the intrinsic/extrinsic
> peak-f_max ratio, 5.892 → 5.899, a 0.13 % change. The error was
> multiplicative in both terms of every *ratio* this section uses, and
> divided out of all of them — which is exactly why it survived five weeks of
> audits that compare ratios against the literature. Second, and less
> comfortably, the corrected model puts **f_max/f_T = 0.65 rather than 0.93**
> at the literature-scale device, i.e. *further* from the f_max/f_T = 1.3–1.4
> that Feijoo *et al.* report for de-embedded devices, not closer. The
> discussion immediately below reads as though the model and that literature
> had been brought into rough agreement; with R_g correct, they have not
> been. The gap is now attributable rather than hidden: at 8 fingers this
> model's R_g = 8.33 Ω is already in the engineered-low range Feijoo *et al.*
> describe, so the remaining shortfall sits in g_ds and in the R_g·C_gd
> feedback term, not in the gate resistance. ~~Quantifying that split is
> left open and is recorded as such in `AUTOMATION_LOG.md`.~~
>
> **Quantified 2026-10-04 in Section 4.6.2, and the sentence immediately
> above is half wrong.** The R_g·C_gd feedback term carries **0.2154 %** of
> the f_max denominator, so it is not a contributing cause at all; and R_g
> is not the binding constraint either, since f_max/f_T is still 0.935276
> in the limit R_g = 0 exactly. The shortfall is g_ds, and g_ds is large
> because the DC model has no current saturation. See Section 4.6.2.
>
> See `graphene_cross_module_delivery_audit.py` (Sections 3, 4 and 7) and
> `notes/2026-10-03-a-frozen-default-in-another-module.md`.

<a id="s4-6-1"></a>

*Moved from Chapter 4, §4.6.1. Original section number kept.*

### 4.6.1 Numerical provenance of g_ds, and a caution about auditing a discretised model

The output conductance g_ds entering the f_max expression above is not a
closed-form derivative. It is a central difference in V_ds of
`transfer_characteristic()`, which is itself a 50-point quadrature of the
local channel resistance over the drain-bias drop
(`V_channel_profile = linspace(0, V_ds, 50)`). Two numerical defaults are
therefore entangled in every f_max number in this chapter: a step size
ΔV_ds = 10⁻³ V and a resolution n_segments = 50. They are not independent,
because the quadrature grid is itself set by V_ds — differencing in V_ds
differences the quadrature error too.

That entanglement was audited on 2026-09-26
(`graphene_gds_quadrature_audit.py`,
`notes/2026-09-26-gds-step-quadrature-audit.md`) because the analogous
absolute step in Chapter 6 (`H_DIFF = 10⁻³`) had been convicted eight days
earlier of moving eleven of Section 6.13's sensitivities by up to 44%.
**Both defaults here are innocent, and the numbers in this section stand
unchanged:**

| quantity | measured effect on g_ds | effect on peak f_max |
|---|---|---|
| step ΔV_ds = 10⁻³ V vs. anchored-plateau step | 1.25 × 10⁻⁶ | — |
| n_segments = 50 vs. n → ∞ (Richardson in 1/(n−1)) | 8.4 × 10⁻⁷ | — |
| both, propagated | — | 1.1 × 10⁻⁷ (0.00001%) |

f_T contains no g_ds and is untouched; the 400-point V_g grid behind
g_m = ∂I_d/∂V_g moves peak f_T by 0.001% under sixteenfold refinement.

The methodological point survives the null result, and is the reason this
subsection exists rather than a one-line footnote. **A step-refinement
study of a discretised function converges to the derivative of the
discretisation that was held fixed, not to the derivative of the model.**
The anchored step criterion adopted in Section 6.13 — accept a step only if
the derivative is unchanged at h/10 *and* h/100 — is therefore blind to a
quadrature bias by construction, at any tolerance. This was verified rather
than argued: the n = 50 bias measured 3.04 × 10⁻⁷ at ΔV_ds = 10⁻³ and
3.04 × 10⁻⁷ at ΔV_ds = 10⁻⁷, drifting 2.3 × 10⁻⁴ across four decades of
step refinement. And the mechanism is not small in general: in a test case
with R(V_ch) = a + b·V_ch², where both the discretised and the continuum
channel average are exact closed forms (the discrete mean of V_ch² over an
endpoint-inclusive grid is V²(2n−1)/(6(n−1)) against a continuum V²/3, so the
bias is exactly V²/(6(n−1))), the same mechanism reaches **1.32%**. Graphene's
channel resistance simply has too little curvature in V_ch across a 50 mV
drop for it to bite here. That is a property of this model, not a general
licence to ignore the effect — a higher-V_ds or a shorter-channel model with
stronger pinch-off curvature would not inherit this exoneration.

Two defects in the surrounding machinery were found by the same audit and
are recorded here because they bear on how the chapter's numbers should be
read. First, the guard `max(V_ds − ΔV_ds, 10⁻⁴)` in `output_conductance()`
moved the evaluation interval without changing the divisor, returning the
true secant slope times exactly (V_ds + ΔV_ds − 10⁻⁴)/(2ΔV_ds) — an 8.42%
silent under-report whenever it fired. It never fired at either operating
point in use (it is reachable only for V_ds ≤ 1.1 mV), so no number in this
chapter was ever affected; it is fixed, with the unaffected path verified
bitwise unchanged. Second, the 2026-09-25 default-scale census scored
ΔV_ds/V_ds as 2 × 10⁻² by reading the function's *signature* default, whereas
at the call site that actually produces this chapter's figures it is
1 × 10⁻². **A default-scale census must read call sites, not signatures**;
the discrepancy is the ratio between the two, here a factor of two and in
general unbounded.

<a id="s4-6-2"></a>

*Moved from Chapter 4, §4.6.2. Original section number kept.*

### 4.6.2 The f_max shortfall, decomposed: one of the two named causes carries 0.22 %, and R_g was never the lever

The annotation above left the residual f_max/f_T shortfall attributed to
"g_ds and the R_g·C_gd feedback term". That attribution was written without
measuring it. `graphene_fmax_shortfall_decomposition.py` measures it, and
the measurement is separable exactly rather than perturbatively, because
the f_max denominator is a **sum**:

```
  f_max/f_T = 1 / (2 sqrt(denom)),   denom = g_ds(R_g+R_s) + 2π f_T C_gd R_g
                                             \_____ A ____/   \_____ B ____/
```

At the literature-scale 40 µm / 8-finger geometry and V_ds = 0.1 V:

| term | value | share of denominator |
|---|---|---|
| A = g_ds(R_g + R_s) | 6.033531×10⁻¹ | **99.7846 %** |
| — of which g_ds·R_g | 3.175543×10⁻¹ | 52.5182 % |
| — of which g_ds·R_s | 2.857988×10⁻¹ | 47.2664 % |
| B = 2π f_T C_gd R_g | 1.302534×10⁻³ | **0.2154 %** |

**Result 1: the feedback term is not a cause.** Deleting C_gd *entirely* —
not improving it, deleting it — moves f_max/f_T from 0.643007 to 0.643701,
an improvement of **0.108 %** against a shortfall of a factor of 2.2. Its
share rises with drain bias but only to 1.99 % at V_ds = 1 V. Section 4.6's
annotation named two residual causes and one of them is not one. The
superseded sentence is struck through in place above rather than rewritten,
per this thesis's standing practice.

**Result 2, which contradicts a thread this chapter has been pursuing for
weeks: R_g is not the binding constraint, even in the limit R_g = 0.** The
same module evaluates each counterfactual through the same code path with
one input deleted:

| counterfactual | f_max/f_T | vs. baseline |
|---|---|---|
| model as committed | 0.643007 | — |
| C_gd = 0 | 0.643701 | +0.108 % |
| 16 fingers (R_g/4) | 0.827025 | +28.618 % |
| 64 fingers (R_g/64) | 0.927229 | +44.202 % |
| **R_g = 0 (perfect gate)** | **0.935276** | +45.453 % |
| R_s = 0 (perfect contacts) | 0.885467 | +37.707 % |
| g_ds ÷ 4.8483 | 1.410000 | +119.282 % |

A zero-resistance gate still lands **28.06 % below the bottom of the
1.3–1.4 band**, and 64 fingers is already within 1 % of that limit, so the
multi-finger mechanism is exhausted at about 16 fingers. The 2026-10-03
session found a real 40× delivery bug in R_g and was right to; the belief
that R_g was the lever on f_max did not survive being measured.

**Why term A is so large: the DC model has no current saturation.**
`transfer_characteristic()` computes I_d = V_ds/(R_channel(V_g,V_ds) + R_c),
a bias-dependent *resistor* — no velocity saturation, no pinch-off, no
drain-field cutoff. Its output conductance is therefore the channel
conductance itself. Measured, g_ds·R_total lies in [1.0026, 1.0487] across
V_ds = 0.05–1 V, i.e. g_ds ≈ 1/R_total to between 0.3 % and 4.9 %. Term A
then collapses to a pure resistance ratio and

```
  f_max/f_T  ≈  (1/2) sqrt( R_total / (R_g + R_s) )
```

which reproduces the exact value to 0.36 % (0.645352 vs 0.643007). Reaching
1.41 requires R_total/(R_g+R_s) ≥ 4·1.41² = 7.95; this device has
26.38 Ω / 15.83 Ω = 1.67. **f_max is primarily a measurement of output
resistance, and this model has none to measure.** The intrinsic voltage gain
g_m/g_ds runs 0.00514–0.09648 over the same bias range, one to two orders of
magnitude below a real GFET's. The g_ds that would put this device at
f_max/f_T = 1.41 at its own R_g, R_s and C_gd is 7.859726×10⁻³ S
(0.1965 mS/µm) against the model's 3.810651×10⁻² S (0.9527 mS/µm): **the
missing physics is worth a factor of 4.8483.**

**Result 3, a caution that applies to f_T and therefore to the ≈20 GHz
figure quoted throughout this chapter.** Because g_m = d(V_ds/R)/dV_g is
nearly proportional to V_ds in a resistor model, peak f_T here scales
nearly linearly with drain bias:

| V_ds (V) | 0.05 | 0.1 | 0.3 | 0.6 | 1.0 |
|---|---|---|---|---|---|
| peak f_T (GHz) | 10.140 | 20.279 | 60.779 | 121.168 | 200.440 |
| f_max/f_T | 0.6442 | 0.6430 | 0.6385 | 0.6312 | 0.6218 |

That is ×19.8 over a ×20 bias range. A saturating device's f_T does not do
this. **So the agreement between this model's ≈20 GHz and the literature's
tens of GHz is a fact about the bias point `plot_fT_fmax()` happens to pass,
not independent corroboration** — the same model "agrees" with 200 GHz
record devices if asked at V_ds = 1 V. None of the numbers in Section 4.6
are withdrawn; what is withdrawn is reading their proximity to the
literature as evidence. Note also that f_max/f_T *falls* slightly with
drain bias here, where a real device's rises as it saturates: the sign of
that trend is itself a saturation diagnostic, and this model has it
backwards.

**What this section does not claim.** It does not claim the model is wrong
for the use Chapter 4 puts it to. A resistor-plus-contact model is the right
instrument for the contact-resistance questions of Sections 4.3–4.5, 4.7 and
4.8, which are its substantive results, and those are untouched. What it
claims is narrower and, for an RF reader, more important: **f_max from this
model is not a prediction of f_max.** It is a restatement of the device's
access-resistance-to-channel-resistance ratio. Making it a prediction needs
a saturation term in `transfer_characteristic()` — a drain-field-dependent
velocity, or an explicit saturation velocity v_sat with the usual
I_d = W·q·n·v_sat ceiling — and that is now the top open item for this
chapter in `AUTOMATION_LOG.md`, with a measured target (4.85× in g_ds)
rather than a direction.

Validations: 30/30 and 6 of 6 mutants killed, with control A and control B.
Five checks are tolerance-free, each naming the operation it is exact
under: bitwise reconstruction of the model's own f_max (exact under
summation *order* — re-associating the identical algebra breaks it, which
mutant M1 confirms), multiplication by zero, the 1/N² finger law at powers
of two, and the overlapping-limit identity where this module's
default-geometry path must agree with `rf_small_signal_model` bitwise. One
further result is worth reading off the mutation report: the 2026-10-03
frozen-default bug itself, reinstated as mutant M4, is killed by **exactly
one** check — and that check is the only one in the battery that retains a
*magnitude* rather than a sign. See
`notes/2026-10-04-fmax-was-a-resistance-ratio.md` and
`fmax_shortfall_decomposition.png`.

<a id="s4-6-3"></a>

*Moved from Chapter 4, §4.6.3. Original section number kept.*

### 4.6.3 A saturation term, and the third exhausted lever: the contacts now cap g_ds

Section 4.6.2 ended with a target that had a tolerance rather than a
direction: the missing physics was worth **×4.8483** in g_ds, and velocity
saturation was the named route. That term has now been built
(`graphene_velocity_saturation_model.py`, 26/26; mutation report 7 of 7 with
control A and control B) and evaluated against all three of the acceptance
criteria §4.6.2 set. **Two are met. The third is missed by 4.3×, and the
reason it is missed is the result of this section.**

#### The model

Soft-saturation (Caughey–Thomas) local drift velocity at β = 1 — the γ = 1
form of Feijoo *et al.* (2019) — with current continuity fixing I_d along the
channel. The channel-length integral then separates, so the drain current
stays closed-form even though v_sat varies with position:

```
  I_d = μWe·Q / (L + μ·S),     Q = ∫ n dV,     S = ∫ dV/v_sat
```

closed on the contacts by V_ds,ch = V_ds − I_d·R_c. Setting S = 0 recovers
drift-diffusion with no saturation along the *same code path*, which is what
makes the overlapping-limit comparison below a comparison and not a second
implementation. The model is added **alongside** `transfer_characteristic()`,
which is left bitwise unchanged; rewiring the repository onto it is a
deliberate act under the 2026-10-02 rule and is not performed here.

#### v_sat is not fitted, and that is load-bearing

Criterion A names a number and v_sat is the only knob that moves it, so
choosing v_sat to satisfy the criterion would be fitting the model to its own
acceptance test. v_sat is therefore taken from optical-phonon emission,
v_sat(n) = (2/π)Ω/√(πn), with ħΩ = 0.10 eV as published by Feijoo *et al.*
(2019). At this device's density that gives 6.63 × 10⁷ cm/s — **above**
Dorgan *et al.*'s measured 1–3 × 10⁷ cm/s band on SiO₂, i.e. the generous
end of the physics.

Inverting the question makes the gap falsifiable rather than rhetorical:
**criterion A requires v_sat = 6.50 × 10⁶ cm/s, below the measured band.** A
sweep across the available phonon energies (0.059 eV and 0.149 eV, the SiO₂
surface polar modes; 0.100 eV, Feijoo's fit; 0.196 eV, graphene's intrinsic
optical phonon) reaches a g_ds factor of **1.2491 at best**. No phonon energy
in the physical range satisfies criterion A. This also discharges the standing
question of whether the remote-polar-phonon cap on SiO₂ matters here: it is
now a measured 1.07–1.25 in g_ds rather than an assertion.

Dorgan *et al.*'s own best fit is β = 2, which saturates *less* hard at these
fields (24.6–26.0 % higher velocity at the largest field reached), so the
β = 1 choice is generous in that direction too and the shortfall is bounded
from both sides.

#### Criterion B — met

| V_ds (V) | f_max/f_T, resistor | f_max/f_T, saturated |
|---|---|---|
| 0.05 | 0.644461 | 0.663187 |
| 0.10 | 0.642966 | 0.683186 |
| 0.20 | 0.640511 | 0.728061 |
| 0.50 | 0.633440 | 0.911727 |
| 1.00 | 0.621775 | 1.384374 |

The slope reverses from −0.0227 to **+0.7212**, 31.8× the magnitude of the
wrong-signed behaviour it replaces. A saturating device's f_max/f_T rises
with drain bias; this one now does, which was the part of §4.6.2's diagnosis
that made a claim about mechanism.

**What this table is not.** The V_ds = 1 V row lands inside Feijoo *et al.*'s
1.3–1.4 band, and that is **explicitly not offered as corroboration**. §4.6.2
recorded why: peak f_T here scales ×19.8 over a ×20 drain-bias range, and a
model with a knob that size will meet a comparably wide literature band
somewhere. The claim is the sign, which is bias-independent; where the ladder
crosses 1.3 is a fact about where the ladder was stopped.

#### Criterion A — missed, and the miss is the finding

| | resistor | saturated | factor |
|---|---|---|---|
| g_ds at the resistor's peak-f_T bias | 3.811e-2 S | 3.355e-2 S | **1.1361** |
| g_ds at the saturated model's peak | 3.767e-2 S | 3.354e-2 S | 1.1230 |
| **required** | | | **4.8483** |

23.4 % of the requirement. The saturation term is not weak: μS/L = 0.3036 at
this bias, which alone would divide the **channel** conductance by 1.699. It
is diluted because saturation acts only on the channel while the contacts sit
in series with it. At the 40 µm / 8-finger geometry R_c,total = 15.0 Ω of a
≈ 26 Ω device, so:

> **Roughly half of g_ds is a contact resistance that no saturation mechanism
> can reach. The contacts, not the channel, now cap g_ds.**

This is the third candidate lever on the f_max shortfall to come back
non-binding. §4.6.2 exhausted the gate resistance (R_g = 0 *exactly* still
falls 28 % short) and the feedback capacitance (C_gd deleted entirely moves
the answer 0.11 %). Velocity saturation now joins them. **All three point at
the access resistance** — which is §4.5, §4.7 and §4.8 of this chapter, the
part of this work with the most content and the one place it carries a
negative result worth having. The RF thread does not end in a tuning problem;
it ends in the contact.

#### Criterion C — met, and it corrects this chapter's own model

The overlapping limit cannot be checked as bitwise agreement, because the
saturated model's S → 0 limit takes the *arithmetic* average of n along the
channel while `transfer_characteristic()` takes a *harmonic* one. Making four
artefacts converge as V_ds → 0 separates three mechanisms, which are
separable precisely because their orders in V_ds differ:

| mechanism | size at V_ds = 0.1 V | order in V_ds |
|---|---|---|
| quadrature rule (50-sample mean vs. trapezoid) | −4.3 × 10⁻⁵ % | 2.04 |
| **profile domain** | **+0.253555 %** | **1.00** |
| Jensen gap (⟨1/n⟩ vs. 1/⟨n⟩) | +3.8 × 10⁻⁴ % | 2.00 |

Two things follow, and both were predicted wrongly before being measured.

First, the `n_segments = 50` discretisation that §4.6.1 flagged and that
§4.6.2 called load-bearing — because g_ds is a V_ds difference of exactly that
quadrature — is **measured at 4.3 × 10⁻⁵ % and is not a problem at the
committed biases.** A long-standing caution closes with a negative result.

Second, the dominant discrepancy is one nobody here was looking for.
`transfer_characteristic()` profiles the local channel resistance over the
**full V_ds**, but over half of V_ds is dropped across the contacts and never
appears across the channel at all, so the model evaluates n(V_ch) over a
potential range about twice too wide. The error is first order in V_ds and
674× the Jensen gap the check was written to find.

**No number in this chapter is withdrawn.** The total is one-signed,
0.2539 % at V_ds = 0.1 V, far below anything §4.6 concludes from, and in the
direction that makes I_d *larger* — so every shortfall against literature
reported in this chapter is, if anything, understated rather than overstated.
What has changed is that the bound is measured, its mechanisms are separated
and ordered, and the larger one was not the one that had been suspected.

#### One unexplained observation, recorded rather than resolved

The peak-f_T bias **moves branch** between the two models: V_g = −1.1128 V
(hole side) for the resistor model, +2.8571 V (electron side) for the
saturated one. Both readings are reported above rather than one being chosen.
A figure of merit whose optimum jumps 4 V between two models of the same
device is either a real electron–hole asymmetry — which this chapter has
independent reason to expect — or an artefact of where g_m peaks, and this
section does not distinguish them.

See `notes/2026-10-05-saturation-arrived-and-the-contacts-ate-it.md`,
`velocity_saturation_model.png` and `velocity_saturation_output.txt`.

<a id="s4-6-4"></a>

*Moved from Chapter 4, §4.6.4. Original section number kept.*

### 4.6.4 The perfect contact: the fourth exhausted lever, and a factor of four the ratio hid

§4.6.3 ended by naming the contacts as what caps `g_ds` and left one
counterfactual unrun: `f_max/f_T` at `R_c = 0` *exactly*, in the saturated
model, in the way §4.6.2 ran `R_g = 0` and `R_s = 0`. This section runs it.

`R_c` has to be removed from **both** places it enters, and that is not a
formality. It is a series resistance inside the self-consistent solve for the
intrinsic channel drop, *and* it is `R_s = R_c/2` in the `f_max` denominator,
reached through a different module that reads `gfet.Rc_total` as a global.
Zeroing only the first — which is what the natural keyword-argument route
does — gives `f_max/f_T` = 0.552349, a **19 % degradation**. Zeroing both
gives 0.759363, an **11 % improvement**. The two answers have opposite signs,
and the misleading one is the plausible one, since a 19 % degradation is what
"the contacts are binding" predicts if read without care. Check C3 of
`graphene_perfect_contact_counterfactual.py` asserts the sign inversion so
that it is a recorded measurement rather than a caveat.

**The result, at the 40 µm / 8-finger geometry and `V_ds` = 0.1 V, at a fixed
gate bias:**

| quantity | baseline | `R_c = 0` coherent | change |
|---|---|---|---|
| `g_ds` [S] | 3.374491e-02 | 5.145076e-02 | ×1.5248 |
| `f_T` [Hz] | 2.07398e+10 | 7.46504e+10 | ×3.5994 |
| `f_max` [Hz] | 1.41692e+10 | 5.66868e+10 | **×4.0007** |
| `f_max/f_T` | 0.683186 | 0.759363 | +11.15 % |

Note that the comparison has to be at a **fixed** bias. With `R_c = 0` the
peak of `f_T(V_g)` leaves the sweep entirely and sits at the grid edge in both
models, because the contact is what produced the turn-over: `I_d` is capped by
`R_c` at large `|V_g − V_dirac|` while `C_gs` keeps growing, so
`f_T = g_m/2πC_gs` turns over. Remove the cap and it does not. Widening the
sweep moves the edge, not the physics. "At the peak-`f_T` bias" is therefore
not an available convention for this counterfactual, and check P1 asserts the
interior/boundary distinction rather than leaving it as a remark.

**The fourth lever is exhausted on the ratio.** 0.759363 is 58.4 % of the
bottom of Feijoo *et al.*'s 1.3–1.4 band. A residual `g_ds` requirement of
**×2.9955** survives the perfect contact, and §4.6.3 bounded what `v_sat` can
supply at ×1.2491 across the whole physically available phonon range. So all
four quantities that have been proposed as the `f_max` shortfall's cause —
`R_g`, `C_gd`, `v_sat`, `R_c` — have now been set to their most favourable
value and none closes the gap.

**And the ratio hid a factor of four.** The same counterfactual is worth
×4.0007 in `f_max` itself. Section 6b of the module sweeps contact resistances
that have been measured rather than assumed:

| `R_c` per contact [Ω·µm] | `f_max` [GHz] | `f_max/f_T` | share of the `R_c = 0` gain |
|---|---|---|---|
| 0 | 56.687 | 0.759363 | 100.0 % |
| 65 (lowest reported) | 39.526 | 0.722202 | 59.6 % |
| 165 (Feijoo *et al.*) | 24.308 | 0.695514 | 23.8 % |
| **300 (this model)** | **14.169** | **0.683186** | 0.0 % |
| 470 (practical Ni flow) | 8.268 | 0.680047 | −13.9 % |
| 4000 (typical literature) | 0.222 | 0.698615 | −32.8 % |

> **ANNOTATION 2026-11-10 — every `f_max` in this table is an *isothermal*
> number, and that is now quantified rather than merely admitted. No value
> above is changed.**
>
> These are computed at `V_ds` = 50 mV. `graphene_self_heating_model.py`
> solves the coupled electro-thermal problem for this device and finds a
> self-consistent temperature rise of **1.0675 K** at that bias, costing
> **1.8 × 10⁻⁶** in drain current. The isothermal assumption behind this
> table is therefore sound *at this bias*, to far better than the precision
> printed, and the six-decimal `f_max/f_T` column is not being flattered by
> a neglected mechanism.
>
> This is the opposite of what the 2026-10-06 reading of Feijoo *et al.*
> implied, and the reason is a bias difference rather than a disagreement:
> their 65 → 40 GHz degradation is measured on a device driven at volts,
> where the same network here gives Δ*T* = 168.7 K and the device leaves
> the regime its transport parameters were fitted in (past the ~200 K
> fit-exclusion ceiling of Dorgan, Bae and Pop at `V_ds` = 1.25 V). **So
> the table is safe and the extrapolation from it is not**: a reader must
> not carry these numbers to a bias where the channel is 170 K hot. See
> §7.8.1c and `self_heating_output.txt`.
>
> One consequence for this table specifically. The `R_c` sweep above holds
> the *channel* fixed and varies the contacts, and self-heating couples to
> the total dissipation including the contacts — so at a bias where heating
> matters, the low-`R_c` rows would heat *more* (they pass more current),
> not less, partially eating the gain the table attributes to them. At
> 50 mV that coupling is 1.8e-06 and ignorable. It is recorded because the
> sign is counter-intuitive and the table invites exactly that
> extrapolation.

Two readings, and the second was not anticipated.

First, **this model's contact is 1.8× worse than that of the devices its
`f_max/f_T` is compared against**, which was not known before this section and
which closes an obvious objection: give this model Feijoo's own 165 Ω·µm and
it reaches 0.695514, not 1.3. The shortfall is not a contact-quality artefact.

Second, **`f_max/f_T` is non-monotonic in `R_c`.** It bottoms out at 0.680047
at 470 Ω·µm and rises again to 0.698615 at 4000 Ω·µm, while `f_max` falls
**37×** across those same two rows. This is §4.6.2's resistance-ratio identity
reasserting itself inside the *saturated* model: once `R_c` dominates
`R_total` faster than it dominates `R_g + R_s`, the ratio is rewarded for a
worse contact. There is a region of this design space where the ratio and the
figure of merit it is supposed to summarise **move in opposite directions**,
and no `f_max/f_T` criterion can be trusted inside it. §7.8.1c draws the
methodological consequence.

**One prediction of §4.6.3 is confirmed, from the other side.** That section
attributed the ×1.1361 `g_ds` factor to dilution by the contacts, computing
from a channel integral that the undiluted factor would be 1.699. That is a
testable prediction at `R_c = 0`, and the measured factor there is
**1.70661** — 100.45 % of it — obtained from a finite difference of a solved
terminal current rather than from the integral that produced the 1.699
(check D1).

**One claim is narrowed, by literature rather than by computation.** Feijoo,
Pasadas, Bonmann *et al.*, *Nanoscale Advances* **2** (2020), report that the
largest `f_max` in their measured GFETs occur "far from the saturated velocity
regime", at ~45 % of `v_sat`, with the **diffusion** contribution to the
current comparable to drift. Equation (4) of §4.6.3 is drift-only. §4.6.3's
criterion A therefore remains correct arithmetic about *this model* and is not
a statement about real devices, and the absence of a diffusion term — not the
absence of saturation — is now the RF thread's binding question. See §7.8.1c
and §7.9 item 7.

See `graphene_perfect_contact_counterfactual.py`,
`perfect_contact_counterfactual.png`,
`perfect_contact_counterfactual_output.txt` and
`notes/2026-10-06-the-perfect-contact-and-the-ratio-that-divided-out-the-prize.md`.

<a id="s4-6-5"></a>

*Moved from Chapter 4, §4.6.5. Original section number kept.*

### 4.6.5 The diffusion term: a boundary term worth ±0.5 %, and a potential this chapter never named

§4.6.4 closed by naming the absence of a diffusion term in Eq. (4) as the RF
thread's binding question, on the strength of Feijoo *et al.*'s measurement
that diffusion is comparable to drift at the peak-`f_max` bias. The term has
now been derived and built (`graphene_diffusion_current_model.py`,
**30 checks, 30 passed**), and the result divides into three parts that do not
agree with each other.

**The derivation, and the `beta = 1` separability it preserves.** Zebrev's
generalised Einstein relation ([arXiv:1102.2348](https://arxiv.org/pdf/1102.2348),
Eq. 20) is `mu = e D/eps_D` with `eps_D = n/(dn/dmu_c)`, and graphene's linear
dispersion gives `eps_D = E_F/2` exactly in the degenerate limit, so
`D = mu E_F/(2e)`. Writing `V_F(n) = E_F/e = A_F sqrt(n)` with
`A_F = hbar v_F sqrt(pi)/e` — this thesis's own Dirac dispersion, with no new
parameter — the diffusion term collapses exactly:

```
    D dn/dx = mu n dV_F/dx
    I_d = W e mu n d(V - lambda V_F)/dx                                  (10)
```

so the drift-diffusion current is transport down the gradient of the
quasi-Fermi potential `Phi = V - V_F`, and the separability that makes Eq. (4)
closed-form survives intact:

```
    I_d = mu W e Q_D / (L + mu S_D)
    Q_D = Q + lambda (A_F/3) (n_s^{3/2} - n_d^{3/2})                     (11)
    S_D = S + int kappa(V)/v_sat(V) dV,     kappa = -lambda dV_F/dV      (12)
```

`lambda = 0` recovers Eq. (4) **bitwise** (check X1, 0.0 ULP at twelve
(`V_g`, `V_ds`) points), which is the overlapping-limit requirement of
§4.6.3's criterion C in its strongest available form.

**Eq. (11) is a boundary term, and that matters more than its size.**
`int n dV_F = A_F int n d(sqrt n) = (A_F/3)[n^{3/2}]` is an exact
antiderivative, so the diffusion contribution to the channel integral depends
only on the **source and drain densities** and not at all on the density
profile between them. §4.6.3's criterion C found that this chapter's
profile-domain handling carried a one-signed error 674× larger than the Jensen
gap it was written to look for; the diffusion term **cannot inherit that
error**. Check X6a asserts this structurally rather than arguing it — a
deliberately lopsided interior quadrature grid with the same endpoints leaves
the term bitwise unchanged — and X6c is the contrast, the quadrature route to
the same number moving by 0.4 % on the same regridding.

**The size, at this chapter's own biases, `V_ds = 0.1 V`:**

| bias | `I_d` (λ = 0) [A] | `I_d` (λ = 1) [A] | share | `Q_d/Q` |
|---|---|---|---|---|
| saturated-peak, `V_g = +2.857143 V` | 9.007791e-05 | 9.052996e-05 | **+0.5018 %** | +1.1621 % |
| resistor-peak, `V_g = −1.112782 V` | 8.910103e-05 | 8.865209e-05 | **−0.5038 %** | −1.1546 % |

**The correction is signed and its sign flips between branches.** It follows
`d|n|/dV`, and `carrier_density()` returns a magnitude
(`sqrt(n_eff² + n_puddle²)`), so on the hole branch `|n|` rises toward the
drain and the term subtracts. A model that has discarded the sign of the
carrier cannot settle whether that is the correct physical sign for holes;
this is recorded as an open item (§7.9 item 9) and not as a result. The share
of `I_d` is about half of `Q_d/Q` because `S_d > 0`: Eq. (12) **opposes**
Eq. (11), with `S_d/S = +1.1621 %`.

**The RF consequence.** At the literature 40 µm / 8-finger geometry:

| quantity | λ = 0 | λ = 1 | change |
|---|---|---|---|
| `f_T` | 2.073963e+10 | 2.085305e+10 | +0.5469 % |
| `f_max` | 1.416902e+10 | 1.421070e+10 | +0.2941 % |
| `f_max/f_T` | 0.683186 | 0.681468 | **−0.2514 %** |
| `g_ds` | 3.374491e-02 | 3.391519e-02 | +0.5046 % |

`f_max` rises while the ratio **falls** — §4.6.4's ratio finding in a third
mechanism, and this time the two move in opposite directions rather than
merely by different amounts. §7.8.1a's verdict is unaffected: 0.681468 is
52.42 % of 1.3 either way.

#### A prediction written down first, and refuted — with its mechanism measured

Zebrev's closed form for the local ratio, `kappa = C_ox/(C_Q + C_it)`
(his Eq. 62), **diverges** as `C_Q → 0` at charge neutrality, and Feijoo
*et al.* put the peak-`f_max` bias near the onset of bipolar conduction. Both
predict the diffusion share **peaks at the Dirac point**, and that prediction
was written into the module before its sweep was run.

Measured: the share **collapses** there — `−0.006925 %` at the nearest grid
point to the Dirac point against `0.5053 %` at the maximum, a factor of 73,
with the maximum two volts away at `V_g = −1.272 V`. The mechanism is this
chapter's own regularisation, measured rather than argued:
`n = sqrt(n_eff² + n_puddle²)` makes `d|n|/dV = (n_eff/n) dn_eff/dV`
**exactly 0.0** at `n_eff = 0`, while `C_Q` stays finite at `n_puddle`
(G4b: `0.000000e+00` against `−2.290589e+15 m⁻²/V` unfloored). **The puddle
floor flattens precisely the region where the capacitor ratio diverges.** The
diffusion-dominated regime is not small in this model — it is absent, and that
is a statement about the regularisation rather than about graphene. §4.6.4 was
right that the peak-`f_max` bias cannot be asked about; the reason is the
puddle floor, not the missing term, and it does not go away by adding terms.

#### Annotation (2026-10-07): §4.3's `C_q` is evaluated at the gate overdrive, not at `E_F/e`

`quantum_capacitance()` takes `(V_g − V_ch − V_dirac)` in the slot the
dispersion wants `E_F/e = V_F`. At the RF bias those are **2.057 V and
0.0977 V**, a ratio of **21.05**, so this chapter's `C_q` overstates the
dispersion-consistent `C_Q = (2e³/πℏ²v_F²)V_F` by that factor and would
understate `kappa` by it. §4.6.5's comparison with Zebrev therefore uses the
dispersion-consistent form. **No number in this chapter is withdrawn**: `C_q`
is used downstream only through the series factor `C_q/(C_q + C_ox)`, where
`C_q ≫ C_ox` makes the factor close to unity under either evaluation. The
figures and tables of §4.3 are annotated here rather than recomputed, per the
2026-10-02 rule that a changed rule makes derived numbers stale until
deliberately re-derived.

#### What this section does not claim

It does not claim that Eq. (4) was missing this term. §7.8.1d sets out why
that question is not decided by the measurement above: the standard GFET
compact model this chapter's Eq. (4) descends from writes its drift expression
in the **quasi-Fermi** potential, under which Eq. (4) is already a complete
drift-diffusion current and `lambda = 1` double-counts. `lambda = 1` is
therefore **not** wired into any other module of this thesis, and the
production path remains `lambda = 0` with ±0.50 % as a stated bound.

See `graphene_diffusion_current_model.py`, `graphene_diffusion_mutation.py`,
`diffusion_current_model.png`, `diffusion_current_output.txt`,
`diffusion_mutation_output.txt` and
`notes/2026-10-07-the-diffusion-term-was-a-boundary-term-and-the-potential-was-never-named.md`.

<a id="s4-6-6"></a>

*Moved from Chapter 4, §4.6.6. Original section number kept.*

### 4.6.6 The potential is named — and the naming question was standing in front of a larger one

§4.6.5 ended with a variable that had two readings and a term that was exactly
the difference between them. This section names the variable, reports what that
does to §4.6.5's result, and reports the thing the naming question had been
hiding.

#### The declaration

> **`V_ch` is the quasi-Fermi (electrochemical) potential of the channel
> carriers, in volts, measured from the source.**

It is declared in `graphene_fet_model.py` itself, as a machine-readable
in-source marker above `carrier_density()`, so that the repository's statement
of what its own variable means travels with the variable rather than living in
a note. Three grounds, none of them new work:

1. The compact model Eq. (4) descends from defines exactly that variable.
   Pasadas and Jiménez (*IEEE TED* 63(7) 2016,
   [arXiv:1605.08235](https://arxiv.org/pdf/1605.08235)) write `v = μF` with
   `F = −dV/dx` and state that *"V(x) is the quasi-Fermi level along the
   graphene channel"*.
2. Classical long-channel theory uses the same variable and states what it
   buys: the channel variable `V(y)` is the quasi-Fermi potential, and the
   current-density equation it feeds is *"both drift and diffusion"*
   ([long-channel MOSFET notes](https://dunham.ece.uw.edu/ee531/Long_Channel_MOSFET.pdf)).
   In the quasi-Fermi variable, **one term is the complete current**.
3. `carrier_density()` already assumes it. The quantum-capacitance series
   factor `C_q/(C_q + C_ox)` is the correction one applies when the channel
   variable is the quasi-Fermi level and the graphene drop `E_F/e` is carried
   separately. Under a strictly electrostatic reading, that factor is itself a
   double count.

**What the literature does not say, recorded because it bears on the thesis.**
A review of GFET compact models (Lu *et al.*,
[arXiv:1703.09759](https://arxiv.org/pdf/1703.09759)) does *not* make this
identification: it calls the imref splitting "a.k.a. channel voltage" in one
section and describes `V(x)` only as "the voltage along the channel" in the
compact-model sections, never connecting them; its Eq. (12) carries a diffusion
term that its Eq. (15) writes with the same expression as drift; and its
compact model uses a drift form while stating only that drift-diffusion "is
assumed". **The ambiguity §4.6.5 found is inherited from the compact-model
literature rather than introduced here** — which does not excuse it in a
thesis, and does raise the value of declaring it explicitly.

#### §4.6.5's term is a double count, and its numbers stand

`graphene_diffusion_current_model.py` writes the current as `d(V − λV_F)/dx`
and names `V − V_F` the quasi-Fermi potential. Under the declaration `V` *is*
that potential, so `λ = 0` is already the complete drift-diffusion current.

**The term built in §4.6.5 is therefore a double count, not missing physics.**
Measured three ways at `V_g` = 2.0 V, `V_ds` = 0.1 V:

| Level | `I_d` (λ = 0) | `I_d` (λ = 1) | difference |
|---|---|---|---|
| intrinsic (fixed channel drop) | 1.432116×10⁻⁴ A | 1.442091×10⁻⁴ A | **+0.6965 %** |
| terminal (contact feedback on) | 8.190217×10⁻⁵ A | 8.225263×10⁻⁵ A | **+0.4279 %** |
| §4.6.5, at its saturated-peak bias | — | — | +0.5018 % |

The terminal figure must be the smaller of the first two, because the series
contact resistance is a negative feedback on `I_d`; that inequality is checked
rather than assumed, and its sign is fixed by the circuit, not by a fit.

**No number in §4.6.5 is deleted or changed.** The derivation is correct, the
boundary-term structure is correct, and the measured size is correct. What
changes is the direction: the production path is `λ = 0`, and `λ` is retained
as the instrument that *measures* the difference between the two readings. The
±0.50 % ambiguity §4.6.5 attached to every `I_d` in this chapter is therefore
**resolved rather than bounded**.

#### And the larger question the naming had been standing in front of

Declaring `V_ch` quasi-Fermi is not free: it promises that the charge relation
is the one that variable belongs to. Under the declaration the gate drive
divides exactly,

```
V_g − V_dirac − V_ch  =  e·n/C_ox  +  E_F(n)/e ,    E_F(n)/e = A_F·√n    (4.28)
```

which is a quadratic in `√n` with a closed-form root and no new parameter —
`A_F` is this thesis's own Dirac dispersion. §4.3's `carrier_density()` instead
uses the linearised series factor

```
n = (C_ox·dV/e) · C_q/(C_q + C_ox)                                      (4.29)
```

(4.29) is the linearisation of (4.28). [**SUPERSEDED 2026-10-10, §4.6.8:**
this sentence is wrong and is kept as written. (4.29) is (4.28)'s *exact
differential relation* used algebraically — a one-point rectangle rule for
∫(4.31) — compounded with `C_q` evaluated at the overdrive rather than at
`E_F/e`. The two forms do not share a leading order as `dV → 0`: (4.28) is
quadratic, (4.29) linear. **The −0.6023 % below is unaffected and was
reproduced independently.**] Running Eq. (4) on (4.28) with nothing
else changed:

| | `I_d` at the RF bias |
|---|---|
| on (4.29), λ = 0 — **the committed model** | 8.190217×10⁻⁵ A |
| on (4.28), λ = 0 — the declaration taken exactly | 8.140885×10⁻⁵ A |
| **cost of the declaration** | **−0.6023 %** |
| the double count it removes | +0.4279 % |

**The charge-model question is 1.4 times larger than the transport term this
chapter spent 2026-10-06 and 2026-10-07 on.** That is the finding of this
section, and it is not a flattering one: two sessions went into a term worth
+0.43 % while a −0.60 % question about the chapter's own charge relation sat
unexamined directly underneath it, reachable in one line of code.

**A number this section declines to quote.** The raw (4.28)-vs-(4.29) charge
disagreement reaches **+68 %** at `dV` = 0.01 V. Quoting that would be wrong:
`n_puddle` = 5×10¹⁵ m⁻² swamps both forms there and no result in this chapter
lives in that regime. In what `carrier_density()` actually returns, the worst
disagreement across the drives this chapter sweeps is **+2.22 %**; in `I_d` it
is **−0.60 %**. The raw figure would have overstated the cost by a factor of
31. This is §4.10's lesson and 2026-10-06's rule applied to this section's own
result.

#### Annotation to §4.3 and §4.6.5, in place

§4.3's `carrier_density()` docstring is annotated in the source rather than
rewritten: the wording "electrostatic charge from `V_g − V_ch`" is retained and
marked loose, since under the declaration the exact relation is (4.28) and the
series factor is its linearisation rather than an additional physical effect.
[**SUPERSEDED 2026-10-10, §4.6.8:** read "its exact differential relation,
used algebraically" for "its linearisation". The statement that it is not an
additional physical effect stands.]
**The numbers that function returns are unchanged — verified bitwise over a
41-point sweep — and what changed is what they are a model of.**

#### A numerical defect found on the way, measured and left open

`quantum_capacitance()` returns a finite value up to `dV` = **18.3493 V** and
`inf`/`NaN` above it, because `log(2(1 + cosh η))` overflows once
`η = dV/(kT/e)` exceeds about 710 (`kT/e` = 25.85 mV at 300 K). The exact
large-`η` limit is `|η| + log 2`. The onset is far outside the ~2.7 V maximum
drive this chapter sweeps, so **no committed number is affected**; it is
recorded as §7.9 item 13 rather than patched here, because that module owns
committed transcripts.

#### What this section does not claim

1. It does **not** claim (4.29) is wrong. (4.28) and (4.29) have been compared
   only to each other, never to measurement. Under the declaration (4.29) is
   non-exact; whether it is less accurate is not established here.
2. It does **not** re-derive this chapter. Every committed number remains a
   (4.29) number. The −0.60 % is a measurement of the gap, obtained by
   replacing one function and re-running, not a new result set.
3. It does **not** touch the hole branch. Both (4.28) and (4.29) carry `n` as a
   magnitude, which §4.6.5 already identified as leaving every branch-asymmetric
   result in this chapter unverified. The declaration does not help, and the
   signed-carrier item (§7.9 item 9) now **blocks** the (4.28) re-derivation
   rather than sitting beside it: re-deriving on a magnitude would rebuild the
   same flaw in a new equation.

Machine output: `potential_declaration_output.txt` (10 checks, 10 passed) and
`potential_census_output.txt` (16 checks, 15 passed, 1 expected failure).
Modules: `graphene_potential_declaration.py`,
`graphene_potential_census_audit.py`. Note:
`notes/2026-10-08-the-declaration-and-the-question-underneath-it.md`.

<a id="s4-6-7"></a>

*Moved from Chapter 4, §4.6.7. Original section number kept.*

### 4.6.7 The sign was computed and discarded — and the disorder floor is 2.07× Chapter 3's

**2026-10-09.** §4.6.6 declared `V_ch` and left the exact charge relation (4.28)
as §7.9's top item, blocked by item 9. Item 9 is now closed, and closing it
answered two other open items and falsified a docstring. Module:
`graphene_fet_signed_carrier_model.py` (19/19 checks). Transcript:
`fet_signed_carrier_output.txt`. Figure: `fet_signed_carrier_model.png`.
Session notes:
`notes/2026-10-09-the-sign-was-computed-then-thrown-away-and-the-floor-was-never-a-floor.md`.

#### Why (4.28) could not be solved before today

Both terms on the right of (4.28) are **odd** in `n`: `e·n/C_ox` trivially, and
`E_F(n)/e = A_F·sign(n)·√|n|` because the Fermi level of a hole gas lies below
the Dirac point. A relation with two odd sides has an odd root, and an odd root
cannot come out of a solver whose first line is `dV = np.abs(dV)` — which is the
first line of the 10-08 solver. That solver returns the correct *magnitude*;
what no object in this repository could return was **which branch the channel is
on at a given point**.

`carrier_density()` is the origin. It *computes* the sign — `n_eff` is odd in
`dV` because `C_q` is even — and then discards it in the quadrature floor
`√(n_eff² + n_puddle²)`. Everything downstream takes `abs(n)`, so for the
conductivity the discard is harmless; for (4.28) it removes the branch.

`n_net_exact()` restores it. Validated against exactly known values, not ranges:
the root is odd **bitwise** over 801 points, and substituting it back into
(4.28) leaves a residual of **7.690e-16** relative over `dV ∈ [−2, 2] V`.
**§7.9 item 12 is unblocked.**

#### The decomposition, and what it reveals `n_puddle` to be

    n_net = n_e − n_h  (electrostatics, signed),  S = √(n_net² + n_puddle²)
    n_e = (S + n_net)/2,    n_h = (S − n_net)/2

- `n_e + n_h` **is** `S`, i.e. exactly what `carrier_density()` returns. So the
  quantity this chapter has been flooring since August is the **total conducting
  density** — not the net density, not a per-species density.
- `n_e·n_h = n_puddle²/4` to 3.603e-16 relative: the August quadrature floor,
  chosen for smoothness, is algebraically a **mass-action law**.
- `n_e(+dV) = n_h(−dV)` bitwise: graphene's ambipolar symmetry, machine-checked.

**Contract: no number in this chapter moves.** The two-carrier conductivity
equals the shipped conductivity to 3 ulp (3.609e-16 relative, against the four
significant figures this chapter quotes), and the shipped sweep is re-run after
the module's monkey-patching with the constants verified restored.

[Wiedmann et al., Phys. Rev. B **84**, 115314 (2011), arXiv:1107.3929] splits
net-from-gate and total-from-disorder in the same way and argues from the
divergence of `1/R_H` at neutrality that `n + p` must remain finite —
independent, measurement-side support for the identification above. It imposes
**no** electron–hole product, so the mass-action law is a *stronger* assumption
than the literature makes and is recorded as a limitation, not a result.

#### §4.9's floor item, open since 2026-09-29 and untouched for eleven sessions

With `n_puddle` named as the total, `dV = 0` gives `n_e = n_h = n_puddle/2`
exactly and the floor is a one-line calculation:

| | value | sheet resistance |
|---|---|---|
| this chapter: `σ_min = e·μ·n_puddle` | 3.204353 × 10⁻⁴ S/sq | **3.1208 kΩ/sq** |
| Chapter 3 §3.3, measured `4q_e²/h` | 1.549618 × 10⁻⁴ S/sq | **6.4532 kΩ/sq** |

**Chapter 4 is 2.0678× more conductive at neutrality than Chapter 3.** The
`n_puddle` that reconciles them is 2.417989 × 10¹⁵ m⁻² = 2.418 × 10¹¹ cm⁻²,
against the shipped 5 × 10¹¹ cm⁻². (§3.3's 6.45 kΩ/sq is reproduced here from
CODATA to 0.0032 kΩ/sq rather than re-typed, so the comparison cannot drift.)

Quoted as a factor first, per 10-03: a PASS/FAIL on "do the two floors agree"
returns FAIL and says nothing about whether the gap is 2× or 200×.

**This is decidable only because item 9 is closed.** A quadrature floor on a
*magnitude* cannot say whether it floors the net density, the total, or each
species — three readings differing by factors of 1, 1 and 2 in `σ_min`, the same
order as the answer. And the convention is not the explanation: under the
per-species reading the total would be `2·n_puddle` and this chapter would be
**4.1357×** too conductive, so that reading moves the discrepancy the **wrong
way**.

**What is not claimed.** `4q_e²/h` is a quantum/ballistic minimum conductivity
observed experimentally; `e·μ·n_puddle` is a diffusive disorder floor. They are
not required to be equal. The defect is that this thesis quotes both, in
chapters that feed each other, and had never compared them — and that the
diffusive one lands **below** the measured one, i.e. this chapter's channel at
neutrality is less resistive than any measured graphene sheet.

#### What the factor costs, and the contacts attenuate it again

At `V_ds` = 0.1 V, on-state `V_g` = 3.5 V, off-state at `V_g = V_dirac`:

| | shipped | Ch.3-matched | on/off multiplier |
|---|---|---|---|
| terminal (`R_c` = 600 Ω) | on/off **1.2408** | 1.7820 | **1.4362×** |
| intrinsic (`R_c` = 0) | on/off **1.6146** | 2.8024 | **1.7357×** |

`I_off` falls 35.19 % at the terminals, 51.58 % intrinsically. The 2.0678× floor
error reaches the on/off ratio as 1.7357× intrinsically and only 1.4362× at the
terminals: **the contacts attenuate it by a further 1.2085×** — §4.6.6's
negative-feedback structure, on a completely different parameter, which makes it
a property of the device topology rather than of the charge relation.

**§4.4's on/off ratio is therefore understated**, by up to 1.44× as measured at
the terminals, *if* §3.3's measured floor is the right one for this device.
§4.4's numbers stand and are not deleted; this is the limitation attached to
them. It does not disturb §7.2's conclusion — an on/off ratio of 1.78 is no more
a logic transistor than one of 1.24.

#### Annotation to §4.3 and §4.4, in place: `n_puddle` is not a Dirac-point term

§7.9 item 11, created 10-07, asked whether `n_puddle` removes other regimes.
It does not remove regimes; **it adds a floor to all of them.** Measured:

| `V_g` (V) | \|n_net\| (m⁻²) | `n_puddle`/\|n_net\| | floor adds to total |
|---|---|---|---|
| 0.90 | 2.357 × 10¹⁴ | 21.22 | +2023.9 % |
| 2.00 | 2.870 × 10¹⁵ | 1.742 | +100.9 % |
| 3.00 | 5.265 × 10¹⁵ | 0.950 | +37.9 % |
| **3.50** | **6.462 × 10¹⁵** | **0.774** | **+26.4 %** |

At the largest overdrive this thesis sweeps, the *residual* puddle density is
**77 % of the gate-induced density** and still adds **+26.4 %** to the total.
Its contribution falls below 1 % only at `dV` = **14.73 V**, 5.46× anything swept
here. `carrier_density()`'s docstring — the floor "regularizes `n` at the Dirac
point" — describes a regime this thesis never enters. The docstring's wording is
annotated in place rather than deleted.

This was found by a **check written to pass**: §5 of the module asserted that
`I_on` is untouched by an off-state parameter, on the docstring's authority, and
it failed at −6.92 %. 09-28 established that prose is a detector; today adds
that so is a sanity check one expects to pass.

#### Three exactness checks failed as first written, and none was a code defect

Recorded because the *repair* is the transferable part.

1. The ambipolar symmetry failed "bitwise" by 2.0 /m². The symmetry is not
   broken — it is bitwise in `dV`. `(V_dirac − x) − V_dirac` is not
   `−((V_dirac + x) − V_dirac)` in binary floating point (worst 2.220e-16 V at
   `x` = −2 V), and the `C_ox/e` lever of 2.395 × 10¹⁵ /(V·m²) turns that into
   2.0 /m². **Every sweep in this thesis is a `linspace` over `V_g` followed by
   a subtraction of `V_dirac`**, so any bitwise symmetry claim must name its
   parameterisation. The residual is accounted for, not thresholded: feeding the
   reconstructed `dV` back into the `dV`-level model reproduces the `V_g`-route
   species bitwise.
2. "`n_e − n_h` recovers `n_net`" failed at 1.772e-14 and passes at 1.5e-16 —
   the same 1-ulp residual, **119× apart**, because `n_net → 0` at neutrality
   while `S` stays at the floor. The denominator decided PASS from FAIL: §7.9's
   10-06 ratio-criterion item, in a check written by the audit itself.
3. The `σ` identity failed "bitwise" at 3 ulp and was **not** repaired by
   widening a ulp budget — 09-29's α-anchor policy — but by stating the identity
   and naming the slack against the precision this chapter commits.

#### A number item 13 produced on its first use outside its own module

§4.6.6 recorded the `cosh` overflow above `dV` = 18.3493 V as item 13, latent,
"no committed number affected". The first draft of the 1 %-crossover bisection
above scored the resulting `NaN` as "the floor dominates" (`nan > 0` is `False`)
and returned its own 200 V bracket ceiling as the answer; the correct value is
14.73 V. The bracket is now capped below the onset and `NaN` is a hard stop.
**Latent and harmless are not the same property.** Item 13's ranking is
unchanged — it is a numerics change inside a module with committed transcripts
across three chapters — but "no committed number affected" is no longer the
whole of its description.

#### What this section does not claim

- It does not claim `n_puddle` should be changed. The 2.418 × 10¹⁵ m⁻² figure is
  what would reconcile two floors that are not required to agree; adopting it
  would move §4.4, §4.7 and §4.10 and is not done here.
- It does not claim the mass-action closure is physical. It is algebraically
  equivalent to a smoothness choice, and the literature it most resembles
  declines to impose it.
- It does not claim the signed model is a new physical model. Its contract is
  the opposite: every quantity it introduces is a decomposition of one already
  being computed, and the conductivity is unmoved to 3 ulp.
- It does not re-derive anything on (4.28). That is item 12, now unblocked and
  still open. [**CLOSED 2026-10-10, §4.6.8.**]

<a id="s4-6-8"></a>

*Moved from Chapter 4, §4.6.8. Original section number kept.*

### 4.6.8 Chapter 4 re-derived on (4.28) — and (4.29) was never a linearisation

§4.6.6 declared `V_ch` the quasi-Fermi potential, made (4.28) the exact charge
relation, measured the gap to (4.29) at −0.6023 % of `I_d`, and left the
re-derivation as §7.9 item 12. §4.6.7 supplied the branched root the item was
blocked on. This section does the re-derivation. **The −0.60 % reproduces
exactly; the finding is that the chapter has been describing the pair of
equations wrongly, in three places, since 2026-10-08.**

Machine output: `exact_charge_rederivation_output.txt` (21 checks, 21 passed).
Module: `graphene_exact_charge_rederivation.py`. Figure:
`exact_charge_rederivation.png`.

#### (4.29) is not the linearisation of (4.28), and the test is one line

A linearisation of a function agrees with it to first order by construction.
These two do not share a *leading* order. In (4.28), `a·n + b·√n = dV` with
`a = e/C_ox` and `b = A_F`; as `dV → 0` the `√n` term dominates and
`n → (dV/b)²` — **quadratic**. In (4.29), `C_q(0)` is finite, so
`n → (C_ox/e)·f(C_q(0))·dV` — **linear**. Measured log-log exponents:

| `dV` (V) | d ln n/d ln dV, (4.28) | d ln n/d ln dV, (4.29) |
|---|---|---|
| 10⁻¹⁰ | **1.99999854** | **1.00000000** |
| 10⁻⁶ | 1.999387 | 1.000000 |
| 10⁻⁴ | 1.943766 | 1.000000 |
| 10⁻² | 1.274493 | 1.002201 |

#### What (4.29) actually is: the exact differential relation, used algebraically

Differentiating (4.28),

```
d(dV)/dn = e/C_ox + (1/e)·dE_F/dn = e/C_ox + e/C_Q(n),
C_Q(n) := e²·dn/dE_F = 2·e·√n / A_F                                  (4.30)

  =>  dn/d(dV) = (C_ox/e) · C_Q(n)/(C_Q(n) + C_ox)                   (4.31)
```

(4.31) **is** (4.29) with `dn/d(dV)` in place of `n/dV`. So the series factor
is not an approximation *of* the charge relation; it is the charge relation's
exact derivative, and (4.29) uses it as a **one-point rectangle rule** for
`∫₀^dV (4.31) dV'`. That is a quadrature error, not a truncated expansion, and
a rectangle rule is worst where the slope varies fastest — near the Dirac
point, which is exactly where §4.6.6's raw disagreement reached +68 %.

(4.30) and (4.31) are verified as identities, not asserted: `e·dn/d(dV)` from
(4.28)'s closed-form root equals the series combination of `C_ox` with `C_Q(n)`
to **1 ulp**, and numerically integrating (4.31) back to the closed-form root
converges at **observed order 3.994** against Simpson's 4, with the `h⁴`
extrapolation landing **1.5×10⁻¹³** relative from it. An error in the prefactor
of `C_Q` — a dropped factor of 2, say — would leave an O(1) plateau and drive
the observed ratio to 1 instead of 16; that is what the order, rather than a
threshold on the error, is testing.

#### Two independent errors, separated for the first time

(4.29) differs from (4.28) in two unrelated ways: the quadrature above, and
`quantum_capacitance()` evaluating `C_q` at the gate overdrive `dV` instead of
at `E_F/e = A_F·√n` (noticed in §4.6.5's 2026-10-07 annotation and never
costed). Running the 2×2 — {rectangle, integral} × {argument `dV`, argument
`E_F/e`} — on the raw charge, with the puddle floor off so that the floor does
none of the work:

| `dV` (V) | total error of (4.29) | quadrature alone | **argument alone** | sum |
|---|---|---|---|---|
| 0.01 | +68.2193 % | +32.1722 % | **+68.0927 %** | +100.2649 % |
| 0.50 | +8.0543 % | +4.1176 % | **+7.2144 %** | +11.3320 % |
| 1.20 (RF bias) | **+5.2064 %** | **+2.6396 %** | **+4.7416 %** | +7.3812 % |
| 2.70 | +3.4729 % | +1.7523 % | **+3.2193 %** | +4.9716 % |

> **SUPERSEDED 2026-10-11, §4.6.9(c). Numbers retained, not rewritten.**
> Every figure in this table is differenced against (4.28), which is a T = 0
> relation, while the model is a 300 K one — so each of them contains a
> temperature difference as well as the knob it names. Re-measured in a 2×2
> at one temperature, the RF-bias row reads **argument +2.6408 %, quadrature
> +0.7603 %, total +3.0963 %**: this table overstates them by 1.80×, 3.47×
> and 1.68× respectively. The *ordering* claimed below survives — the
> argument error is still the larger contributor — and the conclusion that
> they are not additive survives, though the departure from additivity is
> 1.109 rather than 1.42 once the reference is at the model's temperature.
> The mistake is kept visible because its shape (a decomposition measured
> against a reference at a different operating point) recurs, and because
> §4.6.9's 300 K argument error of +2.6408 % is numerically within 0.0012 %
> of this table's quadrature figure, which reads as a reproduction and is not
> one.

Three things follow. **The argument error is the larger single contributor**
(+4.74 % against +2.64 % at the RF bias) — which no earlier session could have
known, because none separated them. **Neither is the explanation on its own.**
And **they are not additive**: the sum overstates the total by 1.42× at the RF
bias and by 1.47× at `dV` = 0.01 V, because they compose multiplicatively in
the series factor. A decomposition reported as a sum would have been wrong in
the direction that flatters it.

The puddle floor then attenuates all of this by **4.3×** at the RF bias
(+5.21 % raw → +1.22 % in what `carrier_density()` returns). That attenuation
is a property of `n_puddle`, which §4.6.7 measured to be 2.0678× off Chapter
3's floor — so it is a reason the re-derivation moves little, not a reason for
confidence in either form.

#### §§4.3–4.4 re-derived

| quantity | on (4.29) | on (4.28) | change |
|---|---|---|---|
| `n` at `V_g` = 3.5 V (m⁻²) | 8.170466×10¹⁵ | 8.000036×10¹⁵ | **−2.0859 %** |
| `R_ch` at `V_g` = 3.5 V (Ω) | 381.9555 | 390.0926 | +2.1304 % |
| `R_ch` at the Dirac point (Ω) | 624.1509 | 624.1509 | **0.0000 %** |
| `I_d`, drift, RF bias (A) | 8.719345×10⁻⁵ | 8.671268×10⁻⁵ | −0.5514 % |
| `I_d`, saturated, λ = 0, RF bias (A) | 8.190217×10⁻⁵ | 8.140885×10⁻⁵ | **−0.6023 %** |
| `L_x`, Pd 110 Ω·µm (nm) | 115.20 | 112.79 | −2.0859 % |
| `L_x`, 300 Ω·µm (nm) | 314.17 | 307.62 | −2.0859 % |
| `L_x`, 500 Ω·µm (nm) | 523.62 | 512.70 | −2.0859 % |
| terminal on/off | 1.240756 | 1.230617 | −0.8171 % |
| intrinsic on/off (`R_c` = 0) | 1.614570 | 1.581261 | −2.0630 % |

The −0.6023 % line **reproduces §4.6.6's published figure to four decimal
places from an independently written patch path**. Per 2026-10-02 that is a
fact about two artefacts rather than a validation; what it establishes is that
the (4.29) column above is the committed model and not a re-run with drift in
it.

At the Dirac point the two forms agree to 0.0000 % because the floor supplies
the entire density there — the one place in this chapter where the charge
relation does not matter at all, and §4.6.7 already showed it is also the
place the chapter's off-state numbers come from.

**Every qualitative claim of §§4.3–4.4 survives.** The best literature contact
still crosses over below the 200 nm channel and the mid- and worst-case
contacts still do not; the on/off ratio is still between 1 and 2, i.e. still
not a switch; the channel is still more resistive at neutrality than in the on
state. §4.4's contact-engineering conclusion and §7.2's "no switch" are
untouched.

**On/off moves less at the terminals (−0.8171 %) than intrinsically
(−2.0630 %).** This is the third parameter in three sessions to show that
structure — the diffusion term (2026-10-08), `n_puddle` (2026-10-09), and now
the charge relation itself. It is therefore a property of the **device
topology**, not of any one parameter: the series contact resistance is a
negative feedback on `I_d` and attenuates whatever moves upstream of it.

#### §4.6 re-derived — and the item was mis-scoped

Under (4.31) the small-signal gate capacitance is **not an independent
ingredient**: `e·dn/d(dV)` *is* `(1/C_ox + 1/C_Q(n))⁻¹`. So a consistent
re-derivation moves `C_gs` as well as `I_d`, and `f_T = g_m/(2πC_gs)` sees
both. Applied separately, at the peak-`f_T` bias, `V_ds` = 0.1 V:

| | on (4.29) | (4.28) charge only | (4.28) charge + `C_Q` |
|---|---|---|---|
| peak-`f_T` bias `V_g` (V) | −1.1115 | −1.1992 | **−1.1867** |
| `f_T` (GHz) | 20.2787 | 19.8734 | 20.2666 |
| `f_max` (GHz) | 18.5566 | 18.1781 | 18.5495 |
| `C_gs` (fF) | 0.0767 | 0.0767 | 0.0752 |
| change in `f_T` | — | **−1.9986 %** | **−0.0594 %** |
| change in `f_max` | — | −2.0399 % | −0.0386 % |

**The combined −0.06 % is a cancellation, not a smallness.** The charge half
is −1.9986 % and the capacitance half +1.9392 % — opposite in sign and within
3 % of each other in magnitude. A re-derivation that changed only the charge
would have reported −2.0 % and been wrong by a factor of 34; one that reported
only the combined number would have called item 12 negligible in `f_T` when it
is in fact two coupled 2 % corrections that happen to oppose. **An item that
names one of two coupled quantities gets the other one whether or not it asks,
and the sum is not the way to find out which mattered.**

The capacitance half is small for a reason worth stating: `C_q` at the
overdrive is 4.500×10⁻¹ F/m² and `C_Q` at `E_F/e` is 1.820×10⁻² F/m², a factor
of **24.7** apart, but both are ≫ `C_ox` = 3.837×10⁻⁴ F/m², so a factor of 25
in `C_q` becomes −1.94 % in `C_gs`. `E_F/e` = 0.0773 V is 4.04 % of the 1.91 V
overdrive, which is the argument error of the previous subsection seen
directly.

**The peak bias is the sensitive quantity, not the peak value.** `f_T` moves
0.0594 % in height and 6.76 % in *location* (−1.1115 → −1.1867 V), a ratio of
**114**. Both peaks sit on the **hole branch** and neither crosses to the
electron branch, so the open "does the peak-`f_T` bias change branch between
the two models" question of §4.6.7 is **answered NO for the (4.28)-vs-(4.29)
pair**; it was asked about the branched-vs-magnitude pair and stays open there.

At the literature-scale geometry (40 µm, 8 fingers) that §4.6.2 and §7.8.1a
quote: peak `f_T` 20.2787 → 20.2666 GHz, `f_max` at that bias 13.0393 →
13.0346 GHz, extrinsic `f_T` 3.4423 → 3.3847 GHz, extrinsic `f_max` 2.2105 →
2.1727 GHz. **§7.8.1a is untouched, and for the right reason**: its verdict is
structural — that `f_max` is a restatement of `R_total/(R_g+R_s)` — so the
quantity to check is the ratio, and `f_max/f_T` moves **+0.0228 %**, from
0.643007 to 0.643154, still below the bottom of Feijoo *et al.*'s 1.3–1.4 band
that 2026-10-04 showed a *perfect* gate cannot reach. Five levers have now
been tried — velocity saturation, the perfect contact, the diffusion term, the
exact charge relation, and its capacitance — and the binding constraint is
still that `transfer_characteristic()` is a resistor with no output resistance.

#### Three exactness checks failed as first written, and two were wrong in their prose

Continuing 2026-10-09's series, with a new instance of it.

**(a) The quadratic-exponent check failed at `dV` = 10⁻⁴ V**, returning 1.9438
against a ±0.01 criterion. It did not measure the exponent. The two terms of
(4.28) are equal at `dV_x = 2·A_F²·C_ox/e` = **6.519 mV**, so 10⁻⁴ V is only
65× below the crossover, where the correction is O(√(dV/dV_x)) ≈ 12 %. **It
measured the width of the crossover**, and the repair is to evaluate the limit
where the limit is and to print `dV_x` so the criterion comes from the model
instead of from a round number.

**(b) The integral identity failed at 1.1×10⁻⁷ with order 3.49 on a uniform
grid — and it measured the same 6.519 mV.** The integrand of (4.31) turns over
inside a layer of width `dV_x` at the lower limit, and a uniform grid over
[0, 1.2 V] puts only 22 nodes in it even at 4001 points. A grid graded about
`dV_x` reaches order 3.994 with **12× fewer** points. **Two exactness checks
written for different purposes, in different sections, both failed by measuring
the same number** — which is the strongest evidence this thesis has that
2026-10-09's rule is a rule and not an anecdote.

**(c) "Every number moves by under 2 %" failed at 2.0859 %, and the bound was
removed rather than widened.** The 2 % was a round number chosen before the
run and derived from nothing, so widening it to 2.5 % would have been exactly
the tuned tolerance §7.9 item 15 exists to find. The magnitudes are now
printed with no threshold on them and the check tests the qualitative claims,
which are what a re-derivation can actually overturn. **Of 2026-10-09's three
available repairs — loosen, fix the model, find out what it measured — this
session used the third twice and the fourth option, delete the criterion, once.**

Two further checks were wrong in their *prose* and failed on it: one asserted
the capacitance half was the *larger* half (it is 0.970× the charge half), and
one asserted `f_T` was "two orders below the 100–300 GHz RF-frontend
requirement" — wrong twice, since 20.3 GHz is about one order below 100 GHz
and §7.8.1a does not rest on a frequency threshold at all. 2026-09-28
established that prose is a detector; here it detected the author's own
checks.

#### A delivery fault found on the way

`contact_resistance_crossover.py` opens with `from graphene_fet_model import
carrier_density`, so it holds its **own** reference, bound at import time.
Rebinding `graphene_fet_model.carrier_density` does not reach it, and a
re-derivation that patched only the owning module would have left §4.4 on
(4.29) while reporting it as re-derived. This is the 2026-10-03
frozen-default fault in its monkey-patch form. An **arrival control**
(2026-09-30's rule) now sits in the module ahead of the §4.4 numbers and fails
if the patch did not reach the module that computes them.

#### What this section does not claim

1. It does **not** claim (4.29) is less accurate than (4.28). Neither has been
   compared to measurement. (4.29) is non-exact under the declaration and is
   now additionally *named* — a quadrature of the exact differential — which
   is a statement about its derivation, not about its accuracy.
2. It does **not** claim (4.28) is exact physics. `E_F = ħv_F√(πn)` is the
   **T = 0** dispersion and carries no thermal broadening; (4.28) is exact
   with respect to the declaration and the Dirac dispersion only. The
   compact-model review [arXiv:2209.00388](https://arxiv.org/pdf/2209.00388)
   uses `C_q ≈ k·c₁·√(1 + (V_c/c₁)²)`, whose `c₁ → 0` limit is exactly (4.30),
   and `c₁` is the thermal rounding (4.28) drops. That regime is far inside
   the puddle floor and is stated as a limitation rather than costed.
3. It does **not** replace Chapter 4's numbers. Every committed number remains
   a (4.29) number, annotated with its (4.28) value above. The chapter's
   production path is unchanged and is now carrying a measured, decomposed
   statement of what that choice costs.
4. It does **not** touch the hole-branch sign of this chapter's
   branch-asymmetric results. §4.6.7 made that expressible; it remains
   unverified.
5. It does **not** change the floor. `n_puddle` is carried into (4.28)
   unchanged so the comparison is apples to apples, which means (4.28) + floor
   is not itself an exact model either.

#### Annotation to §4.6.6 and §4.3, in place

§4.6.6's sentence "(4.29) is the linearisation of (4.28)" and the matching
wording in `carrier_density()`'s 2026-10-08 annotation are **marked
superseded where they stand, not rewritten**: (4.29) is the exact differential
relation of (4.28) used algebraically, compounded with an argument error. No
number in §4.6.6 changes — the −0.6023 % it published is reproduced here — and
the superseded wording is kept visible because the shape of the mistake
(describing a quadrature as a truncation) is the kind that recurs.

<a id="s4-6-9"></a>

*Moved from Chapter 4, §4.6.9. Original section number kept.*

### 4.6.9 `C_q`'s argument: the fix is a self-consistent solve, the obvious version of it deletes the temperature, and §4.6.8's decomposition was measured against the wrong reference

§4.6.8 named a second error alongside the quadrature error of (4.29) —
`quantum_capacitance()` is evaluated at the gate overdrive `dV` while the
dispersion it descends from requires `E_F/e = A_F·√n` — sized it at
**+4.7416 %** of the raw charge at the RF bias, called it the largest cost
attached to anything then open on this chapter, and made it §7.9 item 16.
This section acts on it. **The item's direction is right and all three of its
numbers are wrong, by factors of 1.68 to 3.47.**

Machine output: `quantum_capacitance_argument_output.txt` (18 checks, 18
passed, five against exactly known values). Module:
`graphene_quantum_capacitance_argument.py`. Figure:
`quantum_capacitance_argument.png`.

#### (a) The argument cannot be fixed by changing an argument

`quantum_capacitance(V)` takes a potential. The quantity it must be evaluated
at is `E_F/e = A_F·√n`, and `n` is what `C_q` is being used to compute. The
correct argument is a **state**, not a voltage, so the fix is a
self-consistent solve

    n = (C_ox·dV/e) · C_q(E_F(n)/e) / (C_q(E_F(n)/e) + C_ox)             (4.32)

and not an edit to a call site. That is the reason the largest open cost on
this chapter sat for a session as a one-line description: the item as written
understated the work, not the size.

At the RF bias, on the committed charge, `E_F/e = 0.062498 V` — **5.2082 %** of
the 1.2 V overdrive — and the two candidate arguments give

| argument | `C_q` (F m⁻²) | `C_q/C_ox` | series factor |
|---|---|---|---|
| `dV` = 1.2 V (committed) | 2.825142e-01 | 736.3 | 0.998644 |
| `E_F/e` = 0.0625 V | 1.575324e-02 | 41.06 | 0.976223 |

A **17.9×** difference in `C_q` is a **2.2966 %** difference in the series
factor. The reason the error is small is that both values are far above
`C_ox`, i.e. **the device is not quantum-capacitance-limited at this bias** —
not that the argument is nearly right. Where it would matter is exactly where
this chapter has least to say: near the Dirac point, where the puddle floor
already dominates.

#### (b) The obvious fix is a T = 0 expression, and it deletes the rounding the function exists for

The dispersion `C_Q(n) = 2e√n/A_F` is a zero-temperature result.
`quantum_capacitance()` is thermally rounded for a stated reason — its own
docstring says *"so C_q does not diverge to zero exactly at the neutrality
point"*. Substituting the dispersion removes that:

| | `C_q` at the Dirac point (F m⁻²) |
|---|---|
| committed thermal, at `dV` = 0 | 8.43739892e-03 |
| T = 0 dispersion, at `n` = 0 | **0** |
| thermal at `E_F/e`, at `n` = 0 | 8.43739892e-03 |

`C_q(0)` has a closed form, `PREFACTOR·(kT/e)·ln 4`, because `ln(2(1+cosh 0)) =
ln 4` exactly; the thermal-at-`E_F/e` form reproduces it to **0 ulp**. The
dispersion gives zero, hence a series factor of zero and an induced charge of
*exactly* zero at the Dirac point — the floor would supply 100 % of the
density there instead of its committed 99.9 %.

A sharper form of the same statement: with the dispersion, `n ≡ 0` solves
`dn/d(dV) = (C_ox/e)·f(C_Q(n))` through the origin, so the initial-value
problem is **degenerate** there. With the thermal form it is not. *Temperature
is what makes the charge relation well-posed at the Dirac point*, and it is
what the obvious fix throws away. §4.6.8 sized item 16 with that substitution,
so **its +4.7416 % is the cost of a fix that is wrong in a second way.**

The form used here, (4.32) with the thermal `C_q` at `E_F/e`, is called `B_T`.
It keeps the temperature, reduces to the dispersion as `T → 0` (checked, §(d)),
and keeps a finite `C_q` at neutrality.

#### (c) The decomposition, re-measured at one temperature

(4.28) is a T = 0 relation and the committed model is a 300 K one, so every
percentage in §4.6.8's 2×2 contains a temperature difference as well as the
knob it names. The repair is a 2×2 at **one** temperature — argument ∈ {`dV`,
`E_F/e`} × quadrature ∈ {rectangle, exact} — with the thermal `C_q`
throughout, and a reference `D_T` obtained by integrating the model's own
differential relation as an initial-value problem. That integrator is
validated against an exactly known value: with the dispersion in place of the
thermal form it must reproduce (4.28)'s closed-form root, and does, to
**−1.5e-14** relative.

Raw charge at `dV` = 1.2 V:

| | `n` (m⁻²) | vs `D_T` |
|---|---|---|
| A `C_q(dV)`, rectangle **(committed)** | 2.86980420e15 | +3.0963 % |
| `B_T` `C_q(E_F/e)`, rectangle | 2.80477739e15 | +0.7603 % |
| `C_T` `C_q(dV)`, exact quadrature | 2.85712585e15 | +2.6408 % |
| `D_T` `C_q(E_F/e)`, exact quadrature | 2.78361495e15 | 0 |

| | §4.6.8, vs T = 0 (4.28) | here, vs 300 K `D_T` | overstated by |
|---|---|---|---|
| argument error | +4.7416 % | **+2.6408 %** | 1.80× |
| quadrature error | +2.6396 % | **+0.7603 %** | 3.47× |
| both wrong | +5.2064 % | **+3.0963 %** | 1.68× |

§4.6.8's three figures are **marked superseded where they stand and not
deleted**, per this chapter's standing practice. What survives is the
*ordering*: the argument error is still the larger of the two, so item 16 is
still correctly ranked ahead of the quadrature error. What does not survive is
the headline number this chapter quotes in two places.

**A trap, recorded because it is the kind that is read as a confirmation.**
Today's argument error at 300 K (+2.6408 %) agrees with §4.6.8's *quadrature*
error measured at T = 0 (+2.6396 %) to 0.0012 %. They are different quantities
against different references. Anyone reading the two sections in sequence
would take the agreement for a reproduction. §4.6.6's rule — that a
coincidence between two artefacts is not a validation — is applied here for
the first time to a **numerical** coincidence rather than to a re-run.

#### (d) Item 13 was ranked last, and it is what blocked this section

The natural exact-limit check on `B_T` is `T → 0`: it must reproduce the
dispersion fixed point. It **could not be run on the committed code**.
`η = eV/kT` grows as `T` falls, `cosh(η)` overflows at **η = 710.4758**,
`ln(2(1+∞)) = ∞`, and the series factor returns `∞/∞ = nan`. That is §7.9 item
13, added 2026-10-08, judged *latent* and ranked **last** because η = 710 needs
**18.35 V** of overdrive at 300 K, which the device never sees. At T = 1 K and
the ordinary RF bias, η = 719.

> **Latency was judged against the inputs the model takes, not against the
> inputs its validations take.** A defect unreachable from a device's
> operating range can still be reached by a check on that device's model, and
> ranking by reachability missed it. The item ranked last gated the item
> ranked as the largest cost.

Item 13 is closed in the same session (§4.6.10). With the overflow-free
spelling in place the limit runs: at T = 10⁻³ K, `B_T` reproduces the
dispersion fixed point to **7e-16** relative, against `nan` before.

#### (e) What actually moves

`B_T` installed, with an arrival control on each patched name:

| | (4.29), `C_q(dV)` | `C_q(E_F/e)` | change |
|---|---|---|---|
| `n` at 3.5 V (m⁻²) | 8.170466e15 | 8.087032e15 | **−1.0212 %** |
| `R_ch` at the Dirac point (Ω) | 624.150907 | 624.150907 | **0.0000 %** |
| `L_x` Pd / 300 / 500 (nm) | 314.17 / 314.17 / 523.62 | 310.96 / 310.96 / 518.27 | −1.0212 % |
| terminal on/off | 1.240756 | 1.235801 | −0.3993 % |
| intrinsic on/off | 1.614570 | 1.598207 | −1.0135 % |
| `I_d` at the RF bias (A) | 8.190217e-05 | 8.169499e-05 | **−0.2530 %** |

The raw charge moves **−2.2659 %** and `n` at 3.5 V moves **−1.0212 %**: the
puddle floor attenuates by a factor **0.451**, which is why item 16 was
survivable as a change that moves committed numbers. `R_ch` at the Dirac point
is unchanged to the printed digit because the floor supplies the whole density
there.

Every qualitative claim of §§4.3–4.4 survives: the on/off ratio stays between
1 and 2, `L_x` at 300 Ω·μm stays above the 200 nm channel (314.17 → 310.96 nm),
and the resistance ordering by metal is unchanged. **The correction moves
numbers and no verdicts.**

RF, per-μm convention at `V_ds` = 0.1 V:

| | `C_q(dV)` | charge only | charge + `C_gs` |
|---|---|---|---|
| peak-`f_T` bias (V) | −1.1115 | −1.1617 | −1.1491 |
| `f_T` (GHz) | 20.2787 | 20.0574 | 20.4422 |
| `C_gs` (fF) | 0.0767 | 0.0767 | 0.0752 |

`f_T`: charge alone **−1.0911 %**, the consistent `C_gs` **+1.8976 %**, total
**+0.8065 %**. Opposite signs again. §4.6.8 found the same structure for item
12 and recorded it as a property of the **device topology** rather than of any
one parameter, on a single instance; it now has a second, independent one,
which is what that claim needed. `f_max/f_T` moves **+0.0009 %** and stays at
0.9151, so §7.8.1a is untouched by a **sixth** lever.

#### (f) The two-knob check's first use, and it fails in the direction the check cannot see

§4.6.8 created a backlog item for a check that can see a cancellation: move
each input separately, and flag the case where the combined move is *smaller*
than either. Its first instance is here, and it has the **opposite** sign.

| reference | argument only | quadrature only | sum | both | both/sum |
|---|---|---|---|---|---|
| 300 K, consistent | −2.2659 % | −0.4418 % | −2.7077 % | −3.0033 % | **1.109** |
| T = 0 (§4.6.8) | −2.4398 % | −0.4418 % | −2.8816 % | −4.9488 % | **1.717** |

The two knobs **reinforce**: the combined move is larger than the sum of the
separate ones, so a check that tests for smallness would not have flagged it.
The item is **amended, not closed** — it must report the ratio and flag any
departure from 1 in either direction.

And a requirement on the item that could only be found by running it: **the
interaction measure inherits the reference error of whatever it is differenced
against.** The T = 0 frame puts the departure from additivity at 0.717 where
the consistent frame puts it at 0.109 — an inflation of **6.6×**. A two-knob
check must be run at the model's own operating point.

#### (h) The literature regularises `C_q` with disorder, not with temperature — and here the two are the same size

Live search reached a **second, independent** compact-model lineage, the
[Southampton bilayer-GFET model](https://eprints.soton.ac.uk/366757/1/GFetIEEEtran.pdf).
It confirms §(a) from a different source than §4.6.8's
[Pasadas/Jiménez review](https://arxiv.org/pdf/2209.00388): its `C_q` is a
function of the **surface potential** `φ_s`, never of the gate overdrive, and
its layer potentials form a **coupled system solved together**. The claim that
item 16 is a self-consistency problem is now held by two independent groups
rather than inferred from one.

But it regularises `C_q` at neutrality **differently**, and the difference
opens a question this chapter has not asked. Temperature never enters its
`C_q` at all; it floors the capacitance with the residual disorder density,
`C_q,min = e²√n₀/(√π·ħ·v_F)`. The puddle density does the job this thesis
gives to `kT`.

| | `C_q` at the Dirac point (F m⁻²) |
|---|---|
| this thesis, thermal rounding | 8.437399e-03 |
| literature form, at `n_puddle` = 5e15 m⁻² | 9.710793e-03 |
| ratio | **1.1509** |

They agree exactly at `n₀ = 3.77e15 m⁻²` (0.755× `n_puddle`) or, holding `n₀`,
at **T = 345.28 K** (1.151× 300 K). The agreement is therefore not a
coincidence of one parameter: **SiO₂-supported graphene at room temperature
sits where disorder and temperature regularise `C_q` equally.**

And §4.3 applies **both**. The thermal rounding is inside `C_q`, and
`n_puddle` is then added in quadrature to `n`. The literature form shows those
are two spellings of one physical effect at this `n₀`, so **using both is a
double count of the same regularisation, worth about a factor 2 in `C_q` at the
Dirac point, and it has never been costed** (§7.9 item 18). It would not be
true for suspended graphene (`n₀` two orders lower) or at 4 K, where the two
conventions diverge — which is why neither is obviously wrong in the literature
and why the overlap here is easy to miss.

This also bounds what §(a)'s measurement can say about neutrality. The region
where the argument error would be *large* is precisely the region where two
overlapping regularisations are fighting over the same physics.

#### (g) What this section does not claim

1. It does **not** install `B_T` as the chapter's charge relation. The default
   remains (4.29); `B_T` is measured, committed and available. Flipping the
   default would move every number in Chapters 4, 5 and 6 and is its own
   session's work, with its own transcript regeneration.
   *[Annotated 2026-10-21, §4.6.17. Two corrections to this sentence, neither
   of which changes the claim that the default is unflipped. First, §4.6.12
   measured that Chapters 5 and 6 never call the charge path at all, so "every
   number in Chapters 4, 5 and 6" is wrong in that direction and "more of
   Chapter 4 than was measured" in the other. Second — and this is what
   §4.6.17 adds — "its own session's work" was read as meaning the work is
   the flip. It is not. The two halves of the flip are each correct alone and
   **do not compose**: the contact replacement reads the shipped
   `gfet.quantum_capacitance` at call time, which is the name the channel half
   rebinds, and installing both at once costs −51 % to −62 % in `n_contact`,
   4.86× to 10.77× the effect being installed. The numbers in this subsection
   stand; what was incomplete was the account of what flipping requires.]*
2. It does **not** claim `D_T` is exact physics. It is the exact integral of
   *this model's* differential relation at 300 K. The relation itself still
   uses a phenomenological rounding rather than Fermi–Dirac integrals, and
   (4.28) remains the T = 0 limit of a standard literature form.
3. It does **not** re-open §4.6.8's −0.6023 % in `I_d`. That number compared
   (4.29) against (4.28) and is reproduced unchanged; what is corrected is the
   *attribution* of the raw-charge error between the two knobs.
4. It does **not** claim the argument error is negligible where it matters.
   2.64 % at the RF bias is a statement about a bias where `C_q ≫ C_ox`. The
   quantum-capacitance-limited regime, where it would be large, is not studied
   in this chapter.

<a id="s4-6-10"></a>

*Moved from Chapter 4, §4.6.10. Original section number kept.*

### 4.6.10 Item 13 closed, and a neutrality pass that measured the wrong thing three times

The rounding factor in `quantum_capacitance()` is now
`|η| + 2·log1p(exp(−|η|))`, an identity with `ln(2(1+cosh η))` because
`2(1+cosh x) = 4cosh²(x/2)`, and finite for every real η. The committed line
is retained as a comment above it.

Machine output: `quantum_capacitance_overflow_output.txt` (5 checks, 5 passed).
Module: `graphene_quantum_capacitance_overflow.py`.

**The result, stated as a measurement rather than as neutrality.** `C_q`
itself moves by **1 ulp at 4 of 400** sampled gate biases, all within 0.093 V
of the Dirac point where η is small and the two spellings round differently in
the last bit. **No committed quantity moves by a single bit**: `n`, `R_ch`,
`I_d` at both `V_ds` of record, `f_T`, `f_max`, `g_m`, `g_ds`, `C_gs` and the
crossover lengths are bit-identical over 3603 values. Saying it that way
rather than "output-neutral" is §4.6.6's prose-is-a-detector rule applied to
this section's own summary: a 1-ulp move in the function is a real difference,
and the claim that survives is about what it reaches.

**The pass failed three times before it measured the spelling.** Each draft
re-typed the old function body to compare against:

1. `e`, `ħ` and `k_B` entered as literals. `graphene_fet_model` takes them
   from `scipy.constants`, whose `ħ = h/2π = 1.05457181764615647e-34`; the
   literal `1.054571817e-34` differs at the tenth digit, squared in the
   prefactor that is ~1.2e-9, and the pass reported **3983 of 4003** values
   differing by up to 1.1e7 ulp. Sixth instance of §4.6.7's rule (find out
   what a failed check measured), and the first in which the substituted
   quantity was a *physical constant* rather than a model parameter.
2. Constants fixed, product regrouped as `prefactor·(kT/e)·log(…)` against the
   module's `smoothed_V = (kT/e)·log(…)` then `prefactor·smoothed_V`. Floating
   multiplication is not associative; the regrouping alone moved **118 of 400**
   `C_q` values by up to 2 ulp.
3. The patch meant to fix the association silently matched nothing, so draft 2
   ran again and was very nearly read as a confirmed result.

> **A function that "restores the old spelling for comparison" is a
> re-implementation, and every detail of it not copied exactly — the
> constants, the grouping, the order of operations — becomes part of the
> difference it reports.** The repair is not a better transcription. It is to
> make the thing that changed *substitutable*: `_thermal_rounding()` now holds
> the single expression item 13 touched, the pass swaps that one function, and
> the rest of the call path is shipped code running unmodified. That form
> reproduces the direct file-against-file result exactly — 4 of 400, 1 ulp —
> where all three transcriptions did not.

**And a correction to this section's own first draft.** It claimed the two
spellings agree "to ≤ 1 ulp over η ∈ [0, 700]". That was measured on eleven
hand-picked values; on 401 it is **2 ulp**. The ≤ 1 criterion was a round
number chosen before the run — precisely what §7.9 item 15 exists to find —
and it is **not widened to 2**. For an *identity*, a ulp count cannot answer
the question that matters, so `mpmath` at 60 digits adjudicates, and the answer
is not the expected one:

| | worst relative error | mean | strictly better | strictly worse |
|---|---|---|---|---|
| committed `ln(2(1+cosh η))` | 1.793e-16 | 4.643e-17 | — | — |
| installed `|η| + 2·log1p(e^−|η|)` | **1.527e-16** | 5.162e-17 | 45 / 401 | **81 / 401** |

The installed spelling is **worse at more samples than it is better**, with a
slightly larger mean error. All it has is the smaller worst case. **Accuracy is
therefore not a reason to install it**; being finite where the other is not is
the only reason, which is the only reason item 13 ever gave. The claim stops
there.

<a id="s4-6-11"></a>

*Moved from Chapter 4, §4.6.11. Original section number kept.*

### 4.6.11 The double regularisation, costed: the two regularisers never meet in (4.29), and where they do they compose max-like

Section 7.9 **item 18** was created 2026-10-11 out of a literature reading,
not out of an internal audit, and it was ranked **first** on the next-run list
because it was the only open item that could change a *physical* claim rather
than a numeric one. Its premise was that this chapter regularises `C_q` at
neutrality with temperature *and* regularises `n` with the disorder density
`n_puddle`, that at this operating point the two are the same size, and that
applying both therefore double-counts one physical effect by about a factor 2.

Costing it closes it, and it closes by having been **wrong twice over**. The
measurement is `graphene_double_regularisation.py`; machine output
`double_regularisation_output.txt` (19 checks, 19 passed); figure
`double_regularisation.png`.

#### (a) In the committed relation the two regularisers are disjoint

The premise "this thesis applies both" is false for (4.29), and the reason is
structural rather than numerical.

`C_q`'s argument in (4.29) is the gate overdrive `dV`, and `dV` is never
floored. `n_puddle` enters only the final quadrature step
`n = √(n_eff² + n_puddle²)`, which happens *after* `C_q` has been evaluated.
So **the disorder regulariser cannot reach `C_q` at all** in the committed
form.

The converse is the sharper half. The induced charge is

```
n_eff = (C_ox·dV/e) · f(C_q),      f(C) = C/(C + C_ox) ∈ (0,1)
```

and `dV` appears as a **factor**. At `dV = 0` the product is zero whatever
`f` returns, because `f` is bounded. Measured on the shipped function with
four different `C_q` conventions substituted — the thermal rounding, the T = 0
dispersion, the Southampton disorder floor, and this thesis's own
`C_Q(n_puddle)` — `n_eff(0)` is **`0.0` bitwise in all four**, and
`carrier_density(V_dirac)` returns **exactly** `n_puddle`
(`5000000000000000.0`, difference `0.0`). So **the thermal regulariser cannot
reach the Dirac-point density at all.**

Both halves are checked against *exactly-known* values — `0.0` and
`n_puddle` — and not against a tolerance. That matters here for the reason
2026-09-17 established and §4.6.10 restated: a plausible-range check on these
quantities would have passed under either premise.

> **What the committed model does is apply exactly one regulariser to each of
> two different quantities. That is a reuse of one parameter for two jobs — a
> modelling choice to defend — and not a double count, which is an arithmetic
> error.**

Item 18 was created by comparing two *literature forms* without first asking
which quantity each term in **this** model reaches. The generalisable form:
a cross-model comparison establishes that two conventions exist; it does not
establish that one model contains both.

#### (b) The spelling that does double-count, and its size

There is one spelling in which both regularisers reach the same quantity.
Item 16's correction (form B_T, §4.6.9) makes `C_q` a function of `E_F(n)`;
place the floor **inside** that argument,

```
C_q = C_q( E_F( √(n² + n_puddle²) ) / e ),
```

and the capacitance now carries the thermal rounding and the disorder floor
at once. This is the most double-counting spelling available, and it is
arguably the *consistent* one: the channel's density of states should be
evaluated at the density the channel has.

Measured at the Dirac point, in form B_T:

| | thermal OFF (T = 0) | thermal ON |
|---|---|---|
| floor outside `C_q` | 0.000000000e+00 | 8.437398908e-03 |
| floor inside `C_q` | 1.942158579e-02 | **1.991221401e-02** |

| composition | value |
|---|---|
| both / max(alone) | **1.025262** |
| both / (sum of increments) | 0.714750 |
| thermal **on top of** the floor | **+2.5262 %** |
| floor **on top of** thermal | +135.9994 % |

The double count is **+2.5262 %** in `C_q(0)`. It is not a factor 2.

> **The "factor 2" was the ratio of the two candidate values. That ratio is
> the cost of *choosing wrongly* between the conventions. It is not the cost
> of *applying both*.**

Those are different quantities and only one of them was the item. Reading one
for the other is the 2026-10-06 backlog item — *census every criterion that is
a ratio* — in its **first instance where the ratio was not a criterion but an
estimate of a magnitude**. The item has been open five sessions as a question
about verdicts; it turns out to have been a question about magnitudes too.

And the reason the two do not add is worth stating as a rule, because it is a
property of regularisers and not of these two:

> **Two regularisers that enter through a monotone function of a quadrature
> sum compose max-like, not additively: the larger swallows the smaller.**

Here `max(8.44e-03, 1.94e-02) = 1.94e-02` and `both = 1.99e-02`, which is
2.53 % above pure max and 28.5 % *below* the additive null.

#### (c) What it costs downstream, at the model's own biases

Item 17's second requirement — run a two-knob check at the model's own
operating point, not at a reference that inflates it — applied here:

| `ΔV` | neither | thermal | disorder | both | both/max |
|---|---|---|---|---|---|
| 0 V | 5.000000e+15 | 5.000000e+15 | 5.000000e+15 | 5.000000e+15 | 1.000000 |
| 1.2 V | 5.730515e+15 | 5.732955e+15 | 5.741244e+15 | 5.741709e+15 | 1.000081 |
| 2.7 V | 8.085749e+15 | 8.087032e+15 | 8.095579e+15 | 8.096222e+15 | 1.000079 |

At the Dirac point **all four conventions agree bitwise**, because every one
of them is overwritten by the quadrature floor on the *result*. The quantity
item 18 was created to protect is the one quantity in the model that cannot
move. `σ_min` is identical to the last bit in all four
(`3.204353268e-04` S/sq): it is `n_puddle·e·μ` and nothing else, so **no
committed conductivity, resistance or on/off number in Chapters 4, 5 or 6 is
exposed to item 18 at neutrality.**

At the RF bias the whole double count is worth **+0.1527 %** in `n`, against
−1.0212 % for item 16 itself and −0.60 % for the exact charge relation. Item
18 is **real but small, and it is smaller at the operating point than at
neutrality** — the opposite of how it was ranked when it was created, and the
second time in two sessions that a cost estimate has been found to price
something other than the repair it named (§4.6.9's lesson, one level up).

#### (d) The floor's own spelling carries an exact factor 2

One finding sharpens the comparison rather than the cost, and it is the reason
yesterday's ratio should not be read as settled. The Southampton floor is

```
C_q,min = e²√n₀ / (√π·ħ·v_F)
```

while graphene's density of states, `D(E) = 2|E|/(π ħ² v_F²)`, carries the
spin–valley degeneracy `g = 4` and gives this thesis's own T = 0 capacitance

```
C_Q(n) = e² D(E_F) = 2e²√n / (√π·ħ·v_F) = 2e√n / A_F.
```

**The literature floor is exactly half this thesis's own `C_Q` at the same
density** — measured at **0.000 ulp from ½** with one consistent `ħ`. It is a
convention difference, not a numerical one. Whatever its origin — a
per-species `n₀`, a bilayer convention, a different definition of the residual
density — it is not a quantity that can be compared against this thesis's `kT`
before the factor is fixed.

The consequence for yesterday's number: the ratio **1.1509** is `kT` against
the *halved* spelling. On this thesis's own degeneracy convention it is
**2.3018**, and the claim "disorder and temperature regularise `C_q` equally
at this operating point" weakens from *the same size* to *the same order*.

**Yesterday's 1.1509 is not rewritten here, and §4.6.9(h) is annotated rather
than corrected.** Both figures are reported, because which one is right
depends on what `n₀` denotes in the source, and that is a question for the
literature rather than for arithmetic. It is now the only open item on this
list whose resolution requires reading rather than computing.

#### (e) The module's first draft failed its own substitutability check

Recorded because the check that caught it was written **before** the
measurement rather than during a neutrality pass afterwards, which is the
first time that ordering has held in this repository.

Draft 1 built the two switches on `graphene_diffusion_current_model`'s literal
`E_CHARGE` and `HBAR` while `quantum_capacitance()` uses `scipy`'s, and
reproduced the shipped `carrier_density(2.0)` only to **2375 ulp**. That is
the **eighth** instance of the 2026-10-09 rule and the second in which the
substituted quantity was a *physical constant*.

The repair is 2026-10-11's repair and not a better transcription: each switch
**substitutes one shipped name** — `gfet._thermal_rounding` for the thermal
regulariser, `gfet.n_puddle` for the disorder one — and leaves every other
operation in the shipped call path. Form A reached through the switches is
then the shipped function, and agrees with it at **0 ulp** by construction
rather than by care. One arrival counter per patched name, each checked at a
call that reads that name (§4.6.10's sixth repair).

One caveat is recorded in the output rather than smoothed over: the T → 0
validation reports a departure of *exactly* zero, which is weaker than
`1e-16` would be. At `η ~ 10⁷` the term `2·log1p(exp(−η))` underflows to
exactly `0.0`, so the two expressions reduce to the same floating-point
operations. The informative version of that check is at moderate `η`, where
§4.6.10 already adjudicates the identity against 60-digit arithmetic.

#### (f) The source was read to settle the factor, and it leaves `n₀` undefined — but its own composition law bounds the whole class of error

The one residue of (d) that arithmetic cannot reach is what `n₀` denotes in
the source. The paper was fetched and read for exactly that, and **it does not
say.** At zero surface potential it states only that "the channel has been
reported to have a charge density, `n₀`", cites two references for it, and
never calls it total, per-species, or a disorder parameter; in its resistance
section the same symbol is "the minimum charge density". **No degeneracy
factor appears anywhere in the paper, and it gives no density of states at
all** — it goes directly from `n` to `C_q`.

> **The comparison is blocked not by disagreement but by the source having
> left a symbol undefined and delegated it to a citation.** There is no stated
> convention to convert, so the factor 2 of (d) cannot be reconciled by
> reading harder; both 1.1509 and 2.3018 stay in the record and neither is
> installed.

Its fitted values are in the `10¹⁶ m⁻²` range — **1.2e16 m⁻²** for its Test A,
**2.4×** this thesis's `n_puddle`. That is a second literature anchor on the
standing `n_puddle` question, and it points the same way as §4.6.7's: this
thesis's floor may be low rather than high.

What the reading *does* settle is the composition law, and it settles it in
favour of (b). Their Eq. (3) is

```
C_q(φ_s) = (C_q,min² + 2h²) / √(C_q,min² + h²),      h = C_q,var/2
```

which returns `C_q,min` **bitwise** at `φ_s = 0`, tends to `C_q,var` to twelve
digits when that dominates, and is therefore a **smooth max and not a sum**.
Its excess over pure max is maximal exactly where the two terms are equal,
and there it is

```
sup  C_q / max(C_q,min, C_q,var)  =  1.5/√1.25  =  3/√5  =  1.341640786
```

measured at `C_var/C_min = 1.000000` over a 200 001-point sweep, against that
closed form.

Two consequences, both for item 18 and not for the model.

**The literature applies `n₀` twice itself.** Eq. (2) floors `C_q` with it and
its Drude relation `n = √(n₀² + n*²)` floors the density with it, with no
temperature anywhere. So the comparison item 18 rested on is not between a
model that double-counts and one that does not; it is between two
**allocations of regularisers to quantities**. This thesis spends `kT` on
`C_q` and `n₀` on `n`; the source spends `n₀` on both.

**And §4.6.9(h)'s intuition was right about where the worst case is and wrong
about how big it is.** A smooth-max composition is furthest from pure max
exactly when its two terms are *equal* — which is precisely the condition that
section identified as holding at this `n₀` and this `T`. That reasoning was
sound. But at that worst case the excess is `3/√5`, **+34.164 %**, and this
thesis sits at **+2.5262 %** because its two terms are 2.30:1 rather than 1:1.

> **A quadrature-composed regulariser pair cannot exceed the larger of its two
> terms by more than 34.2 %. The "factor 2" was unreachable in principle, not
> merely unmeasured.**

An upper bound on the entire class of error the item named, derived from the
source's own composition law rather than from this thesis's. It is quoted as
context for the item's estimate and **not** as a bound on this thesis's own
arithmetic: the nested spelling of (b) is a different function, and its
+2.5262 % is measured rather than bounded.

#### (g) What survives

Not a double count. What survives is the **reuse of one parameter for two
jobs**: `n₀` buys the density-of-states broadening in the literature and the
carrier floor here, and this thesis has never stated which of the two it means
`n_puddle` to be. That is the standing `n_puddle` question — §4.6.7 reached it
from the disorder floor being 2.07× Chapter 3's, and item 18 reaches it from
the opposite direction — and it is still deliberately not an item, because
changing `n_puddle` moves committed numbers in three chapters and the case for
changing it has not been made.

<a id="s4-6-12"></a>

*Moved from Chapter 4, §4.6.12. Original section number kept.*

### 4.6.12 Item 16's blast radius, measured: the contact model was outside the harness, Chapters 5 and 6 were never inside it, and Pd's residual does not survive

Section 7.9 item 16 — "flip `quantum_capacitance()`'s argument from the gate
overdrive to `E_F/e`" — has been ranked first since 2026-10-12, carrying the
note that flipping it "moves every committed number in Chapters 4, 5 and 6 and
needs its own session". This is that session, and it does not flip the
default. It measures the blast radius first, and the blast radius is not what
that sentence says. It is **smaller across the thesis and larger inside
Chapter 4**, and inside Chapter 4 it reaches a verdict.

`graphene_contact_argument_blast_radius.py`, 24/24,
`contact_argument_blast_radius_output.txt`,
`contact_argument_blast_radius.png`.

#### (a) The harness that sized item 16 could not reach half of Chapter 4

`graphene_quantum_capacitance_argument.on_consistent_argument()` rebinds three
names: `gfet.carrier_density`, `gfet.quantum_capacitance`, and
`crx.carrier_density` — the third added 2026-10-10 after exactly this class of
fault. It does not rebind `graphene_contact_doping_model.quantum_capacitance`,
which that module bound at import time on its line 50

```python
from graphene_fet_model import quantum_capacitance, sheet_conductivity, ...
```

and which `contact_edge_carrier_density()` calls on its line 105 with

```python
C_q = quantum_capacitance(dV, T=T)        # dV = -(W_metal - W_graphene)
```

— **exactly the defect item 16 names**, in a second module, against a
work-function-induced potential instead of a gate overdrive.

Measured rather than read from the source: under `on_consistent_argument()`
that binding is reached **11 times**, and **exactly zero** of those calls
arrive at the corrected function. Every cost ever quoted for item 16 —
−2.2659 % in the raw charge, −1.0212 % in `n` at 3.5 V, −0.2530 % in `I_d`
(Section 4.6.9(e)) — was therefore measured with Sections 4.5, 4.7 and 4.8
**switched out of the circuit**.

> **Nothing errored and no number moved, which is the whole difficulty: a cost
> table built through this harness reports a model that was absent as a model
> that is insensitive.**

#### (b) And Chapters 5 and 6 were never inside it either

The same counters, on nine committed summary routines —
`graphene_interconnect_model.summary_numbers()` and
`self_check_parallel_limit()`, the collection model's `summary_numbers()` and
`modeled_eqe_by_metal()`, the signed-carrier model's three tables, the
two-contact model's two — record **0 charge-path calls between them**.
Chapter 5 imports nothing from the charge model. Chapter 6 imports `mu`,
`METAL_WORK_FUNCTIONS`, `W_GRAPHENE` and `lambda_decay`: parameters, not the
charge path.

> **Item 16 is a Chapter 4 operation that reaches two models, not a
> three-chapter one that reaches one.**

#### (c) The correction is 4.7× to 12.2× larger at the contact

The contact geometry is where the argument error is least forgiving, for two
compounding reasons.

`C_interface = ε₀/3.5 Å = 2.529768e-02 F m⁻²` is **65.93×** `C_ox`, and `C_q`
is *comparable* to it rather than small against it: the series factor
`C_q/(C_q + C_s)` runs **0.555–0.806** at the contact against **0.999** in the
channel at `V_g = 3.5 V`. The channel's near-unity factor is what attenuated
item 16 by 0.451 (Section 4.6.9, §7); the contact has no such attenuation.
And the potential is a work-function difference of 0.15–1.15 V rather than a
few hundred mV of overdrive, so the gap between `ΔV` and `E_F(n)/e` is wider
to begin with — `E_F/e` lands at 0.315 V where `ΔV` is 0.620 V for Pd.

| Metal | W_F (eV) | ΔV (V) | n_contact, form A (m⁻²) | form B_T (m⁻²) | change | E_F/e (V) |
|---|---|---|---|---|---|---|
| Ti | 4.33 | +0.170 | +1.644935e16 | +1.541997e16 | **−6.2578 %** | 0.1449 |
| Cr | 4.50 | −0.000 | −0.000000e+00 | −0.000000e+00 | — (ΔW_F = 0 exactly) | 0.0000 |
| Cu | 4.65 | −0.150 | −1.380517e16 | −1.314672e16 | **−4.7696 %** | 0.1338 |
| Ni | 5.04 | −0.540 | −7.111301e16 | −6.227844e16 | **−12.4233 %** | 0.2911 |
| Au | 5.10 | −0.600 | −8.034794e16 | −7.031424e16 | **−12.4878 %** | 0.3094 |
| Pd | 5.12 | −0.620 | −8.343500e16 | −7.300861e16 | **−12.4964 %** | 0.3152 |
| Pt | 5.65 | −1.150 | −1.660634e17 | −1.463454e17 | **−11.8738 %** | 0.4463 |

Against **−1.0212 %** for the same correction in the channel.

#### (d) And it changes a verdict — the one Section 4.10 said was fragile

Section 4.7's Table reports `R_transmission = Rc,measured − R_extra` for four
metals and finds three negative, i.e. unphysical, with **only Pd positive** at
+50.949 Ω·µm. Section 4.10's condition-number audit then found that this one
positive row — the only evidence the additive decomposition might survive —
**flips sign on a 9.6 % error in `R_extra`**, while Cu's would need 81.1 %.

The argument correction moves Pd's `R_extra` by **+47.6925 %**, five times
that margin.

| Metal | Rc, measured | R_extra, form A | R_extra, form B_T | change | R_tr, A | R_tr, B_T |
|---|---|---|---|---|---|---|
| Cu | 184.0 | 971.461 | 959.236 | −1.2583 % | −787.461 | −775.236 |
| Ni | 470.0 | 830.462 | 1003.490 | +20.8352 % | −360.462 | −533.490 |
| Au | 519.0 | 609.159 | 848.120 | +39.2281 % | −90.159 | −329.120 |
| Pd | 584.0 | 533.051 | 787.276 | **+47.6925 %** | **+50.949** | **−203.276** |

> **Pd's residual does not survive the correction. Under form B_T the additive
> decomposition fails on 4 of 4 metals rather than 3 of 4, and Section 4.10's
> 9.6 % margin — computed as a tolerance to INPUT error — is consumed by a
> MODELLING error that was already in the committed code.**

Two things that do *not* move, and both were checked rather than assumed: the
metal ordering by residual (Cu < Ni < Au < Pd under both forms) and the
n-type/p-type assignment of Section 4.5 for all seven metals.

Note the direction. `R_extra` goes **up** because `n_contact` goes **down**:
a less strongly doped contact edge is a more resistive one, and the
doping-profile integral of Section 4.5 integrates `1/σ(n(x))`. Cu is the
exception at −1.2583 % because Cu's contact density (−1.38e16 m⁻²) is below
the on-state bulk reference (2.0e16 m⁻²), so reducing it moves the profile
*toward* the bulk over part of the junction.

#### (e) Why the default is still not flipped

Nothing here is an argument for keeping form A. It is an argument that the
flip is a larger operation than the item said, and that two of its
preconditions were discovered today rather than before:

1. `contact_edge_carrier_density()` must be flipped **in the same commit** as
   `carrier_density()`, or Chapter 4 will hold a corrected channel against an
   uncorrected contact and Section 4.7 will mix the two.
2. Section 4.7's text and Section 4.10's annotation both have to be rewritten,
   not re-run: the prose reasons explicitly from "only Pd gives a small,
   plausible positive residual", and that sentence becomes false.

#### (f) What this section does not claim

It does not claim form B_T is right. `E_F/e` is the argument the dispersion
requires and the argument two independent literature lineages use (Sections
4.6.9(h) and 4.6.11(f)), but the contact model's *other* approximations —
a vacuum-like 3.5 Å gap with no interface dipole layer, a single
`lambda_decay` for every metal — are untouched here and are plausibly larger
than 12 %. What is claimed is narrower and is measured: **the two forms
differ by 4.7–12.5 % at the contact, the harness that sized item 16 could not
see that difference, and the difference is large enough to move a verdict
Section 4.10 had already flagged as the fragile one.**

It also does not claim the Pd flip is the *end* of Table 4.7's story. Under
form B_T all four rows are negative, which strengthens rather than weakens
Section 4.7's stated explanation — that TLM-extracted `Rc,measured` already
absorbs part of the same near-contact region this model integrates
separately. A four-of-four failure is more consistent with a double count
than a three-of-four one was.

#### (g) Two instrument defects found on the way, both kept in place

**(i) This section's own first draft counted the wrong arrival.** It counted
calls to the *name* `cd.quantum_capacitance` and asserted they would be zero
under the patch. They are 11, whether or not anything was substituted into
that name. The question is not whether the binding is reached but which
function is behind it, and the counter has to live inside the *replacement*.
**Ninth** instance of the 2026-10-10 rule, and the second in its "counts a
different arrival from the one it gates" direction.

> **A control that counts calls to a name cannot see whether the right
> function is bound to it.**

**(ii) Item 13's overflow is closed in the shipped function and open in the
harness.** `qca.cq_thermal` is still the committed `ln(2(1+cosh η))` spelling,
deliberately retained there, and still overflows at `η = 710.4758`. A fixed
point built on it takes the *electrostatic* density as its first iterate, and
that iterate grows with the series capacitance — so the solver returns `nan`
for any `C_s` above ≈ 100 F m⁻². `C_interface` is 2.53e-02 and the worst
first-iterate `η` in this session is **19.2**, so no number here is affected.
But it is item 13's lesson a third time: *a defect unreachable from the
device's operating range is reachable from a check on the device's model.*
Form B_T is built here on the **shipped** `gfet.quantum_capacitance`
accordingly — which also removes a **1.2254e-09** constants gap (prefactor
ratio `qca/gfet = 1.000000001225439`) from the A→B_T comparison, so the two
forms differ in the argument and in nothing else by construction.

**(iii) And one tolerance chosen before the run failed.** The
large-`C_s` limit check demanded `C_s`-independence to `1e-9` and measured
`5.254e-04`. The limit is approached at **first order** in `C_q/C_s`, so no
fixed tolerance states what it claims. Replaced by the rate — the departure
must fall by ten per decade of `C_s`, measured **9.9933** and **9.9993**.
**Item 15's eighth worked instance**, and the repair is 2026-10-11's: let the
structure set the criterion.

<a id="s4-6-13"></a>

*Moved from Chapter 4, §4.6.13. Original section number kept.*

### 4.6.13 The rebind set is a runtime object: 25 stale sites, 1 that can fire, and a repair that blinded its own audit

§4.6.12 deferred item 16's flip for a stated reason — nobody knew what a
rebinding of the charge path reaches — and closed by noting that the
`cd.quantum_capacitance` miss found there had been found **by hand**, as had
the two before it. Three hand-found instances in ten days is a rate. This
section measures the exposure that rate implies, and then removes the
mechanism that produced it.

**(a) What a hand-maintained list leaves behind.** An AST pass over all 53
modules finds **149** module-level `from`-imports of repository names across
**73** distinct `(owner, name)` pairs. Against the repository's **14**
hand-listed rebinding sites (`alias.name = ...`), **25** of those bindings are
left stale. Three had been found by hand; **22 had not.** The identity classes
that matter are small and the two that item 16 needs are the smallest:

| name | sites a single `setattr` fails to reach |
|---|---|
| `gfet.carrier_density` | 1 (`crx.carrier_density`) |
| `gfet.quantum_capacitance` | 1 (`cd.quantum_capacitance`) |
| `gfet.mu` | 4 (`cd`, and Chapter 6's three photodetector modules) |
| `gfet.W` | 2 (`crx`, `cd`) |
| `gfet.n_puddle` | 2 (`cd.n_puddle`, `dr.N0_SHIPPED`) |
| `gfet.C_ox` | 2 (`cd.C_ox`, `dr.C_OX`) |
| `gfet.sheet_conductivity` | 2 (`crx`, `cd`) |

**(b) Which of them are USED, measured without entering any manager.**
Perturb **only** the holder's own copy and see whether the holder's own
numbers move. This is the decomposition that matters, and it needs no
manager — which is what makes it applicable to the ten inline `try/finally`
rebinding sites that cannot be entered at all.

| stale binding | max relative change in its holder's numbers | verdict |
|---|---|---|
| `cd.quantum_capacitance` | 9.8781 (988 %) | USED |
| `cd.n_puddle` | 35.764 % | USED |
| `crx.carrier_density` | 27.007 % | USED |
| `crx.W` | 27.007 % | USED |
| `pcm.mu` | 9.1743 % | USED |
| `psc.mu` | 5.3067 % | USED |
| `ptc.mu` | 3.9039 % | USED |
| `cd.mu`, `cd.C_ox`, `cd.W` | exactly 0 | imported, never reaches a number |

**(c) The repair, and why it is a computed class and not a longer list.** A
list is exactly what was wrong all three times, so the mechanism had to stop
being one. `from M import n` copies a **reference** at import time, so the set
of places to rebind is the set of names that at runtime hold the same
**object** — discoverable by identity, with no import graph, no AST pass and
no list. `graphene_binding_class.rebind()` computes it at entry and is now
installed at all three charge-path sites (`on_consistent_argument`,
`on_exact`, and `graphene_potential_declaration.main`, which had no second
line at all). Behaviour is unchanged, and that is **verified rather than
asserted**: all three modules' committed transcripts regenerate
**byte-identically**, and each was confirmed byte-reproducible *before* the
change so that the diff is a control and not a hope.

**(d) The class is necessary and not sufficient — item 16 needs two
replacements, not one rebinding.** `cd.quantum_capacitance` is deliberately
**not** rebound. It is reached at a contact where the series capacitance is
`C_interface = 2.529768e-02 F m⁻²`, **65.93×** `C_ox`, and the series factor
runs **0.555–0.806** against **0.999** in the channel (§4.6.12). The object
the owner holds after item 16's flip solves for `n` against `C_ox`; copying it
to the contact installs a **channel solver at the contact**, which is a
different wrong answer rather than a repair. The contact-correct form is the
one §4.6.12 measured separately (**−12.4964 %** in `n` for Pd). The exclusion
lives once, in code, with its reason attached, and `rebind()` **refuses an
exclusion without a stated reason** — because an unstated exclusion is how all
three 2026-10 faults were introduced.

**(e) And the static count over-reports, which retires a claim made earlier
the same day.** A miss needs the stale holder to be **imported in the patching
process**, and an import graph over the repository cannot see that. Measured
per rebinding site in a **fresh interpreter**, so that the audit's own imports
cannot supply the answer: of the 22 remaining static misses the holder is
**PRESENT for 1** and **ABSENT for 21**. `import
graphene_velocity_saturation_model` brings in three repository modules and
none of the four `mu` holders is among them. `graphene_potential_declaration`
does not import `contact_resistance_crossover` at all — four repository
modules enter its process and `crx` is not one of them — so the reading of
that site as "the 2026-10-10 fault unrepaired in a third module" was **wrong
about this repository while right about the shape**: there was no stale
binding to leave behind and no wrong number. The single survivor is
`blast_radius.series_capacitance` leaving `cd.C_ox` stale, and (b) already
measures `cd.C_ox` as UNUSED, so it is cosmetic too.

> **The 25 is retained in place as an upper bound rather than restated. The
> honest standing is three independent axes — PRESENT, USED, REACHED — and
> nothing currently in the repository clears all three. What the sweep found
> is a landmine count, not a defect count.**

**(f) Four instrument faults, all four kept in place.**

**(i) A repair that blinds its own audit is not a repair** — §7.5.4's "created
by a fix" class in its purest form yet. Retrofitting `rebind()` made all three
charge-path sites **invisible** to the sweep, whose only notion of a rebinding
site was `alias.name = ...`, and the validation that exists to keep
2026-10-13's finding on the record **failed while nothing was wrong**: the
manager it names had stopped being a rebinding site in the only sense the
instrument understood. The sweep now recognises both forms and reports a
class-complete site as covered by construction with its exclusions listed. The
validation is rewritten to assert the same fact in the form the repair left it
in — zero *silent* misses plus one exclusion named with a reason of substance —
and its original no-silent-miss clause is **kept**, not replaced.

**(ii) The planted-binding control failed on the instrument.** The module
scanner skipped every file whose name begins with an underscore, so the sweep
reported **zero** bindings in a file it had itself just created with one. An
underscore-prefixed module is as importable as any other; the filter was a
blind spot and not a convention. **Fourth** occurrence in this repository of a
control failing on the instrument rather than on the thing measured.

**(iii) The identity control found a blind spot in the AST rule, and the rule
was extended rather than the control dropped.** Twelve second-order
**renaming** aliases (`dr.C_OX is gfet.C_ox`, `dr.N0_SHIPPED is
gfet.n_puddle`, `psc.W_CROSS_VACUUM is cd.W_GRAPHENE`) are invisible to any
rule that reads import statements. Each is a real site. What identity yields
is not a directed edge but the **equivalence class** of names sharing one
object, and for the question "what must a patch rebind" the class is the safe
criterion and the import graph the unsafe one.

**(iv) The arrival counter carries a positive control, so that LATENT is a
measurement.** Under their own workloads both charge-path managers reach `cd`
**zero** times, which agrees with §4.6.12's finding that the contact model
sits outside the harness — but a counter that reads zero because it is broken
reads the same. The same counter reads **exactly 7** when `cd` is exercised
over its seven metals.

**(g) Validations against exactly known values, six, all passing.** The
2026-10-10 repair must show 0 misses (0, and both managers now cover it by
construction); 2026-10-13's finding must survive as a reasoned exclusion with
no silent miss (it does); a planted binding must be found exactly once by each
rule (1, 1, and 2 by identity — the import and the rename); identity-patch
neutrality over **98** wrapped sites at **0 ulp**; the usedness instrument's
own null, a perturbation by a factor of exactly 1.0, at **exactly 0**; and the
arrival counter's positive control at **7 of 7**.

What stays open is named and small: `graphene_contact_doping_model` needs its
own contact-correct replacement before item 16 can flip, and the eleven
surviving hand-listed sites outside the charge path are **not** converted,
because each is latent by import today and converting one changes nothing
until it is not.

<a id="s4-6-14"></a>

*Moved from Chapter 4, §4.6.14. Original section number kept.*

### 4.6.14 The contact gets a replacement instead of an exclusion, and the identity class turns out to span two operating points

Section 7.9 item 16 — flip `quantum_capacitance()`'s argument from the gate
overdrive to `E_F/e` — has been the first-ranked open item on this chapter
since 2026-10-10. §4.6.13 reduced its blockers to exactly two, of which the
first was:

> `graphene_contact_doping_model` needs its **own** contact-correct
> replacement, because the contact cannot take the channel's object. The flip
> is two replacements, not one rebinding.

**That blocker is now discharged.** What follows is the discharge, and the
structural finding that came out of attempting it. Neither flips the default:
item 16 remains open, and what has been removed is the reason it could not be
flipped.

#### (a) A measurement harness is not a replacement

The contact-correct numbers were not unknown. §4.6.12 measured them a
fortnight before, on 2026-10-13: the correction at the contact edge runs
**−4.7696 %** (Cu) to **−12.4964 %** (Pd) in `n`, against **−1.0212 %** for
the same correction in the channel at `V_g = 3.5 V`. The contact is **4.7× to
12.2× less forgiving**, because its series factor `C_q/(C_q + C_s)` runs
**0.555–0.806** where the channel's is **0.999**, and it has a wider `dV` to
be wrong about.

Having the numbers was not having the object, for three reasons, and only the
third was a surprise:

1. **The object belonged to the audit.** It lived inside
   `graphene_contact_argument_blast_radius`, so nothing outside that module
   could install it. Item 16's flip had no second object to install even
   though its size was known to four decimal places.
2. **The harness rebound by hand.** A single assignment,
   `cd.quantum_capacitance = ...` — precisely the mechanism §4.6.13 retired
   for the channel after measuring that hand-listed rebinding leaves **25**
   stale sites at **14** listed sites, **22** of which no audit had ever
   found. The contact half was still on the retired mechanism.
3. **And the replacement mechanism could not express what the contact needed.**

#### (b) Identity answers *which names move together* and does not answer *to what*

`bc.rebind()` — §4.6.13's repair — computes a name's runtime **identity
class**, the set of names that hold the same object, and gives every member
**the same replacement**. Its one escape hatch is `CONTACT_OPERATING_POINT`:
leave this member alone, with a reason recorded. So its vocabulary is
`{same replacement, left alone}`, and it has no word for *rebind this member
to a **different** object* — which is exactly what the contact requires the
moment it has an object of its own.

Measured at call time, the class is:

| member | operating point | series capacitance | series factor |
|---|---|---|---|
| `graphene_fet_model.quantum_capacitance` | channel | `C_ox` | 0.999 |
| `graphene_contact_doping_model.quantum_capacitance` | contact | `C_interface` = 65.93 × `C_ox` | 0.555–0.806 |

> **Two names can hold one object and still sit at two operating points. An
> operating point is not something identity can see, so the class must be
> PARTITIONED by one — and the partition has to come from the caller.**

This is a limitation of §4.6.13's repair discovered by trying to use it, not a
defect in it. The repair's claim — that the set of names to rebind is
computable and needs no list — survives intact; what fails is the implicit
assumption that one class takes one replacement.

`bc.rebind_by_operating_point()` is that generalisation. `bc.rebind()` is
**left untouched**, and all three transcripts it owns regenerate
byte-identically. The new function takes a map from class member to
`(replacement, reason)`, where `replacement = None` expresses §4.6.13's
exclusion in the same vocabulary rather than as a separate mechanism, and the
`reason` is mandatory in both cases.

The refusal carries across in the form the same lesson takes here. §4.6.13
refuses an **exclusion without a stated reason**, because an unstated
exclusion is how all three of the 2026-10 delivery faults were introduced.
This refuses an **unassigned class member**:

> **Silently leaving a member at the shipped object is the identical fault,
> reached by omission from a map instead of omission from a list. The whole
> value of computing the class is that the omission is detectable, so it is
> detected: the manager raises and rebinds nothing.**

#### (c) The emptiness that discharges the precondition

`on_consistent_argument_everywhere()` installs the charge relation and **both**
capacitance operating points together, and its left-alone set is **empty** —
four sites installed, zero left at the shipped object:

```
graphene_fet_model.carrier_density                     channel
contact_resistance_crossover.carrier_density           channel
graphene_fet_model.quantum_capacitance                 channel  (C_ox)
graphene_contact_doping_model.quantum_capacitance      contact  (C_interface)
```

The emptiness is **asserted by the mechanism, not by this section**: the
manager refuses to open if any member of either class is unassigned, and a
negative control plants a missing member and requires the refusal to fire.
Without that control, an empty left-alone set and a mechanism that never
checks would read identically — the 2026-10-14 rule that a null result is a
measurement only when the instrument carries a control.

#### (d) Validations against exactly known values

Nine checks, 9/9, five of them against values known exactly rather than
plausible ranges:

| check | value required | result |
|---|---|---|
| `ΔWF = 0` → `n = 0` under **both** forms | exactly 0 | `−0.0` / `−0.0` |
| form `B_T` odd in `dV` | **0 ulp** | 0 ulp over 5 offsets |
| new machinery fed the **OLD** argument vs shipped | bitwise | **0 ulp**, all 7 metals |
| quantum limit `n = 4(dV/A_F)²`, error order in `1/C_s` | ratio exactly 100 | **100.000** at `C_s × 100` |
| today's owned object vs §4.6.12's independent construction | bitwise | **0 ulp**, all 7 metals |
| planted unassigned member | must raise | raises |
| three arrival counters | positive control on each | 1 / 19 / 1 |

The third is the one that makes the rest mean anything. Feeding the new code
path a `C_q` that ignores `n` and returns the shipped function at the **old**
argument reproduces the shipped `contact_edge_carrier_density()` **bitwise**
over all seven metals. That is what establishes that form A and form `B_T` at
the contact differ in `C_q`'s **argument** and in nothing else — no second
solver, no second set of constants, no altered series factor. Without it,
every per-metal number in §4.6.12 and here could be the residue of some other
difference between two implementations. The fourth is the honest form of a
closed-form check: at finite `C_s` the solve differs from the limit by
`O(C_q/C_s)`, so a single threshold would be a plausible-range check wearing a
closed form's clothes, and the exact statement is the **convergence order**.

#### (e) And the first draft of one check misattributed a real defect

The oddness check (second row) was first written through the shipped
`contact_edge_carrier_density()` rather than on the solver, and **failed at 3
ulp**. It then failed at **2–3 ulp under form A as well**, which is what
identified the cause:

> **The shipped function's own `dWF = work_function_metal - W_GRAPHENE` is not
> antisymmetric in double precision.** `(4.5 + 0.05) − 4.5 =
> 0.050000000000000266`, while `(4.5 − 0.05) − 4.5 = −0.04999999999999982`.
> The solver is exactly odd; the work-function subtraction is not.

So a pre-existing 2–3 ulp asymmetry in the shipped contact model was briefly
read as a defect introduced by the new solver. It is kept as its own
validation row rather than quietly dropped, because this is §4.6.6's rule for
the second time in two sessions: **a check that reads two knobs at once
attributes the effect to the wrong one.** The defect is cosmetic at 2–3 ulp
against a 12 % physical effect, and it is now on the record instead of being
unknown.

#### (f) What this does and does not settle

| | standing |
|---|---|
| Precondition (i) — the contact's own object | **DISCHARGED** |
| Precondition (ii) — §4.7's prose and §4.10's annotation | **OPEN**; writing, not plumbing |
| Item 16 itself — flipping the default | **OPEN**; nothing here flips it |

§4.7's prose becomes false under the flip — "only Pd gives a small, plausible
positive residual" cannot survive a 12.4964 % correction to Pd's `n_contact`
— and §4.10's conditioning annotation inherits that. Both are rewritten under
item 16 proper, not here, and §4.6.12's numbers stay annotated in place rather
than replaced.

Module: `graphene_contact_consistent_argument.py`. Machine output:
`contact_consistent_argument_output.txt`. Figure:
`contact_consistent_argument.png`. Mechanism:
`graphene_binding_class.rebind_by_operating_point()`.

<a id="s4-6-15"></a>

*Moved from Chapter 4, §4.6.15. Original section number kept.*

### 4.6.15 The inert-knob census was counting a syntactic pattern for a semantic fault, and three instruments had been dark for twelve days

Sections 4.6.13 and 4.6.14 built the machinery item 16's flip needs: an
identity class computed at run time, and a partition of that class by
operating point. Both deliver a change the same way everything in this
repository does — by **rebinding a name**. That makes one question
load-bearing for every result in Chapters 4 to 6 that was ever established
by mutation: *which readers does a rebinding actually reach?*

Since 2026-09-28 the known answer has had one documented exception, the
**inert knob**:

```
K = 1.0
def f(x=K): ...          # x's default is bound when `def` executes
```

A rebinding of `K` never reaches `f`, so a mutation audit that rebinds `K`
and watches `f` measures nothing and reports a pass. Three separate
instruments census the pattern —
`graphene_figure_provenance_audit.census()` (Section C),
`graphene_mutation_arrival_probe.repo_wide_frozen_census()` (Section 7b),
and `graphene_cross_module_delivery_audit` for the cross-module subset.

#### (a) All three transcripts were written on the same day, and none of them has been read since

The committed transcripts read **0**, **0** (an empty census behind a
`[PASS]` on "no module-level constant is captured as a default anywhere")
and **8**. All three were last written **2026-10-04**. Run today they read
**21**, **20** and **16**.

> **A fault class that 2026-10-04 drove to zero regrew over twelve days in
> front of three censuses, none of which was being read — and the fourth
> instrument whose only purpose is to notice exactly that
> (`graphene_pristine_transcript_audit.py`, "A COMMITTED TRANSCRIPT IS A
> CLAIM ABOUT THE CODE BESIDE IT") was last run on the same day.** Asked
> today, unprompted, it reported four of its subjects STALE or ERROR.

This is §7.7's unequal-scrutiny habit one level out. The repository has
spent three weeks learning that a guard which is not run is not a guard;
what it had not noticed is that the rule applies to the guard that checks
the guards.

#### (b) 21 against 20 is not a rounding difference

Over one path set both censuses name the **same five files**, and **four of
the five agree exactly**. They differ by one occurrence, in
`graphene_velocity_saturation_model.py` at line 589 —
`_R_harmonic(V_g, V_hi, n_quad=N_QUAD)`, defined **inside the body of
another function**. Measured rather than asserted, in all three
configurations, against exactly known values:

| where the `def` is | when the `def` statement executes | captured default |
|---|---|---|
| module body | once, at import | **FROZEN** |
| class body | once, at import | **FROZEN** |
| function body | on **every call** of that function | **LIVE** |

With `K` rebound 1.0 → 99.0: the module-level arm returns exactly 1.0, the
class-body arm exactly 1.0, the function-nested arm exactly 99.0.

> **THE SYNTACTIC PATTERN `def f(x=K)` IS TWO CLASSES WITH OPPOSITE
> BEHAVIOUR UNDER THE ONE OPERATION THIS REPOSITORY USES TO DELIVER
> CHANGES. A FUNCTION-NESTED CAPTURE IS NOT AN INERT KNOB. IT IS LIVE.**

So neither census is simply right:

- **figure_provenance's 21** counts the syntax correctly and over-counts the
  **fault** by exactly the live captures. Its own historical figures — "49
  occurrences across 14 files, 14 of them inside the audit modules … 37
  remain, and the 13 in audit modules are the next run's top item" — are
  syntactic counts of a class that is not uniformly the fault, and have
  never been split into the two.
- **mutation_arrival_probe's 20** counts the fault correctly **for the wrong
  reason**. Its `classify()` confirmed each statically-found site with
  `getattr(module, fname)` — added correctly, to drop a function a later
  redefinition or decorator had replaced — and a function-nested function is
  not a module attribute, so the site was dropped **silently**: no row, no
  class, no reason. The right verdict by an accident of resolution.

The accident does not generalise. The same line also drops a **class-body**
capture, which is *frozen* — so the next capture of that shape would have
been **undercounted**, and undercounting is the dangerous direction here: an
instrument that silently omits a true site has the identical signature to the
fault being censused. The repository contains no class-body capture today;
that is luck, not a property of the mechanism.

#### (c) Reconciled to zero residue

```
21 syntactic  =  20 FROZEN (the genuine inert knobs)  +  1 LIVE
```

Four checks, because a reconciliation that leaves a remainder means a third
cause is present and unaccounted for: the arbiter reproduces
figure_provenance's total exactly (R1), mutation_arrival_probe's confirmed
total exactly (R2), accounts for the whole gap with nothing left over (R3),
and agrees **per `(file, constant)`** rather than only in total (R4 — added
because two totals can agree while the rows behind them disagree in
compensating directions; it is the only check that catches a single dropped
row).

Per file: `graphene_velocity_saturation_model` 9 = 8 + 1,
`graphene_diffusion_current_model` 5 = 5 + 0,
`graphene_exact_charge_rederivation` 3 = 3 + 0,
`graphene_perfect_contact_counterfactual` 2 = 2 + 0,
`graphene_second_anchor_bound` 2 = 2 + 0 — the last of these added by the
session immediately before this one.

#### (d) The repairs, and the one deliberately not made

`classify()` now takes the verdict from the **scope** instead of discarding
what it cannot resolve: function-enclosed → `live_by_scope`, otherwise →
`frozen_unresolved`, which still counts as frozen because the static evidence
stands and `getattr` merely cannot reach it. The fix is a verdict from the
scope and **not a widening of the resolver**, because resolvability was never
the question being asked. The amendment is **count-neutral** — the frozen
total is 20 before and after, and the suite's own output differs in exactly
the two lines that now carry `1 rescoped-live`. The `getattr` confirmation's
original purpose is preserved, not removed: a module attribute rebound after
import genuinely no longer carries the captured default (V6).

The pristine audit needed its own repair, created by this session's decision
rather than inherited. Every suite here ends `return 0 if failed == 0 else
1`, so **a red guard is a nonzero exit by design** — and `run_suite`
discarded stdout on any nonzero exit and reported `ERROR` without comparing
the transcript at all. Tolerable while every suite was green; not tolerable
once this repository began committing red transcripts in preference to stale
green ones, since two suites would then have had transcripts nothing checks,
permanently. ERROR is now reserved for a run producing no usable stdout, and
a new check requires a transcript's own `TOTAL: N passed, M failed` line to
agree with its suite's exit code — the two being independent claims about one
run, which contradicted each other for the whole of the twelve days.

**Not done here, and the reason matters more than the omission.** The 20
genuine inert knobs are not converted. Converting 21 captures as one class
would have rewritten a live capture into a `None`-sentinel for no reason and
recorded it as a fix; the work is ranked first for the next run and is now
**specified correctly for the first time** — 20, not 21, with the live one
named so it is not converted by mistake.

#### (e) A defect in this section's own instrument, kept on the record

The first draft of `graphene_frozen_default_semantics.py` classified by
**nesting depth**: depth 0 frozen, depth ≥ 1 live. Its own control planted a
method in a class body, called it depth 1, and asserted *in prose* that it
was therefore not an inert knob. A class body executes once at import, so
that is wrong — and the repository contains no class-body capture, so the
draft's count would have been right **by luck**, which is the same luck
`classify()` was running on. That is what made it worth recording rather than
quietly fixing. Controls C3/C4/C5 now pin all three cases, including a class
defined *inside* a function, where a depth rule and a class rule both get it
wrong.

#### (f) What this means for item 16

The flip's blast radius is established by regenerating committed transcripts
and reading what moved. The instrument of record for that is the pristine
transcript audit — which, at the start of this session, reported four of its
own subjects STALE or ERROR.

> **A flip measured against that baseline could not have distinguished its own
> blast radius from twelve days of accumulated drift. That is a precondition
> (iii) for item 16 that no entry had listed, and it is discharged here: 13
> IDENTICAL, 0 STALE, 0 NUMERIC_DRIFT, 0 ERROR.**

Module: `graphene_frozen_default_semantics.py` (**24/24**). Machine output:
`frozen_default_semantics_output.txt`. Amended:
`graphene_mutation_arrival_probe.py`,
`graphene_pristine_transcript_audit.py`. Mutation-tested outside the
repository, so the instrument and the measurement do not move together: the
depth rule restored fails C1–C4; every capture called frozen fails C1, C2,
C5, R2, R3, R4; one frozen row dropped fails only R2 and R4; `run_suite`
reverted fails C7; and the 2026-10-04 transcript restored verbatim over
today's red guard fails the new consistency check **by name**.

<a id="s4-6-16"></a>

*Moved from Chapter 4, §4.6.16. Original section number kept.*

### 4.6.16 The list that ranks this work had never been audited, and two of its items asserted the opposite of the tree

Every session in this project chooses what to do from the "Not yet covered
(candidates for future runs)" list at the foot of the preceding
`AUTOMATION_LOG.md` entry. §4.6.12 audited a harness, §4.6.13 a rebind set,
§4.6.15 three censuses. Nothing had ever audited **the list**, and a claim
nothing checks is the standing subject of this chapter's instrument
sections — §4.6.15's finding was that three instruments had been dark for
twelve days, and the list is the artefact with the widest reach of all of
them, because it is upstream of which instrument gets read at all.

Two of the 2026-10-17 list's twenty-five items assert a state of this tree
that the tree denies.

**A2, item 16's precondition (ii).** Recorded 2026-10-17 as "§4.7's prose and
§4.10's annotation — is still untouched and is still writing, not plumbing",
and recorded in §7.9 on 2026-10-14 as "the only precondition still open".
Commit `7604c05` of **2026-10-13** is titled *"Section 4.6.12, and the two
sentences in Sections 4.7 and 4.10 that the measurement makes false, annotated
in place"*. Both annotations are in this chapter, above §4.7's table and below
it, carrying −1.2583 % (Cu), +20.8352 % (Ni), +39.2281 % (Au), +47.6925 %
(Pd), and `+51.0 → −203.276`. The §7.9 sentence was therefore written **one
day after** the work it describes as outstanding, and carried for four days
after that.

What is genuinely open is narrower than the item said and **cannot be done
ahead of the flip**, which is the part worth recording. Both annotations are
written in the subjunctive: §4.7's says the sentence it annotates is "true of
the committed model and false of the corrected one" and leaves it standing
"pending the item 16 flip". The writing that remains is the writing only the
flip makes possible — the point at which "pending the item 16 flip" is itself
stale, the committed table becomes the superseded column rather than the live
one, and §4.10's margin-against-input-error becomes a statement about a margin
already consumed. A session that believes (ii) is missing can discharge it by
writing a sentence that asserts a state of the code that does not exist: the
§7.5.4 "created by a fix" class, in prose. **Item 16's remaining content is
the flip and nothing else.**

**A1, and it is a revival rather than a staleness.** The 2026-10-17 list
carries *"CHAPTERS 2 AND 3 ARE STILL UNDRAFTED — their computational results
have been complete for weeks. This is the largest remaining block of pure
writing and has now gone three sessions without being touched"*, and closes it
*"Worth ranking against item 16 deliberately rather than by default"*. Chapter
2 was drafted 2026-09-28 (`5bbd5f3`, eleven sections) and Chapter 3 on
2026-09-29 (`c1d5abc`, ten sections); the tree holds 375 and 330 non-blank
lines and 11 and 9 `##` sections. So the item offered as a competitor for
first rank was describing writing committed three weeks earlier.

And the list had already got this right. Measured across every dated entry:

| wording | entries |
|---|---|
| "Chapters 2–3 remain/are still undrafted" | 2026-09-17 … 2026-09-27 (12), **and 2026-10-17** |
| "Chapters 2 and 3 **are drafted**, but §4.8.2's form of error is not audited for elsewhere" | 2026-10-01 … 2026-10-06 (6) |
| neither | 2026-10-07 … 2026-10-16 (10) |

Every "undrafted" entry up to 09-27 was **true when written**. The list then
tracked the drafting correctly for six entries, carried nothing for ten, and
the entry after those ten revived the retired wording. The two failure modes
have different cures and are worth separating: *staleness* is cured by
re-checking an item, *revival* by making retirement checkable. An item absent
long enough stops being remembered as retired and becomes available to be
re-derived.

**The cause is outside the repository.** These sessions are driven by a
scheduled brief carrying a frozen description of repository state dated
2026-09-17, which contains the A1 sentence almost verbatim — eleven entries
after the list retired it. No audit inside the tree can read that brief, so
the cause is recorded in `graphene_standing_list_audit.py` as a
`NOT_MECHANISABLE` row with its reason rather than silently omitted, and the
compensating rule is procedural: **a claim about the state of the tree that
arrives from outside the tree is re-verified against the tree before it ranks
anything.**

**The instrument, and why its retired rows are the load-bearing half.**
`graphene_standing_list_audit.py` parses the latest entry's list, matches each
baseline assertion, and evaluates a predicate over the tree — so a verdict is
gated on the mechanism it reports rather than on text similarity, which is
Design-Verification-Roadmap's 2026-10-16 rule. The trap this file had to avoid
is §4.6.15's: a suite that erases its own finding the moment the finding is
acted on. Once A1 and A2 are retired from the list, "absent" must not read as
a pass, so a **retired** row asserts two things — that the wording stays gone
(returning it reports `REVIVED`, which is a fault class, and is the check that
would have caught 2026-10-17 on the day) and that the tree-side fact which
justified the retirement still holds (`RETIREMENT-STALE` otherwise, so
truncating Chapter 2 later would re-open the item rather than leave a silent
retirement). `ABSENT` is also a fault class, so a list the parser cannot read
cannot pass as clean.

Two by-products. **The module's first draft failed its own plant, on its own
subject matter.** `p_item16_ii_untouched` tested `phrase in file.read()` for a
sentence this file *wraps*, and for numbers written here with a Unicode minus;
it therefore reported §4.10 "untouched" about the file carrying the amendment
— it **confirmed the item it was written to contradict**, by answering "are
these exact bytes present" instead of "is this claim present". That is the
same defect as §4.6.15's two docstrings, which passed an audit by spelling
while being the thing it audited for, and as the 2026-10-17 grep for a
`setattr` route over computed names. Three for three: this repository's
instruments fail first by answering a weaker question than the one asked, and
the fix each time is to normalise or to put the question to an oracle rather
than to a string.

**And the instrument already existed, in the other repository.**
Design-Verification-Roadmap built `tools/open_item_staleness_audit.py` for
this failure on 2026-10-13, and on 2026-10-17 used it to retire an item for
the identical reason — *"it names work finished before the item started being
carried"*. That repository caught its instance on the same day this one
revived its instance. The difference between the two days was one module; the
generalisation is that an instrument built in one place does not protect the
other, which is §7.5.2 across repositories rather than across modules.

<a id="s4-6-17"></a>

*Moved from Chapter 4, §4.6.17. Original section number kept.*

### 4.6.17 The two halves of item 16 do not compose: a replacement that reads a name the same flip rebinds

§4.6.9 built form `B_T` in the channel. §4.6.14 gave the contact a
replacement of its own instead of an exclusion. §4.6.12 measured what each
costs. The log of 2026-10-20 recorded, on that basis, that item 16's
"remaining content is the flip and nothing else".

**It was not.** Composition was left as well, and composition fails.

#### (a) The measurement

`graphene_default_flip.py` installs both halves at once and asks what
`n_contact` becomes. Against the committed contact-only form `B_T`:

| Metal | `n_contact`, form A (m⁻²) | committed `B_T` (m⁻²) | both halves, unrepaired (m⁻²) | `B_T` vs A | unrepaired vs `B_T` | ratio |
|---|---|---|---|---|---|---|
| Ti | 1.644935e16 | 1.541997e16 | 7.281739e15 | −6.2578 % | −52.7772 % | 8.43× |
| Cr | −0.0 | −0.0 | −0.0 | `n = 0` | `n = 0` | — |
| Cu | −1.380517e16 | −1.314672e16 | −6.391282e15 | −4.7696 % | −51.3850 % | 10.77× |
| Ni | −7.111301e16 | −6.227844e16 | −2.466546e16 | −12.4233 % | −60.3949 % | 4.86× |
| Au | −8.034794e16 | −7.031424e16 | −2.760960e16 | −12.4878 % | −60.7340 % | 4.86× |
| Pd | −8.343500e16 | −7.300861e16 | −2.859729e16 | −12.4964 % | −60.8302 % | 4.87× |
| Pt | −1.660634e17 | −1.463454e17 | −5.573671e16 | −11.8738 % | −61.9143 % | 5.21× |

The unrepaired composition overshoots the effect item 16 exists to install by
**4.86× to 10.77×, metal by metal, and in the same direction** — so it would
not have presented as a sign error, and a session that landed the flip and
regenerated transcripts would have committed it as the answer.

> The per-metal ratio is the honest statement, and the first draft of this
> table did not give it. It divided the smallest-magnitude unrepaired
> deviation by the **largest**-magnitude committed effect and the other way
> round, quoting "4.1× to 13.0×" from two pairs of numbers belonging to
> different metals. That is this chapter's recorded dominant failure mode — a
> claim answering a weaker question than the one it was written to ask — and
> it was caught by *reading the transcript*, not by any check in the file.

#### (b) The cause is one attribute read

The contact replacement is deliberately built on the **shipped**
`gfet.quantum_capacitance`, for the two reasons §4.6.12 records: the
alternative carries its own constants and differs by 9.48e6 ulp, and the
committed `log(2(1+cosh η))` spelling overflows where the fixed point's first
iterate lands. It reads that name by attribute lookup at call time:

```python
return gfet.quantum_capacitance(qca.A_FERMI * np.sqrt(nn), T=T)
```

and the channel half of the flip rebinds exactly that name. With both halves
installed, the contact's self-consistency is solved against a `C_q` that has
itself already been solved — once, against the wrong series capacitance.

#### (c) The general form, which no mechanism in this thesis covers

§4.6.13 established that the set a patch must reach is the name's runtime
**identity class**, so that no binding is left stale. §4.6.14 established that
the class must be **partitioned** by operating point, so that no member is
silently left at the shipped object. Both are statements about *which names
move*.

**Neither asks whether a replacement READS a name the same flip rebinds.**

A replacement that does is neither stale nor unassigned. It is correctly
installed, computed by the right mechanism, reported as installed by that
mechanism's own report — and it reads a different object than the one it was
validated against. This is a third failure mode beside the two §4.6.13 and
§4.6.14 were built for, and it is invisible to both. The audit that would
catch it without running anything is a one-line question asked of each
replacement: *which module-level names does this function read at call time,
and does the flip rebind any of them?*

> **Annotated 2026-10-22 (§4.6.18). The question as worded above answers
> *no* for this very replacement.** It says "*the* flip" — the flip that
> installs the replacement. `blast.contact_on_consistent_argument()`
> installs `_cq_consistent_contact` and rebinds `cd.quantum_capacitance`;
> the replacement reads `gfet.quantum_capacitance`. The intersection is
> **empty**. The defect is **cross-flip**: the name it reads is rebound by a
> *different* manager, which is only installed at the same time because
> item 16 composes them. The audit of §4.6.18 measures both forms and
> reports 1-of-4 under the composition against 0-of-4 under the wording
> above. The sentence is left unrevised because it is the record of what was
> asked; the correct question is *every flip installable concurrently*.

#### (d) Two mechanisms, and they partially cancel

The failure has two parts, not one. **(i)** the inner `C_q` is the channel's
form-`B_T` solve rather than the shipped form-A `C_q`; **(ii)** that inner
solve reaches `qca.series_factor`, which reads `gfet.C_ox` at call time, and
the contact solve has `C_ox` substituted by `C_interface` for its duration —
so the inner channel solve runs against the **contact's** series capacitance
too. Isolating (i) by restoring `C_ox` for the inner call:

| Metal | (i) alone | (i)+(ii) | (ii) contributes |
|---|---|---|---|
| Ti | −53.8545 % | −52.7772 % | +1.0772 |
| Cu | −52.4054 % | −51.3850 % | +1.0204 |
| Ni | −62.1287 % | −60.3949 % | +1.7338 |
| Au | −62.5403 % | −60.7340 % | +1.8063 |
| Pd | −62.6595 % | −60.8302 % | +1.8293 |
| Pt | −64.2107 % | −61.9143 % | +2.2965 |

**(ii) offsets (i) for every metal.** Repairing (ii) alone — which is the more
obviously wrong-looking of the two, since a global being mutated across a
nested call is visible on inspection — would have made the error *larger*.
The repair has to be made at (i).

> **Annotated 2026-10-22 (§4.6.18). (i) and (ii) are not two defects.** A
> direct count of how often the channel replacement is reached during the
> contact census, and of the value `gfet.C_ox` holds at each arrival, gives
> **192 arrivals under the unrepaired flip, 181 of them with `C_ox` already
> substituted by `C_interface`** — and **0 arrivals under the repaired
> flip**. (ii) cannot occur without (i), because (ii) is a property of the
> context (i)'s arrivals land in. One cause, two consequences, 181/192 of
> the arrivals carrying both at once. The percentages in the table are
> unrevised and remain the correct decomposition *of the error*; what is
> corrected is the reading that they are two independent faults. It also
> explains why a single change at (i) repairs both, which the table above
> could not have predicted.

#### (e) The repair, and what it does not do

Bind the contact replacement to the shipped object **captured at manager
entry** rather than to the live name. `on_default_flip()` then installs the
whole flip through the two managers §4.6.14 named — `bc.rebind()` for the
charge, whose class is `{gfet, crx}` and both members in the channel, and
`bc.rebind_by_operating_point()` for the capacitance, whose class spans the
two operating points. This is the first caller to give that manager a
**non-`None` contact assignment**; until §4.6.14 the contact had no
replacement of its own and the exclusion was correct exactly as long as that
was true.

Validated against exactly known values, not plausible ranges: the identity
assignment (every member → its own shipped object) is **bitwise** neutral on
all nineteen committed numbers; the flip's contact half reproduces the
committed contact-only form `B_T` **bitwise** on all seven metals; its channel
half reproduces `on_consistent_argument()` **bitwise** on all eight census
entries; `n_contact(ΔW_F = 0)` is exactly `0.0`; the electron/hole antisymmetry
is bitwise; an incomplete assignment raises and rebinds **nothing**; and a
positive control asserts the repair is load-bearing, so that if the two ever
agree this section's premise is known to have expired. 18 passed, 0 failed;
transcript byte-identical across two consecutive runs.

**With both halves in the circuit for the first time, every cost quoted in
Chapter 7 item 16 is reproduced exactly** — `n` −1.0212 %, `I_d` −0.2530 %,
`n_contact` −4.7696 % to −12.4964 %, §4.7 `R_extra` −1.2583 % to +47.6925 %,
Pd's implied residual +50.9489 → −203.2764 Ω·µm, additive decomposition
4-of-4 negative. That reproduction is the evidence that the defect lives in
the composition and in neither half.

**This section does not flip the default.** `graphene_fet_model.py` is
unedited and the chapter's charge relation is still (4.29). What has changed
is that the flip is a single installable, validated, costed object, so what
remains is a decision about what this thesis claims rather than a build.

> **Annotated 2026-10-22 (§4.6.18). The repair is correct and its
> correctness is contingent.** Binding to the object captured at entry
> protects the **name** the cell holds. It does not protect the names the
> **held object** reads at call time. It is sufficient here only because the
> shipped `gfet.quantum_capacitance` reads nothing in the rebind set — its
> own live reads are `{gfet._thermal_rounding}`. Its neighbour
> `gfet.carrier_density` reads `gfet.quantum_capacitance` at call time, so
> the identical construction applied to *it* is **not** immune, which
> §4.6.18 demonstrates with a planted control rather than arguing. Had the
> contact needed the shipped charge relation instead of the shipped
> capacitance, the same repair would have been written, would have looked
> right, and would have left the defect in place. Separately, the live read
> is still **in** the function on the `live=True` branch, and
> `make_contact_cq` is public: this is a default, not a removal.

<a id="s4-6-18"></a>

*Moved from Chapter 4, §4.6.18. Original section number kept.*

### 4.6.18 Auditing what a replacement reads: three readers, a cross-flip relation, and the depth of a capture

§4.6.17(c) named the audit that would have caught its own defect without
running anything, and the log of 2026-10-21 ranked it first precisely because
that defect had been found by accident. This section builds it. The result is
that the exposure §4.6.17 found is, as framed there, **complete** — and that
the instrument built to confirm that found two things §4.6.17 did not ask
about, one of which bears on its own repair.

#### (a) Three readers, chosen so that agreement is evidence

A replacement's dependence on a name is a dependence on the **live binding**,
so the question is which `(module, attribute)` pairs it resolves at call time.

**Reader A — static.** A bytecode walk of the transitive closure of live name
lookups: every `LOAD_GLOBAL` of a repository module followed by a `LOAD_ATTR`,
every repository function reached through such a name or through a global,
followed into recursively. Closure cells are resolved from the **real function
object**, so a *captured* object is distinguished from a *live name read* —
which is the entire content of §4.6.17's repair. Branch-**insensitive**, hence
an over-approximation; §(d) turns that into a result rather than a caveat.

**Reader B — sentinels.** An identity counting wrapper is installed at each
member of the flip's rebind set, the replacement is called, and the counters
that fired are the names it read. Object substitution, so it sees a read
whether it is spelled `mod.attr` or held as a module global by
`from mod import attr` — the spelling behind all three of the October delivery
faults.

**Reader C — module-type interception.** Each module's `__class__` is swapped
for a `ModuleType` subclass whose `__getattribute__` and `__setattr__` record
the access and delegate. It sees two things B cannot: reads of
**non-callables** (`gfet.C_ox` is a float and cannot be wrapped) and
**writes**.

B and C are blind in **opposite** directions — a read made from inside a
name's own module is a globals-dict lookup and never reaches the module type —
so their agreeing is a test and not a restatement. They agree exactly on all
five replacements. The planted controls confirm the blind spot is real: three
of them are flagged by B and invisible to C. Rebind sets are obtained by
entering each manager under reader C and collecting its own `setattr` calls,
so a manager that grows a name is reported with it.

Both instruments are bitwise neutral on the seven committed contact densities,
checked before any count is used.

#### (b) The answer, and the question's wording

Of the four replacements §4.6.17 names, **exactly one** reads a name that is
rebound while it is installed: `blast._cq_consistent_contact`. The accident of
2026-10-21 found the whole of the exposure as it framed it.

Asked in §4.6.17(c)'s own words, however — *does the flip **that installs
it*** rebind any of them — the answer for that same replacement is **empty**.
The defect is cross-flip, and the table makes the two readings explicit:

| Replacement | own flip | the composition |
|---|---|---|
| `qca._carrier_density_consistent` | — | — |
| `qca._quantum_capacitance_consistent` | — | — |
| `blast._cq_consistent_contact` | **—** | **`gfet.quantum_capacitance`** |
| `flip.make_contact_cq(live=False)` | — | — |
| `flip.make_contact_cq(live=True)` | — | `gfet.quantum_capacitance` |

A guard written to §4.6.17(c)'s wording would have passed the composed flip.

#### (c) What the question did not ask: a write/read relation

Reader C reports that **all five** replacements read `gfet.C_ox` at call time
and that **three** of them write it. `gfet.C_ox` is in nobody's identity
class, so it is invisible to §4.6.13's mechanism, to §4.6.14's, and to the
question of §4.6.17(c), all three of which are about names a flip *rebinds*.

The relation that governs composition is therefore not a set of names but a
**bipartite relation between what replacements write and what replacements
read**, and §4.6.17(d)'s mechanism (ii) is one edge of it rather than a quirk
of one function. Its reachability under the repaired flip is currently zero —
but that zero is a property of the capture, not of the write, and it returns
the moment the capture is undone.

#### (d) Capture is one level deep

A cell protects the **name** it holds. It does not protect the names the held
**object** reads at call time:

| Captured object | its own live reads | in the rebind set? |
|---|---|---|
| `gfet.quantum_capacitance` | `{gfet._thermal_rounding}` | **no** — safe to capture |
| `gfet.carrier_density` | `{gfet._thermal_rounding, gfet.quantum_capacitance}` | **yes** — not safe |

> **Annotated 2026-10-24 (§4.6.20). Every entry in this table is correct and
> the conclusion drawn in its third column is not.** "In the rebind set?" asks
> whether a flip over the name exists *today*; "safe to capture" was read off
> it as if it answered whether the name can move the result. §4.6.20 measures
> the second question directly and finds that a rebind of
> `gfet._thermal_rounding` — the single live read of the captured object, the
> `no` row — moves the repaired contact replacement by **+58.3421 %**, to
> bitwise the same value a rebind of `gfet.quantum_capacitance` moves the
> *unrepaired* one to. Not being in a rebind set is not safety; it is the
> absence of a mechanism that would notice. The numbers and the read sets
> above are unrevised.

Both captures are the same construction and are indistinguishable at the call
site. §4.6.17's repair is safe because of a contingent fact about the one
function it happens to capture. A planted control holding `gfet.carrier_
density` by the identical construction is still flagged, which demonstrates
the limit rather than asserting it.

Stated in the form the other members of this family take:

- §4.6.13 — a repair to one binding is not a repair to the class.
- §4.6.14 — a class can span two operating points, so it must be partitioned.
- §4.6.17 — a replacement that *reads* a rebound name is neither stale nor
  unassigned.
- **§4.6.18 — a capture that fixes one name does not fix the names the
  captured object reads. The safe criterion is the transitive live-read set
  of the captured object.**

#### (e) What this section does not claim

It does not claim the repository is clean. It audits the four replacements
§4.6.17 names plus the two variants of the contact replacement; the eleven
surviving hand-listed rebinding sites are not replacements and are not
covered. It does not flip the default. And reader A is an
over-approximation by construction: it reports the `live=True` branch for
both variants of `make_contact_cq`, which is correct — the read is in the
function — and is why the three readers are reported separately rather than
reduced to one verdict.

Machine output: `replacement_read_audit_output.txt` (30/30 validations,
transcript byte-identical across two runs). Figure:
`replacement_read_audit.png`. Notes:
`notes/2026-10-22-what-a-replacement-reads.md`.

<a id="s4-6-19"></a>

*Moved from Chapter 4, §4.6.19. Original section number kept.*

### 4.6.19 The write/read relation as a guard: the two composition defects are two edges of one relation

Section 4.6.17 reported a replacement reading a name its own flip rebinds.
Section 4.6.18 §4 reported `gfet.C_ox`, written by three of the five
replacements and read by all five, in no identity class and therefore
invisible to `bc.rebind()`, to `bc.rebind_by_operating_point()` and to
4.6.18's own question — all three of which are about names a flip *rebinds*.
It closed by stating that what governs composition is not a set of names but
a bipartite relation between what replacements write and what they read.

**Those are two sources of writes into one relation, not two findings.** A
name is written either by a replacement or by the flip itself, which is every
member of the rebind set. Composing the two makes the composed column —
the one 4.6.18 identified as the column that answers the question —
computable instead of argued. With reads traced and each edge derived from
the read:

| composition | observed `(writer, name, reader)` edges |
|---|---|
| unrepaired (`live_contact=True`) | `THE FLIP --gfet.quantum_capacitance--> make_contact_cq`; `make_contact_cq --gfet.C_ox--> qca._quantum_capacitance_consistent`; `make_contact_cq --gfet.C_ox-->` itself |
| repaired (`live_contact=False`) | `make_contact_cq --gfet.C_ox-->` itself |

The first unrepaired edge **is** 4.6.17's defect; the second **is**
4.6.18 §§4–5. The survivor is the contact substituting `gfet.C_ox` for
`cd.C_interface` and reading it back inside its own call — the same operating
point, restored before return — so it is the substitution and not a leak.
**One change removes both cross-replacement edges**, which is 4.6.18 §5's
"one cause, two consequences" as a measurement rather than as a reading of a
coincidence. No number in 4.6.18 is revised; what is revised is the status of
the relation, from a census to a guard.

**The guard refuses, and it refuses in both directions.** Each edge must
carry a verdict from `{GOVERNED, IMPOSSIBLE, OCCURS}` and a mandatory reason,
in the vocabulary `bc.rebind_by_operating_point()` already uses for
assignments; an unadjudicated edge raises before anything runs. What is new
is that completeness is static but a verdict is a claim about
**reachability**, which no static set contains — so every verdict is checked
against an order probe. `GOVERNED` is refused if the name is in no rebind
set, `IMPOSSIBLE` if the probe observes the edge, and `OCCURS` **if the probe
never observes it**. The last is the half a census cannot have: a guard that
only catches a wrongly-permissive claim answers a weaker question than it was
written to ask, and an `OCCURS` reason describing an edge that has stopped
occurring stays true-looking indefinitely because nothing asks it the
opposite question.

The operational result is that the repaired composition's adjudication — the
same stated reasons, unchanged — is **accepted** for `live_contact=False` and
**refused** for `live_contact=True`, naming three disagreeing verdicts. The
distinction between the two compositions is now a raise rather than a table.

Validation is anchored to exactly-known values: the 4.6.18 arrival counts are
**reproduced** rather than cited (0 arrivals repaired, 192 unrepaired, exactly
181 of those with `C_ox` already substituted), the repaired composition has
**exactly zero** cross-replacement edges, and the deriver itself is put under
synthetic traces whose answers are known without it — a deriver that returns
few edges and a repository that has few are otherwise indistinguishable.

Two defects in the work that produced this section are recorded because both
generalise. The first probe **did not trace reads**: it closed an edge on
every *arrival* following a write, which asks whether a replacement *ran*
after a name moved rather than whether it *read* it, and reported twelve
edges of which four were on names the reader does not read at all. Same shape
as the dominant failure mode, reached through the choice of which events to
instrument. The second is sharper: a **planted control found a crash in a
refusal path that nothing else reaches**, so without it the guard would have
shipped with one of its three refusals unable to report itself while looking
correct, because the two refusals the live subject exercises both work. The
generalisation, and the fourth member of the 4.6.14 / 4.6.17 / 4.6.18 family:
**a refusal that nothing reaches is not a refusal**, so every verdict in a
vocabulary needs its own planted control and not only the ones the subject
happens to exercise.

Machine output: `write_read_guard_output.txt` (21/21). Module:
`graphene_write_read_guard.py`. Notes:
`notes/2026-10-23-two-findings-one-relation.md`.

<a id="s4-6-20"></a>

*Moved from Chapter 4, §4.6.20. Original section number kept.*

### 4.6.20 A capture is adjudicated by a counterfactual, and the repaired contact's exact twin was open all along

§4.6.17 repaired a composition defect by binding the contact replacement to
the shipped `gfet.quantum_capacitance` **captured at manager entry**, so that
a flip over that name could not reach it. §4.6.18 (d) recorded the limit of
that construction — a cell protects the name it holds, not the names the held
object reads — and judged the repair safe on the ground that the held
object's own live reads, `{gfet._thermal_rounding}`, are in no rebind set.
§4.6.19 made the write/read relation a guard. This section asks what the
capture is actually worth, and answers it with a measurement rather than with
a read set.

#### (a) The instrument used to make the claim cannot check it

§4.6.17's argument and §4.6.18's table are both written in the vocabulary of
`rra.static_reads`. On the two objects that the entire repair is the
difference between, that reader returns:

| | static live-read set |
|---|---|
| `flip.make_contact_cq(shipped_cq, live=False)` — repaired | 8 names |
| `flip.make_contact_cq(shipped_cq, live=True)` — unrepaired | the same 8 |

with symmetric difference **zero**, both containing `gfet.quantum_capacitance`
as a *live* read. `cq_of_n` carries both branches in its bytecode and a
bytecode walk cannot know which runs, so the reader is sound and
over-approximate exactly as designed — and therefore **cannot distinguish the
repaired construction from the unrepaired one at all**. The correctness
argument was never checkable where it was made.

#### (b) The verdict is a counterfactual, and it is exact

For each name on the capture surface — the holder's own live reads together
with the live reads of everything it holds — the name is rebound to a
deliberately wrong stand-in scaled by ×1.5, the holder is called, and the two
results are compared **bitwise**. Nothing is compared to a tolerance: *does a
rebind of this name reach this object* has an exact answer. An **identity**
rebind of all nine names moves the result by exactly zero, checked before any
verdict is used.

Contact replacement at `V = 0.5 V`, base `C_q = 6.555495551420368 × 10⁻²`
F/m² (`0x1.0c835a6063492p-4`):

| Name on the capture surface | verdict | repaired | unrepaired |
|---|---|---|---|
| `gfet._thermal_rounding` | EXPOSED | **+58.3421 %** | +58.3421 % |
| `qca.A_FERMI` | EXPOSED | +58.3417 % | +58.3417 % |
| `qca.series_factor` | EXPOSED | +26.1720 % | +26.1720 % |
| `qca._fixed_point` | EXPOSED | +22.4740 % | +22.4740 % |
| `cd.C_interface` | EXPOSED | +18.1154 % | +18.1154 % |
| `gfet.T` | EXPOSED | +0.0034 % | +0.0034 % |
| `cd.quantum_capacitance` | UNREACHED | exact 0 | exact 0 |
| `gfet.C_ox` | PROTECTED | **exact 0** | exact 0 |
| `gfet.quantum_capacitance` | PROTECTED | **exact 0** | **+58.3421 %** |

The two replacements differ on **exactly one** of the nine names, and it is
the one §4.6.17 is about. That column is the repair, measured.

#### (c) The result: the repair has an exact twin it does not cover

A rebind of `gfet._thermal_rounding` moves the **repaired** replacement to
`0x1.a92b527d6852cp-4` — **bitwise** the value a rebind of
`gfet.quantum_capacitance` moves the **unrepaired** one to.
`gfet._thermal_rounding` is the *sole* live read of the captured object, so
the route is the capture itself. §4.6.17 closed one path to that number and
left an exactly equally powerful one open **through the object it captured**.

It is in no flip's rebind set and is written by no replacement, so
`bc.rebind_by_operating_point`, §4.6.19's write/read guard and §4.6.17's
argument all reach it equally: not at all. **A name outside every rebind set
is not thereby harmless — it is only unwatched**, and §4.6.18 (d) is
annotated in place to say so.

#### (d) A second result, in the opposite direction

Walking the holder's closure cells on **object identity** recovers not one
held name but the held object's whole identity class,
`{gfet.quantum_capacitance, cd.quantum_capacitance}`, reproducing
`bc.identity_class(gfet, 'quantum_capacitance')` exactly and independently —
two derivations of one set, which is this section's validation against the
existing implementation in the overlapping limit. The capture therefore
neutralises a rebind of **both** members, including
`blast.contact_on_consistent_argument()`'s, which is **stronger** than
§4.6.17 claimed. But only one member is ever looked up, so calling both
PROTECTED would credit the capture with work it did not do — §4.6.18 (d)'s
conflation with the sign reversed. `UNREACHED` keeps the two apart, and
separating them needs the one judgement `static_reads` **can** make soundly:
an over-approximate read set has reliable *negatives*, so *not on the read
surface* is trustworthy even though *on it* is not.

> **Annotated 2026-10-25 (§4.6.21). The sentence this paragraph ends on is
> false.** `rra.static_reads` records a `(module, attr)` pair only when a
> module object is loaded and an attribute taken from it, or when a bare
> global resolves to a repo *function*. A bare global resolving to a
> module-level **scalar**, read from inside the module that defines it —
> which is how every constant here is read by its own module's functions —
> is recorded **nowhere**, so a negative means *not read in a spelling the
> walker can see*, not *not read*. The demonstration is a single name:
> `qca.A_FERMI` is **off** the read surface of both channel replacements
> while a rebind of it moves one of them by **+29.8900 %**, and it is **on**
> the surface of both contact holders, where `blast` spells it qualified.
> §4.6.21 measures this. **Nothing in (d) above is withdrawn**: the identity
> class, the two independent derivations of it and the `PROTECTED`/`UNREACHED`
> distinction all stand, and the `UNREACHED` verdicts of §4.6.20 are
> re-measured in §4.6.21 and still hold. What is withdrawn is the *ground*
> offered for them here. `UNREACHED` is sound when it rests on a measured
> bitwise zero, and unsound when it rests on an absence.

#### (e) A device-physics consequence, not only a software one

`gfet.C_ox` is **exactly** protected for both variants. §4.6.19 left one
`OCCURS` edge standing — the contact replacement's write-then-read of
`gfet.C_ox` inside its own call — adjudicated *the substitution itself and
not a leak* on the argument that the write and the read are at the same
operating point. The counterfactual converts that argument into a measured
zero: the **external** value of `gfet.C_ox` reaches the contact
replacement's output not at all, because the replacement substitutes
`cd.C_interface` for the whole solve and restores the old value before
returning. **The contact self-consistency is referenced to `C_interface` and
is exactly independent of the oxide capacitance**, so no §4.6 contact number
inherits an oxide-thickness error through that path. The three names that
*do* carry the contact `C_q` are `gfet._thermal_rounding` and `qca.A_FERMI`
at ≈ +58 %, and `cd.C_interface` itself at +18 %, which is the sensitivity
ranking a reader of §4.5 should have.

> **Annotated 2026-10-25 (§4.6.21). The `gfet.C_ox` result stands; the
> sensitivity ranking in the last sentence is incomplete and its "three
> names" is wrong.** The exact zero on `gfet.C_ox` is re-measured in
> §4.6.21 on `blast._cq_consistent_contact`, which has **no capture at
> all**, so the zero is a property of the contact model and not of §4.6.17's
> repair — a strengthening of this paragraph, not a retraction. The ranking
> is another matter. Five further names reach the contact `C_q` and were
> invisible to the surface this paragraph was computed from:
>
> | name | leverage on contact `C_q` |
> |---|---|
> | `gfet.e` | **+278.4911 %** |
> | `gfet._thermal_rounding` | +58.3421 % |
> | `qca.A_FERMI` | +58.3417 % |
> | `gfet.hbar` | −63.5273 % |
> | `gfet.v_F` | −63.5273 % |
> | `qca.E_CHARGE` | −21.2518 % |
> | `cd.C_interface` | +18.1154 % |
>
> The true ranking is headed by `gfet.e`, and `gfet.hbar` and `gfet.v_F`
> are interleaved among the three names named above. The numbers quoted in
> this paragraph are all correct; the claim that they are *the* ranking is
> not, and §4.6.21 (c) draws the device-physics consequence for `v_F`.

#### (f) It is a guard, not a census

`graphene_capture_surface_guard.py` raises `UnadjudicatedCaptureSurface`
before measuring anything if a surface name carries no `(verdict, reason)`,
and `FalseCaptureVerdict` after measuring if a verdict disagrees with the
counterfactual **in either direction**: PROTECTED on a name that moves,
EXPOSED on a name that does not, UNREACHED on a name that *is* read,
PROTECTED on a name that is read by nothing, GOVERNED on a name in no rebind
set. One refusal is §4.6.18 (d)'s error verbatim — a PROTECTED verdict on a
held name obliges every live read of the object held there to carry its own
verdict, and any that moves must be adjudicated EXPOSED. The repaired
capture's adjudication is **refused** for the unrepaired replacement, which
(a) showed the static reader cannot tell apart.

29/29 checks, with a planted control for every verdict and every refusal
branch — the §4.6.19 lesson that a refusal nothing reaches is not a refusal,
applied inside the module that states it. Seven defects injected one at a
time into a copy of the repository are all detected (7/7); undoing §4.6.17's
repair drops the score to 23/29, and replacing the bitwise comparison with a
1 × 10⁻³ relative tolerance is caught too, so the result depends on exactness
rather than on a threshold.

Stated in the form the other members of this family take:

- §4.6.13 — a repair to one binding is not a repair to the class.
- §4.6.14 — a class can span two operating points, so it must be partitioned.
- §4.6.17 — a replacement that *reads* a rebound name is neither stale nor
  unassigned.
- §4.6.18 — a capture that fixes one name does not fix the names the captured
  object reads.
- §4.6.19 — a refusal that nothing reaches is not a refusal.
- §4.6.20 — a name outside every rebind set is not thereby harmless, only
  unwatched. Whether a name can move a result is a counterfactual, and no
  read set answers it.
- **§4.6.21 — a name's absence from a read surface is a fact about the
  surface's vocabulary, not about the name. Two names for one parameter are
  not one knob, and a precomputed copy is unreachable from the name it was
  computed from.**


**Annotated 2026-10-26 (§4.6.22).** This section's `gfet.C_ox` **PROTECTED**
verdict on `blast._cq_consistent_contact`, and the number behind it, are
**correct and are not withdrawn** — but every leverage figure in this section
and in §4.6.21 was measured with the replacement's **flip not entered**, and
§4.6.22 shows that for this one holder the two configurations disagree about
eleven of fifteen names. The stated reason here generalised from one
configuration; §4.6.22 §(c) states the scope and corrects the wording in the
code, keeping the superseded text. The §4.6.20 counterfactual instrument
itself is vindicated: it is bitwise identical to the as-installed instrument
for five of the six objects the flips install.

<a id="s4-6-21"></a>

*Moved from Chapter 4, §4.6.21. Original section number kept.*

### 4.6.21 The negatives were never sound, and `v_F` was never one name

§4.6.20 adjudicated one capture and closed with a sentence it did not test:
that an over-approximate read set has reliable *negatives*. This section
tests it on the three remaining replacements, finds it false, and follows the
consequence into a statement about the Fermi velocity that is a device-physics
result rather than a software one.

#### (a) The blind spot, and one name that exhibits it

`rra.static_reads` records a name in two situations only: a module object
loaded and an attribute taken from it (`qca.A_FERMI`), or a bare global
resolving to a repo **function** (`cq_thermal`). A bare global resolving to a
module-level **scalar** — the spelling every constant in this repository is
read by inside its own module — matches neither, and is recorded nowhere.

The cleanest exhibit is one name read two ways:

| holder | on the §4.6.20 surface? | leverage of `qca.A_FERMI` |
|---|---|---|
| `qca._carrier_density_consistent` | **no** | +0.0121 % |
| `qca._quantum_capacitance_consistent` | **no** | **+29.8900 %** |
| `blast._cq_consistent_contact` | yes | +58.3417 % |
| `flip.make_contact_cq(live=False)` | yes | +58.3417 % |

One name, one value, one mechanism; `qca` reads it bare and `blast` spells it
`qca.A_FERMI`, and that spelling is the whole of the difference. This is an
exact binary prediction of the mechanism on a single name, and it is checked
as such.

#### (b) The answer to §4.6.20's own question, in both directions

§4.6.20 asked whether `gfet._thermal_rounding` and `qca.A_FERMI` carry its
≈ +58 % onto the other replacements, which would make them repository-wide
levers rather than contact facts. The two halves answer oppositely.
`gfet._thermal_rounding` is **bitwise exactly zero** on both channel
replacements — a contact fact, and the zero is a *measured counterfactual*,
not an absence from a surface, which after (a) would prove nothing.
`qca.A_FERMI` is repository-wide, and off-surface in both places where it
matters.

Completing the surface on the four holders adds **18 (holder, name) pairs**,
8 distinct names, and **every one moves the value**; the largest is `gfet.e`
at +278.4911 % and the largest derived model prefactor is `qca.CQ_PREFACTOR`
at +50.4145 %. None had ever carried a verdict, because a guard can only
demand one for a name its surface names — and the §4.6.20 guard is shown to
**accept** an adjudication of all three replacements with every one of them
omitted. Repository-wide, **75 of the 90 module-level scalars (84 %)** are in
the blind class; `gfet.t_ox` and `gfet.epsilon_0` are scalars outside it,
because they are read only at import time to build `C_ox`, which pins the
census from the other side.

#### (c) The device-physics consequence: `v_F` is two names, not one

`C_q` carries the Fermi velocity twice, as `CQ_PREFACTOR = 2e³/(π ħ² v_F²)`
and as `A_FERMI = ħ v_F √π / e`, so a consistent change of `v_F` by *k* must
move one as *k*⁻² and the other as *k*⁺¹. The two families do not share a
name for it:

| holder | `gfet.v_F` | `qca.CQ_PREFACTOR` | `qca.A_FERMI` |
|---|---|---|---|
| `qca._carrier_density_consistent` | **exact 0** | +0.000671 | +0.000247 |
| `qca._quantum_capacitance_consistent` | **exact 0** | +1.008232 | +0.519524 |
| `blast._cq_consistent_contact` | −2.323463 | **exact 0** | +1.161680 |

*(d ln value / d ln name; `exact 0` = a ×1.5 rebind leaves the value bitwise
unchanged.)* The contact replacement reads `gfet`'s **live** prefactor;
both channel replacements read `qca.CQ_PREFACTOR` and `qca.A_FERMI`, which
`graphene_exact_charge_rederivation` computed from `gfet.v_F` **at import
time** and copied in. So `gfet.v_F` is bitwise unreachable from every
Chapter 4 channel number, and `qca.CQ_PREFACTOR` is bitwise unreachable from
the contact number. **A `v_F` sensitivity study performed by editing
`gfet.v_F` would have reported channel insensitivity to the Fermi velocity
as a result.** This is the §4.6.13 delivery fault in a new place: the name is
neither stale nor unassigned but **precomputed**, and `bc.identity_class`
cannot see it because the two values are not the same object — they are a
float and a float derived from it.

#### (d) And the perturbation is the physical case

`v_F` is not a constant of nature. Hwang *et al.* measure, by ARPES,
1.15 ± 0.02, 1.49 ± 0.08 and 2.49 ± 0.30 (×10⁶ m s⁻¹) for graphene on
SiC(000-1) (ε = 7.26), BN (ε = 4.22) and quartz (ε = 1.80) respectively,
against a fully screened *ab initio* value of 0.85 × 10⁶ m s⁻¹, with a 1/ε
dependence at fixed carrier concentration [Sci. Rep. **2**, 590 (2012),
doi:10.1038/srep00590].

The value committed throughout this chapter, `v_F = 1.00 × 10⁶ m s⁻¹`, sits
just above the **bare** value. The 90 nm SiO₂ back gate modelled in §4.2 has
ε ≈ 3.9, between the BN and quartz cases, so the substrate-appropriate value
is near 1.5 × 10⁶ m s⁻¹ — which is exactly the ×1.5 the probe applies by
default. The **−63.5273 %** it reports for the contact `C_q` is therefore not
a synthetic perturbation: it is the SiO₂-screened Fermi velocity, and this
chapter quotes the contact `C_q` at a `v_F` its own substrate does not have.

This is left as an open question rather than resolved here, and §7.9 records
it as such. Using the bare value consistently is defensible — it is what most
compact device models do, and the renormalisation is itself carrier-density
dependent, so any single number is already an approximation. Quoting it as
*the* contact `C_q` while the same chapter models a specific SiO₂ stack is
weaker. What is not defensible is the present state, in which the question
cannot be *asked* by changing one name.

#### (e) Anchored to an exactly known symmetry

All of (c) and (d) rest on a probe, so the probe is anchored to a value known
in advance: `gfet.hbar` and `gfet.v_F` must have **bitwise identical**
leverage, because the prefactor depends on them only through `(ħ v_F)²`. They
do — both move the contact `C_q` to `0x1.87bc5bc0c26d0p-6`. A probe reaching
the value by any other route would not reproduce that coincidence. 36/36
checks, with 9 defects injected one at a time into a copy of the repository
and all 9 detected; the mutation that hides this section's two findings behind
the measuring module's own apparatus carve-out scored **31/31** against the
first draft of the check meant to prevent exactly that, which is recorded in
`notes/2026-10-25-the-negatives-were-never-sound.md` §8 rather than quietly
fixed.

**Annotated 2026-10-26 (§4.6.22).** Every number in (a)–(d) above is measured
on the replacement **in isolation**, with its own flip not entered, and
§4.6.22 shows that this is bitwise equivalent to the as-installed measurement
for both channel holders and for `flip.make_contact_cq(live=False)` — so
(a)–(c) stand unrevised, and the ×1.5 `v_F` figure of **−63.5273 %** holds in
*every* configuration on the holder the default flip actually installs. The
one exception is `blast._cq_consistent_contact`, on which `gfet.v_F` is
**bitwise dead** once either other flip is entered. The 18 added pairs and
the 84 % blind-class census are unaffected: both are statements about
spelling, not about configuration.

<a id="s4-6-22"></a>

*Moved from Chapter 4, §4.6.22. Original section number kept.*

### 4.6.22 A verdict measured in isolation is a verdict about a configuration the replacement never runs in

§4.6.20 and §4.6.21 adjudicate four *holders* — the replacement functions the
flip managers of §4.6.16 install. Both measure them by calling them
**directly**, with no manager entered. A replacement only ever runs with its
manager entered, and all three instruments used here are relative to the
namespace at the moment of use: `csg.capture_surface` and `vrs.value_reads`
are static walks that resolve each name to the object *bound now*, and
`csg.leverage` calls the holder in the current namespace. So the question
§4.6.21 left untouched is whether isolation and installation agree. Until this
section, nothing had asked.

#### (a) The population, derived rather than listed

`rra.FLIPS` × `rra.rebind_set` gives **8 (flip, site) pairs and 6 distinct
installed objects**. §4.6.20 and §4.6.21 between them audited four. The two
that had never had a read surface computed at all are the arrival-counting
wrappers `flip._counting` builds — the objects `gfet.carrier_density`,
`gfet.quantum_capacitance` and `crx.carrier_density` *actually become* under
the default flip.

#### (b) The answer: the instruments are vindicated on five of six objects

Measured over all 6 × 3 cells, with the dependence of each (holder, name)
pair classified on a four-valued axis (`BITWISE_SAME`, `BOTH_LIVE`,
`DEAD_AS_INSTALLED`, `LIVE_ONLY_AS_INSTALLED`):

| installed object | under its own flip | under the other two |
|---|---|---|
| `qca._carrier_density_consistent` | 14/14 `BITWISE_SAME` | invariant |
| `qca._quantum_capacitance_consistent` | 11/11 | invariant |
| `_counting(_carrier_density_consistent)` | 15/15 | invariant |
| `_counting(_quantum_capacitance_consistent)` | 12/12 | invariant |
| `flip.make_contact_cq(live=False)` | 19/19 | invariant |
| **`blast._cq_consistent_contact`** | **15/15** | **not invariant under either** |

**Every one of the six is bitwise flip-invariant under its own flip** — every
name on its surface, every counterfactual, and the holder's own value. That
proposition is what licenses every leverage figure quoted in §4.6.20,
§4.6.21 and §4.6.18, and it had never been stated, let alone checked. It is
now an assertion rather than an assumption.

A **surface** change is also not a **verdict** change, and the measurement
separates them: four of the six objects gain names on their surface in some
configuration — the counting wrapper for `carrier_density` goes 15 → 17 under
`qca.on_consistent_argument()` — while every one of their counterfactuals
stays bitwise identical.

#### (c) The one row that moves, and what it does and does not revise

`blast._cq_consistent_contact` under `flip.on_default_flip(live_contact=False)`
returns `0x1.4ef081d146707p-7` instead of `0x1.0c835a6063492p-4`, a factor of
**6.4134**, and **eleven of its fifteen names cross the exact-zero
boundary**: five go bitwise dead (`gfet.e` +278.4911 %, `gfet.hbar` and
`gfet.v_F` −63.5273 %, `gfet.kB` +0.0242 %, and `gfet._thermal_rounding`
+58.3421 % — the whole of §4.6.20's finding) and six come alive from bitwise
zero (`qca.cq_thermal`, `qca.cq_thermal_of_n` and `qca.CQ_PREFACTOR` at
+52.4103 %, `qca.KB` +35.8863 %, `qca.n_arg_EF_thermal` +9.0738 %).

**The value change is not a new result.** It is §4.6.17's — *the two halves do
not compose* — which measured this exact composition and costed it at −51 %
to −62 % in `n_contact`. What is new is that the same non-composition appears
as a change in the **read surface**, which is the object §4.6.19's and
§4.6.20's guards demand adjudication over. §7.9's standing item on the
unadjudicated compositions is therefore partly discharged: the relation is
measured over all 6 × 3 cells and exactly one row is non-composing.

What it revises is one **reason**, not one number. §4.6.21's dictionary gives
`gfet.C_ox` on this holder the verdict PROTECTED with the reason *"the contact
solve closes against `C_interface`, so the oxide capacitance reaches this
value not at all … the zero is a property of the contact model and not of the
2026-10-21 repair."* The first clause is right. The generalisation is not: the
same holder gives `gfet.C_ox` **+7.8110 %** under the default flip, so the
exact zero is a property of the contact model **composed with nothing else**.
The reason is corrected in place in the code and the superseded wording is
kept there.

#### (d) And the vocabulary cannot express the state of that name

Exhaustively over the four verdicts, for `gfet.C_ox` on this holder:

| verdict | in isolation | as installed |
|---|---|---|
| `PROTECTED` | admissible | refused — it moves |
| `UNREACHED` | refused — it is read | refused — it moves |
| `EXPOSED` | refused — bitwise zero | admissible |
| `GOVERNED` | refused — in no rebind set | refused |

Exactly one verdict is admissible in each configuration **and they are
different ones**. The four-verdict vocabulary of §4.6.20 is therefore not
merely incomplete here; non-composition is **not dischargeable by
adjudication at all**, which is the structural reason that item has survived
three sessions of being noted. A fifth value — a verdict indexed by
configuration — is the repair, and it is not attempted here.

#### (e) A second finding: `flip._counting` is a capture

Each counting wrapper's read surface contains exactly one name the
replacement it wraps does not — that replacement's **own name** — and a
rebind of it leaves the wrapper's value **bitwise unchanged**, because
`_counting` holds the function object in a closure cell. This is a **third**
capture in the model code, after the §4.6.17 contact repair and the
`bc.rebind` family, and nothing had noticed it; its consequence is that the
channel replacement *name* is immune under the default flip. The verdict is
`UNREACHED` and not `PROTECTED`, by §4.6.20's own argument for the other
identity-class member: the wrapper never *looks the name up*, so the bitwise
zero is not a capture neutralising a read. This section's first draft said
`PROTECTED` and the §4.6.20 guard refused it — the fourth time a rule planted
in an earlier section has caught a claim written in a later one.

#### (f) Exact anchors, and two of this section's own checks

The result rests on probes, so the probes are anchored to values known in
advance. `gfet.hbar` and `gfet.v_F` move `flip.make_contact_cq(live=False)` to
the **bitwise identical** `0x1.87bc5bc0c26d0p-6` in all four configurations,
because the `C_q` prefactor depends on them only through `(ħ v_F)²`; the
counting wrapper returns a value bitwise equal to the replacement it wraps;
an identity control flip leaves value and surface bitwise unchanged; and —
the sharpest of the four — **the repaired and unrepaired contact replacements
are bitwise indistinguishable in isolation**, both `0x1.0c835a6063492p-4`, so
the entire content of the §4.6.17 repair appears *only* as installed. Five
sections measured both holders in the one configuration in which they cannot
be told apart.

Nine defects were injected one at a time into a copy of the repository, each
required to be caught **by a named check** rather than by a non-zero exit, and
**two of this section's own checks failed that bar on the first run**: the
positive control for the new PROTECTED-as-installed refusal planted its false
verdict on a holder whose adjudication was incomplete in that configuration,
so the missing-name refusal fired first and the branch under test was never
reached (26/26, undetected); and the flip-invariance check asserted two hex
values but never the invariance *flag* computed from them, so a relative
tolerance could report a factor of 6.4 as invariant (27/28, undetected). Both
are repaired. Final: **28/28 checks, 9/9 defects detected on their own terms,
control clean.**

#### (g) What this section does not claim

It measures three configurations, not all of them: `live_contact=True` and the
*simultaneous* composition of two managers are not covered. "Invariant under
its own flip" is a bitwise statement at one probe point with one perturbation,
not a theorem. And `vrs.value_reads` remains over-approximate with unsound
negatives (§4.6.21 (a)); nothing here repairs that.

Machine output: `installed_surface_audit_output.txt`,
`installed_surface_mutation_output.txt`. Working notes:
`notes/2026-10-26-the-configuration-the-replacement-runs-in.md`.

<a id="s4-6-23"></a>

*Moved from Chapter 4, §4.6.23. Original section number kept.*

### 4.6.23 A repaired instrument's own repair was mis-predicted twice, and the second prediction was the informative one

Every adjudication in §4.6.16 through §4.6.22 is taken with two constants:
`csg.PROBE_V`, the gate voltage the holder is evaluated at, and
`csg.PERTURBATION`, the factor the counterfactual multiplies by. Until
2026-10-27 both were captured as function defaults rather than read at call
time, which put them in the delivery class `MIXED` — some readers live, some
frozen — and `MIXED` is not an intermediate state. A rebind of a `MIXED` name
reaches the live readers and not the frozen ones, so if the live readers are
an audit's expected-value expressions and the frozen ones are its measuring
functions, the audit reports a disagreement its own mutation manufactured.
These two names are not incidental to this chapter: **every counterfactual in
it is taken with them**, so the split was a defect in the instrument and not
in any model the instrument measures. Three further names were `MIXED` for the same
reason in `graphene_chord_decay_length.py`, the module behind §7.9 item 2's
2026-10-20 annotation — the one that shows `γ` cancels from the chord/local
decay-length gap and that `ℓ` is the local decay length at no metal in the
table. Eleven of the repository's
twelve capture sites are converted; the twelfth is discussed below.

**No number in this chapter moves.** The two affected modules' machine
transcripts are bitwise unchanged by the conversion, which is the standard the
2026-09-30 and 2026-10-01 conversions were held to and the only available
evidence that a repair to a measuring instrument did not also move the
measurement. What changes is what a *future* mutation of those names will do.

**The informative part is that the repair was mis-predicted, twice.** The
suite written to measure it predicted the textbook `MIXED` signature on
`D1_CHORD`: the live reader `L'(d)` moves under a rebind and the captured
reader `⟨L'⟩` freezes. It does not freeze. `⟨L'⟩` captured the chord endpoint
`d₁` as a default, but the function it calls to evaluate `P` reads `D1_CHORD`
*live*, to re-centre the polynomial at `d − d₁`. So under a rebind the
pre-repair code averaged `L'` over the **old** interval while evaluating a
polynomial centred on the **new** `d₁`. The resulting number corresponds to no
consistent choice of `d₁` at all. That is strictly worse than the frozen case
the prediction assumed, because a frozen number is at least one the model once
produced, and the quantity carrying it was Eq. (D) — the γ-free discrepancy
this chapter quotes.

Restating the check as *coherence* rather than *movement* then caught a second
wrong draft of the same kind. The reference value for "the old interval" was
first taken **with the rebind still in place**, which made the reference the
same incoherent quantity it was meant to be distinguished from; it agreed with
the measurement to `0.000000e+00` and the check failed for the right reason.
A reference value has to be taken in the state it claims to describe. Under a
rebind of `d₁` from 3.3 Å to 3.9 Å the pre-repair default path lands
`1.617099e-02` from the new-interval value and `3.838368e-02` from the
pre-rebind one — between the two and equal to neither. After the repair it
equals the new-interval value **bitwise**.

**The residual, and a conflict this chapter does not resolve.** One capture
site is left: `n_quad=N_QUAD` in the velocity-saturation module. Its `def` sits
inside a function *body*, so the default is re-evaluated on every call of the
enclosing function and a rebind does arrive — it is live by scope, not frozen,
and is measured as live. It is also the only *real* site the
frozen-default-semantics suite has for demonstrating that a static pass can
find a capture its confirmation step then silently discards. So the standing
guard "no module-level constant is captured as a default anywhere" cannot go
green while that demonstration keeps its site, while its sibling "and in
particular none is `MIXED`" now can, for the first time since it was written
on 2026-10-01. The conflict is in the guard's **vocabulary** — it counts
syntactic captures and names the class it cares about only in the second check
— and not in the tree. Resolving it means either splitting the guard, or
trading a real demonstration site for a synthetic one, which is the trade that
suite exists to argue against.

Machine output: `mixed_capture_repair_output.txt` (33/33),
`mutation_arrival_probe_output.txt` (Section 7b, 32/2 → 33/1). Notes:
`notes/2026-10-27-the-repair-that-refuted-its-own-prediction.md`.

**Resolved 2026-10-28** by splitting the guard — see §4.6.24. The paragraph
above is left as written because it is the diagnosis the split was built from.

<a id="s4-6-24"></a>

*Moved from Chapter 4, §4.6.24. Original section number kept.*

### 4.6.24 A guard stated so that it cannot be true, and a guard blind to its own source

The conflict §4.6.23 ends on is not about this chapter's physics. It is about
a sentence. Section 7b asserted *"no module-level constant is captured as a
default anywhere"* — a count of **syntactic** sites — and separately *"none is
`MIXED`"* — a statement about a **delivery class**. By 2026-10-27 the second
was true and the first could not become true, because the one remaining
syntactic site is live by scope and is the only real demonstration site the
frozen-default-semantics suite has. A check that cannot go green is not a
guard; and its report, `1 site-groups remaining`, cannot distinguish the known
site from a new one appearing while the known one was removed.

Split into three, with no coverage given up:

| | statement | exemptions | state |
|---|---|---|---|
| 7b-i | no capture site anywhere **freezes a read** (`FROZEN`, `MIXED`, or an unimportable module) | **none, by construction** | GREEN |
| 7b-ii | no capture site in the tree is **unregistered** | one, with its class and its reason | GREEN |
| 7b-iii | every registry entry still **describes the tree** (`VANISHED` / `DRIFTED`) | — | GREEN |

The load-bearing control of the eight planted ones is the seventh: the
*registered* site, planted as `FROZEN`, must still fire 7b-i. A registry that
could silence the semantic check would make the split a weakening rather than
a restatement, so `census_semantic_faults()` does not take the registry as a
parameter at all and no caller can supply one. Suite 33/1 → 43/0.

**Then the split committed the fault it guards against.** The three new
predicates were first written with their tables as parameter defaults —
`faults=SEMANTIC_CAPTURE_FAULTS`, `registry=SYNTACTIC_CAPTURE_EXEMPTIONS` —
module-level constants captured at `def` time, in the module that defines the
class. `graphene_figure_provenance_audit.py`'s repository-wide census of the
pattern went from 1 occurrence in 1 file to 4 in 2 and that is how it was
found. Converted to the 2026-09-30 late-bound form; transcript bitwise
unchanged.

Section 7b did **not** catch it and cannot: `repo_wide_frozen_census()` skips
`graphene_mutation_arrival_probe` by name, because that module is full of
capture patterns planted deliberately for its own controls. The exclusion is
defensible; its consequence had never been written down. **A guard excluded
from its own census is blind in exactly the place where being wrong is most
likely**, since the module that defines a fault class is the one most often
editing that class's shape. This is the same lesson as §4.6.8 — knowing a
failure mode is not the same as not committing it — and a sharper instance of
it, because here the committing module is the defining module.

For this chapter's numbers nothing changes: the residual is still one site and
still the same site, and `mixed_capture_repair_output.txt` reports 33/33
before and after. What changes is that the guard now says which site, in which
vocabulary, and why it is permitted.

Machine output: `mutation_arrival_probe_output.txt` (Section 7b split, 43/0),
`figure_provenance_audit_output.txt` (Section C, 4 occurrences → 1 after the
conversion), `pristine_transcript_audit_output.txt` (19/0). Notes:
`notes/2026-10-28-three-unchecked-claims-and-the-guard-blind-to-itself.md`.

<a id="s4-7-superseded"></a>

*Superseded 2026-11-16 by the crossover correction ([A.1](#corrections)). Chapter 4 §4.7 as it stood from 2026-08-31, kept as written, and the row of §4.9's summary table that reported it. Every number here was computed with the contact doping measured from graphene's own 4.5 eV rather than from the 5.4 eV n/p crossover. The two annotations that follow annotate this table.*

> ## 4.7 Per-metal Rc recalibration: a genuine negative-residual result
>
> Section 4.5's contact-doping model computes a work-function-dependent
> "extra" resistance, R_extra, contributed by the extended in-plane doping
> junction beyond a contact's edge. Section 4.5 used R_extra only
> *diagnostically* (comparing metals to each other and to the model's own
> qualitative predictions). This section attempts to go one step further —
> recalibrate the model against real, individually-cited, per-metal
> measured Rc values from the literature — and reports a genuine negative
> result rather than a forced fit.
>
> Four literature Rc values were sourced this session (see
> `notes/2026-08-31-per-metal-contact-resistance-literature-and-
> recalibration.md` for full citations and search-access notes; Ti and Cr
> could not be recalibrated because no single clean per-metal figure was
> obtained this session — two ResearchGate fetches were rate-limited),
> all for top (surface) contacts at room temperature, gate-biased
> (on-state): Cu 184 Ω·µm, Ni 470 Ω·µm, Au 519 Ω·µm, Pd 584 Ω·µm.
>
> The natural hypothesis going in was an additive decomposition,
> Rc,measured = R_extra (Section 4.5's doping-junction term) +
> R_transmission (an unmodeled interface-tunneling/mode-limited-injection
> term, concentrated at the contact itself). Computing R_extra at the same
> on-state bulk channel density used throughout this chapter
> (`n_bulk = 2.0e16 m^-2`) and subtracting it from each metal's measured
> Rc (`graphene_contact_doping_model.recalibrate_metal_rc()`,
> `rc_recalibration.png`):
>
> | Metal | Rc, measured (Ω·µm) | R_extra, computed (Ω·µm) | Naive R_transmission (Ω·µm) |
> |---|---|---|---|
> | Cu | 184.0 | 971.5 | **−787.5** |
> | Ni | 470.0 | 830.5 | **−360.5** |
> | Au | 519.0 | 609.2 | **−90.2** |
> | Pd | 584.0 | 533.1 | +51.0 (~9% of total) |
>
> *Two annotations to this table (2026-09-22, 2026-10-13) are in [Appendix A, §4.7 notes](appendix-A-verification-record.md#s4-7-note-a): the negative rows are the robust ones, and under the C_q-argument correction Pd's +51.0 becomes −203.276.*
>
>
> For three of the four metals, the computed doping-junction contribution
> **alone exceeds the entire measured lumped contact resistance**, which
> would require a negative "transmission resistance" to balance the
> books — unphysical, since a resistance cannot be negative. Only Pd gives
> a small, plausible positive residual. Rather than adjusting
> `lambda_decay` or the bulk-density reference ad hoc until the four data
> points come out positive — which would just be curve-fitting a
> parameter to this session's four measurements — this is reported as-is:
> **the simple additive decomposition does not hold**, most plausibly
> because the transfer length method (TLM) used to extract all four
> literature Rc values fits total device resistance vs. contact spacing
> back to zero spacing, and its "contact resistance" term can already
> partially absorb the same near-contact doping-gradient region this
> model integrates separately — i.e. R_extra and Rc,measured likely
> double-count part of the same physical region rather than summing as
> independent series terms. A second, non-exclusive possibility is that
> `lambda_decay` = 250 nm (Section 4.5, from Khomyakov et al. 2010's
> doped-graphene asymptotic regime) over-integrates the extra-resistance
> contribution for these specific literature devices' actual channel
> geometry. Distinguishing between these needs either the source papers'
> extracted transfer length L_T (not available from the excerpts obtained
> this session) or an explicit simulation of the TLM extraction procedure
> on top of the doping profile — both left as open items rather than
> resolved here.
>
> *Section 4.9 summary row, as written:*
>
> | Using the doping-profile model to *recalibrate* Rc per metal (vs. using it only diagnostically) | Attempted (Section 4.7) — additive decomposition found NOT to hold for 3 of 4 metals (Cu, Ni, Au); only Pd gives a physically plausible residual. Root cause (TLM double-counting vs. `lambda_decay` mismatch) not yet isolated; Ti and Cr not recalibrated (no literature Rc sourced this session) |

<a id="s4-7-note-a"></a>

*Annotation moved from the body of Chapter 4, §4.7.*

> **Annotation added 2026-10-13 (Section 4.6.12).** No number in the table
> above is changed or retracted, and the four `R_extra` values remain the
> committed ones. They are computed through
> `contact_edge_carrier_density()`, which evaluates `C_q` at the
> work-function potential rather than at `E_F/e` — the Section 7.9 item 16
> defect, in a second module, and one that no measurement of item 16 had
> reached because the harness rebinds `gfet.quantum_capacitance` while this
> module holds its own import-time reference. Corrected, `R_extra` moves
> −1.2583 % (Cu), +20.8352 % (Ni), +39.2281 % (Au), **+47.6925 % (Pd)**, and
> **Pd's +51.0 becomes −203.276**: all four residuals are then negative. The
> sentence below, "Only Pd gives a small, plausible positive residual", is
> therefore **true of the committed model and false of the corrected one**,
> and is left standing with this annotation rather than rewritten, pending
> the item 16 flip. The ordering of the four metals and the n/p type
> assignments are unchanged. The section's *conclusion* — that the additive
> decomposition does not hold — is strengthened, not weakened: a
> four-of-four failure is more consistent with the TLM double-count
> explanation given below than a three-of-four one was.

<a id="s4-7-note-b"></a>

*Annotation moved from the body of Chapter 4, §4.7.*

> **Annotation added 2026-09-22 (Section 4.10).** No number in the table
> above is changed or retracted. A condition-number audit
> (`graphene_sensitivity_audit.py`) finds that the four rows' **sensitivity
> to input error runs in the opposite order to their apparent
> plausibility**: Pd's `+51.0` — the one row read below as "a small,
> plausible positive residual", i.e. the only evidence that the additive
> decomposition might survive — flips sign on a **9.6 %** error in
> `R_extra`, while Cu's `−787.5` would need **81.1 %**. The negative results
> are therefore the *robust* rows of this table and Pd is the fragile one.
> Section 4.10 gives the full margins and their consequence: the additive
> decomposition fails robustly, and amplified input error is eliminated as
> the cause of the negatives.
>
> *Amended 2026-10-13 (Section 4.6.12). The 9.6 % figure stands as computed
> and is not revised. What this annotation did not know is that it is a
> margin against INPUT error only, and that a known MODELLING error in the
> same quantity already exceeds it by 5x: correcting `C_q`'s argument in
> `contact_edge_carrier_density()` moves Pd's `R_extra` by +47.6925 %. The
> sentence "Pd is the fragile one" was right; the fragility was realised by
> the model and not by its inputs.*

<a id="s4-8-2"></a>

*Moved from Chapter 4, §4.8.2. Original section number kept.*

### 4.8.2 Annotation (2026-10-01): the ~11x was compared against the wrong column, and the model is 14% off the right one

**No number in 4.8 or 4.8.1 above is withdrawn.** What changes is which
measured number the model's 2.59x should have been held against, and the
reason is that Passi et al.'s table has two columns and this chapter had
only ever read one.

Their TLM table reports contact resistance at **two gate biases**: the
on-state (V_BG = −40 V) and the **Dirac point**. The repository ingested
both arrays in 2026-08 and read only the on-state one; the Dirac-point
array was a cited, unread constant until `graphene_dead_name_sweep.py`
found it (`PASSI_RC_DIRAC_OHM_UM`, dead for ten sessions). Reading it:

| hole D (nm) | Rc on-state (Ω·µm) | Rc at Dirac (Ω·µm) | ratio |
|---|---|---|---|
| 0 (unpatterned) | 519 | 1372 | 2.64 |
| 50 | 212 | 620 | 2.92 |
| 100 | 352 | 732 | 2.08 |
| **200** | **45** | **456** | **10.13** |
| 500 | 410 | 1354 | 3.30 |
| 1000 | 560 | 1590 | 2.84 |

The ratio is the **gate-tunable multiple** of the contact resistance. It
sits between 2.08 and 3.30 at five of the six diameters and jumps to
**10.13** at D = 200 nm — the one diameter whose 11x this chapter quotes,
and **3.7x the mean of the other five**.

**So the headline reduction is a gate-bias statement, not a geometry
statement.** D = 0 → D = 200 nm is a **11.53x** reduction in the on state
and a **3.01x** reduction at the Dirac point, for the same two devices and
the same etched geometry.

**And 3.01x is the comparison this model's own physics selects.**
`R_extra` is computed from the *contact-induced* carrier density alone:
the metal's Fermi-level shift sets it, and the gate does not enter. At the
Dirac point the gate contributes no channel carriers, so the measured
access resistance is dominated by exactly the contact-doping-limited term
the model isolates; in the on state a −40 V back gate floods the channel
with holes and the measured Rc reflects a different balance entirely.
Holding a gate-independent model against an on-state measurement is a
category error, and it is the error 4.8 above makes.

Against the right column the model's **2.59x** (982.1 → 379.0 Ω·µm) is
within **14%** of the measured **3.01x**, where against the on-state
column it was short by a factor of **4.5**.

**Two caveats, stated rather than absorbed.** (i) The model's 2.59x is a
*pure-mode* ratio — fully top-mode against fully edge-mode — while the
D = 200 nm device is a mixture at a fill fraction Passi et al. do not
report (4.8.1). Treating the two as comparable assumes the patterned
contact is essentially edge-dominated at that diameter, which the
ligament-overlap criterion of 4.8.1 supports but does not measure.
(ii) The agreement is a *ratio* agreement; Section 4.7's finding that
`R_extra` runs higher than measured Rc in absolute terms is untouched by
this and still open.

**Consequence for 4.8.1's open item.** The small-D upturn the geometry
model fails to reproduce was attributed above to a fabrication effect
"outside this compact model's scope". That attribution is now much better
supported, and more specifically: whatever produces the D = 200 nm optimum
is **3.5x more of an effect in the on state than at the Dirac point**
(4.71x against 1.36x, taking D = 50 nm as the reference), so it scales with
*carrier density*, not with interface area. A model built entirely from
interface geometry — which 4.8.1's is, edge length per unit area at fixed
fill fraction — cannot produce a minimum whose depth depends on the gate.
Its failure was not a missing refinement; it was the wrong class of model
for that feature, and the evidence was in the column the chapter cited and
did not read.

Source: `graphene_edge_contact_model.py`,
`gate_tunable_fraction_of_contact_resistance()` and `dirac_point_summary()`;
sweep `graphene_dead_name_sweep.py`, log `dead_name_sweep_output.txt`.

<a id="s4-10"></a>

*Moved from Chapter 4, §4.10. Original section number kept.*

### 4.10 Conditioning: which Chapter 4 numbers inherit their inputs' errors, and by how much

This section applies to Chapter 4 a diagnostic developed in Chapter 6. It
adds no physics and changes no input; it asks only how each of this
chapter's computed numbers would respond if one of its inputs were wrong.

<a id="s4-10-1"></a>

#### 4.10.1 Why the question is not rhetorical

Section 6.12.5 measured a **65-fold** difference in how two classes of
quantity respond to the *same* 17.12 % input perturbation: a same-sign
(near-cancelling) contact pair moved 77.94 %, a straddling (reinforcing)
pair 1.21 %. The physics of the two cases is identical. What differs is that
the first is a *difference of two comparable numbers*, so its fractional
error is the inputs' fractional error multiplied by the ratio of input
magnitude to output magnitude.

That is arithmetic, not photodetector physics, and Chapter 4's central
recalibration result is a difference of two comparable numbers:

    R_transmission(metal) = R_c,measured(metal) − R_extra(metal)          (4.21)

<a id="s4-10-2"></a>

#### 4.10.2 The diagnostic

For a quantity written as a sum of additive terms `Q = Σ_i T_i`, the signed
logarithmic sensitivity of `Q` to term `T_i` is exactly `S_i = T_i / Q`,
read as *a 1 % change in `T_i` produces an `S_i` % change in `Q`*, and

    Σ_i S_i = 1        exactly                                           (4.22)

by Euler's homogeneous-function theorem. The **condition number** is
`κ(Q) = max_i |S_i|`, and the classification used here is `κ ≤ 1`
reinforcing, `1 < κ < 3` mildly ill-conditioned, `κ ≥ 3` near-cancelling
(Chapter 6's same-sign pairs sit at 4.55).

Equation (4.22) is the implementation's exact validation: it held to
`1.8 × 10⁻¹⁵` over all fifteen decompositions audited across Chapters 4, 5
and 6. Three further exact checks are reported in Section 11 of the
2026-09-22 note, and a fourth — three independent exact zeros — is described
in Chapter 5, Section 5.5.

<a id="s4-10-3"></a>

#### 4.10.3 Result: Section 4.7's four metals, conditioned

| metal | residual (Ω·µm) | `S(R_c)` | `S(R_extra)` | `κ` | class |
|---|---|---|---|---|---|
| Pd | +50.9 | +11.462 | −10.462 | **11.46** | near-cancelling |
| Au | −90.2 | −5.756 | +6.756 | **6.76** | near-cancelling |
| Ni | −360.5 | −1.304 | +2.304 | 2.30 | mildly ill-cond. |
| Cu | −787.5 | −0.234 | +1.234 | 1.23 | mildly ill-cond. |

A prediction written into the note before the audit ran — that **at least
three** of the four would be near-cancelling — is **falsified**: only two
are. The hand arithmetic behind it had been done for Pd and Au and
generalised from those two; Ni and Cu were never computed. This is the third
consecutive session in which a claim generalised from a partial enumeration
failed on the full one, and the first in which the failure was in a
*prediction about* the model rather than in the model.

<a id="s4-10-4"></a>

#### 4.10.4 Result: κ and sign-robustness run in opposite order

`κ` says how an input error is amplified. The question Section 4.7 actually
needs answered is the inverse and sharper one: **how wrong would `R_extra`
have to be to flip the *sign* of the residual?** For (4.21) that fraction is

    f = (R_extra − R_c) / R_extra = −R_transmission / R_extra             (4.23)

| metal | residual | `κ` | `R_extra` must move by | sign is |
|---|---|---|---|---|
| Pd | +50.9 | 11.46 | **−9.6 %** | **fragile** |
| Au | −90.2 | 6.76 | +14.8 % | intermediate |
| Ni | −360.5 | 2.30 | +43.4 % | intermediate |
| Cu | −787.5 | 1.23 | **+81.1 %** | **robust** |

The two orderings are **exactly reversed**, and that reversal carries the
section's two conclusions.

**First, the additive decomposition fails robustly rather than marginally.**
Section 4.7 reported three failures and one survivor and hedged
accordingly. The survivor is the weakest row in the table: a 9.6 % error in
a quantity computed from a model with an admittedly uncertain `λ_decay`
erases Pd's positive residual entirely. Read with its conditioning, the
table does not contain one plausible case and three implausible ones; it
contains one case that establishes nothing and three that establish
something, the strongest being the one that looked worst.

**Second, amplified input error is eliminated as the cause of the negative
residuals** — one branch of an item open since 2026-08-31. Had the negatives
been an artefact of a large `κ` acting on a mis-specified input, the largest
negatives would sit at the largest `κ`. They sit at the smallest: Cu is both
the largest violation and the best-conditioned row. The cause must therefore
lie in the model's structure, which leaves Section 4.7's own two hypotheses
(TLM double-counting the near-contact region; `λ_decay = 250 nm`
over-integrating) standing and removes a third that had not been named. This
is a candidate eliminated by measurement rather than by argument, and it is
the first narrowing of that open item since it was opened.

<a id="s4-10-5"></a>

#### 4.10.5 What this section does not claim

1. It does **not** claim any Chapter 4 number is wrong. `κ` converts an
   input error into an output error and is silent on whether an input error
   exists.
2. The 5 % perturbation used for the illustrative error bars is
   **hypothetical**: none of the four literature `R_c` values is quoted with
   an uncertainty. The resulting bars (Pd 57 %, Au 34 %, Ni 12 %, Cu 6 %)
   inherit that status.
3. A second pre-registered prediction — that some *published* Chapter 4 or 5
   number would carry a >100 % bar under a 5 % input error — is also
   **falsified**. The worst is Pd's 57 %. `κ` is large exactly where the
   output is small, so large `κ` here threatens signs, not orders of
   magnitude; the sign-flip margin of (4.23), not the error bar, is the right
   instrument and it was not the one predicted.
4. Sections 4.5, 4.6 and 4.8 are **not** audited here. Only the Section 4.7
   decomposition and Chapter 5's two families were, and Section 4.8's
   geometry model in particular contains a ratio of two computed resistances
   that has never been conditioned.

Figure: `sensitivity_audit.png`, panel (a). Machine output:
`audit_output.txt`. Pre-registration and outcome:
`notes/2026-09-22-condition-number-audit-of-chapters-4-and-5.md`.

**Annotated 2026-10-23.** Section 4.6.18's closing claim — that what governs
composition is a bipartite write/read relation rather than a set of names —
is **correct and was a census**. It enumerated the relation; it did not
refuse a composition containing an edge. Section 4.6.19 makes it a guard and,
in doing so, shows that 4.6.17's defect and 4.6.18 §4's `gfet.C_ox` are **two
edges of the same relation** removed by **one** change. The numbers in
4.6.18 stand unrevised; its status does not.

---

<a id="s4-11-5a"></a>

*Moved from Chapter 4, §4.11.5a. Original section number kept.*

#### 4.11.5a Annotation (2026-11-03): the crossover quoted above is a *hole*-branch number, and "dies off as it screens" is conditional on an assumption this chapter never stated

Two things in §4.11.5 need correcting in place. Neither changes a committed
number in this chapter — §4.11.5 already records that the shipped model sits on
the `r = 0` diagonal — but both change what the surrounding prose is entitled to
claim.

**(i) `p ≈ 3 × 10¹² cm⁻²` is not a material constant, and it is 13 % too
high.** §4.11.5 quotes it as though it were a property of graphene. It is not.
`graphene_impurity_limited_crossover.py` (2026-11-02) showed that the quantity
invariant under charged-impurity density is the *product* `A = n_c·μ̄_L`, so
`n_c` is proportional to the charged-impurity density and moves with the gate
stack. `graphene_kappa_impurity_bound.py` (2026-11-03) then bounded the
proportionality constant from Zion *et al.*'s own two numbers — their crossover
density and the charged-acceptor density `N_A^(S) = C_ox·V_NP/e` they estimate
for the same sample — giving `κ ≤ 2.089` against the `κ = 2.4` the 11-02 module
had obtained by pairing Zion's crossover with *this* repository's mobility.
`κ = 2.4` is therefore **excluded**, and every `n_c` derived from it, including
the `3 × 10¹² cm⁻²` quoted above, is high by at least 13 %. The superseded
number is left in the sentence above rather than edited out; the corrected
reference value is `2.61 × 10¹² cm⁻²`.

**(ii) The clause "dies off as the long-range impurity potential screens" is
conditional, and the condition has never been checked.** Zion read their
`3 × 10¹² cm⁻²` off the **hole** branch (their Section III c). The
Matthiessen form this repository uses applies a *single* `n_c` to both
branches. Writing it per branch,

&nbsp;&nbsp;&nbsp;&nbsp;`1/μ_b(n) = 1/μ_{L,b} + |n|/(n_{c,b}·μ̄_L)`

and parameterising a split as `n_{c,e/h} = n_c(1 ± s/2)`, there is an exact
identity at `s = 2r₀`: the common factor `(1 + |n|/n_c)` cancels out of
`r = (μ_e − μ_h)/(μ_e + μ_h)`, and

&nbsp;&nbsp;&nbsp;&nbsp;**if `n_{c,b} ∝ μ_{L,b}` then `r(n) = r₀` at every density, exactly, and the
asymmetry does not screen at all.**

So the sentence above does not state a consequence of the model; it states a
consequence of the model *plus* the assumption that the neutral-defect channel
is branch-**symmetric** while the charged channel is not. That assumption was
never argued — it was made by writing one `n_c` — and it is not the obviously
correct one: `r₀` is derived from a beyond-Born cross-section ratio, which is a
property of scattering off a potential of a given sign, and vacancies and grain
boundaries are also potentials.

The magnitude is not marginal. At the top-gated operating point
(`|n|/n_c = 5.29`) the shipped `s = 0` answer is `r = 0.0160`, a 6.2× screening
of `r₀ = 0.10`; at `s = 0.10` it is `0.0580`; at `s = 2r₀` it is `0.1000`
exactly. A 20 % split — modest for a quantity that has never been measured —
moves `r` by 524 %, and the endpoint is not a perturbation of the shipped
answer but its negation.

**What this costs the chapter.** Nothing in §4.3–§4.10 moves, for the reason
§4.11.5 already gives: the shipped model is the `r = 0` case. What moves is
§4.11.5's own framing. "A density-dependent pair is logged as an open item
rather than guessed at" is still the right call, but the open item is now
strictly larger than it was written to be: it is not only *which* pair, it is
whether the density dependence has the sign and the asymptote this chapter
attributes to it. Discharging it needs the **electron**-branch transition
density, which is in Zion Fig. 5 next to the hole one and has never been read.

Source: `graphene_kappa_impurity_bound.py`, 16/16 checks green, six against
exactly known values. Figure: `kappa_impurity_bound.png` (panel (c) is the
identity above). Machine output: `kappa_impurity_bound_output.txt`.
Derivation and literature: `notes/2026-11-03-the-bound-is-in-the-papers-own-numbers.md`.

## A.6 Chapter 5: Interconnects

<a id="s5-5"></a>

*Moved from Chapter 5, §5.5. Original section number kept.*

### 5.5 Conditioning: the calibration is the worst-conditioned step in this thesis

Section 4.10 introduces the condition-number diagnostic and its exact
validation, Eq. (4.22). This section applies it to Chapter 5. It changes no
number and recalibrates nothing; the results below are error bars that were
previously unstated, not corrections.

<a id="s5-5-1"></a>

#### 5.5.1 The formula is reinforcing; the calibration behind it is not

`resistivity_vs_width()` is a Matthiessen sum of three strictly positive
terms, and audited against its own arguments it behaves exactly as such:

| W (nm) | ρ (µΩ·cm) | S(bulk) | S(edge) | S(impurity) | κ |
|---|---|---|---|---|---|
| 18 | 4.093 | 0.293 | 0.662 | 0.045 | 0.662 |
| 22 | 3.600 | 0.333 | 0.616 | 0.051 | 0.616 |
| 30 | 3.009 | 0.399 | 0.541 | 0.061 | 0.541 |
| 40 | 2.602 | 0.461 | 0.469 | 0.070 | 0.469 |
| 52 | 2.321 | 0.517 | 0.404 | 0.079 | 0.517 |

Every `κ ≤ 0.66`: reinforcing at every width, as a sum of positive terms
must be. Read at that level, Chapter 5 contains no near-cancellation at all.

**But `λ_impurity` is not an argument of that formula — it is *solved* from a
single calibration point, and the solve is a subtraction:**

    residual = ρ_cal/ρ_bulk − 1 − [(1−p)/(1+p)]·(λ_bulk/W_cal)
             = 3.0000 − 1.0000 − 1.8478  =  0.1522                       (5.9)
    λ_impurity = λ_bulk / residual = 361.4 nm

Three terms of order unity producing a result of order 0.15 is the textbook
signature, and the audit measures

    κ(residual) = 19.71                                                  (5.10)

**larger than anything in Chapter 4 (11.46) or Chapter 6 (4.55). Equation
(5.9) is the worst-conditioned step in this thesis.** The parameter-level
sensitivities follow: `d ln λ_imp / d ln ρ_cal = −19.71`,
`d ln λ_imp / d ln λ_bulk = +13.14`, `d ln λ_imp / d ln W_cal = −12.14`,
`d ln λ_imp / d ln p = −3.73`. **A 5 % error in the single calibration datum
(ρ = 3.6 µΩ·cm at W = 22 nm, "best GNR" from Murali et al.) moves the
extracted impurity mean free path by 99 %.**

<a id="s5-5-2"></a>

#### 5.5.2 Result: the near-cancellation propagates, damped but not removed

Re-solving `λ_impurity` inside the audit — the chapter's *actual* dependency
chain, rather than the formula read in isolation — gives:

| W (nm) | S(ρ_cal) | S(ρ_bulk) | S(p) | S(W_cal) | S(λ_bulk) | κ |
|---|---|---|---|---|---|---|
| 18 | 0.880 | +0.120 | −0.037 | 0.542 | +0.120 | 0.880 |
| 22 | **1.000** | **0.000** | **0.000** | 0.616 | **0.000** | 1.000 |
| 30 | 1.197 | −0.197 | +0.060 | 0.737 | −0.197 | 1.197 |
| 40 | 1.384 | −0.384 | +0.118 | 0.852 | −0.384 | 1.384 |
| 52 | 1.551 | −0.551 | +0.169 | 0.956 | −0.551 | **1.551** |

The `κ ≈ 20` of (5.9) is damped by the small weight the impurity term
carries in the total ratio, but not removed: `κ(ρ)` crosses 1 at the
calibration width and reaches **1.551** at 52 nm. A prediction registered
before the audit — that `κ` would exceed 1 but stay below the
near-cancelling threshold of 3 — holds.

The methodological point is general and worth stating in Chapter 7's voice:
**a quantity can be well-conditioned in its own arguments and
ill-conditioned in the quantities those arguments were derived from.**
Auditing the visible formula alone would have returned `κ ≤ 0.66` and a
clean bill of health.

<a id="s5-5-3"></a>

#### 5.5.3 An exact result: the calibration point pins the model, and S(ρ_bulk) changes sign through it

The bold row above is not rounding. At `W = W_calibration = 22 nm` the model
reproduces `ρ_calibration = 3.6 µΩ·cm` **bitwise**, and

    S(ρ_bulk) = S(p) = S(λ_bulk) = 0.0        exactly                    (5.11)

— three independent exact zeros, verified as such and not merely small
(`conv = 0.0e+00`). The reason is structural: at the calibration width the
impurity term is defined as whatever makes the total come out to the
calibration datum, so it absorbs any change in the other three inputs
exactly. **At 22 nm this chapter's prediction carries no information from
`ρ_bulk`, `p` or `λ_bulk` whatsoever.**

`S(ρ_bulk)` then changes sign through that width — `+0.120` at 18 nm,
`0.000` at 22 nm, `−0.551` at 52 nm — because on the wide side an increase in
`ρ_bulk` raises the modelled bulk term but lowers the fitted impurity term by
more. A one-point-calibrated model inherits a sensitivity structure that
pivots about its calibration point, and **the further a prediction is made
from 22 nm, the more of its content comes from the single calibration datum
and the less from the physical inputs** — which is the opposite of how a
calibrated model is usually read. Section 5.3's headline comparisons are
drawn at the narrow end (18–22 nm), where this effect is weakest; the
wide-end numbers carry it.

<a id="s5-5-4"></a>

#### 5.5.4 Result: the Cu liner model contains an unremarked subtraction

`W_eff = W − 2t` (Section 5.2) is a difference of comparable quantities:

| W (nm) | W_eff (nm) | S(W) | S(2t) | κ(W_eff) | \|S_t(ρ_eff)\| |
|---|---|---|---|---|---|
| 18 | 12.0 | 1.500 | −0.500 | **1.500** | **1.167** |
| 22 | 16.0 | 1.375 | −0.375 | 1.375 | 0.852 |
| 30 | 24.0 | 1.250 | −0.250 | 1.250 | 0.550 |
| 40 | 34.0 | 1.176 | −0.176 | 1.176 | 0.379 |
| 52 | 46.0 | 1.130 | −0.130 | 1.130 | 0.276 |

`κ` grows as `W → 2t` and is mildly ill-conditioned across the whole
modelled range. At the narrow end a 1 % error in the liner thickness is a
1.17 % error in the reported effective resistivity. Section 5.2 takes
`t = 3.0 nm` from two sources that themselves differ — 3 nm for a Cu/TaN-Co
stack (*Nanomaterials* **12**(10), 1760 (2022)) against a "2–3 nm functional
floor" (Domenichini *et al.*, arXiv:2406.09106) — a spread of roughly 17 %,
which (5.12)

    δρ_eff/ρ_eff ≈ |S_t| · δt/t                                          (5.12)

converts into a **~20 % bar on ρ_eff at 18 nm**, previously unstated. This
is the one input in either chapter whose uncertainty is available from the
literature rather than hypothesised, and it is also the width range where
Chapter 5's graphene-vs-Cu crossover argument is made.

<a id="s5-5-5"></a>

#### 5.5.5 What this section does not claim

1. No number in Sections 5.2–5.4 is changed, retracted or recalibrated. The
   calibration datum is still `3.6 µΩ·cm at 22 nm` and `λ_impurity` is still
   361.4 nm.
2. `κ(ρ_GNR)` grows monotonically with `W` over the audited 18–52 nm, so a
   wider wire would be *worse* conditioned. Beyond 52 nm this is untested.
3. The parallel-conduction refinement (Section 5.3, thin Ru/Co liners) and
   the `cu_resistivity_vs_width()` Fuchs–Sondheimer baseline are **not**
   audited here. The baseline in particular contains its own `p_cu = 0.6`
   with no stated uncertainty.
4. `κ ≈ 20` on (5.9) is not evidence that `λ_impurity = 361 nm` is wrong. It
   is evidence that the value is determined by one datum to a precision
   twenty times worse than that datum's own — a statement about what the
   chapter is entitled to claim, not about what is true.

Figure: `sensitivity_audit.png`, panels (b) and (c). Machine output:
`audit_output.txt`. Pre-registration and outcome:
`notes/2026-09-22-condition-number-audit-of-chapters-4-and-5.md`.

## A.7 Chapter 6: Photodetectors

<a id="s6-8-note"></a>

*Annotation moved from the body of Chapter 6, §6.8.*

> **Annotation added 2026-09-21 (Section 6.12).** The single
> `W_CROSS_CHEM = 5.4 eV` used throughout this section and Sections
> 6.9–6.11 is now known to be the crossover *at the physisorbed separation*
> rather than a constant shared by all metals. Nothing in this section is
> retracted and no number below is changed, but every ΔW here carries an
> uncertainty of up to **17%** (and every same-sign pair up to **78%**) from
> that assumption alone — quantified, for the three metals it can be
> quantified for, in Section 6.12.3 and 6.12.5.

<a id="s6-9-note"></a>

*Correction moved from the body of Chapter 6, §6.9.5.*

> **RETRACTED 2026-09-20 — see Section 6.11.2.** The two numbers above are
> correct. The generalisation from them is not. Enumerating all 21
> asymmetric pairs under this same model gives **14 violations of the
> "under 1.02" claim, the worst being Au/Pd at 22.4**. The sample of two
> was drawn entirely from pairs that *straddle* the 5.4 eV crossover, and
> every straddling pair does satisfy the bound; same-sign pairs need not,
> because max|k|/|N_uniform| diverges as N_uniform → 0. The final sentence
> — "masks are for symmetric devices only" — is therefore **false**: a
> perfect mask gains 8.4× on Au/Pd and more than 2× on five pairs, per
> incident photon. What survives is the *practical* recommendation, on
> different grounds, in Section 6.11.2. The original numbers are left in
> place above rather than rewritten.

<a id="s6-11-note"></a>

*Annotation moved from the body of Chapter 6, §6.11.2.*

> **Annotation added 2026-09-22 (Section 6.13.7).** The count *five* holds
> only at exactly `w_cross = 5.4 eV`. Re-scored across the scalar-crossover
> band `|δ| ≤ 0.2 eV` it runs **3 → 4 → 5 → 5 → 6**: at `δ = −0.2 eV` Ni/Au
> and Ni/Pd fall below 2×, and at `δ = +0.2 eV` Ti/Cu rises above it. The
> numbers in this table are **not** retracted — they are correct at `δ = 0` —
> but the *count* is, because it is a threshold evaluated at one point of a
> parameter both source papers give as "~5.4 eV". The qualitative claim below
> survives at every offset.

<a id="s6-12-note"></a>

*Annotation moved from the body of Chapter 6, §6.12.3.*

> **Audited 2026-09-24 -- the number stands, and it nearly did not deserve to.**
> The bisection that produced this boundary was, from the day it was written
> until 2026-09-24, one of the *unguarded* root-finders whose fault class
> Section 6.13.6 put on record: it ran 200 halvings with no check that its
> interval contained a sign change. It has now been audited in full
> (`graphene_rootfinder_audit.py`, seven validations, 7/7). At the values
> actually used here the interval **is** bracketed -- `f(lo) = +1.371`,
> `f(hi) = -0.090` -- so 0.4997 A is a genuine root, with residual
> `worst(root) - tol = -9.6e-16`, and no number in this section moves. The
> fault was latent rather than realised: on a deliberately root-free interval
> the same loop returns 5.00 A, or one ulp above 0.05 A, with equal composure.
> The audit also converted this section's monotonicity assumption -- that the
> |dW| shift decreases in ell, without which a bracketed bisection could still
> return the wrong one of several roots -- from an assertion into a
> term-by-term proof, with Pt's identically zero term (the paragraph below)
> as its exact test. And it found that the *fix* initially degraded this
> number by ten significant digits, because the shared guard's default
> tolerance is absolute and had been calibrated for a variable of order 1 eV
> rather than 1e-10 m; see
> `notes/2026-09-24-root-finder-bracket-audit.md`, Section 4. **The honest
> reading is that this section's headline survived an audit it had a real
> chance of failing, and that the audit's own repair had to be audited.**

<a id="s6-6"></a>

*Moved from Chapter 6, §6.6. Original section number kept.*

### 6.6 Spatially resolved, metal-dependent photocarrier collection

Section 6.4's model applies a single literature midpoint, EQE_bare =
0.15%, regardless of contact metal. But the "collection bottleneck"
described in Section 6.2 -- only photocarriers generated within a
built-in-field region near a contact survive graphene's ~1 ps carrier
lifetime to be collected -- is exactly the same contact-induced-doping
physics already modeled quantitatively for Chapter 4's contact
resistance work (`graphene_contact_doping_model.py`, Section 4.5). This
section closes the gap flagged since this chapter was first drafted
(2026-08-24): integrating that spatial doping-profile machinery into the
photodetector side (full literature review and design derivation in
`notes/2026-09-07-spatial-photocarrier-collection-model.md`).

**Model.** `graphene_photodetector_collection_model.py` superposes the
contact-doping-induced field (differentiating the same saturating
profile used for `graphene_contact_doping_model.py`'s density profile)
on the bare-device model's uniform bias field E_bias = V_bias/L_channel
= 5.0e5 V/m, then integrates the resulting position-dependent drift
velocity to get each photocarrier's transit time to the contact, and
from that an exponential-lifetime survival (collection) probability
p(x) = exp(-t(x)/tau_carrier). Averaging p(x) over the channel and
normalizing against a bias-field-only baseline (physically, a
work-function-matched contact -- already present in this thesis's metal
table as Cr, W = 4.50 eV vs. graphene's 4.5 eV) gives a metal-dependent
**collection enhancement factor**, plotted alongside the resulting
metal-resolved modeled EQE in
`photodetector_collection_efficiency_by_metal.png`.

**Mechanism support from the literature.** Xia et al. (*Nature
Nanotechnology* 4, 839 (2009)) demonstrated zero-external-bias graphene
photodetection using two *different*-work-function metal contacts (Ti,
Pd), which is only possible if contact-metal work function sets a real,
usable built-in field in the channel -- direct experimental support for
this section's core mechanism. Weiss & Duan (*NPG Asia Materials* 5, e74
(2013)) state the mechanism explicitly ("contacting graphene with metals
of variable work function... creates a potential offset... that
facilitate[s] the separation and transport of photocarriers") and make
an important scope-limiting point this model respects: identical
contacts produce "equal and opposing" built-in potentials that partially
cancel in the net terminal photocurrent. Mueller, Xia & Avouris (*Nature
Photonics* 4, 297 (2010)) independently estimate the built-in-field
collection region at ~100-200 nm, consistent in order of magnitude with
(though not identical to) the lambda_decay = 250 nm decay length reused
here from `graphene_contact_doping_model.py` (itself sourced from
Khomyakov et al. 2010's DFT-derived doping-decay length) -- a useful
independent cross-check, not a perfect match, and reported as such.

**Results.** Sweeping the seven contact metals already in this thesis's
`METAL_WORK_FUNCTIONS` table, the modeled collection-efficiency
enhancement ranges from 1.00x (Cr, the baseline) to 1.47x (Pt, the
largest work-function mismatch, |dV| = 1.15 V), i.e. a modest,
physically plausible boost -- not a multiple-order-of-magnitude effect,
since the doping field only dominates within roughly one lambda_decay of
the contact and the 200 nm channel modeled here is comparable in length
to that decay length. Commonly used real contact metals fall in between:
Ti and Cu (both common adhesion/contact metals) give ~1.21-1.23x; Ni,
Au, and Pd (higher work function, used specifically for their
historically low contact resistance -- Section 4.5/4.7) give ~1.39-1.41x.
This offers an additional physical rationale, alongside the
injection-resistance/quantum-capacitance reasoning of Chapter 4, for why
the high-work-function metals already favored for low contact resistance
are also a reasonable choice for photodetector contacts specifically.

**Scope limitation (stated explicitly).** This model represents only the
single contact whose doping field reinforces the bias-driven carrier
sweep direction; the real device's second contact, where the doping
field would partially oppose that sweep (per Weiss & Duan's "equal and
opposing" result, generalized to the biased case), is not modeled. A
fully self-consistent two-contact treatment is left as an open item
below rather than approximated away silently.

> **Resolved, and it overturns this section's ranking (2026-09-17).**
> Section 6.7 now models both contacts self-consistently. The result is
> not a small correction to the numbers below: for a *symmetric*
> two-terminal device the metal ranking is essentially reversed, and Pt
> -- this section's best metal at 1.47x -- becomes the worst, overstated
> here by a factor of 5.5. **The enhancement factors reported in this
> section remain valid only for a single junction in isolation, which is
> not the device geometry this chapter is otherwise about.** They are
> retained rather than deleted because the single-junction calculation
> is the correct building block and the limiting case that validates the
> two-contact model, but any device-level conclusion should be taken
> from Section 6.7.

**What this result can and cannot claim.** The relative, metal-to-metal
ordering (larger work-function mismatch -> shorter transit time -> higher
survival probability -> higher effective EQE) follows directly from the
mechanism above and is this section's main result. The *absolute*
EQE_bare = 0.15% literature midpoint this factor scales is not
attributed to any specific contact metal in the material reviewed for
Section 6.2, so the resulting absolute EQE_model(metal) numbers should
be read as relative-to-that-midpoint predictions, not validated absolute
values for each named metal.

A further correction, added 2026-09-17: the relative metal-to-metal
ordering described in the paragraph above -- the part called "this
section's main result" -- is precisely the part Section 6.7 overturns for
a two-terminal device. The chain "larger work-function mismatch ->
shorter transit time -> higher survival probability" is sound for one
junction, but in a symmetric device the same larger mismatch also
produces a stronger *opposing* field at the other contact, and the
second effect wins.

<a id="s6-7"></a>

*Moved from Chapter 6, §6.7. Original section number kept.*

### 6.7 Self-consistent two-contact collection, and the reversal of Section 6.6's metal ranking

*(Heading retained for the record. The "reversal" it names was
downgraded to crossover-dependent by Section 6.8 on 2026-09-18; this
section's Result 2, by contrast, was confirmed and strengthened there.)*

Section 6.6 closed one gap and opened another, which it named: it models
a single contact and assumes every photocarrier is collected there. The
geometry makes that assumption untenable rather than merely approximate.
The modeled channel is L = 200 nm while the contact-doping decay length
reused from Chapter 4 is lambda_decay = 250 nm, so neither contact's
field has decayed appreciably by mid-channel -- the second contact is not
a perturbation on the first, it is comparable to it everywhere.

**Model.** `graphene_photodetector_two_contact_model.py` places contact
A (metal A) at x = 0 and contact B (metal B) at x = L, and defines a
signed total field, positive meaning "sweeps carriers toward contact A":

    F(x) = E_bias + g_A(x) - g_B(x)

where g_A is byte-for-byte the doping-field profile of Section 6.6 and
g_B is its mirror image about mid-channel. The minus sign on g_B is the
entire physical content that Section 6.6 was missing. Because F can
change sign inside the channel, the model assumes no destination: a
carrier drifts along F and is collected at A (signed +1) only if F > 0
along its whole path to x = 0, or at B (signed -1) only if F < 0 along
its whole path to x = L. A carrier whose path crosses a sign change
drifts into a **stagnation point** (F = 0) rather than a contact and is
counted as uncollected -- a real feature of the two-contact geometry,
since with two opposing doping fields and a modest bias there is
generically an interior null. Survival to the contact uses the same
tau = 1 ps photocarrier lifetime as Section 6.6. The figure of merit is
the signed, channel-averaged net response N(A,B), directly comparable to
Section 6.6's collection efficiency.

**Two validations.** Both are computed and printed, not asserted. First,
switching off contact B's doping field (by giving it graphene's own work
function) must collapse the model onto Section 6.6's result, since then
F > 0 everywhere and every carrier reaches A: it does, for all seven
metals, to a worst relative deviation of 7.8e-08. Second, at zero bias
two *identical* contacts must give exactly zero net response, because
F(x) = g(x) - g(L-x) is antisymmetric about mid-channel -- this is Weiss
& Duan's "equal positive and negative flow with a net zero
photocurrent", and Shimomura et al.'s observation that in symmetric devices
"the polarities of the photovoltages at each graphene/electrode
interface ... are canceled out under macroscopic light irradiation". The
model reproduces it to |N| <= 6.6e-17, i.e. machine precision, as an
emergent consequence rather than an input. The antisymmetry of the
metal-pair matrix under contact swap is a third consistency check.

**Result 1: the reversal.** At the working bias V_bias = 0.1 V
(E_bias = 0.5 MV/m), comparing Section 6.6's single-contact collection
efficiency against the symmetric two-contact device's actual net
response:

| metal (both contacts) | W (eV) | single-contact (6.6) | symmetric two-contact (6.7) | ratio |
|---|---|---|---|---|
| Ti | 4.33 | 0.776 | 0.672 | 0.87 |
| Cr | 4.50 | 0.632 | 0.632 | 1.00 |
| Cu | 4.65 | 0.766 | 0.673 | 0.88 |
| Ni | 5.04 | 0.877 | 0.321 | 0.37 |
| Au | 5.10 | 0.885 | 0.295 | 0.33 |
| Pd | 5.12 | 0.888 | 0.287 | 0.32 |
| Pt | 5.65 | 0.929 | 0.169 | 0.18 |

> **Downgraded by Section 6.8 (2026-09-18).** The reversal below is
> **crossover-dependent and is no longer claimed as a result.** Section
> 6.8's signed, carrier-resolved model reproduces this section exactly
> under its own restrictions (holes only, |ΔW|, crossover at 4.5 eV --
> agreement to 0.000e+00 over 98 comparisons), so the arithmetic here is
> right. But placing the n/p crossover at its physical value of 5.4 eV
> instead of graphene's 4.5 eV puts Pt **back at the top** (0.556, best)
> rather than the bottom. Three sections have now given three different
> answers to this one question, and the honest reading is that the
> symmetric metal ranking is not an established result of this thesis.
> The numbers below are retained because the calculation is correct for
> the convention it states; the *ranking* is not. What survives is
> Section 6.8.5's structural statement: in a symmetric device contact
> doping is purely parasitic, so the best metal is the one that dopes
> graphene least.

Section 6.6 ranked the metals Pt > Pd > Au > Ni > Ti > Cu > Cr. For a
symmetric device the ranking is essentially reversed: Cu ~ Ti > Cr > Ni >
Au > Pd > Pt. The mechanism is clear once the sign is right: a strongly
doping contact sweeps carriers toward *itself*, so two of them facing
each other produce a large opposing-field region and a mid-channel
stagnation point that the 0.5 MV/m bias cannot overcome against Pt's
~4.5 MV/m near-contact field. **Strong contact doping helps an isolated
junction and hurts a symmetric two-terminal device.** Cr is the one
metal unaffected (ratio exactly 1.00), because its work function matches
graphene's to within 0.01 eV and so contributes essentially no field to
cancel, leaving the bias to sweep the channel unopposed.

> **Restated by Section 6.8.6 (2026-09-18).** The tension as framed
> below -- that the low-contact-resistance metals are the worst for
> photoresponse -- does not survive the signed model: at the physical
> crossover Pd sits near the top of the symmetric ranking, not the
> bottom. The durable tension is narrower and different in kind:
> Chapter 4 optimises a *single* junction, whereas the two-terminal
> photoresponse depends on the *difference between two* junctions and is
> to first order blind to either one's own quality. Chapter 7 should
> synthesise that version.

This puts two of this thesis's own device-design arguments in tension,
which is worth stating plainly rather than smoothing over. Chapter 4
favours high-work-function metals (Ni, Au, Pd) for low contact
resistance, and Section 6.6 appeared to agree for photocarrier
collection. Section 6.7 shows the opposite for the photodetector: the
metals that are best for contact resistance are the worst for symmetric
two-terminal photoresponse. That is a real trade-off in device design,
not a modelling artefact, and a device-oriented treatment should make
the choice explicit rather than inherit it from the contact-resistance
chapter.

**Result 2: asymmetric contacts, at zero bias.** Because the symmetric
diagonal is exactly zero at zero bias, every nonzero entry isolates the
work-function-asymmetry mechanism:

| A \ B | Ti | Cr | Cu | Ni | Au | Pd | Pt |
|---|---|---|---|---|---|---|---|
| Ti | 0.000 | 0.596 | 0.067 | -0.753 | -0.794 | -0.802 | -0.903 |
| Cr | -0.596 | 0.000 | -0.563 | -0.835 | -0.849 | -0.853 | -0.917 |
| Cu | -0.067 | 0.563 | 0.000 | -0.778 | -0.805 | -0.812 | -0.905 |
| Ni | 0.753 | 0.835 | 0.778 | 0.000 | -0.078 | -0.103 | -0.581 |
| Au | 0.794 | 0.849 | 0.805 | 0.078 | 0.000 | -0.025 | -0.504 |
| Pd | 0.802 | 0.853 | 0.812 | 0.103 | 0.025 | -0.000 | -0.480 |
| Pt | 0.903 | 0.917 | 0.905 | 0.581 | 0.504 | 0.480 | 0.000 |

This is the quantitative version of the Xia et al. (2009) Ti/Pd
zero-bias device that Section 6.6 already cites as experimental support
for the mechanism: a work-function-asymmetric contact pair produces a
net zero-bias photoresponse where a symmetric pair produces none.

One non-obvious feature, worth flagging because it distinguishes a
mechanism model from a monotonic fit in |W_A - W_B|: the largest
response is Cr/Pt (0.917), not Ti/Pt (0.903), despite Ti/Pt having the
larger work-function difference (1.32 eV vs 1.15 eV). Ti carries its own
doping field (0.17 eV of mismatch) that sweeps carriers back toward Ti
and partially opposes Pt's, whereas Cr contributes no opposing field at
all. Under this model, the best zero-bias pairing is one strongly doping
contact against a *work-function-matched* one.

**The simplification that could overturn Result 2.** Both doping fields
are taken as sweeping carriers toward their own contact, through
|W_metal - W_graphene| -- the same magnitude convention Section 6.6 uses,
and the convention that makes identical contacts cancel exactly. But it
erases the n/p distinction: Ti and Cu (W < W_graphene) dope graphene
n-type while Pt, Pd and Au (W > W_graphene) dope it p-type, and a signed
treatment tracking electrons and holes separately would let an n/p pair
such as Ti/Pt *add* rather than partially cancel for one carrier
species. That would plausibly make Ti/Pt the strongest pairing and
invert the Cr/Pt-vs-Ti/Pt conclusion above. **Confirmed in Section 6.8.4
(2026-09-18): Ti/Pt is the strongest pair, |N| = 1.832, roughly double
the 0.917 quoted here, and the result is insensitive to where the
crossover is placed. The "Result 1 does not depend on it" clause that
follows, however, turned out to be wrong -- see the callout at the head
of Result 1 above.** Result 2's ordering should
therefore be read as provisional; Result 1's reversal does not depend on
it, since it concerns identical contacts, where the two conventions
agree. Full discussion in
`notes/2026-09-17-two-contact-self-consistent-collection.md`, Section 4.

Figure: `photodetector_two_contact_net_response.png` -- (a) the two
opposing fields and the Pt/Pt stagnation point at x = 100 nm, (b) the
single-contact vs symmetric two-contact comparison showing the reversal,
(c) zero-bias net response against W_A - W_B.

<a id="s6-14-7"></a>

*Moved from Chapter 6, §6.14.7. Original section number kept.*

### 6.14.7 What this section does and does not change

**Changes.**
- Section 6.13.5's sign-robustness result is **strengthened in kind**: it was
  a 12621-point scan and is now a parity identity, valid without an interval.
- Section 6.13.10's closing remark — *"a per-metal crossover difference of
  0.02 eV would flip Au/Pd"* — is confirmed **exactly**, and generalised: the
  flipping offset is `W_B − W_A` for every pair.
- The chapter now has an explicit statement of which of its sign claims are
  safe (all of them against common offsets; three of twenty-one against
  differential ones at the modelled magnitudes) and which are unevaluable
  (the four smallest margins).

**Does not change.** No number in Sections 6.8–6.13 is retracted or altered.
Predictions P2 and P4 of this session's own note are falsified; they were
predictions about this model, not results of an earlier section.

<a id="s6-10"></a>

*Moved from Chapter 6, §6.10. Original section number kept.*

### 6.10 Further follow-on work

- ~~Model both contacts self-consistently~~ -- **done 2026-09-17,
  Section 6.7**, and it reversed Section 6.6's metal ranking for a
  symmetric device.
- Replace Section 6.7's |W_metal - W_graphene| magnitude convention with
  a signed, carrier-resolved treatment that distinguishes n-type
  (Ti, Cu) from p-type (Pt, Pd, Au) contacts and tracks electrons and
  holes separately. This is now the single most consequential open item
  in this chapter: it could invert Section 6.7's Result 2 ordering by
  making an n/p pair such as Ti/Pt add rather than partially cancel.
- ~~Model non-uniform illumination (a generation weight g(x))~~ --
  **done 2026-09-19, Section 6.9**, against Shimomura et al.'s shadow-mask
  device. It required no new transport physics (illumination is a weight on
  a fixed kernel) and it produced the chapter's first *bound*, max|k|,
  rather than another ranking.
- ~~Replace the linear dW -> doping profile relation~~ -- **done
  2026-09-20, Section 6.11**, and it both confirmed Section 6.8's Ti/Pt
  headline (4.9% compression) and retracted Section 6.9's ceiling claim and
  its "masks are for symmetric devices only" conclusion.
- **Implement a per-metal crossover w_cross = W_graphene + Delta_c(d)**, now
  the top open item in this chapter. Section 6.11.3 item 3: holding the
  crossover at 5.4 eV for all seven metals is the last unexamined constant
  in the chain, and unlike the magnitude nonlinearity it would move every
  dW in the table, not just compress them.
- **Find a description of Ti, Ni and Pd that does not go through their work
  function.** Section 6.11.1 excludes four of seven metals from Eq. 7, and
  Ti is contact A of both headline pairs -- so the chapter's central numbers
  now rest on the one assumption the literature most explicitly disowns.
- Reconcile Chapter 4's contact-metal recommendation with Section 6.7's:
  the two chapters now point in opposite directions for the same metals,
  and Chapter 7 (discussion/outlook) should resolve this explicitly
  rather than leaving the reader to notice it.
- Integrate Section 6.5's plasmonic near-field picture with the
  spatially-resolved contact-doping machinery now used in both Section
  6.6 and Chapter 4 (`graphene_contact_doping_model.py`,
  `graphene_edge_contact_model.py`) -- a real device's plasmonic hot-spot
  and contact depletion region sit at related but distinct locations,
  which Section 6.5's spatially-uniform EQE multiplier does not capture.
- Model the photo-bolometric mechanism (dominant in the telecom bowtie/
  waveguide device reviewed in Section 6.5) as its own responsivity
  model -- dR/dT of graphene's resistance vs. absorbed power under bias
  -- rather than treating photogating as the only gain mechanism this
  thesis represents.
- Find or derive a measured (rather than assumed) plasmon resonance
  quality factor/FWHM for either Section 6.5 design, to replace the
  assumed Q = 6-8 values with literature-anchored ones.
- Revisit Park, Ahn et al.'s scanning-photocurrent-microscopy sign-map
  result (PubMed 19326919), which WebFetch could not retrieve this
  session (HTTP 429), as a more direct cross-check of Section 6.6's
  two-contact sign-reversal reasoning than the two sources it currently
  relies on (Xia et al. 2009, Weiss & Duan 2013).

<a id="s6-13"></a>

*Moved from Chapter 6, §6.13. Original section number kept.*

### 6.13 How far can the crossover move? A sensitivity band for Sections 6.8–6.12

Section 6.12 replaced the single 5.4 eV crossover with a per-metal one, and
in doing so raised a question it could not answer: **what does the rest of
the chapter look like if 5.4 eV is simply wrong by some amount?**

That question has to be asked separately, because the per-metal machinery
**cannot reach the pairs that matter**. It needs the metal–graphene
separation `d_eq`, which Khomyakov *et al.* tabulate for three of this
thesis's seven metals (Cu, Au, Pt). Ti, Ni and Pd are chemisorbed, where
Section 6.12.4 showed the anchored exponential diverges; Cr has no tabulated
separation at all. And every one of the five headline pairs of Section
6.11.2's retraction table — Au/Pd, Ti/Cr, Ni/Au, Cr/Cu, Ni/Pd — contains at
least one metal outside that set of three.

<a id="s6-13-1"></a>

#### 6.13.1 The perturbation, and why a scalar one is honest

Both papers the 5.4 eV comes from write it with a tilde:

> "the crossover from p-type to n-type doping occurs for a metal work
> function of **~5.4 eV**"
> — Giovannetti *et al.*, PRL **101**, 026803 (2008); repeated verbatim in
> Khomyakov *et al.*, PRB **79**, 195425 (2009).

Sections 6.8–6.12 have read that tilde as three significant figures. This
section puts an unknown scalar offset `δ` on the crossover, `w_cross = 5.4 +
δ`, and sweeps it. The offset has the *wrong shape* — it is not per-metal —
but the *right reach*: it needs no `d_eq`, so it applies to all seven metals.
The two treatments are complementary and neither supersedes the other.

It is also algebraically identical to a **common error in every entry of the
work-function table**, since only `dW = W − w_cross` enters. That is not a
hypothetical: Ni is entered at 5.04 eV from a range commonly quoted as
4.9–5.35 eV, and the 0.12 eV potential step Mueller *et al.* measure has been
an unreconciled discrepancy in this thesis since Section 6.8. A nominal band
of `|δ| ≤ 0.2 eV` is adopted, with a ±0.5 eV outer sweep; no result below
depends on that choice, because every result is reported as a function of δ.

<a id="s6-13-2"></a>

#### 6.13.2 Validation against five exactly known values

Four were planned; the fifth was found while trying to score a prediction and
turned out to be the section's main result.

| # | check | outcome |
|---|---|---|
| 1 | shift invariance: `N(W; 5.4+δ) == N(W−δ; 5.4)` | 215/231 bitwise, worst **1.78e-15** |
| 2 | `δ = 0` reproduces the Section 6.8 model | **bitwise**, 21/21 pairs |
| 3 | symmetric pair → `N == 0` at every δ | worst **1.3e-16**, 77 cases |
| 4 | charge conjugation `N(−dW) == −N(dW)` at every δ | **exactly 0.000e+00**, 539/539 |
| 5 | a straddling pair's `⟨\|dW\|⟩` is independent of `w_cross` | **exactly 0.000e+00**, 126 cases |

Check 1 is exact in real arithmetic but not in floating point, and the note
predicted that in advance — along with a bound of 1e-15, which the measured
1.78e-15 **misses by a factor of 1.8**. The bound was the wrong shape: the
right statement is a few ulp of `N` (1.78e-15 is 8 ulp of 0.593), not a fixed
absolute number. It is recorded because a pre-registration practice is worth
nothing if only the comfortable misses are recorded.

<a id="s6-13-3"></a>

#### 6.13.3 Result: the amplification ratio of Section 6.12 is not a property of a pair

Section 6.12.5's headline was a **65-fold** difference in how hard a
same-sign pair and a straddling pair amplify an input uncertainty, measured
on three pairs. Repeating that measurement here gives, for all six straddling
pairs, an amplification of **infinity** — not large, exactly undefined.

Check 5 is why. For a pair straddling the crossover,

    ⟨|dW|⟩ = [ (w_cross − W_A) + (W_B − w_cross) ] / 2 = (W_B − W_A) / 2

in which `w_cross` has cancelled identically. Moving the crossover changes
such a pair's response (`S` below is small but nonzero) while changing the
input measure by **exactly nothing**, so the ratio divides by zero.

**The amplification ratio is a property of the (pair, perturbation) couple,
not of the pair.** Section 6.12.5's number should not be quoted as a device
property, and this chapter came close to doing so.

<a id="s6-13-4"></a>

#### 6.13.4 Result: the dichotomy is real, at 30× rather than 65×

The statistic that survives is the logarithmic sensitivity
`S = d ln|N| / dδ`, which is well defined for every pair because it does not
divide by the input:

| group | range of `\|S\|` (eV⁻¹) | pairs |
|---|---|---|
| straddling | 0.0077 – 0.0218 | Ti/Pt, Cr/Pt, Cu/Pt, Ni/Pt, Au/Pt, Pd/Pt |
| same-sign | 0.6625 – 2.5628 | the other 15 |

`min(same-sign) / max(straddling) = **30.4**`, with no overlap anywhere in
the table. So the effect Section 6.12 saw on three pairs under a *per-metal*
perturbation reappears on all twenty-one under a *scalar* one — genuinely
independent evidence that it belongs to the pair's sign structure — but the
magnitude comes down from 65× to 30× when the sample is widened from three
pairs to twenty-one. That is the third time in this chapter that widening a
sample has shrunk a claim.

The most sensitive pair is **Au/Pd** (`|S| = 2.56 eV⁻¹`), whose contacts are
0.02 eV apart; the least is Pd/Pt (0.0077). This was predicted in advance and
held.

> **Annotated 2026-09-25 — the numbers in this subsection were computed with a
> finite-difference step 3.5 decades too coarse; see Section 6.13.11.** The
> table's numbers are left as computed. Converged, the same-sign range is
> **0.5817 – 2.5662** rather than 0.6625 – 2.5628, and the dichotomy ratio is
> **26.7×** rather than 30.4×. The straddling column does not move at all
> (all six pairs agree to better than 1 part in 10⁵). The result of this
> subsection — a clean separation with no overlap, an order of magnitude wide
> — is unaffected, and so is the claim in this paragraph: Au/Pd remains the
> most sensitive pair at `δ = 0` converged, at `|S| = 2.5662`. What does not
> survive is the *prediction* P2, which was about `δ = +0.1 eV`, and it is the
> subject of 6.13.11.

<a id="s6-13-5"></a>

#### 6.13.5 Result: no scalar offset can reverse a response sign

Write `dW_A = m − s` and `dW_B = m + s`, with `s = (W_B − W_A)/2` and
`m = (W_A + W_B)/2 − w_cross`. A scalar offset moves `m` and leaves `s`
**exactly** fixed, so the contact field separates as

    E(x) = m [ f(L−x) − f(x) ] + s [ f(x) + f(L−x) ],    f(x) = λ⁻¹(1+x/λ)⁻²

into an antisymmetric part carrying `m`, which contributes exactly zero to
`N` by the same argument as Validation 3, and a symmetric part carrying `s`,
which the perturbation cannot touch. **`sign(N)` is fixed by which contact
has the higher work function.** A scan over `|δ| ≤ 3 eV` (601 points × 21
pairs = 12621 evaluations) finds zero sign changes; the smallest `|N|`
anywhere is 5.1e-03, approached and never crossed. `m` only rescales `|N|`,
vanishing as `|m| → ∞` — and that limit *is* the near-cancellation of Section
6.11.2, now with an explicit mechanism.

This is a **stronger robustness statement than any in Sections 6.8–6.12**:
the *direction* of a two-terminal photoresponse, which is what a design rule
actually asserts, is immune to the entire "~" on 5.4 eV, and indeed to a
scalar error six times its size.

It is also the opposite of what was predicted. The pre-registered derivation
D5 held that a sign flip occurs whenever `w_cross` moves between a pair's two
work functions, and tabulated 21 flip windows accordingly. The premise is
wrong: `sign(dW_A · dW_B)` is not what sets `sign(N)`.

<a id="s6-13-6"></a>

#### 6.13.6 A methodological finding: the five exact checks could not see the bug

The first implementation of the sign-flip scan **reported 21 roots.** It
bisected between `δ = 0` and the midpoint of D5's window without ever
checking that a root was bracketed; with no sign change present, such a loop
walks its lower bound up to its upper bound and returns the **endpoint**.
Every printed root equalled `(W_A + W_B)/2 − 5.4`, every one was physically
plausible, and its headline — *"nearest sign flip: Pd/Pt at −0.0150 eV"* —
is exactly the kind of number that reaches a thesis unchallenged.

All five exact validations passed while it did so, because they validate the
*model* and the fault was in the *analysis wrapped around it*. This is a
third distinct failure class, and the chapter should keep all three separate:

| class | example | caught by |
|---|---|---|
| implementation error in the model | the `inf − inf` transit integral, 2026-09-17 | an exactly known zero |
| a claim quantified on an unrepresentative sample | the `< 1.02` ceiling, Section 6.9 → 6.11 | enumerating all 21 pairs |
| **error in the analysis layer** | this section's unbracketed bisection | **neither** — hand-checking one row |

<a id="s6-13-7"></a>

#### 6.13.7 Result: Section 6.11.2's "five pairs" is retracted as a count

Section 6.11.2's table is re-scored at each offset (best perfect mask, per
*incident* photon, divided by `|N_uniform|`):

| pair | δ=−0.20 | δ=−0.10 | δ=0 | δ=+0.10 | δ=+0.20 |
|---|---|---|---|---|---|
| Au/Pd | 2.985 | 5.741 | **8.421** | 11.055 | 13.663 |
| Ti/Cr | 2.735 | 3.066 | 3.396 | 3.727 | 4.056 |
| Ni/Au | 1.534 | **2.442** | 3.357 | 4.277 | 5.195 |
| Cr/Cu | 2.496 | 2.868 | 3.242 | 3.614 | 3.989 |
| Ni/Pd | 1.122 | 1.801 | **2.486** | 3.176 | 3.869 |
| Ti/Cu | 1.337 | 1.513 | 1.689 | 1.865 | **2.041** |

The count above 2× runs **3, 4, 5, 5, 6** across the band. The qualitative
replacement claim — masks help whenever `|N_uniform|` is small — holds at
every offset, and so does the design recommendation, since Au/Pd's best
masked response is 0.377 per incident photon against unmasked Ti/Pt's 1.832.
**What does not survive is the number five**, which is an artefact of
evaluating a threshold at one point of a parameter the literature gives with
a tilde. Section 6.11.2 is annotated in place rather than rewritten.

The Section 6.9.2 ceiling `|N[g]| ≤ max|k|` held at every offset: 0
violations in 210 checks.

<a id="s6-13-8"></a>

#### 6.13.8 Result: Ti/Pt is now the most-tested claim in the chapter

Over the whole nominal band Ti/Pt's `|N|` ranges 1.8241–1.8367 — a **0.45%**
worst-case change, against 4.9% under Section 6.11's nonlinearity and
1.2–1.8% under Section 6.12's per-metal crossover. Three independent
refinements have now failed to move it. Whatever the risk to Chapter 6's
headline is, it is not the value of `w_cross`.

<a id="s6-13-9"></a>

#### 6.13.9 What this section does and does not change

**Changes:** Section 6.12.5's 65× amplification is restated as 30× and
reclassified as perturbation-dependent rather than a device property; Section
6.11.2's "five pairs" becomes "three to six across the band"; Sections
6.8–6.12 gain a stated scalar-crossover band, which for straddling pairs is
**under 0.5%** and for same-sign pairs reaches **51%** at `|δ| = 0.2 eV`.

**Does not change:** every doping-type assignment, every response sign, the
Ti/Pt headline, the `max|k|` ceiling identity, and the design recommendation
to choose metals before illumination.

<a id="s6-13-10"></a>

#### 6.13.10 What this section does not claim

1. It does **not** determine the crossover. It measures sensitivity to an
   unknown offset; the offset remains unknown.
2. A scalar offset is **not** the per-metal structure of Section 6.12, and
   the robustness result of 6.13.5 does not transfer to it. A scalar offset
   cannot change `s = (W_B − W_A)/2`; a *differential* one can. **A per-metal
   crossover difference of 0.02 eV would flip Au/Pd's response sign, where a
   3 eV scalar offset cannot.** That is now the chapter's sharpest open
   question, and Section 6.12 can answer it for exactly three metals.
3. `λ = 250 nm` is still one value for all metals, and the
   photo-thermoelectric term is still absent.
4. The band is a sensitivity band, not an error bar on a measurement: it
   propagates one parameter's uncertainty through a model whose other
   approximations are unchanged.

Figure: `crossover_sensitivity.png` — (a) sensitivity of all 21 pairs on a
log axis, straddling in blue and same-sign in red; (b) `N(δ)` over ±0.5 eV
with the nominal band shaded, showing curves approaching but never crossing
zero; (c) the Section 6.11.2 ">2×" claim crossing its own threshold twice
inside ±0.2 eV.

<a id="s6-13-11"></a>

#### 6.13.11 Annotation (2026-09-25): the step size, and the one prediction it overturned

Every logarithmic sensitivity in Section 6.13 was computed by a central
difference with the step `H_DIFF = 10⁻³ eV`. A numeric-default audit on
2026-09-25 measured, for the first time, where that step actually sits relative
to the plateau of its own derivative. The plateau — the range of `h` over which
`S(h)` is stationary in `h` — runs from about `10⁻¹⁰` to `3 × 10⁻⁷ eV`. The
shipped step is **3.5 decades above its top.**

The consequences are uneven and the unevenness is the whole story. Six pairs
move not at all: every straddling pair agrees to better than 1 part in 10⁵,
because `|S|` there is small and smooth. Four same-sign pairs move by under
0.4%. The remaining eleven move by 5.9% to **44.4%**:

| pair | `S` at `h = 10⁻³` | `S` converged | shift |
|---|---|---|---|
| Ti/Ni | −1.4037 | −0.9723 | 44.4% |
| Cr/Ni | −1.4790 | −2.4208 | 38.9% |
| Cr/Pd | −1.6942 | −1.2836 | 32.0% |
| Cr/Au | −1.6374 | −1.3223 | 23.8% |
| Cu/Au | −1.7729 | −2.2627 | 21.6% |
| Ti/Cu | −0.9503 | −0.7962 | 19.4% |
| Ti/Cr | −0.6625 | −0.5817 | 13.9% |
| Cu/Ni | −1.6256 | −1.4757 | 10.2% |
| Cu/Pd | −1.8114 | −1.6673 | 8.7% |
| Cr/Cu | −1.1687 | −1.0836 | 7.9% |
| Ni/Au | −2.3306 | −2.4768 | 5.9% |

Note that the shift is not even one-signed: Cr/Ni and Cu/Au were
*under*-estimated, the rest over-estimated. A uniform correction factor would
not have helped.

**What survives.** The dichotomy of 6.13.4 — the sole load-bearing result of
this section — survives cleanly, at 26.7× instead of 30.4×, with no overlap.
Au/Pd remains the most sensitive pair at `δ = 0`. Section 6.13.5's sign
robustness does not involve a derivative at all. Section 6.13.8's Ti/Pt band is
a ratio of `|N|` values, not a derivative, and is untouched. Nothing in
Sections 6.8–6.12 or 6.14 reads a number off `S`.

**What does not survive: prediction P2.** P2, pre-registered in the 2026-09-22
note and recorded there as **HELD**, asserted that the largest `|S|` at
`δ = +0.1 eV` belongs to **Au/Pd**, on the rationale that amplification is
near-cancellation and Au/Pd is the nearest-cancelling pair in the table
(0.02 eV apart). Converged, the ranking at `δ = +0.1 eV` is

| rank | pair | `\|S\|` (eV⁻¹) | contact separation |
|---|---|---|---|
| 1 | **Ti/Cu** | 2.828 | 0.32 eV |
| 2 | Cr/Au | 2.518 | 0.60 eV |
| 3 | Ni/Pd | 2.291 | 0.08 eV |
| … | Au/Pd | 2.016 | 0.02 eV |

and that ordering is identical at every step inside the plateau, across four
decades. **P2 is falsified, and it is falsified in exactly the way its own
rationale said it could be.** The 2026-09-22 note wrote: *"near-cancellation
could be non-monotone in separation once the kernel's spatial structure enters,
in which case some other pair wins."* Some other pair wins. The two closest
pairs in the table, Au/Pd at 0.02 eV and Ni/Pd at 0.08 eV, place fourth and
third; a pair 16× further apart places first. Near-cancellation is not
monotone in contact separation at finite `δ`, and the `δ = 0` extremum
(Au/Pd, which does hold) was giving a misleading impression that it was.

The mechanism of the artefact is worth one sentence, because it generalises
past this table. Au/Pd's `|S|` is nearly independent of `h` — 2.0162 across
five decades, 2.0609 at the shipped step, a 2% drift — while Ti/Cu's collapses
from 2.8282 to 0.8725, a 69% error. Au/Pd did not lead at `h = 10⁻³` because it
was the most sensitive pair; it led because it was the pair whose derivative
that step happened to get right, while its competitors' were destroyed. **A
ranking taken at a step outside the plateau ranks the entries by how well the
step suits them.** No amount of care about the physics would have caught that;
only a sweep over `h` does.

The shipped `H_DIFF` is deliberately *not* changed — these numbers stay
reproducible and comparable — and a second accessor,
`sensitivity_converged()`, supplies the plateau-anchored derivative with a
runtime check that refuses to return an unconverged one. Full measurement in
`default_scale_audit_output.txt` and `notes/2026-09-25-*.md`; the figure is
`default_scale_audit.png`.

## A.8 Chapter 7: Discussion and outlook

<a id="s7-8-1b"></a>

*Moved from Chapter 7, §7.8.1b. Original section number kept.*

#### 7.8.1b Annotation (2026-10-05): the saturation term was built, and it did not change the verdict — it sharpened the reason

Section 7.8.1a left the RF case *unaddressed* rather than refuted, and named
the one repair that could change that: a saturation term worth ×4.8483 in
g_ds. **That term has now been built and the verdict stands — but the reason
it stands is no longer "the model has no output resistance".**

What changed. With velocity saturation present and v_sat taken from
optical-phonon physics rather than fitted, the trend §7.8.1a called "the first
thing any such term has to fix" **is fixed**: f_max/f_T now rises with drain
bias (0.6632 → 1.3844 over V_ds = 0.05–1 V) where it previously fell. The
model is no longer a pure resistor, and the sign of its central RF diagnostic
is now right.

What did not change, and why it is a sharper statement than before. The same
term delivers only **1.1361** of the required 4.8483 — 23.4 % — and the
shortfall is not a weakness in the saturation law. At the literature-scale
geometry, R_c,total is 15.0 Ω of a ≈ 26 Ω device, so **roughly half of g_ds is
a contact resistance that no channel mechanism can reach.** Inverting the
question: criterion A would need v_sat = 6.50 × 10⁶ cm/s, *below* Dorgan
*et al.*'s measured 1–3 × 10⁷ cm/s band on SiO₂, and no phonon energy in the
physically available range (0.059–0.196 eV) gets past a factor of 1.2491.

**Three levers, three non-binding results.** The gate resistance is exhausted
(R_g = 0 exactly still falls 28 % short). The feedback capacitance is
negligible (C_gd deleted entirely moves the answer 0.11 %). Velocity
saturation is diluted by the contacts. Each was proposed in turn as the
missing piece of the f_max shortfall, each was built or bounded, and **all
three point at the access resistance.**

This converts §7.8.1a's verdict from a confession into a finding. The earlier
statement was that this thesis contributes nothing to the RF case because its
instrument measured the wrong quantity — a fact about the instrument. The
statement now available is a fact about the device: **in this model, at this
geometry, the RF figures of merit are contact-limited, and the channel physics
that the GFET literature treats as the interesting part is not where the
constraint sits.** That is a claim this thesis is entitled to make, because
it was reached by building the mechanism and measuring it rather than by
assuming it.

It also makes Section 7.8.2 load-bearing for an application it did not
previously claim. §7.8.1a already noted that R_s enters the f_max denominator
and carries 47 % of it, so the contact work was "not wasted" on RF. The
stronger version is now available: the contact work is the *only* lever left
standing for RF in this model, which means §7.8.2's negative result — that
work function alone does not describe Ti, Ni and Pd — is a limitation on the
RF case as much as on the contact case.

**What would overturn this annotation.** A single counterfactual, and it has
not been run: f_max/f_T at R_c = 0 *exactly*, in the saturated model, in the
way §4.6.2 ran R_g = 0 and R_s = 0. If a perfect contact also falls short of
1.3, then all four levers are exhausted and §7.8.1a's structural verdict
becomes final rather than provisional. If a perfect contact reaches the band,
then the RF case for graphene in this model is a contact-engineering problem
with a quantified target, which is a far more useful conclusion than either
of the two this chapter currently carries. **That counterfactual is the top
open item created by this annotation**, and it is deliberately not asserted
either way here.

See Chapter 4 §4.6.3 and
`notes/2026-10-05-saturation-arrived-and-the-contacts-ate-it.md`.

<a id="s7-8-1c"></a>

*Moved from Chapter 7, §7.8.1c. Original section number kept.*

#### 7.8.1c Annotation (2026-10-06): the counterfactual was run, and it answered both branches at once

§7.8.1b closed by naming one counterfactual as the thing that would overturn
it — `f_max/f_T` at `R_c = 0` *exactly*, in the saturated model — and stated
its two branches in advance. **It has now been run, and the answer is not
either branch. It is both, because the question was asked about a ratio.**

At the 40 µm / 8-finger geometry and `V_ds = 0.1 V`, with the contact removed
exactly and coherently:

| quantity | baseline | `R_c = 0` | change |
|---|---|---|---|
| `f_max/f_T` | 0.683186 | 0.759363 | +11.15 % — **58.4 % of the 1.3 floor** |
| `f_T` | 20.740 GHz | 74.650 GHz | **×3.5994** |
| `f_max` | 14.169 GHz | 56.687 GHz | **×4.0007** |

**On the ratio, §7.8.1b's first branch holds and the verdict is now final.** A
*perfect* contact — not a better one — leaves `f_max/f_T` 41.6 % below the
bottom of Feijoo *et al.*'s band. Four levers have now been built or bounded
and all four are non-binding: `R_g = 0` exactly gives 0.935276, deleting `C_gd`
entirely moves the answer 0.11 %, velocity saturation from optical-phonon
physics delivers ×1.1361 of a required ×4.8483, and `R_c = 0` exactly gives
0.759363. §7.8.1a's structural verdict on the *ratio* is no longer
provisional. A residual `g_ds` requirement of **×2.9955** survives the perfect
contact, and §7.8.1b already showed no phonon energy in the physically
available 0.059–0.196 eV reaches even ×1.2491 of it.

**On `f_max` itself, §7.8.1b's second branch also holds, and this is the part
an industrial reader acts on.** The same counterfactual is worth a factor of
**four** in `f_max`, and most of it is reachable: at the lowest contact
resistance reported in the literature (65 Ω·µm per contact, against this
model's 300) `f_max` rises ×2.79, which is **59.6 %** of the perfect-contact
gain. So the RF case in this model *is* a contact-engineering problem with a
quantified target — ×4.0 available, ×2.8 reachable — and that statement
coexists with the final verdict on the ratio rather than replacing it.

**Why both branches could be true, and it is not a drafting failure.**
`f_max/f_T` divides out exactly the quantity the contact dominates: the
contact limits `I_d`, `I_d` sets `g_m`, `g_m` sets `f_T`, and the same contact
sets `R_s` in the `f_max` denominator. Removing it moves numerator and
denominator together and the ratio keeps only the residue. **A criterion
written on a ratio is silent about any mechanism that scales both of its
arguments.** This is checkable rather than rhetorical: `f_max/f_T` is
**non-monotonic** in `R_c`, bottoming at 0.680047 at 470 Ω·µm and rising again
to 0.698615 at 4000 Ω·µm, while `f_max` falls **37×** across those same two
rows. There is a region of this design space where satisfying a ratio
criterion means losing on the figure of merit.

**One objection is now closed.** This model's contact is 300 Ω·µm per contact,
**1.8× the 165 Ω·µm** Feijoo *et al.* report for the very devices whose
1.3–1.4 band Chapter 4 is measured against. The whole `f_max` thread could
therefore have been nothing but a worse contact. Given Feijoo's own contact
resistance, this model reaches **0.695514**, not 1.3.

**Two findings consolidate, and one claim is narrowed.** §7.8.1b's dilution
*explanation* — that the contacts absorb the saturation — made an untested
prediction, that the `g_ds` saturation factor at `R_c = 0` must be the
undiluted 1.699. Measured: **1.70661**, i.e. 100.45 %, from a derivative of a
solved terminal current rather than from the channel integral that produced
the 1.699. And §7.8.1a's resistance-ratio identity acquires a sign test: the
*superseded* resistor model is made **worse** by the same perfect contact
(0.642901 → 0.581707), because that identity carries `R_total` in the
numerator and therefore rewards a worse contact. A figure of merit that
improves when the device gets worse is not measuring the device, which is what
§7.8.1a asserted and now has evidence for.

**What this annotation does NOT claim, and the reason is a literature result
that cuts against §7.8.1b's framing.** Feijoo, Pasadas, Bonmann *et al.*,
*Nanoscale Advances* **2** (2020) — a paper whose title asks "Does carrier
velocity saturation help to enhance `f_max` in graphene field-effect
transistors?" — answer their own question in the negative: "the largest
`f_max` are located at biases close to the onset of bipolar conduction and
**far from the saturated velocity regime**", with drift at only ~45 % of
`v_sat` and the **diffusion** contribution to the current comparable to drift
at the peak-`f_max` bias, concluding that "our results do not support that
operating in the regime of velocity saturation results in the highest
`f_max`". Three things follow, and none of them is comfortable:

1. **§7.8.1b's criterion A was constructed from this model's algebra, not from
   device physics.** Its arithmetic stands — ×4.8483 required, ×1.1361
   delivered — but it licenses a statement about *this model* and not about
   real GFETs, and that narrowing is recorded here rather than left implied.
2. **This model has no diffusion current at all.** Eq. (4) of the saturated
   model is drift-only. If diffusion is ~40 % of the current near the drain at
   the peak-`f_max` bias in a measured device, a drift-only model cannot be
   asked about that bias, and the `V_ds` ladders of §7.8.1b and of this
   annotation are a drift-only slice of a two-mechanism problem. **This is now
   the top open item for the RF thread**, ahead of anything else in it.
3. **This model has no self-heating**, a mechanism the same paper measures at
   65 → 40 GHz in `f_max`. Today's ×4.0007 is a ratio between two runs of one
   model and is largely insulated from this; the absolute 56.687 GHz is not.

   > **ANNOTATION 2026-11-10 — this worry is now MEASURED, and it is
   > NARROWED rather than confirmed. No number above is changed.**
   >
   > `graphene_self_heating_model.py` builds the missing term from the
   > thermal network that Dorgan, Bae and Pop (*APL* 97, 082112 (2010))
   > carry alongside the mobility and saturation-velocity fits this chapter
   > already uses, so `R_th` is checkable rather than asserted. The
   > implementation reproduces that paper's own term shares
   > (3.60 / 84.49 / 11.91 % against its stated 4 / 84 / 12 %, total
   > 9915.6 K/W against its ~1e4 K/W) and its footnote-18 branch-cut
   > densities (9.4043e11 and 1.1113e11 cm⁻² against 9.4e11 and 1.1e11).
   >
   > **At the bias this chapter's `f_max` is computed at — `V_ds` = 50 mV —
   > the self-consistent temperature rise is 1.0675 K and the current loss
   > is 1.8 × 10⁻⁶.** The absolute 56.687 GHz is therefore exposed *in
   > principle* and negligibly exposed *in fact*, at this bias. The clause
   > above — "the absolute 56.687 GHz is not" — stands as written, because
   > it is true of the model's silence; what is now known is the size of
   > what the silence was hiding, and it is small here.
   >
   > **Where it is not small is the bias a real RF GFET is driven at.** At
   > `V_ds` = 1 V the same network gives Δ*T* = 168.7 K, and the sweep
   > passes Dorgan's own ~200 K fit-exclusion ceiling at `V_ds` = 1.25 V —
   > i.e. this device's thermal environment puts it outside the regime its
   > own transport parameters were measured in before it reaches the bias
   > where the literature's peak `f_max` is found. That, and not a
   > correction factor, is the usable result.
   >
   > Three caveats that keep this an upper bound rather than an estimate.
   > (i) Eq. (5) assumes `W, L >> t_ox`; this device has `L/t_ox` = 2.22,
   > so lateral spreading in the oxide is neglected and `R_ox` is
   > **over**-estimated — the error is directional, so every Δ*T* here is
   > a worst case. (ii) The boundary conductance `h` = 1e8 W m⁻² K⁻¹ is a
   > **secondary citation** (Chen *et al.* 2009 *via* Dorgan). (iii) The
   > 65 → 40 GHz figure that created this item is Feijoo *et al.*'s
   > measurement on *their* device, and nothing here reproduces it; what is
   > reproduced is the thermal *model*, on the device this thesis actually
   > models.
   >
   > **The 160 meV branch could not be used at all**, which is a result and
   > not an omission: at `V_g` = 3 V this device's density is 7.26e15 m⁻²
   > against a branch cut of 9.40e15 m⁻², so Dorgan's two phonon energies
   > **cannot bracket** `v_sat(T)` here. `V_g` = 4.1277 V would be needed.
   > The remaining admissible spread is 113.97 K (55 meV) to 168.72 K
   > (82 meV) at `V_ds` = 1 V, so the phonon energy — not the mobility-fit
   > ambiguity, which moves Δ*T* by 0.004 K — is the dominant uncertainty.

So the RF verdict of §7.8.1a — *unaddressed*, not refuted — survives in a
third form. The ratio is settled and contact-limited only mildly; `f_max`
itself is contact-limited by a factor of four; and the reason this thesis
still cannot speak to the RF case is no longer the absence of saturation, but
the absence of diffusion and of self-heating in the model that was built to
replace the resistor.

**AMENDED 2026-11-10.** Self-heating is no longer absent; diffusion still is.
The sentence above is kept as written because it records what the chapter
could say on 2026-10-06, and the correction is that the RF thread's
remaining blocker is now **diffusion alone** — see the annotation to item 2
above, which is unchanged, and the annotation to item 3, which is the work
that removed it from the list. The self-heating half of the verdict is
discharged in the direction of *the model was conservative here*, which was
not the expected direction.

See `graphene_perfect_contact_counterfactual.py`,
`perfect_contact_counterfactual_output.txt` and
`notes/2026-10-06-the-perfect-contact-and-the-ratio-that-divided-out-the-prize.md`.

<a id="s7-8-1d"></a>

*Moved from Chapter 7, §7.8.1d. Original section number kept.*

#### 7.8.1d Annotation (2026-10-07): the diffusion term was built, and the item's premise did not survive it

§7.9 item 7, created 2026-10-06 and named that session's top item, said that
Eq. (4) is drift-only and that every `V_ds` ladder in §7.8.1b and §7.8.1c is
therefore a drift-only slice. The term has been derived and built
(§4.6.5, `graphene_diffusion_current_model.py`, 30 checks, 30 passed). Three
parts of the item are now settled and they do not settle the same way.

**It is worth ±0.50 % of `I_d`, and the sign flips between branches.**
+0.5018 % at the saturated-peak bias and −0.5038 % at the resistor-peak bias.
`f_max` gains +0.2941 % while `f_max/f_T` *falls* 0.2514 %, from 0.683186 to
0.681468. §7.8.1a's structural verdict is untouched — 52.42 % of 1.3 either
way — and the ratio moving against `f_max` is §7.8.1c's finding in a third
mechanism.

**The premise is the part that moved.** Pasadas and Jiménez
(*IEEE TED* 2016, [arXiv:1605.08235](https://arxiv.org/pdf/1605.08235)), from
the same group as the Feijoo *et al.* paper that created this item, write the
GFET channel current as `I = −W Q_tot v` with `v = μF`, `F = −dV/dx`, and
state explicitly that **"V(x) is the quasi-Fermi level along the graphene
channel"**. A drift expression whose driving potential is the quasi-Fermi
level is *already* the complete drift-diffusion current. Under that reading,
the term added in §4.6.5 is a double count, not missing physics.

`carrier_density(V_g, V_ch)` contains evidence for both readings and settles
neither. `V_ch` enters an electrostatic charge relation
`C_ox(V_g − V_ch − V_dirac)`, which reads electrostatic; it is then corrected
by a quantum-capacitance series factor, which is exactly what one applies when
`V_ch` is the quasi-Fermi level and the graphene drop `E_F/e` is a separate
voltage. **The two readings differ by exactly `λV_F`** — the term itself. So
±0.50 % is simultaneously the diffusion current under one reading and the size
of an ambiguity that was already present in every `I_d` this thesis has
committed.

This is §7.5's recurring class in a new location. 2026-10-04 recorded that a
quantity can be computed correctly under a name that claims more than it
measures; here the name is a **variable**, and what it fails to pin down is
which of two potentials it is. No value is wrong. The equation is not
determined.

**What this changes in this chapter.** The claim that §7.8.1b's and §7.8.1c's
ladders are *drift-only* slices is withdrawn: they are slices of a model whose
driving potential was never named, bounded at ±0.50 %. The repair is not the
one item 7 asked for — **it is to name the potential, not to add a term** —
and the recommendation is to declare `V_ch` the quasi-Fermi potential, keep
`λ = 0` as the production path, and keep §4.6.5 as the bound on the choice.

**And item 7's other claim is upheld for a different reason.** The
peak-`f_max` bias remains a bias this thesis cannot be asked about. Not for
the missing term: §4.6.5 shows the diffusion-dominated region near charge
neutrality is *absent* from this model, because the `n_puddle` regularisation
makes `d|n|/dV` exactly zero where Zebrev's `C_ox/C_Q` diverges. A
regularisation chosen to keep the conductivity minimum finite removed the one
regime the literature says matters, and no amount of care about the transport
term recovers it.

See §4.6.5 and
`notes/2026-10-07-the-diffusion-term-was-a-boundary-term-and-the-potential-was-never-named.md`.

<a id="s7-8-1e"></a>

*Moved from Chapter 7, §7.8.1e. Original section number kept.*

#### 7.8.1e Annotation (2026-10-08): the potential is named, and the naming question was the smaller half

§7.8.1d recommended declaring `V_ch` the quasi-Fermi potential, keeping
`λ = 0` as the production path, and keeping §4.6.5 as the bound on the choice.
**All three have been done** (§4.6.6). The declaration lives in
`graphene_fet_model.py` as a machine-readable in-source marker, so the
statement of what this thesis's own variable means travels with the variable.

**What §7.8.1d got right.** The term is a double count, not missing physics:
+0.6965 % intrinsic and **+0.4279 % at the terminals** at `V_g` = 2.0 V,
`V_ds` = 0.1 V, bracketing §4.6.5's +0.5018 % at its own bias. The claim
withdrawn in §7.8.1d — that §7.8.1b's and §7.8.1c's `V_ds` ladders are
*drift-only* slices — can now be restated positively: **they are slices of the
complete drift-diffusion current in the quasi-Fermi variable**, with no
±0.50 % ambiguity attached. That ambiguity is resolved rather than bounded, and
§7.8.1a's structural verdict is untouched for the third consecutive session.

**What §7.8.1d did not see.** A declaration is a promise that the rest of the
model belongs to the declared variable, and that promise has a price. Under the
declaration the exact charge relation is
`V_g − V_dirac − V_ch = e·n/C_ox + E_F(n)/e`; §4.3 computes its linearisation.
Running Eq. (4) on the exact form moves `I_d` at the RF bias by **−0.6023 %**
— **1.4 times the term §7.9 items 7 and 8 were competing over.**

This is the uncomfortable part, and it belongs in this chapter rather than only
in Chapter 4. §7.5.6 identified the recurring class and §7.7 identified unequal
scrutiny as its cause. Here is a clean instance of the latter: **two sessions
(2026-10-06 and 2026-10-07) went into a transport term worth +0.43 %, while a
−0.60 % question about the same chapter's charge relation sat directly
underneath it, reachable by replacing one function and re-running.** Nothing
protected the charge model from scrutiny except that no item pointed at it, and
no item pointed at it because the attention had been captured by a literature
paper's framing (Feijoo *et al.*'s diffusion measurement) rather than by an
audit of this thesis's own equations. §7.7's thesis survives contact with a
third example.

**What it does not change.** Every number in Chapters 4 and 5 remains a
linearised-series-factor number, annotated in place and not deleted, per the
standing rule. The −0.60 % is a *measurement of a gap*, not a new result set,
and the two charge forms have been compared only to each other — never to
measurement — so nothing here licenses calling the committed model wrong, only
non-exact under its own declaration. §7.8.5's summary judgement is unchanged.

**And one ranking change that follows.** §7.9 item 9 (a signed-carrier FET
model, created 2026-10-07) is promoted from a sibling of item 7 to the
**prerequisite** of the new top item 12. Both charge forms carry `n` as a
magnitude, so re-deriving Chapter 4 on the exact relation before a
signed-carrier treatment exists would rebuild the branch-sign flaw in a new
equation rather than fix it. The hole branch is now blocking, not adjacent.

See §4.6.6 and
`notes/2026-10-08-the-declaration-and-the-question-underneath-it.md`.

<a id="s7-8-1f"></a>

*Moved from Chapter 7, §7.8.1f. Original section number kept.*

#### 7.8.1f Annotation (2026-10-09): the prerequisite is cleared, and the off-state number this thread assumed was 2.07x off

§7.8.1e ended by promoting §7.9 item 9 to the prerequisite of the RF thread's
top item. **Item 9 is closed** (§4.6.7), item 12 is unblocked, and for the first
time in three sessions the top item has nothing in front of it. That is the
whole of this annotation's RF content: no `f_T` or `f_max` number moved today,
and the saturated-model queue — item 12, then the §4.6.5 rewiring, then
self-heating — is unchanged in order and one position shorter at the front.

**The non-RF finding is the one that bears on §7.8.1 rather than §7.8.1a.** This
thread's verdict that graphene is not a logic transistor rests on the off-state
floor, and §4.6.7 finds that Chapter 4's floor and Chapter 3's measured floor
disagree by **2.0678×**, with Chapter 4 on the more conductive side. Correcting
to Chapter 3's floor raises the on/off ratio from **1.2408 to 1.7820** at
`V_ds` = 0.1 V. **The verdict is undisturbed and this is worth saying plainly
rather than quietly**: an on/off ratio of 1.78 is no more a switch than one of
1.24, and §7.8.1's reason — that the floor is physics and not a process
limitation — is if anything strengthened by the discovery that the floor was
being *understated*. A finding that moves a number by a factor of two without
moving the conclusion it supports is the useful kind of robustness result, and
§7.6 is where it belongs.

But it is also an instance of §7.7. The floor comparison was a one-line
calculation flagged in §4.9 on 2026-09-29 and left for **eleven sessions**,
while two consecutive sessions went into a ±0.5 % transport term in the RF
thread. §7.7's unequal-scrutiny finding now has a fourth instance, and the
pattern across all four is the same: **the cheap check on the quantity the
conclusion actually rests on loses to the expensive check on the quantity the
current thread is about.**

<a id="s7-8-1g"></a>

*Moved from Chapter 7, §7.8.1g. Original section number kept.*

#### 7.8.1g Annotation (2026-10-10): the fifth lever, and it was two levers that cancelled

§7.9 item 12 — re-deriving Chapter 4 on the exact charge relation — is closed
(§4.6.8). It is the **fifth** lever this thread has tried against the RF
verdict, after velocity saturation (10-05), the perfect contact (10-06), the
diffusion term (10-07/10-08) and the charge relation's prerequisite (10-09).
The verdict is unchanged and the *reason* it is unchanged has changed shape
twice.

**What moved.** `f_max/f_T` moves **+0.0228 %**, from 0.643007 to 0.643154 at
the 40 µm 8-finger geometry. Peak `f_T` moves −0.0594 %, peak `f_max`
−0.0366 %, extrinsic `f_T` −1.6720 %, extrinsic `f_max` −1.7091 %. The
structural statement of §7.8.1a — that this model's `f_max` is a restatement
of `R_total/(R_g+R_s)` and not an estimate of `f_max` — is therefore
untouched, which is the correct thing to have checked. An earlier draft of
that check tested a *frequency* instead (that `f_T` sat "two orders below the
100–300 GHz RF-frontend requirement") and was wrong twice: the gap is about
one order, and §7.8.1a never rested on a frequency threshold. §7.5's taxonomy
gains a small class (v) instance in the audit layer rather than in the model.

**What the small number conceals, and this is the new thing.** The −0.0594 %
in `f_T` is **a cancellation, not a smallness**. The charge change alone is
−1.9986 %; the gate-capacitance change that the same differential identity
forces is +1.9392 %; the two are within 3 % of each other in magnitude and
opposite in sign. Four of the five levers in this thread were reported as
single numbers. This one could not have been: reporting the charge alone would
have overstated the effect by 34×, and reporting only the total would have
called a pair of 2 % corrections negligible. **§7.4's distinction between a
bound and a ranking has a third member — a number that is small because two
larger numbers oppose — and nothing in this thesis's audit record would detect
it, because every check that has ever been written here asks about one
quantity at a time.**

**And the peak bias is where the sensitivity lives.** Peak `f_T` moves
0.0594 % in height and **6.76 %** in location, a ratio of 114. Both peaks stay
on the hole branch, so §4.6.7's open "does the peak-`f_T` bias change branch
between the two models" question is answered **NO** for the
(4.28)-vs-(4.29) pair, and stays open for the branched-vs-magnitude pair it
was actually asked about. The observation that this thread's figures of merit
are far more sensitive in *where* they peak than in *how high* they peak is new
and is not something any of the five levers was designed to reveal.

**What this does to the ledger.** The RF thread has now exhausted its
model-side levers; what remains (item 16, the `C_q` argument) is the first
open item with a cost larger than any lever tried, and it is a defect rather
than a missing physical effect. §7.8.1a's conclusion is unchanged and better
supported than it was, because the levers that could have moved it have been
tried and measured rather than assumed to be small.

<a id="s7-0-note"></a>

*Drafting note moved from the head of Chapter 7.*

> **Drafting note (2026-09-27).** This chapter's material has been complete for
> some weeks; it was planned and displaced on seven consecutive working
> sessions, each time by a numerical result that was worth the displacement.
> The ten synthesis threads it was handed are listed in Chapter 1's status
> table. This draft discharges them, and adds two that arose on the day it was
> finally written.

<a id="s7-8-2-superseded"></a>

*Superseded 2026-11-16 by the crossover correction ([A.1](#corrections)). The §7.8.2 bullet as written; its result was an artefact of measuring the contact doping from 4.5 eV.*

> * **The additive decomposition `R_c = R_transmission + R_extra` fails for three
>   of four metals, and it fails robustly.** The negative residuals are not
>   amplified input error: Cu's would need an 81 % error in `R_extra` to flip
>   sign (Section 4.10.3). This is a real modelling result and it says the
>   two-term picture is incomplete, not mis-parameterised.

<a id="s7-4-3-superseded"></a>

*Superseded 2026-11-16, same cause. The §7.4.3 bullet as written:*

> * **Chapter 4 §4.7's additive decomposition.** Its failure is the cleanest
>   candidate for a bound: see §7.5.2 below, where the conditioning result turns a
>   negative result into a robust one.

<a id="s7-9"></a>

*Moved from Chapter 7, §7.9. Original section and item numbers kept.*

### 7.9 What would change these conclusions

Stated as falsifiable items, in order of how much they would move:

1. **A second anchor for `Δ_c` at any separation other than 3.3 Å.** Every
   "reachable" number in Section 6.12's margin table is one anchored
   exponential with a swept decay length. A second anchor would make Section
   6.12's per-metal crossover a measurement instead of an interpolation, and it
   is the top *physics* item in this thesis — open since 2026-09-21 and
   untouched for seven consecutive working sessions while the numerical thread
   ran.
   **ACTED ON 2026-10-15, after twenty-four consecutive sessions untouched.
   The item NARROWS and does not close, and it narrows in a way the item did
   not anticipate: the second anchor is in the source, stated in words, as an
   INEQUALITY rather than a number.** Khomyakov et al. introduce the
   short-range term `Δ_c(d) = e^(−γd)(a₀ + a₁d + a₂d²)` as their Eq. 4 and
   **never evaluate it numerically in the text** — no `γ`, no `aᵢ`, no `Δ_c`
   tabulated at any separation — so twenty-four sessions of looking for a
   second `(d, Δ_c)` pair were looking for something the paper does not
   contain. What it does contain is a **second separation**: its
   charge-transfer calculation is run at `d = 5.0 Å` as well as at
   `d_eq ≈ 3.3 Å`, and of the far one it says the crossover "decreases to its
   expected value, `W_M ∼ W_G`, for large separations". Through
   `W_0(d) = W_G + Δ_c(d)` that is exactly `Δ_c(5.0 Å) ≤ ε`, which inverts to
   a bound on the swept parameter of (N1):
   `ℓ ≤ (5.0 Å − 3.3 Å)/ln(0.9 eV/ε)`. At the only small energy scale the
   paper attaches to "small" nearby — its 0.08 eV — this gives
   **`ℓ ≤ 0.7024 Å`, excluding 66.5 % of the committed sweep range
   [0.3, 1.5] Å**;
   *[Annotated 2026-10-18. Both numbers stand and neither is revised. The
   hedge above — "the only small energy scale the paper attaches to 'small'
   nearby" — understates the borrowing in two ways the source's verbatim
   sentence makes plain: the `∼0.08 eV` is the deviation at `d ≈ 3.3 Å`, the
   separation about which that same sentence says Eq. (1) **holds**, and the
   paper traces it to a perturbation "that cannot be described as a rigid
   shift" — while `Δ_c` **is** a rigid shift. So it is the wrong quantity at
   the wrong separation, and the two defects compound rather than cancel.
   **There is therefore no privileged row:** the family over
   `ε ∈ [0.02, 0.20] eV` is the result, and no single `ε` is "the paper's".
   The sentence is kept as written because the bound it reports is correct
   as arithmetic; what was wrong was the attribution, and the transcript —
   which until today carried the attribution and none of the caveat — now
   prints its own provenance. See §4.6.16's sibling finding: a qualification
   with no output behind it is 2026-10-01's dead-name class inverted.]* the top of that range would put the 5.0 Å crossover
   **0.2898 eV above `W_G`**, which is not "≈ `W_G`" under any reading. The
   committed default `ℓ = 0.6 Å` **survives**, and that was not guaranteed:
   it was chosen as "the midpoint of the swept range", a choice with no
   physics in it, and it lands inside the surviving window by luck. **No
   doping verdict moves** — the sign of `W_M − W_0` is identical at the bound
   and at the old top of the range for all six metals with a tabulated
   `d_eq` — so this is a **tightening of the error bar on `W_0`, not a
   correction**, and Section 6.12's margin table is unchanged. What stays
   open is a **two-sided** constraint: nothing here bounds `ℓ` from below,
   and the cheapest remaining route is digitising the paper's Fig. 4 upper
   point set, which the fetch could not retrieve (it truncated three times at
   the opening of Sec. III.3). Six exact validations, 6/6, including a
   planted-`ℓ` round trip at 2 ulp and a positive control on the exclusion
   counter. See `notes/2026-10-15-the-second-anchor-exists-and-it-is-an-inequality.md`,
   `graphene_second_anchor_bound.py`.
2. **A description of Ti, Ni and Pd that does not go through work function.**
   Three independent failures are on record for these chemisorbed metals: the
   nonlinear doping relation excludes them (`d_eq < d_0`), the anchored
   exponential diverges on them, and they carry the four smallest reversal
   margins. Ti has the largest `|ΔW|` in the table (−1.07 eV) and is contact A
   of *both* Chapter 6 headline pairs. The chapter's central numbers rest on the
   one assumption the literature most explicitly disowns.
   *[Annotated 2026-10-20. A **fourth** independent failure is now on record,
   it is upstream of the three above, and — unlike them — it is **not confined
   to the chemisorbed metals.** No number in this item is revised.
   `graphene_chord_decay_length.py` costs the question of whether `ℓ` is the
   right parameter at all. Writing Eq. (4) as `e^(−γd) P(d)` and `L = ln P`,
   the chord and local decay lengths are `1/ℓ_chord = γ − ⟨L'⟩` and
   `1/ℓ_local(d) = γ − L'(d)`, so their difference is `L'(d) − ⟨L'⟩` and **`γ`
   cancels exactly** — as does `a₀`, since `L' = P'/P` is scale-invariant.
   Two of the four coefficients the paper never reports leave the question,
   and the gap is computable from `P`'s two shape ratios alone.
   By the mean value theorem `ℓ_chord` equals `ℓ_local(ξ)` at one interior
   separation `ξ`. Measured over 204 admissible quadratic shapes,
   **`ξ ∈ [3.8482, 4.4952] Å`, above Al's 3.41 Å — the largest tabulated
   `d_eq` — by `+0.4382 Å`.** So the swept parameter is the local decay length
   at **no metal in the table, under no admissible shape**, and the four
   physisorbed metals that sit *inside* the chord interval are not spared:
   the gap vanishes only at `ξ`, not across the interval. At `ℓ = 0.6 Å` the
   ratio `ℓ_local/ℓ` runs 0.89–1.95 for the physisorbed metals and 0.79–2.73
   for Ni and Co; at the top of the committed sweep, `ℓ = 1.5 Å`, **92
   metal-shape combinations have no representable local length at all** —
   `ℓ_local` infinite or negative, which no single exponential reproduces at
   any `ℓ` — and those reach the physisorbed metals too.
   This **compounds with** the `d_eq < d_0` exclusion rather than replacing
   it: Ti, Ni, Co and Pd now fail for two independent reasons, a regime
   boundary and a parameter read off an interval they are not in. It does
   **not** move item 1's bound, which is arithmetic about the chord and
   remains correct as arithmetic; what it says is that the bounded quantity
   is not the quantity the per-metal evaluation uses. See
   `notes/2026-10-20-gamma-cancels-and-the-parameter-is-right-at-no-metal.md`.]*
3. **Reconciling the 0.12 eV measured step with the 0.25–1.07 eV assumed
   offsets** (§7.8.2). Bounded in the common channel by Section 6.13's
   sensitivity table; **not** bounded in the differential channel, which is the
   one that sets polarity.
4. **A photo-thermoelectric term.** The parity of Section 6.14 is derived for a
   drift-only model with `μ_e = μ_h`. Whether it survives a
   photo-thermoelectric contribution is the single most consequential open
   question about the thesis's one symmetry-derived result — because if it does
   not, §7.2.3's polarity statement loses its exactness.
5. **Asking the differential question of Chapter 5** (§7.3's prediction).
6. ~~**Re-running Chapter 4's contact-resistance results at the 5.4 eV
   crossover**~~ rather than at graphene's 4.5 eV, for consistency with Chapter 6.
   — **CLOSED 2026-11-16** (§4.7, [A.1](#corrections)). The contact model now measures
   the doping from `W_CROSSOVER = 5.4 eV`; `tests/test_contact_crossover.py` pins
   n_contact(W₀) = 0 exactly. The §4.7 negative residuals did not survive, and the
   section's conclusion is withdrawn.
7. ~~**Naming the potential in Eq. (4)**~~ — **CLOSED 2026-10-08** (§4.6.6,
   §7.8.1e). `V_ch` is declared the **quasi-Fermi potential**, in the source
   where it is defined. The consequence is that 2026-10-07's diffusion term is
   a **double count** worth +0.4279 % of `I_d` at the terminals, not a missing
   0.50 %, so the committed `λ = 0` currents stand and `λ` is retained as the
   instrument that measures the difference between readings. The item's
   history — created 2026-10-06 as *a diffusion term in the saturated
   transport model*, rewritten 2026-10-07 when the term turned out to be the
   difference between two readings of a variable — is kept because the shape of
   it recurs. **It is replaced, not succeeded, by item 12.**
8. ~~**A self-heating term.**~~ The same paper measures 65 → 40 GHz in `f_max`
   from self-heating alone. This bounds every absolute frequency in Chapter 4;
   the *ratios* between model runs, which is what §7.8.1c reports, are much
   less exposed. — **CLOSED 2026-11-10**,
   `graphene_self_heating_model.py`, 41/41 checks, 10/10 mutants detected,
   both of Dorgan's published benchmarks reproduced. The item stood for 34
   sessions and the obstacle was never the physics: it was that a thermal
   resistance with no provenance would have been a number that looked like a
   measurement. It became answerable when the primary source for this
   chapter's *transport* parameters turned out to carry its own *thermal*
   budget, with the geometry and the per-term shares, making `R_th`
   checkable. **The expected conclusion was wrong in a useful direction**:
   the bound on Chapter 4's absolute frequencies is 1.8e-06 in current at
   the bias those frequencies are computed at, not the large correction the
   65 → 40 GHz figure suggested — because that figure belongs to a device
   driven at volts, and this one is characterised at 50 mV. What the model
   does find is a *regime boundary* at `V_ds` = 1.25 V, past Dorgan's own
   fit-exclusion ceiling. Two of the session's own errors were caught by
   the suite before commit (a tolerance that tested a precision the source
   never claimed, and a monotonicity claim that was simply false — the
   isothermal current turns over too, at 1.11 V), and both are recorded in
   place rather than quietly fixed.
9. ~~**A signed-carrier FET model**~~ (added 2026-10-07) — **CLOSED
   2026-10-09**, `graphene_fet_signed_carrier_model.py`, 19/19 checks, §4.6.7.
   The sharper statement of why it blocked item 12 than 10-08 gave: both terms
   on the right of (4.28) are **odd** in `n`, so the relation has an odd root,
   and an odd root cannot come out of a solver whose first line is
   `dV = np.abs(dV)`. `carrier_density()` *computes* the sign — `n_eff` is odd
   in `dV` because `C_q` is even — and discards it in the quadrature floor.
   `n_net_exact()` carries it: odd **bitwise**, satisfying (4.28) to 7.690e-16
   relative. The decomposition `n_e,h = (S ± n_net)/2` reproduces the shipped
   magnitude as an identity, yields `n_e·n_h = n_puddle²/4` (so the August
   quadrature floor is algebraically a mass-action law), and makes graphene's
   ambipolar symmetry a bitwise check. **No number in Chapter 4 moved** — 3 ulp
   on `σ`. The two sub-questions this item carried are now separated: the
   hole-branch sign of Chapter 4's branch-asymmetric results is *expressible*
   but still unverified, and the peak-`f_T` branch-crossing question is
   untouched; both are retained in the backlog under their own entries rather
   than inside a closed item.
10. ~~**A census of every potential-like variable in this thesis**~~ —
   **CLOSED 2026-10-08** as a standing instrument,
   `graphene_potential_census_audit.py`. 86 names across 25 modules, each
   carrying an explicit adjudication; a name present in the repository and
   absent from the registry is reported as a fault, so the census cannot go
   stale silently — and it caught its own staleness on the very next file
   added. **One defect in 86 names.** All three siblings this item predicted
   were checked and none is one, for three different reasons: Chapter 5's
   interconnect model contains *no potential-like variable at all*; Chapter 6's
   two-contact `V_bias` is imported and printed and enters no computation in
   its own module; §4.4's `V_channel_shift` appears only in a docstring
   derivation the evaluated closed form never uses. The first of those is kept
   as a check that is *expected to fail*, because a pass there would mean a
   potential had appeared in Chapter 5. Seven names remain `UNADJUDICATED` and
   are recorded as open work; none is channel-swept, so none can carry the
   defect. **A predicted class is not a measured class**, and the measured rate
   is 1 in 86.
11. **Whether the `n_puddle` regularisation removes other regimes, not only
   the diffusion-dominated one** (added 2026-10-07). §4.6.5 found that
   `d|n|/dV` is *exactly* zero at `n_eff = 0`, which deletes the
   charge-neutrality physics that Zebrev's `C_ox/C_Q` and Feijoo *et al.*'s
   peak-`f_max` bias both live in. A floor added to keep one quantity finite
   silenced a different one. Which other near-Dirac results in Chapters 4
   and 6 rest on a derivative the floor has set to zero is not known, and this
   is a sibling of the §4.9 item on whether `n_puddle` reproduces §3.3's
   6.45 kΩ/sq floor — open since 2026-09-29 and untouched for ten sessions.
   **2026-10-08 adds a reason to prefer this item over item 13:** §4.6.6's
   exact charge relation (4.28) and the series factor (4.29) disagree by +68 %
   raw at `dV` = 0.01 V and by 0.001 % once the floor is applied, so the floor
   is now known to be the dominant term in the near-Dirac charge of *both*
   forms. Whatever it has silenced, it has silenced in the replacement model
   too.
   **2026-10-09 — PARTIALLY ANSWERED, and the premise of the question was
   wrong.** The question asked which *other near-Dirac regimes* the floor had
   removed. §4.6.7 measures that the floor is not a near-Dirac term at all: at
   `V_g` = 3.5 V, the largest overdrive swept anywhere in this thesis, it is
   **0.774** of the gate-induced net density and adds **+26.4 %** to the total,
   and its contribution falls below 1 % only at `dV` = **14.73 V**, 5.46× the
   maximum swept. It does not remove regimes; it adds a parallel conducting
   density to all of them, and `carrier_density()`'s docstring calling it a
   Dirac-point regularizer describes a regime this thesis never enters. It was
   found by a check written to PASS on that docstring's authority, which failed
   at −6.92 % in `I_on` — an allegedly off-state-only parameter moving the on
   state. **What remains open is the original `d|n|/dV` = 0 question**, which
   is about the derivative and is not answered by this. The sibling §4.9 floor
   item is **CLOSED** (see item 11a).
11a. ~~**Whether Chapter 4's `n_puddle` reproduces §3.3's 6.45 kΩ/sq floor**~~
   (§4.9, added 2026-09-29, **untouched for eleven sessions**) — **CLOSED
   2026-10-09**, §4.6.7 and §3.3.1. It does not. `σ_min = e·μ·n_puddle` =
   3.204353e-04 S/sq = **3.1208 kΩ/sq** against §3.3's measured `4q_e²/h` =
   **6.4532 kΩ/sq**: Chapter 4's channel at neutrality is **2.0678× more
   conductive** than the floor Chapter 3 calls measured, and the reconciling
   `n_puddle` is 2.418e11 cm⁻² against the shipped 5e11 cm⁻². It stayed open
   eleven sessions not because the calculation is hard — it is one line, as the
   item itself always said — but because `n_puddle` is a floor on a
   **magnitude**, and a magnitude cannot say whether it floors the net density,
   the total, or each species: three readings differing by factors of 1, 1 and 2
   in `σ_min`, the same order as the answer. Item 9 settles it as the total, and
   the per-species reading moves the discrepancy the **wrong way** (4.1357×), so
   the convention is not the explanation. Cost: §4.4's on/off ratio is
   understated by 1.7357× intrinsically and 1.4362× at the terminals, the
   contacts attenuating the error by a further 1.2085× — §4.6.6's
   negative-feedback structure on a different parameter, hence a property of the
   device topology. §7.2 is undisturbed: 1.78 is no more a switch than 1.24.
   **The generalisable finding is the eleven sessions**, not the factor: a
   one-line comparison was not blocked by its difficulty but by an undeclared
   signature in one of its two inputs.
12. **Re-deriving Chapter 4 on the exact charge relation (4.28)** (added
   2026-10-08, and **the new top item of the RF thread**). The declaration in
   item 7 makes
   `V_g − V_dirac − V_ch = e·n/C_ox + E_F(n)/e` the exact relation; §4.3
   computes its linearisation. Running Eq. (4) on the exact form moves `I_d` at
   the RF bias by **−0.6023 %**, which is **1.4 times** the transport term
   items 7 and 8 were competing over. Two sessions went into a +0.43 % term
   while this sat underneath it, reachable in one line. The item is to re-derive
   §§4.3–4.6 on (4.28) and re-report, keeping the (4.29) numbers annotated in
   place. **It is blocked by item 9**: (4.28) and (4.29) both carry `n` as a
   magnitude, so re-deriving before the signed-carrier treatment exists would
   rebuild the branch-sign flaw in a new equation. Item 9 is therefore promoted
   from a sibling of item 7 to the **prerequisite** of the top item.
   **2026-10-09: UNBLOCKED AND UNTOUCHED.** Item 9 is closed and
   `n_net_exact()` supplies the branched root (odd bitwise, residual 7.690e-16).
   This item has no prerequisite left and is the **top item with nothing in
   front of it** — the first time in three sessions that is true. It is also now
   the only remaining reason the §4.6.5 rewiring item is deferred, so it blocks
   work rather than merely ranking above it.
   **2026-10-10 — CLOSED**, §4.6.8,
   `graphene_exact_charge_rederivation.py`, 21/21 checks. The −0.6023 % is
   reproduced to four decimal places from an independently written patch path,
   and **the item's own description of its two equations was wrong**. (4.29) is
   not the linearisation of (4.28); the two do not share a leading order as
   `dV → 0` (exponents 1.99999854 and 1.00000000 at `dV` = 10⁻¹⁰ V, quadratic
   against linear). (4.29) is (4.28)'s **exact differential relation used
   algebraically** — a one-point rectangle rule for the integral of
   `dn/d(dV) = (C_ox/e)·C_Q/(C_Q+C_ox)` with `C_Q(n) = 2e√n/A_F`, verified an
   identity to 1 ulp and by integrating it back to the closed-form root at
   observed order 3.994 — compounded with a **second, independent** error:
   `C_q` evaluated at the overdrive rather than at `E_F/e`. Separated for the
   first time, the argument error (+4.7416 % of the charge at the RF bias) is
   **larger** than the quadrature error (+2.6396 %), and the two are not
   additive. *[Both figures superseded 2026-10-11 — see item 16 and §4.6.9(c):
   +2.6408 % and +0.7603 % against a reference at the model's own temperature.
   The "larger" and "not additive" claims both survive; the magnitudes do
   not.]* §§4.3–4.4 move by at most 2.0859 % with every qualitative claim
   intact. **The item was also mis-scoped**: under the differential identity
   the gate capacitance *is* the charge relation's derivative, so `C_gs` moves
   too, and in `f_T` the charge half (−1.9986 %) and the capacitance half
   (+1.9392 %) nearly cancel to −0.0594 % — two coupled 2 % corrections, not a
   negligible effect. The three sentences in this thesis that call (4.29) a
   linearisation (§4.6.6, the §4.3 annotation, and this item) are marked
   superseded in place, not deleted. The §4.6.5 rewiring item is therefore
   **unblocked**, and the argument error this exposed becomes item 16.
13. **The `cosh` overflow in `quantum_capacitance()`** (added 2026-10-08).
   `log(2(1 + cosh η))` returns `inf`/`NaN` once `η = dV/(kT/e)` passes ~710,
   i.e. above `dV` = **18.3493 V** at 300 K, measured by bisection. The exact
   large-`η` limit is `|η| + log 2`. No committed number is affected — the
   maximum drive swept anywhere in this thesis is ~2.7 V — so this is a
   *latent* defect and is ranked last deliberately. It is recorded rather than
   patched because `graphene_fet_model.py` owns committed transcripts across
   Chapters 4, 5 and 6, and a numerics change there needs its own
   output-neutrality pass of the kind 2026-10-07 ran for the f-string fix.
   **2026-10-09 amends the "no committed number affected" clause, and the
   ranking stands.** On its first use from outside the module that owns it, the
   overflow produced a wrong number: above the onset `C_q` is `inf`, so
   `C_q/(C_q + C_ox)` is `inf/inf = nan`, `nan > 0` evaluates `False`, and the
   first draft of §4.6.7's 1 %-crossover bisection scored the `NaN` branch as
   "the floor dominates", walked away from the root and returned its own 200 V
   bracket ceiling as a result. The correct value is 14.73 V. Caught before
   commit because the printed figure was 74× the swept range; nothing reached a
   transcript. **Latent and harmless are not the same property** — the defect is
   still latent with respect to committed numbers, and is now known to be a trap
   for every new caller. The ranking is unchanged for the output-neutrality
   reason above, but the *mitigation* changed: callers outside this module must
   treat non-finite `C_q` as a hard stop rather than as a magnitude.
   **CLOSED 2026-10-11 (§4.6.10), and the ranking was wrong for a reason worth
   keeping.** This item was ranked last twice on the grounds that η = 710 needs
   18.35 V of overdrive, which the device never sees. But η also grows as `T`
   falls, and the `T → 0` exact-limit check that validates item 16 drives η to
   7e5 at the ordinary RF bias and returned `nan`: **the item ranked last was
   gating the item ranked as the largest open cost.** Latency had been judged
   against the inputs the *model* takes and not against the inputs its
   *validations* take, and that distinction is now a standing caution on any
   future latency claim and part of item 15's remit. The fix is the identity
   `ln(2(1+cosh η)) = |η| + 2·log1p(exp(−|η|))`, from
   `2(1+cosh x) = 4cosh²(x/2)`, installed in `_thermal_rounding()`. The
   output-neutrality pass this item had been waiting for is in
   `graphene_quantum_capacitance_overflow.py`: `C_q` itself moves by 1 ulp at 4
   of 400 biases near the Dirac point and **no committed quantity moves by a
   single bit** over 3603 values. Two by-products are recorded in §4.6.10: the
   pass measured the wrong thing in each of its first three drafts (literal
   constants, then product regrouping, then a patch that silently matched
   nothing), which is why the fix is a *substitutable* one-function kernel
   rather than a transcription; and the installed spelling is **not** more
   accurate — adjudicated at 60 digits it is strictly worse at 81 of 401
   samples and better at 45, with only the worst case improved — so being
   finite is the whole of the case for it.
14. **Measuring what a declaration costs, as a standing habit** (added
   2026-10-08; the mechanisable form of §4.6.6's second half). A declaration
   that resolves an ambiguity is a promise that the rest of the model belongs
   to the declared variable, and the size of that promise is computable: here
   it was one monkey-patched function and one re-run. Nothing in this thesis
   asks it of the other declarations it has made. Items 2, 3 and 6 all turn on
   a named quantity, and none of them has been costed this way.
15. **An audit of every tolerance in this repository, recording what fixed it
   at that value** (added 2026-10-09; the mechanisable form of §4.6.7's last
   subsection). Three exactness checks failed as first written today, and all
   three were measurements of something other than what they were pointed at:
   one measured the floating-point behaviour of this thesis's `V_g` sweep
   parameterisation (the ambipolar symmetry is bitwise in `dV`, and
   `(V_dirac − x) − V_dirac ≠ −((V_dirac + x) − V_dirac)` to 2.22e-16 V, which
   the `C_ox/e` lever turns into 2.0 /m²); one measured its own denominator
   (PASS and FAIL 119× apart on a 1-ulp residual, because `n_net → 0` where `S`
   does not); one measured four significant figures against 53 bits. **Three
   repairs were available at the moment of failure — loosen the tolerance, fix
   the model, or find out what it measured — and only the third produced the
   findings.** The first would have silently accepted all three. The item is to
   sweep every numerical tolerance in the repository and flag any that was set
   *after* a failure rather than derived from the precision being claimed, which
   is the only way to tell a justified tolerance from a tuned one after the
   fact. Sibling of item 6's magnitude census and of 10-06's ratio-criterion
   item, and the natural instrument for both: a tolerance set after a failure
   and a criterion written as a ratio are the two ways a check stops being a
   check without anyone editing it.
   **2026-10-10 — a near-miss of exactly this kind, caught and recorded.**
   §4.6.8's first draft required "every number moves by under 2 %" and failed
   at 2.0859 %. Widening it to 2.5 % would have produced a passing criterion
   that measured nothing but the run that set it. The bound was **deleted**
   instead: the magnitudes are printed with no threshold and the check tests
   the qualitative claims. The same session replaced a guessed 10⁻¹⁰ error
   bound on a quadrature identity with the **observed convergence order and
   the absence of a plateau**, which the method supplies and no one has to
   choose. The item now has two patterns to look for and one worked repair for
   each: *delete the bound and check the claim*, and *let the method set the
   criterion*. It is still unbuilt, and now **five sessions overdue**.

16. **Fixing `quantum_capacitance()`'s argument** (added 2026-10-10, and the
   first item on this list that arrives already costed). `C_q` is evaluated at
   the gate overdrive `dV`; the dispersion it is derived from requires
   `E_F/e = A_F·√n`, which at the RF bias is **4.04 %** of the overdrive.
   §4.6.8 separates this from the quadrature error of (4.29) and sizes it at
   **+4.7416 %** of the raw charge at that bias — **larger than the quadrature
   error (+2.6396 %), and the largest cost attached to anything currently open
   on Chapter 4**. The compact-model review
   [arXiv:2209.00388](https://arxiv.org/pdf/2209.00388) evaluates its `C_q` at
   the chemical-potential shift and never at the overdrive, so this is a
   deviation from the literature's own convention rather than a modelling
   choice this thesis is entitled to. It is ranked here rather than acted on
   for the same process reason as item 13 — `graphene_fet_model.py` owns
   committed transcripts across Chapters 4, 5 and 6, and a numerics change
   there needs its own output-neutrality pass — but unlike item 13 it is **not
   latent**: it is active in every committed number in those chapters. The
   honest statement is that it has the largest known cost on this list and is
   ranked low only for process reasons. Acting on it would be the first change
   in this thesis to **move** committed Chapter 4 numbers rather than annotate
   them, and the puddle floor's 4.3× attenuation of charge errors (§4.6.8) is
   the reason that is survivable.
   **ACTED ON 2026-10-11 (§4.6.9). The direction was right; all three numbers
   above are wrong, and the item understated the work rather than the size.**
   **PRECONDITION (i) DISCHARGED 2026-10-15 (§4.6.14).** The contact now owns
   a contact-correct replacement instead of being an exclusion, so the flip's
   two replacements both exist. Getting there required generalising §4.6.13's
   mechanism: `quantum_capacitance`'s runtime identity class spans **two
   operating points** — the channel's `C_ox` and the contact's `C_interface`
   (65.93×) — and identity answers *which names move together* without
   answering *to what*, so the class has to be partitioned by operating point
   and the partition has to come from the caller. The composed manager now
   leaves **nothing** at the shipped object (four sites installed, zero left
   alone, asserted by the mechanism and backed by a negative control).
   **What remains is exactly precondition (ii) — §4.7's prose and §4.10's
   annotation — which is writing, not plumbing.** The item is still OPEN and
   the default is still unflipped.
   **A PRECONDITION (iii) FOUND AND DISCHARGED 2026-10-16 (§4.6.15), WHICH NO
   ENTRY ABOVE HAD LISTED.** The flip's blast radius is established by
   regenerating committed transcripts and reading what moved, and the
   instrument of record for that is `graphene_pristine_transcript_audit.py`.
   At the start of 2026-10-16 that instrument — unrun since 2026-10-04 —
   reported **four of its own subjects STALE or ERROR**. A flip measured
   against that baseline could not have distinguished its own blast radius
   from twelve days of accumulated drift, so the word "exactly" in the
   sentence above was wrong when it was written: two preconditions were known
   and a third was in force unnoticed. It is discharged — 13 IDENTICAL, 0
   STALE, 0 NUMERIC_DRIFT, 0 ERROR — but the baseline now carries **two
   deliberately red suites** (item 20's 20 inert knobs and four dead names),
   and a flip run over a baseline with two failing guards would again mix its
   own effect with theirs. **This is why item 20 is ranked ahead of this
   item** despite this one carrying the larger physics cost: the ordering is a
   measurement constraint, not a change of priority. Preconditions (i) and
   (ii) are unaffected.
   (i) The argument cannot be fixed by changing an argument: `E_F/e` depends on
   `n`, which is what `C_q` is being used to compute, so the fix is the
   self-consistent solve (4.32). (ii) The +4.7416 % was obtained by
   substituting the T = 0 dispersion `C_Q(n) = 2e√n/A_F`, which deletes the
   thermal rounding `quantum_capacitance()` exists for: it gives `C_q(0) = 0`
   and an induced charge of *exactly* zero at the Dirac point, and makes
   `n ≡ 0` a solution through the origin — so temperature is what makes the
   relation well-posed there. (iii) Measured in a 2×2 at **one** temperature,
   against an initial-value reference validated to −1.5e-14 against (4.28)'s
   closed form, the argument error is **+2.6408 %** (not +4.7416 %), the
   quadrature error **+0.7603 %** (not +2.6396 %) and the total **+3.0963 %**
   (not +5.2064 %): overstated by 1.80×, 3.47× and 1.68×, because a 300 K
   model had been differenced against a T = 0 relation. The **ordering
   survives**, so this item remains correctly ranked ahead of the quadrature
   error. Installed, it moves the raw charge −2.2659 %, `n`(3.5 V) −1.0212 %
   after the floor's 0.451 attenuation, `I_d` at the RF bias −0.2530 %, `R_ch`
   at the Dirac point 0.0000 %, and **no verdict at all**; `f_T` moves
   +0.8065 % as another near-cancellation of a −1.0911 % charge half against a
   +1.8976 % capacitance half, and `f_max/f_T` +0.0009 %, leaving §7.8.1a
   untouched by a sixth lever. **The item does not close.** The default in
   `graphene_fet_model.py` is still (4.29); flipping it moves every committed
   number in Chapters 4, 5 and 6 and needs its own session with its own
   transcript regeneration. What changes today is that the flip is now *costed
   on the right form* and the blocking defect (item 13) is gone. One trap is
   recorded: today's 300 K argument error (+2.6408 %) sits within 0.0012 % of
   2026-10-10's T = 0 quadrature error (+2.6396 %), two different quantities
   that read as a reproduction.
   **PRECONDITIONS DISCHARGED 2026-10-14 (§4.6.13), and one of the two was
   stated wrongly on 2026-10-13.** The flip was deferred on 2026-10-13 for a
   named reason — nobody knew what a rebinding of the charge path reaches —
   and that reason is now gone, replaced by a mechanism rather than by a
   list. (i) The set a patch must reach is the name's **runtime identity
   class**: every `(module, attribute)` holding the same object, which is
   discoverable by identity and needs no import graph. `carrier_density`'s
   class has **2** members, `quantum_capacitance`'s **2**, and all three
   charge-path rebinding sites in the repository now compute it at entry
   instead of listing it (`graphene_binding_class.rebind`), verified
   behaviour-neutral by **byte-identical** regeneration of all three
   committed transcripts. (ii) **But the flip cannot be a rebinding
   everywhere.** `graphene_contact_doping_model.quantum_capacitance` is
   reached at a contact whose series capacitance is `C_interface`, 65.93×
   `C_ox`, where the series factor runs 0.555–0.806 against 0.999 in the
   channel; the object the owner would hold after the flip solves for `n`
   against `C_ox`, so copying it to the contact installs a **channel solver
   at the contact** — a different wrong answer, not a repair. The flip
   therefore needs **two** replacements, not one rebinding, and the
   contact-correct form is the one measured on 2026-10-13 (−12.4964 % in `n`
   for Pd). The exclusion is recorded in code with its reason attached, and
   the rebinder **refuses an exclusion without one**, because an unstated
   exclusion is how all three 2026-10 delivery faults were introduced.
   The only precondition still open is the one that was never numerical:
   §4.7's prose and §4.10's annotation have to be **rewritten**, since "only
   Pd gives a small, plausible positive residual" becomes false.
   **ANNOTATION 2026-10-18 (§4.6.16). The sentence above is kept verbatim and
   is false, and it was false on the day it was written.** Commit `7604c05`,
   dated **2026-10-13** — the day *before* this annotation's own 10-14
   paragraph — is titled "the two sentences in Sections 4.7 and 4.10 that the
   measurement makes false, annotated in place", and both annotations are in
   Chapter 4 with their numbers: §4.7's names −1.2583 % (Cu), +20.8352 % (Ni),
   +39.2281 % (Au), **+47.6925 % (Pd)** and `+51.0 → −203.276`, states that
   "Only Pd gives a small, plausible positive residual" is "true of the
   committed model and false of the corrected one", and leaves the sentence
   standing "**pending the item 16 flip**"; §4.10's amendment says the 9.6 %
   margin "stands as computed" and is a margin against *input* error that a
   known *modelling* error already exceeds fivefold. So a session recorded a
   writing precondition as open while its own tree carried the writing, and
   the 2026-10-17 entry then carried it forward as "still untouched and still
   writing, not plumbing" for a fourth day.
   **What is actually open is narrower and cannot be done first.** The
   annotations are written *in the subjunctive* — they describe what the
   corrected model would say and explicitly defer to the flip. The writing
   that remains is the writing that only the flip makes possible: when
   "pending the item 16 flip" itself goes stale, the committed table becomes
   the superseded column rather than the live one, and §4.10's "margin against
   input error" becomes a statement about a margin that was consumed. That is
   why 2026-10-17's operational advice — precondition (ii) and the flip in the
   same session — is right while its reason is wrong, and the difference is
   not pedantic: a session that believes (ii) is missing can discharge it by
   writing a sentence asserting a state of the code that does not exist, which
   is the §7.5.4 "created by a fix" class applied to prose. **Item 16's
   remaining content is therefore the flip and nothing else**, and it is
   ranked first for the next run on that basis.
   **AMENDED 2026-10-21 (§4.6.17). "The flip and nothing else" was wrong in
   one specific way, and the way is worth more than the flip.** Composition
   was left as well. The two halves of item 16 are each committed, each
   validated and each correct alone, and **installing both at once does not
   give form `B_T` everywhere**: `n_contact` lands **4.86× (Ni) to 10.77×
   (Cu)** past the effect this item exists to install, in the **same
   direction**, so it would not have read as a sign error. Measured:
   −51.3850 % to −61.9143 % against the committed contact-only form `B_T`,
   where the item's own effect is −4.7696 % to −12.4964 %. The cause is one
   attribute read. The contact replacement is deliberately built on the
   **shipped** `gfet.quantum_capacitance` — for two good reasons recorded in
   §4.6.12 (a 9.48e6 ulp constants gap in the alternative, and the overflow
   item 13 closed) — and reads that name by attribute lookup at call time,
   which is exactly the name the channel half rebinds. The contact's
   self-consistency is then solved against a `C_q` that has itself already
   been solved, against the wrong series capacitance. **No number above is
   revised**: every cost this entry quotes is reproduced exactly by the
   composed flip, which is the evidence that the defect is in the composition
   and not in either half.
   The general form is the part that outlives item 16, and **no guard in this
   repository currently covers it** [*as of 2026-10-21; `graphene_replacement_
   read_audit.py` covers it as of 2026-10-22, in the broadened cross-flip form
   — see §4.6.18*]. `bc.rebind()` asks whether all bindings
   of a name move together. `bc.rebind_by_operating_point()` asks whether any
   member of the class is unassigned. Neither asks **whether a replacement
   reads a name the same flip rebinds**. Such a replacement is neither stale
   nor unassigned — it is correctly installed, and reads a different object
   than the one it was validated against. That is a third failure mode beside
   the two the 2026-10-14 and 2026-10-17 mechanisms were built for, and it is
   invisible to both.
   Two further findings. (a) The failure decomposes into **(i)** the inner
   `C_q` and **(ii)** the inner series capacitance — `qca.series_factor`
   reads `gfet.C_ox` at call time and the contact solve has it substituted —
   and the two **partially cancel**, by +1.02 to +2.30 percentage points, so
   repairing (ii) alone makes the error *larger*. (b) The repair is to bind
   the contact replacement to the object **captured at entry** rather than to
   the live name; it reproduces the committed contact-only form `B_T`
   **bitwise** on all seven metals, and the channel half reproduces
   `on_consistent_argument()` **bitwise** on all eight census entries.
   **Status unchanged: the item is still OPEN and the default is still
   (4.29).** What changed is that the flip is now built, validated and costed
   as a single installable object (`graphene_default_flip.on_default_flip()`),
   so what remains is editing `graphene_fet_model.py` — a decision about what
   the thesis claims, not a build. With both halves in the circuit for the
   first time, Pd's implied transmission residual moves **+50.9489 →
   −203.2764 Ω·µm** and the additive decomposition fails **4 of 4**, which is
   §4.10's margin consumed exactly as the paragraph above predicted.
   **AMENDED 2026-10-22 (§4.6.18). The guard now exists, and it reports two
   things the 2026-10-21 amendment above did not ask for.** (a) The audit
   confirms the exposure as framed: of the four replacements, exactly one —
   `blast._cq_consistent_contact` — reads a name rebound while it is
   installed. **No number above is revised.** (b) But the question as the
   amendment above words it, "*whether a replacement reads a name the same
   flip rebinds*", answers **empty** for that very replacement: the name it
   reads is rebound by a *different* manager, and the two are concurrent only
   because this item composes them. A guard written to that wording would
   have passed the composed flip. The relation is replacement × every flip
   installable **concurrently**. (c) A second exposure, outside the wording
   entirely: all five replacements **read** `gfet.C_ox` at call time and
   three of them **write** it, and `gfet.C_ox` is in no identity class, so
   neither existing mechanism nor this item's own question can see it.
   (d) The decomposition into (i) and (ii) is not two defects — **192**
   channel arrivals under the unrepaired flip, **181** of them with `C_ox`
   already substituted, **0** under the repaired flip — one cause with two
   consequences, which is why one change repairs both. (e) The repair is
   **contingently** correct: a capture protects the name it holds, not the
   names the held object reads, and it suffices here only because the shipped
   `gfet.quantum_capacitance` reads nothing rebound. **Status unchanged: the
   item is still OPEN and the default is still (4.29)**, and the ranking
   consequence is unchanged too — the flip should not be landed in the same
   session as this audit, and the audit is now done.
17. **A two-knob check that reports the ratio, not smallness** (added
   2026-10-10 as "a check that can see a cancellation", **amended 2026-10-11
   on its first use**, §4.6.9(f)). Every check in this repository asks about
   one quantity, so a number that is small because two larger numbers oppose is
   invisible to all of them. The prescribed form was: move each input
   separately and flag the case where the combined move is *smaller* than
   either. Its first instance has the **opposite sign** — the argument and
   quadrature knobs **reinforce**, combined/(sum of separate) = **1.109** — so
   the check as specified would not have flagged it. It must report the ratio
   and flag any departure from 1 in either direction. And a requirement that
   could only be found by running it: **the interaction measure inherits the
   reference error of whatever it is differenced against.** The same two knobs
   give a departure from additivity of 0.717 against the T = 0 relation and
   0.109 against a 300 K reference — an inflation of **6.6×** — so a two-knob
   check must be run at the model's own operating point. Still unbuilt as a
   sweep over every committed number; it now has one worked instance and two
   requirements.
   **AMENDED AGAIN 2026-10-12 (§4.6.11(b)), on its second use, and the
   amendment is to the null rather than to the statistic.** The prescribed
   form reports combined/(sum of separate) and flags any departure from 1.
   Today's pair — the thermal rounding and the disorder floor, both
   regularising `C_q` at neutrality — gives **0.7148** against the additive
   null, a departure of 28.5 %, which the check would flag as a large
   interaction. It is not one. Regularisers entering through a monotone
   function of a quadrature sum **compose max-like**, so the additive null is
   simply the wrong null for them: against `max(alone)` the same measurement
   is **1.0253**, a departure of 2.53 %, which is the honest size.
   **A two-knob check must therefore state which composition law it is testing
   against, because the same numbers read as a 28.5 % interaction or a 2.5 %
   one depending on a null that was never chosen explicitly.** Three
   requirements now, and still unbuilt as a sweep: report the ratio, run it at
   the model's own operating point, and name the null.
18. **Costing the double regularisation of `C_q` at neutrality** (added
   2026-10-11, §4.6.9(h); the first item on this list to arrive from a
   *literature* source rather than from an internal audit). A second,
   independent compact-model lineage — the Southampton bilayer-GFET model —
   confirms that `C_q`'s argument is the surface potential and that the
   potentials are solved as a coupled system, which is item 16's claim held by
   two groups rather than inferred from one. But it regularises `C_q` at the
   Dirac point with the **residual disorder density** and not with temperature:
   `C_q,min = e²√n₀/(√π·ħ·v_F)`, with temperature never entering its `C_q`
   at all. **This thesis applies both** — a thermal rounding inside `C_q`, and
   `n_puddle` added in quadrature to `n` — and at `n_puddle = 5e15 m⁻²` the two
   are the same size: 9.710793e-03 against 8.437399e-03 F m⁻², a ratio of
   **1.1509**. They would agree exactly at `n₀ = 3.77e15 m⁻²` (0.755×
   `n_puddle`) or at **T = 345 K** (1.151× 300 K), so the agreement is not an
   accident of one parameter: SiO₂-supported graphene at room temperature sits
   where disorder and temperature regularise `C_q` equally. Using both is
   therefore a **double count of one physical effect**, worth about a factor 2
   in `C_q(0)`, and it has never been costed. It would not be true for
   suspended graphene (`n₀` two orders lower) or at 4 K, which is why neither
   convention is obviously wrong in the literature and why the overlap here is
   easy to miss. Related to, and partly the other side of, the uncosted third
   knob §7.10 now names: the rounding *form* itself has never been compared
   against Fermi–Dirac integrals.

   **COSTED AND CLOSED 2026-10-12 (§4.6.11). The item was wrong twice over,
   and the premise above is annotated rather than deleted.**
   (i) **The committed model does not double-count at all.** The two
   regularisers reach *disjoint* quantities. `C_q`'s argument in (4.29) is
   `dV`, and `dV` is never floored, so `n_puddle` cannot reach `C_q`;
   conversely `n_eff = (C_ox·dV/e)·f(C_q)` has `dV` as a **factor** and
   `f ∈ (0,1)`, so `n_eff(0) = 0.0` **bitwise** under all four `C_q`
   conventions tried and `carrier_density(V_dirac)` is **exactly**
   `n_puddle` — the thermal rounding cannot reach the Dirac-point density.
   Applying one regulariser to each of two quantities is a **reuse of one
   parameter for two jobs**, a modelling choice to defend; it is not a double
   count, which is an arithmetic error. The item was created by comparing two
   *literature forms* without first asking which quantity each term in **this**
   model reaches.
   (ii) **In the one spelling that does double-count the cost is +2.5262 %,
   not a factor 2.** Form B_T with the floor inside `C_q`'s argument gives
   `C_q(0) = 1.991221e-02` against 1.942159e-02 for disorder alone and
   8.437399e-03 for thermal alone: **both/max = 1.025262**, both/sum =
   0.714750. **Two regularisers entering through a monotone function of a
   quadrature sum compose max-like, not additively** — the larger swallows the
   smaller. The "factor 2" above was the **ratio of the two alternatives**,
   which prices *choosing wrongly between* them and not *applying both*. First
   instance of the 2026-10-06 ratio item in which the ratio was not a
   criterion but an estimate of a magnitude. Downstream: **+0.1527 %** in `n`
   at the RF bias, **0.0000 %** at the Dirac point (bitwise in all four
   conventions), and `σ_min` identical to the last bit — so **no committed
   conductivity, resistance or on/off number in Chapters 4, 5 or 6 is exposed
   to this item at neutrality.** It is real but small, and smaller at the
   operating point than at neutrality, which is the opposite of how it was
   ranked.
   (iii) **And the floor carries an exact factor 2 of its own.** The
   Southampton `C_q,min = e²√n₀/(√π·ħ·v_F)` is **exactly half** this thesis's
   own DOS-exact `C_Q(n₀) = 2e²√n₀/(√π·ħ·v_F) = 2e√n₀/A_F` — **0.000 ulp from
   ½** with one consistent `ħ` — because graphene's `D(E) = 2|E|/(π ħ² v_F²)`
   carries the `g = 4` spin–valley degeneracy and the literature form does
   not. So the ratio **1.1509** above is measured against the halved
   spelling; on this thesis's own convention it is **2.3018** and "the same
   size" weakens to **"the same order"**. The 1.1509 is **not** rewritten,
   because which spelling is right depends on what `n₀` denotes in the source
   — making this the only item on this list whose residue is resolved by
   *reading* rather than by computing.
   What survives is not a double count but the standing `n_puddle` question,
   now reached from two opposite directions (§4.6.7 from the floor being
   2.07× Chapter 3's, item 18 from the literature's use of the same symbol for
   DOS broadening).
   **AND A BOUND, FROM THE SOURCE'S OWN COMPOSITION LAW** (§4.6.11(f), added
   the same session after fetching the paper to settle what `n₀` denotes — it
   does not: the paper cites `n₀` to two references, never calls it total or
   per-species, carries no degeneracy factor and gives no DOS at all, so the
   factor 2 of (iii) cannot be reconciled by reading and both figures stay in
   the record). Its Eq. (3) composes the floor with the gate-driven part as
   `(C_min² + 2h²)/√(C_min² + h²)`, `h = C_var/2` — a **smooth max, not a
   sum**: `C_min` bitwise at zero surface potential, `C_var` to twelve digits
   when that dominates, and an excess over pure max that is maximal exactly
   where the two terms are **equal**, where it equals `3/√5 = 1.341640786`.
   So **the literature applies `n₀` twice itself** — Eq. (2) floors `C_q`,
   `n = √(n₀² + n*²)` floors the density, no temperature anywhere — making
   this a comparison between two *allocations* of regularisers to quantities
   rather than between a model that double-counts and one that does not. And
   the item's intuition was **right about where the worst case is and wrong
   about how big**: equality is indeed the worst case, and there the excess is
   **+34.164 %**, against the **+2.5262 %** this thesis incurs at 2.30:1.
   **A quadrature-composed regulariser pair cannot exceed the larger of its
   two terms by more than 34.2 %, so the "factor 2" was unreachable in
   principle and not merely unmeasured.** Also: the source's fitted
   `n₀ = 1.2e16 m⁻²` is **2.4×** this thesis's `n_puddle`, a second literature
   anchor on the standing question and pointing the same way as §4.6.7's.

19. **Making the rebind set a computed object rather than a remembered
   list** (added and largely closed 2026-10-14, §4.6.13; the second item on
   this list to be created by a *rate* rather than by a single finding).
   Three times in ten days a monkey-patch in this repository reported a model
   it was not running, each time because a second module held its own
   import-time reference to the patched name: the frozen-default fault in its
   monkey-patch form (2026-10-03), `contact_resistance_crossover.carrier_density`
   (2026-10-10) and `graphene_contact_doping_model.quantum_capacitance`
   (2026-10-13). Each was repaired by adding one line to one manager, and
   2026-10-13 stated why that is insufficient — **a repair to one binding is
   not a repair to the class** — without measuring the exposure. Measured:
   **149** module-level `from`-imports of repository names across **73**
   distinct `(owner, name)` pairs, of which **25** sites were left stale by
   the repository's **14** hand-listed rebinding sites. **Twenty-two of those
   had never been found by hand.** Of the ten distinct stale bindings, **seven
   are measured USED** by their own holder — the holder's numbers move when
   only its local copy is perturbed, which needs no manager to be entered and
   so covers the ten inline `try/finally` sites that cannot be entered at all:
   `cd.quantum_capacitance` at 9.88 (988 %), `cd.n_puddle` at 35.8 %,
   `crx.carrier_density` and `crx.W` at 27.0 %, and the three Chapter 6
   photodetector modules' `mu` at 9.2 %, 5.3 % and 3.9 %. Three are cosmetic:
   `cd.mu`, `cd.C_ox` and `cd.W` are imported and never reach a computed
   number.

   **And the import graph is the wrong instrument, in both directions.**
   (a) It under-reports. A runtime identity control found **12** second-order
   *renaming* aliases invisible to any rule that reads import statements
   (`double_regularisation.C_OX is gfet.C_ox`; `N0_SHIPPED is gfet.n_puddle`;
   `W_CROSS_VACUUM is W_GRAPHENE`), each a real site a patch of either the
   owner or the first-order holder leaves untouched. Identity sees them for
   free, because what it yields is not a directed edge but the **equivalence
   class** of names sharing one object — and for the question "what must a
   patch rebind", the class is the safe criterion and the import graph is the
   unsafe one. (b) It over-reports, and this correction retires a claim made
   earlier the same day. A miss needs the stale holder to be **imported in the
   patching process**. Measured per site in a fresh interpreter: of the 22
   remaining static misses the holder is **PRESENT for 1** and **ABSENT for
   21**. `import graphene_velocity_saturation_model` brings in three
   repository modules and none of the four `mu` holders is among them;
   `graphene_potential_declaration` does not import
   `contact_resistance_crossover` at all. So the morning's reading of that
   last site as "the 2026-10-10 fault unrepaired in a third module" was
   **wrong about this repository while right about the shape**: there was no
   stale binding to leave behind and no wrong number. The 25 is retained in
   place as an **upper bound** rather than restated, and the honest standing
   is **three independent axes** — PRESENT, USED, REACHED — of which nothing
   currently in the repository clears all three. The exposure is a **landmine
   count, not a defect count**, which is a weaker claim than the one the
   morning's commit made and the only one the measurement supports.

   **A third instrument fault, found and fixed in the same session, and it is
   the §7.5.4 class ("created by a fix") in its purest form yet.** Retrofitting
   the class-complete rebinder made all three charge-path sites **invisible to
   the audit**, whose only notion of a rebinding site was `alias.name = ...`;
   the validation that exists to keep 2026-10-13's finding on the record
   failed because the manager it names had stopped being a rebinding site in
   the only sense the instrument understood. **A repair that blinds the audit
   is not a repair.** The audit now recognises both forms; the validation is
   rewritten to assert the same fact in the form the repair left it in — zero
   *silent* misses, plus one exclusion named with a reason of substance — and
   its original no-silent-miss clause is kept rather than replaced. A fourth
   instrument fault was caught by the audit's own planted-binding control,
   which failed while nothing was wrong with the repository: the module
   scanner skipped every file whose name begins with an underscore, so the
   sweep reported zero bindings in a file it had itself just created with
   one. What remains open is small and named: `graphene_contact_doping_model`
   needs its own contact-correct replacement (item 16(ii)), and the eleven
   surviving hand-listed rebinding sites outside the charge path — `mu`, `W`,
   `n_puddle`, `C_ox` — have not been converted, because each is latent by
   import today and converting a site changes nothing until it is not.
   **ANNOTATION 2026-10-15 (§4.6.14). The first of those two is now done, and
   doing it found a limitation of this item's own repair.** The repair gives
   every member of an identity class the same replacement, with one escape
   hatch that leaves a member alone for a stated reason; its vocabulary is
   therefore `{same replacement, left alone}` and it cannot say *rebind this
   member to a different object*. That is what a contact with a replacement of
   its own needs, so `bc.rebind_by_operating_point()` was added alongside
   `bc.rebind()` (which is **untouched**, and whose three transcripts
   regenerate byte-identically). **The claim this item makes — that the
   rebind set is computable and needs no list — survives intact; what fails is
   the unstated assumption that one class takes one replacement.** The
   refusal carries across in the form the lesson takes here: an **unassigned
   class member** is refused exactly as an unstated exclusion is, because
   silently leaving a member at the shipped object is the identical fault
   reached by omission from a map instead of from a list. The eleven
   hand-listed sites outside the charge path are still unconverted.

20. **Converting the 20 inert knobs, and the four dead names behind them**
   (added 2026-10-16, §4.6.15; the third item on this list created by a
   *rate* rather than by a finding, and the first created by three
   instruments agreeing that they had not been read). On 2026-10-04 the
   inert-knob census — a module-level constant captured as a function
   parameter default, where the default binds once when the `def` executes,
   so a rebinding never arrives and a mutation audit reports a pass while
   measuring nothing — was driven to **zero**. Twelve days later it reads
   **21 syntactic captures**, in front of three separate censuses, none of
   which was being read; the committed transcripts of all three were last
   written on that same 2026-10-04, and so was the transcript of the
   instrument whose only purpose is to notice that
   (`graphene_pristine_transcript_audit.py`). That instrument, asked today,
   reported four of its own subjects STALE or ERROR without being prompted.
   The 21 is **not** the count of the fault. Measured against exactly known
   values, a `def` in the module body or a **class** body executes once at
   import and its capture is FROZEN, while a `def` in a **function** body
   re-executes on every call of that function and its capture is **LIVE** —
   so the syntactic pattern is two classes with opposite behaviour under the
   one operation this repository uses to deliver every change. The split is
   **20 FROZEN + 1 LIVE**, reconciled to zero residue against both standing
   censuses and per `(file, constant)` as well as in total. What remains:
   **convert the 20**, which is mechanical and was deliberately not done in
   the session that specified it, because converting 21 as one class would
   have rewritten a live capture for no reason and recorded it as a fix. In
   the same twelve days `graphene_dead_name_sweep.py` went red on **one
   genuinely DEAD name** (`graphene_potential_census_audit.CHANNEL_POTENTIAL_READINGS`)
   and **three hidden dead names** — CROSS_MODULE resolutions resolved on
   *spelling* whose claimed readers do not import the module at all
   (`graphene_contact_consistent_argument.E_CHARGE`, added by the preceding
   session; `graphene_diffusion_current_model.HBAR`;
   `graphene_quantum_capacitance_argument.A_OX`), over 12 unsupported rows.
   The 2026-10-01 precedent for dispositioning a dead name — give it a
   reader, or declare it exempt with a stated reason — applies to all four.
   **Ranked first for the next run**, ahead of item 16, for a reason
   internal to item 16: see its precondition (iii).
21. **The standing list is not a closed system, and two of its items asserted
   the opposite of what the tree contained** (added 2026-10-18, §4.6.16; the
   first item on this list created by auditing *the list itself* rather than
   the code or the physics). Every session here chooses its work from the
   "Not yet covered" list of the preceding `AUTOMATION_LOG.md` entry. Nothing
   had ever checked a claim that list makes. Two of the 2026-10-17 list's
   twenty-five items assert a state of the tree that the tree denies: item
   16's precondition (ii), above, and **"CHAPTERS 2 AND 3 ARE STILL
   UNDRAFTED … the largest remaining block of pure writing"**, when Chapter 2
   was drafted 2026-09-28 and Chapter 3 2026-09-29 and both have been in the
   tree for nineteen days.
   The second is the more instructive, because **the list had already got it
   right**. From 2026-10-01 to 2026-10-06 it carried the corrected form
   ("Chapters 2 and 3 *are drafted* but §4.8.2's form of error is not audited
   for elsewhere"), carried nothing for the ten entries after that, and
   returned on 2026-10-17 to the pre-09-28 wording — and returned it to
   *first-rank contention*, since the entry closes it "worth ranking against
   item 16 deliberately rather than by default". So the one item offered as a
   competitor to the top measured item was describing writing committed three
   weeks earlier. **This is not a stale item carried past its repair; it is a
   retired wording revived**, and the distinction matters because the two have
   different fixes: staleness is cured by re-checking an item, revival is cured
   by making retirement checkable.
   **The cause is outside the repository, which is why no audit here could
   have found it by looking inward.** These sessions are driven by a scheduled
   brief that carries a frozen description of repository state dated
   2026-09-17, including the A1 sentence almost verbatim. A session that reads
   its own brief as current state will re-derive items the list has retired,
   and will do so with the authority of an instruction rather than the
   tentativeness of a memory. The standing rule this produces is checkable
   only in its second half, and both halves are recorded: **(a) procedural —
   a claim about the state of the tree that arrives from outside the tree is
   re-verified against the tree before it ranks anything; (b) mechanical —
   `graphene_standing_list_audit.py`**, which re-checks every standing item
   that asserts a tree state, and whose *retired* rows assert both that a
   retired wording stays gone and that the fact which justified retiring it
   still holds, so the suite cannot go quiet the moment its finding is acted
   on.
   Two by-products are worth keeping. The module's **first draft failed its
   own plant on its own subject matter** — it tested for a sentence by exact
   bytes, the markdown file wraps that sentence, and it therefore reported
   §4.10 "untouched" about a file carrying the amendment, i.e. it *confirmed*
   the item it was written to contradict. Resolving a claim on spelling is the
   same defect as the two docstrings of §4.6.15 and the 2026-10-17 grep, and
   it is now three for three that this repository's instruments fail first by
   answering a weaker question than the one asked. And
   **Design-Verification-Roadmap built `tools/open_item_staleness_audit.py`
   for this exact failure on 2026-10-13, and on 2026-10-17 used it to retire
   an item of its own for the identical reason** ("it names work finished
   before the item started being carried"). That repository caught its
   instance on the same day this one revived its instance. The difference
   between the two days was one module, and the generalisation — *an
   instrument built in one place does not protect the other* — is the
   cross-repository form of §7.5.2.
22. **Which `v_F` this thesis is quoting, and the fact that it cannot
   currently say** (added 2026-10-25, §4.6.21 (c)–(d); the first item on this
   list that is a question about the *physics* produced by an audit of the
   *code*). Two findings, and the second is a precondition for even asking
   the first.
   **(i) The Fermi velocity is two names, not one knob.** `C_q` carries it as
   `CQ_PREFACTOR = 2e³/(π ħ² v_F²)` and `A_FERMI = ħ v_F √π / e`. The contact
   replacement follows `gfet.v_F` with d ln `C_q`/d ln `v_F` = −2.323463;
   both channel replacements are **bitwise unreachable** from `gfet.v_F` and
   follow `qca.CQ_PREFACTOR` (+1.008232) and `qca.A_FERMI` (+0.519524), which
   `graphene_exact_charge_rederivation` computed from `gfet.v_F` at import
   time and copied in. A `v_F` sensitivity study done by editing `gfet.v_F`
   would have reported **channel insensitivity to the Fermi velocity as a
   result**. The repair takes three names moved together, one as *k*⁻²; it is
   a mechanical change and is not blocked on anything.
   **(ii) And the committed value is the bare one, on a substrate that
   renormalises it.** Hwang *et al.* measure 1.15 ± 0.02 (SiC, ε = 7.26),
   1.49 ± 0.08 (BN, ε = 4.22) and 2.49 ± 0.30 (quartz, ε = 1.80), ×10⁶ m s⁻¹,
   against a fully screened *ab initio* 0.85 × 10⁶ [Sci. Rep. **2**, 590
   (2012)]. This thesis commits 1.00 × 10⁶ throughout, just above the bare
   value, while §4.2 models a 90 nm SiO₂ stack whose ε ≈ 3.9 puts the
   substrate-appropriate value near 1.5 × 10⁶ — and ×1.5 is exactly the
   perturbation that moves the contact `C_q` by **−63.5273 %**.
   **What this item asks for is a decision, not a build.** Using the bare
   value consistently is defensible: most compact device models do, and the
   renormalisation is carrier-density dependent, so one number is already an
   approximation. The weak position is the present one, which quotes the
   contact `C_q` at a `v_F` the modelled substrate does not have *and* cannot
   express the alternative by changing one name. **(i) should be done first
   regardless of how (ii) is decided**, because (i) is what makes (ii)
   measurable, and because until it is done any number this thesis quotes for
   the `v_F` sensitivity of a channel result is zero for the wrong reason.
   **Annotated 2026-10-26 (§4.6.22).** The −63.5273 % of (ii) was measured on
   the contact replacement **in isolation**. It **survives** for the holder the
   default flip actually installs — on `flip.make_contact_cq(live=False)`,
   `gfet.v_F` carries −63.5273 % in *every* configuration measured, and
   `gfet.hbar` carries the bitwise identical figure, as `(ħ v_F)²` requires. It
   does **not** survive on `blast._cq_consistent_contact` composed with either
   other flip, where `gfet.v_F` is bitwise dead. (i) and (ii) are unchanged in
   substance; what is added is that the number has a stated configuration.
23. **The configurations §4.6.22 did not measure, and a verdict indexed by
   configuration** (added 2026-10-26, §4.6.22 (d) and (g)). §4.6.22 measures
   six installed objects against three flips and isolation, and finds exactly
   one non-composing row. Three things are open, in descending order of cost.
   **(i) The missing cells.** `live_contact=True` is the unrepaired contact
   replacement *as the default flip installs it*, which is the configuration
   §4.6.17 actually costed, and it is not in the 6 × 3 table. The
   *simultaneous* entry of two managers is not either, and it is the
   configuration a future chapter that flips both halves would run in. Both
   are mechanical: the table is built from `rra.FLIPS`, so the work is adding
   managers to it and supplying the verdicts.
   **(ii) A verdict indexed by configuration.** §4.6.22 (d) shows
   exhaustively that no single member of {PROTECTED, UNREACHED, EXPOSED,
   GOVERNED} is true in both configurations for `gfet.C_ox` on the live
   contact replacement. So the long-standing item asking for the
   non-composing compositions to be *adjudicated* is not merely unfinished —
   it is **not satisfiable in the present vocabulary**, and that is why three
   sessions have noted it without discharging it. The repair is a verdict
   that names its configuration, with the guard requiring one per
   (holder, flip) cell rather than one per (holder, name) pair. This is a
   change to §4.6.20's and §4.6.21's guards, not an addition to them, and it
   should be costed before being built.
   **(iii) The probe point.** "Invariant under its own flip" is asserted
   bitwise at a single `csg.PROBE_V` with a single `csg.PERTURBATION`. A
   corner sweep over both would turn the strongest claim of §4.6.22 (b) from
   a measurement at a point into a measurement over a box, and it shares its
   cheapest entry with the standing *is the admissible shape box the right
   one?* item.
24. **Does `n_c` differ between the electron and hole branches, and does the
   asymmetry therefore screen at all?** (added 2026-11-03,
   `graphene_kappa_impurity_bound.py` §7, and ranked **first** of the open
   items in this chapter.) Every other item here moves a number. This one
   decides whether a *qualitative* conclusion of Chapters 3 and 4 survives.
   §4.11.5 says the electron–hole asymmetry "dies off as the long-range
   impurity potential screens". That is a consequence of the Matthiessen form
   **plus** a single `n_c` shared between the branches, and the second part
   was never argued. Written per branch with a split `n_{c,e/h} = n_c(1±s/2)`,
   there is an exact identity at `s = 2r₀`: the `(1+|n|/n_c)` factor cancels
   out of the asymmetry ratio and `r(n) = r₀` at **every** density, with no
   screening whatever. At the top-gate operating point a 20 % split moves `r`
   by 524 %, so the range between "screened 6.2×" and "not screened at all" is
   spanned by an assumption nobody has tested.
   **Why the symmetric case is not the safe default.** `r₀` is derived from a
   beyond-Born cross-section ratio — a property of scattering off a potential
   of a given sign. Neutral defects are also potentials. There is no argument
   anywhere in this thesis for the asymmetry vanishing for vacancies and grain
   boundaries while surviving for charged impurities; it was assumed by
   writing one symbol.
   **What would settle it, and it is cheap.** Zion *et al.* read their
   `3 × 10¹² cm⁻²` off the **hole** branch (their Section III c). The electron
   branch's transition density is in the same figure. Reading both off Zion
   Fig. 5 gives `s` for one sample directly — and that is the **same picture**
   the standing *lower bound on ℓ* and *`Δ_c` form for the chemisorbed metals*
   items are blocked on, so one figure read discharges three items. This is
   now the cheapest high-value unblocking action available to a human reader,
   displacing the full-text fetch of arXiv:1801.03565, which was attempted on
   2026-11-03 and is **closed**: the text contains no mobility anywhere, so
   the pairing it was meant to resolve is unresolvable from that source and
   was instead resolved by a different route (item 25).
25. **The `A` pairing, closed 2026-11-03, and what replaced it.** Recorded here
   because it is the first item on this list to be discharged by finding that
   the datum it wanted does not exist. The 11-02 module obtained `κ = A/C`
   by pairing Zion's crossover density with *this* repository's mobility, and
   logged "read Zion's own mobility" as the fix. The full text has no mobility
   in it. But it has Zion's own charged-acceptor density,
   `N_A^(S) = C_ox·V_NP/e`, for the same sample — so `κ = n_c/N_imp` is
   readable off **one** device from **two** numbers in **one** paper, with no
   mobility and no cross-paper pairing. Because `N_A^(S)` is a net charge
   inferred from a compensating voltage, because Zion's own `N_t` oxide-defect
   density is acknowledged and unquantified, and because the film is
   transferred CVD, `N_imp ≥ N_A^(S)` and therefore `κ ≤ 2.089`. The 11-02
   value `κ = 2.4` is **excluded** by 14.9 %. Every `n_c` and `ΔV_c` in this
   thesis derived from it falls 13.0 %; every operating margin and impurity
   budget rises 14.9 %. No verdict reverses — the thinnest margin, `1.10×`,
   becomes `1.26×`, which is the useful direction. As a by-product the
   mobility the paper never prints is bounded anyway:
   `μ_{L,Zion} ≤ 0.348 m²/Vs`. **The methodological point worth carrying
   forward** is that the blocked item was discharged by asking what the source
   *does* contain rather than by re-attempting the fetch for what it was
   assumed to contain, and §7.5's taxonomy has no class for this: it is not a
   claim that failed, it is a *question* that was mis-posed.


---
