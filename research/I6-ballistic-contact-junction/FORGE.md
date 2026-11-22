# I6 — FORGE 2026-11-22: ballistic (Landauer) transport across the metal-doped contact junction vs the repo's diffusive R_extra

Seed S1 (forge brief). This idea uses I1's doping profile but changes the transport model on it. src's R_extra integrates a local diffusive σ(x) along the doping tail. I1 (11-21) found that clean-Au data need λ ≈ 5–15 nm, which is shorter than the typical mean free path in SiO₂-supported graphene (tens of nm). On that length scale the junction is ballistic, so a diffusive R_extra is the wrong model in exactly the regime the data point to.

**Candidate claim (to be tested):** with λ ≈ 5–15 nm, the per-metal junction resistance has to be computed as Landauer transmission through the p–p′/p–n potential step, including Klein angular filtering. This gives a different size, and possibly a different gate-voltage asymmetry, from the src diffusive R_extra. The proposal is to compare both with measured R_c(V_g) for clean Au (Malec 2011, Passi 2018), so the model choice is settled by data rather than by assumption.

## Stage 1 — Prior art

Searches: one standard search, "metal graphene contact resistance NEGF Landauer ballistic transmission doping profile". Two pages fetched (both OK, partly summarised): arXiv:1507.08412 (ar5iv) and Sci. Rep. 7, 5109 (PMC). The rows marked *not fetched today* are cited from the literature already in the repo or from standard references. They must be checked before any novelty claim.

| who | what they did | outcome | gap left | what we can reuse |
|---|---|---|---|---|
| Chaves, …, Neumaier, arXiv:[1507.08412](https://arxiv.org/abs/1507.08412) (2015) | analytic: Bardeen transfer for metal→graphene (R_mg) plus a Landauer step (Cayssol) for the lateral R_gg; 8 metals (Table 1) | R_c = R_mg + R_gg, with a secondary peak at V_D that broadening suppresses; Ni/Ti model 150–600 Ω·µm against their own TLM data | the step is abrupt (no λ-tail), and the parameters are fitted per metal (Table 1 asterisks); no physisorbed-metal clean data | the Cayssol step transmission formula; Figs 5–7 R_c(V_g) for digitising |
| Cusati, …, Iannaccone, Sci. Rep. 7, 5109 (2017), doi:[10.1038/s41598-017-05069-7](https://doi.org/10.1038/s41598-017-05069-7) | DFT (QE/PWCOND) + Landauer–Büttiker for Ni, Cu, Pd, Pt; own TLM devices | ballistic R_c: Ni 360/120, Cu 1153/374, Pt 209/123, Pd 2899/2076 Ω·µm (Table 2); fair agreement except Pd | they state the Dirac-point screening length is 1–4 nm; no Au and no λ ≈ 5–15 nm tail | Table 2/3 numbers as reference values; their statement that ballistic R still depends on the Dirac-point energy |
| Xia, Perebeinos, Lin, Wu, Avouris, Nat. Nanotechnol. 6, 179 (2011), doi:[10.1038/nnano.2011.6](https://doi.org/10.1038/nnano.2011.6) (*not fetched today*) | Pd contacts: R_c(V_g) data plus a model with a junction potential profile and ballistic transmission | known as the standard ballistic contact-junction treatment | Pd is chemisorbed (I3); no per-metal physisorbed sweep tied to a W → doping rule | **closest prior art**; the wedge must not be "ballistic junction transport" by itself |
| Cayssol, Huard, Goldhaber-Gordon, PRB 79, 075428 (2009), arXiv:[0810.4568](https://arxiv.org/abs/0810.4568) (*not fetched today*) | contact resistance and shot noise of a graphene p–n step | analytic transmission through a smooth/abrupt step | — | smooth-step transmission formula as a second check |
| Khomyakov et al., PRB 82, 115437 (2010), arXiv:[0911.2027](https://arxiv.org/abs/0911.2027) (already in repo, I1) | nonlinear screening of contact doping, with the gate distance as decay scale | doping tail shape vs gate | transport not computed in that paper | the tail shape for I1's gate-dependent profile |
| Malec & Davidović, PRB 84, 033407 (2011), arXiv:[1108.3238](https://arxiv.org/abs/1108.3238) (in repo, I1) | clean Au on SLG, R_c max at +120 V | data | — | the R_c(V_g) peak position as the comparison target |
| DFT–NEGF metal–graphene contact model, PCCP 15 (2013), [C3CP52589A](https://pubs.rsc.org/en/Content/ArticleLanding/CP/2013/C3CP52589A) (title only) | DFT–NEGF R_c | — | — | read only if the wedge is challenged |

**Wedge (narrow):** the same per-metal W → ΔE_F rule (I2, 5/5 signs) and the same tail profile (I1), solved ballistically and diffusively side by side. The question is which one matches the clean-Au R_c(V_g) data, and by how much the src R_extra is wrong in the λ ≈ 5–15 nm regime. This is a **method-consistency result** that supports I1 and the thesis. It is likely a section of the I1 paper rather than a standalone paper.

## Stage 2 — Mechanisms

1. **λ vs mean free path ℓ.** If λ ≪ ℓ, the junction is ballistic (Landauer); if λ ≫ ℓ, it is diffusive (src). I1 says λ ≈ 5–15 nm, against ℓ ≈ 20–100 nm on SiO₂.
2. **Klein angular filtering in a smooth p–n junction**, with transmission ∝ exp(−π ħ v_F k_F² / |e dV/dx|) for oblique modes. This gives the n/p asymmetry of R_c(V_g).
3. **Mode count in the channel (Sharvin limit).** R_min = h/(4e²) · π/(k_F W) sets a floor that the diffusive model does not have.
4. **Under-contact density (gate-tuned vs pinned; I1 11-21).** It sets the step height.
5. **Puddle broadening near the CNP.** It smears the Klein filter.

## Stage 3 — Options

- **Ambitious:** k_y-resolved Kwant or transfer-matrix calculation through the self-consistent gate-dependent profile for 5 physisorbed metals; R_c(V_g) vs Malec, Passi and Nouchi; ballistic–diffusive crossover vs ℓ.
- **MVP (≤ 7 EXECUTE days, can share code with I5's normal-state Kwant set-up):** 1D Dirac transfer matrix (or a Kwant ribbon) through f(x) = 1/(1+x/λ) at λ ∈ {5, 15, 250 nm} for Au (calc W). Compute R_bal(V_g) against the src R_extra(V_g) and compare the peak position and size with Malec 2011 (+120 V) and Passi 2018 (519 Ω·µm at −40 V).

## Stage 4 — Decision block

| field | decision |
|---|---|
| **Pick** | MVP, Au only first; physisorbed metals (Al, Ag, Cu, Pt) after |
| **7-day test** | (d1) check Xia 2011 for the exact claim; (d2) transfer-matrix code + validation; (d3) R_bal vs R_extra at 3 λ; (d4–5) R_c(V_g) vs Malec/Passi; (d6–7) RESULT + decide whether to fold into I1 |
| **Success metric** | at the λ the data prefer (5–15 nm), R_bal and src R_extra differ by > 2× **and** the ballistic one is closer to the measured Au R_c (Passi 519 Ω·µm, and the sign/position of the Malec peak) |
| **Kill condition** | R_bal / R_extra stays within 0.5–2× for λ = 5–250 nm (the model choice does not matter), or Xia 2011 / Cusati 2017 already give this per-metal physisorbed comparison → KILLED (or folded into I1 as a robustness check) |
| **Risks** | strong prior art (Xia 2011) means the wedge is narrow; I1's profile is itself being replaced (gate-dependent), so the I6 numbers may need a rerun; device-VM network is too slow for pip, so the computation runs in the cloud workspace |
| **First experiment** | k_y-resolved Dirac transfer matrix through f(x) for Au, with exactly ONE validation: flat profile (n_c = n_bulk) gives T(k_y) = 1 for every propagating mode, so R equals the Sharvin value h/(4e²)·π/(k_F W) |
