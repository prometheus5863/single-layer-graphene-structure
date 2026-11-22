# I5 — FORGE 2026-11-22: Andreev reflection at a graphene–superconductor contact, with the metal-induced doping profile

Seed S3 (forge brief). It serves the registered thesis title directly: quantum transport in twisted bilayer graphene–superconductor junctions. It also **supersedes I4**, because I4's next step (listing what Ch 4's contact model must supply to a BTK interface: barrier Z, doping, mismatch) is this idea's step 1.

**Candidate claim (to be tested):** most theory on graphene/SC Andreev reflection uses a sharp step U0 at the interface, or a CNP that varies on its own. A metal superconductor (Nb, MoRe, NbN, Al) deposited on graphene also dopes the graphene next to it: the I1/I2 contact-doping profile, with ΔE_F ≈ 0.1–0.35 eV and a tail λ ≈ 5–15 nm (I1, 11-21). Under the contact, the specular interband AR window needs |E_F| < Δ ≈ 1 meV, and this doping keeps that window closed. The proposed result is where (in gate voltage, bias and λ) a specular-AR signature survives, computed with Kwant BdG and compared with Efetov 2016's measured dI/dV maps.

## Stage 1 — Prior art

Searches (3, all standard): "twisted bilayer graphene superconductor junction Andreev reflection tight-binding", "graphene superconductor junction Kwant BTK retro specular crossover", "specular Andreev reflection bilayer graphene experiment". Two pages were fetched: arXiv:1605.04254 (abstract page only, the full text was not read) and arXiv:2211.11389 (ar5iv, about the first 40k of 246k characters verbatim, the rest summarised). Every other entry below is **title-only** from search results and must be read before any novelty claim.

| who | what they did | outcome | gap left | what we can reuse |
|---|---|---|---|---|
| Beenakker, PRL 97, 067007 (2006), arXiv:[cond-mat/0604594](https://arxiv.org/abs/cond-mat/0604594) | Dirac-BdG theory of a graphene NS interface | specular AR when E_F < Δ; retro AR when E_F ≫ Δ | ideal step interface; no contact doping profile | the analytic G_NS(eV, E_F) to validate Kwant against |
| D. K. Efetov & K. B. Efetov, PRB 94, 075403 (2016), arXiv:[1605.04254](https://arxiv.org/abs/1605.04254) (abstract read) | BdG + BTK for a bilayer graphene/SC interface | retro→specular crossover; conductance suppressed at \|eV\| = \|E_F\| < Δ; dip sharpness depends strongly on the step U0 | U0 is a free step, not derived from a contact model | the U0-sensitivity statement is exactly the lever our doping profile sets |
| Efetov et al., Nat. Phys. 12, 328 (2016), arXiv:[1505.04812](https://arxiv.org/abs/1505.04812) (cited in 1605.04254's abstract; not fetched today) | experiment: BLG/NbSe2 van der Waals interface | first observation of specular interband AR | vdW interface with no metal doping; deposited-metal SCs untested | **the data to digitise** (dI/dV vs V_bias and V_g); figure numbers to be confirmed |
| Sainz-Cruz, Pantaleón, Phong, Jimeno-Pozo, Guinea, arXiv:[2211.11389](https://arxiv.org/abs/2211.11389) (2022; partly read) | tight-binding BTK for STM-tip/superconducting-TBG junctions and TBG Josephson junctions, s- vs f-wave | f-wave gives a sub-gap peak; compared with Oh 2021 and Kim 2022 STM data | the normal side is an STM tip, so no graphene/metal contact doping; untwisted graphene junctions are called a direction for future work (summary only, unverified) | the TBG step: their three-channel BTK-as-TB set-up for the later TBG stage |
| arXiv:[1706.03457](https://arxiv.org/abs/1706.03457), "Andreev reflection in a bilayer graphene junction: role of spatial variation of the charge neutrality point" (title only) | presumably a CNP that varies in space at an NS junction | **unknown, closest threat** | may already contain the "doping profile kills specular AR" result | must be read on day 1; if it uses a metal-contact doping profile and compares with data, I5 is KILLED (prior art) |
| arXiv:[1801.00620](https://arxiv.org/abs/1801.00620), "Enhanced specular Andreev reflection in bilayer graphene" (title only) | — | — | — | read on day 1 |
| arXiv:[2306.00529](https://arxiv.org/abs/2306.00529), "Andreev and normal reflections in gapped bilayer graphene–superconductor junctions" (title only) | — | — | — | read if 1706.03457 is not decisive |
| arXiv:[2304.10741](https://arxiv.org/abs/2304.10741), "Revisiting Andreev processes in SGS Josephson junctions: comparison with experimental results" (title only) | — | — | — | possible second data source |

**Wedge:** derive the interface potential from the I1/I2 contact model (sign and magnitude of ΔE_F from W, tail λ from the clean-Au fit) instead of a free U0. Then ask a falsifiable question: for which superconductors (by W) and which λ does specular AR survive within experimental reach, and does the vdW NbSe2 case (with little doping) come out as the special one? Novelty is **unconfirmed**, because the 1706.03457 check is pending.

## Stage 2 — Mechanisms (the levers that decide the outcome)

1. **Contact doping ΔE_F under the SC** (W_SC − W_G and d_eq; I2). If |E_F| ≫ Δ under and next to the contact, AR there is retro whatever the gate does.
2. **Tail length λ.** A tail shorter than the coherence length ξ acts as a step U0 (Efetov & Efetov's dip sharpness). A tail longer than the Fermi wavelength at E_F ≈ Δ forms a p–n region that Klein-filters oblique modes.
3. **Electron–hole puddles δE_F ≈ 5–50 meV, against Δ ≈ 1 meV.** These limit any specular window. They are why the hBN-encapsulated data are needed.
4. **Interface barrier Z / Fermi-velocity mismatch** (Sainz-Cruz Fig. 1(c), velocity mismatch case).
5. **Monolayer vs bilayer vs TBG dispersion.** In the BLG parabolic band, E_F < Δ is easier to reach at the same density; the TBG flat bands are the later stage.

## Stage 3 — Options

- **Ambitious:** TBG (or BLG) / metal-SC junction in Kwant with a self-consistent contact doping profile. Map G_NS(V_bias, V_g) and compare with Efetov 2016 BLG/NbSe2 data and a deposited-metal-SC dataset. This becomes a thesis chapter plus a PRB/PRApplied paper.
- **MVP (7 EXECUTE days):** monolayer graphene NS junction in Kwant (translation-invariant along y, k_y-resolved), with a BdG on-site potential following the I1 profile (ΔE_F, λ). Step 1: validate against Beenakker. Step 2: sweep λ ∈ {0 (step), 5, 15, 250 nm} and ΔE_F from I2 for Nb/Al/MoRe, and report the specular-dip visibility. Step 3: digitise one Efetov 2016 dI/dV line cut and set it next to the λ→0, ΔE_F→0 (vdW) case.

## Stage 4 — Decision block

| field | decision |
|---|---|
| **Pick** | MVP above. Monolayer first, BLG second. TBG is out of scope for the first 7 days |
| **7-day test** | (d1) read 1706.03457 / 1801.00620 for prior art; (d2) Kwant NS validation; (d3–4) doping-profile sweep; (d5–6) digitise Efetov 2016 and compare; (d7) write RESULT + decide |
| **Success metric** | (a) Kwant reproduces Beenakker G_NS within 2 % at E_F = 0 and E_F ≫ Δ; (b) the doping-profile model gives a dip-visibility vs (ΔE_F, λ) map that puts NbSe2-vdW inside and Nb/Al-deposited contacts outside the specular window; (c) one digitised Efetov 2016 curve is reproduced in dip position (\|eV\| = \|E_F\|) within the data's resolution |
| **Kill condition** | 1706.03457 (or any paper found on d1) already computes AR with a metal-contact doping profile and compares with data → KILLED (prior art), cited. Or the doping profile changes dip visibility by < 10 % for every λ in 5–250 nm, so it is not a lever |
| **Risks** | Kwant needs Python 3.11 (`uv venv -p 3.11`, `uv pip install --no-build-isolation kwant` built 1.5.0 in the cloud workspace on 11-22; no MUMPS, so SciPy fallback is slow). The device VM network is ~6 KiB/s, so computation runs in the cloud. Graphene lattice + Δ ≈ 1 meV needs ξ ≫ a; use a scaled Δ (e.g. 10–50 meV) and check the scaling. Efetov figure data may be image-only |
| **First experiment** | read arXiv:1706.03457 and arXiv:1801.00620 (abstract + model section) and record whether either uses a metal-contact doping profile; then build the k_y-resolved Kwant NS junction with exactly ONE validation: normal incidence (k_y = 0), clean interface, eV < Δ, heavily doped SC → Andreev probability R_A = 1 (Klein tunnelling + perfect AR, so G = 2G_0 for that mode) |
