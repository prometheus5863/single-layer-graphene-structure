# Research log

One dated entry per daily session: mode, idea, what was done, result in one line, next experiment.

## 2026-11-17 (Tue) — EXECUTE — I2
- **Mode decision:** IDEAS.md uses its own status words (new/checking/running/parked/killed/done) and no row says LIVE. Rows I1, I2 and I4 (not parked) were read as live, which is 3, so today is EXECUTE, not FORGE. Status mapping used from now on: LIVE ≈ new/checking/running, PAPER-READY ≈ done. The forge seeds S1–S4 wait for Sunday 2026-11-22.
- **Idea choice:** none of the three had progress. I2 was picked because its next experiment can be run directly on existing `src/` code.
- **Done:** digitised Giovannetti et al. PRL 101, 026803 (2008) Table I and reran the §4.7 contact-doping model (src unchanged) with facet-matched (111) W. Compared the contact ΔE_F with DFT. Validation: n_contact(W0) = 0 passes.
- **Result:** Au's doping flips n→p (R_extra −161 → +947 Ω·µm). The model matches the DFT doping sign for 5/5 physisorbed metals, with RMS ΔE_F error 0.074 eV. The corrected §4.7 claim "R_extra < 0 for all physisorbed contacts" holds only for Al/Ag once W and W0 come from one source.
- **Next experiment:** I2: second reference (Khomyakov 2009 tabulated ΔE_F) + one experimental Au-contact doping sign, then close I2 and feed the doping into I1/I4.
- **Environment note:** the device VM's /sessions disk was full (0 B free), so the repo was cloned to /tmp/work and scipy/matplotlib were installed with `pip --target /tmp/pylib`. A scratch copy of the figure was left in the connected folder at `_tmp_view/` (deletion is not permitted there).

## 2026-11-18 (Wed) — EXECUTE — I2
- **Mode decision:** not Sunday. Under the status mapping set on 11-17, I1, I2 and I4 count as live (3), so today is EXECUTE. I2 was picked because it has the most progress.
- **Done:** digitised Khomyakov et al. PRB 79, 195425 (2009), Table I, using the printed ΔE_F plus the experimental W_M column. Added one experimental Au doping sign, from STS on polycrystalline Au (Sławińska et al. PRB 85, 235430 (2012), p-type, +0.25 to +0.55 eV). Reran the §4.7 model (src unchanged). Validation: ΔE_F(W0) = 0, which passes.
- **Result:** against the printed DFT values the model gets 5/5 signs, RMS 0.053 eV. Yesterday's W_cov − W_G proxy was itself 0.089 eV off. Au p-type is confirmed by experiment; the model's magnitude of 0.13 eV is below the STS range. With experimental W_M, Au wrongly goes n-type, so W and W0 must come from one source. I2 is closed as **parked (consistency fix done, not PAPER-READY)**. Its doping output feeds I1 and I4.
- **Next experiment:** I1, the TLM fit of total R vs contact spacing on the §4.5 profile, using I2's calc (111) W. The forge on Sunday 2026-11-22 should add seed S1 or S3, because only 2 live ideas remain.
- **Environment note:** /sessions on the device VM is still full (0 B free, shared with other sessions). The repo was cloned to /tmp/hwork instead of $HOME/work, and the /tmp/pylib packages from 11-17 were reused.

## 2026-11-19 (Thu) — EXECUTE — I1
- **Mode decision:** not Sunday. Under the 11-17 status mapping, I1 and I4 are live (2), so today is EXECUTE. I1 was picked because its next experiment was concrete and was set on 11-18. I4 has none yet.
- **Done:** a two-contact TLM simulation on the src doping profile (src unchanged) with calc (111) W for Al, Ag, Cu, Au and Pt. It sweeps the channel length, fits straight lines over three TLM windows, and compares R_c,TLM with the src R_extra. Validation: with n_c = n_bulk, the intercept is 0 and the slope equals R_sheet, which passes.
- **Result:** the kill condition is not hit, because R_c,TLM = 0.48–1.13× R_extra for Al–Au. However, R_extra grows like ln L_j on the 1/(1+x/λ) tail, so the src 1 µm value is an arbitrary cut, and R_c,TLM shifts by ×2–3 with the fit window. The Pt p–n crossing sits at 1.25 µm, outside the src window (+173 Ω·µm at 1 µm vs +4439 at 10 µm). The model predicts visible R_tot(L) curvature for Au, Al and Pt.
- **Next experiment:** I1, digitise a published Au TLM plot (Passi et al. arXiv:1807.04772) and test the predicted curvature against the data scatter. The Sunday 2026-11-22 forge should add S3 (BTK in Kwant) or S1.
- **Environment note:** /sessions is still full (0 B free). /tmp/work from earlier runs belongs to another user, so today's clone is in a fresh mktemp directory under /tmp. The /tmp/pylib packages were reused.

### 2026-11-19 (Thu), continued on Harsh's request — EXECUTE — I1
- **Push:** retried from the cloud workspace. It was refused because the repo is not in that sandbox's authorized repository set, and the device VM blocks git push. The commits stay in `_to_push/2026-11-19/`. The token copy was shredded.
- **Done:** compared the Au TLM simulation with Passi et al. arXiv:1807.04772 (Table 1 R_c, R² > 0.984). No R(L) points were accessible, because the figures, the supplement and arxiv.org were all blocked. Validation: R² = 1 at n_c = n_bulk, which passes.
- **Result:** the comparison is consistent but weak. In a p-channel at −40 V, R(L) stays linear (R² ≥ 0.988), and the junction gives 24–55 % of the measured 519 Ω·µm. The curvature prediction needs n-branch data. The Dirac-point comparison flags the src puddle-σ model.
- **Next experiment:** I1, an n-branch or ambipolar TLM dataset for Au or Pt contacts (Xia 2011 / Nagashio 2009–2010, to be confirmed). The Sunday forge should still add S3.

## 2026-11-20 (Fri) — EXECUTE — I1
- **Mode decision:** not Sunday. Under the 11-17 status mapping, I1 and I4 are live (2), so today is EXECUTE. I1 was picked because it has the most progress; I4 still has no concrete experiment.
- **Done:** found a test of the contact doping tail that the model can fail. Nouchi et al. arXiv:1102.2954 (2011) measure V_NP vs channel length L = 0.4–2.5 µm for Au, Cu and Ag contacts; only the text values could be read. I computed V_NP(L) on the src profile, with src unchanged and calc (111) W. Validation: with n_c = 0 the shift is 0 at every L, which passes. On the way, I1's two-contact superposition was found to go above 1 for L ≲ 2λ, which inverts the gate's leverage mid-channel. It was replaced by a bounded combination.
- **Result:** **0/3 doping signs match** (Au, Cu, Ag). The Au shift of −28 V is reproduced with n_c = +4.5×10¹¹ cm⁻² and λ = 250 nm, so the tail length is plausible and the W → sign rule is what fails here. Nouchi's Cr adhesion layer under Au and the Au-capped Cu/Ag films may explain the reversal, so this is not yet a kill.
- **Next experiment:** I1, digitisable V_NP(L) or R_c(V_g) data for clean contacts (arXiv:1108.3238, clean Au), to settle the sign question. The Sunday 2026-11-22 forge should add S3 (BTK in Kwant).
- **Environment note:** /sessions on the device VM is still full (0 B free), so the clone is in /tmp/hv-work and the /tmp/pylib packages were reused. The subagent fetches of arXiv PDF text and ar5iv worked today. Also read: Huard et al. PRB 78, 121402 (2008), which reports Ti/Au contacts p-doping the channel (R_odd > 0), against Nouchi's Au result. It gives no numbers that can be digitised.

## 2026-11-21 (Sat) — EXECUTE — I1
- **Mode decision:** Saturday, not Sunday. Under the 11-17 status mapping, I1 and I4 are live (2), so today is EXECUTE. I1 was picked because it has the most progress and a concrete next experiment from 11-20.
- **Done:** read the clean-contact dataset named on 11-20, Malec & Davidović PRB 84, 033407 (2011), arXiv:1108.3238. Their Au grains sit directly on SLG with no resist and no adhesion layer. Text values only: p-doping, bulk R max at +85 V, R_c max at +120 V (9.3e12 cm⁻², ΔE_F = 0.35 eV), C_g = 12.4 nF/cm², 35 nm gaps between grains. Ran the 11-20 gap model (src unchanged) with λ as a free parameter. Validation: with n_c = 0, V_NP = 0 at λ = 5/50/250 nm, which passes.
- **Result:** the doping **sign matches** (p), so clean Au now agrees with the model 2/2 and Nouchi's n-type Au points to its Cr underlayer. The contact density is 7× low (model 0.13 eV vs 0.35 eV). **src λ = 250 nm fails:** it pins the 35 nm gap, giving gap neutrality near +29 kV against +85 V measured. The data need λ ≈ 5–15 nm, and with λ that short, Nouchi's −28 V needs an implausible n_c = +4.4e13 cm⁻². The pinned, gate-independent edge in the src profile is the weak point.
- **Next experiment:** I1, a gate-dependent contact profile (C_q-limited, gate-tuned under-contact density; tail screened by the gate per Khomyakov 2010, arXiv:0911.2027), in research/ only. Fit Malec and Nouchi with one parameter set, then rerun the 11-19 TLM numbers. The forge on Sunday 2026-11-22 should add S3 (BTK in Kwant).
- **Environment note:** /sessions on the device VM is still full (0 B free), so the clone is in /tmp/hv-mnhj (shallow, blob-less, sparse: research/ + src/). The device network was too slow to pip-install scipy (37 MB wheel timed out). The calculation therefore ran in the cloud workspace on raw copies of the src files (unchanged), and the outputs were copied to the device via the connected folder's `_tmp_transfer/2026-11-21/` (deletion is not permitted there, so the copy remains).

## 2026-11-22 (Sun) — FORGE — S3 → I5, S1 → I6
- **Mode decision:** Sunday (IST), so today is FORGE. Two candidates were forged from the seed list: S3, because it serves the registered TBG–SC title, and S1, because I1's 11-21 finding (λ ≈ 5–15 nm, shorter than the mean free path) makes the src diffusive R_extra suspect. Adding both would make 4 live ideas (I1, I4, I5, I6), so **I4 was parked as superseded by I5**. I5's step 1 is I4's next experiment. S2 (photo-thermoelectric) and S4 (failure taxonomy) wait for a later forge.
- **Done:** Stages 1–4 for both candidates. There were 4 web searches and 4 pages fetched (2 per candidate). All four fetches worked. arXiv:1605.04254 was abstract only, and arXiv:2211.11389 was about the first 16 % verbatim with the rest summarised. Several prior-art rows are title-only and are marked that way in the FORGE files. Kwant 1.5.0 was confirmed to build in the cloud workspace (Python 3.11 venv via uv, `--no-build-isolation`, no MUMPS). It does not build on Python 3.13, and the device VM cannot install it because its network runs at about 6 KiB/s.
- **Result:** **I5** (Andreev reflection with the metal-contact doping profile replacing a free step U0) is LIVE. Its novelty is unconfirmed: arXiv:1706.03457 ("role of spatial variation of the CNP") is the closest threat and has not been read yet. **I6** (ballistic vs diffusive junction R) is LIVE with a narrow wedge (Xia 2011 is the closest prior art) and may fold into I1.
- **Next experiment:** EXECUTE Monday. I1 keeps the gate-dependent profile (the most progress). I5's first experiment is the 1706.03457/1801.00620 prior-art read followed by the Kwant NS validation (normal-incidence R_A = 1).
- **Environment note:** /sessions on the device VM is still full (0 B free), and the device network is about 6 KiB/s, so a git clone there did not finish. The repo was cloned and edited in the cloud workspace. The commits were moved to the device as a git bundle through the connected folder (`_tmp_transfer/2026-11-22/`). **The push from the device was refused** (Cowork HTTP 403: only clone and fetch are allowed, as on 11-19), so the commits are in `_to_push/2026-11-22/` (bundle, format-patch and README).

## 2026-11-23 (Mon) — EXECUTE — I1
- **Mode decision:** Monday (IST), and I1, I5 and I6 are live (3), so today is EXECUTE. I1 was picked because it has the most progress and an unblocked next experiment set on 11-21. I5 and I6 have 0 days.
- **Done:** replaced src's pinned contact edge, in research/ only (src unchanged), with a gate-tuned edge. The back-gate charge enters src's own C_q-in-series-with-C_I contact equation, E_F(n)/e + e(n − n_g)/C_I = −(W_M − W0)/e. Fitted Malec 2011 clean Au: the +120 V R_c maximum fixes W − W0, and the +85 V gap maximum fixes λ. Then applied the same parameters to Nouchi 2011 (Au, −28 V). Validation: C_I → 0 gives n_c = n_g (residual 5×10² cm⁻²), which passes after a metric bug was fixed (relative error divided by 1 m⁻² at n_g = 0).
- **Result:** the gate-tuned edge fits clean Au with **W − W0 = 0.59 eV and λ = 27 nm**. src's 0.14 eV puts the R_c maximum at +29 V instead of +120 V, and src's λ = 250 nm now gives +116 V instead of +29 kV. Out of sample, Nouchi's Au magnitude is −15.6 V against −28 V measured (11-21: 4.7× off in n_c); the sign is still opposite (Cr underlayer). Side finding: src's linearised n_c overestimates the exact series solution by 1.8×.
- **Next experiment:** I1, rerun the 11-19 TLM extraction on the gate-tuned edge (λ = 27 nm, Au W − W0 = 0.59 eV) at −40 V and on the n-branch, to retest the kill condition.
- **Environment note:** /sessions on the device VM is still full (0 B free). A git clone into /tmp on the device stalled (no data after 150 s; codeload timed out). The repo was cloned and run in the cloud workspace with src unchanged.

## 2026-11-24 (Tue) — EXECUTE — I1
- **Mode decision:** Tuesday (IST), and I1, I5 and I6 are live (3), so today is EXECUTE. I1 was picked because it has the most progress (4 days) and its next experiment, set on 11-23, was unblocked. I5 and I6 are still at 0 days.
- **Done:** reran the 11-19 TLM extraction on the 11-23 gate-tuned contact edge (Au W − W0 = 0.59 eV, λ = 27 nm). The p-branch at Passi's −40 V was bracketed with n = −2 and −5 ×10¹² cm⁻², and the n-branch used +2 and +5 ×10¹². Three fit windows were used, with src unchanged. The 11-23 functions were copied, not imported, so its outputs are not regenerated. Validation: with C_I → 0 (metal decoupled), R_c,TLM = 1.7×10⁻⁸ Ω·µm (exact value 0) and the slope equals R_sheet. It passes.
- **Result:** the kill condition is **not hit**: R_c,TLM = 0.99–1.11× R_extra, with R² ≥ 0.9996. The 11-19 window and curvature artefact was caused by λ = 250 nm and is withdrawn. What survives is the junction p/n asymmetry, −112 vs +306 Ω·µm at |n| = 2×10¹² cm⁻². On this model, Passi's −40 V R_c is all interface resistance.
- **Next experiment:** I1, digitise a published ambipolar Au TLM R_c(V_g) and compare R_c(n) − R_c(p) with the model's ≈ +420 Ω·µm.
- **Environment note:** /sessions on the device VM is full (0 B free, used by other data). A device clone into /tmp ran at about 10 KB/s and was killed with its shell. The work was done in the cloud workspace; the push goes through a small bundle applied to a blobless shallow clone on the device.

## 2026-11-25 (Wed) — EXECUTE — I1
- **Mode decision:** it is Wednesday (IST) and three ideas are live (I1, I5, I6), so today is EXECUTE. I1 was picked because it has the most progress (5 days) and an unblocked next experiment set on 11-24.
- **Done:** searched for a published ambipolar Au TLM R_c(V_g). Two searches, three fetches.
  - Huard et al., PRB 78, 121402 (2008), arXiv:0804.2040, measures exactly the odd part R(n) − R(−n) at Ti/Au and Pd contacts. Only the text was readable. arxiv.org PDFs are proxy-blocked (403) from both the cloud workspace and the device VM, so the figures were not digitised, and `huard2008_text_constraints.csv` holds text values only.
  - Nagashio arXiv:1004.2074 uses Cr/Au and gives per-device R_c only, so it is not usable.
  - Ran a density sweep of ΔR = R_c,TLM(n) − R_c,TLM(p) on the 11-24 gate-tuned Au edge (src unchanged). Validation: C_I → 0 gives ΔR = 1.3×10⁻⁸ Ω·µm (exact value 0). Pass.
- **Result:** the sign matches Huard (electron side more resistive). The density dependence does not. The model falls 11.7× over 1.2–7×10¹² cm⁻² (slope −1.42) because the gate closes the Au p–n junction at 9.3×10¹², whereas Huard report a saturated R_odd there. Caveat: Huard's metal is Ti/Au, and the comparison uses text only. At 8×10¹² the diffusive ΔR reaches the ballistic sharp-junction bound (1.05×), which supports I6.
- **Next experiment:** I6. Compute the ballistic (Klein-filtered, k_y-resolved) junction R on the same gate-tuned edge and test it for saturation over 1.2–7×10¹² cm⁻² (validation: flat profile gives T = 1, the Sharvin limit). I1 resumes once the Huard PDF is available for digitising.
- **Environment note:** /sessions on the device VM is full (0 B free, used by other data). A device clone into /tmp ran at about 15 KB/s and timed out. The work was done in the cloud workspace, and the push goes through a small bundle applied to a blobless shallow clone on the device, as on 11-24.

## 2026-11-26 (Thu) — EXECUTE — I6
- **Mode decision:** it is Thursday (IST) and three ideas are live (I1, I5, I6), so today is EXECUTE. I1 is blocked: it needs the Huard 2008 PDF, which is proxy-blocked. Its own 11-25 entry named I6 as the next experiment, so I6 was picked. I5's first step is a prior-art read with 0 days of progress, so it waits.
- **Done:** wrote a k_y-resolved Dirac transfer-matrix (Landauer) calculation through I1's gate-tuned Au contact edge (W − W0 = 0.59 eV; src unchanged, imported) at λ = 5, 15, 27 and 250 nm. Both branches were run over |n| = 0.5–8 ×10¹² cm⁻², alongside src's diffusive R_extra on the same tail. A first forward transfer-matrix version overflowed on the p–n branch and was replaced, before any number was used, by backward propagation of the transmitted state with log-norm renormalisation. Validation: a flat profile gives T = 1 and R = the Sharvin value (relative error 1e-13). PASS.
- **Result:** the transport model matters.
  - On the p–n (n-) branch, R_J/R_extra ranges from 1.4–10× at λ = 5 nm to 0.17–0.62× at 250 nm.
  - On the p–p′ branch, ballistic R_J ≈ 0 while diffusive R_extra < 0, a sign Landauer cannot give.
  - The ballistic p/n asymmetry falls with slope −0.4 to −0.6, against −1.5 diffusive. That is closer to Huard's "saturated" R_odd but not saturated.
  - Neither junction model gives Passi's 519 Ω·µm at −40 V, which points to metal injection.
  - The kill condition was not hit, and success-metric part 1 was met.
- **Next experiment:** I6 day 2. Add R_mg(n_c) using the Chaves 2015 form, with one parameter fixed by Passi, to get a full ambipolar R_c(V_g). Compare it with Malec 2011 and Huard with one parameter set, and read Xia 2011 for prior art first.
- **Environment note:** /sessions on the device VM is still full (0 B free, other data). A device clone into /tmp/hw-work worked quickly today, but the device has no scipy. The calculation ran in the cloud workspace on a fresh clone with src unchanged.

## 2026-11-27 (Fri) — EXECUTE — I6
- **Mode decision:** it is Friday (IST) and three ideas are live (I1, I5, I6), so today is EXECUTE.
  - I1 is still blocked on the Huard PDF.
  - I6 has the most recent progress and an unblocked next step, so it was picked.
  - I5 waits, at 0 days.
- **Done:**
  - Added metal injection R_mg as a transmission-line ladder (ρ_c with the src sheet resistance under the metal at n_c). Then R_c = R_mg + ballistic R_J on the gate-tuned Au edge.
  - Fixed the one parameter ρ_c to Passi 2018 at −40 V (519 Ω·µm), then predicted Passi's Dirac-point R_c and Malec 2011's gate ratio with no refit.
  - Read the Xia 2011 abstract.
  - Validation: the uniform ladder equals √(ρ_c R_s)coth(L_c/L_T) to 2×10⁻⁸. PASS.
- **Result:** with one parameter, the predictions split.
  - Malec R_c(0)/R_c(−100 V) comes out at 1.51 at λ = 5 nm (1.76 at 15 nm, 1.93 at 27 nm) against the measured 1.35.
  - Passi CNP/−40 V comes out at 1.80–2.09 against the measured 2.64.
  - Passi's p-side R_c is all injection (R_J < 0.1 Ω·µm).
  - ρ_c is degenerate with the unknown contact length (L_T > L_c).
  - The Xia abstract is Pd only, so the kill condition is not hit.
- **Next experiment:** I6 day 3. Run the diffusive head-to-head: the same fit with src R_extra in place of R_J, and both ratios ballistic vs diffusive at λ = 5/15/27 nm.
- **Environment note:**
  - /sessions on the device VM is full (0 B free, other sessions' data) and /tmp there is not writable, so the clone and calculation ran in the cloud workspace with src unchanged.
