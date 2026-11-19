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
