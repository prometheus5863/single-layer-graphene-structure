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
