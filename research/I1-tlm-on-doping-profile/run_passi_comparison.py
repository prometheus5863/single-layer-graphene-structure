"""I1, experiment 2 (2026-11-19, continued on Harsh's request): compare the I1 TLM
simulation for Au with the only published Au TLM numbers accessible today, Passi et al.
arXiv:1807.04772 (Table 1 R_C, and "R^2 > 0.984 for all TLM fits"). No R_tot(L) points
could be digitised (figure images and Supplementary Fig. S3 not reachable), so this is a
weak test: (a) does the model's predicted curvature for Au keep a linear fit at R^2 > 0.984,
and (b) how much of the measured R_C can the junction term explain, at the two bias points
Passi used (Dirac point, and V_bg = -40 V = hole-doped channel, same type as the p-doped
Au contact, so no p-n junction)?
The channel densities are not given in the paper; the model is run at n_bulk in
{+2e16 (n, I1 baseline), ~0 (Dirac: puddle-limited), -2e16, -5e16 m^-2 (p)}; the
V_bg = -40 V case is bracketed by the two p values (78 V overdrive on a typical 90-300 nm SiO2).
Imports src/ unchanged. Validation (exact limit): R^2 = 1 when n_c = n_bulk.
"""
import os, sys, csv
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import run_tlm_extraction as T  # reuses r_tot / tlm_fit (re-runs and rewrites experiment 1 outputs identically)
cdm = T.cdm

def fit(n_c, n_b, Lmin, Lmax, k=8):
    Ls = np.linspace(Lmin, Lmax, k); R = np.array([T.r_tot(L, n_c, n_b) for L in Ls])
    p = np.polyfit(Ls * 1e6, R, 1); res = R - np.polyval(p, Ls * 1e6)
    r2 = 1 - np.sum(res**2) / np.sum((R - R.mean())**2)
    return p[1] / 2, p[0], r2

W_AU = 5.54  # calc (111), Khomyakov 2009, validated in I2
n_c = cdm.contact_edge_carrier_density(W_AU)
dEF = -np.sign(n_c) * cdm.hbar * cdm.v_F * np.sqrt(np.pi * abs(n_c)) / cdm.e if hasattr(cdm, "v_F") else None
L = []; P = L.append
_, _, r2v = fit(2e16, 2e16, 1e-6, 10e-6)
P(f"VALIDATION n_c = n_bulk: R^2 = {r2v:.12f} (exact 1) -> {'PASS' if abs(1 - r2v) < 1e-12 else 'FAIL'}")
P(f"Au calc (111) W = {W_AU} eV: n_c = {n_c:+.3e} m^-2 (p-type); Passi's DFT surface value E_D-E_F = +0.14 eV vs model +0.13 eV (I2)")
P("")
P("n_bulk (m^-2) | window | R_c,TLM junction (Ohm.um) | R^2 | R^2 > 0.984? | implied R_int to match Passi (Dirac 1372 / -40V 519)")
rows = []
cases = [(+2e16, "n, I1 baseline"), (+1e14, "Dirac (puddle-limited)"), (-2e16, "p (-40 V, low)"), (-5e16, "p (-40 V, high)")]
for nb, lab in cases:
    target = 1372 if "Dirac" in lab else (519 if "p (" in lab else float("nan"))
    for a, b in [(1e-6, 10e-6), (2e-6, 20e-6), (0.5e-6, 5e-6)]:
        rc, s, r2 = fit(n_c, nb, a, b)
        rows.append(dict(n_bulk=nb, case=lab, window_um=f"{a*1e6:g}-{b*1e6:g}", Rc_TLM_junction=rc, slope=s, R2=r2,
                         R2_ok=r2 > 0.984, Passi_Rc=target, implied_R_int=target - rc))
        P(f"{nb:+.0e} {lab:24s} | {a*1e6:g}-{b*1e6:g} um | {rc:+8.1f} | {r2:.4f} | {r2>0.984} | {target - rc:+.0f}")
P("")
au_n = [r for r in rows if r["n_bulk"] == 2e16]
P("Au, n-type channel (I1 baseline): curvature visible to R^2 in a linear fit: min R^2 = %.4f" % min(r["R2"] for r in au_n))
p_rows = [r for r in rows if r["n_bulk"] < 0]
P("Au, p-type channel (Passi -40 V): junction R_c,TLM in [%.0f, %.0f] Ohm.um -> fraction of measured 519: [%.2f, %.2f]"
  % (min(r["Rc_TLM_junction"] for r in p_rows), max(r["Rc_TLM_junction"] for r in p_rows),
     min(r["Rc_TLM_junction"] for r in p_rows) / 519, max(r["Rc_TLM_junction"] for r in p_rows) / 519))
with open(os.path.join(HERE, "passi_comparison_results.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
open(os.path.join(HERE, "passi_comparison_transcript.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
