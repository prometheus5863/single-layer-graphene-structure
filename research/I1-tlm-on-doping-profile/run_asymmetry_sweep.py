"""I1, experiment 7 (2026-11-25): density sweep of the p/n asymmetry of the TLM junction term for Au,
tested against the text-level constraints of Huard et al. PRB 78, 121402 (2008), arXiv:0804.2040.

Model: unchanged from 11-24 (run_tlm_gate_tuned.py). The contact edge is gate-tuned (Au W - W0 = 0.59 eV,
lambda = 27 nm, Malec 2011 fit), src is unchanged, and the TLM window is 1-10 um with R_int = 0.
For |n| = 0.5 ... 8 x 10^12 cm^-2 on both branches the script computes
    dR(|n|) = R_c,TLM(+|n|) - R_c,TLM(-|n|)    (per contact, Ohm.um)  = 2 R_odd * w in Huard's Eq. (1)-(2)
            (assumes rho_e = rho_h, as Huard measure 1 +/- 0.03, and one invasive interface per R_c).

Text-level data from Huard 2008 (figures not readable here; see RESULT-2026-11-25.md):
  H1  sign: R_odd > 0 (electron side more resistive) <=> p-type contact doping.
  H2  for n > ~1.2e12 cm^-2, R_odd "saturates to a finite value" (Ti/Au and Pd, data up to ~7e12).
  H3  Eq. (3): sharp p-n interface lower bound r > 2.06 (h/4e^2) n^-1/2 (ballistic, n_c -> infinity).
Huard's metal is Ti(5 nm)/Au(25 nm), not clean Au, so the comparison is a shape and sign test, not a fit.

Validation (exact limit): with C_I -> 1e-12 F/m^2 the profile is flat on both branches, so dR = 0 exactly.
"""
import os, sys, csv
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.constants import e, hbar, h
import graphene_contact_doping_model as cdm
from graphene_fet_model import v_F
sig = cdm._sheet_conductivity_puddle_regularized

# copied verbatim from run_tlm_gate_tuned.py (11-24); not imported, because importing re-runs it and rewrites its outputs
def EF(n): return np.sign(n) * hbar * v_F * np.sqrt(np.pi * np.abs(n)) / e
def n_edge(dV, n_g, CI=cdm.C_interface):
    gfun = lambda n: EF(n) + e * (n - n_g) / CI - dV
    span = 1e18 + abs(n_g) + abs(CI * dV / e) * 2
    return brentq(gfun, -span, span, xtol=1e6, rtol=1e-13)
def S(x, L, lam):
    a, b = 1 / (1 + x / lam), 1 / (1 + (L - x) / lam)
    return a + b - a * b
DV_AU, LAM = -0.59, 27e-9
def r_tot(L, n_b, n_c, lam=LAM, npts=4000):
    hh = np.geomspace(1e-13, L / 2, npts); x = np.concatenate([[0.0], hh])
    n = n_b + (n_c - n_b) * S(x, L, lam)
    return 2 * np.trapezoid(1 / sig(n), x) * 1e6
def tlm(n_b, n_c, a=1e-6, b=10e-6, lam=LAM, k=8):
    Ls = np.linspace(a, b, k); R = np.array([r_tot(L, n_b, n_c, lam) for L in Ls])
    p = np.polyfit(Ls * 1e6, R, 1); fit = np.polyval(p, Ls * 1e6)
    return p[1] / 2, 1 - np.sum((R - fit) ** 2) / np.sum((R - R.mean()) ** 2)

T = []; P = T.append
# --- validation
nb = 3e16
d0 = tlm(nb, n_edge(DV_AU, nb, CI=1e-12))[0] - tlm(-nb, n_edge(DV_AU, -nb, CI=1e-12))[0]
P(f"VALIDATION C_I -> 1e-12 F/m^2 at |n| = 3e12 cm^-2: dR = {d0:.2e} Ohm.um (exact 0) -> {'PASS' if abs(d0) < 1e-3 else 'FAIL'}")
ng_star = cdm.C_interface * abs(DV_AU) / e**2 * e   # gate density at which the Au contact edge reaches neutrality (C_I dV / e)
P(f"Au p-n junction closes at n_g* = C_I|dV|/e = {ng_star*1e-4:.2e} cm^-2 (n_c = 0 there exactly)")
P("")
P("|n| (cm^-2) | n_c(p) n_c(n) | R_c,TLM p | R_c,TLM n | dR = n - p | Huard Eq.3 bound | dR / bound | min R^2")
rows = []
for nabs in np.array([0.5, 0.8, 1.2, 1.6, 2, 2.5, 3, 4, 5, 6, 7, 8]) * 1e16:
    ncp, ncn = n_edge(DV_AU, -nabs), n_edge(DV_AU, nabs)
    (rp, r2p), (rn, r2n) = tlm(-nabs, ncp), tlm(nabs, ncn)
    bound = 2.06 * h / (4 * e**2) / np.sqrt(nabs) * 1e6       # Ohm.um
    d = rn - rp
    rows.append(dict(n_abs_cm2=nabs * 1e-4, n_c_p_cm2=ncp * 1e-4, n_c_n_cm2=ncn * 1e-4, Rc_TLM_p=rp, Rc_TLM_n=rn,
                     dR_n_minus_p=d, huard_eq3_bound=bound, ratio_dR_bound=d / bound, min_R2=min(r2p, r2n)))
    P(f"{nabs*1e-4:.1e} | {ncp*1e-4:+.2e} {ncn*1e-4:+.2e} | {rp:+8.1f} | {rn:+8.1f} | {d:+8.1f} | {bound:7.1f} | {d/bound:5.2f} | {min(r2p, r2n):.5f}")
P("")
dR = np.array([r["dR_n_minus_p"] for r in rows]); na = np.array([r["n_abs_cm2"] for r in rows])
sel = (na >= 1.2e12) & (na <= 7e12)
P(f"H1 sign: dR > 0 for all |n| in the sweep: {bool(np.all(dR > 0))}  (Huard: R_odd > 0, p-type contacts)")
P(f"H2 saturation over 1.2-7e12 cm^-2: dR falls from {dR[sel][0]:.0f} to {dR[sel][-1]:.0f} Ohm.um "
  f"(ratio {dR[sel][0]/dR[sel][-1]:.1f}x); a 'saturated' R_odd would keep this ratio near 1")
fit = np.polyfit(np.log(na[sel]), np.log(dR[sel]), 1)
P(f"   power-law slope d ln dR / d ln n over 1.2-7e12 = {fit[0]:.2f} (Eq. 3 bound: -0.50; saturation: 0)")
bnd = np.array([r["ratio_dR_bound"] for r in rows])
P(f"H3 ballistic sharp-junction bound: dR / bound = {bnd.min():.2f}-{bnd.max():.2f}, min at |n| = {na[bnd.argmin()]:.1e} cm^-2; "
  f"{'never below 1 in the sweep' if bnd.min() >= 1 else 'below 1 somewhere: diffusive term under the ballistic minimum'}")
with open(os.path.join(HERE, "asymmetry_sweep_results.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader()
    for r in rows: w.writerow({k: f"{v:.6g}" for k, v in r.items()})
open(os.path.join(HERE, "asymmetry_sweep_transcript.txt"), "w").write("\n".join(T) + "\n")
print("\n".join(T))

fig, ax = plt.subplots(figsize=(6.4, 4.4))
ax.plot(na, dR, "C0o-", label="model ΔR = R_c,TLM(n) − R_c,TLM(p), Au, λ = 27 nm")
ax.plot(na, [r["huard_eq3_bound"] for r in rows], "k--", label="Huard 2008 Eq. 3: sharp p–n ballistic lower bound")
ax.axvspan(0, 1.2e12, color="gray", alpha=.15, label="Huard: R_odd fluctuates (n < 1.2×10¹²)")
ax.axvline(ng_star * 1e-4, c="C3", ls=":", label=f"Au junction closes, n_g* = {ng_star*1e-4:.1e}")
ax.set_xlabel("|n| (cm⁻²)"); ax.set_ylabel("Ω·µm per contact (= 2 R_odd · w)")
ax.set_ylim(0, None); ax.legend(fontsize=7)
ax.set_title("I1: p/n asymmetry vs density; text constraints from Huard 2008 (arXiv:0804.2040)\n"
             "Huard figure points not digitised (images inaccessible); Huard contacts are Ti/Au", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "asymmetry_sweep.png"), dpi=150)
