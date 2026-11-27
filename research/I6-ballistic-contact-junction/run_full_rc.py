"""I6, experiment 2 (2026-11-27): full ambipolar contact resistance R_c(V_g) = R_mg(n_c) + R_J on I1's gate-tuned
Au edge, with ONE parameter (specific contact resistivity rho_c) fixed by Passi 2018 at -40 V, then compared, with no
further fitting, to (a) Passi's Dirac-point R_c, (b) Malec 2011's R_c(0 V)/R_c(-100 V) ratio and R_c peak position.

R_mg: metal -> graphene injection over the contact length L_c, transmission-line model (Berger 1972; Schroder,
Semiconductor Material and Device Characterization, Ch. 3): sheet resistance under the metal R_s = 1/sigma(n_c) (src's
puddle-regularised sigma), interface resistivity rho_c,  R_mg = sqrt(rho_c R_s) coth(L_c / L_T),  L_T = sqrt(rho_c/R_s).
Solved here numerically as a resistor ladder (so it can later take a non-uniform n under the metal); the closed form is
the validation.
R_J: the 11-26 ballistic (k_y-resolved Dirac transfer-matrix) junction excess over the channel Sharvin term, imported
from run_ballistic_junction.py's functions (that script is executed once on import; its outputs are only re-written).

Validation (exactly ONE): uniform ladder vs closed-form sqrt(rho_c R_s) coth(L_c/L_T).
src is imported and not edited.
"""
import os, sys, csv, io, contextlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.constants import e

# import 11-26 functions; suppress its re-run printout. It rewrites its own outputs with identical content (deterministic).
with contextlib.redirect_stdout(io.StringIO()):
    import run_ballistic_junction as bj
sig, n_edge, DV_AU, landauer, sharvin = bj.sig, bj.n_edge, bj.DV_AU, bj.landauer, bj.sharvin
n_puddle = bj.n_puddle

def r_mg_ladder(rho_c, Rs, Lc, N=4000):
    """Contact of length Lc (m), width 1 m: sheet segments Rs*dx in series, interface conductance dx/rho_c to the metal
    (grounded) at each node; current injected at the contact front edge x = 0, back edge x = Lc open. Returns Ohm.m."""
    dx = Lc / N
    g_sh = 1 / (Rs * dx); g_if = dx / rho_c
    # tridiagonal nodal analysis, nodes 0..N, interface conductance halved at the two end nodes (trapezoid)
    gi = np.full(N + 1, g_if); gi[0] *= 0.5; gi[-1] *= 0.5
    diag = gi.copy(); diag[:-1] += g_sh; diag[1:] += g_sh
    off = -g_sh * np.ones(N)
    from scipy.linalg import solve_banded
    ab = np.zeros((3, N + 1)); ab[0, 1:] = off; ab[1] = diag; ab[2, :-1] = off
    I = np.zeros(N + 1); I[0] = 1.0
    V = solve_banded((1, 1), ab, I)
    return V[0]

def r_mg_closed(rho_c, Rs, Lc):
    LT = np.sqrt(rho_c / Rs); return np.sqrt(rho_c * Rs) / np.tanh(Lc / LT)

T = []; P = T.append
# ---- validation
rho_v, Rs_v, Lc_v = 1e-11, 1 / sig(-5e16), 0.3e-6     # 1e-7 Ohm cm^2, Rs at -5e12, Lc ~ L_T
num, ex = r_mg_ladder(rho_v, Rs_v, Lc_v), r_mg_closed(rho_v, Rs_v, Lc_v)
P(f"VALIDATION uniform TLM ladder vs closed form sqrt(rho_c R_s) coth(L_c/L_T): {num*1e6:.5f} vs {ex*1e6:.5f} Ohm.um "
  f"(rel. err {abs(num/ex-1):.1e}; L_T = {np.sqrt(rho_v/Rs_v)*1e9:.0f} nm) -> {'PASS' if abs(num/ex-1) < 1e-4 else 'FAIL'}")

LAM = 27e-9            # I1 11-23 Malec fit; 5 and 15 nm reported as sensitivity
LCS = [0.5e-6, 1e-6, 10e-6]
def rc_parts(ng, rho_c, lam=LAM, Lc=1e-6):
    nc = n_edge(DV_AU, ng)
    Rs = 1 / sig(nc)
    Rmg = r_mg_ladder(rho_c, Rs, Lc) * 1e6
    RL, _ = landauer(nc, ng, lam); RJ = RL - sharvin(ng)
    return nc, Rs, Rmg, RJ

# ---- fix rho_c from Passi -40 V (519 Ohm.um), at both ends of I1's density bracket
PASSI_M40, PASSI_DIRAC = 519.0, 1372.0
P(""); P("Step 1: rho_c fixed so that R_c = R_mg + R_J = 519 Ohm.um (Passi 2018, Au, -40 V)")
fits = {}
for Lc in LCS:
    for ngp in (-2e16, -5e16):
        f = lambda lr: rc_parts(ngp, 10**lr, Lc=Lc)[2] + rc_parts(ngp, 10**lr, Lc=Lc)[3] - PASSI_M40
        lr = brentq(f, -14, -6, xtol=1e-6); fits[(Lc, ngp)] = 10**lr
        nc, Rs, Rmg, RJ = rc_parts(ngp, 10**lr, Lc=Lc)
        P(f"  L_c = {Lc*1e6:4.1f} um, n_g = {ngp*1e-4:+.0e}: rho_c = {10**lr*1e4:.3e} Ohm.cm^2, L_T = {np.sqrt(10**lr/Rs)*1e9:.0f} nm, "
          f"n_c = {nc*1e-4:+.2e}, R_s(n_c) = {Rs:.0f} Ohm/sq, R_mg = {Rmg:.1f}, R_J = {RJ:.2f} Ohm.um")

# ---- prediction 1: Passi Dirac point, channel at the puddle density on either side
P(""); P(f"Step 2 (prediction, no refit): Passi Dirac point R_c (measured {PASSI_DIRAC:.0f}; ratio to -40 V = {PASSI_DIRAC/PASSI_M40:.2f})")
pred_dirac = []
for (Lc, ngp), rho in fits.items():
    vals = []
    for ng0 in (-n_puddle, +n_puddle):
        nc, Rs, Rmg, RJ = rc_parts(ng0, rho, Lc=Lc); vals.append(Rmg + RJ)
    rc0 = np.mean(vals); pred_dirac.append(rc0 / PASSI_M40)
    P(f"  L_c = {Lc*1e6:4.1f} um, rho_c from n_g = {ngp*1e-4:+.0e}: R_c(CNP, n_g = -/+ n_puddle) = {vals[0]:.0f} / {vals[1]:.0f}, mean {rc0:.0f} Ohm.um "
      f"-> ratio {rc0/PASSI_M40:.2f} (measured {PASSI_DIRAC/PASSI_M40:.2f})")

# ---- gate sweep on Malec's oxide (C_g = 12.4 nF/cm^2), Malec V_g = 0 assumed neutral (see malec2011_clean_au.csv)
CG = 12.4e-9 * 1e4      # F/m^2
VG = np.concatenate([np.arange(-100, 131, 10.0)])
rho_mid = fits[(1e-6, -2e16)]
rows = []
P(""); P(f"Step 3 (prediction): R_c(V_g) on Malec's gate (C_g = 12.4 nF/cm^2), L_c = 1 um, rho_c = {rho_mid*1e4:.3e} Ohm.cm^2")
P("  V_g | n_g | n_c | R_mg | R_J | R_c (Ohm.um), for lambda = 5 / 15 / 27 nm")
for vg in VG:
    ng = CG * vg / e
    if vg != 0 and abs(ng) < n_puddle: ng = np.sign(vg) * n_puddle
    rec = dict(Vg=vg, n_g_cm2=ng * 1e-4)
    # at V_g = 0 the channel is at the CNP: average the two puddle-density sides, as in Step 2
    ngs = (-n_puddle, +n_puddle) if vg == 0 else (ng,)
    if vg == 0: rec["n_g_cm2"] = 0.0
    for lam in (5e-9, 15e-9, 27e-9):
        parts = [rc_parts(g, rho_mid, lam=lam) for g in ngs]
        nc = np.mean([p[0] for p in parts]); Rmg = np.mean([p[2] for p in parts]); RJ = np.mean([p[3] for p in parts])
        rec[f"n_c_cm2"] = nc * 1e-4; rec["R_mg"] = Rmg
        rec[f"R_J_{lam*1e9:.0f}"] = RJ; rec[f"R_c_{lam*1e9:.0f}"] = Rmg + RJ
    rows.append(rec)
    P(f"  {vg:5.0f} | {ng*1e-4:+.2e} | {rec['n_c_cm2']:+.2e} | {rec['R_mg']:7.1f} | {rec['R_J_5']:6.1f}/{rec['R_J_15']:6.1f}/{rec['R_J_27']:6.1f} | "
      f"{rec['R_c_5']:7.1f}/{rec['R_c_15']:7.1f}/{rec['R_c_27']:7.1f}")
def at(v, k): return [r[k] for r in rows if r["Vg"] == v][0]
P("")
for lam in (5, 15, 27):
    k = f"R_c_{lam}"; vpk = rows[int(np.argmax([r[k] for r in rows]))]["Vg"]
    P(f"  lambda = {lam:2d} nm: R_c(0)/R_c(-100) = {at(0,k)/at(-100,k):.2f} (Malec 128/95 = {128/95:.2f}); peak at V_g = {vpk:.0f} V (Malec ~120 V); "
      f"R_c(0) = {at(0,k):.0f} (Malec <= 128), R_c(-100) = {at(-100,k):.0f} (Malec 95) Ohm.um")
P("  NB: rho_c is Passi's (different devices, Ti-free Au on SiO2 in ambient); Malec's grains are a different interface, so only"
  " ratios and the peak position are compared, not absolute values.")
P("")
P("Xia et al. Nat. Nanotechnol. 6, 179 (2011), doi:10.1038/nnano.2011.6, abstract read 2026-11-27 (full text not fetched):"
  " Pd only; R_c attributed to transport in graphene under the metal (coupling length vs mean free path), ~75 % transfer"
  " efficiency at low T, 110+/-20 Ohm.um at 6 K. Abstract names no physisorbed metal and no p-n junction term.")

with open(os.path.join(HERE, "full_rc_results.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader()
    for r in rows: w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
open(os.path.join(HERE, "full_rc_transcript.txt"), "w").write("\n".join(T) + "\n")
print("\n".join(T))

fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
vg = np.array([r["Vg"] for r in rows])
axs[0].plot(vg, [r["R_mg"] for r in rows], "k-", label="R_mg (injection, TLM)")
for i, lam in enumerate((5, 15, 27)):
    axs[0].plot(vg, [r[f"R_c_{lam}"] for r in rows], f"C{i}o-", ms=3, label=f"R_c = R_mg + R_J, λ = {lam} nm")
axs[0].plot([0, -100], [128, 95], "rs", ms=7, label="Malec 2011 (text; different interface)")
axs[0].axvline(120, color="r", ls=":", label="Malec R_c max ~120 V")
axs[0].set_xlabel("V_g (V), Malec gate C_g = 12.4 nF/cm²"); axs[0].set_ylabel("R_c (Ω·µm)"); axs[0].set_yscale("log")
axs[0].legend(fontsize=6); axs[0].set_title("Predicted ambipolar R_c(V_g), ρ_c fixed by Passi −40 V", fontsize=9)
labs = [f"L_c={k[0]*1e6:g}µm\nn={k[1]*1e-4:.0e}" for k in fits]
axs[1].bar(range(len(pred_dirac)), pred_dirac, color="C0", label="model R_c(CNP)/R_c(−40 V)")
axs[1].axhline(PASSI_DIRAC / PASSI_M40, color="r", label="Passi 2018 measured 1372/519")
axs[1].set_xticks(range(len(labs))); axs[1].set_xticklabels(labs, fontsize=6)
axs[1].set_ylabel("ratio"); axs[1].legend(fontsize=7); axs[1].set_title("Prediction vs Passi Dirac-point R_c", fontsize=9)
fig.suptitle("I6 2026-11-27: full R_c = R_mg + R_J on the gate-tuned Au edge (one parameter)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "full_rc.png"), dpi=150)
