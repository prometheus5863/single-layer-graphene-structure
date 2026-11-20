"""I1, experiment 3 (2026-11-20): channel-length dependence of the charge-neutrality gate
voltage V_NP(L) predicted by the src contact-doping profile (src/graphene_contact_doping_model.py,
unchanged), compared with Nouchi, Saito, Tanigaki arXiv:1102.2954 (2011): Au, Cu, Ag contacts,
L = 0.4-2.5 um, 300 nm SiO2.

Device model: run_tlm_extraction.py superposes the two contacts as S = f(x) + f(L-x),
f = 1/(1+x/lambda), lambda = src 250 nm. For L < ~2 lambda that S exceeds 1 mid-channel, and then
the gate's leverage there, n = n_b (1-S) + n_c S, changes sign. That is unphysical and gave spurious
V_NP shifts in a first pass today. The main model here is therefore the bounded combination
    S_or = f(x) + f(L-x) - f(x) f(L-x)   (0 <= S_or <= 1, S_or = 1 exactly at each contact),
    n(x) = n_b + (n_c - n_b) S_or(x).
The plain sum is kept as RULE="sum" for comparison only. The choice of S_or is an ansatz;
src gives a single-contact profile only.
    R(L, n_b) = int_0^L dx / sigma(n(x)),   sigma = src puddle-regularised Drude
n_b = C_ox (V_g - V_0)/e (geometric capacitance only). V_NP(L) = gate voltage of the maximum
two-terminal resistance (Nouchi plot two-terminal conductivity). Shifts are quoted relative to
L = 2.5 um, the longest Nouchi device, because the background V_0 (-43 V for Au) is substrate doping
the model does not contain.
Validation (exact limit): n_c = 0 -> n(x) = n_b (1 - f - f) vanishes everywhere at n_b = 0, so
V_NP(L) = V_0 for every L, shift exactly 0.
"""
import os, sys, csv
import numpy as np
from scipy.optimize import minimize_scalar, brentq
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import graphene_contact_doping_model as cdm

lam = cdm.lambda_decay
sig = cdm._sheet_conductivity_puddle_regularized
f = lambda x: 1.0 / (1.0 + x / lam)
e = 1.602176634e-19; eps0 = 8.8541878128e-12
C_OX = 3.9 * eps0 / 300e-9            # F/m^2, Nouchi 300 nm SiO2
N_PER_V = C_OX / e                     # m^-2 per volt

def S(x, L, rule):
    a, b = f(x), f(L - x)
    return a + b - a * b if rule == "or" else a + b

def r_tot(L, n_c, n_b, npts=3000, rule="or"):
    h = np.geomspace(1e-12, L / 2, npts); x = np.concatenate([[0.0], h])
    n = n_b + (n_c - n_b) * S(x, L, rule)
    return 2 * np.trapezoid(1.0 / sig(n), x)

def v_np(L, n_c, rule="or"):
    """gate voltage (relative to V_0) of max R; bounded search over +-300 V."""
    res = minimize_scalar(lambda v: -r_tot(L, n_c, v * N_PER_V, rule=rule), bounds=(-300, 300),
                          method="bounded", options=dict(xatol=1e-4))
    return res.x

def shift(L, n_c, Lref=2.5e-6, rule="or"): return v_np(L, n_c, rule) - v_np(Lref, n_c, rule)

T = []; P = T.append
# --- the one validation
v0 = [v_np(L, 0.0) for L in (0.4e-6, 1e-6, 2.5e-6)]
ok = max(abs(v) for v in v0) < 1e-3
P(f"VALIDATION n_c = 0: V_NP - V_0 at L = 0.4/1/2.5 um = {', '.join(f'{v:+.1e}' for v in v0)} V (exact 0) -> {'PASS' if ok else 'FAIL'}")
P(f"C_ox = {C_OX*1e5:.3f} nF/cm^2 ; {N_PER_V*1e-4:.3e} cm^-2 per V ; lambda = {lam*1e9:.0f} nm")
P("")
K = {r["metal"]: r for r in csv.DictReader(l for l in open(os.path.join(HERE, "..", "I2-au-facet-work-function", "khomyakov2009_table1.csv")) if not l.startswith("#"))}
data_sign = {"Au": "neg", "Cu": "pos", "Ag": "pos"}
Lg = np.array([0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5]) * 1e-6
rows = []
P("metal | W111 calc (eV) | n_c model (cm^-2) | model doping | dV_NP(0.4 um) vs 2.5 um, S_or (V) | [plain-sum rule, artefact] | model sign | Nouchi sign | agree")
for m in ("Au", "Cu", "Ag"):
    W = float(K[m]["W_M_111_calc_eV"])
    n_c = cdm.contact_edge_carrier_density(W)
    sh = [shift(L, n_c) for L in Lg]
    msign = "neg" if sh[0] < 0 else "pos"
    agree = msign == data_sign[m]
    P(f"{m} | {W:.2f} | {n_c*1e-4:+.2e} | {'n' if n_c > 0 else 'p'} | {sh[0]:+7.2f} | [{shift(0.4e-6, n_c, rule='sum'):+7.2f}] | {msign} | {data_sign[m]} | {'yes' if agree else 'NO'}")
    for L, s in zip(Lg, sh): rows.append(dict(metal=m, source="model_calc111W", n_c_cm2=n_c * 1e-4, L_um=L * 1e6, dV_NP_V=s))
P("")
# --- inverse: what contact density reproduces Nouchi's Au -28 V between L = 2.5 and 0.4 um?
target = -71 - (-43)
g = lambda nc: shift(0.4e-6, nc) - target
nc_fit = brentq(g, 1e14, 1e17, xtol=1e12)  # upper end stays inside the +-300 V search
# equivalent W in the src model (monotone n_c(W))
W_fit = brentq(lambda W: cdm.contact_edge_carrier_density(W) - nc_fit, 3.0, 5.4)
nc_au = cdm.contact_edge_carrier_density(float(K["Au"]["W_M_111_calc_eV"]))
P(f"INVERSE Au: n_c reproducing dV_NP(0.4 vs 2.5 um) = {target} V : n_c = {nc_fit*1e-4:+.2e} cm^-2 (n-type), "
  f"src-equivalent W = {W_fit:.2f} eV (calc 111 W = {K['Au']['W_M_111_calc_eV']} eV, model n_c = {nc_au*1e-4:+.2e} cm^-2)")
sh_fit = [shift(L, nc_fit) for L in Lg]
for L, s in zip(Lg, sh_fit): rows.append(dict(metal="Au", source="model_fit_to_Nouchi", n_c_cm2=nc_fit * 1e-4, L_um=L * 1e6, dV_NP_V=s))
P("Au fitted-n_c shape dV_NP(L) [V]: " + ", ".join(f"{L*1e6:.1f}um {s:+.1f}" for L, s in zip(Lg, sh_fit)))
nsig = sum(1 for m in ("Au", "Cu", "Ag") if (shift(0.4e-6, cdm.contact_edge_carrier_density(float(K[m]['W_M_111_calc_eV']))) < 0) == (data_sign[m] == "neg"))
P(f"SIGN SCORE (model calc-111 W vs Nouchi 2011): {nsig}/3")
with open(os.path.join(HERE, "nouchi_vnp_results.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
open(os.path.join(HERE, "nouchi_vnp_transcript.txt"), "w").write("\n".join(T) + "\n")
print("\n".join(T))

fig, ax = plt.subplots(figsize=(6.2, 4.4))
Lf = np.linspace(0.4e-6, 2.5e-6, 30)
for i, m in enumerate(("Au", "Cu", "Ag")):
    n_c = cdm.contact_edge_carrier_density(float(K[m]["W_M_111_calc_eV"]))
    ax.plot(Lf * 1e6, [shift(L, n_c) for L in Lf], f"C{i}-", label=f"model {m} (calc 111 W, n_c={n_c*1e-4:+.1e} cm⁻²)")
ax.plot(Lf * 1e6, [shift(L, nc_fit) for L in Lf], "C0--", label=f"model Au, n_c fitted = {nc_fit*1e-4:+.1e} cm⁻²")
ax.plot([2.5, 0.4], [0, target], "ks", ms=7, label="Nouchi 2011 Au (text endpoints)")
ax.annotate("Nouchi: Cu, Ag shift +\n(figure only)", xy=(0.55, 4), fontsize=8)
ax.axhline(0, c="gray", lw=.6); ax.set_xlabel("channel length L (µm)")
ax.set_ylabel(r"$V_{NP}(L) - V_{NP}(2.5\,\mu m)$ (V)")
ax.set_title("I1: CNP shift vs L, src doping profile vs Nouchi 2011 (arXiv:1102.2954)", fontsize=9)
ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(HERE, "nouchi_vnp_shift.png"), dpi=150)
