"""I1, experiment 4 (2026-11-21): clean (adhesion-free) Au contacts, Malec & Davidovic,
PRB 84, 033407 (2011), arXiv:1108.3238 -- the dataset named as the next experiment on 11-20.

Question: does the src contact-doping model (src/graphene_contact_doping_model.py, unchanged)
get the doping SIGN and size right for a clean Au contact, and is its tail length lambda = 250 nm
compatible with the observation that the 35 nm graphene gaps between Au grains are still
gate-tuned through neutrality at V_g = +85 V?

Device model for one gap (same construction as run_nouchi_vnp_shift.py, but lambda is a parameter
here instead of src's fixed value, and the gate capacitance is Malec's measured C_g):
    f(x) = 1/(1+x/lam),  S_or = f(x) + f(L-x) - f(x) f(L-x),
    n(x) = n_b + (n_c - n_b) S_or(x),   n_b = C_g V_g / e,
    R(V_g) = int_0^L dx / sigma(n(x)),  sigma = src puddle-regularised Drude.
V_NP(gap) = gate voltage of maximum R. n_c is pinned (gate-independent), as in src.
Malec's numbers imply a gap-to-contact density ratio at V_g = 0 of 85 V / 120 V (same gate leverage
in both regions); the model's ratio is 1 - (1-f(L/2))^2 for the pinned-edge profile.
Validation (exact limit): n_c = 0 -> V_NP(gap) = 0 for every lambda.
"""
import os, sys, csv
import numpy as np
from scipy.optimize import minimize_scalar, brentq
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import graphene_contact_doping_model as cdm
from graphene_fet_model import v_F
from scipy.constants import e, hbar, epsilon_0

sig = cdm._sheet_conductivity_puddle_regularized
D = {r["quantity"]: r for r in csv.DictReader(l for l in open(os.path.join(HERE, "malec2011_clean_au.csv")) if not l.startswith("#"))}
C_G = float(D["C_g"]["value"]) * 1e-9 / 1e-4          # F/m^2
NPV = C_G / e                                          # m^-2 per V
L_GAP = float(D["gap_channel_length_L"]["value"]) * 1e-9
V_BULK = float(D["Vg_max_bulk_high_coverage"]["value"])
V_RC = float(D["Vg_max_contact_resistance"]["value"])
N_MALEC = -float(D["n_added_at_120V"]["value"]) * 1e4  # signed, p-type, m^-2
DEF_MALEC = float(D["dEF_contact"]["value"])

def S(x, L, lam):
    a, b = 1 / (1 + x / lam), 1 / (1 + (L - x) / lam)
    return a + b - a * b

def r_gap(n_c, v, lam, L=L_GAP, npts=2000):
    x = np.linspace(0, L, npts)
    n = v * NPV + (n_c - v * NPV) * S(x, L, lam)
    return np.trapezoid(1 / sig(n), x)

def v_np(n_c, lam, L=L_GAP, vmax=5e4):
    # coarse log-spaced scan then bounded refine; the peak can sit far outside +-300 V for long tails
    vs = np.concatenate([-np.geomspace(vmax, 1e-3, 400), [0.0], np.geomspace(1e-3, vmax, 400)])
    R = [r_gap(n_c, v, lam, L) for v in vs]
    i = int(np.argmax(R)); lo, hi = vs[max(i - 1, 0)], vs[min(i + 1, len(vs) - 1)]
    return minimize_scalar(lambda v: -r_gap(n_c, v, lam, L), bounds=(lo, hi), method="bounded",
                           options=dict(xatol=1e-6 * max(1, abs(vs[i])))).x

def dEF(n): return np.sign(-n) * hbar * v_F * np.sqrt(np.pi * abs(n)) / e   # +: p-type (E_F below Dirac)

T = []; P = T.append
v0 = [v_np(0.0, lam) for lam in (5e-9, 50e-9, 250e-9)]
P(f"VALIDATION n_c = 0: V_NP(gap) at lambda = 5/50/250 nm = {', '.join(f'{v:+.1e}' for v in v0)} V (exact 0) -> "
  f"{'PASS' if max(map(abs, v0)) < 1e-2 else 'FAIL'}")
P(f"C_g = {C_G*1e5:.2f} nF/cm^2 ({NPV*1e-4:.3e} cm^-2/V), gap L = {L_GAP*1e9:.0f} nm, src lambda = {cdm.lambda_decay*1e9:.0f} nm, src v_F = {v_F:.1e} m/s")
P("")
K = {r["metal"]: r for r in csv.DictReader(l for l in open(os.path.join(HERE, "..", "I2-au-facet-work-function", "khomyakov2009_table1.csv")) if not l.startswith("#"))}
W = float(K["Au"]["W_M_111_calc_eV"]); n_model = cdm.contact_edge_carrier_density(W)
P("1) SIGN (clean Au)")
P(f"   model: W(111, calc) = {W} eV -> n_c = {n_model*1e-4:+.2e} cm^-2 -> {'p' if n_model < 0 else 'n'}-type ; Malec 2011: {D['doping_type']['value']}-type"
  f" -> {'MATCH' if (n_model < 0) == (D['doping_type']['value'] == 'p') else 'MISMATCH'}")
P("")
P("2) SIZE of the contact doping")
P(f"   model dE_F = {dEF(n_model):+.3f} eV ; Malec {DEF_MALEC:+.2f} eV (|n| = {abs(N_MALEC)*1e-4:.1e} cm^-2) ; "
  f"Khomyakov 2009 DFT {float(K['Au']['dEF_printed_eV']):+.2f} eV ; ratio model/Malec n = {n_model/N_MALEC:.2f}")
W_fit = brentq(lambda w: cdm.contact_edge_carrier_density(w) - N_MALEC, 5.4, 9.0)
P(f"   src-equivalent W that gives Malec's density: {W_fit:.2f} eV (vs 5.54 calc(111), 5.31 exp)")
P("")
P("3) TAIL LENGTH: can a pinned edge with lambda = 250 nm leave a 35 nm gap gate-tunable to neutrality at +85 V?")
ratio_data = V_BULK / V_RC
P(f"   data gap/contact density ratio at V_g = 0: {V_BULK:.0f}/{V_RC:.0f} = {ratio_data:.2f}")
rows = []
P("   lambda (nm) | S_or(L/2) | V_NP(gap), n_c = Malec (V) | V_NP(gap), n_c = model (V)")
for lam in (2e-9, 5e-9, 10e-9, 15e-9, 25e-9, 50e-9, 100e-9, 250e-9):
    s_mid = S(L_GAP / 2, L_GAP, lam); va = v_np(N_MALEC, lam); vb = v_np(n_model, lam)
    rows.append(dict(lambda_nm=lam * 1e9, S_or_mid=s_mid, VNP_gap_nc_malec_V=va, VNP_gap_nc_model_V=vb))
    P(f"   {lam*1e9:6.0f} | {s_mid:.3f} | {va:+10.1f} | {vb:+10.1f}")
lam_ratio = brentq(lambda lam: S(L_GAP / 2, L_GAP, lam) - ratio_data, 1e-10, 1e-6)
lam_85 = brentq(lambda lam: v_np(N_MALEC, lam) - V_BULK, 1e-10, 2e-7, xtol=1e-12)
P(f"   lambda matching the density ratio {ratio_data:.2f}: {lam_ratio*1e9:.1f} nm ; lambda putting V_NP(gap) at +{V_BULK:.0f} V with Malec n_c: {lam_85*1e9:.1f} nm")
P(f"   src lambda = 250 nm puts the gap neutrality at {rows[-1]['VNP_gap_nc_malec_V']:+.0f} V (Malec n_c) -- "
  f"{rows[-1]['VNP_gap_nc_malec_V']/V_BULK:.0f}x the measured +{V_BULK:.0f} V")
P("")
P("4) CROSS-CHECK against 11-20 (Nouchi 2011, Au, dV_NP(0.4 vs 2.5 um) = -28 V on 300 nm SiO2, C_ox = 11.51 nF/cm^2)")
C_N = 3.9 * epsilon_0 / 300e-9; NPV_N = C_N / e
def r_ch(n_c, v, lam, L, npts=3000):
    h = np.geomspace(1e-12, L / 2, npts); x = np.concatenate([[0.0], h])
    n = v * NPV_N + (n_c - v * NPV_N) * S(x, L, lam)
    return 2 * np.trapezoid(1 / sig(n), x)
def vnp_ch(n_c, lam, L):
    return minimize_scalar(lambda v: -r_ch(n_c, v, lam, L), bounds=(-300, 300), method="bounded", options=dict(xatol=1e-4)).x
def sh(n_c, lam): return vnp_ch(n_c, lam, 0.4e-6) - vnp_ch(n_c, lam, 2.5e-6)
for lam in (lam_85, 250e-9):
    grid = np.geomspace(1e14, 5e17, 40); g = [sh(nc, lam) + 28 for nc in grid]   # first sign change (the +-300 V window clips large n_c)
    k = next((i for i in range(len(grid) - 1) if g[i] * g[i + 1] < 0), None)
    if k is None:
        P(f"   lambda = {lam*1e9:5.1f} nm: no n_c up to 5e13 cm^-2 reproduces Nouchi -28 V")
    else:
        nc = brentq(lambda nc: sh(nc, lam) + 28, grid[k], grid[k + 1], xtol=1e12)
        P(f"   lambda = {lam*1e9:5.1f} nm: n_c reproducing Nouchi -28 V = {nc*1e-4:+.2e} cm^-2 (|n_c|/|n_Malec| = {abs(nc/N_MALEC):.2f})")
with open(os.path.join(HERE, "malec_clean_au_results.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
open(os.path.join(HERE, "malec_clean_au_transcript.txt"), "w").write("\n".join(T) + "\n")
print("\n".join(T))

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
lams = np.geomspace(1e-9, 300e-9, 40)
ax[0].loglog(lams * 1e9, [v_np(N_MALEC, l) for l in lams], "C0-", label="model gap V_NP, n_c = Malec (−9.3e12)")
ax[0].loglog(lams * 1e9, [v_np(n_model, l) for l in lams], "C1-", label="model gap V_NP, n_c = src Au(111) (−1.25e12)")
ax[0].axhline(V_BULK, c="k", ls="--", lw=1, label="Malec 2011: +85 V (text)")
ax[0].axvline(250, c="gray", ls=":", lw=1); ax[0].annotate("src λ", (255, 2), fontsize=8)
ax[0].set_xlabel("tail length λ (nm)"); ax[0].set_ylabel("gap V_NP (V)"); ax[0].legend(fontsize=7)
ax[0].set_title("35 nm gap between clean Au grains", fontsize=9)
x = np.linspace(0, L_GAP, 200)
for lam, c in ((lam_85, "C0"), (250e-9, "C3")):
    ax[1].plot(x * 1e9, S(x, L_GAP, lam), c, label=f"S_or, λ = {lam*1e9:.0f} nm")
ax[1].axhline(ratio_data, c="k", ls="--", lw=1, label=f"Malec gap/contact ratio {ratio_data:.2f}")
ax[1].set_xlabel("x in gap (nm)"); ax[1].set_ylabel("contact share of doping S"); ax[1].set_ylim(0, 1.05); ax[1].legend(fontsize=7)
ax[1].set_title("doping profile across the gap", fontsize=9)
fig.suptitle("I1: src contact-doping model vs clean Au (Malec & Davidović 2011, arXiv:1108.3238)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "malec_clean_au.png"), dpi=150)
