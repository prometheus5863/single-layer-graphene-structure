"""I1, experiment 5 (2026-11-23): gate-tuned contact edge (C_q-limited) instead of src's pinned edge.

11-21 found that src's pinned, gate-independent edge density n_c is the weak point: Malec & Davidovic
(PRB 84, 033407 (2011), arXiv:1108.3238) see the contact-resistance maximum at V_g = +120 V, so the
graphene UNDER the Au is gate-tuned. src already models the contact as graphene's quantum capacitance
C_q in series with C_I = eps0/d_interface (src/graphene_contact_doping_model.py, unchanged). Here the
back-gate charge n_g = C_g V_g / e is added to the same electrostatics. Charge neutrality of
gate + graphene + metal, with the metal equipotential at the source, gives (signed n, + = electrons):

    E_F(n)/e + e (n - n_g) / C_I = dV,      E_F(n) = sgn(n) hbar v_F sqrt(pi |n|),   dV = -(W_M - W0)/e

At n_g = 0 this is the exact (non-linearised) form of src's contact_edge_carrier_density. The gate
leverage under the metal is dn/dn_g = C_q/(C_q + C_I): zero at neutrality, rising with |n|.
Consequence used below: the contact is neutral (R_c maximal) exactly at n_g* = -C_I dV / e.

Tail: n(x) = n_g + (n_c(V_g) - n_g) S_or(x), with the 11-20 bounded two-contact S_or and f = 1/(1+x/lam);
lam is a free parameter (gate screening of the tail, Khomyakov et al. PRB 82, 115437 (2010),
arXiv:0911.2027, is NOT implemented here -- only the edge is changed today).

Tests: (A) Malec clean Au: dV from the +120 V contact maximum; lam from the +85 V gap maximum.
       (B) Nouchi 2011 Au (Cr/Au leads): -28 V shift between L = 0.4 and 2.5 um, with the same lam.
Validation (exact limit): C_I -> 0 (metal decoupled) must give n_c(V_g) = n_g exactly.
"""
import os, sys, csv
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.constants import e, hbar, epsilon_0
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import graphene_contact_doping_model as cdm
from graphene_fet_model import v_F

sig = cdm._sheet_conductivity_puddle_regularized
C_I = cdm.C_interface
D = {r["quantity"]: r for r in csv.DictReader(l for l in open(os.path.join(HERE, "malec2011_clean_au.csv")) if not l.startswith("#"))}
C_M = float(D["C_g"]["value"]) * 1e-5            # F/m^2 (nF/cm^2 -> F/m^2)
L_GAP = float(D["gap_channel_length_L"]["value"]) * 1e-9
V_BULK = float(D["Vg_max_bulk_high_coverage"]["value"]); V_RC = float(D["Vg_max_contact_resistance"]["value"])
C_N = 3.9 * epsilon_0 / 300e-9                   # Nouchi 2011, 300 nm SiO2
NOUCHI_SHIFT = -71.0 - (-43.0)                   # V_NP(0.4 um) - V_NP(2.5 um), text values, -28 V

def EF(n): return np.sign(n) * hbar * v_F * np.sqrt(np.pi * np.abs(n)) / e   # eV, + = electrons

def n_edge(dV, n_g, CI=C_I):
    g = lambda n: EF(n) + e * (n - n_g) / CI - dV
    span = 1e18 + abs(n_g) + abs(CI * dV / e) * 2
    return brentq(g, -span, span, xtol=1e6, rtol=1e-13)

def S(x, L, lam):
    a, b = 1 / (1 + x / lam), 1 / (1 + (L - x) / lam)
    return a + b - a * b

def R_strip(dV, Vg, lam, L, Cg, npts=3000):
    n_g = Cg * Vg / e; nc = n_edge(dV, n_g)
    h = np.geomspace(1e-12, L / 2, npts); x = np.concatenate([[0.0], h])
    n = n_g + (nc - n_g) * S(x, L, lam)
    return 2 * np.trapezoid(1 / sig(n), x)

def vnp(dV, lam, L, Cg, lo=-400, hi=400):
    vs = np.linspace(lo, hi, 161); R = [R_strip(dV, v, lam, L, Cg) for v in vs]
    i = int(np.argmax(R))
    return minimize_scalar(lambda v: -R_strip(dV, v, lam, L, Cg), bounds=(vs[max(i-1, 0)], vs[min(i+1, 160)]),
                           method="bounded", options=dict(xatol=1e-3)).x

T = []; P = T.append
# --- validation: decoupled metal, C_I -> 0 => n_c = n_g exactly
errs = [abs(n_edge(-0.5, ng, CI=1e-12) - ng) * 1e-4 for ng in (-5e16, -1e15, 0.0, 2e15, 7e16)]   # cm^-2
P(f"VALIDATION C_I -> 1e-12 F/m^2 (metal decoupled): max |n_c - n_g| over n_g = -5e12..+7e12 cm^-2 = {max(errs):.1e} cm^-2 (exact 0) -> {'PASS' if max(errs) < 1e3 else 'FAIL'}")
nc_src = cdm.contact_edge_carrier_density(5.54)
P(f"consistency (not a validation): n_g = 0, dV = -(5.54-5.40): exact series eq. n_c = {n_edge(-0.14, 0.0)*1e-4:+.3e} cm^-2 vs src linearised {nc_src*1e-4:+.3e}")
P(f"C_I = eps0/{cdm.d_interface*1e10:.1f} A = {C_I:.3e} F/m^2 ; Malec C_g = {C_M*1e5:.1f} nF/cm^2 ; C_I/C_g = {C_I/C_M:.0f}")
P("")
# --- (A1) contact maximum -> dV
P("(A1) Malec clean Au: R_c maximum at V_g = +120 V  ->  n_c = 0 there  ->  dV = -e n_g*/C_I")
rows = []
for d in (3.3e-10, 3.5e-10, 4.0e-10):
    CI = epsilon_0 / d; dV = -C_M * V_RC / CI   # from E_F(0) = 0: -e n_g*/C_I = dV
    nc0 = n_edge(dV, 0.0, CI=CI)
    rows.append(dict(d_A=d * 1e10, dV_V=dV, W_minus_W0_eV=-dV, nc0_cm2=nc0 * 1e-4, dEF0_eV=-EF(nc0)))
    P(f"   d_I = {d*1e10:.1f} A: W_M - W0 = {-dV:+.3f} eV -> n_c(V_g=0) = {nc0*1e-4:+.2e} cm^-2, E_D - E_F = {-EF(nc0):+.3f} eV")
dV_M = -C_M * V_RC / C_I
P(f"   src model Au(111): W - W0 = +0.14 eV -> predicted R_c maximum at V_g = {C_I*0.14/C_M:+.0f} V (measured +{V_RC:.0f} V)")
P(f"   Malec's own reading (gate fully reaches graphene under Au): 9.3e12 cm^-2, 0.35 eV; series-capacitor reading at d_I = 3.5 A: "
  f"{rows[1]['nc0_cm2']:+.2e} cm^-2, {rows[1]['dEF0_eV']:.2f} eV")
P("")
# --- (A2) gap maximum -> lambda, with the gate-tuned edge
P("(A2) 35 nm gap between grains, gate-tuned edge (dV from A1, d_I = 3.5 A): gap V_NP vs lambda (measured +85 V)")
lams = [2e-9, 5e-9, 10e-9, 15e-9, 25e-9, 50e-9, 100e-9, 250e-9]
gap = []
for lam in lams:
    v = vnp(dV_M, lam, L_GAP, C_M); gap.append(v); P(f"   lambda = {lam*1e9:5.0f} nm: gap V_NP = {v:+7.1f} V")
lam_fit = brentq(lambda l: vnp(dV_M, l, L_GAP, C_M) - V_BULK, 1e-10, 250e-9, xtol=1e-11)
P(f"   lambda fitting +85 V: {lam_fit*1e9:.1f} nm  (11-21, pinned edge: 5 nm; src: 250 nm)")
P("")
# --- (B) Nouchi with the same lambda
P("(B) Nouchi 2011 Au (Cr/Au leads, 300 nm SiO2): dV_NP(0.4 um) - dV_NP(2.5 um), measured -28 V (n-type)")
def shift(dV, lam): return vnp(dV, lam, 0.4e-6, C_N, -300, 300) - vnp(dV, lam, 2.5e-6, C_N, -300, 300)
for lam in (lam_fit, 250e-9):
    s_same = shift(-dV_M, lam)   # Malec |dV| with the sign flipped to n-type
    P(f"   lambda = {lam*1e9:5.1f} nm, |W - W0| = {abs(dV_M):.3f} eV (n-type): shift = {s_same:+.2f} V")
dvs = np.geomspace(0.05, 20, 30); sh = [shift(d, lam_fit) for d in dvs]
k = next((i for i in range(len(dvs) - 1) if (sh[i] - NOUCHI_SHIFT) * (sh[i+1] - NOUCHI_SHIFT) < 0), None)
if k is None:
    P(f"   lambda = {lam_fit*1e9:.1f} nm: no |W - W0| in 0.05-20 eV reaches -28 V (largest |shift| = {max(map(abs, sh)):.2f} V)")
    dv_need = np.nan
else:
    dv_need = brentq(lambda d: shift(d, lam_fit) - NOUCHI_SHIFT, dvs[k], dvs[k+1], xtol=1e-3)
    P(f"   lambda = {lam_fit*1e9:.1f} nm: -28 V needs |W - W0| = {dv_need:.2f} eV (Malec: {abs(dV_M):.2f} eV)")
P(f"   gate-tuned edge, saturation: as |W - W0| -> large the shift tends to {sh[-1]:+.2f} V at lambda = {lam_fit*1e9:.1f} nm")

with open(os.path.join(HERE, "gate_tuned_edge_results.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["test", "parameter", "value", "unit", "note"])
    for r in rows: w.writerow(["A1", f"d_I={r['d_A']:.1f}A", f"{r['W_minus_W0_eV']:.4f}", "eV", f"n_c0={r['nc0_cm2']:.3e} cm^-2; E_D-E_F={r['dEF0_eV']:.3f} eV"])
    for l, v in zip(lams, gap): w.writerow(["A2", f"lambda={l*1e9:.0f}nm", f"{v:.2f}", "V", "gap V_NP, gate-tuned edge"])
    w.writerow(["A2", "lambda_fit", f"{lam_fit*1e9:.2f}", "nm", "gap V_NP = +85 V"])
    w.writerow(["B", "|W-W0| needed for Nouchi -28 V", f"{dv_need:.3f}", "eV", f"lambda={lam_fit*1e9:.1f} nm"])
    for d, s in zip(dvs, sh): w.writerow(["B", f"|W-W0|={d:.3f}eV", f"{s:.3f}", "V", f"Nouchi shift, lambda={lam_fit*1e9:.1f} nm"])
open(os.path.join(HERE, "gate_tuned_edge_transcript.txt"), "w").write("\n".join(T) + "\n")
print("\n".join(T))

fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
vg = np.linspace(-50, 250, 301)
for dv, lab, c in ((dV_M, f"fit to Malec, W−W0 = {-dV_M:.2f} eV", "C0"), (-0.14, "src Au(111), W−W0 = 0.14 eV", "C1")):
    ax[0].plot(vg, [n_edge(dv, C_M * v / e) * 1e-16 for v in vg], c, label=lab)
ax[0].plot(vg, C_M * vg / e * 1e-16 - 9.3, "k:", lw=1, label="Malec reading (full gate leverage)")
ax[0].axvline(V_RC, c="k", ls="--", lw=1); ax[0].axhline(0, c="gray", lw=0.5)
ax[0].set_xlabel("V_g (V)"); ax[0].set_ylabel("n under Au (10¹² cm⁻²)"); ax[0].legend(fontsize=7); ax[0].set_title("gate-tuned edge (C_q in series with C_I)", fontsize=9)
ax[1].semilogx(np.array(lams) * 1e9, gap, "C0o-", label="gap V_NP, gate-tuned edge")
ax[1].axhline(V_BULK, c="k", ls="--", lw=1, label="Malec 2011: +85 V (text)")
ax[1].axvline(lam_fit * 1e9, c="C0", ls=":", lw=1)
ax[1].set_xlabel("λ (nm)"); ax[1].set_ylabel("gap V_NP (V)"); ax[1].legend(fontsize=7); ax[1].set_title("35 nm gap between clean Au grains", fontsize=9)
ax[2].semilogx(dvs, sh, "C2o-", ms=3, label=f"model shift, λ = {lam_fit*1e9:.0f} nm")
ax[2].axhline(NOUCHI_SHIFT, c="k", ls="--", lw=1, label="Nouchi 2011 Au: −28 V (text)")
ax[2].axvline(abs(dV_M), c="C0", ls=":", lw=1, label="Malec |W−W0|")
ax[2].set_xlabel("|W − W0| (eV), n-type"); ax[2].set_ylabel("V_NP(0.4 µm) − V_NP(2.5 µm) (V)"); ax[2].legend(fontsize=7); ax[2].set_title("same λ applied to Nouchi", fontsize=9)
fig.suptitle("I1: gate-tuned contact edge vs Malec & Davidović 2011 (arXiv:1108.3238) and Nouchi 2011 (arXiv:1102.2954); text values only", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "gate_tuned_edge.png"), dpi=150)
