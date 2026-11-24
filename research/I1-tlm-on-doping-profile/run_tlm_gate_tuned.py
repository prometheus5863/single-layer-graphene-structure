"""I1, experiment 6 (2026-11-24): rerun the 11-19 TLM extraction on the 11-23 gate-tuned contact edge.

Profile (from run_gate_tuned_edge.py, imported unchanged): the graphene under Au obeys
    E_F(n_c)/e + e (n_c - n_g)/C_I = dV,  dV = -(W_M - W0)/e = -0.59 eV (Malec 2011 fit, d_I = 3.5 A)
and the channel follows n(x) = n_g + (n_c - n_g) S_or(x, L, lam), lam = 27 nm (Malec gap fit).
src (graphene_contact_doping_model.py) is unchanged; its puddle-regularised Drude sigma is used.

Assumption: the gate-induced density n_g equals the channel density n_bulk, i.e. densities are counted
from the channel's own neutrality point (Passi 2018 V_D = 38 V; the -40 V point is p-type). Passi gives
no oxide thickness, so the p-branch is bracketed with n_bulk = -2e16 and -5e16 m^-2 (as on 11-19b),
and the n-branch uses +2e16 and +5e16 m^-2.

TLM: R_tot(L) W = int_0^L dx / sigma(n(x)) (R_int = 0), fitted linearly over 8 lengths in three windows
(1-10, 2-20, 0.5-5 um); R_c,TLM = intercept / 2. R_extra (src definition, on this profile) is the
single-contact excess int_0^Lj [1/sigma(n_g + (n_c - n_g) f) - 1/sigma(n_g)] dx at Lj = 1 and 10 um.
Kill condition (IDEAS.md): R_c,TLM absorbs < 10 % of R_extra for every physisorbed metal.

Validation (exact limit): C_I -> 1e-12 F/m^2 (metal decoupled) gives n_c = n_g, a flat profile, so
R_c,TLM = 0 and slope = R_sheet exactly.
"""
import os, sys, csv
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.constants import e, hbar
import graphene_contact_doping_model as cdm
from graphene_fet_model import v_F
sig = cdm._sheet_conductivity_puddle_regularized

# copied verbatim from run_gate_tuned_edge.py (11-23); not imported, because importing re-runs it and rewrites its outputs
def EF(n): return np.sign(n) * hbar * v_F * np.sqrt(np.pi * np.abs(n)) / e   # eV, + = electrons
def n_edge(dV, n_g, CI=cdm.C_interface):
    gfun = lambda n: EF(n) + e * (n - n_g) / CI - dV
    span = 1e18 + abs(n_g) + abs(CI * dV / e) * 2
    return brentq(gfun, -span, span, xtol=1e6, rtol=1e-13)
def S(x, L, lam):
    a, b = 1 / (1 + x / lam), 1 / (1 + (L - x) / lam)
    return a + b - a * b
class g: C_I = cdm.C_interface; n_edge = staticmethod(n_edge)

DV_AU = -0.59          # eV, 11-23 Malec fit at d_I = 3.5 A (gate_tuned_edge_results.csv: 0.5897)
LAM = 27e-9            # m, 11-23 gap fit (gate_tuned_edge_results.csv)
PASSI_RC_P = 519.0     # Ohm.um at -40 V, passi2018_au_tlm.csv
WINDOWS = [(1e-6, 10e-6), (2e-6, 20e-6), (0.5e-6, 5e-6)]

def nc_of(n_b, dV=DV_AU, CI=g.C_I): return g.n_edge(dV, n_b, CI=CI)

def r_tot(L, n_b, n_c, lam=LAM, npts=4000):
    h = np.geomspace(1e-13, L / 2, npts); x = np.concatenate([[0.0], h])
    n = n_b + (n_c - n_b) * S(x, L, lam)
    return 2 * np.trapezoid(1 / sig(n), x) * 1e6           # Ohm.um

def r_extra(n_b, n_c, lam=LAM, Lj=1e-6, npts=20000):
    h = np.geomspace(1e-13, Lj, npts); x = np.concatenate([[0.0], h])
    n = n_b + (n_c - n_b) / (1 + x / lam)
    return np.trapezoid(1 / sig(n) - 1 / sig(n_b), x) * 1e6

def tlm(n_b, n_c, a, b, lam=LAM, k=8):
    Ls = np.linspace(a, b, k); R = np.array([r_tot(L, n_b, n_c, lam) for L in Ls])
    p = np.polyfit(Ls * 1e6, R, 1); fit = np.polyval(p, Ls * 1e6)
    r2 = 1 - np.sum((R - fit) ** 2) / np.sum((R - R.mean()) ** 2)
    return p[1] / 2, p[0], r2, np.max(np.abs(R - fit))

T = []; P = T.append
# --- validation
nb = 2e16; nc0 = nc_of(nb, CI=1e-12)
rc, sl, r2, _ = tlm(nb, nc0, 1e-6, 10e-6)
Rsh = 1 / sig(nb)
ok = abs(rc) < 1e-6 * Rsh and abs(sl - Rsh) / Rsh < 1e-6
P(f"VALIDATION C_I -> 1e-12 F/m^2: R_c,TLM = {rc:.2e} Ohm.um (exact 0), slope = {sl:.6f} vs R_sheet = {Rsh:.6f} Ohm/sq -> {'PASS' if ok else 'FAIL'}")
P(f"profile: Au W - W0 = {-DV_AU:.2f} eV, lambda = {LAM*1e9:.0f} nm, C_I = {g.C_I:.3e} F/m^2 (src d_I = 3.5 A); R_int = 0")
P("")
cases = [("p (-40 V, low)", -2e16), ("p (-40 V, high)", -5e16), ("n (low)", 2e16), ("n (high)", 5e16)]
P("channel | n_b (cm^-2) | n_c (cm^-2) | R_extra 1um | R_extra 10um | R_c,TLM [1-10] [2-20] [0.5-5] | R^2 [1-10] [2-20] [0.5-5] | TLM/R_extra(1um) | max resid [1-10] | apparent R_sh vs true")
rows = []
for lam in (LAM, 250e-9):
    P(f"-- lambda = {lam*1e9:.0f} nm" + ("  (11-23 fit)" if lam == LAM else "  (src value, for contrast)"))
    for lab, n_b in cases:
        n_c = nc_of(n_b)
        Rx1, Rx10 = r_extra(n_b, n_c, lam), r_extra(n_b, n_c, lam, Lj=10e-6, npts=60000)
        F = [tlm(n_b, n_c, a, b, lam) for a, b in WINDOWS]
        rat = F[0][0] / Rx1
        rows.append(dict(lambda_nm=lam * 1e9, channel=lab, n_b_cm2=n_b * 1e-4, n_c_cm2=n_c * 1e-4, R_extra_1um=Rx1, R_extra_10um=Rx10,
                         Rc_TLM_1_10=F[0][0], Rc_TLM_2_20=F[1][0], Rc_TLM_05_5=F[2][0], R2_1_10=F[0][2], R2_2_20=F[1][2], R2_05_5=F[2][2],
                         ratio_TLM_Rextra=rat, max_resid_1_10=F[0][3], slope_1_10=F[0][1], R_sheet_true=1 / sig(n_b)))
        P(f"{lab:16s} | {n_b*1e-4:+.1e} | {n_c*1e-4:+.2e} | {Rx1:+8.1f} | {Rx10:+8.1f} | {F[0][0]:+8.1f} {F[1][0]:+8.1f} {F[2][0]:+8.1f} | "
          f"{F[0][2]:.5f} {F[1][2]:.5f} {F[2][2]:.5f} | {rat:5.2f} | {F[0][3]:7.1f} | {F[0][1]:.0f} vs {1/sig(n_b):.0f}")
P("")
main = [r for r in rows if r["lambda_nm"] == LAM * 1e9]
absorbed = [r for r in main if abs(r["ratio_TLM_Rextra"]) >= 0.10]
P(f"kill condition for Au on the gate-tuned profile (lambda = 27 nm): {'NOT hit' if absorbed else 'HIT'} "
  f"(channels with |TLM/R_extra| >= 10 %: {[r['channel'] for r in absorbed]}); other metals not refitted on this edge")
for r in main[:2]:
    P(f"Passi 2018 -40 V ({r['channel']}): junction R_c,TLM = {r['Rc_TLM_05_5']:+.0f} to {r['Rc_TLM_1_10']:+.0f} Ohm.um over windows -> "
      f"{min(r['Rc_TLM_1_10'], r['Rc_TLM_2_20'], r['Rc_TLM_05_5'])/PASSI_RC_P*100:.0f}-{max(r['Rc_TLM_1_10'], r['Rc_TLM_2_20'], r['Rc_TLM_05_5'])/PASSI_RC_P*100:.0f} % of measured 519; "
      f"min R^2 = {min(r['R2_1_10'], r['R2_2_20'], r['R2_05_5']):.4f} (Passi: > 0.984)")
with open(os.path.join(HERE, "tlm_gate_tuned_results.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader()
    for r in rows: w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
open(os.path.join(HERE, "tlm_gate_tuned_transcript.txt"), "w").write("\n".join(T) + "\n")
print("\n".join(T))

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
Lg = np.linspace(0.2e-6, 20e-6, 70)
for i, (lab, n_b) in enumerate(cases):
    n_c = nc_of(n_b)
    R = np.array([r_tot(L, n_b, n_c) for L in Lg]) - Lg * 1e6 / sig(n_b)
    ax[0].plot(Lg * 1e6, R, f"C{i}-", label=f"{lab}, n = {n_b*1e-4:+.0e} cm⁻²")
ax[0].axhline(PASSI_RC_P * 2, c="k", ls="--", lw=1, label="2 × Passi R_c at −40 V (519 Ω·µm)")
ax[0].axvspan(1, 10, color="gray", alpha=.1); ax[0].axhline(0, c="gray", lw=.5)
ax[0].set_xlabel("L (µm)"); ax[0].set_ylabel("R_tot·W − L·R_sheet (Ω·µm), R_int = 0"); ax[0].legend(fontsize=7)
ax[0].set_title("two-contact junction excess, gate-tuned edge, λ = 27 nm", fontsize=9)
xs = np.arange(len(main)); wd = .2
for j, (k, lab) in enumerate([("R_extra_1um", "R_extra (1 µm)"), ("R_extra_10um", "R_extra (10 µm)"),
                              ("Rc_TLM_1_10", "R_c,TLM [1–10]"), ("Rc_TLM_2_20", "R_c,TLM [2–20]")]):
    ax[1].bar(xs + (j - 1.5) * wd, [r[k] for r in main], wd, label=lab)
ax[1].axhline(PASSI_RC_P, c="k", ls="--", lw=1, label="Passi −40 V (total R_c)")
ax[1].set_xticks(xs); ax[1].set_xticklabels([r["channel"] for r in main], fontsize=7); ax[1].axhline(0, c="k", lw=.5)
ax[1].set_ylabel("Ω·µm per contact"); ax[1].legend(fontsize=7); ax[1].set_title("Au: TLM sees vs R_extra", fontsize=9)
fig.suptitle("I1: TLM extraction on the gate-tuned contact edge (W−W0 = 0.59 eV, λ = 27 nm); Passi 2018 arXiv:1807.04772", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "tlm_gate_tuned.png"), dpi=150)
