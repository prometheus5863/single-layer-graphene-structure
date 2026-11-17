"""I2, experiment 1 (2026-11-17): rerun the Section 4.7 contact-doping model with
facet-matched DFT (111) work functions (Giovannetti et al. PRL 101, 026803 (2008),
Table I) in place of the polycrystalline table, and compare the model's contact
Fermi-level shift with the DFT shift from the same table.

Imports src/ unchanged. Writes outputs next to this file.
"""
import os, sys, csv
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.constants import hbar, e
import graphene_contact_doping_model as cdm
from graphene_fet_model import v_F

N_BULK = 2.0e16  # 1/m^2, same on-state density as recalibrate_metal_rc()

def ef_shift_eV(n_signed):
    """Fermi-level shift for a signed density (n>0 electron-doped), in the
    paper's sign convention (n-type negative): -sign(n) hbar v_F sqrt(pi|n|)."""
    return -np.sign(n_signed) * hbar * v_F * np.sqrt(np.pi * abs(n_signed)) / e

def giovannetti_eq1(WM, d, WG=4.5):
    """Giovannetti 2008 Eq. (1), universal parameters as printed."""
    alpha, D0, d0 = 34.93, 0.09, 2.4
    kap, a0, a1, a2 = 1.6443, -2048.56, 1363.87, -205.737
    dc = np.exp(-kap * d) * (a0 + a1 * d + a2 * d * d)
    x = WM - WG - dc
    k = alpha * D0 * (d - d0)
    return np.sign(x) * (np.sqrt(1 + 2 * k * abs(x)) - 1) / k, dc

rows = list(csv.DictReader(l for l in open(os.path.join(HERE, "giovannetti2008_table1.csv")) if not l.startswith("#")))
out = []
lines = []
P = lines.append

# --- the one validation: at W = W0 the model must leave graphene undoped ---
n0 = cdm.contact_edge_carrier_density(cdm.W_CROSSOVER)
P(f"VALIDATION n_contact(W0={cdm.W_CROSSOVER}) = {n0:.3e} m^-2 (exact limit: 0) -> {'PASS' if n0 == 0 else 'FAIL'}")

P("")
P("metal | W_poly | W_111 | sign(poly) | sign(111) | dEF_model(111) | dEF_Eq1(d_eq) | dEF_DFT | R_extra_poly | R_extra_111 (Ohm.um)")
for r in rows:
    m = r["metal"]; W111 = float(r["W_M_111_eV"]); d = float(r["d_eq_A"]); dft = float(r["dEF_DFT_eV"])
    Wp = cdm.METAL_WORK_FUNCTIONS.get(m)
    n111 = cdm.contact_edge_carrier_density(W111)
    R111 = cdm.junction_extra_resistance(W111, N_BULK)[0]
    if Wp is not None:
        npoly = cdm.contact_edge_carrier_density(Wp); Rp = cdm.junction_extra_resistance(Wp, N_BULK)[0]
        sp = "n" if npoly > 0 else "p"
    else:
        npoly = Rp = np.nan; sp = "-"
    s111 = "n" if n111 > 0 else "p"
    efm = ef_shift_eV(n111); efq, dc = giovannetti_eq1(W111, d)
    out.append(dict(metal=m, W_poly=Wp, W_111=W111, sign_poly=sp, sign_111=s111,
                    n_contact_111_m2=n111, dEF_model_111_eV=efm, dEF_Eq1_eV=efq,
                    dEF_DFT_eV=dft, R_extra_poly_ohm_um=Rp, R_extra_111_ohm_um=R111))
    P(f"{m:2s} | {Wp if Wp else '-':>5} | {W111:.2f} | {sp} | {s111} | {efm:+.3f} | {efq:+.3f} | {dft:+.2f} | {Rp:9.1f} | {R111:9.1f}")

ef_m = np.array([o["dEF_model_111_eV"] for o in out]); ef_d = np.array([o["dEF_DFT_eV"] for o in out])
ef_q = np.array([o["dEF_Eq1_eV"] for o in out])
P("")
P(f"sign agreement model(111) vs DFT: {int(np.sum(np.sign(ef_m)==np.sign(ef_d)))}/{len(out)}")
P(f"RMS |dEF| error model vs DFT: {np.sqrt(np.mean((ef_m-ef_d)**2)):.3f} eV; Eq.1 vs DFT: {np.sqrt(np.mean((ef_q-ef_d)**2)):.3f} eV")
P(f"ratio model/DFT per metal: " + ", ".join(f"{o['metal']} {o['dEF_model_111_eV']/o['dEF_DFT_eV']:.2f}" for o in out))
flips = [o["metal"] for o in out if o["sign_poly"] not in ("-", o["sign_111"])]
P(f"physisorbed metals whose doping sign flips poly -> (111): {flips or 'none'}")

with open(os.path.join(HERE, "facet_matched_results.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
open(os.path.join(HERE, "facet_matched_transcript.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
Ws = np.linspace(4.0, 6.3, 200)
ax[0].plot(Ws, [ef_shift_eV(cdm.contact_edge_carrier_density(W)) for W in Ws], "C0-", label="this repo's model (src/, unchanged)")
ax[0].plot(Ws, [giovannetti_eq1(W, 3.3)[0] for W in Ws], "C1--", label="Giovannetti 2008 Eq. (1), d = 3.3 Å")
ax[0].plot([o["W_111"] for o in out], ef_d, "ks", label="DFT, Giovannetti 2008 Table I")
for o in out: ax[0].annotate(o["metal"], (o["W_111"], o["dEF_DFT_eV"]), xytext=(4, -10), textcoords="offset points")
ax[0].axhline(0, c="gray", lw=.6); ax[0].axvline(cdm.W_CROSSOVER, c="gray", ls=":", lw=.8)
ax[0].set_xlabel(r"metal work function $W_M$ (eV)"); ax[0].set_ylabel(r"$\Delta E_F$ (eV), n-type < 0")
ax[0].legend(fontsize=8); ax[0].set_title("Contact Fermi-level shift vs DFT")
ms = [o["metal"] for o in out if o["W_poly"]]
xi = np.arange(len(ms))
ax[1].bar(xi-.2, [o["R_extra_poly_ohm_um"] for o in out if o["W_poly"]], .4, label="polycrystalline W (src table)")
ax[1].bar(xi+.2, [o["R_extra_111_ohm_um"] for o in out if o["W_poly"]], .4, label="(111) W, Giovannetti 2008")
ax[1].set_xticks(xi, ms); ax[1].axhline(0, c="gray", lw=.6)
ax[1].set_ylabel(r"$R_{extra}$ ($\Omega\,\mu$m), $n_{bulk}=2\times10^{16}$ m$^{-2}$")
ax[1].legend(fontsize=8); ax[1].set_title(r"§4.7 $R_{extra}$, physisorbed metals")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "facet_matched_W.png"), dpi=150)
