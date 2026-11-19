"""I1, experiment 1 (2026-11-19): simulate a transfer-length-method (TLM) extraction on
top of the repo's contact-doping profile (src/graphene_contact_doping_model.py, unchanged),
so that the TLM-extracted contact resistance and the model's R_extra are defined on the
same device. Metal work functions: calculated (111) values from Khomyakov et al. PRB 79,
195425 (2009) Table I (validated in I2), with W0 = 5.4 eV from src.

Device: two identical contacts at x=0 and x=L, channel gated to n_bulk. The excess
carrier density from each contact is superposed linearly (linear electrostatic response):
    n(x) = n_bulk + (n_c - n_bulk) [f(x) + f(L-x)],  f = 1/(1+x/lambda)   (src doping_profile)
Total width-normalised resistance R_tot(L) W = 2 R_int + int_0^L dx / sigma(n(x)), with the
src puddle-regularised Drude sheet conductivity. R_int (interface transmission) is unknown
and only shifts the intercept, so it is set to 0: R_c,TLM below is the *junction* part of
what a TLM fit would report. A straight line is fitted to R_tot(L) over a TLM length window;
R_c,TLM = intercept/2.
Validation (exact limit): n_c = n_bulk -> R_tot = L/sigma_bulk exactly, intercept 0, slope R_sheet.
"""
import os, sys, csv
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import graphene_contact_doping_model as cdm

N_BULK = 2.0e16  # 1/m^2, same on-state density src uses for R_extra / recalibration
lam = cdm.lambda_decay
sig = cdm._sheet_conductivity_puddle_regularized
f = lambda x: 1.0 / (1.0 + x / lam)

def r_tot(L, n_c, n_b=N_BULK, npts=4000):
    """Width-normalised junction+channel resistance (Ohm.um) of a device of length L (m)."""
    h = np.geomspace(1e-12, L / 2, npts); x = np.concatenate([[0.0], h])  # dense near each contact
    n = n_b + (n_c - n_b) * (f(x) + f(L - x))
    half = np.trapezoid(1.0 / sig(n), x)  # profile symmetric about L/2
    return 2 * half * 1e6

def tlm_fit(n_c, Lmin, Lmax, k=8):
    Ls = np.linspace(Lmin, Lmax, k)
    R = np.array([r_tot(L, n_c) for L in Ls])
    slope, icpt = np.polyfit(Ls * 1e6, R, 1)
    resid = R - (slope * Ls * 1e6 + icpt)
    return icpt / 2, slope, np.max(np.abs(resid)), Ls, R

def rd(fn): return list(csv.DictReader(l for l in open(fn) if not l.startswith("#")))
K = rd(os.path.join(HERE, "..", "I2-au-facet-work-function", "khomyakov2009_table1.csv"))
L_ = []; P = L_.append
Rsh = 1e6 / sig(N_BULK) / 1e6  # Ohm/sq ; slope in Ohm.um per um = Ohm/sq

# --- the one validation (exact limit)
rc0, s0, _, _, _ = tlm_fit(N_BULK, 1e-6, 10e-6)
ok = abs(rc0) < 1e-6 * Rsh and abs(s0 - Rsh) / Rsh < 1e-9
P(f"VALIDATION n_c = n_bulk: R_c,TLM = {rc0:.2e} Ohm.um (exact 0), slope = {s0:.6f} vs R_sheet = {Rsh:.6f} Ohm/sq -> {'PASS' if ok else 'FAIL'}")
P(f"n_bulk = {N_BULK:.1e} m^-2, lambda = {lam*1e9:.0f} nm, R_sheet = {Rsh:.1f} Ohm/sq")
P("")
windows = [(1e-6, 10e-6), (2e-6, 20e-6), (0.5e-6, 5e-6)]
P("metal | W111 | n_c (m^-2) | R_extra src(L_j=1um) | R_extra(L_j=10um) | R_extra(L_j=100um) | R_c,TLM [1-10um] | [2-20] | [0.5-5] | ratio TLM[1-10]/R_extra(1um) | max fit resid [1-10]")
out = []
for r in K:
    m = r["metal"]; W = float(r["W_M_111_calc_eV"])
    Rex1, _, _, n_c = cdm.junction_extra_resistance(W, N_BULK)  # src, L_junction = 1 um
    Rex10 = cdm.junction_extra_resistance(W, N_BULK, L_junction=10e-6, n_points=20000)[0]
    Rex100 = cdm.junction_extra_resistance(W, N_BULK, L_junction=100e-6, n_points=200000)[0]
    fits = [tlm_fit(n_c, a, b) for a, b in windows]
    ratio = fits[0][0] / Rex1
    o = dict(metal=m, W111_eV=W, n_c_m2=n_c, R_extra_src_1um=Rex1, R_extra_10um=Rex10, R_extra_100um=Rex100,
             Rc_TLM_1_10um=fits[0][0], Rc_TLM_2_20um=fits[1][0], Rc_TLM_05_5um=fits[2][0],
             ratio_TLM_over_Rextra=ratio, slope_1_10=fits[0][1], max_resid_1_10=fits[0][2])
    out.append(o)
    P(f"{m:2s} | {W:.2f} | {n_c:+.2e} | {Rex1:+8.1f} | {Rex10:+8.1f} | {Rex100:+8.1f} | {fits[0][0]:+8.1f} | {fits[1][0]:+8.1f} | {fits[2][0]:+8.1f} | {ratio:5.2f} | {fits[0][2]:.2f}")
P("")
absorbed = [o for o in out if abs(o["ratio_TLM_over_Rextra"]) >= 0.10]
P(f"kill condition (R_c,TLM absorbs < 10% of R_extra for every physisorbed metal): {'HIT' if not absorbed else 'NOT hit'}; metals with >= 10%: {[o['metal'] for o in absorbed]}")
P(f"apparent sheet resistance from [1-10um] slope vs true R_sheet {Rsh:.1f}: " + ", ".join(f"{o['metal']} {o['slope_1_10']:.1f}" for o in out))
with open(os.path.join(HERE, "tlm_results.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
open(os.path.join(HERE, "tlm_transcript.txt"), "w").write("\n".join(L_) + "\n")
print("\n".join(L_))

fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
Lg = np.linspace(0.3e-6, 20e-6, 60)
for i, o in enumerate(out):
    R = np.array([r_tot(L, o["n_c_m2"]) for L in Lg]) - Lg * 1e6 * Rsh
    ax[0].plot(Lg * 1e6, R, f"C{i}-", label=o["metal"])
ax[0].axhline(0, c="gray", lw=.6); ax[0].axvspan(1, 10, color="gray", alpha=.1, label="TLM window 1–10 µm")
ax[0].set_xlabel("channel length L (µm)"); ax[0].set_ylabel(r"$R_{tot}W - L R_{sheet}$ (Ω·µm)")
ax[0].set_title("junction excess in a two-contact device"); ax[0].legend(fontsize=8)
xs = np.arange(len(out)); wdt = .2
for j, (k, lab) in enumerate([("R_extra_src_1um", "R_extra (src, 1 µm)"), ("R_extra_10um", "R_extra (10 µm)"),
                              ("Rc_TLM_1_10um", "R_c,TLM [1–10 µm]"), ("Rc_TLM_2_20um", "R_c,TLM [2–20 µm]")]):
    ax[1].bar(xs + (j - 1.5) * wdt, [o[k] for o in out], wdt, label=lab)
ax[1].set_xticks(xs); ax[1].set_xticklabels([o["metal"] for o in out]); ax[1].axhline(0, c="k", lw=.6)
ax[1].set_ylabel("Ω·µm (per contact)"); ax[1].set_title("what TLM sees vs the model's R_extra"); ax[1].legend(fontsize=7.5)
fig.suptitle("I1: TLM extraction on the §4.5 doping profile, calc (111) W (Khomyakov 2009)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "tlm_extraction.png"), dpi=150)
