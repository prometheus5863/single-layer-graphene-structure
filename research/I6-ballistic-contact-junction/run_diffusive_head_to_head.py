"""I6, experiment 3 (2026-11-28): diffusive vs ballistic junction term in the one-parameter full R_c, head to head.

Same procedure as run_full_rc.py (11-27): R_c = R_mg(n_c) + R_junction, the single parameter rho_c fixed so that
R_c = 519 Ohm.um at -40 V (Passi 2018, Au), then two no-refit predictions:
  (a) Passi Dirac-point ratio R_c(CNP)/R_c(-40 V), measured 1372/519 = 2.64
  (b) Malec 2011 ratio R_c(0 V)/R_c(-100 V) on C_g = 12.4 nF/cm^2, measured 128/95 = 1.35
R_junction is either
  ballistic : R_J = R_Landauer - R_Sharvin(channel)        (11-26 k_y-resolved Dirac transfer matrix)
  diffusive : R_extra = int [1/sigma(n(x)) - 1/sigma(n_b)] dx  (src's puddle-regularised sigma, 11-26 definition)
Difference from 11-27: rho_c is refitted for EACH model and EACH lambda (5, 15, 27 nm), so both models get the same
one-parameter freedom. 11-27 used the lambda = 27 nm fit for all lambda in the Malec sweep.
At the CNP both predictions average the n_g = -/+ n_puddle sides, as on 11-27.

Validation (exactly ONE, new code path = the diffusive junction leg inside the fit): flat profile n_c = n_b gives
R_extra = 0 exactly, so R_c reduces to the closed-form TLM sqrt(rho_c R_s) coth(L_c/L_T).
src is imported and not edited. run_full_rc.py is imported (it re-runs once and rewrites its outputs identically).
"""
import os, sys, csv, io, contextlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.constants import e

with contextlib.redirect_stdout(io.StringIO()):
    import run_full_rc as fr
bj = fr.bj
sig, n_edge, DV_AU, landauer, sharvin, r_extra_diff = bj.sig, bj.n_edge, bj.DV_AU, bj.landauer, bj.sharvin, bj.r_extra_diff
n_puddle = bj.n_puddle
r_mg_ladder, r_mg_closed = fr.r_mg_ladder, fr.r_mg_closed

PASSI_M40, PASSI_DIRAC = 519.0, 1372.0
MALEC_0, MALEC_M100 = 128.0, 95.0
CG = 12.4e-9 * 1e4  # F/m^2
LAMS = (5e-9, 15e-9, 27e-9)
LCS = (0.5e-6, 1e-6, 10e-6)
NGP = (-2e16, -5e16)

_cache = {}
def junction(model, ng, lam):
    key = (model, ng, lam)
    if key not in _cache:
        nc = n_edge(DV_AU, ng)
        if model == "ballistic":
            RL, _ = landauer(nc, ng, lam); _cache[key] = RL - sharvin(ng)
        else:
            _cache[key] = r_extra_diff(nc, ng, lam)
    return _cache[key]

def rc(model, ng, rho, lam, Lc):
    nc = n_edge(DV_AU, ng)
    return r_mg_ladder(rho, 1 / sig(nc), Lc) * 1e6 + junction(model, ng, lam)

def rc_cnp(model, rho, lam, Lc):
    return np.mean([rc(model, s * n_puddle, rho, lam, Lc) for s in (-1, 1)])

T = []; P = T.append
# ---- validation: flat profile
nb = -2e16
Rflat = r_extra_diff(nb, nb, 15e-9)
Rs = 1 / sig(nb); rho_v, Lc_v = 1e-10, 1e-6
tot = r_mg_ladder(rho_v, Rs, Lc_v) * 1e6 + Rflat; ex = r_mg_closed(rho_v, Rs, Lc_v) * 1e6
ok = Rflat == 0.0 and abs(tot / ex - 1) < 1e-4
P(f"VALIDATION flat profile (n_c = n_b = {nb*1e-4:.0e} cm^-2): diffusive R_extra = {Rflat:.3g} Ohm.um (exact 0); "
  f"R_c = {tot:.4f} vs closed-form TLM {ex:.4f} Ohm.um (rel. err {abs(tot/ex-1):.1e}) -> {'PASS' if ok else 'FAIL'}")

ng_m100 = CG * -100 / e
P(""); P(f"Junction terms (Ohm.um) at the fit and prediction densities; n_g(-100 V, Malec) = {ng_m100*1e-4:+.2e} cm^-2")
for lam in LAMS:
    s = []
    for lab, g in (("-40V/-2e12", -2e16), ("-5e12", -5e16), ("CNP p", -n_puddle), ("CNP n", n_puddle), ("Malec -100V", ng_m100)):
        s.append(f"{lab}: bal {junction('ballistic', g, lam):7.1f} / dif {junction('diffusive', g, lam):7.1f}")
    P(f"  lambda = {lam*1e9:2.0f} nm | " + " | ".join(s))

rows = []
P(""); P("model | lambda | L_c | n_g(fit) | rho_c (Ohm.cm^2) | Passi CNP/-40V (meas 2.64) | Malec 0/-100V (meas 1.35)")
for model in ("ballistic", "diffusive"):
    for lam in LAMS:
        for Lc in LCS:
            for ngp in NGP:
                f = lambda lr: rc(model, ngp, 10**lr, lam, Lc) - PASSI_M40
                try:
                    rho = 10**brentq(f, -16, -5, xtol=1e-7)
                except ValueError:
                    P(f"  {model:9s} | {lam*1e9:2.0f} nm | {Lc*1e6:4.1f} um | {ngp*1e-4:+.0e} | no rho_c > 0 reaches 519 Ohm.um"); continue
                rP = rc_cnp(model, rho, lam, Lc) / PASSI_M40
                rM = rc_cnp(model, rho, lam, Lc) / rc(model, ng_m100, rho, lam, Lc)
                rows.append(dict(model=model, lam_nm=lam * 1e9, Lc_um=Lc * 1e6, ng_fit_cm2=ngp * 1e-4, rho_c_Ohm_cm2=rho * 1e4,
                                 passi_ratio=rP, malec_ratio=rM))
                P(f"  {model:9s} | {lam*1e9:2.0f} nm | {Lc*1e6:4.1f} um | {ngp*1e-4:+.0e} | {rho*1e4:.3e} | {rP:5.2f} | {rM:5.2f}")

P(""); P("Summary per model and lambda: range over L_c and fit density; error = |log(model/measured)| (geometric-mean of the two)")
summ = {}
for model in ("ballistic", "diffusive"):
    for lam in LAMS:
        sel = [r for r in rows if r["model"] == model and r["lam_nm"] == lam * 1e9]
        if not sel: continue
        pr = np.array([r["passi_ratio"] for r in sel]); mr = np.array([r["malec_ratio"] for r in sel])
        err = np.sqrt(np.mean(np.log(pr / (PASSI_DIRAC / PASSI_M40)))**2 + np.mean(np.log(mr / (MALEC_0 / MALEC_M100)))**2)
        summ[(model, lam)] = (pr.min(), pr.max(), mr.min(), mr.max(), err)
        P(f"  {model:9s} lambda = {lam*1e9:2.0f} nm: Passi ratio {pr.min():.2f}-{pr.max():.2f} (2.64), Malec ratio {mr.min():.2f}-{mr.max():.2f} (1.35), "
          f"combined log-error {err:.3f}")
P("")
for lam in LAMS:
    if ("ballistic", lam) in summ and ("diffusive", lam) in summ:
        b, d = summ[("ballistic", lam)][4], summ[("diffusive", lam)][4]
        P(f"  lambda = {lam*1e9:2.0f} nm: closer to data = {'BALLISTIC' if b < d else 'DIFFUSIVE'} (log-error {b:.3f} vs {d:.3f})")

with open(os.path.join(HERE, "diffusive_head_to_head_results.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader()
    for r in rows: w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
open(os.path.join(HERE, "diffusive_head_to_head_transcript.txt"), "w").write("\n".join(T) + "\n")
print("\n".join(T))

fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for ax, key, meas, ttl in ((axs[0], "passi_ratio", PASSI_DIRAC / PASSI_M40, "Passi 2018: R_c(CNP)/R_c(−40 V)"),
                           (axs[1], "malec_ratio", MALEC_0 / MALEC_M100, "Malec 2011: R_c(0 V)/R_c(−100 V)")):
    for j, model in enumerate(("ballistic", "diffusive")):
        for i, lam in enumerate(LAMS):
            v = [r[key] for r in rows if r["model"] == model and r["lam_nm"] == lam * 1e9]
            x = i + (j - 0.5) * 0.3
            ax.plot([x] * len(v), v, "o" if j == 0 else "s", color=f"C{j}", ms=4,
                    label=model if i == 0 else None)
    ax.axhline(meas, color="r", ls="--", label=f"measured {meas:.2f}")
    ax.set_xticks(range(len(LAMS))); ax.set_xticklabels([f"λ = {l*1e9:.0f} nm" for l in LAMS])
    ax.set_ylabel("ratio"); ax.set_title(ttl, fontsize=9); ax.legend(fontsize=7)
fig.suptitle("I6 2026-11-28: one-parameter R_c = R_mg + junction term, ballistic vs src diffusive (ρ_c refitted per model, λ)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "diffusive_head_to_head.png"), dpi=150)
