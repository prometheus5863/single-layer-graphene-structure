"""I6, experiment 1 (2026-11-26): ballistic (Landauer, k_y-resolved Dirac transfer matrix) resistance of the
Au contact junction on I1's gate-tuned edge, against the src diffusive R_extra on the same profile, and the
saturation test of the p/n asymmetry against Huard et al. PRB 78, 121402 (2008), arXiv:0804.2040 (text values).

Profile (I1, 11-23/11-24, unchanged): contact-edge density n_c = n_edge(W - W0 = 0.59 eV, n_g) from src's
C_q-in-series-with-C_I equation with the gate charge included; channel density n_b = n_g; one tail
n(x) = n_b + (n_c - n_b) / (1 + x / lambda) for x >= 0. Graphene under the metal (x < 0) is the left lead at n_c.
The tail is truncated at X = max(4 um, 200 lambda) and joined abruptly to the n_b lead (residual step reported).

Transport: massless Dirac equation at the Fermi level, k_y conserved, local Fermi wavenumber
eps(x) = sign(n) sqrt(pi |n|) (eps = (E_F - U)/hbar v_F). d psi/dx = A psi with A = i eps sigma_x + k_y sigma_z,
integrated exactly slice by slice (piecewise constant, exp(A dx) = cosh(kap dx) I + sinh(kap dx)/kap A,
kap^2 = k_y^2 - eps^2); this has no singularity at the p-n neutrality point. Lead modes are the eigenvectors of A,
classified and normalised by the flux psi^dag sigma_x psi.
  G/W = (4 e^2/h) int dk_y/(2 pi) T(k_y),   R_L * W = 1 / (G/W)
  R_J = R_L - R_Sh(n_b),  R_Sh = (h/4e^2) pi / k_F(n_b)   (junction excess over the channel's own Sharvin term)
Diffusive comparison: R_extra = int_0^X [1/sigma(n(x)) - 1/sigma(n_b)] dx with src's puddle-regularised sigma.

Validation (exactly ONE, exact limit): flat profile n_c = n_b gives T(k_y) = 1 for every propagating mode and
R_L * W = Sharvin value (h/4e^2) pi / k_F.
src is imported and not edited.
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
from graphene_fet_model import v_F, n_puddle, mu
sig = cdm._sheet_conductivity_puddle_regularized
RQ = h / (4 * e**2)

# copied verbatim from I1/run_tlm_gate_tuned.py (11-24); not imported, because importing re-runs it
def EF(n): return np.sign(n) * hbar * v_F * np.sqrt(np.pi * np.abs(n)) / e
def n_edge(dV, n_g, CI=cdm.C_interface):
    gfun = lambda n: EF(n) + e * (n - n_g) / CI - dV
    span = 1e18 + abs(n_g) + abs(CI * dV / e) * 2
    return brentq(gfun, -span, span, xtol=1e6, rtol=1e-13)
DV_AU = -0.59
X_END = 4e-6

def keps(n): return np.sign(n) * np.sqrt(np.pi * np.abs(n))
SX = np.array([[0, 1], [1, 0]], complex); SZ = np.array([[1, 0], [0, -1]], complex)

def xend(lam): return max(X_END, 200 * lam)

def grid(lam, npts=2500):
    return np.concatenate([[0.0], np.geomspace(lam * 1e-3, xend(lam), npts)])

def lead_modes(eps, ky):
    A = 1j * eps * SX + ky * SZ
    w, V = np.linalg.eig(A)
    out = {}
    for i in range(2):
        v = V[:, i]; fl = np.real(np.conj(v) @ SX @ v)
        out["R" if fl > 0 else "L"] = v / np.sqrt(abs(fl))
    return out

def transmission(nprof, xs, kys, epsL, epsR):
    """T(k_y), batched over k_y. The purely transmitted state t*u_R(right lead) is propagated backwards from X to 0
    with the exact slice propagators (det = 1), renormalised every slice with the log-norm accumulated, so evanescent
    growth cannot overflow. At x = 0 it is split into lead modes, psi(0) = a u_R + b u_L, and T = |t|^2 = 1/|a|^2."""
    kys = np.atleast_1d(np.asarray(kys, float))
    epsm = keps(0.5 * (nprof[1:] + nprof[:-1])); dx = np.diff(xs)
    psi = np.array([lead_modes(epsR, k)["R"] for k in kys])          # (nk, 2)
    logs = np.zeros(len(kys))
    for ep, d in zip(epsm[::-1], dx[::-1]):
        A = 1j * ep * SX[None] + kys[:, None, None] * SZ[None]
        kap = np.sqrt((kys**2 - ep**2).astype(complex))
        kd = kap * d
        s = np.where(np.abs(kd) > 1e-12, np.sinh(kd) / np.where(kap == 0, 1, kap), d)
        Pinv = np.cosh(kd)[:, None, None] * np.eye(2)[None] - s[:, None, None] * A   # exp(-A d)
        psi = np.einsum("kij,kj->ki", Pinv, psi)
        nrm = np.linalg.norm(psi, axis=1); psi /= nrm[:, None]; logs += np.log(nrm)
    Ts = np.empty(len(kys))
    for i, k in enumerate(kys):
        mL = lead_modes(epsL, k)
        a_, b_ = np.linalg.solve(np.column_stack([mL["R"], mL["L"]]), psi[i])
        Ts[i] = np.exp(-2 * logs[i]) / abs(a_) ** 2
    return Ts

def landauer(n_c, n_b, lam, nk=161):
    xs = grid(lam); nprof = n_b + (n_c - n_b) / (1 + xs / lam)
    epsL, epsR = keps(n_c), keps(n_b)
    kmax = min(abs(epsL), abs(epsR))
    if kmax == 0: return np.inf, 0.0
    # Gauss-Legendre in phi with ky = kmax sin(phi), to resolve the band edges
    ph, wt = np.polynomial.legendre.leggauss(nk); ph = ph * np.pi / 2; wt = wt * np.pi / 2
    ky = kmax * np.sin(ph) * (1 - 1e-12); jac = kmax * np.cos(ph)
    Ts = transmission(nprof, xs, ky, epsL, epsR)
    GW = (1 / RQ) * np.sum(wt * jac * Ts) / (2 * np.pi)
    resid = abs(nprof[-1] - n_b) / abs(n_b)
    return 1 / GW * 1e6, resid          # Ohm.um

def r_extra_diff(n_c, n_b, lam, npts=6000):
    xs = np.concatenate([[0.0], np.geomspace(1e-13, xend(lam), npts)])
    nprof = n_b + (n_c - n_b) / (1 + xs / lam)
    return np.trapezoid(1 / sig(nprof) - 1 / sig(n_b), xs) * 1e6

def sharvin(n_b): return RQ * np.pi / np.sqrt(np.pi * abs(n_b)) * 1e6

T = []; P = T.append
# ---- validation: flat profile
nb = 3e16
xs = grid(27e-9); flat = np.full_like(xs, nb)
kmax = keps(nb)
Tmin = transmission(flat, xs, np.linspace(0, 0.999 * kmax, 25), kmax, kmax).min()
RL = landauer(nb, nb, 27e-9)[0]
P(f"VALIDATION flat profile n = 3e12 cm^-2: min T(k_y) = {Tmin:.12f} (exact 1); R_L W = {RL:.4f} Ohm.um vs "
  f"Sharvin {sharvin(nb):.4f} (rel. err {abs(RL/sharvin(nb)-1):.1e}) -> {'PASS' if abs(Tmin-1)<1e-9 and abs(RL/sharvin(nb)-1)<1e-3 else 'FAIL'}")
l_mfp = hbar * np.sqrt(np.pi * 2e16) * mu / e
P(f"src mean free path l = hbar k_F mu / e at 2e12 cm^-2 = {l_mfp*1e9:.0f} nm (src mu = {mu*1e4:.0f} cm^2/Vs); "
  f"src puddle density = {n_puddle*1e-4:.1e} cm^-2")
P("")

LAMS = [5e-9, 15e-9, 27e-9, 250e-9]
NS = np.array([0.5, 0.8, 1.2, 1.6, 2, 2.5, 3, 4, 5, 6, 7, 8]) * 1e16
rows = []
P("lambda | |n| | branch | n_c (cm^-2) | R_L | R_Sh | R_J = R_L - R_Sh | R_extra diffusive | R_J / R_extra | tail residual")
for lam in LAMS:
    for nabs in NS:
        for br, sgn in (("p", -1), ("n", +1)):
            ng = sgn * nabs; nc = n_edge(DV_AU, ng)
            RL, res = landauer(nc, ng, lam); RS = sharvin(ng); RJ = RL - RS
            RD = r_extra_diff(nc, ng, lam)
            rows.append(dict(lam_nm=lam * 1e9, n_abs_cm2=nabs * 1e-4, branch=br, n_c_cm2=nc * 1e-4, R_L=RL, R_Sh=RS,
                             R_J=RJ, R_extra_diff=RD, ratio_RJ_Rextra=RJ / RD if RD else np.nan, tail_resid=res))
            P(f"{lam*1e9:5.0f} | {nabs*1e-4:.1e} | {br} | {nc*1e-4:+.2e} | {RL:8.1f} | {RS:6.1f} | {RJ:8.1f} | {RD:8.1f} | "
              f"{RJ/RD if RD else float('nan'):6.2f} | {res:.3f}")
P("")

def get(lam, br, key): return np.array([r[key] for r in rows if r["lam_nm"] == lam * 1e9 and r["branch"] == br])
na = NS * 1e-4; sel = (na >= 1.2e12) & (na <= 7e12)
P("Saturation test (Huard 2008 H2: R_odd saturates for n > 1.2e12 up to ~7e12; slope 0 = saturated):")
summ = []
for lam in LAMS:
    dB = get(lam, "n", "R_L") - get(lam, "p", "R_L"); dD = get(lam, "n", "R_extra_diff") - get(lam, "p", "R_extra_diff")
    sB = np.polyfit(np.log(na[sel]), np.log(np.abs(dB[sel])), 1)[0] if np.all(dB[sel] > 0) else np.nan
    sD = np.polyfit(np.log(na[sel]), np.log(np.abs(dD[sel])), 1)[0] if np.all(dD[sel] > 0) else np.nan
    rJ = get(lam, "n", "ratio_RJ_Rextra")[sel]     # p-n branch, Huard window; p-p' branch reported separately
    summ.append((lam, sB, sD, np.nanmin(rJ), np.nanmax(rJ), dB, dD))
    P(f"  lambda = {lam*1e9:3.0f} nm: ballistic dR = R_L(n) - R_L(p): {dB[sel][0]:.0f} -> {dB[sel][-1]:.0f} Ohm.um "
      f"(sign > 0 everywhere in 1.2-7e12: {bool(np.all(dB[sel] > 0))}), slope {sB:.2f}; diffusive slope {sD:.2f}; "
      f"n-branch (p-n) R_J/R_extra over 1.2-7e12: {np.nanmin(rJ):.2f}-{np.nanmax(rJ):.2f}; "
      f"p-branch (p-p'): R_J = {get(lam,'p','R_J')[sel].min():.2f}...{get(lam,'p','R_J')[sel].max():.2f}, diffusive R_extra = "
      f"{get(lam,'p','R_extra_diff')[sel].min():.0f}...{get(lam,'p','R_extra_diff')[sel].max():.0f} Ohm.um")
P("")
p40 = [r for r in rows if r["branch"] == "p" and r["lam_nm"] == 27 and r["n_abs_cm2"] in (2e12, 5e12)]
P("Passi 2018 (arXiv:1807.04772, Table 1, Au, -40 V): R_c = 519 Ohm.um; I1 bracket n = -2 ... -5e12 cm^-2 at lambda = 27 nm:")
for r in p40:
    P(f"  n = -{r['n_abs_cm2']:.0e}: ballistic R_J = {r['R_J']:.0f} (R_L = {r['R_L']:.0f}), diffusive R_extra = {r['R_extra_diff']:.0f} Ohm.um")
P("")
P("Metric definition: R_J/R_extra on the n-branch (p-n junction, both terms > 0), |n| = 1.2-7e12 cm^-2. On the p-branch the ratio is"
  " not defined: src's diffusive R_extra < 0 (doped tail conducts better) while Landauer R_J >= 0 (channel Sharvin is the bottleneck).")
kill = all(0.5 <= s[3] and s[4] <= 2 for s in summ)
succ1 = any(s[3] < 0.5 or s[4] > 2 for s in summ if s[0] <= 15e-9)
P(f"I6 kill condition (R_J/R_extra within 0.5-2x for all lambda 5-250 nm): {'HIT' if kill else 'not hit'}")
P(f"I6 success metric part 1 (>2x difference at lambda 5-15 nm): {'met' if succ1 else 'not met'}")

with open(os.path.join(HERE, "ballistic_junction_results.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader()
    for r in rows: w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
open(os.path.join(HERE, "ballistic_junction_transcript.txt"), "w").write("\n".join(T) + "\n")
print("\n".join(T))

fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
for i, (lam, sB, sD, lo, hi, dB, dD) in enumerate(summ):
    c = f"C{i}"
    axs[0].plot(na, get(lam, "n", "R_J"), c + "o-", ms=3, label=f"ballistic R_J, λ={lam*1e9:.0f} nm")
    axs[0].plot(na, get(lam, "n", "R_extra_diff"), c + "--", label=f"diffusive R_extra, λ={lam*1e9:.0f} nm")
    axs[1].plot(na, dB, c + "o-", ms=3, label=f"ballistic ΔR, λ={lam*1e9:.0f} nm (slope {sB:.2f})")
    axs[1].plot(na, dD, c + "--", label=f"diffusive ΔR, λ={lam*1e9:.0f} nm (slope {sD:.2f})")
axs[0].set_yscale("log"); axs[0].set_xlabel("n (cm⁻²), electron branch (p–n junction at the Au edge)"); axs[0].set_ylabel("Ω·µm per contact")
axs[0].legend(fontsize=6); axs[0].set_title("Junction term, n-branch (p-branch: ballistic R_J < 0.2 Ω·µm, diffusive < 0)", fontsize=8)
axs[1].plot(na, 2.06 * RQ / np.sqrt(NS) * 1e6, "k:", label="Huard Eq. 3 sharp p–n bound")
axs[1].axvspan(1.2e12, 7e12, color="gray", alpha=.12, label="Huard: R_odd saturated here")
axs[1].set_yscale("log"); axs[1].set_xlabel("|n| (cm⁻²)"); axs[1].set_ylabel("ΔR = R(n) − R(p), Ω·µm")
axs[1].legend(fontsize=6); axs[1].set_title("p/n asymmetry vs Huard 2008 (text constraints)", fontsize=9)
fig.suptitle("I6 2026-11-26: ballistic (Dirac transfer matrix) vs src diffusive, Au W−W0 = 0.59 eV", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "ballistic_junction.png"), dpi=150)
