"""I2, experiment 2 (2026-11-18): check the 2026-11-17 facet-matched result against
(a) a second DFT reference with *printed* Fermi shifts (Khomyakov et al. PRB 79,
195425 (2009), Table I) instead of the W_covered - W_G proxy used yesterday, and
(b) one experimental Au doping sign (Slawinska et al. PRB 85, 235430 (2012), STS).
Also tests the sign with the table's *experimental* metal work functions, to see
whether the sign depends on W and W0 coming from one source.
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

def rd(fn): return list(csv.DictReader(l for l in open(os.path.join(HERE, fn)) if not l.startswith("#")))
def ef(n): return -np.sign(n) * hbar * v_F * np.sqrt(np.pi * abs(n)) / e
def s(x): return "n" if x < 0 else ("p" if x > 0 else "0")

K = rd("khomyakov2009_table1.csv"); G = {r["metal"]: r for r in rd("giovannetti2008_table1.csv")}
X = rd("experimental_au_doping.csv")[0]
L = []; P = L.append

# --- the one validation (exact limit): graphene on a metal with W = W0 is undoped,
# so the model Fermi shift must be exactly 0 there.
v = ef(cdm.contact_edge_carrier_density(cdm.W_CROSSOVER))
P(f"VALIDATION dEF_model(W0={cdm.W_CROSSOVER} eV) = {v:.3e} eV (exact: 0) -> {'PASS' if v == 0 else 'FAIL'}")
P("")
P("metal | W111calc | W_exp | dEF_model(calc W) | dEF_Khomyakov printed | proxy W_cov-W_G (11-17) | dEF_model(exp W) | sign calc/exp/DFT")
out = []
for r in K:
    m = r["metal"]; Wc = float(r["W_M_111_calc_eV"]); We = float(r["W_M_exp_eV"]); ref = float(r["dEF_printed_eV"])
    prox = float(G[m]["dEF_DFT_eV"])
    mc = ef(cdm.contact_edge_carrier_density(Wc)); me = ef(cdm.contact_edge_carrier_density(We))
    out.append(dict(metal=m, W_calc=Wc, W_exp=We, dEF_model_calcW=mc, dEF_printed=ref, dEF_proxy=prox, dEF_model_expW=me,
                    R_extra_expW_ohm_um=cdm.junction_extra_resistance(We, 2.0e16)[0]))
    P(f"{m:2s} | {Wc:.2f} | {We:.2f} | {mc:+.3f} | {ref:+.2f} | {prox:+.2f} | {me:+.3f} | {s(mc)}/{s(me)}/{s(ref)}")
a = lambda k: np.array([o[k] for o in out])
mc, ref, prox, me = a("dEF_model_calcW"), a("dEF_printed"), a("dEF_proxy"), a("dEF_model_expW")
rms = lambda x, y: np.sqrt(np.mean((x - y) ** 2))
P("")
P(f"signs model(calc W) vs Khomyakov printed: {int(np.sum(np.sign(mc)==np.sign(ref)))}/5; RMS {rms(mc,ref):.3f} eV")
P(f"yesterday's proxy vs printed: RMS {rms(prox,ref):.3f} eV (proxy is not the Fermi shift: W_cov-W_G includes the interface dipole)")
P(f"model(calc W) vs proxy (yesterday's number): RMS {rms(mc,prox):.3f} eV")
P(f"signs model(exp W, W0=5.4 calc) vs printed: {int(np.sum(np.sign(me)==np.sign(ref)))}/5; RMS {rms(me,ref):.3f} eV")
au = [o for o in out if o["metal"] == "Au"][0]
lo, hi = float(X["dEF_min_eV"]), float(X["dEF_max_eV"])
P(f"Au experiment (STS, Slawinska 2012): dEF in [{lo:+.2f}, {hi:+.2f}] eV, p-type")
P(f"  model calc W: {au['dEF_model_calcW']:+.3f} eV sign {s(au['dEF_model_calcW'])} -> sign {'AGREES' if au['dEF_model_calcW']>0 else 'DISAGREES'}; inside range: {lo<=au['dEF_model_calcW']<=hi}")
P(f"  model exp  W: {au['dEF_model_expW']:+.3f} eV sign {s(au['dEF_model_expW'])} -> sign {'AGREES' if au['dEF_model_expW']>0 else 'DISAGREES'}")
P(f"  Khomyakov DFT +0.19 inside range: {lo<=0.19<=hi}")
with open(os.path.join(HERE, "second_reference_results.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
open(os.path.join(HERE, "second_reference_transcript.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))

fig, ax = plt.subplots(figsize=(6.4, 4.6))
Ws = np.linspace(4.0, 6.3, 200)
ax.plot(Ws, [ef(cdm.contact_edge_carrier_density(W)) for W in Ws], "C0-", label="repo model (src/, unchanged), W0 = 5.4 eV")
ax.plot(a("W_calc"), ref, "ks", label="Khomyakov 2009 Table I, printed ΔE_F (DFT)")
ax.plot(a("W_calc"), prox, "o", mfc="none", c="gray", label="W_cov − W_G proxy (used 2026-11-17)")
ax.plot(a("W_exp"), me, "C3^", mfc="none", label="model at experimental W_M (same table)")
ax.fill_between([5.25, 5.6], lo, hi, color="C2", alpha=.25, label="Au STS, Sławińska 2012 (poly Au)")
for o in out: ax.annotate(o["metal"], (o["W_calc"], o["dEF_printed"]), xytext=(4, 4), textcoords="offset points")
ax.axhline(0, c="gray", lw=.6); ax.axvline(cdm.W_CROSSOVER, c="gray", ls=":", lw=.8)
ax.set_xlabel(r"metal work function $W_M$ (eV)"); ax.set_ylabel(r"$\Delta E_F$ (eV), n-type < 0")
ax.legend(fontsize=7.5, loc="upper left"); ax.set_title("I2: second DFT reference + experimental Au sign")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "second_reference.png"), dpi=150)
