"""
graphene_tolerance_census_figure.py

What the source licenses, against what the check tests.

Three panels, one per adjudicated row of graphene_tolerance_census_audit.py on
2026-11-11.  Each shows, on the quantity's own axis:

  * the LICENSED interval -- the set of values that round to the number the
    source actually printed, at the number of figures the source printed it to;
  * the TESTED interval -- the set of values the repository's check admits;
  * the computed value.

A check is sound when the tested interval is the licensed one.  Wider and it
admits values that round to a different published figure (VACUOUS); narrower
and it can reject a value the source's own rounding allows (OVER_TIGHT), which
is the defect found on 2026-11-10.

PANEL (b) IS DRAWN AS THE CENSUS FOUND IT AND IS NOT THE FINAL STATE.  Reading
Dorgan et al. afterwards showed that 1.2794 W/m.K is not a value the paper
prints at any precision: supp. Sec. 6 gives the fit ln(T_ox^0.52) - 1.687, and
the only numerical conductivity in the paper is supp. Sec. 5's "kappa_ox ~ 1.3
Wm-1K-1 ... (at 300 K)".  Two significant figures license +/-0.05, a band 1000
times wider than the +/-5e-5 the phantom five-figure value appeared to license
and far too wide to draw on the same axis, so it is annotated instead.  The
computed 1.2789669 rounds to 1.3, so the transcription was sound and only its
warrant was not.  The panel is kept in the as-found state because the defect it
shows -- a tolerance adjudicated against a precision no source ever claimed --
is the finding.

Run: python3 graphene_tolerance_census_figure.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import graphene_fet_signed_carrier_model as fsc
import graphene_self_heating_model as sh
from graphene_tolerance_census_audit import half_width

OUT = "tolerance_census_licensed_vs_tested.png"

# --- the three rows, with their computed values recomputed here ------------
R_sheet = 1e-3 / fsc.SIGMA_MIN_CH3                       # kOhm/sq
kappa = sh.kappa_ox(300.0)                               # W/m.K
R_th = sh.thermal_resistances(W=7e-6, L=4e-6, t_ox=300e-9,
                              T=sh.T_AMBIENT, T0=sh.T_AMBIENT)["R_th"]

ROWS = (
    dict(title="(a) Chapter 3's sheet resistance floor",
         module="graphene_fet_signed_carrier_model:648",
         token="6.45", quoted=6.45, got=R_sheet, tol=0.01, rel=False,
         unit=r"$R_{\rm sheet}$  (k$\Omega$/sq)", verdict="VACUOUS"),
    dict(title="(b) oxide thermal conductivity vs a value with no source",
         module="graphene_self_heating_model:A1, as found",
         token="1.2794", quoted=1.2794, got=kappa, tol=5e-4, rel=False,
         unit=r"$\kappa_{\rm ox}(300\,$K$)$  (W/m$\cdot$K)",
         verdict="VACUOUS"),
    # PANEL (c) UPDATED 2026-11-12.  A2's tolerance was +-5 % when this
    # figure was drawn; it is now +-1.7857 %, tightened rather than widened,
    # because the paper's SAME SENTENCE prints the fact to two figures as
    # ~2.8e-7 m^2.K/W per unit area.  The panel is kept because the census
    # still returns OVER_TIGHT on it -- but that verdict is now a FALSE
    # POSITIVE (census check 9l): the instrument reads the one-figure centre
    # token "1e4" and cannot see that the width's warrant is the two-figure
    # areal number.  The panel therefore illustrates a limitation of the
    # instrument, which is a different thing from what it illustrated
    # yesterday, and its title says so.
    dict(title="(c) Dorgan's total thermal resistance -- the tested band is "
               "now the 2 s.f. areal one; OVER_TIGHT here is a FALSE POSITIVE",
         module="graphene_self_heating_model:A2, repaired 2026-11-12",
         token="1e4", quoted=1.0e4, got=R_th, tol=0.05 / 2.8, rel=True,
         unit=r"$R_{\rm th}$  (K/W)", verdict="OVER_TIGHT (false positive)"),
)

C_LIC = "#1b6ca8"
C_TEST = "#c1433c"

fig, axes = plt.subplots(3, 1, figsize=(9.4, 8.2))

for ax, row in zip(axes, ROWS):
    q = row["quoted"]
    h_abs, h_rel, sf = half_width(row["token"])
    lic_lo, lic_hi = q - h_abs, q + h_abs
    t_abs = row["tol"] * q if row["rel"] else row["tol"]
    test_lo, test_hi = q - t_abs, q + t_abs

    span = max(h_abs, t_abs)
    ax.set_xlim(q - 1.45 * span, q + 1.45 * span)
    ax.set_ylim(-0.95, 1.35)

    ax.broken_barh([(test_lo, test_hi - test_lo)], (0.05, 0.42),
                   facecolors=C_TEST, alpha=0.22, edgecolors=C_TEST,
                   linewidth=1.4)
    ax.broken_barh([(lic_lo, lic_hi - lic_lo)], (-0.52, 0.42),
                   facecolors=C_LIC, alpha=0.26, edgecolors=C_LIC,
                   linewidth=1.4)

    ax.axvline(q, color="0.25", linewidth=1.1, linestyle="--", zorder=1)
    ax.plot([row["got"]], [-0.31], marker="o", markersize=7.5,
            color="black", zorder=5)

    ax.text(test_hi, 0.26, "  tested:  $\\pm$%s%s" % (
        ("%g" % (row["tol"] * 100)) + " %" if row["rel"]
        else "%g" % row["tol"],
        "" if row["rel"] else ""),
        va="center", ha="left", fontsize=8.6, color=C_TEST)
    ax.text(lic_hi, -0.31, "  licensed by %d s.f.:  $\\pm$%.3g"
            % (sf, h_abs), va="center", ha="left", fontsize=8.6, color=C_LIC)
    ax.text(q, 1.08, "label quotes %s" % row["token"], ha="center",
            fontsize=8.6, color="0.25")
    if row["module"].startswith("graphene_self_heating_model:A1"):
        ax.annotate("the paper prints only ~1.3 W/m$\\cdot$K (supp. Sec. 5):\n"
                    "2 s.f. license $\\pm$0.05, a band 1000x wider, off scale.\n"
                    "1.2794 is in no part of the paper at any precision.",
                    xy=(0.015, 0.055), xycoords="axes fraction", ha="left",
                    va="center", fontsize=8.2, color="#6a4c93")

    ratio = h_abs / t_abs
    ax.set_title("%s   --   %s, %s by %.3gx"
                 % (row["title"], row["module"], row["verdict"],
                    (1.0 / ratio) if ratio < 1 else ratio),
                 fontsize=9.6, loc="left")
    ax.set_xlabel(row["unit"], fontsize=9)
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x", labelsize=8.4)

axes[0].text(0.995, 1.18, "computed value", transform=axes[0].transAxes,
             ha="right", fontsize=8.4, color="black")

fig.suptitle("What the source licenses, against what the check tests\n"
             "graphene_tolerance_census_audit.py, 2026-11-11",
             fontsize=11.5)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(OUT, dpi=160)
print("wrote %s" % OUT)
print("  (a) R_sheet  computed %.6f, licensed [%.4f, %.4f], tested [%.4f, %.4f]"
      % (R_sheet, 6.45 - 0.005, 6.45 + 0.005, 6.44, 6.46))
print("  (b) kappa_ox computed %.6f, licensed [%.5f, %.5f], tested [%.4f, %.4f]"
      % (kappa, 1.2794 - 5e-5, 1.2794 + 5e-5, 1.2794 - 5e-4, 1.2794 + 5e-4))
print("  (c) R_th     computed %.1f, licensed [%.0f, %.0f] from the "
      "one-figure ~1e4, tested [%.1f, %.1f] from the two-figure ~2.8e-7"
      % (R_th, 5e3, 1.5e4, 1.0e4 * (1 - 0.05 / 2.8), 1.0e4 * (1 + 0.05 / 2.8)))
print("      the OVER_TIGHT verdict on this row is a FALSE POSITIVE as of")
print("      2026-11-12: centre and width come from different printed")
print("      numbers in one sentence, which the census cannot represent.")
