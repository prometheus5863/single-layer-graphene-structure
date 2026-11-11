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
    dict(title="(b) thermal conductivity of the oxide",
         module="graphene_self_heating_model:327",
         token="1.2794", quoted=1.2794, got=kappa, tol=5e-4, rel=False,
         unit=r"$\kappa_{\rm ox}(300\,$K$)$  (W/m$\cdot$K)",
         verdict="VACUOUS"),
    dict(title="(c) Dorgan's total thermal resistance",
         module="graphene_self_heating_model:330",
         token="1e4", quoted=1.0e4, got=R_th, tol=0.05, rel=True,
         unit=r"$R_{\rm th}$  (K/W)", verdict="OVER_TIGHT"),
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
    ax.text(q, 1.08, "source prints %s" % row["token"], ha="center",
            fontsize=8.6, color="0.25")

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
print("  (c) R_th     computed %.1f, licensed [%.0f, %.0f], tested [%.0f, %.0f]"
      % (R_th, 5e3, 1.5e4, 9.5e3, 1.05e4))
