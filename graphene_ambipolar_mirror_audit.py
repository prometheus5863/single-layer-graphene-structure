"""
graphene_ambipolar_mirror_audit.py

THE SYMMETRY WAS NEVER LOST AT FINITE DRAIN BIAS.  ITS CENTRE MOVED, BY
EXACTLY V_ds/2, AND TWO PLACES IN THIS REPOSITORY STILL READ IT AT V_dirac.

Provenance.  The 2026-11-15 entry ranked second an item created by its own
Section 2: "AUDIT THE REPOSITORY FOR OTHER UNIPOLAR FORMS TRANSCRIBED TO AN
AMBIPOLAR DEVICE", and recorded a prediction to be scored --

    "The |g_m| audit of Section 2 will find at least one more instance, and
     it will be in a figure-of-merit expression rather than in a current
     model.  Reason: current models in this repository take abs(n) early
     (sheet_conductivity does), while figures of merit are transcribed from
     unipolar RF literature where the sign never arises."

This module is that audit.  It does not find a missing `abs`.  It finds the
weaker and more interesting form of the same fault: two quantities that are
each individually even about a reflection centre, divided by one another,
with DIFFERENT centres -- and the mismatch is exactly the drain bias that a
unipolar treatment never has to reflect about.

The law, stated first so the measurement can refute it.
----------------------------------------------------------------------------
`gfet.transfer_characteristic(Vg, Vds)` averages the local channel
resistance over `V_channel_profile = np.linspace(0, Vds, 50)`:

    R_bar(V_g) = mean_{V_ch in profile} R_ch(V_g - V_ch - V_dirac)

and `R_ch` depends on its argument only through `|.|` (carrier_density
returns sqrt(n_eff^2 + n_puddle^2); sheet_conductivity takes abs).  Write
u = V_g - V_dirac.  Then

    R_bar(u) = mean_{V_ch} f(|u - V_ch|),   f even,

and the profile is itself symmetric about V_ds/2 (linspace with both
endpoints included).  Substituting V_ch -> V_ds - V_ch maps
u - V_ch -> u - V_ds + V_ch, so

    R_bar(V_ds - u) = R_bar(u)          i.e.   R_bar is EVEN about u = V_ds/2

    ==>  V_MIRROR = V_dirac + V_ds/2.

Consequences, all exact rather than asymptotic:
  * I_d is EVEN about V_MIRROR at every V_ds, not only as V_ds -> 0;
  * g_m = dI_d/dV_g is ODD about V_MIRROR;
  * the gate voltage of the current minimum is V_dirac + V_ds/2, which is
    the gate voltage an experimentalist reports as "the Dirac point";
  * `rf_small_signal_model.gate_capacitance(Vg, Vch=0.0)` is even about
    V_dirac -- BITWISE, since its argument is literally Vg - Vch - V_dirac
    passed through the even `quantum_capacitance` -- so f_T = |g_m|/(2 pi
    C_gs) is even about NEITHER centre.

What this corrects, and what it does not.
----------------------------------------------------------------------------
This is a strengthening of a result already in the tree, not a refutation of
one.  `graphene_low_frequency_noise_model.py` (2026-11-15) states in its V2
that "at finite Vds the asymmetry is the drain bias itself", and obtains its
exact symmetry check by evaluating a SEPARATE V_ds -> 0 limit function,
`id_zero_vds_limit`.  That reading is correct and its check is sound.  What
it does not say -- and what is measured here -- is that nothing is lost at
finite V_ds: the exact evenness survives, displaced.  A check that needed a
limit can be taken at the operating bias instead.  No number in that module
is withdrawn; its V2b convergence rate (asymmetry linear in V_ds, measured
ABOUT V_dirac) remains exactly right about the quantity it measures.

Literature.  The half-drain-bias centre is not this model's artefact.
Jimenez, "Explicit drain current, charge and capacitance model of graphene
field-effect transistors", arXiv:1105.2734 (DOI 10.1109/TED.2011.2168960),
names

    V_gs,D = V_gs0 + V_ds/2

the "Dirac voltage", and says of it that "electron and hole charges are
exactly balanced over the channel length" -- a statement about the channel
INTEGRAL, which is why the centre is the midpoint of the channel potential
and not the source end.  The same paper takes the gate charge as
Q_g = -W * integral_0^L Q_c(y) dy, i.e. the gate capacitance of a GFET at
finite drain bias is a channel-integrated quantity; evaluating it at V_ch=0
is the approximation this module measures.  See
notes/2026-11-16-the-mirror-centre-of-an-ambipolar-transfer-curve.md.

Outputs: ambipolar_mirror_audit_output.txt, ambipolar_mirror_audit.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

from bracketed_root import bracketed_bisect

import graphene_fet_model as gfet
import rf_small_signal_model as rf

VDS_RF = 0.1          # the bias Chapter 4's RF numbers are quoted at
N_SEGMENTS = 50       # must match transfer_characteristic's literal
VDS_DECADES = (0.001, 0.01, 0.05, 0.1, 0.2, 0.5)


# ------------------------------------------------------------------ the law
def v_mirror(Vds):
    """The exact reflection centre of I_d(V_g) at drain bias Vds."""
    return gfet.V_dirac + 0.5 * Vds


def channel_profile(Vds, n_segments=N_SEGMENTS):
    """The same grid transfer_characteristic() integrates over."""
    return np.linspace(0.0, Vds, n_segments)


def gate_capacitance_channel_averaged(Vg_range, Vds=VDS_RF,
                                      n_segments=N_SEGMENTS):
    """C_gs averaged over the channel potential, on the SAME grid the drain
    current is integrated over.

    This is the repair, and it is a repair of BOOKKEEPING rather than of
    physics: `rf.gate_capacitance` is called unchanged, once per segment,
    and the segment values are averaged.  Two properties make it the right
    object to compare g_m against:

      (i)  it is even about V_MIRROR, because the profile is symmetric and
           rf.gate_capacitance is even in its own argument -- so f_T built
           from it is even, as an ambipolar device's f_T must be;
      (ii) at Vds = 0 it reduces to rf.gate_capacitance(Vg, Vch=0.0)
           BITWISE: the profile collapses to 50 copies of 0.0, and the mean
           of N identical floats is that float.  V8 pins this.

    `rf.gate_capacitance` is deliberately NOT changed, so every committed
    number stays reproducible and the artefact stays measurable.
    """
    prof = channel_profile(Vds, n_segments)
    return np.mean([rf.gate_capacitance(Vg_range, Vch=v) for v in prof],
                   axis=0)


def _mirror_pairs(Vds, centre, d):
    """I_d at centre +- d, as a pair of arrays."""
    return (gfet.transfer_characteristic(centre + d, Vds=Vds),
            gfet.transfer_characteristic(centre - d, Vds=Vds))


def mirror_asymmetry(Vds, centre, d):
    """max |I(c+d) - I(c-d)| / mean, the relative mirror asymmetry."""
    hi, lo = _mirror_pairs(Vds, centre, d)
    return float(np.max(np.abs(hi - lo) / (0.5 * np.abs(hi + lo))))


def locate_current_minimum(Vds, span=0.45):
    """The gate voltage of the I_d minimum, by Brent rather than by argmax
    over a grid, so the answer is not quantised by the sweep.

    Reported for the record; the LAW is tested by locate_mirror_centre
    below instead, because a minimiser localises an extremum only to
    ~sqrt(eps) of the curvature scale and that floor (measured at 1.1e-6 in
    V_g) is the minimiser's, not the model's.  Stating a 1e-6 agreement as
    though it were the law's accuracy would be reporting a tolerance as a
    result.
    """
    r = minimize_scalar(
        lambda v: float(gfet.transfer_characteristic(np.array([v]),
                                                     Vds=Vds)[0]),
        bracket=(gfet.V_dirac - span, gfet.V_dirac + span + Vds))
    return float(r.x)


def locate_mirror_centre(Vds, h=0.3, span=0.6):
    """The reflection centre of I_d, located as the ROOT of the odd function

        F(c) = I_d(c + h) - I_d(c - h)

    rather than as the minimum of an even one.  A root of a function with
    non-zero slope is localisable to machine precision; an extremum is not.
    This is the same distinction graphene_low_frequency_noise_model.py's
    threshold bisection ran into from the other side (there the predicate
    was degenerate, so a root-find was NOT available and a predicate bisect
    was used instead).
    """
    def F(c):
        hi = float(gfet.transfer_characteristic(np.array([c + h]),
                                                Vds=Vds)[0])
        lo = float(gfet.transfer_characteristic(np.array([c - h]),
                                                Vds=Vds)[0])
        return hi - lo
    lo_b, hi_b = gfet.V_dirac - span, gfet.V_dirac + span + Vds
    return bracketed_bisect(F, lo_b, hi_b, tol=1e-13)


# ------------------------------------------------- figures of merit at bias
def fT_both_ways(Vg, Vds=VDS_RF):
    """f_T from the committed source-end C_gs and from the channel-averaged
    one.  Everything else (g_m, the estimator, the grid) is identical, so
    the two differ ONLY in which reflection centre the capacitance has."""
    gm, Id = rf.transconductance(Vg, Vds=Vds)
    C_src = rf.gate_capacitance(Vg)
    C_avg = gate_capacitance_channel_averaged(Vg, Vds=Vds)
    return (np.abs(gm) / (2 * np.pi * C_src),
            np.abs(gm) / (2 * np.pi * C_avg), gm, Id, C_src, C_avg)


def branch_peaks(Vg, f, centre):
    """(V_g, value) of the peak on each side of `centre`."""
    hi = Vg > centre
    e_v, e_f = Vg[hi][int(np.argmax(f[hi]))], float(np.max(f[hi]))
    h_v, h_f = Vg[~hi][int(np.argmax(f[~hi]))], float(np.max(f[~hi]))
    return (float(e_v), e_f), (float(h_v), h_f)


def mirror_profile(Vg_half, f_fn, centre):
    """Relative mirror asymmetry of a callable f_fn(Vg) about `centre`,
    as a function of distance from it."""
    hi, lo = f_fn(centre + Vg_half), f_fn(centre - Vg_half)
    return np.abs(hi - lo) / (0.5 * np.abs(hi + lo))


# ------------------------------------------------------------- validation
def validate():
    """V1, V3, V4, V5, V8, V9 are against EXACTLY known values (an exact
    reflection symmetry, a bitwise identity, an exact ratio of 1, an exact
    linear law in Vds).  V2 and V7 are EXISTENCE checks on the artefact --
    they must FAIL to hold, i.e. they assert a departure is measurable, and
    are labelled so.  V6 is a convergence rate, which is weaker than both.
    """
    out, ok = [], True

    def check(name, detail, passed, exact=True):
        nonlocal ok
        ok &= bool(passed)
        out.append((name, detail, bool(passed), exact))

    d = np.linspace(0.02, 2.0, 41)

    # V1  EXACT.  I_d is even about V_dirac + Vds/2 at EVERY drain bias.
    worst, worst_vds = 0.0, None
    for Vds in VDS_DECADES:
        a = mirror_asymmetry(Vds, v_mirror(Vds), d)
        if a > worst:
            worst, worst_vds = a, Vds
    check("V1 EXACT  I_d is even about V_dirac+Vds/2 at all six drain biases",
          "worst rel asymmetry %.3e (at Vds = %g V), over Vds spanning "
          "a factor 500" % (worst, worst_vds), worst < 1e-13)

    # V2  THE ARTEFACT EXISTS.  About V_dirac the same curve is NOT even,
    #     and the departure is large compared with V1's floor.
    a_dirac = mirror_asymmetry(VDS_RF, gfet.V_dirac, d)
    check("V2 ARTEFACT  the same curve is NOT even about V_dirac",
          "rel asymmetry %.4e at Vds = %g V, i.e. %.3g x V1's floor"
          % (a_dirac, VDS_RF, a_dirac / max(worst, 1e-300)),
          a_dirac > 1e-3, exact=False)

    # V3  EXACT, BITWISE.  rf.gate_capacitance's argument is Vg-Vch-V_dirac
    #     passed through an even function, so reflection about V_dirac is
    #     an exact identity on floats -- not an approximation to 1e-16.
    cp, cm = rf.gate_capacitance(gfet.V_dirac + d), \
             rf.gate_capacitance(gfet.V_dirac - d)
    bitwise = bool(np.all(cp == cm))
    check("V3 EXACT  C_gs is BITWISE even about V_dirac (wrong centre)",
          "max |diff| %.3e, bitwise-equal on all %d pairs: %s"
          % (float(np.max(np.abs(cp - cm))), d.size, bitwise), bitwise)

    # V4  EXACT.  The channel-averaged C_gs is even about V_MIRROR.
    ap = gate_capacitance_channel_averaged(v_mirror(VDS_RF) + d, Vds=VDS_RF)
    am = gate_capacitance_channel_averaged(v_mirror(VDS_RF) - d, Vds=VDS_RF)
    a4 = float(np.max(np.abs(ap - am) / (0.5 * (ap + am))))
    check("V4 EXACT  the channel-averaged C_gs IS even about V_dirac+Vds/2",
          "max rel asymmetry %.3e" % a4, a4 < 1e-13)

    # V5  EXACT RATIO 1.  f_T built on the channel-averaged C_gs has
    #     electron- and hole-branch peaks of EQUAL height, because |g_m| and
    #     C_avg now share a centre.  Nothing else in the expression changed.
    Vg = np.linspace(-2.6, 4.3, 691)
    fT_src, fT_avg, gm, Id, C_src, C_avg = fT_both_ways(Vg, Vds=VDS_RF)
    (ev_a, ef_a), (hv_a, hf_a) = branch_peaks(Vg, fT_avg, v_mirror(VDS_RF))
    (ev_s, ef_s), (hv_s, hf_s) = branch_peaks(Vg, fT_src, v_mirror(VDS_RF))
    r_avg, r_src = ef_a / hf_a, ef_s / hf_s
    check("V5 EXACT  repaired f_T: electron/hole peak ratio is 1",
          "ratio %.14f (|ratio-1| = %.3e)" % (r_avg, abs(r_avg - 1.0)),
          abs(r_avg - 1.0) < 1e-9)

    # V6  CONVERGENCE (weaker).  The source-end/channel-averaged C_gs
    #     discrepancy must vanish LINEARLY in Vds, since the average is a
    #     midpoint rule over an interval of length Vds.
    #     The gate voltage must be held FIXED while Vds halves: the first
    #     draft evaluated at v_mirror(Vds), which moves the evaluation point
    #     too, and measured ratios 1.407 / 2.428 / 3.392 -- a mixture of the
    #     interval shrinking and the point sliding along a strongly curved
    #     C_gs.
    #
    #     THE SECOND DRAFT'S PREMISE WAS ALSO WRONG, and is recorded here
    #     rather than quietly replaced.  It asserted O(Vds^2) "because the
    #     average is a midpoint rule", and predicted a ratio near 4; held at
    #     a fixed V_g it measured 3.494 / 2.394 / 2.159, i.e. heading to 2.
    #     The premise confused two different comparisons.  The average is
    #     compared here against the SOURCE-END value C(u), not against the
    #     exact channel integral, and
    #
    #         C_avg - C(u) = (1/Vds) int_0^Vds [C(u - V) - C(u)] dV + O(h^2)
    #                      = -C'(u) * Vds/2 + O(Vds^2)
    #
    #     so the discrepancy is FIRST order in Vds with a known coefficient.
    #     The quadrature error of the 50-point equal-weight rule is the
    #     O(h^2) term and is not what dominates.  The ratio must therefore
    #     tend to 2, and the leading coefficient is itself checkable (V6b).
    V_FIX = np.array([gfet.V_dirac + 0.25])
    c0 = float(rf.gate_capacitance(V_FIX)[0])
    rel, signed = [], []
    for Vds in (0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125):
        c1 = float(gate_capacitance_channel_averaged(V_FIX, Vds=Vds)[0])
        rel.append(abs(c1 - c0) / c0)
        signed.append((c1 - c0) / Vds)
    ratios = [rel[i] / rel[i + 1] for i in range(len(rel) - 1)]
    check("V6 RATE  C_gs discrepancy at FIXED V_g is O(Vds) (ratio -> 2)",
          "rel discrepancies %s; successive ratios %s"
          % (", ".join("%.3e" % r for r in rel),
             ", ".join("%.3f" % r for r in ratios)),
          abs(ratios[-1] - 2.0) < 0.25, exact=False)

    # V6b LIMIT, against an independently computed derivative.  The leading
    #     coefficient of (C_avg - C_src)/Vds must be -C'(u)/2, where C' is
    #     the slope of rf.gate_capacitance in V_g.  Richardson-extrapolate
    #     the two smallest drain biases and compare.
    hh = 1e-5
    dC = float((rf.gate_capacitance(V_FIX + hh)[0]
                - rf.gate_capacitance(V_FIX - hh)[0]) / (2 * hh))
    predicted = -0.5 * dC
    extrap = 2 * signed[-1] - signed[-2]
    err6b = abs(extrap - predicted) / abs(predicted)
    check("V6b LIMIT  the leading coefficient is -C'(V_g)/2",
          "Richardson limit %.6e F/V against -C'/2 = %.6e F/V, rel err "
          "%.2e; Richardson sequence %s"
          % (extrap, predicted, err6b,
             ", ".join("%.5e" % (2 * signed[i + 1] - signed[i])
                       for i in range(len(signed) - 1))),
          err6b < 1e-3, exact=False)

    # V7  THE ARTEFACT EXISTS IN THE FIGURE OF MERIT.  The committed f_T's
    #     two branch peaks are NOT equal, and the inequality is physical
    #     nonsense on a device whose I_d is exactly mirror-symmetric.
    check("V7 ARTEFACT  committed f_T: electron/hole peak ratio is NOT 1",
          "ratio %.9f (|ratio-1| = %.3e), peaks %.4f GHz at %+.3f V and "
          "%.4f GHz at %+.3f V" % (r_src, abs(r_src - 1.0), ef_s / 1e9,
                                   ev_s, hf_s / 1e9, hv_s),
          abs(r_src - 1.0) > 1e-7, exact=False)

    # V8  EXACT, BITWISE.  At Vds = 0 the repair is the committed function.
    ref = rf.gate_capacitance(Vg)
    z1 = gate_capacitance_channel_averaged(Vg, Vds=0.0, n_segments=1)
    check("V8a EXACT  the repair is a BITWISE no-op at Vds=0, n_segments=1",
          "bitwise-equal to rf.gate_capacitance on all %d points: %s"
          % (Vg.size, bool(np.all(z1 == ref))), bool(np.all(z1 == ref)))

    # V8b  NOT EXACT, and the reason is worth stating.  With the default
    #      50 segments the profile is 50 copies of 0.0 and the average is
    #      the float mean of 50 IDENTICAL values -- which is not the
    #      identity, because summing x fifty times rounds.  The first draft
    #      of this check asserted bitwise equality here and FAILED; the
    #      measured distance below is floating-point summation, not a model
    #      difference, and calling it "bitwise" would have been wrong.
    z50 = gate_capacitance_channel_averaged(Vg, Vds=0.0)
    ulp = np.abs(z50 - ref) / np.spacing(np.abs(ref))
    #      The bound asserted is the TEXTBOOK one for summing n floats,
    #      (n-1) ULP = 49 here; the first draft asserted 1 ULP from no
    #      argument at all and measured 8.  The gap between 8 and 49 is
    #      numpy's pairwise summation, which is why the measured number is
    #      printed next to the bound rather than the bound being tuned to it.
    check("V8b NEAR   with 50 segments at Vds=0 the mean of identical "
          "floats is within the (n-1) ULP summation bound",
          "max distance %.1f ULP against the n-1 = %d bound, exactly equal "
          "on %d of %d points"
          % (float(np.max(ulp)), N_SEGMENTS - 1,
             int(np.sum(z50 == ref)), Vg.size),
          float(np.max(ulp)) <= N_SEGMENTS - 1, exact=False)

    # V9  EXACT LAW.  The located current minimum is V_dirac + Vds/2, so the
    #     offset divided by Vds is 1/2 for every Vds.
    halves, mins = [], []
    for Vds in VDS_DECADES[1:]:
        halves.append((locate_mirror_centre(Vds) - gfet.V_dirac) / Vds)
        mins.append((locate_current_minimum(Vds) - gfet.V_dirac) / Vds)
    err9 = float(np.max(np.abs(np.array(halves) - 0.5)))
    err9m = float(np.max(np.abs(np.array(mins) - 0.5)))
    check("V9 EXACT  (reflection centre - V_dirac)/Vds = 1/2 for every Vds",
          "max |ratio - 0.5| = %.2e by root-find, against %.2e for the "
          "same quantity located as an extremum (the minimiser's floor, "
          "not the law's) -- over Vds = %s"
          % (err9, err9m, ", ".join("%g" % v for v in VDS_DECADES[1:])),
          err9 < 1e-11)

    return out, ok, dict(Vg=Vg, fT_src=fT_src, fT_avg=fT_avg, gm=gm, Id=Id,
                         C_src=C_src, C_avg=C_avg, r_src=r_src, r_avg=r_avg,
                         peaks_src=((ev_s, ef_s), (hv_s, hf_s)),
                         peaks_avg=((ev_a, ef_a), (hv_a, hf_a)),
                         a_dirac=a_dirac, a_mirror=worst)


# ---------------------------------------------------------------- census
CENSUS_SITES = (
    ("rf_small_signal_model.gate_capacitance",
     "line 224: dV = Vg_range - Vch - gfet.V_dirac, with Vch defaulting to "
     "0.0, so C_gs is the SOURCE-END capacitance while g_m belongs to a "
     "channel-averaged device at finite Vds",
     "MEASURED BELOW -- this is the figure-of-merit instance prediction 1 "
     "named"),
    ("graphene_diffusion_current_model.py:693  j_dirac",
     "j_dirac = int(np.argmin(np.abs(Vg - gfet.V_dirac))), used to report "
     "the diffusion share 'at the nearest grid point to the Dirac point' "
     "and as the left-hand side of check G4a",
     "MEASURED BELOW -- the grid point it selects is Vds/2 = 50 mV away "
     "from the bias at which this model's channel is actually neutral"),
    ("graphene_low_frequency_noise_model.py  V2",
     "takes a SEPARATE Vds -> 0 limit function (id_zero_vds_limit) to get "
     "an exact evenness check about V_dirac",
     "CORRECT AS WRITTEN, and strengthened rather than refuted: the exact "
     "check is available at finite Vds about V_dirac+Vds/2.  No number "
     "withdrawn."),
)


def measure_jdirac_site():
    """What the diffusion module's j_dirac selects, and what it would select
    at the model's actual neutrality bias.  Two gate points only -- the full
    sweep is that module's job, not this one's."""
    import graphene_diffusion_current_model as dd
    Vds = dd.VDS_REF
    Vg_grid = np.linspace(-2.0, 4.0, dd.N_VG)
    j_dirac = int(np.argmin(np.abs(Vg_grid - gfet.V_dirac)))
    j_mirror = int(np.argmin(np.abs(Vg_grid - v_mirror(Vds))))
    pts = np.array(sorted({Vg_grid[j_dirac], Vg_grid[j_mirror]}))
    I0 = dd.transfer_characteristic_dd(pts, Vds=Vds, lam=0.0)
    I1 = dd.transfer_characteristic_dd(pts, Vds=Vds, lam=1.0)
    share = (I1 / I0 - 1.0) * 100.0
    by_v = dict(zip(pts.tolist(), share.tolist()))
    return dict(Vds=Vds, n_vg=dd.N_VG,
                v_dirac_grid=float(Vg_grid[j_dirac]),
                v_mirror_grid=float(Vg_grid[j_mirror]),
                same_point=(j_dirac == j_mirror),
                share_at_dirac=by_v[float(Vg_grid[j_dirac])],
                share_at_mirror=by_v[float(Vg_grid[j_mirror])],
                grid_step=float(Vg_grid[1] - Vg_grid[0]))


# ---------------------------------------------------------------- figure
def make_figure(r, fname="ambipolar_mirror_audit.png"):
    Vg, Vm = r["Vg"], v_mirror(VDS_RF)
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6))

    d = np.linspace(0.02, 1.9, 220)
    for centre, lab, c in ((Vm, r"about $V_{dirac}+V_{ds}/2$", "tab:green"),
                           (gfet.V_dirac, r"about $V_{dirac}$", "tab:red")):
        hi, lo = _mirror_pairs(VDS_RF, centre, d)
        ax[0].semilogy(d, np.abs(hi - lo) / (0.5 * np.abs(hi + lo)),
                       color=c, lw=1.6, label=lab)
    ax[0].axhline(2.2e-16, color="k", ls=":", lw=1,
                  label="double precision")
    ax[0].set_xlabel(r"$|V_g - $ centre$|$  (V)")
    ax[0].set_ylabel(r"relative mirror asymmetry of $I_d$")
    ax[0].set_title(r"(a) $I_d$ at $V_{ds}=%g$ V is even about"
                    "\n" r"$V_{dirac}+V_{ds}/2$, not $V_{dirac}$" % VDS_RF)
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=.3)

    ax[1].plot(Vg, r["C_src"] * 1e18, lw=1.5, color="tab:red",
               label=r"committed $C_{gs}(V_{ch}=0)$")
    ax[1].plot(Vg, r["C_avg"] * 1e18, lw=1.5, color="tab:green", ls="--",
               label="channel-averaged")
    ax[1].axvline(gfet.V_dirac, color="tab:red", lw=.8, ls=":")
    ax[1].axvline(Vm, color="tab:green", lw=.8, ls=":")
    ax[1].set_xlabel(r"$V_g$  (V)")
    ax[1].set_ylabel(r"$C_{gs}$  (aF)")
    ax[1].set_title("(b) the two capacitances, and their\n"
                    "two reflection centres (dotted)")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=.3)

    dd2 = np.linspace(0.02, 1.9, 220)
    ax[2].semilogy(dd2, mirror_profile(dd2, lambda v:
                   np.abs(rf.transconductance(v, Vds=VDS_RF)[0])
                   / (2 * np.pi * rf.gate_capacitance(v)), Vm),
                   color="tab:red", lw=1.6, label="committed $f_T$")
    ax[2].semilogy(dd2, mirror_profile(dd2, lambda v:
                   np.abs(rf.transconductance(v, Vds=VDS_RF)[0])
                   / (2 * np.pi * gate_capacitance_channel_averaged(
                       v, Vds=VDS_RF)), Vm),
                   color="tab:green", lw=1.6, label="repaired $f_T$")
    ax[2].set_xlabel(r"$|V_g - (V_{dirac}+V_{ds}/2)|$  (V)")
    ax[2].set_ylabel(r"relative mirror asymmetry of $f_T$")
    ax[2].set_title("(c) the figure of merit inherits the\n"
                    "centre mismatch; the repair removes it")
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=.3)

    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------ main
def main():
    print("=" * 78)
    print("AMBIPOLAR MIRROR AUDIT -- where is the reflection centre, and")
    print("which expressions in this repository disagree about it?")
    print("=" * 78)
    print("  Scoring the 2026-11-15 prediction 1: 'the |g_m| audit will find")
    print("  at least one more instance, and it will be in a figure-of-merit")
    print("  expression rather than in a current model.'")
    print()

    print("-" * 78)
    print("SECTION 1.  THE LAW, measured over three decades of drain bias")
    print("-" * 78)
    d = np.linspace(0.02, 2.0, 41)
    print("  %8s %14s %16s %16s" % ("Vds [V]", "V_mirror [V]",
                                    "asym about V_m", "asym about V_dirac"))
    for Vds in VDS_DECADES:
        print("  %8.4g %14.6f %16.3e %16.3e"
              % (Vds, v_mirror(Vds), mirror_asymmetry(Vds, v_mirror(Vds), d),
                 mirror_asymmetry(Vds, gfet.V_dirac, d)))
    print()
    print("  and the LOCATED current minimum, by Brent (not grid argmin):")
    print("  %8s %18s %18s %12s" % ("Vds [V]", "V_g at min I_d",
                                    "V_dirac + Vds/2", "offset/Vds"))
    for Vds in VDS_DECADES[1:]:
        vmin = locate_current_minimum(Vds)
        print("  %8.4g %18.9f %18.9f %12.8f"
              % (Vds, vmin, v_mirror(Vds), (vmin - gfet.V_dirac) / Vds))
    print()
    print("  So the gate voltage this model's channel is actually neutral at")
    print("  is V_dirac + Vds/2 = %.3f V at the RF bias, not the %.3f V that")
    print("  `gfet.V_dirac` names.  The 50 mV is not an error in V_dirac --")
    print("  V_dirac is a FLAT-BAND parameter and is exactly what it says.")
    print("  It is the CURVE's centre that moves, and nothing in this")
    print("  repository had written the relation down." )

    print()
    print("-" * 78)
    print("SECTION 2.  CENSUS -- every site that treats a centre as known")
    print("-" * 78)
    for name, what, verdict in CENSUS_SITES:
        print("  SITE  %s" % name)
        print("        %s" % what)
        print("        -> %s" % verdict)
        print()

    print("-" * 78)
    print("SECTION 3.  THE FIGURE OF MERIT -- f_T's two centres")
    print("-" * 78)
    out, ok, r = validate()
    (ev_s, ef_s), (hv_s, hf_s) = r["peaks_src"]
    (ev_a, ef_a), (hv_a, hf_a) = r["peaks_avg"]
    print("  C_gs is BITWISE even about V_dirac; |g_m| is even about")
    print("  V_dirac + Vds/2.  f_T = |g_m|/(2 pi C_gs) is therefore even")
    print("  about neither, and the two branch peaks of a device whose I_d")
    print("  is exactly symmetric come out UNEQUAL:")
    print()
    print("  %-26s %14s %14s %14s" % ("", "electron peak", "hole peak",
                                      "ratio e/h"))
    print("  %-26s %9.5f GHz %9.5f GHz %14.10f"
          % ("committed C_gs(Vch=0)", ef_s / 1e9, hf_s / 1e9, r["r_src"]))
    print("  %-26s %9.5f GHz %9.5f GHz %14.10f"
          % ("channel-averaged C_gs", ef_a / 1e9, hf_a / 1e9, r["r_avg"]))
    print("  %-26s %9.3f V   %9.3f V" % ("peak bias, committed", ev_s, hv_s))
    print()
    print("  Magnitude, honestly: at the PEAK the two capacitances barely")
    print("  differ, because there C_gs is C_ox-dominated and flat -- the")
    print("  peak-f_T asymmetry is %.3e, i.e. the committed peak numbers in"
          % abs(r["r_src"] - 1.0))
    print("  Chapter 4 are NOT materially wrong and none is withdrawn.  The")
    print("  mismatch is large where C_q matters, i.e. near neutrality:")
    Vm = v_mirror(VDS_RF)
    for off in (0.05, 0.1, 0.2, 0.5, 1.0):
        cs = float(rf.gate_capacitance(np.array([Vm + off]))[0])
        cs2 = float(rf.gate_capacitance(np.array([Vm - off]))[0])
        print("    |Vg - V_m| = %4.2f V : committed C_gs differs between the "
              "two branches by %7.4f %%" % (off, 100 * abs(cs - cs2)
                                            / (0.5 * (cs + cs2))))
    print()
    print("  The claim is therefore about the SHAPE of f_T near the Dirac")
    print("  point, which is the region the repository's own RF chapter")
    print("  says it cannot be asked about for a different reason (the")
    print("  puddle floor, graphene_diffusion_current_model Section 3).")

    print()
    print("-" * 78)
    print("SECTION 4.  THE j_dirac SITE, measured")
    print("-" * 78)
    try:
        j = measure_jdirac_site()
        print("  graphene_diffusion_current_model.py Section 3 reports the")
        print("  diffusion share 'at the nearest grid point to the Dirac")
        print("  point'.  On its own %d-point grid (step %.6f V):"
              % (j["n_vg"], j["grid_step"]))
        print("    nearest grid point to V_dirac        = %+.6f V, "
              "share %+.6f %%" % (j["v_dirac_grid"], j["share_at_dirac"]))
        print("    nearest grid point to V_dirac+Vds/2  = %+.6f V, "
              "share %+.6f %%" % (j["v_mirror_grid"], j["share_at_mirror"]))
        print("    the two grid points are the same point: %s"
              % j["same_point"])
        print()
        if j["same_point"]:
            print("  VERDICT: the grid is coarse enough (%.3f V step) that"
                  % j["grid_step"])
            print("  both centres land on the SAME grid point, so check G4a's")
            print("  NUMBER is unaffected.  Its stated justification -- 'at")
            print("  charge neutrality' -- is nonetheless off by Vds/2 = %.3f"
                  % (0.5 * j["Vds"]))
            print("  V, and on a finer grid the number would move.  This is a")
            print("  LATENT defect, and it is annotated in place rather than")
            print("  silently repaired.")
        else:
            print("  VERDICT: the two grid points DIFFER, so the reported")
            print("  number is the share at the wrong bias.  Both are printed")
            print("  above; neither is deleted.")
    except Exception as exc:                       # pragma: no cover
        print("  NOT MEASURED: importing graphene_diffusion_current_model "
              "raised %s: %s" % (type(exc).__name__, exc))
        print("  Recorded as a failure rather than worked around.")

    print()
    print("-" * 78)
    print("SECTION 5.  VALIDATION")
    print("-" * 78)
    for name, detail, passed, exact in out:
        print("  [%s] %-62s" % ("PASS" if passed else "FAIL", name))
        print("         %s%s" % ("" if exact else "(not exact) ", detail))
    print()

    print("-" * 78)
    print("SECTION 6.  WHAT THIS DOES NOT ESTABLISH")
    print("-" * 78)
    print("  DOES NOT: the channel-averaged C_gs is a MIDPOINT average of")
    print("  the source-end formula, not a derived non-reciprocal")
    print("  capacitance.  Jimenez arXiv:1105.2734 obtains C_gs by")
    print("  differentiating channel-integrated charges with Ward-Dutton")
    print("  partitioning; that is a different and better object, and it is")
    print("  not built here.  What is claimed is only that the averaged form")
    print("  shares a reflection centre with g_m and the source-end one")
    print("  does not.")
    print("  DOES NOT: `rf.gate_capacitance` is left UNCHANGED, so every")
    print("  committed f_T/f_max number in this repository is reproducible")
    print("  bit for bit.  Nothing in Chapter 4 is withdrawn by this module.")
    print("  DOES NOT: the law V_m = V_dirac + Vds/2 is a property of THIS")
    print("  model's symmetric channel-potential grid and of mu_e = mu_h.")
    print("  graphene_branch_asymmetric_mobility.py breaks the second")
    print("  premise on purpose; the centre is not expected to survive it,")
    print("  and that case is not measured here.")
    print("  DOES NOT: velocity saturation and diffusion (vsat, dd modules)")
    print("  integrate over the same profile but with an asymmetric")
    print("  velocity field; whether the centre survives there is OPEN.")
    print()
    make_figure(r)
    print("  figure written: ambipolar_mirror_audit.png")
    print("=" * 78)
    print("VERDICT: %s over V1..V9 (11 checks)" % ("PASS" if ok else "FAIL"))
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
