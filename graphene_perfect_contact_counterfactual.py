"""
graphene_perfect_contact_counterfactual.py

Closes the TOP open item created 2026-10-05: *f_max/f_T at R_c = 0 EXACTLY,
in the saturated model.*  That item was written as a dichotomy, and it is
quoted here verbatim because the finding is that the dichotomy is not
exhaustive:

    "If a perfect contact also falls short of 1.3, four levers are
     exhausted and Section 7.8.1a's structural verdict is final; if it
     reaches the band, the RF case becomes a contact-engineering problem
     with a quantified target.  Either outcome is more useful than what
     this thesis currently carries."

RESULT, stated up front.  **BOTH halves are true at once, because the
item named a RATIO and the decision turns on its NUMERATOR.**  At the
literature-scale 40 um / 8-finger geometry and V_ds = 0.1 V, removing the
contact resistance exactly:

    f_max/f_T   0.683186  ->  0.759363   (+11.15 %; 58.41 % of 1.3,
                                            short by 0.5406)
    f_T         20.740 GHz ->  74.650 GHz  (x3.5994)
    f_max       14.169 GHz ->  56.687 GHz  (x4.0007)

So the *ratio* falls short -- four levers are exhausted and Section 7.8.1a's
structural verdict on f_max/f_T is final -- AND the RF case becomes a
contact-engineering problem with a quantified target anyway, because a
perfect contact is worth a factor of **4.0 in f_max itself**.  The item
could not see this because f_max/f_T divides out exactly the quantity the
contact dominates.

THE SECOND RESULT, which is about the counterfactual and not about the
device.  R_c enters this model in TWO places: as a series resistance in the
drain-current solve (`transfer_characteristic_saturated`'s bisection) and
as R_s = R_c/2 in the f_max denominator (`rf.source_access_resistance`).
Those are different code paths reading different names, and zeroing only
the first -- which is what passing `Rc_total=0.0` through `**kw` does, since
`fT_fmax_saturated` calls `rf.source_access_resistance()` which reads the
MODULE GLOBAL -- gives

    f_max/f_T   0.683186  ->  0.552349   (-19.15 %)

**The half-removed perfect contact answers the question with the opposite
sign.**  Both numbers are reported below, labelled, because the one that
answers the item is the coherent one and the plausible one is the
incoherent one: a 19 % degradation is exactly what a reader expects from
"the contacts are binding" read carelessly, and it is an artefact of
removing R_c from one of the two places it enters.

THE THIRD RESULT, a sign-level disagreement between the two models.  In the
superseded resistor model the SAME counterfactual, done coherently, gives

    f_max/f_T   0.642901  ->  0.581707   (-9.52 %)

i.e. **the resistor model is made WORSE by a perfect contact, and the
saturated model is made BETTER.**  Both columns compare the SAME
counterfactual -- the contact removed from the current path and from R_s
together -- so the two signs are not an artefact of different conventions.
This is not a discrepancy in magnitude,
it is a disagreement in sign, and it is stated rather than reconciled.  It
also has a diagnostic reading, which is the strongest single sentence this
module can offer: Section 7.8.1a's identity

    f_max/f_T  ~  (1/2) sqrt( R_total / (R_g + R_s) )

carries R_total in the NUMERATOR, so a model obeying it is rewarded for a
worse contact.  A figure of merit that improves when the device gets worse
is not measuring the device.  The resistor model's f_max was already known
(2026-10-04) to be a resistance ratio wearing f_max's name; today that
statement acquires a falsifiable consequence, and the consequence holds.

----------------------------------------------------------------------------
WHAT CONFIRMS THE 2026-10-05 MECHANISM, AND IT IS A PREDICTION, NOT A FIT
----------------------------------------------------------------------------
2026-10-05 measured a g_ds saturation factor of 1.1361 against a required
4.8483, and explained the shortfall mechanistically: mu*S/L = 0.3036 at that
bias "alone would divide the CHANNEL conductance by 1.699", and the measured
1.1361 is that 1.699 diluted by a contact resistance in series with the
channel.  That explanation makes a quantitative prediction it did not test:
**at R_c = 0 the dilution vanishes and the factor must be 1.699.**

This module measures it: **1.70661** at the item's own bias and **1.69621**
at the other reference bias, i.e. 100.45 % and 99.84 % of the predicted
1.699.  The 2026-10-05 mechanism is therefore confirmed to 0.45 % by a
measurement
taken from the other side of the dilution, which is a stronger form of
evidence than the agreement-between-two-artefacts the 2026-10-02 entry
warned about -- the 1.699 was computed from a channel integral and the
1.70661 from a finite-difference derivative of a self-consistent solve.

----------------------------------------------------------------------------
WHY "AT THE PEAK-f_T BIAS" IS NOT AVAILABLE HERE, AND WHAT IS USED INSTEAD
----------------------------------------------------------------------------
Every comparison in this module is at a FIXED gate bias, because with
R_c = 0 the peak of f_T(V_g) is **not interior to the sweep**: f_T rises
monotonically to the grid edge at V_g = +4 V in both the resistor and the
saturated model, so `argmax` returns the boundary and "peak f_T" would be a
property of where the sweep was stopped.  Check P1 below asserts this
explicitly -- interior at committed R_c, at the boundary at R_c = 0 -- so
the reason the fixed-bias convention is used is a measured fact in the
transcript rather than a choice in a docstring.  Widening the sweep is not
the repair: the peak moves with the sweep because the contact was what
produced it.

----------------------------------------------------------------------------
EXACTNESS
----------------------------------------------------------------------------
Per the 2026-10-03 rule that an exactness claim must name the operation it
is exact under:

  X1  exact under ABSENCE OF A SOLVE.  At R_c = 0 the bisection branch in
      `transfer_characteristic_saturated` is skipped entirely, so the result
      is bit-identical under any change of `tol` and `max_iter`.  Asserted
      at 0.0 ULP, and the committed-R_c case is asserted to MOVE under the
      same perturbation with its displacement printed (the 2026-10-03
      magnitude rule), so X1 is a two-sided claim about the code path.

  X2  exact under IDENTITY OF CODE PATH, against the COMMITTED TRANSCRIPT.
      Running this module's machinery at the committed R_c must reproduce
      velocity_saturation_output.txt's Section 4/5 numbers.  The four pinned
      values are listed in PINNED below and compared bitwise at the printed
      precision; this is the overlapping-limit validation required of any
      module that supersedes another (rule set 2026-10-03).

  X3  exact against a CLOSED FORM, i.e. against a known value and not
      against another artefact.  With n frozen and v_sat fixed, Eq. (4) at
      R_c = 0 has the analytic solution
          I_d / (W e n v_sat)  =  1 / (1 + L v_sat / (mu V_ds))
      which is checked to <= 2 ULP.  This is the X3b anchor of 2026-10-05
      re-taken on the R_c = 0 path specifically, so the branch this module
      lives on is anchored to arithmetic rather than to a sibling branch.

  X4  exact under MULTIPLICATION BY ZERO.  R_s = 0 and R_g = 0 must each
      delete their own addend from the f_max denominator exactly, with the
      reconstruction bitwise in the order the model sums them.

MAGNITUDES.  Every MUST_CHANGE check prints the size of the response it
measured, per the standing top methodological item (2026-10-03, restated
2026-10-04 and 2026-10-05).

Run:    python3 graphene_perfect_contact_counterfactual.py
Writes: perfect_contact_counterfactual_output.txt
        perfect_contact_counterfactual.png
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import graphene_fet_model as gfet
import rf_small_signal_model as rf
import graphene_velocity_saturation_model as vsat

VDS_REF = 0.1
VDS_LADDER = (0.05, 0.1, 0.2, 0.5, 1.0)
FEIJOO_LO, FEIJOO_HI = 1.3, 1.4

# The counterfactual's zero, named rather than written twice as a literal.
# X1a's whole claim is that this is EXACTLY zero and therefore skips the
# bisection branch in transfer_characteristic_saturated; a near-zero would
# restore the solve and make "exactly" false while changing no digit a
# tolerance-based check would notice.  Mutant M3 of the harness attacks this
# name for that reason.
RC_ZERO = 0.0

# Values pinned from the COMMITTED transcript velocity_saturation_output.txt
# (2026-10-05).  X2 compares against these, not against a re-run.
PINNED = {
    "gds_sat_at_resistor_peak": 3.354728e-02,
    "ratio_sat_at_resistor_peak": 6.851967e-01,
    "gds_sat_at_saturated_peak": 3.374491e-02,
    "ratio_sat_at_saturated_peak": 0.683186,
    # Added after the mutation harness's FIRST run, which found that M5 -- a
    # one-sided g_m stencil, i.e. an error of order dVg in f_T -- SURVIVED.
    # It survived because f_max/f_T is almost independent of f_T: term B is
    # 1.11 % of the f_max denominator here, so a 1e-5 relative shift in f_T
    # moves the ratio by ~5e-8 and no check that reads only the ratio can see
    # it.  That is a real property of this geometry and not a tolerance
    # mistake, so the repair is a pin on f_T ITSELF, which nothing had.
    "fT_sat_at_resistor_peak": 2.058532e+10,
    "fT_sat_at_saturated_peak": 2.0740e+10,
}
# The two f_T pins are REPORTED, not asserted equal, and the reason is a
# finding of its own.  The committed transcript's f_T comes from np.gradient
# over the 400-point gate grid (spacing 0.01504 V); this module differences a
# 3-point stencil at dVg = 1e-4.  Those are different ESTIMATORS of the same
# derivative, and the 2026-10-05 module's own rule -- "the comparison is
# between models and not between estimators" -- says a pin across them is not
# an exactness claim.  Measured deviation is ~1.4e-5 relative, i.e. the
# O(h^2) gate-grid discretisation of np.gradient, and it had never been
# quantified.  X2c reports it with a loose bound; X2b instead asserts the
# ORDER of this module's own stencil, which is what a one-sided g_m would
# break (mutant M5).
FT_ESTIMATOR_BOUND = 1e-4
# How many significant figures velocity_saturation_output.txt actually PRINTED
# for each pin.  Stated per pin rather than assumed uniform: that transcript
# prints its Section 4 table in %e with 6 decimals (7 figures), its Section 5
# ladder to 6 figures, and the saturated peak-f_T value in its Section 4
# header to 5.  Comparing all six at the tightest of those would fail on the
# formatting of the source rather than on the arithmetic -- the 2026-10-01
# fault, a reference that is not content-pinned, in the direction of a false
# alarm rather than a false pass.
PINNED_SIGFIGS = {
    "gds_sat_at_resistor_peak": 7,
    "ratio_sat_at_resistor_peak": 7,
    "gds_sat_at_saturated_peak": 7,
    "ratio_sat_at_saturated_peak": 6,
}
# Per-contact width-specific contact resistance, Ohm*um, for the
# reachability sweep of Section 6b.  gfet.Rc_per_width_ohm_um is 300.0; the
# others are literature values, each attached to its source in
# notes/2026-10-06-the-perfect-contact-and-the-ratio-that-divided-out-the-prize.md:
#   65   Liu et al. 2019, bottom-contact, e-beam -- the lowest reported
#   165  Feijoo et al., Nanoscale Adv. 2, 2020 (Rc*Wg/2), the device family
#        whose f_max/f_T = 1.3-1.4 band this thesis is measured against
#   470  Khosravi Rad et al., Sci. Rep. 14, 9190 (2024), "two-in-one" Ni
#        process -- a good result from a practical photolithographic flow
#   4000 the same paper's statement of what most of the literature exceeds
RC_PER_WIDTH_LADDER = (0.0, 65.0, 165.0, 300.0, 470.0, 4000.0)

# The 2026-10-05 channel-integral estimate of the UNDILUTED g_ds factor.
# Quoted from that transcript's Section 4 and used as a PREDICTION here.
PREDICTED_UNDILUTED_FACTOR = 1.699

_FAST = os.environ.get("PCC_FAST") == "1"
N_VG = 80 if _FAST else 400

_FAIL = []
_LINES = []


def say(s=""):
    _LINES.append(s)
    print(s)


def check(tag, ok, label, detail=""):
    say("  [%s] %-5s%s%s" % ("PASS" if ok else "FAIL", tag, label,
                             ("  -- " + detail) if detail else ""))
    if not ok:
        _FAIL.append(tag)
    return ok


def ulps(a, b):
    if a == b:
        return 0.0
    return abs(a - b) / np.spacing(max(abs(a), abs(b)))


# ---------------------------------------------------------------------------
# Geometry context.  Late binding only -- no module-level name is captured as
# a default (the 2026-10-03 frozen-default fault and its 2026-10-04 repeat).
# ---------------------------------------------------------------------------
def at_literature_geometry(fn, W=None, N_fingers=None, rc_per_width=None):
    """rc_per_width overrides the per-contact Ohm*um figure COHERENTLY: it
    sets gfet.Rc_total, which is read both by the drain-current solve and by
    rf.source_access_resistance(), so R_s moves with it.  Section 4's whole
    second result is that those two are separate paths and that moving one
    alone inverts the answer, so the sweep in Section 6b moves both."""
    W = rf.W_RF if W is None else W
    N_fingers = rf.N_FINGERS_RF if N_fingers is None else N_fingers
    rc_pw = (gfet.Rc_per_width_ohm_um if rc_per_width is None
             else rc_per_width)
    W_saved, Rc_saved = gfet.W, gfet.Rc_total
    gfet.W = W
    gfet.Rc_total = 2 * (rc_pw * 1e-6) / W
    try:
        return fn(N_fingers)
    finally:
        gfet.W, gfet.Rc_total = W_saved, Rc_saved


# ---------------------------------------------------------------------------
# The four configurations.  `rc0` zeroes the contact in the CURRENT path;
# `rs0` zeroes it in the f_max DENOMINATOR.  They are independent flags on
# purpose, because the whole second result of this module is that they are
# separate code paths and that setting only one of them inverts the answer.
# ---------------------------------------------------------------------------
def small_signal_at(V_g, N_fingers, Vds=VDS_REF, saturate=True,
                    rc0=False, rs0=False, dVg=1e-4, **kw):
    """f_T, f_max, g_m, g_ds, C_gs at ONE bias, at whatever geometry is set.

    g_m is a central difference on a 3-point stencil rather than
    np.gradient on a sweep, so a single-bias call cannot raise IndexError
    (the crash 2026-10-04 and 2026-10-05 both recorded in reporting paths).
    """
    if rc0:
        kw["Rc_total"] = RC_ZERO
    Vb = np.array([V_g - dVg, V_g, V_g + dVg])
    Id = vsat.transfer_characteristic_saturated(Vb, Vds=Vds,
                                                saturate=saturate, **kw)
    gm = (Id[2] - Id[0]) / (2.0 * dVg)
    gds = vsat.gds_at_bias(V_g, Vds=Vds, saturate=saturate, **kw)
    Cgs = float(rf.gate_capacitance(np.array([float(V_g)]))[0])
    Cgd = rf.Cgd_over_Cgs * Cgs
    Rg = rf.gate_resistance(N_fingers=N_fingers)
    Rs = 0.0 if rs0 else rf.source_access_resistance()
    fT = abs(gm) / (2 * np.pi * Cgs)
    termA = gds * (Rg + Rs)
    termB = 2 * np.pi * fT * Cgd * Rg
    denom = termA + termB
    fmax = fT / (2 * np.sqrt(denom)) if denom > 0 else np.inf
    return dict(V_g=float(V_g), gm=gm, gds=gds, Cgs=Cgs, Cgd=Cgd,
                Rg=Rg, Rs=Rs, fT=fT, fmax=fmax, ratio=fmax / fT,
                termA=termA, termB=termB, denom=denom,
                Id=float(Id[1]))


def sweep(Vg, N_fingers, Vds=VDS_REF, saturate=True, rc0=False, **kw):
    if rc0:
        kw["Rc_total"] = RC_ZERO
    gm, Id = vsat.gm_saturated(Vg, Vds=Vds, saturate=saturate, **kw)
    gds = vsat.gds_saturated(Vg, Vds=Vds, saturate=saturate, **kw)
    Cgs = rf.gate_capacitance(Vg)
    fT = np.abs(gm) / (2 * np.pi * Cgs)
    return fT, gds, Cgs


# ===========================================================================
def main():
    Vg_grid = np.linspace(-2.0, 4.0, N_VG)
    if _FAST:
        say("  [PCC_FAST=1: gate grid %d pts -- mutation-harness mode]" % N_VG)

    say("the perfect contact: f_max/f_T at R_c = 0 EXACTLY in the saturated model")
    say("(the 2026-10-05 top item; its dichotomy is quoted in the module docstring)")
    say()

    # The two committed reference biases, taken from the 400-point grid the
    # 2026-10-05 module used, so the fixed-bias comparisons land on the
    # numbers that transcript printed.
    Vg400 = np.linspace(-2.0, 4.0, 400)
    BIAS = (("resistor-peak", float(Vg400[59])),
            ("saturated-peak", float(Vg400[323])))

    say("=" * 78)
    say("SECTION 1.  Why every comparison below is at a FIXED bias")
    say("=" * 78)

    def sec1(N):
        out = {}
        for lbl, rc0 in (("committed R_c", False), ("R_c = 0", True)):
            for mlbl, sat in (("saturated", True), ("resistor", False)):
                fT, gds, Cgs = sweep(Vg_grid, N, saturate=sat, rc0=rc0)
                i = int(np.argmax(fT))
                interior = 0 < i < len(fT) - 1
                out[(lbl, mlbl)] = (Vg_grid[i], fT[i], interior)
        return out

    s1 = at_literature_geometry(sec1)
    say("  peak of f_T(V_g) over the committed sweep V_g in [-2, +4] V:")
    say("    %-15s %-10s %12s %14s %s"
        % ("R_c", "model", "V_g at peak", "f_T [Hz]", "interior?"))
    for k in sorted(s1):
        V, f, interior = s1[k]
        say("    %-15s %-10s %+12.6f %14.5e %s"
            % (k[0], k[1], V, f, "yes" if interior else "NO (grid edge)"))
    interior_cmt = all(s1[k][2] for k in s1 if k[0] == "committed R_c")
    edge_rc0 = all(not s1[k][2] for k in s1 if k[0] == "R_c = 0")
    check("P1", interior_cmt and edge_rc0,
          "the peak-f_T bias is INTERIOR at committed R_c and AT THE GRID EDGE "
          "at R_c = 0",
          "so 'at the peak-f_T bias' is not a well-defined comparison for the "
          "counterfactual; every table below is at a fixed bias")
    say("  The peak exists at committed R_c BECAUSE of the contact: I_d is")
    say("  capped by R_c at high |V_g - V_dirac| while C_gs keeps growing, so")
    say("  f_T = g_m/(2 pi C_gs) turns over.  Remove the cap and it does not.")
    say("  Widening the sweep moves the edge, not the physics.")
    say()

    say("=" * 78)
    say("SECTION 2.  X2: the overlapping limit against the COMMITTED transcript")
    say("=" * 78)

    def sec2(N):
        return {lbl: small_signal_at(V, N, saturate=True) for lbl, V in BIAS}

    s2 = at_literature_geometry(sec2)
    say("  quantity                                    committed       this module")
    pairs = (("gds_sat_at_resistor_peak", "resistor-peak", "gds",
              "g_ds, saturated, at the resistor peak [S]"),
             ("ratio_sat_at_resistor_peak", "resistor-peak", "ratio",
              "f_max/f_T, saturated, at the resistor peak"),
             ("gds_sat_at_saturated_peak", "saturated-peak", "gds",
              "g_ds, saturated, at the saturated peak [S]"),
             ("ratio_sat_at_saturated_peak", "saturated-peak", "ratio",
              "f_max/f_T, saturated, at the saturated peak"),
             ("fT_sat_at_resistor_peak", "resistor-peak", "fT",
              "f_T, saturated, at the resistor peak [Hz]"),
             ("fT_sat_at_saturated_peak", "saturated-peak", "fT",
              "f_T, saturated, at the saturated peak [Hz]"))
    x2_ok, x2b_ok, worst = True, True, 0.0
    for key, bias, field, label in pairs:
        got = s2[bias][field]
        want = PINNED[key]
        # Compare at the precision the committed transcript printed, which
        # for the two f_T pins is 5 significant figures and for the rest is 7.
        sig = PINNED_SIGFIGS.get(key)
        tol = (0.5 * 10.0 ** (np.floor(np.log10(abs(want))) - (sig - 1))
               if sig else np.inf)
        agree = abs(got - want) <= tol
        rel = abs(got - want) / abs(want)
        worst = max(worst, rel)
        if key.startswith("fT_sat"):
            x2b_ok = x2b_ok and (rel < FT_ESTIMATOR_BOUND)
            say("  %-42s %13.6e %13.6e   (rel %.3e, estimator difference)"
                % (label, want, got, rel))
            continue
        x2_ok = x2_ok and agree
        say("  %-42s %13.6e %13.6e%s"
            % (label, want, got, "" if agree else "   <-- DIFFERS"))
    check("X2", x2_ok,
          "exact under IDENTITY OF CODE PATH: the four committed g_ds and "
          "f_max/f_T values reproduced at printed precision",
          "pinned from velocity_saturation_output.txt, 2026-10-05, each at "
          "the figure count that transcript printed (see PINNED_SIGFIGS); "
          "worst relative deviation over all six pins %.3e" % worst)
    check("X2c", x2b_ok,
          "the committed f_T values are reproduced to within the ESTIMATOR "
          "difference, which is reported rather than asserted away",
          "bound %.0e; the committed numbers come from np.gradient on the "
          "400-point gate grid and this module differences a 3-point stencil "
          "at dVg = 1e-4, so the ~1.4e-5 gap is the O(h^2) gate-grid error of "
          "np.gradient and is quantified here for the first time"
          % FT_ESTIMATOR_BOUND)

    # X2b: the ORDER of this module's own g_m stencil.  Added after the
    # mutation harness's first run, in which M5 -- a one-sided g_m -- SURVIVED
    # every check, because f_max/f_T is almost independent of f_T here (term B
    # is 1.11 % of the denominator, so a 1e-5 shift in f_T moves the ratio by
    # ~5e-8) and because a pin against the committed f_T is a pin across two
    # estimators.  A Richardson ratio is neither: a central difference has
    # error O(h^2) and halving h divides it by 4; a one-sided difference has
    # error O(h) and halving h divides it by 2.  The check therefore reads the
    # ORDER of the stencil and not the value of anything.
    def sec2b(N):
        Vb = float(BIAS[1][1])
        g = [small_signal_at(Vb, N, saturate=True, dVg=h)["gm"]
             for h in (8e-3, 4e-3, 2e-3)]
        return g

    g8, g4, g2 = at_literature_geometry(sec2b)
    d1, d2 = g8 - g4, g4 - g2
    rich = d1 / d2 if d2 != 0 else np.inf
    say("  Richardson ratio of this module's own g_m stencil, at dVg = "
        "8e-3 / 4e-3 / 2e-3 V:")
    say("    g_m = %.10e, %.10e, %.10e" % (g8, g4, g2))
    say("    (g8-g4)/(g4-g2) = %.4f   [4 = central/O(h^2), 2 = one-sided/O(h)]"
        % rich)
    check("X2b", 3.2 <= rich <= 4.8,
          "ORDER: this module's g_m is a SECOND-order central difference, "
          "which is what a one-sided stencil would break",
          "Richardson ratio %.4f against 4 for O(h^2) and 2 for O(h); "
          "distance to the one-sided value %.4f" % (rich, abs(rich - 2.0)))
    say()

    say("=" * 78)
    say("SECTION 3.  X1/X3/X4: the R_c = 0 branch is exact, and anchored")
    say("=" * 78)

    def sec3(N):
        Vb = float(BIAS[1][1])
        out = {}
        # X1: no solve at R_c = 0, so tol/max_iter cannot matter.
        # RC_ZERO, NOT a literal 0.0.  As first written this check passed a
        # literal, and the mutation harness's first run showed M3 -- RC_ZERO
        # set to 1e-12 -- SURVIVING: the exactness claim was being made about
        # a value the counterfactual did not use.  That is the 2026-09-30
        # arrival fault (a mutation that does not arrive is indistinguishable
        # from a system that does not respond) inside an exactness check, and
        # it is the third instance this one harness run produced.
        a = vsat.transfer_characteristic_saturated(
            np.array([Vb]), Vds=VDS_REF, Rc_total=RC_ZERO, tol=1e-15,
            max_iter=200)[0]
        b = vsat.transfer_characteristic_saturated(
            np.array([Vb]), Vds=VDS_REF, Rc_total=RC_ZERO, tol=1e-3,
            max_iter=3)[0]
        out["x1_rc0"] = (a, b, ulps(a, b))
        c = vsat.transfer_characteristic_saturated(
            np.array([Vb]), Vds=VDS_REF, tol=1e-15, max_iter=200)[0]
        d = vsat.transfer_characteristic_saturated(
            np.array([Vb]), Vds=VDS_REF, tol=1e-3, max_iter=3)[0]
        out["x1_cmt"] = (c, d, ulps(c, d), abs(d - c) / abs(c))
        # X3: closed form at R_c = 0 with n frozen and v_sat fixed.
        v_const = 5.0e5
        Id, Q, S, n, vs_ = vsat._Id_given_Vds_ch(
            Vb, VDS_REF, constant_n=True, v_sat_const=v_const)
        n0 = float(n[0])
        lhs = Id / (gfet.W * vsat.E_CHARGE * n0 * v_const)
        rhs = 1.0 / (1.0 + gfet.L * v_const / (gfet.mu * VDS_REF))
        out["x3"] = (lhs, rhs, ulps(lhs, rhs))
        # X4: multiplication by zero in the denominator.
        full = small_signal_at(Vb, N, saturate=True)
        nors = small_signal_at(Vb, N, saturate=True, rs0=True)
        recon = full["gds"] * (full["Rg"] + 0.0) + full["termB"]
        out["x4"] = (nors["termA"] + nors["termB"], recon,
                     ulps(nors["denom"], recon))
        return out

    s3 = at_literature_geometry(sec3)
    a, b, u = s3["x1_rc0"]
    check("X1a", u == 0.0,
          "exact under ABSENCE OF A SOLVE: RC_ZERO (the constant the "
          "counterfactual actually uses) is bit-identical under tol "
          "1e-15 -> 1e-3, max_iter 200 -> 3",
          "%.1f ULP; I_d = %.15e both ways" % (u, a))
    c, d, u2, rel = s3["x1_cmt"]
    check("X1b", u2 > 0.0,
          "MUST_CHANGE: the committed-R_c case DOES move under the same "
          "perturbation, so X1a is a fact about the branch and not about the "
          "perturbation",
          "%.3e ULP, relative displacement %.3e (%.4f %%)"
          % (u2, rel, 100 * rel))
    lhs, rhs, u3 = s3["x3"]
    check("X3", u3 <= 2.0,
          "exact against a CLOSED FORM on the R_c = 0 path: "
          "I_d/(W e n v_sat) == 1/(1 + L v_sat/(mu V_ds))",
          "%.1f ULP; %.12f vs %.12f" % (u3, lhs, rhs))
    dn, recon, u4 = s3["x4"]
    check("X4", u4 == 0.0,
          "exact under MULTIPLICATION BY ZERO: R_s = 0 deletes its own addend "
          "bitwise, in the order the model sums them",
          "%.1f ULP; denom = %.15e" % (u4, dn))
    say()

    say("=" * 78)
    say("SECTION 4.  The counterfactual, and the sign that depends on coherence")
    say("=" * 78)

    def sec4(N):
        out = {}
        for blab, V in BIAS:
            for mlab, sat in (("resistor", False), ("saturated", True)):
                for clab, rc0, rs0 in (
                        ("baseline", False, False),
                        ("R_c=0 current path only", True, False),
                        ("R_c=0 COHERENT (R_s=0 too)", True, True),
                        ("R_s=0 only", False, True)):
                    out[(blab, mlab, clab)] = small_signal_at(
                        V, N, saturate=sat, rc0=rc0, rs0=rs0)
        return out

    s4 = at_literature_geometry(sec4)
    for blab, V in BIAS:
        say("  --- fixed bias: %s, V_g = %+0.6f V, V_ds = %.2f V"
            % (blab, V, VDS_REF))
        say("    %-10s %-27s %12s %12s %12s %10s"
            % ("model", "contact counterfactual", "g_ds [S]", "f_T [Hz]",
               "f_max [Hz]", "f_max/f_T"))
        for mlab in ("resistor", "saturated"):
            for clab in ("baseline", "R_c=0 current path only",
                         "R_c=0 COHERENT (R_s=0 too)", "R_s=0 only"):
                r = s4[(blab, mlab, clab)]
                say("    %-10s %-27s %12.6e %12.5e %12.5e %10.6f"
                    % (mlab, clab, r["gds"], r["fT"], r["fmax"], r["ratio"]))
        say()

    key = BIAS[1][0]           # the saturated-peak bias: the item's own bias
    base = s4[(key, "saturated", "baseline")]
    coh = s4[(key, "saturated", "R_c=0 COHERENT (R_s=0 too)")]
    half = s4[(key, "saturated", "R_c=0 current path only")]
    rbase = s4[(key, "resistor", "baseline")]
    rcoh = s4[(key, "resistor", "R_c=0 COHERENT (R_s=0 too)")]

    say("  THE ANSWER TO THE 2026-10-05 ITEM, at its own bias (%s):" % key)
    say("    f_max/f_T   %.6f -> %.6f   (%+.2f %%)"
        % (base["ratio"], coh["ratio"],
           100 * (coh["ratio"] / base["ratio"] - 1)))
    say("    f_T  [Hz]   %.5e -> %.5e   (x%.4f)"
        % (base["fT"], coh["fT"], coh["fT"] / base["fT"]))
    say("    f_max[Hz]   %.5e -> %.5e   (x%.4f)"
        % (base["fmax"], coh["fmax"], coh["fmax"] / base["fmax"]))
    say("    distance to the bottom of the Feijoo band (1.3): %.4f -> %.4f"
        % (FEIJOO_LO - base["ratio"], FEIJOO_LO - coh["ratio"]))
    say("    fraction of 1.3 reached              : %.2f %% -> %.2f %%"
        % (100 * base["ratio"] / FEIJOO_LO, 100 * coh["ratio"] / FEIJOO_LO))
    say()
    check("C1", coh["ratio"] < FEIJOO_LO,
          "MAGNITUDE: a PERFECT contact still falls short of the Feijoo band, "
          "so the FOURTH lever is exhausted and Section 7.8.1a's verdict on "
          "the RATIO is final",
          "ratio = %.6f, i.e. %.2f %% of 1.3; short by %.4f"
          % (coh["ratio"], 100 * coh["ratio"] / FEIJOO_LO,
             FEIJOO_LO - coh["ratio"]))
    check("C2", coh["fmax"] / base["fmax"] > 2.0,
          "MAGNITUDE: and the SAME counterfactual is worth more than x2 in "
          "f_max itself, so the item's dichotomy was not exhaustive",
          "f_max x%.4f, f_T x%.4f, ratio only %+.2f %% -- the prize is in the "
          "numerator the ratio divides out"
          % (coh["fmax"] / base["fmax"], coh["fT"] / base["fT"],
             100 * (coh["ratio"] / base["ratio"] - 1)))
    check("C3", (half["ratio"] - base["ratio"]) *
          (coh["ratio"] - base["ratio"]) < 0.0,
          "MUST_CHANGE, SIGN: the half-removed contact (current path only, "
          "R_s left in place) answers the item with the OPPOSITE sign",
          "coherent %+.2f %%, half-removed %+.2f %% -- removing R_c from one "
          "of the two places it enters inverts the verdict"
          % (100 * (coh["ratio"] / base["ratio"] - 1),
             100 * (half["ratio"] / base["ratio"] - 1)))
    check("C4", (rcoh["ratio"] - rbase["ratio"]) *
          (coh["ratio"] - base["ratio"]) < 0.0,
          "MUST_CHANGE, SIGN: the two MODELS disagree on the sign of the "
          "contact lever -- the resistor model is made WORSE by a perfect "
          "contact",
          "resistor %.6f -> %.6f (%+.2f %%), saturated %.6f -> %.6f "
          "(%+.2f %%)"
          % (rbase["ratio"], rcoh["ratio"],
             100 * (rcoh["ratio"] / rbase["ratio"] - 1),
             base["ratio"], coh["ratio"],
             100 * (coh["ratio"] / base["ratio"] - 1)))
    say("  WHY C4 IS A DIAGNOSTIC AND NOT A DISCREPANCY.  Section 7.8.1a's")
    say("  identity f_max/f_T ~ (1/2) sqrt(R_total/(R_g+R_s)) carries R_total")
    say("  in the NUMERATOR, so any model obeying it is REWARDED for a worse")
    say("  contact.  The identity was established 2026-10-04 as a description;")
    say("  C4 is its falsifiable consequence, and it holds.  A figure of merit")
    say("  that improves when the device gets worse is not measuring the")
    say("  device -- which is what 2026-10-04 claimed about this quantity, now")
    say("  with a sign test behind it rather than an algebraic reading.")
    say()

    say("=" * 78)
    say("SECTION 5.  The 2026-10-05 dilution mechanism, tested as a prediction")
    say("=" * 78)

    def sec5(N):
        out = {}
        for blab, V in BIAS:
            for clab, rc0 in (("committed R_c", False), ("R_c = 0", True)):
                r = small_signal_at(V, N, saturate=False, rc0=rc0)
                s = small_signal_at(V, N, saturate=True, rc0=rc0)
                out[(blab, clab)] = (r["gds"], s["gds"], r["gds"] / s["gds"])
        return out

    s5 = at_literature_geometry(sec5)
    say("  g_ds saturation factor = g_ds(no saturation) / g_ds(saturated)")
    say("    %-16s %-15s %14s %14s %9s"
        % ("bias", "R_c", "no-sat [S]", "saturated [S]", "factor"))
    for k in sorted(s5):
        a_, b_, f_ = s5[k]
        say("    %-16s %-15s %14.6e %14.6e %9.5f" % (k[0], k[1], a_, b_, f_))
    f_cmt = s5[(key, "committed R_c")][2]
    f_rc0 = s5[(key, "R_c = 0")][2]
    pred = PREDICTED_UNDILUTED_FACTOR
    say()
    say("  2026-10-05 PREDICTED, from mu*S/L = 0.3036 and a channel integral:")
    say("    the UNDILUTED factor is %.3f, and 1.1361 is that number diluted"
        % pred)
    say("    by a contact resistance in series with the channel.")
    say("  MEASURED here at R_c = 0, from a finite difference of a")
    say("    self-consistent solve:                                 %.5f"
        % f_rc0)
    say("    agreement with the prediction:                         %.2f %%"
        % (100 * f_rc0 / pred))
    say("    dilution actually attributable to the contact:         x%.5f"
        % (f_rc0 / f_cmt))
    check("D1", abs(f_rc0 / pred - 1.0) < 0.02,
          "MAGNITUDE: the R_c = 0 saturation factor reproduces the 2026-10-05 "
          "channel-integral prediction of 1.699 to better than 2 %",
          "measured %.5f vs predicted %.3f, ratio %.5f (%.2f %% error) -- the "
          "mechanism was stated 2026-10-05 and is confirmed today from the "
          "other side of the dilution"
          % (f_rc0, pred, f_rc0 / pred, 100 * abs(f_rc0 / pred - 1.0)))
    check("D2", f_rc0 > f_cmt,
          "MUST_CHANGE: removing the contact INCREASES the saturation factor, "
          "which is the direction the dilution explanation requires",
          "%.5f -> %.5f, i.e. the contacts were absorbing a factor of %.5f "
          "of the available saturation" % (f_cmt, f_rc0, f_rc0 / f_cmt))
    say()

    say("=" * 78)
    say("SECTION 6.  What is left to ask of g_ds after a perfect contact")
    say("=" * 78)

    def sec6(N):
        Vb = float(BIAS[1][1])
        c = small_signal_at(Vb, N, saturate=True, rc0=True, rs0=True)
        need_denom = 1.0 / (4.0 * FEIJOO_LO ** 2)
        # Only term A scales with g_ds; term B scales with f_T, which is held.
        have = c["denom"]
        residual_gds_factor = (c["termA"] / (need_denom - c["termB"])
                               if need_denom > c["termB"] else np.inf)
        # INDEPENDENT recomputation of the same factor, by bisection on the
        # full f_max expression rather than by rearranging it.  Added after
        # the mutation harness's first run, in which M6 -- dropping term B
        # from the closed form above -- SURVIVED, because R1 asserted only
        # that the residual was finite and greater than one and nothing read
        # its value.  The two routes share no algebra: this one never forms
        # need_denom at all, it only evaluates f_max/f_T and compares to 1.3.
        def ratio_at(F):
            return 1.0 / (2.0 * np.sqrt(c["termA"] / F + c["termB"]))
        lo, hi = 1.0, 1e9
        for _ in range(300):
            mid = 0.5 * (lo + hi)
            if ratio_at(mid) < FEIJOO_LO:
                lo = mid
            else:
                hi = mid
        root = 0.5 * (lo + hi)
        return c, need_denom, have, residual_gds_factor, root

    c, need_denom, have, res_factor, res_root = at_literature_geometry(sec6)
    say("  at R_c = 0 COHERENT, at the saturated-peak bias:")
    say("    term A = g_ds (R_g + R_s) = %.6e  (%.2f %% of the denominator)"
        % (c["termA"], 100 * c["termA"] / c["denom"]))
    say("    term B = 2 pi f_T C_gd R_g = %.6e  (%.2f %%)"
        % (c["termB"], 100 * c["termB"] / c["denom"]))
    say("    denominator now                  = %.6e" % have)
    say("    denominator needed for f_max/f_T = 1.3 : %.6e" % need_denom)
    say("    so g_ds must come down a FURTHER factor of  %.4f" % res_factor)
    say("    (and 2026-10-05 showed no phonon energy in 0.059-0.196 eV buys")
    say("     more than 1.2491, so this residual is not available from v_sat)")
    say("    same factor by BISECTION on f_max/f_T, sharing no algebra "
        "with the")
    say("    closed form above                        :  %.4f" % res_root)
    check("R1b", abs(res_root / res_factor - 1.0) < 1e-9,
          "ROUND TRIP: the residual factor from the closed form and from an "
          "independent bisection on f_max/f_T agree",
          "closed form %.10f, bisection %.10f, relative difference %.3e -- "
          "a magnitude pin, which R1 alone was not"
          % (res_factor, res_root, abs(res_root / res_factor - 1.0)))
    check("R1", np.isfinite(res_factor) and res_factor > 1.0,
          "MAGNITUDE: a residual g_ds requirement survives the perfect "
          "contact, and it is finite -- so the shortfall is still a g_ds "
          "shortfall and not a parasitic one",
          "residual factor %.4f after R_c = 0; term B is only %.2f %% of the "
          "denominator, so term A is still what has to move"
          % (res_factor, 100 * c["termB"] / c["denom"]))
    say()

    say("=" * 78)
    say("SECTION 6b.  How much of the perfect contact is REACHABLE")
    say("=" * 78)
    say("  Section 4 prices a perfect contact at x%.4f in f_max.  A perfect"
        % (coh["fmax"] / base["fmax"]))
    say("  contact does not exist, so the sweep below asks the same question")
    say("  at contact resistances that have been MEASURED.  Each row sets the")
    say("  per-contact Ohm*um figure coherently -- the drain-current solve and")
    say("  R_s together -- and the sources are in RC_PER_WIDTH_LADDER's")
    say("  comment.  The fraction captured is measured against the R_c = 0")
    say("  row of this same sweep, not against Section 4's counterfactual, so")
    say("  the comparison is within one code path.")
    say()

    Vb_r = float(BIAS[1][1])
    ladder = []
    for rc_pw in RC_PER_WIDTH_LADDER:
        r = at_literature_geometry(
            lambda n: small_signal_at(Vb_r, n, saturate=True),
            rc_per_width=rc_pw)
        ladder.append((rc_pw, r))
    fmax0 = ladder[0][1]["fmax"]
    fmax_300 = [r for pw, r in ladder if pw == 300.0][0]["fmax"]
    say("    %10s %12s %13s %12s %11s %13s"
        % ("Rc [O*um]", "R_c,tot [O]", "g_ds [S]", "f_T [Hz]", "f_max/f_T",
           "f_max [GHz]"))
    for rc_pw, r in ladder:
        say("    %10.0f %12.4f %13.6e %12.5e %11.6f %13.4f"
            % (rc_pw, 2 * (rc_pw * 1e-6) / rf.W_RF, r["gds"], r["fT"],
               r["ratio"], r["fmax"] / 1e9))
    say()
    say("    %10s %14s %22s" % ("Rc [O*um]", "f_max/f_max(300)",
                                "fraction of the R_c=0 gain captured"))
    for rc_pw, r in ladder:
        gain = r["fmax"] / fmax_300
        frac = ((r["fmax"] - fmax_300) / (fmax0 - fmax_300)
                if fmax0 != fmax_300 else float("nan"))
        say("    %10.0f %14.4f %21.1f %%" % (rc_pw, gain, 100 * frac))
    r65 = [r for pw, r in ladder if pw == 65.0][0]
    r165 = [r for pw, r in ladder if pw == 165.0][0]
    check("Q1", r65["fmax"] > fmax_300 and r65["ratio"] > base["ratio"],
          "MAGNITUDE: the LOWEST reported graphene contact resistance "
          "(65 Ohm*um) captures a measurable share of the perfect contact",
          "f_max x%.4f against this model's 300 Ohm*um, i.e. %.1f %% of the "
          "R_c = 0 gain; f_max/f_T %.6f -> %.6f, still %.2f %% of 1.3"
          % (r65["fmax"] / fmax_300,
             100 * (r65["fmax"] - fmax_300) / (fmax0 - fmax_300),
             base["ratio"], r65["ratio"], 100 * r65["ratio"] / FEIJOO_LO))
    check("Q2", r165["ratio"] < FEIJOO_LO,
          "and at the contact resistance of the device family whose "
          "1.3-1.4 band this thesis is measured against (165 Ohm*um), the "
          "model is STILL below that band",
          "f_max/f_T = %.6f, i.e. %.2f %% of 1.3 -- so the shortfall is not "
          "explained by this model having a worse contact than Feijoo et al."
          % (r165["ratio"], 100 * r165["ratio"] / FEIJOO_LO))
    ratios = [r["ratio"] for _, r in ladder]
    imin = int(np.argmin(ratios))
    interior_min = 0 < imin < len(ratios) - 1
    check("Q3", interior_min,
          "the ratio f_max/f_T is NON-MONOTONIC in R_c, with an INTERIOR "
          "minimum -- so beyond a point the saturated model is rewarded for "
          "a worse contact too",
          "minimum %.6f at %.0f Ohm*um, rising to %.6f at %.0f Ohm*um; this "
          "is C4's resistance-ratio reward reasserting itself in the "
          "SATURATED model once R_c dominates the channel again"
          % (ratios[imin], RC_PER_WIDTH_LADDER[imin], ratios[-1],
             RC_PER_WIDTH_LADDER[-1]))
    say("  Q3 is a NEW result and it sharpens C4 rather than contradicting")
    say("  it.  C4 says the resistor model is rewarded for a bad contact")
    say("  because its f_max/f_T is a resistance ratio.  Q3 says the")
    say("  SATURATED model is too, beyond %.0f Ohm*um: as R_c grows it"
        % RC_PER_WIDTH_LADDER[imin])
    say("  eventually dominates R_total faster than it dominates R_g+R_s, and")
    say("  the ratio turns back up while f_max itself falls by %.0fx across"
        % (ladder[imin][1]["fmax"] / ladder[-1][1]["fmax"]))
    say("  the same two rows.  The ratio and the figure of merit it is")
    say("  supposed to summarise move in OPPOSITE directions over part of")
    say("  the design space, which is the clearest statement available of")
    say("  why 2026-10-05's item could be satisfied on the ratio and still")
    say("  miss the finding.")
    say()
    say("  READ Q2 CAREFULLY, because it is the load-bearing row.  This")
    say("  model's 300 Ohm*um is 1.8x the 165 Ohm*um Feijoo et al. report for")
    say("  the devices whose f_max/f_T band Chapter 4 is compared against.  A")
    say("  natural objection to the whole f_max thread is therefore that the")
    say("  shortfall is just a worse contact.  Q2 answers it: give this model")
    say("  Feijoo's own contact and it reaches %.6f, not 1.3." % r165["ratio"])
    say()

    say("=" * 78)
    say("SECTION 7.  The V_ds ladder, with and without the contact")
    say("=" * 78)

    def sec7(N):
        Vb = float(BIAS[1][1])
        rows = []
        for Vds in VDS_LADDER:
            b = small_signal_at(Vb, N, Vds=Vds, saturate=True)
            cc = small_signal_at(Vb, N, Vds=Vds, saturate=True,
                                 rc0=True, rs0=True)
            rows.append((Vds, b["ratio"], cc["ratio"],
                         b["fmax"], cc["fmax"]))
        return rows

    rows = at_literature_geometry(sec7)
    say("    %6s %14s %16s %13s %15s"
        % ("V_ds", "ratio baseline", "ratio R_c=0 coh.", "f_max base",
           "f_max R_c=0"))
    for Vds, rb, rc, fb, fc in rows:
        say("    %6.2f %14.6f %16.6f %13.4e %15.4e" % (Vds, rb, rc, fb, fc))
    d_base = rows[-1][1] - rows[0][1]
    d_coh = rows[-1][2] - rows[0][2]
    check("L1", d_coh > 0 and d_base > 0,
          "the rising f_max/f_T trend that velocity saturation bought "
          "(2026-10-05 criterion B) SURVIVES the perfect contact",
          "baseline %+.6f over V_ds = 0.05-1 V, R_c = 0 coherent %+.6f"
          % (d_base, d_coh))
    crossings = [Vds for Vds, rb, rc, _, _ in rows
                 if rc >= FEIJOO_LO and rb < FEIJOO_LO]
    say("    V_ds rows where the perfect contact enters the band and the")
    say("    baseline does not: %s"
        % (", ".join("%.2f V" % v for v in crossings) if crossings else "none"))
    say("    NOT offered as agreement with Feijoo et al.: 2026-10-04 and")
    say("    2026-10-05 both recorded that a x20 bias knob against a x20")
    say("    literature band will coincide somewhere.  The ladder is here for")
    say("    the SIGN and for whether the contact changes it.")
    say()

    say("=" * 78)
    say("SECTION 8.  Controls")
    say("=" * 78)

    def sec8(N):
        Vb = float(BIAS[1][1])
        a = small_signal_at(Vb, N, saturate=True)
        b = small_signal_at(Vb, N, saturate=True, rc0=True, rs0=True)
        # Control: R_c = 0 is not the same operation as a large W.  Both make
        # Rc_total smaller, but W also rescales g_m, g_ds and C_gs, so if the
        # two were aliased the counterfactual would be a geometry change.
        W4 = at_literature_geometry(
            lambda n: small_signal_at(Vb, n, saturate=True), W=4 * rf.W_RF,
            N_fingers=rf.N_FINGERS_RF)
        return a, b, W4

    s8a, s8b, s8W = at_literature_geometry(sec8)
    r_rc0 = s8b["ratio"] / s8a["ratio"]
    r_W = s8W["ratio"] / s8a["ratio"]
    check("K1", abs(r_rc0 - r_W) / max(r_rc0, r_W) > 0.05,
          "CONTROL: R_c = 0 and a 4x wider device are NOT aliased responses, "
          "so the counterfactual is a contact change and not a geometry change",
          "R_c=0 response x%.5f, 4x-W response x%.5f, ratio %.5f"
          % (r_rc0, r_W, r_rc0 / r_W))
    check("K2",
          gfet.W == 1e-6 and gfet.Rc_total == 600.0 and gfet.mu == 0.4
          and gfet.L == 200e-9,
          "CONTROL: gfet's module globals are restored to their committed "
          "values after every geometry context in this module",
          "W = %g, Rc_total = %g, mu = %g, L = %g"
          % (gfet.W, gfet.Rc_total, gfet.mu, gfet.L))
    rerun = at_literature_geometry(
        lambda n: small_signal_at(float(BIAS[1][1]), n, saturate=True))
    check("K3", ulps(rerun["ratio"], s8a["ratio"]) == 0.0,
          "CONTROL: the baseline is bitwise reproducible after every "
          "counterfactual in this module has run",
          "%.1f ULP" % ulps(rerun["ratio"], s8a["ratio"]))
    say()

    say("=" * 78)
    say("SECTION 9.  Figure")
    say("=" * 78)

    def sec9(N):
        fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))

        # (a) f_T sweeps, showing the peak leaving the grid at R_c = 0.
        for lbl, rc0, sat, st in (("saturated, committed R_c", False, True, "-"),
                                  ("saturated, R_c = 0", True, True, "--"),
                                  ("resistor, committed R_c", False, False, "-"),
                                  ("resistor, R_c = 0", True, False, "--")):
            fT, gds, Cgs = sweep(Vg_grid, N, saturate=sat, rc0=rc0)
            ax[0].semilogy(Vg_grid, fT / 1e9, st, lw=1.6, label=lbl)
        for _, V in BIAS:
            ax[0].axvline(V, color="0.6", lw=0.8, ls=":")
        ax[0].set_xlabel("$V_g$ (V)")
        ax[0].set_ylabel("$f_T$ (GHz)")
        ax[0].set_title("(a) the peak leaves the sweep when $R_c=0$")
        ax[0].legend(fontsize=7.5)
        ax[0].grid(alpha=0.3, which="both")

        # (b) the four counterfactuals at the saturated-peak bias.
        labels = ["baseline", "$R_c=0$\ncurrent only", "$R_c=0$\ncoherent",
                  "$R_s=0$\nonly"]
        clabs = ("baseline", "R_c=0 current path only",
                 "R_c=0 COHERENT (R_s=0 too)", "R_s=0 only")
        x = np.arange(4)
        for off, mlab, col in ((-0.19, "resistor", "tab:orange"),
                               (0.19, "saturated", "tab:blue")):
            vals = [s4[(key, mlab, c_)]["ratio"] for c_ in clabs]
            ax[1].bar(x + off, vals, 0.36, label=mlab, color=col)
        ax[1].axhspan(FEIJOO_LO, FEIJOO_HI, color="green", alpha=0.15,
                      label="Feijoo et al. 1.3-1.4")
        ax[1].set_xticks(x)
        ax[1].set_xticklabels(labels, fontsize=8)
        ax[1].set_ylabel("$f_{max}/f_T$")
        ax[1].set_title("(b) the two models disagree on the SIGN")
        ax[1].legend(fontsize=7.5)
        ax[1].grid(alpha=0.3, axis="y")

        # (c) f_max itself: where the prize is.
        vals_b = [r[3] / 1e9 for r in rows]
        vals_c = [r[4] / 1e9 for r in rows]
        ax[2].loglog([r[0] for r in rows], vals_b, "o-", label="baseline")
        ax[2].loglog([r[0] for r in rows], vals_c, "s--",
                     label="$R_c=0$ coherent")
        ax[2].set_xlabel("$V_{ds}$ (V)")
        ax[2].set_ylabel("$f_{max}$ (GHz)")
        ax[2].set_title("(c) $f_{max}$ gains $\\times%.2f$ where the ratio "
                        "gains %.0f %%"
                        % (coh["fmax"] / base["fmax"],
                           100 * (coh["ratio"] / base["ratio"] - 1)))
        ax[2].legend(fontsize=8)
        ax[2].grid(alpha=0.3, which="both")

        fig.suptitle("The perfect contact: $R_c=0$ exactly, saturated model, "
                     "40 um / 8 fingers, $V_{ds}=0.1$ V unless shown",
                     fontsize=10)
        fig.tight_layout()
        fig.savefig("perfect_contact_counterfactual.png", dpi=150)
        plt.close(fig)

    at_literature_geometry(sec9)
    say("  wrote perfect_contact_counterfactual.png")
    say()

    say("=" * 78)
    npass = sum(1 for l in _LINES if "[PASS]" in l)
    nfail = sum(1 for l in _LINES if "[FAIL]" in l)
    say("RESULT: %d lines, %d passed, %d failed" % (len(_LINES), npass, nfail))
    say()
    say("  The 2026-10-05 item, answered     : f_max/f_T = %.6f at R_c = 0 "
        "exactly" % coh["ratio"])
    say("  Verdict on the RATIO              : SHORT of 1.3 by %.4f "
        "(%.2f %% of the band floor) -- FOURTH lever exhausted"
        % (FEIJOO_LO - coh["ratio"], 100 * coh["ratio"] / FEIJOO_LO))
    say("  Verdict on f_max ITSELF           : x%.4f -- a contact-engineering "
        "target the ratio divided out" % (coh["fmax"] / base["fmax"]))
    say("  2026-10-05's dilution mechanism   : CONFIRMED to %.2f %% "
        "(%.5f measured vs 1.699 predicted)"
        % (100 * abs(f_rc0 / pred - 1.0), f_rc0))
    say("  The coherence trap                : half-removal inverts the sign "
        "(%+.2f %% vs %+.2f %%)"
        % (100 * (half["ratio"] / base["ratio"] - 1),
           100 * (coh["ratio"] / base["ratio"] - 1)))
    say("=" * 78)

    with open("perfect_contact_counterfactual_output.txt", "w") as fh:
        fh.write("\n".join(_LINES) + "\n")

    return 1 if _FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
