"""
graphene_screening_overdrive_invariance.py

THE GATE STACK DOES NOT DECIDE WHETHER BRANCH ASYMMETRY IS SCREENED, AND THE
56x SPREAD THAT SAID IT DID CAME FROM HOLDING n_c FIXED ACROSS STACKS THAT DO
NOT SHARE AN n_c.

Created 2026-11-07 to close the open item logged at the end of
notes/2026-10-30-density-dependent-mobility-and-the-gate-stack.md:

    "This is a new open item, not something the committed module handles:
     n_c is ... a charged-impurity density that the gate stack determines."

WHAT 10-30 CONCLUDED, AND THE ONE QUANTITY IT HELD FIXED
--------------------------------------------------------
Section 3 of graphene_density_dependent_mobility.py turns the crossover
density into a crossover gate overdrive,

    dV_c = e * n_c / C_ox                                                (A1)

and evaluates it for two stacks with ONE value of n_c, N_CROSSOVER = 3e16
m^-2, for both rows:

    90 nm SiO2 back gate   C_ox = 3.84e-4 F/m^2   dV_c ~ 12.5 V
    10 nm HfO2 top gate    C_ox = 2.21e-2 F/m^2   dV_c ~ 0.22 V

and concludes that whether the electron/hole asymmetry is a first-order effect
"is a property of the GATE STACK, not of the graphene", being screened away in
the thin-oxide stacks a semiconductor process would actually build.

THE CONTRADICTION THIS MODULE REPORTS
-------------------------------------
n_c is not a constant across those two rows.  Two committed results say so:

  * graphene_impurity_limited_crossover.py (2026-11-02): the invariant under
    charged-impurity density is the PRODUCT A = n_c * mu_L_bar, so
    n_c = KAPPA * N_imp -- n_c is proportional to the charged-impurity
    density.
  * graphene_kappa_impurity_bound.py (2026-11-03): the charged-impurity
    density is bounded BELOW by the net charge inferred from the
    neutrality-point offset, N_imp >= N_A^(S) = C_ox * V_NP / e, which is how
    Zion et al. estimate it for their own sample.

Substituting N_imp = f * N_A^(S) with f >= 1 into (A1), the C_ox in n_c
cancels the C_ox in the denominator, EXACTLY and not approximately:

    dV_c = e * KAPPA * f * (C_ox * V_NP / e) / C_ox  =  KAPPA * f * V_NP   (A2)

    x == |n| / n_c = dV / dV_c = dV / (KAPPA * f * V_NP)                   (A3)

So the crossover overdrive is independent of oxide capacitance -- of oxide
thickness and of oxide permittivity alike -- and the 56.3x spread between the
two rows of 10-30's table is an artefact of using the back gate's n_c for the
top-gated device.  Under (A2) the two stacks need the SAME overdrive to reach
n_c, and the difference is exactly zero rather than small.

This does not delete 10-30's claim; it names its condition.  "A property of
the gate stack" survives, but the stack enters through V_NP -- how much net
charge its deposition leaves behind -- and NOT through C_ox.  A thinner or
higher-permittivity dielectric buys no screening at all.  The 11-02
annotation already flagged that "0.22 V uses the REFERENCE n_c"; this module
is the quantitative version of that flag and it is larger than a correction to
one cell, because the cell that moves is the one the conclusion rested on.

THE DIRECTION OF THE INEQUALITY, WHICH IS WHAT MAKES (A3) USABLE
----------------------------------------------------------------
N_A^(S) is a NET charge and therefore a LOWER bound on the impurity count
(graphene_kappa_impurity_bound.py section B, (i)-(ii)): compensating donors
and acceptors both scatter and partially cancel in the voltage offset.  So
f >= 1, and from (A2)-(A3):

    dV_c >= KAPPA * V_NP        x <= dV / (KAPPA * V_NP)                   (A4)

x is bounded ABOVE, so the screening of r is bounded above and the SURVIVING
asymmetry is bounded BELOW.  That is the useful direction: the bound cannot
promise a device more screening than it gets.

VALIDATION
----------
  E1  BITWISE C_ox INVARIANCE OF dV_c.  Over C_ox spanning six decades,
      dV_c from (A1) with n_c = KAPPA*f*N_A^(S)(C_ox) must equal the exactly
      known KAPPA*f*V_NP.  Not "to a tolerance" -- to a counted ulp floor.
  E2  BITWISE C_ox INVARIANCE OF x, the same statement one algebraic step on.
  E3  SYMMETRY THAT MUST GIVE EXACTLY ZERO: at dV = 0 the overdrive is the
      neutrality point, so x == 0.0 exactly and r(0) == r0 exactly, taken
      through the COMMITTED r_closed rather than a local copy.
  E4  THE TWO GATE STACKS, BOTH WAYS, REPORTED NOT ASSERTED.  Under fixed
      n_c reproduce 10-30's two numbers and their 1/C_ox ratio; under (A2)
      report the new ratio.  The old model is not claimed to be wrong
      everywhere -- it is reproduced exactly in its own regime and the
      disagreement is printed as a number.
  E5  THE OVERLAPPING LIMIT IS A REAL LIMIT.  At the C_ox for which
      N_A^(S)(C_ox) == N_CROSSOVER/KAPPA the two models must agree, and agree
      BITWISE, because there the fixed n_c IS the scaled one.  A model that
      disagreed with its predecessor everywhere would be a different model,
      not a correction.
  E6  MONOTONICITY IN f, STRICTLY: more impurities than the net charge
      implies means a larger dV_c, a smaller x and a LARGER surviving r.
  E7  MUTATION: break the cancellation by taking N_A^(S) proportional to
      C_ox**2 instead of C_ox.  E1 and E2 must both go RED.  A cancellation
      check that cannot fail is not evidence of a cancellation.
"""

import importlib
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ddm = importlib.import_module("graphene_density_dependent_mobility")
ilc = importlib.import_module("graphene_impurity_limited_crossover")
kib = importlib.import_module("graphene_kappa_impurity_bound")
bam = importlib.import_module("graphene_branch_asymmetric_mobility")
gfet = importlib.import_module("graphene_fet_model")

E_CHARGE = ddm.E_CHARGE
EPSILON_0 = ddm.EPSILON_0

# KAPPA, AND WHY THE EXCLUDED VALUE IS NOT THE ONE USED.  ilc.KAPPA is 2.4,
# which graphene_kappa_impurity_bound.py (2026-11-03) EXCLUDES: the bound from
# Zion's own two numbers is KAPPA <= 2.0889.  Section 2 below shows the
# cancellation is independent of KAPPA altogether, so nothing structural here
# turns on the choice; but the numbers printed should not quote an excluded
# value as if it stood, so the bound is used and 2.4 is reported beside it.
KAPPA = kib.kappa_bound()
KAPPA_EXCLUDED = ilc.KAPPA

# The shipped electron/hole pair, at the r0 = 0.10 of Chapter 4 Section 4.11.5.
R0_SHIPPED = 0.10
MU_L_E, MU_L_H = bam.mu_pair((1.0 + R0_SHIPPED) / (1.0 - R0_SHIPPED))

# Zion et al. Section III a, the same two inputs the 11-03 bound uses.
V_NP_ZION = 20.0                      # V
C_OX_ZION = EPSILON_0 * 3.9 / 300e-9  # F/m^2, 300 nm SiO2

C_OX_BACKGATE = gfet.C_ox if hasattr(gfet, "C_ox") else 3.84e-4
C_OX_TOPGATE = ddm.C_OX_TOPGATE_HFO2

passed = 0
failed = 0


def check(label, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
        print("  [PASS] %s" % label)
    else:
        failed += 1
        print("  [FAIL] %s" % label)
    if detail:
        print("         %s" % detail)


def ulp(a, b):
    return ddm.ulp_distance(a, b)


# ---------------------------------------------------------------------------
# The two readings of n_c
# ---------------------------------------------------------------------------
def n_a_surface(c_ox, v_np, power=1.0):
    """Net charged-acceptor density inferred from the neutrality offset.

    `power` exists only for E7's mutation: at power != 1 the C_ox dependence
    no longer cancels and the invariance checks must fail.
    """
    return (c_ox ** power) * v_np / E_CHARGE


def n_c_scaled(c_ox, v_np, f=1.0, power=1.0):
    """n_c = KAPPA * N_imp with N_imp = f * N_A^(S): the gate stack's own."""
    return KAPPA * f * n_a_surface(c_ox, v_np, power=power)


def dv_crossover(c_ox, n_c):
    """(A1): the gate overdrive at which |n| reaches n_c."""
    return E_CHARGE * n_c / c_ox


def dv_crossover_scaled(c_ox, v_np, f=1.0, power=1.0):
    return dv_crossover(c_ox, n_c_scaled(c_ox, v_np, f=f, power=power))


def x_of_overdrive(dv, c_ox, v_np, f=1.0, power=1.0):
    """x = |n|/n_c, with |n| = C_ox*dV/e."""
    n_abs = c_ox * dv / E_CHARGE
    return n_abs / n_c_scaled(c_ox, v_np, f=f, power=power)


def main():
    print("=" * 74)
    print("SCREENING OVERDRIVE INVARIANCE -- does the gate stack decide")
    print("whether electron/hole asymmetry is screened?")
    print("2026-11-07.  Closes the open item logged 2026-10-30.")
    print("=" * 74)

    print("\n" + "-" * 74)
    print("SECTION 1 -- inputs, restated from the committed modules")
    print("-" * 74)
    print("  KAPPA, bounded (graphene_kappa_impurity_bound)   : %.6f" % KAPPA)
    print("  KAPPA, the 11-02 value it EXCLUDES               : %.6f" % KAPPA_EXCLUDED)
    print("  r0, shipped pair (Section 4.11.5)                : %.6f"
          % ((MU_L_E - MU_L_H) / (MU_L_E + MU_L_H)))
    print("  V_NP  (Zion Section III a)                       : %.1f V" % V_NP_ZION)
    print("  C_ox  back gate, 90 nm SiO2                      : %.4e F/m^2" % C_OX_BACKGATE)
    print("  C_ox  top gate, 10 nm HfO2                       : %.4e F/m^2" % C_OX_TOPGATE)
    print("  C_ox ratio, top/back                             : %.2f" % (C_OX_TOPGATE / C_OX_BACKGATE))
    print("  N_CROSSOVER, the ONE value 10-30 used for both   : %.3e m^-2" % ddm.N_CROSSOVER)

    # -----------------------------------------------------------------
    print("\n" + "-" * 74)
    print("SECTION 2 -- E1/E2: the cancellation, bitwise, over six decades")
    print("-" * 74)
    c_grid = np.logspace(np.log10(C_OX_ZION) - 3, np.log10(C_OX_ZION) + 3, 601)
    exact = KAPPA * 1.0 * V_NP_ZION          # the exactly known value, (A2)
    exact_bounded = exact
    dvc = dv_crossover_scaled(c_grid, V_NP_ZION)
    d1 = ulp(dvc, np.full_like(c_grid, exact))
    print("  exactly known dV_c = KAPPA * V_NP = %.12f V" % exact)
    print("  max |dV_c - KAPPA*V_NP| over the grid : %.3e V" % np.max(np.abs(dvc - exact)))
    print("  max ulp distance                      : %d" % int(d1.max()))
    check("E1 dV_c is independent of C_ox to a counted ulp floor",
          int(d1.max()) <= 2,
          "601 points over six decades of C_ox, max %d ulp from the exactly "
          "known KAPPA*V_NP.  A 1/C_ox dependence would be 6 decades wide "
          "here, not 2 ulp." % int(d1.max()))

    dv_probe = 1.0                            # V of overdrive
    xs = x_of_overdrive(dv_probe, c_grid, V_NP_ZION)
    x_exact = dv_probe / (KAPPA * V_NP_ZION)  # (A3)
    d2 = ulp(xs, np.full_like(c_grid, x_exact))
    print("  exactly known x at dV = 1 V : %.12f" % x_exact)
    print("  max ulp distance            : %d" % int(d2.max()))
    check("E2 x = |n|/n_c is independent of C_ox to a counted ulp floor",
          int(d2.max()) <= 2,
          "the same statement one algebraic step on: C_ox enters |n| and "
          "n_c with the same power, so it cannot appear in their ratio")

    # E2b -- the cancellation does not depend on KAPPA, so the result does not
    # inherit the 11-02/11-03 dispute about KAPPA's value.  This matters for
    # what the conclusion is entitled to say: the 56x spread is an artefact
    # whatever KAPPA turns out to be.
    k_grid = np.logspace(-2, 2, 401)
    ratio = np.array([dv_crossover_scaled(C_OX_TOPGATE, V_NP_ZION)
                      / dv_crossover_scaled(C_OX_BACKGATE, V_NP_ZION)])
    kk = np.array([(E_CHARGE * (kk_ * 1.0 * n_a_surface(C_OX_BACKGATE, V_NP_ZION))
                    / C_OX_BACKGATE) / (kk_ * V_NP_ZION) for kk_ in k_grid])
    d2b = ulp(kk, np.ones_like(kk))
    print("  dV_c/(KAPPA*V_NP) over KAPPA in [1e-2, 1e2], max ulp from 1 : %d"
          % int(d2b.max()))
    check("E2b the cancellation is independent of KAPPA, so the conclusion "
          "does not inherit the dispute about KAPPA's value",
          int(d2b.max()) <= 2 and ulp(ratio, np.ones(1)).max() <= 2,
          "401 values of KAPPA over four decades; dV_c/(KAPPA*V_NP) is 1 to "
          "%d ulp, and the top/back stack ratio is 1 to %d ulp.  Whether "
          "KAPPA is 2.4 or the bounded 2.0889 or neither, the two stacks "
          "still coincide." % (int(d2b.max()), int(ulp(ratio, np.ones(1)).max())))

    # -----------------------------------------------------------------
    print("\n" + "-" * 74)
    print("SECTION 3 -- E3: the symmetry that must give exactly zero")
    print("-" * 74)
    x0 = x_of_overdrive(0.0, C_OX_BACKGATE, V_NP_ZION)
    print("  x at dV = 0                : %r" % x0)
    check("E3a at the neutrality point x is EXACTLY zero", x0 == 0.0,
          "not 'small' -- the float is 0.0, because |n| is 0.0 and n_c is "
          "finite and nonzero")
    global r0_ref
    r_at_zero = ddm.r_closed(0.0, MU_L_E, MU_L_H,
                             n_c=n_c_scaled(C_OX_BACKGATE, V_NP_ZION))
    r0_ref = (MU_L_E - MU_L_H) / (MU_L_E + MU_L_H)
    print("  r(0) through the COMMITTED r_closed : %.17g" % r_at_zero)
    print("  r0 computed directly from the pair  : %.17g" % r0_ref)
    check("E3b r(0) == r0 exactly, through the committed closed form",
          r_at_zero == r0_ref,
          "the screening factor is 1/(1+(1-r0^2)x) and x is exactly 0 here, "
          "so this is an identity and not a numerical agreement.  Taken "
          "through graphene_density_dependent_mobility.r_closed, not a local "
          "copy of it")

    # -----------------------------------------------------------------
    print("\n" + "-" * 74)
    print("SECTION 4 -- E4: the two gate stacks, under BOTH readings")
    print("-" * 74)
    fixed_back = dv_crossover(C_OX_BACKGATE, ddm.N_CROSSOVER)
    fixed_top = dv_crossover(C_OX_TOPGATE, ddm.N_CROSSOVER)
    scal_back = dv_crossover_scaled(C_OX_BACKGATE, V_NP_ZION)
    scal_top = dv_crossover_scaled(C_OX_TOPGATE, V_NP_ZION)
    print("  %-26s %12s %12s" % ("", "dV_c fixed n_c", "dV_c scaled"))
    print("  %-26s %12.4f %12.4f" % ("90 nm SiO2 back gate", fixed_back, scal_back))
    print("  %-26s %12.4f %12.4f" % ("10 nm HfO2 top gate", fixed_top, scal_top))
    print("  ratio back/top             %12.2f %12.6f"
          % (fixed_back / fixed_top, scal_back / scal_top))
    print()
    print("  REPORTED, NOT ASSERTED: the fixed-n_c reading spreads the two")
    print("  stacks by %.1fx; the scaled reading spreads them by %.6fx."
          % (fixed_back / fixed_top, scal_back / scal_top))
    check("E4a the fixed-n_c model reproduces 10-30's two numbers",
          abs(fixed_back - 12.5) < 0.3 and abs(fixed_top - 0.22) < 0.02,
          "back %.3f V against 10-30's ~12.5 V, top %.4f V against its "
          "~0.22 V -- the predecessor is reproduced in its own terms before "
          "being contradicted" % (fixed_back, fixed_top))
    check("E4b under the scaled reading the two stacks are EXACTLY equal",
          scal_back == scal_top,
          "difference %r V.  Both equal KAPPA*V_NP = %.6f V, so the "
          "difference is 0.0 and not merely small"
          % (scal_back - scal_top, exact))
    check("E4c and the fixed-n_c ratio is exactly the inverse C_ox ratio",
          ulp(np.array([fixed_back / fixed_top]),
              np.array([C_OX_TOPGATE / C_OX_BACKGATE])).max() <= 2,
          "%.4f against %.4f -- which is what identifies the 56x as a "
          "statement about the oxide and not about the graphene"
          % (fixed_back / fixed_top, C_OX_TOPGATE / C_OX_BACKGATE))

    # -----------------------------------------------------------------
    print("\n" + "-" * 74)
    print("SECTION 5 -- E5: the overlapping limit, and it is bitwise")
    print("-" * 74)
    # Where does the scaled n_c equal the fixed one?  KAPPA*C_ox*V_NP/e = n_c
    c_ox_cross = ddm.N_CROSSOVER * E_CHARGE / (KAPPA * V_NP_ZION)
    print("  C_ox at which the scaled n_c EQUALS N_CROSSOVER : %.6e F/m^2" % c_ox_cross)
    a = dv_crossover(c_ox_cross, ddm.N_CROSSOVER)
    b = dv_crossover_scaled(c_ox_cross, V_NP_ZION)
    print("  dV_c, fixed-n_c model   : %.17g V" % a)
    print("  dV_c, scaled model      : %.17g V" % b)
    print("  ulp distance            : %d" % int(ulp(np.array([a]), np.array([b])).max()))
    check("E5 the two models agree BITWISE where their n_c coincide",
          int(ulp(np.array([a]), np.array([b])).max()) <= 1,
          "the correction contains its predecessor: at this one C_ox the "
          "fixed n_c IS the scaled n_c, so disagreement here would mean the "
          "new model is a different model rather than a correction")

    # -----------------------------------------------------------------
    print("\n" + "-" * 74)
    print("SECTION 6 -- E6: the direction of the inequality, f >= 1")
    print("-" * 74)
    fs = np.array([1.0, 1.5, 2.0, 4.0])
    print("  %6s %12s %12s %12s" % ("f", "dV_c (V)", "x at 1 V", "r at 1 V"))
    rows = []
    for f in fs:
        dvc_f = dv_crossover_scaled(C_OX_BACKGATE, V_NP_ZION, f=f)
        x_f = x_of_overdrive(1.0, C_OX_BACKGATE, V_NP_ZION, f=f)
        n_abs = C_OX_BACKGATE * 1.0 / E_CHARGE
        r_f = ddm.r_closed(n_abs, MU_L_E, MU_L_H,
                           n_c=n_c_scaled(C_OX_BACKGATE, V_NP_ZION, f=f))
        rows.append((dvc_f, x_f, r_f))
        print("  %6.2f %12.4f %12.6f %12.6f" % (f, dvc_f, x_f, r_f))
    check("E6 dV_c strictly increases, x strictly decreases and the "
          "surviving r strictly increases with f",
          all(rows[i][0] < rows[i + 1][0] for i in range(len(rows) - 1)) and
          all(rows[i][1] > rows[i + 1][1] for i in range(len(rows) - 1)) and
          all(rows[i][2] < rows[i + 1][2] for i in range(len(rows) - 1)),
          "f = 1 is the LOWER bound on impurity count, so it is the UPPER "
          "bound on screening and the LOWER bound on surviving asymmetry -- "
          "the safe direction for a device claim")

    # -----------------------------------------------------------------
    print("\n" + "-" * 74)
    print("SECTION 7 -- E7: mutation.  Break the cancellation.")
    print("-" * 74)
    dvc_m = dv_crossover_scaled(c_grid, V_NP_ZION, power=2.0)
    dm1 = ulp(dvc_m, np.full_like(c_grid, dvc_m[0]))
    xs_m = x_of_overdrive(dv_probe, c_grid, V_NP_ZION, power=2.0)
    dm2 = ulp(xs_m, np.full_like(c_grid, xs_m[0]))
    print("  with N_A^(S) ~ C_ox**2, dV_c spans %.3e .. %.3e V"
          % (dvc_m.min(), dvc_m.max()))
    print("  max ulp spread, dV_c : %d" % int(dm1.max()))
    print("  max ulp spread, x    : %d" % int(dm2.max()))
    check("E7 the mutant is DETECTED -- both invariance checks go red",
          int(dm1.max()) > 2 and int(dm2.max()) > 2,
          "a cancellation check that passes when the cancellation is removed "
          "would be evidence of nothing.  Spreads %d and %d ulp."
          % (int(dm1.max()), int(dm2.max())))

    # -----------------------------------------------------------------
    print("\n" + "-" * 74)
    print("SECTION 8 -- the device consequence, which is a REVERSAL")
    print("-" * 74)
    # 10-30's conclusion was operational: the thin-oxide top gate reaches n_c
    # at ~0.22 V, so a process device sits far past the crossover and the
    # asymmetry is screened away.  Under (A2) the crossover overdrive is
    # KAPPA*V_NP for BOTH stacks, and that is larger than either device's
    # operating overdrive -- so neither stack screens.
    dv_op = 1.0
    for name, c_ox in (("90 nm SiO2 back gate", C_OX_BACKGATE),
                       ("10 nm HfO2 top gate", C_OX_TOPGATE)):
        n_abs = c_ox * dv_op / E_CHARGE
        x_fixed = n_abs / ddm.N_CROSSOVER
        x_scal = x_of_overdrive(dv_op, c_ox, V_NP_ZION)
        r_fixed = ddm.r_closed(n_abs, MU_L_E, MU_L_H, n_c=ddm.N_CROSSOVER)
        r_scal = ddm.r_closed(n_abs, MU_L_E, MU_L_H,
                              n_c=n_c_scaled(c_ox, V_NP_ZION))
        print("  %-22s at %.1f V overdrive" % (name, dv_op))
        print("      x  fixed n_c %10.4f   scaled %10.4f   factor %8.1f"
              % (x_fixed, x_scal, x_fixed / x_scal))
        print("      r  fixed n_c %10.6f   scaled %10.6f   of r0 = %.2f"
              % (r_fixed, r_scal, r0_ref))
    n_top = C_OX_TOPGATE * dv_op / E_CHARGE
    r_top_scal = ddm.r_closed(n_top, MU_L_E, MU_L_H,
                              n_c=n_c_scaled(C_OX_TOPGATE, V_NP_ZION))
    check("E8 the top-gated stack does NOT screen the asymmetry away, "
          "contradicting the 10-30 conclusion",
          r_top_scal > 0.9 * r0_ref,
          "r = %.6f against r0 = %.6f -- %.1f%% of the unscreened value, "
          "where the fixed-n_c reading gives %.6f.  10-30 concluded the "
          "asymmetry is 'small in the thin-oxide top-gated devices a "
          "semiconductor process would build'; under the scaled n_c it is "
          "essentially unscreened there.  This is a reversal of that "
          "conclusion, not a refinement of its magnitude."
          % (r_top_scal, r0_ref, 100.0 * r_top_scal / r0_ref,
             ddm.r_closed(n_top, MU_L_E, MU_L_H, n_c=ddm.N_CROSSOVER)))
    check("E8b and the crossover overdrive exceeds BOTH stacks' operating "
          "range, which is why",
          exact_bounded > 10.0,
          "KAPPA*V_NP = %.2f V of overdrive is needed to reach n_c, for "
          "either stack.  10-30's table put the top gate at 0.22 V, so the "
          "crossover looked reachable; it is not." % exact_bounded)

    print("\n" + "-" * 74)
    print("SECTION 9 -- the figure")
    print("-" * 74)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].loglog(c_grid, dv_crossover(c_grid, ddm.N_CROSSOVER),
                 label=r"fixed $n_c$ (10-30): $\propto 1/C_{ox}$")
    ax[0].loglog(c_grid, dv_crossover_scaled(c_grid, V_NP_ZION),
                 label=r"scaled $n_c=\kappa N_A^{(S)}$: $\kappa V_{NP}$")
    for c, lab in ((C_OX_BACKGATE, "90 nm SiO$_2$"), (C_OX_TOPGATE, "10 nm HfO$_2$")):
        ax[0].axvline(c, ls=":", lw=0.8, color="grey")
        ax[0].text(c, dv_crossover(c_grid, ddm.N_CROSSOVER).max(), lab,
                   rotation=90, va="top", ha="right", fontsize=7)
    ax[0].set_xlabel(r"$C_{ox}$ (F m$^{-2}$)")
    ax[0].set_ylabel(r"crossover overdrive $\Delta V_c$ (V)")
    ax[0].set_title("(a) the 56x spread is the fixed-$n_c$ assumption")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3, which="both")

    dv_axis = np.linspace(0.0, 1.25 * exact, 400)   # past KAPPA*V_NP, so the marker is in-plot
    n_axis = C_OX_BACKGATE * dv_axis / E_CHARGE
    for f, ls in ((1.0, "-"), (2.0, "--"), (4.0, ":")):
        r_curve = ddm.r_closed(n_axis, MU_L_E, MU_L_H,
                               n_c=n_c_scaled(C_OX_BACKGATE, V_NP_ZION, f=f))
        ax[1].plot(dv_axis, r_curve, ls, label=r"$f=%.0f$" % f)
    ax[1].axvline(exact, color="k", lw=0.9, ls="-.")
    ax[1].annotate(r"$\kappa V_{NP}=%.1f$ V" % exact, xy=(exact, R0_SHIPPED * 0.62),
                   xytext=(-6, 0), textcoords="offset points",
                   rotation=90, va="center", ha="right", fontsize=8)
    ax[1].axhline(R0_SHIPPED, color="grey", lw=0.7, ls=":")
    ax[1].text(0.5, R0_SHIPPED, r" $r_0$", fontsize=8, va="bottom")
    ax[1].set_xlabel(r"gate overdrive $\Delta V$ (V)")
    ax[1].set_ylabel(r"surviving asymmetry $r$")
    ax[1].set_title(r"(b) $f=1$ bounds screening ABOVE, so $r$ BELOW")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("screening_overdrive_invariance.png", dpi=150)
    print("  wrote screening_overdrive_invariance.png")

    print("\n" + "=" * 74)
    print("TOTAL: %d passed, %d failed" % (passed, failed))
    print("=" * 74)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
