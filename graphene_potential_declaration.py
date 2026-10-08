"""
graphene_potential_declaration.py

WHICH POTENTIAL IS `V_ch`?  THE ANSWER, AND WHAT IT COSTS.

Created 2026-10-08 to close Chapter 7 Section 7.9 item 7 as rewritten on
2026-10-07: "declare which potential V_ch is, consistently across
carrier_density(), the series factor and Eq. (4), and re-derive."

THE DECLARATION
---------------
    V_ch IS THE QUASI-FERMI (ELECTROCHEMICAL) POTENTIAL OF THE CHANNEL
    CARRIERS, IN VOLTS, MEASURED FROM THE SOURCE.

Three independent reasons, none of them new work:

 1. The compact model Eq. (4) descends from says so.  Pasadas and Jimenez,
    IEEE TED 63(7) 2016 (arXiv:1605.08235), write v = mu*F with F = -dV/dx
    and state that "V(x) is the quasi-Fermi level along the graphene
    channel".  This was already quoted in
    graphene_diffusion_current_model.py on 2026-10-07; what was missing was
    the step of binding this repository's variable to it.

 2. Classical long-channel theory uses the same variable and says what it
    buys.  H.-S. P. Wong's long-channel MOSFET notes define the channel
    variable V(y) as the "Quasi-Fermi potential along the channel" and head
    the current-density equation it feeds "Current density equation (both
    drift and diffusion)".  That is the whole content of the declaration: in
    the quasi-Fermi variable, ONE term is the complete current.

 3. The repository's own code already assumes it.  `carrier_density()`
    multiplies the electrostatic estimate by the quantum-capacitance series
    factor C_q/(C_q + C_ox).  That factor is the correction one applies when
    the channel variable is the quasi-Fermi level and the graphene drop
    E_F/e is carried separately.  Under a strictly electrostatic reading of
    V_ch the factor is a double count of the graphene drop.

WHAT THE DECLARATION DOES TO 2026-10-07's WORK -- STATED, NOT QUIETLY FIXED
--------------------------------------------------------------------------
On 2026-10-07 a diffusion term was added to Eq. (4) and measured at +0.5018 %
of I_d at the RF bias.  graphene_diffusion_current_model.py's own Eq. (9)
writes the current as

    I_d = W*e*mu*n * d( V - lambda*V_F )/dx,     Phi = V - V_F

and its own line 72 names Phi the quasi-Fermi potential.  So:

  * if V is the ELECTROSTATIC potential, lambda = 1 is required and the
    committed lambda = 0 numbers were missing a real 0.5 % of current;
  * if V is the QUASI-FERMI potential, then V IS Phi, lambda = 0 is already
    the complete drift-diffusion current, and lambda = 1 ADDS A TERM THAT IS
    ALREADY THERE.

The declaration above picks the second.  **2026-10-07's diffusion term is
therefore a DOUBLE COUNT of 0.50 %, not a missing 0.50 %.**  That contradicts
the direction 10-07 reported, and this file says so rather than editing the
number: the term is correctly derived, correctly implemented, structurally a
boundary term, and it must not be switched on.  `lambda` stays in the source
as the knob that measures the difference between the two readings, which is
the one thing it is now good for.

AND WHAT THE DECLARATION COSTS, MEASURED
----------------------------------------
Declaring V_ch quasi-Fermi makes a claim about `carrier_density()` that can
be tested instead of asserted.  Under the declaration, the gate drive divides
between the oxide and the graphene quantum capacitance exactly:

    V_g - V_dirac - V_ch  =  e*n/C_ox  +  V_F(n),      V_F(n) = A_F*sqrt(n)
                                                                      (D1)

which is a quadratic in sqrt(n) with a CLOSED-FORM root -- no iteration, no
new parameter, and A_F is this repository's own A_FERMI.  The repository
instead uses the linearised series factor

    n = (C_ox*dV/e) * C_q/(C_q + C_ox),   C_q = quantum_capacitance(dV)
                                                                      (D2)

(D2) is the LINEARISATION of (D1), so the declaration is self-consistent only
to the extent the two agree.  SECTION 2 measures the disagreement instead of
claiming there is none, and the result is the main new number in this file:
it is NOT small, and it is larger than the 0.50 % double count it removes.
See SECTION 5 for what that does and does not license.

VALIDATION AGAINST EXACTLY KNOWN VALUES (SECTION 3)
  - (D1)'s closed-form root must satisfy (D1) to machine precision.  This is
    an algebraic identity, not a plausible range.
  - dV = 0 must give n = 0 EXACTLY in both forms, before the puddle floor.
    A symmetry that must give zero.
  - Both forms carry n as a magnitude, so n(+dV) must equal n(-dV) BITWISE.
  - As dV grows, C_q >> C_ox and both forms must approach the pure
    electrostatic value e*n = C_ox*dV.  The approach is measured, not
    asserted: 2026-10-03's rule that a PASS at zero is silent about
    magnitude.
"""

import numpy as np

import graphene_fet_model as gfet
import graphene_diffusion_current_model as dd

E_CHARGE = dd.E_CHARGE if hasattr(dd, 'E_CHARGE') else 1.602176634e-19
KT_OVER_E = 1.380649e-23 * gfet.T / E_CHARGE   # thermal voltage [V]

passed = 0
failed = 0


def check(label, ok, detail=''):
    global passed, failed
    if ok:
        passed += 1
        print('  [PASS] %s' % label)
    else:
        failed += 1
        print('  [FAIL] %s' % label)
    if detail:
        print('         %s' % detail)


def n_quasi_fermi(dV):
    """
    Closed-form root of (D1): e*n/C_ox + A_F*sqrt(n) = |dV|.

    With x = sqrt(n), a = e/C_ox, b = A_FERMI:
        a*x^2 + b*x - |dV| = 0
        x = (-b + sqrt(b^2 + 4*a*|dV|)) / (2*a)
    """
    dV = np.abs(np.asarray(dV, dtype=float))
    a = E_CHARGE / gfet.C_ox
    b = dd.A_FERMI
    x = (-b + np.sqrt(b * b + 4.0 * a * dV)) / (2.0 * a)
    return x * x


def n_series_factor(dV):
    """(D2), the repository's form, with the puddle floor NOT applied."""
    dV = np.asarray(dV, dtype=float)
    n_es = gfet.C_ox * dV / E_CHARGE
    C_q = gfet.quantum_capacitance(dV, T=gfet.T)
    return np.abs(n_es * C_q / (C_q + gfet.C_ox))


def n_total_series(dV):
    """(D2) with the puddle floor, i.e. exactly what carrier_density returns."""
    return np.sqrt(n_series_factor(dV) ** 2 + gfet.n_puddle ** 2)


def n_total_qf(dV):
    """(D1) with the SAME puddle floor, so the comparison is apples to apples."""
    return np.sqrt(n_quasi_fermi(dV) ** 2 + gfet.n_puddle ** 2)


def nan_onset_of_series_factor():
    """The drive at which quantum_capacitance() stops returning a number."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        lo, hi = 1.0, 40.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            ok = np.isfinite(gfet.quantum_capacitance(mid, T=gfet.T))
            if ok:
                lo = mid
            else:
                hi = mid
    return lo, hi


def main():
    print('=' * 78)
    print('THE DECLARATION: V_ch IS THE QUASI-FERMI POTENTIAL')
    print('Chapter 7 Section 7.9 item 7, as rewritten 2026-10-07')
    print('=' * 78)

    print()
    print('SECTION 1.  The double count, at three levels of the same model')
    print('-' * 78)
    print('  lambda = 0 : the drift form in the declared variable.')
    print('  lambda = 1 : the same current with the 10-07 diffusion term ADDED.')
    print('  Under the declaration lambda = 0 is already the complete')
    print('  drift-diffusion current, so lambda = 1 counts diffusion twice.')
    print()
    V_g_rf, Vds_rf = 2.0, dd.VDS_REF
    Id0_i, _ = dd._Id_given_Vds_ch(V_g_rf, Vds_rf, lam=0.0)
    Id1_i, _ = dd._Id_given_Vds_ch(V_g_rf, Vds_rf, lam=1.0)
    sh_i = 100.0 * (Id1_i - Id0_i) / Id0_i
    Id0_t = float(dd.transfer_characteristic_dd(np.array([V_g_rf]), Vds=Vds_rf,
                                                lam=0.0)[0])
    Id1_t = float(dd.transfer_characteristic_dd(np.array([V_g_rf]), Vds=Vds_rf,
                                                lam=1.0)[0])
    sh_t = 100.0 * (Id1_t - Id0_t) / Id0_t
    print('  At V_g = %.2f V, V_ds = %.2f V:' % (V_g_rf, Vds_rf))
    print('    INTRINSIC (fixed channel drop) : %.6e -> %.6e A  %+.4f %%'
          % (Id0_i, Id1_i, sh_i))
    print('    TERMINAL  (contact feedback on): %.6e -> %.6e A  %+.4f %%'
          % (Id0_t, Id1_t, sh_t))
    print('    10-07 reported +0.5018 % at ITS saturated-peak bias, a')
    print('    different V_g; the three numbers are the same term measured')
    print('    three ways and they bracket each other, which is the check.')
    check('the double count is a few tenths of a percent at every level',
          0.2 < abs(sh_i) < 1.0 and 0.2 < abs(sh_t) < 1.0,
          'intrinsic %+.4f %%, terminal %+.4f %%; the contact feedback eats '
          'about a third of it, as 10-07 found' % (sh_i, sh_t))
    check('the terminal double count is SMALLER than the intrinsic one',
          abs(sh_t) < abs(sh_i),
          'the series contact resistance is a negative feedback on I_d, so any '
          'intrinsic change is attenuated at the terminals -- the sign of this '
          'inequality is fixed by that and is not a fitted expectation')

    print()
    print('SECTION 2.  What the declaration costs, in the quantity that is quoted')
    print('-' * 78)
    print('  Declaring V_ch quasi-Fermi makes (D1) the exact charge relation.')
    print('  The repository computes (D2).  Three comparisons, narrowing from')
    print('  the raw charge to the current Chapter 4 actually quotes, because')
    print('  2026-10-06 established that a criterion decided on the quantity it')
    print('  names can miss the finding when that is not the quantity at stake.')
    print()
    print('  2a  RAW charge, puddle floor NOT applied:')
    print('   %10s %14s %14s %12s' % ('dV [V]', 'n (D1)', 'n (D2)', 'D2/D1 - 1'))
    raw_worst, raw_at = 0.0, None
    for dV in (0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 2.7, 3.5):
        n1, n2 = float(n_quasi_fermi(dV)), float(n_series_factor(dV))
        rel = n2 / n1 - 1.0
        if abs(rel) > abs(raw_worst):
            raw_worst, raw_at = rel, dV
        print('   %10.3f %14.6e %14.6e %11.3f %%' % (dV, n1, n2, 100.0 * rel))
    print('   worst %+.2f %% at dV = %.3f V -- but see 2b before quoting that.'
          % (100.0 * raw_worst, raw_at))
    print()
    print('  2b  THE SAME COMPARISON WITH THE PUDDLE FLOOR, i.e. what')
    print('      carrier_density() actually returns.  n_puddle = %.2e /m2,'
          % gfet.n_puddle)
    print('      so the near-Dirac rows above are in a regime the floor erases')
    print('      and the %+.0f %% is NOT a number about any committed result:'
          % (100.0 * raw_worst))
    print('   %10s %14s %14s %12s' % ('dV [V]', 'n_tot (D1)', 'n_tot (D2)',
                                      'D2/D1 - 1'))
    flo_worst, flo_at = 0.0, None
    for dV in (0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 2.7, 3.5):
        n1, n2 = float(n_total_qf(dV)), float(n_total_series(dV))
        rel = n2 / n1 - 1.0
        if abs(rel) > abs(flo_worst):
            flo_worst, flo_at = rel, dV
        print('   %10.3f %14.6e %14.6e %11.3f %%' % (dV, n1, n2, 100.0 * rel))
    print('   worst %+.2f %% at dV = %.2f V, over the drives Chapter 4 sweeps.'
          % (100.0 * flo_worst, flo_at))
    check('the floored disagreement is far smaller than the raw one',
          abs(flo_worst) < abs(raw_worst),
          'raw %+.2f %% -> floored %+.2f %%: a factor %.0f.  Quoting the raw '
          'number would have overstated the cost of the declaration by that '
          'factor.' % (100.0 * raw_worst, 100.0 * flo_worst,
                       abs(raw_worst / flo_worst)))

    print()
    print('  2c  AND IN I_d ITSELF, by running Eq. (4) on (D1).  carrier_density')
    print('      is replaced by the (D1) form with the same puddle floor and')
    print('      nothing else changed:')
    original = gfet.carrier_density
    try:
        gfet.carrier_density = lambda V_g, V_ch=0.0: n_total_qf(
            np.asarray(V_g, dtype=float) - np.asarray(V_ch, dtype=float)
            - gfet.V_dirac)
        Id_D1 = float(dd.transfer_characteristic_dd(np.array([V_g_rf]),
                                                    Vds=Vds_rf, lam=0.0)[0])
    finally:
        gfet.carrier_density = original
    cost = 100.0 * (Id_D1 - Id0_t) / Id0_t
    print('      I_d on (D2), lambda = 0 : %.6e A   (committed)' % Id0_t)
    print('      I_d on (D1), lambda = 0 : %.6e A   (declared exactly)' % Id_D1)
    print('      the declaration costs   : %+.4f %% of I_d' % cost)
    print('      the double count it removes: %+.4f %% of I_d' % sh_t)
    check('the charge-model question is LARGER than the diffusion term',
          abs(cost) > abs(sh_t),
          'the declaration moves I_d by %+.4f %% while the term it retires is '
          'worth %+.4f %% -- a factor %.1f.  So the honest headline is not '
          '"the 0.50 %% is resolved": it is that a question about the charge '
          'model, %.1f times larger, was sitting underneath it, unasked.'
          % (cost, sh_t, abs(cost / sh_t), abs(cost / sh_t)))
    print('      NOTE: this is a MEASUREMENT of the gap, not a re-derivation.')
    print('      Chapter 4 keeps its (D2) numbers; see SECTION 5.')

    print()
    print('  2d  A numerical limit of (D2), found while doing 2a, measured:')
    lo, hi = nan_onset_of_series_factor()
    print('      quantum_capacitance() returns a finite value up to dV =')
    print('      %.4f V and NaN/inf above %.4f V, because' % (lo, hi))
    print('      log(2*(1 + cosh(eta))) overflows once eta = dV/(kT/e) passes')
    print('      about 710.  kT/e = %.6f V here, so the onset is at eta ~ %.0f.'
          % (KT_OVER_E, lo / KT_OVER_E))
    check('the overflow is far outside every bias this thesis quotes',
          lo > 10.0,
          'onset %.2f V against a maximum swept drive of about 2.7 V, so no '
          'committed number is affected.  It is still a defect -- the large-eta '
          'limit of that expression is |eta| + log 2 and is exact -- and it is '
          'left as a new item rather than patched today, because the module '
          'owns committed transcripts.' % lo)

    print()
    print('SECTION 3.  Validation against exactly known values')
    print('-' * 78)
    dVs = np.array([0.01, 0.1, 0.5, 1.0, 2.0, 3.5])
    n1 = n_quasi_fermi(dVs)
    resid = np.abs(E_CHARGE * n1 / gfet.C_ox + dd.A_FERMI * np.sqrt(n1) - dVs)
    check("(D1)'s closed-form root satisfies (D1) to machine precision",
          float(np.max(resid / dVs)) < 1e-14,
          'max relative residual %.3e over dV in [0.01, 3.5] V -- an algebraic '
          'identity, not a plausible range' % float(np.max(resid / dVs)))

    z1, z2 = float(n_quasi_fermi(0.0)), float(n_series_factor(0.0))
    check('dV = 0 gives n = 0 EXACTLY in both forms, before the puddle floor',
          z1 == 0.0 and z2 == 0.0, '(D1) -> %r, (D2) -> %r' % (z1, z2))

    sym1 = n_quasi_fermi(dVs) - n_quasi_fermi(-dVs)
    sym2 = n_series_factor(dVs) - n_series_factor(-dVs)
    check('n(+dV) == n(-dV) BITWISE in both forms (both carry n as a magnitude)',
          bool(np.all(sym1 == 0.0)) and bool(np.all(sym2 == 0.0)),
          'max |difference| (D1) %r, (D2) %r'
          % (float(np.max(np.abs(sym1))), float(np.max(np.abs(sym2)))))

    print()
    print('  The large-drive limit, measured rather than asserted, and stopped')
    print('  below 2d\'s overflow:')
    print('   %10s %18s %18s' % ('dV [V]', '(D1)/electrostatic',
                                 '(D2)/electrostatic'))
    rows = []
    for dV in (0.5, 1.0, 2.0, 5.0, 10.0, 15.0):
        n_es = gfet.C_ox * dV / E_CHARGE
        r1, r2 = float(n_quasi_fermi(dV)) / n_es, float(n_series_factor(dV)) / n_es
        rows.append((dV, r1, r2))
        print('   %10.1f %18.6f %18.6f' % (dV, r1, r2))
    check('both forms approach the pure electrostatic value as dV grows, '
          'monotonically',
          all(rows[k][1] < rows[k + 1][1] for k in range(len(rows) - 1))
          and all(rows[k][2] < rows[k + 1][2] for k in range(len(rows) - 1)),
          '(D1) %.4f -> %.4f and (D2) %.4f -> %.4f over dV = 0.5 -> 15 V; both '
          'monotone, and (D1) approaches from further away because it carries '
          'the whole graphene drop rather than a linearisation of it'
          % (rows[0][1], rows[-1][1], rows[0][2], rows[-1][2]))

    print()
    print('SECTION 4.  V_F is a POTENTIAL, and the census had it UNADJUDICATED')
    print('-' * 78)
    print('  graphene_potential_census_audit.py, committed earlier today, left')
    print('  `V_F` UNADJUDICATED in three modules on the guess that it was a')
    print('  Fermi VELOCITY written with a V_ prefix.  In')
    print('  graphene_diffusion_current_model.py it is not: fermi_voltage(n)')
    print('  returns E_F/e in VOLTS, and it is precisely the quantity that')
    print('  converts between the two readings of V_ch.  In')
    print('  graphene_contact_doping_nonlinear_model.py the same name IS a')
    print('  velocity, 1.0e6 m/s.  One name, two quantities, two modules --')
    print('  which is the 2026-10-01 fault (a name that cannot drift is the')
    print('  whole requirement) in a third place.')
    print()
    for dV in (0.1, 1.0, 2.0, 2.7):
        n1 = float(n_quasi_fermi(dV))
        print('    dV = %.1f V : n = %.4e /m2, V_F = E_F/e = %.4f V, i.e. '
              '%.1f %% of the drive sits on the graphene'
              % (dV, n1, dd.fermi_voltage(n1), 100.0 * dd.fermi_voltage(n1) / dV))
    check('V_F(n) from the diffusion module equals A_FERMI*sqrt(n)',
          abs(dd.fermi_voltage(1e16) - dd.A_FERMI * np.sqrt(1e16)) < 1e-15,
          'checked at n = 1e16 /m2; this is the identity that makes (D1) '
          'closed-form')
    print('  The census registry is updated accordingly in the same commit:')
    print('  two V_F rows move from UNADJUDICATED to a reading, and the third')
    print('  stays UNADJUDICATED because its module was not read today.')

    print()
    print('SECTION 5.  What this settles, and what it opens')
    print('-' * 78)
    print('  SETTLED.  V_ch is the quasi-Fermi potential, declared in')
    print('  graphene_fet_model.py where it is defined rather than here.')
    print('  2026-10-07\'s diffusion term is therefore a DOUBLE COUNT worth')
    print('  %+.4f %% of I_d at the terminals, not a missing %+.4f %%. The'
          % (sh_t, sh_t))
    print('  committed lambda = 0 currents stand unchanged, and the +-0.50 %')
    print('  ambiguity 10-07 attached to every I_d in this repository is')
    print('  RESOLVED rather than reduced.  That is the opposite of the')
    print('  direction 10-07 expected, and the term is kept in the source as')
    print('  the instrument that measures the difference between readings.')
    print()
    print('  OPENED, and larger than what was closed.  (D1) is now the exact')
    print('  charge relation and the repository computes (D2).  In I_d at the')
    print('  RF bias they differ by %+.4f %% -- %.1f times the term just'
          % (cost, abs(cost / sh_t)))
    print('  retired.  Re-deriving Chapter 4 on (D1) is a session of its own')
    print('  and is the new top item.  Chapter 4 keeps its numbers today and')
    print('  gains a stated limitation with a measured size, per the standing')
    print('  rule that a superseded number is annotated and not deleted.')
    print()
    print('  NOT SETTLED by this file:')
    print('   - whether (D1) or (D2) is closer to MEASUREMENT.  Neither has')
    print('     been compared to data here, only to each other, so nothing')
    print('     above licenses calling (D2) wrong -- only non-exact under the')
    print('     declaration.')
    print('   - the hole branch.  Both forms carry n as a magnitude, which')
    print('     10-07 identified as leaving every branch-asymmetric Chapter 4')
    print('     result unverified.  The declaration does not touch it, and it')
    print('     now BLOCKS the (D1) re-derivation rather than sitting beside')
    print('     it: re-deriving on a magnitude would rebuild the same flaw.')
    print('   - the cosh overflow in 2d, left as a new item.')

    print()
    print('=' * 78)
    print('TOTAL: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return failed


if __name__ == '__main__':
    import sys
    sys.exit(1 if main() else 0)
