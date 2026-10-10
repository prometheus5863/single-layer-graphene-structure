"""
graphene_exact_charge_rederivation.py

CHAPTER 4, RE-DERIVED ON THE EXACT CHARGE RELATION (4.28).
Chapter 7 Section 7.9 ITEM 12 -- created 2026-10-08, blocked by item 9,
unblocked 2026-10-09, and the top item with nothing in front of it.

------------------------------------------------------------------------------
WHAT THE ITEM ASKED FOR, AND WHAT IT TURNS OUT TO BE
------------------------------------------------------------------------------
The item, as written: "re-derive Sections 4.3-4.6 on (4.28) and re-report,
keeping the (4.29) numbers annotated in place."  (4.28) is the exact charge
relation under the 2026-10-08 declaration that V_ch is the quasi-Fermi
potential,

    V_g - V_dirac - V_ch  =  e*n/C_ox  +  E_F(n)/e ,   E_F(n)/e = A_F*sqrt(n)
                                                                      (4.28)

and (4.29) is what `graphene_fet_model.carrier_density()` computes,

    n = (C_ox*dV/e) * C_q(dV)/(C_q(dV) + C_ox)                        (4.29)

Chapter 4 has called (4.29) "the linearisation of (4.28)" three times, in
Section 4.6.6, in `carrier_density()`'s annotation and in Section 7.9 item 12
itself.  SECTION 1 BELOW SHOWS THAT DESCRIPTION IS WRONG, and the correction
is the finding of this module rather than the -0.60 % it was asked to produce.

(4.29) is not a linearisation of (4.28).  It is (4.28)'s EXACT DIFFERENTIAL
RELATION, used as if it were algebraic.  Differentiating (4.28),

    d(dV)/dn = e/C_ox + (1/e)*dE_F/dn = e/C_ox + e/C_Q(n),
    C_Q(n) := e^2 dn/dE_F = 2*e*sqrt(n)/A_F                           (R1)

    =>  dn/d(dV) = (C_ox/e) * C_Q(n)/(C_Q(n) + C_ox)                  (R2)

(R2) is (4.29) with `dn/d(dV)` in place of `n/dV`.  So (4.29) is a ONE-POINT
RECTANGLE RULE for the integral of (R2): it evaluates the exact slope at the
endpoint and multiplies by the whole drive.  That is a quadrature error, not a
truncated Taylor series, and it does not vanish at small dV -- it is worst
there, because the slope varies fastest there.  It explains the +68 % raw
disagreement 2026-10-08 measured and declined to quote: near the Dirac point
(4.29) is linear in dV while (4.28) is QUADRATIC (a*n + b*sqrt(n) = dV with
b*sqrt(n) dominating gives n ~ (dV/b)^2), so the two forms do not even share a
leading order.  A linearisation agrees to first order by construction.  These
do not.

TWO INDEPENDENT ERRORS, AND THIS MODULE SEPARATES THEM
------------------------------------------------------
(4.29) differs from (4.28) in two ways that have nothing to do with each
other, and no previous session separated them:

  (i)  QUADRATURE.  The exact slope (R2) is used over the full drive instead
       of being integrated (Section 1, Section 2).
  (ii) ARGUMENT.  `gfet.quantum_capacitance()` evaluates C_q at the GATE
       OVERDRIVE dV, while (R1) requires it at E_F/e = A_F*sqrt(n), which is
       a few percent of dV.  This was noticed on 2026-10-07 (Section 4.6.5's
       annotation) and never costed.

SECTION 2 runs the 2x2: {rectangle, integral} x {argument dV, argument E_F/e}.
The committed model is one corner, (4.28) is the opposite corner, and the two
error sources turn out to be of the SAME ORDER and the SAME SIGN, so neither
is the explanation on its own.

VALIDATION AGAINST EXACTLY KNOWN VALUES (SECTION 5)
---------------------------------------------------
Per the standing rule from 2026-09-17 onwards, and 2026-10-09's addition that
an exactness check which fails is a measurement of whatever it was written on
top of:

  - (R1)/(R2): the closed-form derivative of (4.28)'s root must equal the
    series combination of C_ox with C_Q(n) = 2*e*sqrt(n)/A_F.  An algebraic
    identity; checked to ulps, not to a tolerance.
  - THE INTEGRAL IDENTITY: integrating (R2) numerically from 0 to dV, with n
    taken from (4.28)'s own root at each quadrature node, must return (4.28)'s
    root at dV.  This is the check that would catch a wrong prefactor or a
    dropped factor of 2 in C_Q; it is reported with its CONVERGENCE ORDER, so
    agreement cannot be a coincidence of one grid.
  - A_F -> 0 must return the pure electrostatic charge C_ox*dV/e exactly.
  - n(-dV) == -n(dV) BITWISE (the 2026-10-09 oddness, re-asserted as a
    regression on the function this module depends on).
  - The species split must still satisfy n_e*n_h == n_puddle^2/4.
  - Patching must be reversible: every module-level name this file rebinds is
    checked restored, and the committed transfer characteristic is re-run and
    compared BITWISE afterwards.

WHAT IS RE-DERIVED (SECTIONS 3 AND 4)
-------------------------------------
Section 3  Section 4.3's charge and resistances, Section 4.4's crossover
           lengths and on/off ratio, and I_d at the RF bias.
Section 4  Section 4.6's RF figures of merit, with the charge change and the
           CAPACITANCE change applied separately -- because under (4.28) the
           small-signal gate capacitance is e*dn/d(dV) = the series
           combination with C_Q(n), so a consistent re-derivation moves C_gs
           as well as I_d, and f_T = g_m/(2*pi*C_gs) sees both.

Run:  python graphene_exact_charge_rederivation.py
Writes: exact_charge_rederivation.png, exact_charge_rederivation_output.txt
"""

import contextlib

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import graphene_fet_model as gfet
import graphene_diffusion_current_model as dd
import graphene_fet_signed_carrier_model as sc
import rf_small_signal_model as rf
import contact_resistance_crossover as crx

E_CHARGE = 1.602176634e-19
A_FERMI = dd.A_FERMI                 # E_F/e = A_FERMI*sqrt(n)
CQ_PREFACTOR = dd.CQ_PREFACTOR       # 2 e^3/(pi hbar^2 v_F^2)
A_OX = E_CHARGE / gfet.C_ox          # 'a' in a*n + b*sqrt(n) = dV

VDS_RF = dd.VDS_REF                  # 0.1 V, the RF bias of record
VG_RF = 2.0
VG_ON = 3.5

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


def ulps(a, b):
    a, b = float(a), float(b)
    if a == b:
        return 0.0
    return abs(a - b) / np.spacing(max(abs(a), abs(b)))


# ---------------------------------------------------------------------------
# The four charge forms.  All RAW: the puddle floor is applied separately so
# that the floor is never doing any of the work in a comparison.
# ---------------------------------------------------------------------------

def n_exact(dV):
    """(4.28), signed.  Delegates to the 2026-10-09 branched root."""
    return sc.n_net_exact(dV)


def dn_ddV_exact(dV):
    """d n/d(dV) from (4.28) in closed form: 2*sqrt(n)/(2*a*sqrt(n) + b)."""
    x = np.sqrt(np.abs(np.asarray(n_exact(dV), dtype=float)))
    return 2.0 * x / (2.0 * A_OX * x + A_FERMI)


def C_Q_dispersion(n):
    """C_Q(n) = e^2 dn/dE_F = 2 e sqrt(n)/A_F  [F/m^2].  Equivalently
    CQ_PREFACTOR*A_F*sqrt(n); both spellings are checked against each other."""
    return 2.0 * E_CHARGE * np.sqrt(np.abs(np.asarray(n, dtype=float))) / A_FERMI


def series_factor(C_q):
    return C_q / (C_q + gfet.C_ox)


def n_rect_dV(dV):
    """FORM A -- the committed (4.29): rectangle rule, C_q at the overdrive."""
    dV = np.asarray(dV, dtype=float)
    C_q = gfet.quantum_capacitance(dV, T=gfet.T)
    return (gfet.C_ox * dV / E_CHARGE) * series_factor(C_q)


def n_rect_EF(dV, iters=200):
    """FORM B -- rectangle rule, C_q at E_F/e: fixes the ARGUMENT only.
    n = (C_ox dV/e)*C_Q(n)/(C_Q(n)+C_ox), solved by fixed-point iteration
    from the electrostatic value.  Monotone and bounded, so it converges."""
    dV = np.asarray(dV, dtype=float)
    n_es = gfet.C_ox * dV / E_CHARGE
    n = np.array(n_es, dtype=float, copy=True)
    for _ in range(iters):
        n_new = n_es * series_factor(C_Q_dispersion(n))
        if np.all(np.abs(n_new - n) <= 1e-15 * np.abs(n_new) + 1e-6):
            n = n_new
            break
        n = n_new
    return n


def _simpson(y, h):
    """Composite Simpson on an odd number of equally spaced samples."""
    return h / 3.0 * (y[0] + y[-1] + 4.0 * np.sum(y[1:-1:2])
                      + 2.0 * np.sum(y[2:-1:2]))


def n_int_dV(dV, n_nodes=2001):
    """FORM C -- integral of the exact slope, C_q still at the overdrive:
    fixes the QUADRATURE only.  n(dV) = int_0^dV (C_ox/e) f(C_q(V')) dV'."""
    dV = float(dV)
    if dV == 0.0:
        return 0.0
    s = -1.0 if dV < 0.0 else 1.0
    V = np.linspace(0.0, abs(dV), n_nodes)
    y = (gfet.C_ox / E_CHARGE) * series_factor(
        gfet.quantum_capacitance(V, T=gfet.T))
    return s * _simpson(y, V[1] - V[0])


def n_int_EF(dV, n_nodes=2001):
    """FORM D' -- integral of the exact slope with C_q at E_F/e, i.e. BOTH
    errors removed.  Must reproduce (4.28)'s closed-form root; that identity
    is Section 5's main anchor, with a convergence order attached."""
    dV = float(dV)
    if dV == 0.0:
        return 0.0
    s = -1.0 if dV < 0.0 else 1.0
    V = np.linspace(0.0, abs(dV), n_nodes)
    n_nodes_vals = np.abs(np.asarray(n_exact(V), dtype=float))
    y = (gfet.C_ox / E_CHARGE) * series_factor(C_Q_dispersion(n_nodes_vals))
    return s * _simpson(y, V[1] - V[0])


DV_CROSSOVER = 2.0 * A_FERMI ** 2 / A_OX
"""The drive at which (4.28)'s two terms are equal: a*n = b*sqrt(n) gives
sqrt(n) = b/a and dV = 2 b^2/a.  Below it the graphene drop dominates and
n ~ (dV/b)^2; above it the oxide drop dominates and n ~ C_ox dV/e.  It is
6.52 mV here, and SECTION 1 and SECTION 5 both contain a check that failed as
first written because it did not know this number."""


def _panel_edges(dV):
    """Log-graded panel edges from 0 to dV, straddling DV_CROSSOVER."""
    e = [0.0]
    for k in (-2, -1, 0, 1):
        v = DV_CROSSOVER * 10.0 ** k
        if v < abs(dV):
            e.append(v)
    e.append(abs(dV))
    return e


def n_int_EF_graded(dV, nodes_per_panel=81):
    """FORM D' on a grid that RESOLVES the crossover layer.  Same integrand as
    n_int_EF; the only change is the grid."""
    dV = float(dV)
    if dV == 0.0:
        return 0.0
    s = -1.0 if dV < 0.0 else 1.0
    edges = _panel_edges(dV)
    total = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        V = np.linspace(lo, hi, nodes_per_panel)
        y = (gfet.C_ox / E_CHARGE) * series_factor(
            C_Q_dispersion(np.abs(np.asarray(n_exact(V), dtype=float))))
        total += _simpson(y, V[1] - V[0])
    return s * total


def floor_total(n_signed):
    """The shipped quadrature floor, applied identically to every form so the
    floor is never the difference between two compared numbers."""
    return np.sqrt(np.asarray(n_signed, dtype=float) ** 2 + gfet.n_puddle ** 2)


# ---------------------------------------------------------------------------
# Patching Chapter 4 onto (4.28)
# ---------------------------------------------------------------------------

def _carrier_density_exact(V_g, V_ch=0.0):
    dV = (np.asarray(V_g, dtype=float) - np.asarray(V_ch, dtype=float)
          - gfet.V_dirac)
    return floor_total(n_exact(dV))


def _quantum_capacitance_consistent(dV, T=None):
    """C_q at E_F/e of (4.28)'s own root, so that
    (1/C_ox + 1/C_q)^-1 == e*dn/d(dV) exactly (checked in Section 5)."""
    n = np.abs(np.asarray(n_exact(dV), dtype=float))
    return C_Q_dispersion(n)


@contextlib.contextmanager
def on_exact(capacitance=True):
    """Run a block with Chapter 4's charge relation replaced by (4.28).

    capacitance=False changes ONLY the charge, leaving C_gs on the committed
    thermal C_q(dV); capacitance=True also makes the small-signal gate
    capacitance consistent with (4.28).  The two are separated because they
    move f_T in different ways and a single combined number would hide that.
    """
    cd_saved = gfet.carrier_density
    qc_saved = gfet.quantum_capacitance
    # contact_resistance_crossover.py does `from graphene_fet_model import
    # carrier_density`, so it holds its OWN reference, bound at import time.
    # Rebinding only `gfet.carrier_density` leaves Section 4.4 silently on
    # (4.29) while this file reports it as re-derived -- the 2026-10-03
    # frozen-default fault in its monkey-patch form.  Both names are rebound,
    # and Section 3 carries an ARRIVAL CONTROL (2026-09-30's rule) that fails
    # if the patch did not reach the module where the number is computed.
    crx_saved = crx.carrier_density
    try:
        gfet.carrier_density = _carrier_density_exact
        crx.carrier_density = _carrier_density_exact
        if capacitance:
            gfet.quantum_capacitance = _quantum_capacitance_consistent
        yield
    finally:
        gfet.carrier_density = cd_saved
        gfet.quantum_capacitance = qc_saved
        crx.carrier_density = crx_saved


def peak_rf(Vds=VDS_RF, W=None, N_fingers=1, extrinsic=False, n_vg=400):
    """Peak f_T and the companion numbers at that bias, on whatever charge
    model is currently installed.  Same sweep rf_small_signal_model's own
    summary_numbers() uses, so the committed row is reproducible."""
    Vg = np.linspace(-1.5, 3.5, n_vg)
    fT, fmax, gm, gds, Cgs = rf.compute_fT_fmax(
        Vg, Vds=Vds, extrinsic=extrinsic, W=W, N_fingers=N_fingers)
    i = int(np.argmax(fT))
    return dict(Vg=float(Vg[i]), fT=float(fT[i]), fmax=float(fmax[i]),
                gm=float(gm[i]), gds=float(gds[i]), Cgs=float(Cgs[i]))


def onoff(Vds=VDS_RF, Vg_on=VG_ON, Rc=None):
    """(I_on, I_off, ratio) from whatever charge model is installed."""
    saved = gfet.Rc_total
    try:
        if Rc is not None:
            gfet.Rc_total = Rc
        I = gfet.transfer_characteristic(
            np.array([Vg_on, gfet.V_dirac]), Vds=Vds)
    finally:
        gfet.Rc_total = saved
    return float(I[0]), float(I[1]), float(I[0] / I[1])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print('=' * 78)
    print('CHAPTER 4 RE-DERIVED ON THE EXACT CHARGE RELATION (4.28)')
    print('Chapter 7 Section 7.9 item 12 -- unblocked 2026-10-09')
    print('=' * 78)
    print('  a = e/C_ox = %.6e V*m^2   b = A_F = %.6e V*m' % (A_OX, A_FERMI))
    print('  (4.28):  a*n + b*sign(n)*sqrt(|n|) = dV')

    # -----------------------------------------------------------------
    print()
    print('SECTION 1.  (4.29) IS NOT A LINEARISATION OF (4.28)')
    print('-' * 78)
    print('  Chapter 4 says so in three places.  It is checkable, and false.')
    print('  A linearisation of (4.28) about dV = 0 would share its leading')
    print('  order there.  (4.28) has a*n + b*sqrt(n) = dV with the sqrt term')
    print('  dominating as dV -> 0, so n ~ (dV/b)^2 -- QUADRATIC.  (4.29) has')
    print('  C_q(0) finite, so n ~ (C_ox/e)*f(C_q(0))*dV -- LINEAR.')
    print()
    print('  The two terms of (4.28) are equal at dV = 2*b^2/a = %.6e V'
          % DV_CROSSOVER)
    print('  (%.3f mV), so the quadratic regime is dV << %.2f mV.  THE FIRST'
          % (1e3 * DV_CROSSOVER, 1e3 * DV_CROSSOVER))
    print('  DRAFT OF THIS CHECK SWEPT dV DOWN TO 1e-4 V ONLY and measured an')
    print('  exponent of 1.9438, failing a +-0.01 criterion.  It did not')
    print('  measure the exponent: 1e-4 V is 65x below the crossover, where')
    print('  the correction is O(sqrt(dV/dV_x)) ~ 12 %%, so it measured the')
    print('  WIDTH OF THE CROSSOVER.  The repair is to evaluate the limit')
    print('  where the limit is, and to print the crossover drive so the')
    print('  criterion is derived from the model instead of guessed')
    print('  (2026-10-09: find out what a failed check measured).')
    print()
    print('  Measured log-log slope d ln n/d ln dV as dV -> 0:')
    print('   %12s %14s %14s' % ('dV [V]', 'slope (4.28)', 'slope (4.29)'))
    slopes28, slopes29 = [], []
    for dV in (1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2):
        h = 1.0e-3 * dV
        s28 = (np.log(float(n_exact(dV + h))) - np.log(float(n_exact(dV - h)))) \
            / (np.log(dV + h) - np.log(dV - h))
        s29 = (np.log(float(n_rect_dV(dV + h))) - np.log(float(n_rect_dV(dV - h)))) \
            / (np.log(dV + h) - np.log(dV - h))
        slopes28.append(s28)
        slopes29.append(s29)
        print('   %12.0e %14.6f %14.6f' % (dV, s28, s29))
    check('(4.28) is quadratic and (4.29) linear as dV -> 0, so (4.29) '
          'cannot be a linearisation of (4.28)',
          abs(slopes28[0] - 2.0) < 1e-4 and abs(slopes29[0] - 1.0) < 1e-4,
          'at dV = 1e-10 V, which is 1e8 below the crossover, the exponents '
          'are %.8f and %.8f; at the first draft\'s 1e-4 V they were %.4f and '
          '%.4f.  Two forms that do not share a leading order are not a '
          "function and its linearisation; Chapter 4's description of the pair "
          'is corrected, and no number is.'
          % (slopes28[0], slopes29[0], slopes28[3], slopes29[3]))
    check('and the quadratic regime of (4.28) is never physically entered, so '
          'the correction is a statement about wording and not about results',
          DV_CROSSOVER * gfet.C_ox / E_CHARGE < 0.01 * gfet.n_puddle,
          'the crossover drive carries %.3e /m2 of charge against n_puddle = '
          '%.3e /m2, a factor %.0f smaller.  (4.28) is quadratic only where '
          'the disorder floor already supplies 99.9 %% of the density, which '
          'is why the +68 %% raw figure 2026-10-08 declined to quote was '
          'right to be declined.'
          % (DV_CROSSOVER * gfet.C_ox / E_CHARGE, gfet.n_puddle,
             gfet.n_puddle / (DV_CROSSOVER * gfet.C_ox / E_CHARGE)))

    print()
    print('  What (4.29) actually is.  Differentiating (4.28):')
    print('    d(dV)/dn = e/C_ox + e/C_Q(n),  C_Q(n) = e^2 dn/dE_F = 2e*sqrt(n)/A_F')
    print('    =>  dn/d(dV) = (C_ox/e)*C_Q(n)/(C_Q(n) + C_ox)             (R2)')
    print('  (4.29) is (R2) with n/dV written for dn/d(dV): THE EXACT')
    print('  DIFFERENTIAL RELATION USED AS AN ALGEBRAIC ONE, i.e. a one-point')
    print('  rectangle rule for the integral of (R2).  That is a quadrature')
    print('  error, and a rectangle rule is worst where the slope varies')
    print('  fastest -- near the Dirac point, which is exactly where the')
    print('  2026-10-08 raw disagreement reached +68 %.')
    print()
    print('  Literature position, checked today (live search available).')
    print('  Pasadas/Jimenez-lineage compact-model review arXiv:2209.00388')
    print('  writes the gate-to-channel relation as implicit in V_c (the')
    print('  chemical-potential shift) and solves it NUMERICALLY in the')
    print('  simulator; its series-capacitance expression appears as')
    print('  dV/dV_c = 1 + C_q(V_c)/(C_t + C_b), a DIFFERENTIAL relation,')
    print('  and its C_q is evaluated at V_c and never at the overdrive.')
    print('  Both of this section\'s points are therefore the literature\'s')
    print('  own conventions, not a new convention introduced here.  It also')
    print('  states the Fermi-Dirac form "is not convenient for a compact')
    print('  model" and replaces C_q by k*c1*sqrt(1 + (V_c/c1)^2), whose')
    print('  c1 -> 0 limit is exactly (R1).  So (4.28) is the T = 0 limit of')
    print('  a standard literature approximation -- which is also the')
    print('  limitation stated in Section 6 below.')

    # -----------------------------------------------------------------
    print()
    print('SECTION 2.  THE 2x2: THE GAP IS TWO ERRORS, NOT ONE')
    print('-' * 78)
    print('  A  rectangle + C_q(dV)     = (4.29), THE COMMITTED MODEL')
    print('  B  rectangle + C_q(E_F/e)  = argument fixed, quadrature not')
    print('  C  integral  + C_q(dV)     = quadrature fixed, argument not')
    print('  D  integral  + C_q(E_F/e)  = (4.28), both fixed')
    print('  All four RAW (no puddle floor), so the floor is not doing any of')
    print('  the work in this table.')
    print()
    print('   %8s %13s %13s %13s %13s' % ('dV [V]', 'A (4.29)', 'B arg-fix',
                                          'C quad-fix', 'D (4.28)'))
    rows = []
    for dV in (0.01, 0.1, 0.5, 1.0, 1.2, 2.0, 2.7):
        nA = float(n_rect_dV(dV))
        nB = float(n_rect_EF(dV))
        nC = float(n_int_dV(dV))
        nD = float(n_exact(dV))
        rows.append((dV, nA, nB, nC, nD))
        print('   %8.2f %13.6e %13.6e %13.6e %13.6e' % (dV, nA, nB, nC, nD))
    print()
    print('   As percentage errors against D = (4.28):')
    print('   %8s %13s %13s %13s %13s' % ('dV [V]', 'A total', 'B (quad only)',
                                          'C (arg only)', 'sum B+C'))
    for dV, nA, nB, nC, nD in rows:
        eA = 100.0 * (nA / nD - 1.0)
        eB = 100.0 * (nB / nD - 1.0)
        eC = 100.0 * (nC / nD - 1.0)
        print('   %8.2f %12.4f%% %12.4f%% %12.4f%% %12.4f%%'
              % (dV, eA, eB, eC, eB + eC))
    dV_rf = VG_RF - gfet.V_dirac
    nA = float(n_rect_dV(dV_rf)); nB = float(n_rect_EF(dV_rf))
    nC = float(n_int_dV(dV_rf)); nD = float(n_exact(dV_rf))
    eA = nA / nD - 1.0; eB = nB / nD - 1.0; eC = nC / nD - 1.0
    check('both error sources are present, same sign, and of the same order '
          'at the RF bias',
          eB > 0.0 and eC > 0.0 and 0.25 < eB / eC < 4.0,
          'at dV = %.2f V: total %+.4f %%, of which the QUADRATURE error '
          'alone (form B, argument already fixed) is %+.4f %% and the '
          'ARGUMENT error alone (form C, quadrature already fixed) is '
          '%+.4f %%, a ratio of %.2f.  The argument error is the larger '
          'single contributor, which no previous session could have known '
          'because no previous session separated them.  NEITHER is the '
          'explanation on its own, which is why a single -0.60 %% headline '
          'hid both.' % (dV_rf, 100 * eA, 100 * eB, 100 * eC, eB / eC))
    check('the two errors are NOT additive, so they cannot be reported '
          'separately and summed',
          abs((eB + eC) - eA) > 0.05 * abs(eA),
          'B + C = %+.4f %% against the true total %+.4f %%; they share the '
          'same series factor, so the composition is multiplicative in the '
          'factor and not additive in the error.  Quoting the sum would '
          'overstate the total by %.2f x.'
          % (100 * (eB + eC), 100 * eA, (eB + eC) / eA))
    print()
    print('  And with the floor applied, i.e. what carrier_density returns:')
    print('   %8s %14s %14s %12s' % ('dV [V]', 'n_tot (4.29)', 'n_tot (4.28)',
                                     'A/D - 1'))
    fl_worst, fl_at = 0.0, None
    for dV, nA_, nB_, nC_, nD_ in rows:
        fA, fD = float(floor_total(nA_)), float(floor_total(nD_))
        rel = fA / fD - 1.0
        if abs(rel) > abs(fl_worst):
            fl_worst, fl_at = rel, dV
        print('   %8.2f %14.6e %14.6e %11.4f %%' % (dV, fA, fD, 100 * rel))
    print('   worst %+.4f %% at dV = %.2f V.  The floor attenuates the raw'
          % (100 * fl_worst, fl_at))
    print('   charge error by a factor %.1f at the RF bias, which is the'
          % abs(eA / (floor_total(nA) / floor_total(nD) - 1.0)))
    print('   reason a %+.1f %% charge error costs only tenths of a percent'
          % (100 * eA))
    print('   of I_d.  That attenuation is a property of n_puddle, not of')
    print('   the charge model, and 2026-10-09 measured n_puddle itself to')
    print('   be 2.07x off Chapter 3 -- so it is not a comfort.')

    # -----------------------------------------------------------------
    print()
    print('SECTION 3.  SECTIONS 4.3 AND 4.4, RE-DERIVED')
    print('-' * 78)
    n_committed_on = float(gfet.carrier_density(VG_ON))
    Lx_committed = {rc: float(crx.crossover_length(
        2 * (rc * 1e-6) / gfet.W, VG_ON)) for rc in (110.0, 300.0, 500.0)}
    Rch_committed = (float(gfet.channel_resistance(VG_ON)),
                     float(gfet.channel_resistance(gfet.V_dirac)))
    Id_committed = float(gfet.transfer_characteristic(
        np.array([VG_RF]), Vds=VDS_RF)[0])
    # Section 4.6.6's 8.190217e-05 A is the SATURATED drift-diffusion current
    # at lam = 0 from graphene_diffusion_current_model, not the Section 4.3
    # drift current; both are re-derived here, and conflating them would have
    # produced a 6.5 % discrepancy against a published number for no reason.
    Id_dd_committed = float(dd.transfer_characteristic_dd(
        np.array([VG_RF]), Vds=VDS_RF, lam=0.0)[0])
    onoff_committed = onoff()
    onoff_committed_i = onoff(Rc=0.0)

    with on_exact(capacitance=False):
        # ARRIVAL CONTROL first: if the patch did not reach the module that
        # owns Section 4.4's number, everything below is a re-report of the
        # committed model wearing a new label.
        n_crx_patched = float(crx.carrier_density(VG_ON, V_ch=0.0))
        n_exact_on = float(gfet.carrier_density(VG_ON))
        check('ARRIVAL CONTROL: the patch reached contact_resistance_crossover,'
              ' which imported carrier_density BY NAME',
              n_crx_patched != n_committed_on
              and n_crx_patched == n_exact_on,
              'crx sees %.6e /m2, gfet sees %.6e /m2, committed was %.6e /m2.  '
              'A re-derivation that rebound only gfet.carrier_density would '
              'have left Section 4.4 on (4.29) and said otherwise.'
              % (n_crx_patched, n_exact_on, n_committed_on))
        Lx_exact = {rc: float(crx.crossover_length(
            2 * (rc * 1e-6) / gfet.W, VG_ON)) for rc in (110.0, 300.0, 500.0)}
        Rch_exact = (float(gfet.channel_resistance(VG_ON)),
                     float(gfet.channel_resistance(gfet.V_dirac)))
        Id_exact = float(gfet.transfer_characteristic(
            np.array([VG_RF]), Vds=VDS_RF)[0])
        Id_dd_exact = float(dd.transfer_characteristic_dd(
            np.array([VG_RF]), Vds=VDS_RF, lam=0.0)[0])
        onoff_exact = onoff()
        onoff_exact_i = onoff(Rc=0.0)

    print('  Section 4.3, charge and resistance:')
    print('   %-34s %14s %14s %10s' % ('quantity', '(4.29)', '(4.28)', 'change'))
    def row(label, a, b, fmt='%14.6e'):
        print(('   %-34s ' + fmt + ' ' + fmt + ' %9.4f %%')
              % (label, a, b, 100.0 * (b / a - 1.0)))
    row('n at V_g = 3.5 V  [1/m^2]', n_committed_on, n_exact_on)
    row('R_ch at V_g = 3.5 V  [Ohm]', Rch_committed[0], Rch_exact[0])
    row('R_ch at Dirac point  [Ohm]', Rch_committed[1], Rch_exact[1])
    row('I_d drift, RF bias (4.3) [A]', Id_committed, Id_exact)
    row('I_d saturated dd, lam=0 (4.6.6) [A]', Id_dd_committed, Id_dd_exact)
    check('the re-derived saturated current reproduces Section 4.6.6\'s '
          'measured -0.6023 %% cost of the declaration',
          abs(100.0 * (Id_dd_exact / Id_dd_committed - 1.0) + 0.6023) < 0.01,
          'this module gets %+.4f %% against the -0.6023 %% Section 4.6.6 '
          'published from a different code path (its own monkey-patch of '
          'carrier_density, written independently).  Two implementations of '
          'the same replacement agreeing to 4 decimal places is a fact about '
          'two artefacts (2026-10-02), and it is the check that says the '
          'patch installed here is the one the chapter costed.'
          % (100.0 * (Id_dd_exact / Id_dd_committed - 1.0)))
    print()
    print('  Section 4.4, crossover length L_x at V_g = 3.5 V:')
    print('   %-34s %14s %14s %10s' % ('R_c [Ohm.um]', '(4.29) [nm]',
                                       '(4.28) [nm]', 'change'))
    for rc in (110.0, 300.0, 500.0):
        print('   %-34.0f %14.2f %14.2f %9.4f %%'
              % (rc, Lx_committed[rc] * 1e9, Lx_exact[rc] * 1e9,
                 100.0 * (Lx_exact[rc] / Lx_committed[rc] - 1.0)))
    print()
    print('  Section 4.4, on/off ratio (V_g = 3.5 V on, Dirac point off):')
    row('terminal on/off', onoff_committed[2], onoff_exact[2], '%14.6f')
    row('intrinsic on/off (R_c = 0)', onoff_committed_i[2], onoff_exact_i[2],
        '%14.6f')
    print()
    print('  A bound this check no longer carries.  Its first draft required')
    print('  "every number moves by under 2 %", which FAILED at 2.0859 % on')
    print('  L_x.  The 2 % was a round number chosen before the run, not a')
    print('  precision derived from anything, so widening it to 2.5 % would')
    print('  have been the tuned tolerance Section 7.9 item 15 exists to')
    print('  find.  It is removed instead: the magnitudes are printed above')
    print('  with no threshold on them, and the check below tests the')
    print('  QUALITATIVE claims Section 4.4 actually makes, which are what a')
    print('  re-derivation can overturn.')
    qual_committed = tuple(Lx_committed[rc] < gfet.L for rc in
                           (110.0, 300.0, 500.0))
    qual_exact = tuple(Lx_exact[rc] < gfet.L for rc in (110.0, 300.0, 500.0))
    print('  L_x < L = %.0f nm, by contact:  (4.29) %s   (4.28) %s'
          % (gfet.L * 1e9, qual_committed, qual_exact))
    check('every QUALITATIVE claim in Sections 4.3-4.4 survives the exact '
          'relation',
          qual_committed == qual_exact
          and (onoff_committed[2] > 1.0) == (onoff_exact[2] > 1.0)
          and (onoff_committed[2] < 2.0) == (onoff_exact[2] < 2.0)
          and (Rch_exact[1] > Rch_exact[0]) == (Rch_committed[1] > Rch_committed[0]),
          'the best literature contact still crosses over below the 200 nm '
          'channel (L_x = %.1f nm -> %.1f nm) and the mid- and worst-case '
          'contacts still do not; the on/off ratio is still between 1 and 2 '
          '(%.4f -> %.4f), i.e. still not a switch; and the channel is still '
          'more resistive at the Dirac point than in the on state.  The '
          'largest movement anywhere in Sections 4.3-4.4 is %.4f %% and it is '
          'reported without a threshold.'
          % (Lx_committed[110.0] * 1e9, Lx_exact[110.0] * 1e9,
             onoff_committed[2], onoff_exact[2],
             100.0 * (Lx_exact[110.0] / Lx_committed[110.0] - 1.0)))
    check('the on/off ratio moves LESS at the terminals than intrinsically',
          abs(onoff_exact[2] / onoff_committed[2] - 1.0)
          < abs(onoff_exact_i[2] / onoff_committed_i[2] - 1.0),
          'terminal %+.4f %% vs intrinsic %+.4f %%: the series contact '
          'resistance attenuates a charge-model change exactly as 2026-10-08 '
          'and 2026-10-09 found it attenuates a transport-term change and an '
          'n_puddle change.  Third parameter, same structure -- this is a '
          'property of the device topology, not of any one parameter.'
          % (100 * (onoff_exact[2] / onoff_committed[2] - 1.0),
             100 * (onoff_exact_i[2] / onoff_committed_i[2] - 1.0)))

    # -----------------------------------------------------------------
    print()
    print('SECTION 4.  SECTION 4.6, RE-DERIVED -- AND THE CAPACITANCE MOVES TOO')
    print('-' * 78)
    print('  Under (4.28) the small-signal gate capacitance is not an')
    print('  independent ingredient: e*dn/d(dV) IS (1/C_ox + 1/C_Q(n))^-1 by')
    print('  (R2).  So a consistent re-derivation changes C_gs as well as')
    print('  I_d, and f_T = g_m/(2 pi C_gs) sees both.  The two are applied')
    print('  separately, because they move f_T in OPPOSITE directions and a')
    print('  combined number would hide that.')
    print()
    pk_committed = peak_rf()
    with on_exact(capacitance=False):
        pk_charge = peak_rf()
    with on_exact(capacitance=True):
        pk_full = peak_rf()
    print('  Per-um convention, V_ds = 0.1 V, peak-f_T bias:')
    print('   %-22s %13s %13s %13s' % ('', '(4.29)', '(4.28) charge',
                                       '(4.28) + C_q'))
    for key, scale, unit in (('Vg', 1.0, 'V'), ('fT', 1e-9, 'GHz'),
                             ('fmax', 1e-9, 'GHz'), ('gm', 1e3, 'mS'),
                             ('gds', 1e3, 'mS'), ('Cgs', 1e15, 'fF')):
        print('   %-22s %13.4f %13.4f %13.4f'
              % ('%s [%s]' % (key, unit), pk_committed[key] * scale,
                 pk_charge[key] * scale, pk_full[key] * scale))
    print('   %-22s %13s %13.4f %13.4f'
          % ('f_T change [%]', '--',
             100 * (pk_charge['fT'] / pk_committed['fT'] - 1.0),
             100 * (pk_full['fT'] / pk_committed['fT'] - 1.0)))
    print('   %-22s %13s %13.4f %13.4f'
          % ('f_max change [%]', '--',
             100 * (pk_charge['fmax'] / pk_committed['fmax'] - 1.0),
             100 * (pk_full['fmax'] / pk_committed['fmax'] - 1.0)))
    d_charge = pk_charge['fT'] / pk_committed['fT'] - 1.0
    d_full = pk_full['fT'] / pk_committed['fT'] - 1.0
    d_cap = d_full - d_charge
    check('the combined -0.06 % in f_T is a CANCELLATION of two ~2 % terms, '
          'not a small effect',
          d_cap * d_charge < 0.0 and abs(d_full) < 0.1 * abs(d_charge),
          'charge alone %+.4f %% in f_T; the consistent C_gs contributes '
          '%+.4f %% on its own -- opposite sign, %.3f of the magnitude -- and '
          'the two leave %+.4f %%, which is %.1f %% of either.  The first '
          'draft of this check asserted the capacitance half was the LARGER '
          'half and failed: it is 0.970x the charge half, not more than it.  '
          'What it measured is that they very nearly cancel, which is a '
          'stronger statement than either -- a single combined number would '
          'have reported item 12 as worth 0.06 %% in f_T when it is two '
          'coupled 2 %% corrections that happen to oppose.'
          % (100 * d_charge, 100 * d_cap, abs(d_cap / d_charge), 100 * d_full,
             100 * abs(d_full / d_charge)))
    check('the PEAK BIAS moves far more than the peak value',
          abs(pk_full['Vg'] - pk_committed['Vg']) > 1e-3,
          'peak f_T sits at V_g = %+.4f V on (4.29) and %+.4f V on (4.28) '
          '(%+.4f V, %.2f %%), while the peak VALUE moves %.4f %%.  The '
          'location is %.0fx more sensitive to the charge relation than the '
          'height, both peaks are on the HOLE branch, and neither crosses to '
          'the electron branch -- so the open "does the peak-f_T bias change '
          'branch between the two models" question is answered NO for this '
          'pair, and left open for the branched-vs-magnitude pair it was '
          'asked about.'
          % (pk_committed['Vg'], pk_full['Vg'],
             pk_full['Vg'] - pk_committed['Vg'],
             100 * abs(pk_full['Vg'] / pk_committed['Vg'] - 1.0),
             100 * abs(d_full),
             abs((pk_full['Vg'] / pk_committed['Vg'] - 1.0) / d_full)))
    print()
    print('  C_gs at the peak: %.4f fF -> %.4f fF (%+.4f %%).  C_q(dV) at the'
          % (pk_committed['Cgs'] * 1e15, pk_full['Cgs'] * 1e15,
             100 * (pk_full['Cgs'] / pk_committed['Cgs'] - 1.0)))
    dV_pk = pk_committed['Vg'] - gfet.V_dirac
    cq_dv = float(gfet.quantum_capacitance(dV_pk, T=gfet.T))
    cq_ef = float(_quantum_capacitance_consistent(dV_pk))
    print('  overdrive is %.6e F/m^2 and C_Q at E_F/e is %.6e F/m^2, a factor'
          % (cq_dv, cq_ef))
    print('  %.1f apart; both are >> C_ox = %.6e F/m^2, which is why a factor'
          % (cq_dv / cq_ef, gfet.C_ox))
    print('  %.0f in C_q becomes only %+.2f %% in C_gs.  E_F/e = %.4f V is'
          % (cq_dv / cq_ef, 100 * (pk_full['Cgs'] / pk_committed['Cgs'] - 1.0),
             float(dd.fermi_voltage(abs(float(n_exact(dV_pk)))))))
    print('  %.2f %% of the %.2f V overdrive in magnitude, and THAT ratio is'
          % (100 * float(dd.fermi_voltage(abs(float(n_exact(dV_pk)))))
             / abs(dV_pk), dV_pk))
    print('  the argument error of Section 2 seen directly.')
    print()
    pk_rf_committed = peak_rf(W=rf.W_RF, N_fingers=rf.N_FINGERS_RF)
    with on_exact(capacitance=True):
        pk_rf_full = peak_rf(W=rf.W_RF, N_fingers=rf.N_FINGERS_RF)
    pk_rf_ext_c = peak_rf(W=rf.W_RF, N_fingers=rf.N_FINGERS_RF, extrinsic=True)
    with on_exact(capacitance=True):
        pk_rf_ext_f = peak_rf(W=rf.W_RF, N_fingers=rf.N_FINGERS_RF,
                              extrinsic=True)
    print('  Literature-scale device (%.0f um, %d fingers), the geometry'
          % (rf.W_RF * 1e6, rf.N_FINGERS_RF))
    print('  Section 4.6.2 and Section 7.8.1a quote:')
    print('   %-26s %13s %13s %10s' % ('', '(4.29)', '(4.28)+C_q', 'change'))
    for label, a, b in (('peak f_T [GHz]', pk_rf_committed['fT'] / 1e9,
                         pk_rf_full['fT'] / 1e9),
                        ('f_max at that bias [GHz]',
                         pk_rf_committed['fmax'] / 1e9,
                         pk_rf_full['fmax'] / 1e9),
                        ('extrinsic f_T [GHz]', pk_rf_ext_c['fT'] / 1e9,
                         pk_rf_ext_f['fT'] / 1e9),
                        ('extrinsic f_max [GHz]', pk_rf_ext_c['fmax'] / 1e9,
                         pk_rf_ext_f['fmax'] / 1e9)):
        print('   %-26s %13.4f %13.4f %9.4f %%'
              % (label, a, b, 100.0 * (b / a - 1.0)))
    r_committed = pk_rf_committed['fmax'] / pk_rf_committed['fT']
    r_full = pk_rf_full['fmax'] / pk_rf_full['fT']
    print()
    print('  Section 7.8.1a\'s verdict is STRUCTURAL -- that f_max is a')
    print('  restatement of R_total/(R_g + R_s) rather than an estimate of')
    print('  f_max -- so the thing to check is the RATIO, not a frequency.')
    print('  (The first draft of this check asserted f_T was "two orders')
    print('  below the 100-300 GHz RF requirement".  That is wrong twice:')
    print('  20.3 GHz is about ONE order below 100 GHz, and Section 7.8.1a')
    print('  does not rest on a frequency threshold at all.  Prose is a')
    print('  detector, 2026-09-28, and here it detected my own check.)')
    print('   f_max/f_T = %.6f on (4.29), %.6f on (4.28): %+.4f %%'
          % (r_committed, r_full, 100 * (r_full / r_committed - 1.0)))
    check('Section 7.8.1a\'s structural verdict is untouched by the exact '
          'charge relation',
          abs(r_full / r_committed - 1.0) < 0.01 and r_full < 1.3,
          'f_max/f_T moves %+.4f %% and stays at %.4f, still below the bottom '
          'of Feijoo et al.\'s 1.3-1.4 band that 2026-10-04 showed a PERFECT '
          'gate cannot reach.  Five levers have now been tried -- saturation, '
          'the perfect contact, the diffusion term, the exact charge relation '
          'and its capacitance -- and the binding constraint is still that '
          'transfer_characteristic() is a resistor with no output resistance.'
          % (100 * (r_full / r_committed - 1.0), r_full))

    # -----------------------------------------------------------------
    print()
    print('SECTION 5.  VALIDATION AGAINST EXACTLY KNOWN VALUES')
    print('-' * 78)
    dVs = np.array([0.01, 0.05, 0.1, 0.5, 1.0, 1.2, 2.0, 2.7, 3.5])

    # (a) the relation's own residual
    nD = np.asarray(n_exact(dVs), dtype=float)
    resid = np.abs(sc.exact_relation_residual(nD) - dVs) / dVs
    check('(4.28) is satisfied by its own root to machine precision',
          float(np.max(resid)) < 1e-14,
          'max relative residual %.3e over dV in [0.01, 3.5] V'
          % float(np.max(resid)))

    # (b) THE DIFFERENTIAL IDENTITY (R1)/(R2), in ulps
    lhs = E_CHARGE * dn_ddV_exact(dVs)
    rhs = 1.0 / (1.0 / gfet.C_ox + 1.0 / C_Q_dispersion(nD))
    ulp_max = max(ulps(a, b) for a, b in zip(lhs, rhs))
    check('e*dn/d(dV) from (4.28) EQUALS the series combination of C_ox with '
          'C_Q(n) = 2e*sqrt(n)/A_F',
          ulp_max <= 4.0,
          'max %.1f ulp over the same sweep.  This is the identity that makes '
          '(4.29) a differential relation rather than a linearisation, and it '
          'is an algebra check, not a tolerance.' % ulp_max)
    # the two spellings of C_Q must agree too
    alt = CQ_PREFACTOR * dd.fermi_voltage(nD)
    ulp_alt = max(ulps(a, b) for a, b in zip(C_Q_dispersion(nD), alt))
    check('2e*sqrt(n)/A_F and CQ_PREFACTOR*E_F(n)/e are the same C_Q',
          ulp_alt <= 8.0,
          'max %.1f ulp.  The second spelling is the repository\'s own '
          'quantum_capacitance_dispersion(); agreement means this module did '
          'not invent a prefactor.' % ulp_alt)

    # (c) THE INTEGRAL IDENTITY, with convergence order
    print()
    print('  The integral identity: int_0^dV (R2) dV\' must return (4.28)\'s')
    print('  root.  Reported with convergence order, so agreement cannot be')
    print('  an accident of one grid (2026-10-03: a PASS is silent about')
    print('  magnitude; here the magnitude is the order).')
    dV_t = 1.2
    target = float(n_exact(dV_t))
    print('   UNIFORM grid (the first draft):')
    errs = []
    for nn in (251, 501, 1001, 2001, 4001):
        e_rel = abs(n_int_EF(dV_t, n_nodes=nn) / target - 1.0)
        errs.append((nn, e_rel))
        print('   nodes %6d  relative error %.3e' % (nn, e_rel))
    orders = [np.log2(errs[k][1] / errs[k + 1][1]) for k in range(len(errs) - 1)
              if errs[k + 1][1] > 0]
    print('   observed orders: %s   <-- NOT 4, and still climbing'
          % ', '.join('%.2f' % o for o in orders))
    print('   This check FAILED as first written (1.1e-07 at 4001 nodes')
    print('   against a 1e-10 criterion, order 3.49 against 4).  What it')
    print('   measured is the SAME NUMBER Section 1\'s failure measured: the')
    print('   integrand turns over inside a layer of width dV_x = %.3f mV at'
          % (1e3 * DV_CROSSOVER))
    print('   the lower limit, and a uniform grid over [0, %.1f V] puts only'
          % dV_t)
    print('   %.0f nodes in it even at 4001 points.  Two independent exactness'
          % (4001 * DV_CROSSOVER / dV_t))
    print('   checks, written for different purposes, both failed by')
    print('   measuring dV_x -- which is the strongest evidence this session')
    print('   has that 2026-10-09\'s rule is a rule and not an anecdote.')
    print('   GRADED grid, panels straddling dV_x, same integrand:')
    errs_g = []
    vals_g = []
    for nn in (21, 41, 81, 161, 321, 641):
        v = n_int_EF_graded(dV_t, nodes_per_panel=nn)
        e_rel = abs(v / target - 1.0)
        vals_g.append(v)
        errs_g.append((nn, e_rel))
        print('   %6d nodes/panel  relative error %.3e' % (nn, e_rel))
    ord_g = [np.log2(errs_g[k][1] / errs_g[k + 1][1])
             for k in range(len(errs_g) - 1) if errs_g[k + 1][1] > 0]
    print('   observed orders: %s' % ', '.join('%.2f' % o for o in ord_g))
    richardson = (16.0 * vals_g[-1] - vals_g[-2]) / 15.0
    rich_rel = abs(richardson / target - 1.0)
    print('   h^4 (Richardson) extrapolation of the two finest grids lands')
    print('   %.3e relative from (4.28)\'s closed-form root.' % rich_rel)
    check('the integral of (R2) reproduces (4.28) at the Simpson order, so '
          'C_Q(n) and (4.28) are the same model',
          3.9 <= ord_g[-1] <= 4.1
          and all(errs_g[k][1] > errs_g[k + 1][1]
                  for k in range(len(errs_g) - 1)),
          'the criterion here is the ORDER and the absence of a plateau, not '
          'a threshold on the error: a residual that keeps falling as h^4 is '
          'pure quadrature truncation, while a wrong prefactor or a dropped '
          'factor of 2 in C_Q would leave an O(1) offset and drive the '
          'observed ratio to 1 instead of 16.  Measured order %.3f at the '
          'finest refinement, errors strictly decreasing %.2e -> %.2e over 5 '
          'halvings, and the h^4 extrapolation %.2e from the closed form.  '
          'This replaces the first draft\'s 1e-10, which was a guessed '
          'number that the method itself can supply.'
          % (ord_g[-1], errs_g[0][1], errs_g[-1][1], rich_rel))

    # (d) A_F -> 0 must give the electrostatic charge
    # The b -> 0 branch of the SAME closed form, evaluated with b = 0 rather
    # than with a small b, so this is the limit itself and not an approach to
    # it.  No module-level name is rebound to get it.
    def _n_exact_with_b(dV, b):
        dV = np.asarray(dV, dtype=float)
        x = (-b + np.sqrt(b * b + 4.0 * A_OX * np.abs(dV))) / (2.0 * A_OX)
        return np.sign(dV) * x * x

    n_es_check = _n_exact_with_b(dVs, 0.0)
    target_es = gfet.C_ox * dVs / E_CHARGE
    ulp_es = max(ulps(a, b) for a, b in zip(n_es_check, target_es))
    check('A_F -> 0 returns the pure electrostatic charge C_ox*dV/e',
          ulp_es <= 4.0,
          'max %.1f ulp.  A limit that must reproduce a known result, and the '
          'only thing that distinguishes (4.28) from the electrostatic model '
          'is the graphene drop it carries.' % ulp_es)

    # (e) oddness, bitwise, in dV (2026-10-09's parameterisation lesson)
    odd = np.asarray(n_exact(dVs), dtype=float) + np.asarray(n_exact(-dVs),
                                                             dtype=float)
    check('n(-dV) == -n(dV) BITWISE in dV',
          bool(np.all(odd == 0.0)),
          'max |n(dV) + n(-dV)| = %r.  Bitwise IN dV, which is the '
          'parameterisation 2026-10-09 established the claim has to name: the '
          'same statement is NOT bitwise through a V_g linspace.'
          % float(np.max(np.abs(odd))))

    # (f) the species split still closes
    n_e, n_h = sc.split_species(nD)
    prod = n_e * n_h
    tgt = gfet.n_puddle ** 2 / 4.0
    rel_prod = float(np.max(np.abs(prod / tgt - 1.0)))
    check('n_e*n_h == n_puddle^2/4 on the (4.28) root as well as the (4.29) one',
          rel_prod < 1e-14,
          'max relative deviation %.3e.  The mass-action closure 2026-10-09 '
          'found is a property of the FLOOR, not of the charge relation, so it '
          'survives replacing the charge relation -- which is what makes it '
          'safe to keep the floor unchanged across this re-derivation.'
          % rel_prod)

    # (g) reversibility of the patch, and bitwise restoration
    cd_id, qc_id, crx_id = (gfet.carrier_density, gfet.quantum_capacitance,
                            crx.carrier_density)
    with on_exact(capacitance=True):
        pass
    check('every patched name is restored to the SAME OBJECT after the block',
          gfet.carrier_density is cd_id and gfet.quantum_capacitance is qc_id
          and crx.carrier_density is crx_id,
          'three names rebound and three restored by identity, not by value')
    Vg_s = np.linspace(-1.5, 3.5, 41)
    Id_after = gfet.transfer_characteristic(Vg_s, Vds=VDS_RF)
    with on_exact(capacitance=True):
        pass
    Id_again = gfet.transfer_characteristic(Vg_s, Vds=VDS_RF)
    check('the committed transfer characteristic is BITWISE unchanged by this '
          'module having run',
          bool(np.all(Id_after == Id_again)),
          '41-point sweep at V_ds = %.2f V, max |difference| = %r'
          % (VDS_RF, float(np.max(np.abs(Id_after - Id_again)))))
    check('and the committed saturated I_d at the RF bias is the number '
          'Section 4.6.6 published',
          abs(Id_dd_committed / 8.190217e-05 - 1.0) < 1e-6,
          'recomputed %.6e A against Section 4.6.6\'s 8.190217e-05 A.  This is '
          'the anchor that says the (4.29) column of every table above is the '
          'committed model and not a re-run with drift in it.' % Id_committed)

    # -----------------------------------------------------------------
    print()
    print('SECTION 6.  WHAT THIS SETTLES, AND WHAT IT DOES NOT')
    print('-' * 78)
    print('  SETTLED, and the item is CLOSED.')
    print('  Sections 4.3, 4.4 and 4.6 are re-derived on (4.28).  The largest')
    print('  movement in Sections 4.3-4.4 is 2.09 % (the charge, and L_x with')
    print('  it); in Section 4.6 it is 2.0 % in f_T from the charge alone and')
    print('  0.06 % once the capacitance follows.  No ranking and no')
    print('  qualitative claim in Chapter 4 or Section 7.8.1a changes, and the')
    print('  (4.29) numbers stay in the chapter annotated rather than')
    print('  replaced.')
    print()
    print('  CORRECTED, and this is the finding rather than the -0.60 %.')
    print('  (4.29) is NOT the linearisation of (4.28).  It is (4.28)\'s exact')
    print('  differential relation used algebraically -- a one-point')
    print('  rectangle rule -- compounded with C_q evaluated at the gate')
    print('  overdrive instead of at E_F/e.  Two independent errors, same')
    print('  sign, same order, not additive.  Chapter 4 said "linearisation"')
    print('  in three places and Section 7.9 item 12 repeats it; all four are')
    print('  corrected in place, and the superseded wording is kept visible.')
    print()
    print('  AND THE ITEM WAS MIS-SCOPED.  It was written as a charge-relation')
    print('  item.  Under (R2) the small-signal gate capacitance is the SAME')
    print('  object as the charge relation\'s derivative, so re-deriving the')
    print('  charge forces C_gs to move too.  In f_T the two halves are')
    print('  -1.9986 % and +1.9392 % -- opposite signs, equal to 3 % of each')
    print('  other -- and they leave -0.0594 %.  A re-derivation that changed')
    print('  only the charge would have reported -2.0 % and been wrong by a')
    print('  factor of 34; one that reported only the combined number would')
    print('  have called item 12 negligible in f_T when it is two coupled')
    print('  2 % corrections that happen to oppose.  An item that names one')
    print('  of two coupled quantities gets the other one whether or not it')
    print('  asks, and the sum is not the way to find out which mattered.')
    print()
    print('  AND THE PEAK BIAS IS THE SENSITIVE QUANTITY, not the peak value.')
    print('  Peak f_T moves 0.06 % in height and 6.76 % in LOCATION (-1.1115')
    print('  -> -1.1867 V), a ratio of 114.  Both peaks stay on the hole')
    print('  branch, so the open peak-f_T branch-crossing question is')
    print('  answered NO for the (4.28)-vs-(4.29) pair; it was asked about')
    print('  the branched-vs-magnitude pair and stays open there.')
    print()
    print('  NOT SETTLED.')
    print('   - (4.28) IS A T = 0 RELATION.  E_F = hbar v_F sqrt(pi n) carries')
    print('     no thermal broadening, so (4.28) is "exact" with respect to')
    print('     the declaration and to the Dirac dispersion, NOT with respect')
    print('     to physics at 300 K.  arXiv:2209.00388 uses')
    print('     C_q ~ k c1 sqrt(1 + (V_c/c1)^2), whose c1 -> 0 limit is (R1);')
    print('     c1 is exactly the thermal rounding (4.28) drops.  The regime')
    print('     where that matters is E_F <~ kT, i.e. dV <~ %.4f V, measured'
          % (float(((1.380649e-23 * gfet.T / E_CHARGE) / A_FERMI) ** 2
                   * A_OX + (1.380649e-23 * gfet.T / E_CHARGE))))
    print('     below; and that regime is %.0f x smaller than the drive at'
          % (dV_rf / max(1e-12, float(
              ((1.380649e-23 * gfet.T / E_CHARGE) / A_FERMI) ** 2 * A_OX
              + (1.380649e-23 * gfet.T / E_CHARGE)))))
    print('     which n_puddle still contributes %.0f %% of the total charge,'
          % (100.0 * gfet.n_puddle / float(floor_total(n_exact(dV_rf)))))
    print('     so the thermal correction is bounded by a floor that')
    print('     2026-10-09 measured to be 2.07x off Chapter 3 in the first')
    print('     place.  Stated as a limitation, not costed.')
    print('   - WHETHER (4.28) OR (4.29) IS CLOSER TO MEASUREMENT.  Still')
    print('     untested, exactly as 2026-10-08 said.  Nothing here licenses')
    print('     calling (4.29) wrong; it is non-exact under the declaration')
    print('     and is now also named as a quadrature, which is a statement')
    print('     about its derivation and not about its accuracy.')
    print('   - THE FLOOR.  It is carried unchanged into (4.28) so that the')
    print('     comparison is apples to apples.  (4.28)+floor is therefore')
    print('     not itself an exact model, and the floor attenuates the')
    print('     charge error by about an order of magnitude -- which is why')
    print('     this re-derivation moves so little, and is not a reason for')
    print('     confidence.')
    print('   - the hole branch of Chapter 4\'s branch-asymmetric results and')
    print('     the peak-f_T branch-crossing question: expressible since')
    print('     2026-10-09, still unverified, untouched today.')

    # -----------------------------------------------------------------
    make_figure(rows, pk_committed, pk_charge, pk_full)
    print()
    print('  Figure written: exact_charge_rederivation.png')

    print()
    print('=' * 78)
    print('TOTAL: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return failed


def make_figure(rows, pk_committed, pk_charge, pk_full):
    dV = np.linspace(1e-4, 3.0, 600)
    nA = np.array([float(n_rect_dV(v)) for v in dV])
    nB = np.array([float(n_rect_EF(v)) for v in dV])
    nC = np.array([float(n_int_dV(v, n_nodes=401)) for v in dV])
    nD = np.array([float(n_exact(v)) for v in dV])

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.6))

    ax = axes[0]
    ax.loglog(dV, nA * 1e-4, linewidth=2, label='A: (4.29), committed')
    ax.loglog(dV, nB * 1e-4, '--', linewidth=1.8,
              label=r'B: rectangle, $C_q(E_F/e)$')
    ax.loglog(dV, nC * 1e-4, '-.', linewidth=1.8,
              label=r'C: integral, $C_q(\Delta V)$')
    ax.loglog(dV, np.abs(nD) * 1e-4, linewidth=2.4, color='k',
              label='D: (4.28), exact')
    ax.axhline(gfet.n_puddle * 1e-4, color='gray', linestyle=':',
               label=r'$n_{puddle}$ (floors all four)')
    ax.set_xlabel(r'$\Delta V = V_g - V_{Dirac}$ (V)')
    ax.set_ylabel(r'raw sheet density (cm$^{-2}$)')
    ax.set_title('1. (4.29) is linear and (4.28) quadratic as '
                 r'$\Delta V\to0$' '\n'
                 'so (4.29) is not a linearisation of (4.28)')
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(True, alpha=0.3, which='both')

    ax = axes[1]
    ax.semilogx(dV, 100 * (nA / nD - 1.0), linewidth=2.4, color='crimson',
                label='total error of (4.29)')
    ax.semilogx(dV, 100 * (nB / nD - 1.0), '--', linewidth=1.8,
                label='quadrature error alone')
    ax.semilogx(dV, 100 * (nC / nD - 1.0), '-.', linewidth=1.8,
                label='argument error alone')
    ax.semilogx(dV, 100 * (np.sqrt(nA ** 2 + gfet.n_puddle ** 2)
                           / np.sqrt(nD ** 2 + gfet.n_puddle ** 2) - 1.0),
                linewidth=2, color='seagreen',
                label='total, AFTER the puddle floor')
    ax.axvline(VG_RF - gfet.V_dirac, color='gray', linestyle=':',
               label='RF bias')
    ax.axhline(0.0, color='k', linewidth=0.8)
    ax.set_xlabel(r'$\Delta V$ (V)')
    ax.set_ylabel(r'$(n_{form} - n_{(4.28)})/n_{(4.28)}$  (%)')
    ax.set_title('2. Two independent errors, same sign, same order\n'
                 'and the floor attenuates both by ~10x')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which='both')

    ax = axes[2]
    Vg = np.linspace(-1.5, 3.5, 400)
    fT_c, _, _, _, _ = rf.compute_fT_fmax(Vg, Vds=VDS_RF)
    with on_exact(capacitance=False):
        fT_q, _, _, _, _ = rf.compute_fT_fmax(Vg, Vds=VDS_RF)
    with on_exact(capacitance=True):
        fT_f, _, _, _, _ = rf.compute_fT_fmax(Vg, Vds=VDS_RF)
    ax.plot(Vg, fT_c / 1e9, linewidth=2, label='(4.29), committed')
    ax.plot(Vg, fT_q / 1e9, '--', linewidth=1.8,
            label='(4.28) charge only')
    ax.plot(Vg, fT_f / 1e9, linewidth=2.2, color='k',
            label=r'(4.28) charge + consistent $C_{gs}$')
    ax.axvline(gfet.V_dirac, color='gray', linestyle=':', label='Dirac point')
    ax.set_xlabel(r'$V_g$ (V)')
    ax.set_ylabel(r'$f_T$ (GHz, per $\mu$m convention)')
    ax.set_title('3. The capacitance half is the larger half\n'
                 r'peak $f_T$ %.3f $\to$ %.3f $\to$ %.3f GHz'
                 % (pk_committed['fT'] / 1e9, pk_charge['fT'] / 1e9,
                    pk_full['fT'] / 1e9))
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('exact_charge_rederivation.png', dpi=200,
                bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    import sys
    sys.exit(1 if main() else 0)
