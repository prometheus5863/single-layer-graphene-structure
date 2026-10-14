"""
graphene_quantum_capacitance_argument.py

Chapter 7 Section 7.9 ITEM 16 -- "fixing quantum_capacitance()'s argument",
created 2026-10-10 as the first item on that list to arrive already costed
(+4.7416 % of the raw charge at the RF bias, the largest cost attached to
anything then open on Chapter 4).

The item as written says: C_q is evaluated at the gate overdrive dV, while the
dispersion it is derived from requires E_F/e = A_F*sqrt(n).  Fix the argument.

Acting on it produced two things the item did not anticipate, and both are the
point of this module:

  (1) THE ARGUMENT CANNOT BE FIXED BY CHANGING AN ARGUMENT.  The correct
      argument is E_F/e = A_F*sqrt(n), and n is what C_q is being used to
      compute.  quantum_capacitance(V) takes a VOLTAGE; the quantity it must
      be evaluated at is a STATE.  Any fix is therefore a self-consistent
      solve, not an edit to a call site -- which is why the item has sat as a
      one-line description for a session while being the largest open cost.

  (2) THE OBVIOUS FIX SILENTLY DELETES THE FINITE-TEMPERATURE PHYSICS.  The
      dispersion form C_Q(n) = 2e*sqrt(n)/A_F is a T = 0 expression.  The
      committed quantum_capacitance() exists in its thermally-rounded form for
      a stated reason -- its own docstring says "so C_q does not diverge to
      zero exactly at the neutrality point".  Substituting C_Q(n) for it
      removes that rounding, and C_q(0) = 0 makes the series factor, and hence
      the raw induced charge, EXACTLY ZERO at the Dirac point.  The
      2026-10-10 sizing of this item used precisely that substitution.  So the
      +4.7416 % it quotes is the cost of a fix that is wrong in a second way,
      and the number has to be re-measured on a form that fixes the argument
      WITHOUT throwing away the temperature.

That form is built here and called (B_T):

      C_q = PREFACTOR * (kT/e) * ln(2*(1 + cosh(e*E_F/kT))),  E_F/e = A_F*sqrt(n)

i.e. the committed thermal function evaluated at the right argument, with n
solved self-consistently.  It reduces to the committed form's own large-|eta|
limit, to the T = 0 dispersion as T -> 0, and keeps a finite C_q at n = 0.
Each of those three is an EXACTLY known value and is checked as one.

Validation policy follows 2026-09-17: at least one check against a value that
is known exactly, not a plausible-range check.  There are five here.
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
import graphene_exact_charge_rederivation as X
import graphene_binding_class as bc

E_CHARGE = 1.602176634e-19
KB = 1.380649e-23
A_FERMI = dd.A_FERMI              # E_F/e = A_FERMI*sqrt(n)   [V m]
CQ_PREFACTOR = dd.CQ_PREFACTOR    # 2 e^3/(pi hbar^2 v_F^2)   [F/(m^2 V)]
A_OX = E_CHARGE / gfet.C_ox       # 'a' in a*n + b*sqrt(n) = dV

VDS_RF = dd.VDS_REF               # 0.1 V, the RF bias of record
VG_RF = 2.0
VG_ON = 3.5
DV_RF = VG_RF - gfet.V_dirac      # 1.2 V

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
# The capacitance forms, as functions of their OWN natural argument
# ---------------------------------------------------------------------------

def cq_thermal(V, T=None):
    """The committed quantum_capacitance(), re-expressed so that the argument
    is explicit: a POTENTIAL in volts.  Whether that potential is the gate
    overdrive or E_F/e is the whole of item 16."""
    T = gfet.T if T is None else T
    kT_over_e = KB * T / E_CHARGE
    eta = np.asarray(V, dtype=float) / kT_over_e
    return CQ_PREFACTOR * kT_over_e * np.log(2.0 * (1.0 + np.cosh(eta)))


def cq_disp(n):
    """C_Q(n) = e^2 dn/dE_F = 2e*sqrt(n)/A_F.  A T = 0 expression."""
    return 2.0 * E_CHARGE * np.sqrt(np.abs(np.asarray(n, dtype=float))) / A_FERMI


def cq_thermal_of_n(n, T=None):
    """FORM B_T's capacitance: the thermal form at the CORRECT argument."""
    return cq_thermal(A_FERMI * np.sqrt(np.abs(np.asarray(n, dtype=float))), T=T)


def series_factor(C_q):
    return C_q / (C_q + gfet.C_ox)


# ---------------------------------------------------------------------------
# The three charge forms that differ ONLY in C_q's argument
# ---------------------------------------------------------------------------

def n_arg_dV(dV):
    """FORM A -- committed (4.29): thermal C_q at the gate overdrive."""
    dV = np.asarray(dV, dtype=float)
    return (gfet.C_ox * dV / E_CHARGE) * series_factor(cq_thermal(dV))


def _fixed_point(dV, cq_of_n, iters=600):
    """Solve n = (C_ox dV/e) * f(C_q(n)) by iteration from the electrostatic
    value.  The map is increasing in n and bounded by n_es, so from above it
    is a monotone decreasing bounded sequence."""
    dV = np.asarray(dV, dtype=float)
    n_es = gfet.C_ox * dV / E_CHARGE
    n = np.array(n_es, dtype=float, copy=True)
    for _ in range(iters):
        n_new = n_es * series_factor(cq_of_n(n))
        if np.all(np.abs(n_new - n) <= 1e-16 * np.abs(n_new) + 1e-9):
            n = n_new
            break
        n = n_new
    return n


def n_arg_EF_T0(dV):
    """FORM B -- argument fixed with the T = 0 dispersion.  This is the form
    2026-10-10 used to size item 16 at +4.7416 %; it is also the form that
    deletes the thermal rounding."""
    return _fixed_point(dV, cq_disp)


def n_arg_EF_thermal(dV):
    """FORM B_T -- argument fixed, temperature KEPT.  The honest fix."""
    return _fixed_point(dV, cq_thermal_of_n)


def fixed_point_residual(dV, nfun, cq_of_n):
    """|n - n_es*f(C_q(n))| / |n|, the solve's own error."""
    n = np.asarray(nfun(dV), dtype=float)
    n_es = gfet.C_ox * np.asarray(dV, dtype=float) / E_CHARGE
    g = n_es * series_factor(cq_of_n(n))
    den = np.where(np.abs(n) > 0, np.abs(n), 1.0)
    return np.max(np.abs(g - n) / den)


# ---------------------------------------------------------------------------
# Installing FORM B_T into Chapter 4
# ---------------------------------------------------------------------------

# ONE COUNTER PER PATCHED NAME.  The first draft of Section 6's arrival
# control had a single counter, incremented only in the capacitance patch, and
# FAILED while nothing was wrong: census() exercises the CHARGE patch and
# never asks for C_gs, so zero capacitance calls is the correct answer there.
# A control that counts a different arrival from the one it gates reports a
# fault that does not exist -- the 2026-10-10 rule ("a control must not read
# the thing it controls for") in its other failure mode.  Each patched name
# now has its own counter and is checked where that name is actually used.
_ARRIVED = {'cq': 0, 'charge': 0}


def _quantum_capacitance_consistent(V_g_minus_Vdirac, T=None):
    """Drop-in for gfet.quantum_capacitance with the argument corrected.

    Same signature, so rf_small_signal_model.gate_capacitance() picks up the
    consistent C_gs without being touched.  It takes the overdrive, solves for
    n, and returns C_q at that state's E_F/e."""
    _ARRIVED['cq'] += 1
    n = n_arg_EF_thermal(V_g_minus_Vdirac)
    return cq_thermal_of_n(n, T=T)


def _carrier_density_consistent(V_g, V_ch=0.0):
    _ARRIVED['charge'] += 1
    dV = (np.asarray(V_g, dtype=float) - np.asarray(V_ch, dtype=float)
          - gfet.V_dirac)
    return X.floor_total(n_arg_EF_thermal(dV))


@contextlib.contextmanager
def on_consistent_argument(capacitance=True):
    """Run a block with C_q's argument corrected.

    capacitance=False changes only the charge; True also makes the
    small-signal gate capacitance consistent.  Separated because 2026-10-10
    found those two halves nearly cancel in f_T, and a single combined number
    would hide that.

    REWRITTEN 2026-10-14.  The two lines this manager used to carry --
    `gfet.carrier_density` and `crx.carrier_density` -- were a hand-maintained
    list, and a hand-maintained list is exactly what was wrong in all three of
    the 2026-10 delivery faults.  The sweep of 2026-10-14 measured what such a
    list leaves behind repository-wide (25 stale sites at 14 rebinding sites,
    22 of them never found by hand) and its finding (i) is that the set to
    rebind is the name's RUNTIME IDENTITY CLASS, which is discoverable and
    needs no list.  graphene_binding_class.rebind() computes it at entry, so
    any binding of these names added in future is rebound automatically
    instead of being found by the next audit.

    The two names this installs are NOT symmetric, and that is finding (ii):
    `carrier_density`'s class is rebound whole, while
    `graphene_contact_doping_model.quantum_capacitance` is excluded BY NAME
    with its reason recorded in bc.CONTACT_OPERATING_POINT -- it is reached at
    a contact whose series capacitance is C_interface, and the replacement
    installed here solves for n against C_ox.  Copying this object there
    would be a different wrong answer, not a repair.

    Behaviour is unchanged by the rewrite: the class of `carrier_density` is
    exactly {gfet, crx} and the class of `quantum_capacitance` minus the
    exclusion is exactly {gfet}, i.e. precisely the three sites the old lines
    named.  That is verified by transcript rather than asserted -- this
    module's committed output reproduces byte-identically across the change.
    """
    with contextlib.ExitStack() as stack:
        stack.enter_context(bc.rebind(gfet, 'carrier_density',
                                      _carrier_density_consistent))
        if capacitance:
            stack.enter_context(
                bc.rebind(gfet, 'quantum_capacitance',
                          _quantum_capacitance_consistent,
                          exclude_reasons=bc.CONTACT_OPERATING_POINT))
        yield


def crossover_lengths(metals=(('Pd', None), ('300', 300.0), ('500', 500.0))):
    """L_x for the per-width contact resistances Chapter 4 Section 4.4 quotes,
    on whatever charge model is currently installed."""
    out = {}
    for name, rc_um in metals:
        rc = gfet.Rc_total if rc_um is None else 2 * (rc_um * 1e-6) / gfet.W
        out[name] = float(crx.crossover_length(rc, VG_ON))
    return out


def census():
    """Every committed number that depends on the charge relation."""
    d = {}
    d['n_3V5'] = float(gfet.carrier_density(VG_ON))
    d['R_dirac'] = float(gfet.channel_resistance(gfet.V_dirac))
    d.update(('L_%s' % k, v) for k, v in crossover_lengths().items())
    ion, ioff, ratio = X.onoff()
    d['onoff_term'] = ratio
    ion0, ioff0, ratio0 = X.onoff(Rc=0.0)
    d['onoff_intr'] = ratio0
    d['Id_rf'] = float(gfet.transfer_characteristic(
        np.array([VG_RF]), Vds=VDS_RF)[0])
    return d


# ---------------------------------------------------------------------------
# The overflow-free spelling of the thermal rounding (Section 7.9 item 13)
# ---------------------------------------------------------------------------

def smoothed_V_committed(eta):
    """ln(2(1+cosh eta)) -- the committed spelling.  cosh overflows near
    eta = 710 and the whole expression returns inf."""
    return np.log(2.0 * (1.0 + np.cosh(np.asarray(eta, dtype=float))))


def smoothed_V_stable(eta):
    """The same function, algebraically rearranged:

        2(1 + cosh eta) = e^eta + 2 + e^-eta = (e^(eta/2) + e^(-eta/2))^2
                        = 4 cosh^2(eta/2)
        => ln(2(1+cosh eta)) = |eta| + 2*log1p(exp(-|eta|))

    An IDENTITY, not an approximation, so it must agree with the committed
    spelling to within rounding wherever that one is finite -- and it is
    finite everywhere.  log1p and exp(-|eta|) are both well-conditioned."""
    a = np.abs(np.asarray(eta, dtype=float))
    return a + 2.0 * np.log1p(np.exp(-a))


def cq_thermal_stable(V, T=None):
    T = gfet.T if T is None else T
    kT_over_e = KB * T / E_CHARGE
    return CQ_PREFACTOR * kT_over_e * smoothed_V_stable(
        np.asarray(V, dtype=float) / kT_over_e)


def _overflow_onset():
    """Largest eta for which the committed spelling is finite, by bisection."""
    lo, hi = 1.0, 1.0e4
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if np.isfinite(smoothed_V_committed(mid)):
            lo = mid
        else:
            hi = mid
    return lo


def n_exact_quad_arg(dV, cq_of_n=None, n0=0.0, V0=0.0, rtol=1e-12, atol=1e-6):
    """FORM D_T -- solve dn/d(dV) = (C_ox/e)*f(C_q(n)) as an initial-value
    problem.  With cq_of_n = cq_disp this must reproduce (4.28)'s closed-form
    root; with cq_thermal_of_n it is the finite-T reference that holds BOTH
    knobs right and holds the temperature fixed at the committed 300 K."""
    from scipy.integrate import solve_ivp
    cq_of_n = cq_thermal_of_n if cq_of_n is None else cq_of_n
    dV = float(dV)
    if dV == V0:
        return float(n0)
    s = -1.0 if dV < 0.0 else 1.0
    rhs = lambda V, n: (gfet.C_ox / E_CHARGE) * series_factor(cq_of_n(n))
    sol = solve_ivp(rhs, [V0, abs(dV)], [float(n0)], rtol=rtol, atol=atol,
                    method='Radau')
    return s * float(sol.y[0, -1])


def main():
    print('=' * 78)
    print("FIXING quantum_capacitance()'s ARGUMENT -- Section 7.9 item 16")
    print('and the overflow that had to be fixed first -- item 13')
    print('=' * 78)
    print('  A_F = %.6e V*m    PREFACTOR = %.6e F/(m^2 V)'
          % (A_FERMI, CQ_PREFACTOR))
    print('  C_ox = %.6e F/m^2   kT/e = %.6e V at T = %.1f K'
          % (gfet.C_ox, KB * gfet.T / E_CHARGE, gfet.T))
    print('  RF bias of record: V_g = %.1f V, V_ds = %.1f V, dV = %.1f V'
          % (VG_RF, VDS_RF, DV_RF))

    print()
    print('SECTION 1.  THE ARGUMENT CANNOT BE FIXED BY CHANGING AN ARGUMENT')
    print('-' * 78)
    print('  quantum_capacitance(V) takes a POTENTIAL.  The dispersion needs')
    print('  E_F/e = A_F*sqrt(n), and n is the quantity C_q is being used to')
    print('  compute.  The correct argument is a STATE, so item 16 is a')
    print('  self-consistency problem, not an edit to a call site.  That is')
    print('  why it sat unacted-on for a session while being the largest')
    print('  open cost on Chapter 4.')
    print()
    print('  At the RF bias, on the committed charge:')
    n_A = float(n_arg_dV(DV_RF))
    EF_A = float(A_FERMI * np.sqrt(abs(n_A)))
    print('    n = %.6e /m^2  ->  E_F/e = %.6f V, which is %.4f %% of the'
          % (n_A, EF_A, 100 * EF_A / DV_RF))
    print('    %.1f V overdrive.  C_q at the two arguments:' % DV_RF)
    print('      C_q(dV)     = %.6e F/m^2' % float(cq_thermal(DV_RF)))
    print('      C_q(E_F/e)  = %.6e F/m^2  (factor %.3f)'
          % (float(cq_thermal(EF_A)),
             float(cq_thermal(DV_RF) / cq_thermal(EF_A))))
    print('      C_ox        = %.6e F/m^2' % gfet.C_ox)
    check('the two candidate arguments differ by a factor of ~20, and BOTH '
          'give C_q >> C_ox, which is why a factor 20 in C_q is only a few '
          'per cent in the charge',
          float(cq_thermal(EF_A)) > 5.0 * gfet.C_ox,
          'C_q(E_F/e)/C_ox = %.2f, so the series factor is %.6f against '
          '%.6f -- a %.4f %% difference in the factor from a %.1fx difference '
          'in C_q.  The error is small because the device is NOT '
          'quantum-capacitance-limited at this bias, not because the argument '
          'is nearly right.'
          % (float(cq_thermal(EF_A)) / gfet.C_ox,
             float(series_factor(cq_thermal(EF_A))),
             float(series_factor(cq_thermal(DV_RF))),
             100 * (float(series_factor(cq_thermal(DV_RF)))
                    / float(series_factor(cq_thermal(EF_A))) - 1.0),
             float(cq_thermal(DV_RF) / cq_thermal(EF_A))))
    r = fixed_point_residual(DV_RF, n_arg_EF_thermal, cq_thermal_of_n)
    check('the self-consistent solve converges to its own fixed point at '
          'machine level, so nothing below is an artefact of the iteration',
          r < 1e-14,
          'max |n - n_es*f(C_q(n))|/|n| = %.3e over the solve at the RF bias.'
          % r)

    print()
    print('SECTION 2.  THE OBVIOUS FIX DELETES THE TEMPERATURE')
    print('-' * 78)
    print('  C_Q(n) = 2e*sqrt(n)/A_F is a T = 0 expression.  The committed')
    print('  function is thermally rounded for a stated reason -- its own')
    print('  docstring: "so C_q does not diverge to zero exactly at the')
    print('  neutrality point".  Substituting the dispersion removes that.')
    print()
    print('   %-34s %16s' % ('C_q at the Dirac point', 'F/m^2'))
    print('   %-34s %16.8e' % ('committed thermal, at dV = 0',
                               float(cq_thermal(0.0))))
    print('   %-34s %16.8e' % ('T = 0 dispersion, at n = 0',
                               float(cq_disp(0.0))))
    print('   %-34s %16.8e' % ('thermal at E_F/e, at n = 0',
                               float(cq_thermal_of_n(0.0))))
    closed = CQ_PREFACTOR * (KB * gfet.T / E_CHARGE) * np.log(4.0)
    check('C_q at the Dirac point has a CLOSED FORM, PREFACTOR*(kT/e)*ln 4, '
          'and the thermal-at-E_F/e form reproduces it to the last bit while '
          'the dispersion gives exactly zero',
          ulps(cq_thermal_of_n(0.0), closed) <= 1.0
          and float(cq_disp(0.0)) == 0.0,
          'ln(2(1+cosh 0)) = ln 4 exactly, so C_q(0) = %.10e F/m^2; the '
          'thermal-at-E_F/e form gives it to %.1f ulp.  The dispersion gives '
          '0, hence a series factor of 0 and an induced charge of EXACTLY 0 '
          'at the Dirac point -- the puddle floor would then be supplying '
          '100 %% of the density there rather than its committed 99.9 %%.  '
          'The 2026-10-10 sizing of item 16 used that substitution, so its '
          'number is the cost of a fix that is wrong in a second way.'
          % (closed, ulps(cq_thermal_of_n(0.0), closed)))
    print()
    print('  A sharper statement of the same thing.  With the dispersion,')
    print('  n == 0 is a solution of dn/d(dV) = (C_ox/e)*f(C_Q(n)) through')
    print('  the origin, because f(C_Q(0)) = 0: the initial-value problem is')
    print('  degenerate there.  With the thermal form it is not.  Temperature')
    print('  is what makes the charge relation well-posed at the Dirac point,')
    print('  and that is the physics the obvious fix throws away.')

    print()
    print('SECTION 3.  ITEM 13 BLOCKS ITEM 16, AND ITEM 13 WAS RANKED LAST')
    print('-' * 78)
    print('  The natural exact-limit check on the thermal-at-E_F/e form is')
    print('  T -> 0: it must reproduce the T = 0 dispersion result.  It')
    print('  cannot be run on the committed spelling.  eta = e*V/kT grows as')
    print('  T falls, cosh(eta) overflows, ln(2(1+inf)) = inf, and the series')
    print('  factor returns inf/inf = nan.')
    onset = _overflow_onset()
    T_onset = E_CHARGE * EF_A / (KB * onset)
    print()
    print('    overflow onset:  eta = %.4f' % onset)
    print('    at the RF bias E_F/e = %.6f V, that is T = %.4f K'
          % (EF_A, T_onset))
    print('    at T = 300 K it needs V = %.4f V of overdrive'
          % (onset * KB * gfet.T / E_CHARGE))
    check('item 13 -- added 2026-10-08, called latent, ranked LAST -- is what '
          'gates the exact-limit validation of item 16, the item ranked as '
          'the largest open cost',
          not np.isfinite(float(smoothed_V_committed(onset * 1.01)))
          and np.isfinite(float(smoothed_V_stable(1.0e6))),
          'the committed spelling returns inf for eta > %.4f and the stable '
          'one is finite at eta = 1e6.  Latency was measured against the '
          'DEVICE bias range, where %.2f V of overdrive is unreachable, and '
          'not against the checks the repository would later want to run.  A '
          'defect that cannot be reached by any input the model is used at '
          'can still be reached by a VALIDATION of the model, and ranking by '
          'reachability missed that.'
          % (onset, onset * KB * gfet.T / E_CHARGE))
    etas = np.array([0.0, 1e-8, 1e-3, 0.5, 1.0, 5.0, 20.0, 46.4, 100.0,
                     300.0, 700.0])
    worst = max(ulps(smoothed_V_committed(x), smoothed_V_stable(x))
                for x in etas)
    check('the stable spelling is an IDENTITY, not an approximation: '
          '|eta| + 2*log1p(exp(-|eta|)) == ln(2(1+cosh eta)) because '
          '2(1+cosh eta) = 4cosh^2(eta/2)',
          worst <= 1.0,
          'worst disagreement over eta in [0, 700] is %.1f ulp.  So '
          'installing it changes no committed number -- the output-neutrality '
          'pass item 13 was waiting for is discharged by the algebra, not by '
          'a tolerance.' % worst)
    T_small = 1.0e-3
    n_cold = _fixed_point(DV_RF, lambda n: cq_thermal_of_n(n, T=T_small))
    n_cold_stable = _fixed_point(
        DV_RF, lambda n: cq_thermal_stable(
            A_FERMI * np.sqrt(np.abs(np.asarray(n, dtype=float))), T=T_small))
    n_T0 = float(n_arg_EF_T0(DV_RF))
    check('EXACT LIMIT, now runnable: as T -> 0 the thermal-at-E_F/e form '
          'must reproduce the T = 0 dispersion result, and on the committed '
          'spelling it returns nan instead',
          not np.isfinite(float(n_cold))
          and abs(float(n_cold_stable) / n_T0 - 1.0) < 1e-12,
          'committed spelling at T = %.0e K: %s.  Stable spelling: %.10e '
          'against the dispersion fixed point %.10e, relative %+.3e.  This is '
          'the validation item 16 needed and item 13 was preventing.'
          % (T_small, float(n_cold), float(n_cold_stable), n_T0,
             float(n_cold_stable) / n_T0 - 1.0))

    print()
    print('SECTION 4.  THE COST, RE-MEASURED AGAINST A REFERENCE AT THE SAME')
    print('            TEMPERATURE AS THE MODEL')
    print('-' * 78)
    print('  2026-10-10 decomposed (4.29) against (4.28).  (4.28) is a T = 0')
    print('  relation and the committed model is a 300 K one, so each of')
    print('  those three percentages contains a temperature difference as')
    print('  well as the knob it names.  The fix is a 2x2 at ONE temperature:')
    print('  argument in {dV, E_F/e} x quadrature in {rectangle, exact},')
    print('  thermal C_q throughout.')
    print()
    V0, n0 = 1.0e-3, float(abs(X.n_exact(1.0e-3)))
    n_T0_ode = n_exact_quad_arg(DV_RF, cq_of_n=cq_disp, n0=n0, V0=V0)
    n_T0_cf = float(abs(X.n_exact(DV_RF)))
    check('the integrator that supplies the reference is validated against an '
          'EXACTLY known value: with the T = 0 dispersion it must reproduce '
          "(4.28)'s closed-form root",
          abs(n_T0_ode / n_T0_cf - 1.0) < 1e-10,
          'ODE %.10e vs closed form %.10e, relative %+.3e.  Started at '
          'dV = %.0e V rather than 0 because n == 0 solves the T = 0 problem '
          'through the origin (Section 2).'
          % (n_T0_ode, n_T0_cf, n_T0_ode / n_T0_cf - 1.0, V0))
    A = float(n_arg_dV(DV_RF))
    BT = float(n_arg_EF_thermal(DV_RF))
    CT = float(X.n_int_dV(DV_RF))
    DT = n_exact_quad_arg(DV_RF)
    print()
    print('  Raw charge at dV = %.1f V, all four cells at T = 300 K:' % DV_RF)
    print('   %-38s %16s %12s' % ('', 'n [1/m^2]', 'vs D_T'))
    for lbl, v in (('A   C_q(dV),    rectangle  [committed]', A),
                   ('B_T C_q(E_F/e), rectangle', BT),
                   ('C_T C_q(dV),    exact quadrature', CT),
                   ('D_T C_q(E_F/e), exact quadrature', DT)):
        print('   %-38s %16.8e %11.4f %%' % (lbl, v, 100 * (v / DT - 1.0)))
    arg_err = 100 * (CT / DT - 1.0)
    quad_err = 100 * (BT / DT - 1.0)
    both_err = 100 * (A / DT - 1.0)
    print()
    print('  Attribution, both knobs at 300 K:')
    print('    argument error   (quadrature already exact) = %+.4f %%'
          % arg_err)
    print('    quadrature error (argument already exact)   = %+.4f %%'
          % quad_err)
    print('    both wrong                                  = %+.4f %%'
          % both_err)
    print()
    print('  2026-10-10, against the T = 0 relation:')
    print('    argument error = +4.7416 %, quadrature = +2.6396 %, '
          'both = +5.2064 %')
    check('every number in the 2026-10-10 decomposition is INFLATED, because '
          'a 300 K model was differenced against a T = 0 reference',
          arg_err < 4.7416 and quad_err < 2.6396 and both_err < 5.2064,
          'argument %+.4f %% against +4.7416 %% (%.2fx), quadrature %+.4f %% '
          'against +2.6396 %% (%.2fx), both %+.4f %% against +5.2064 %% '
          '(%.2fx).  The ORDERING survives -- the argument is still the '
          'larger error and item 16 is still ahead of the quadrature -- but '
          "item 16's headline +4.7416 %% does not, and the thesis quotes it.  "
          'Superseded in place, not deleted.'
          % (arg_err, 4.7416 / arg_err, quad_err, 2.6396 / quad_err,
             both_err, 5.2064 / both_err))
    check('AND A TRAP: the argument error at 300 K (%.4f %%) agrees with the '
          '2026-10-10 QUADRATURE error (2.6396 %%) to %.4f %%, and they are '
          'different quantities' % (arg_err, abs(arg_err - 2.6396)),
          abs(arg_err - 2.6396) < 0.01,
          'one is the argument error measured at the model temperature, the '
          'other the quadrature error measured at T = 0.  Anyone reading the '
          'two sessions in sequence would take the agreement for a '
          'reproduction.  By 2026-10-02 a coincidence of two artefacts is '
          'not a validation; this is the first time that rule has had to be '
          'applied to a NUMERICAL coincidence rather than to a re-run.')

    print()
    print('SECTION 5.  THE TWO-KNOB CHECK (the 2026-10-10 item), FIRST USE')
    print('-' * 78)
    print('  The item was created because a number can be small as a result')
    print('  of two larger numbers opposing, and every check here asks about')
    print('  one quantity at a time.  Its prescribed form: move each input')
    print('  separately, and flag the case where the combined move is SMALLER')
    print("  than either.  Today's instance has the OTHER sign.")
    print()
    ma0 = 100 * (float(X.n_rect_EF(DV_RF)) / float(X.n_rect_dV(DV_RF)) - 1.0)
    mb0 = 100 * (float(X.n_int_dV(DV_RF)) / float(X.n_rect_dV(DV_RF)) - 1.0)
    mboth0 = 100 * (float(abs(X.n_exact(DV_RF)))
                    / float(X.n_rect_dV(DV_RF)) - 1.0)
    ma = 100 * (BT / A - 1.0)
    mb = 100 * (CT / A - 1.0)
    mboth = 100 * (DT / A - 1.0)
    for frame, x1, x2, x3 in (('300 K, consistent', ma, mb, mboth),
                              ('T = 0 (2026-10-10)', ma0, mb0, mboth0)):
        print('   %-20s arg %+8.4f %%  quad %+8.4f %%  sum %+8.4f %%  '
              'both %+8.4f %%  ratio %6.3f'
              % (frame, x1, x2, x1 + x2, x3, x3 / (x1 + x2)))
    ratio_300 = mboth / (ma + mb)
    ratio_0 = mboth0 / (ma0 + mb0)
    check('the two knobs REINFORCE: the combined move is LARGER than the sum '
          'of the separate moves, so the check as specified -- flag when the '
          'combined move is smaller than either -- would not have flagged it',
          ratio_300 > 1.0 and abs(ratio_300 - 1.0) > 0.05,
          'combined/(sum of separate) = %.4f at 300 K.  The repair is to '
          'report the RATIO and flag any departure from 1 in either '
          'direction, not to test for smallness.  The item is amended rather '
          'than closed.' % ratio_300)
    check('and the INTERACTION MEASURE ITSELF depends on the reference: the '
          'T = 0 frame inflates the departure from additivity by %.1fx'
          % (abs(ratio_0 - 1.0) / abs(ratio_300 - 1.0)),
          abs(ratio_0 - 1.0) > 2.0 * abs(ratio_300 - 1.0),
          'ratio %.4f at T = 0 against %.4f at 300 K -- departures from 1 of '
          '%.4f and %.4f.  So a two-knob check inherits the reference error '
          'of whatever it is differenced against, and must be run at the '
          "model's own operating point.  That is a requirement on the item "
          'that could only be found by running it.'
          % (ratio_0, ratio_300, abs(ratio_0 - 1.0), abs(ratio_300 - 1.0)))

    print()
    print('SECTION 6.  WHAT MOVES IN CHAPTERS 4, 5 AND 6')
    print('-' * 78)
    print('  FORM B_T installed, with the arrival control of 2026-09-30 and')
    print("  2026-10-10's both-names rebinding.  Raw charge moves %+.4f %%;"
          % ma)
    print('  what reaches the committed numbers is this:')
    print()
    base = census()
    _ARRIVED['charge'] = 0
    _ARRIVED['cq'] = 0
    with on_consistent_argument(capacitance=True):
        moved = census()
        hits_charge, hits_cq = _ARRIVED['charge'], _ARRIVED['cq']
        probe = float(gfet.carrier_density(VG_ON))
    check('ARRIVAL CONTROL on the CHARGE patch, which is the name Section 6 '
          'gates: it was called inside the block, and it reached '
          'contact_resistance_crossover, which holds its own import-time '
          'binding',
          hits_charge > 0 and abs(probe / base['n_3V5'] - 1.0) > 1e-12
          and abs(moved['L_Pd'] / base['L_Pd'] - 1.0) > 1e-12,
          '%d calls to the corrected charge inside the block; n(3.5 V) moved '
          '%+.4f %% and L_x(Pd) moved %+.4f %%.  A patch that had not arrived '
          'would have left both at exactly zero, which is what the '
          '2026-10-10 fault looked like.'
          % (hits_charge, 100 * (probe / base['n_3V5'] - 1.0),
             100 * (moved['L_Pd'] / base['L_Pd'] - 1.0)))
    check("and the CAPACITANCE patch is correctly NOT called here, which is "
          'why the single-counter draft of this control failed on a healthy '
          'run',
          hits_cq == 0,
          'census() computes charge, resistances, currents and ratios and '
          'never asks for C_gs, so %d calls to the corrected C_q is the right '
          'answer in this block.  The first draft counted only this name and '
          'reported the charge patch as not arrived.  A control must be '
          'instrumented on the name it gates; Section 7, which does use '
          'C_gs, carries the counter for this one.' % hits_cq)
    print()
    labels = [('n_3V5', 'n at 3.5 V [1/m^2]', 1.0),
              ('R_dirac', 'R_ch at the Dirac point [Ohm]', 1.0),
              ('L_Pd', 'L_x Pd [nm]', 1e9),
              ('L_300', 'L_x 300 Ohm.um [nm]', 1e9),
              ('L_500', 'L_x 500 Ohm.um [nm]', 1e9),
              ('onoff_term', 'terminal on/off', 1.0),
              ('onoff_intr', 'intrinsic on/off', 1.0),
              ('Id_rf', 'I_d at the RF bias [A]', 1.0)]
    print('   %-32s %15s %15s %10s' % ('', '(4.29) C_q(dV)', 'C_q(E_F/e)',
                                       'change'))
    for key, lbl, sc_ in labels:
        a, b = base[key] * sc_, moved[key] * sc_
        print('   %-32s %15.6f %15.6f %9.4f %%'
              % (lbl, a, b, 100 * (b / a - 1.0)))
    atten = abs(100 * (moved['n_3V5'] / base['n_3V5'] - 1.0)) / abs(ma)
    check('the puddle floor attenuates the raw charge move before it reaches '
          'any committed number, which is why item 16 was survivable',
          atten < 1.0,
          'raw %+.4f %% becomes %+.4f %% in n(3.5 V), a factor %.3f.  R_ch at '
          'the Dirac point moves %+.6f %% because the floor supplies the whole '
          'density there.'
          % (ma, 100 * (moved['n_3V5'] / base['n_3V5'] - 1.0), atten,
             100 * (moved['R_dirac'] / base['R_dirac'] - 1.0)))
    survives = (1.0 < moved['onoff_term'] < 2.0
                and moved['L_300'] > 200e-9 and base['L_300'] > 200e-9)
    check('every qualitative claim of Sections 4.3-4.4 survives the '
          'correction',
          survives,
          'the on/off ratio stays between 1 and 2 (%.6f -> %.6f), the 200 nm '
          'crossover verdict is unchanged (L_x at 300 Ohm.um is %.2f nm -> '
          '%.2f nm, both above the 200 nm channel), and the resistance '
          'ordering by metal is unchanged.  The correction moves numbers and '
          'no verdicts.'
          % (base['onoff_term'], moved['onoff_term'], base['L_300'] * 1e9,
             moved['L_300'] * 1e9))

    print()
    print('SECTION 7.  RF FIGURES OF MERIT, AND THE TWO HALVES SEPARATED')
    print('-' * 78)
    pk_A = X.peak_rf()
    with on_consistent_argument(capacitance=False):
        pk_q = X.peak_rf()
    _ARRIVED['cq'] = 0
    with on_consistent_argument(capacitance=True):
        pk_f = X.peak_rf()
        hits_cq_rf = _ARRIVED['cq']
    check('ARRIVAL CONTROL on the CAPACITANCE patch, checked where C_gs is '
          'actually used: rf_small_signal_model.gate_capacitance() picked up '
          'the corrected C_q without being edited',
          hits_cq_rf > 0 and abs(pk_f['Cgs'] / pk_A['Cgs'] - 1.0) > 1e-12,
          '%d calls to the corrected C_q inside the RF sweep, and C_gs at the '
          'peak moved %+.4f %%.  The signature was kept identical '
          '(overdrive in, F/m^2 out) precisely so that this module would not '
          'have to be touched; the control is what confirms that worked.'
          % (hits_cq_rf, 100 * (pk_f['Cgs'] / pk_A['Cgs'] - 1.0)))
    print('  Per-um convention, V_ds = %.1f V, at the peak-f_T bias:' % VDS_RF)
    print('   %-20s %13s %13s %13s' % ('', 'C_q(dV)', 'charge only',
                                       'charge + C_gs'))
    for key, scale, unit in (('Vg', 1.0, 'V'), ('fT', 1e-9, 'GHz'),
                             ('fmax', 1e-9, 'GHz'), ('gm', 1e3, 'mS'),
                             ('Cgs', 1e15, 'fF')):
        print('   %-20s %13.4f %13.4f %13.4f'
              % ('%s [%s]' % (key, unit), pk_A[key] * scale,
                 pk_q[key] * scale, pk_f[key] * scale))
    d_q = pk_q['fT'] / pk_A['fT'] - 1.0
    d_f = pk_f['fT'] / pk_A['fT'] - 1.0
    d_c = d_f - d_q
    print()
    print('   f_T: charge alone %+.4f %%, the consistent C_gs adds %+.4f %%, '
          'total %+.4f %%' % (100 * d_q, 100 * d_c, 100 * d_f))
    check('the charge half and the capacitance half of item 16 have OPPOSITE '
          'signs in f_T, the same structure 2026-10-10 found for item 12',
          d_q * d_c < 0.0,
          'the two halves are %+.4f %% and %+.4f %% and leave %+.4f %%, which '
          'is %.1f %% of the charge half.  2026-10-10 recorded that structure '
          'as a property of the DEVICE TOPOLOGY rather than of any one '
          'parameter, on one instance.  It now has a second, independent '
          'instance, which is what that claim needed.'
          % (100 * d_q, 100 * d_c, 100 * d_f,
             100 * abs(d_f / d_q) if d_q else 0.0))
    r_A = pk_A['fmax'] / pk_A['fT']
    r_f = pk_f['fmax'] / pk_f['fT']
    print()
    print('   f_max/f_T = %.6f -> %.6f (%+.4f %%)'
          % (r_A, r_f, 100 * (r_f / r_A - 1.0)))
    check("Section 7.8.1a's structural verdict is untouched by item 16 -- the "
          'SIXTH lever tried against it',
          abs(r_f / r_A - 1.0) < 0.01 and r_f < 1.3,
          'f_max/f_T moves %+.4f %% and stays at %.4f, still below the '
          '1.3-1.4 band a perfect gate cannot reach.  Saturation, the '
          'perfect contact, the diffusion term, the exact charge relation, '
          "its capacitance, and now C_q's argument.  The binding constraint "
          'is still that transfer_characteristic() is a resistor with no '
          'output resistance.' % (100 * (r_f / r_A - 1.0), r_f))

    print()
    print('SECTION 9.  THE LITERATURE REGULARISES C_q WITH DISORDER, NOT WITH')
    print('            TEMPERATURE -- AND HERE THE TWO ARE THE SAME SIZE')
    print('-' * 78)
    print('  Live search today reached a SECOND, independent compact-model')
    print('  lineage: the Southampton bilayer-GFET model')
    print('  (eprints.soton.ac.uk/366757, GFetIEEEtran.pdf).  It confirms')
    print("  Section 1 from a different source than 2026-10-10's")
    print('  Pasadas/Jimenez review: its C_q is a function of the SURFACE')
    print('  POTENTIAL phi_s, never of the gate overdrive, and the layer')
    print('  potentials are a COUPLED SYSTEM solved together -- which is')
    print('  exactly the claim that item 16 is a self-consistency problem.')
    print()
    print('  But it regularises C_q at neutrality differently.  It has no')
    print('  thermal rounding at all: temperature never enters its C_q.  It')
    print('  floors the capacitance with the RESIDUAL DISORDER DENSITY n0,')
    print('    C_q,min = e^2*sqrt(n0)/(sqrt(pi)*hbar*v_F)')
    print('  i.e. the puddle density does the job this thesis gives to kT.')
    print()
    cq_floor = (E_CHARGE ** 2 * np.sqrt(gfet.n_puddle)
                / (np.sqrt(np.pi) * 1.054571817e-34 * gfet.v_F))
    cq_therm = float(cq_thermal(0.0))
    print('   this thesis, thermal:        C_q(0)   = %.6e F/m^2' % cq_therm)
    print('   literature form, at n_puddle: C_q,min = %.6e F/m^2' % cq_floor)
    print('   ratio = %.4f' % (cq_floor / cq_therm))
    check('the two regularisations are the SAME SIZE at this disorder level, '
          'so applying both -- as this thesis does -- overcounts C_q at the '
          'Dirac point by about a factor 2',
          0.5 < cq_floor / cq_therm < 2.0,
          'disorder floor %.6e against thermal %.6e F/m^2, ratio %.4f.  This '
          'model applies the thermal rounding INSIDE C_q and then adds '
          'n_puddle in quadrature to n, and the literature form shows those '
          'are two spellings of one physical effect at this n0.  Neither is '
          'wrong on its own; using both is a double count of the same '
          'regularisation, and it has never been costed.  NEW BACKLOG ITEM '
          '(Section 7.9 item 18).'
          % (cq_floor, cq_therm, cq_floor / cq_therm))
    n0_equal = (cq_therm * np.sqrt(np.pi) * 1.054571817e-34
                * gfet.v_F / E_CHARGE ** 2) ** 2
    T_equal = cq_floor / (CQ_PREFACTOR * np.log(4.0)) * E_CHARGE / KB
    print()
    print('   the two agree exactly at n0 = %.4e /m^2 (%.3f x n_puddle)'
          % (n0_equal, n0_equal / gfet.n_puddle))
    print('   or, holding n0 = n_puddle, at T = %.2f K (%.3f x 300 K)'
          % (T_equal, T_equal / gfet.T))
    check('and the coincidence is not a numerical accident of one parameter: '
          'the crossover density is within a factor 2 of n_puddle and the '
          'crossover temperature within a factor 2 of 300 K, so SiO2-supported '
          'graphene at room temperature sits where disorder and temperature '
          'regularise C_q equally',
          0.3 < n0_equal / gfet.n_puddle < 3.0
          and 0.3 < T_equal / gfet.T < 3.0,
          'n0* = %.4e /m^2 = %.3f n_puddle, and T* = %.2f K = %.3f x 300 K.  '
          'That is why neither convention is obviously wrong in the '
          'literature and why the double count here is easy to miss: the two '
          'terms it adds are the same order at the operating point this '
          'thesis studies.  It would NOT be true for suspended graphene '
          '(n0 two orders lower) or at 4 K, where the conventions diverge.'
          % (n0_equal, n0_equal / gfet.n_puddle, T_equal, T_equal / gfet.T))

    # ANNOTATION 2026-10-12, Section 7.9 item 18 COSTED AND CLOSED
    # (graphene_double_regularisation.py, 19/19; thesis Section 4.6.11).
    # Every number printed above stands and none of it is deleted, but TWO
    # of its readings do not survive being costed.
    #
    #  (a) 'applying both -- as this thesis does -- overcounts C_q at the
    #      Dirac point by about a factor 2' is FALSE for the committed
    #      model.  The two regularisers reach DISJOINT quantities: C_q's
    #      argument in (4.29) is dV and dV is never floored, so n_puddle
    #      cannot reach C_q; and n_eff = (C_ox*dV/e)*f(C_q) has dV as a
    #      FACTOR, so n_eff(0) = 0.0 bitwise under every C_q convention and
    #      the thermal rounding cannot reach the Dirac-point density.  In the
    #      one spelling where both DO reach C_q -- form B_T with the floor
    #      inside the argument -- the cost is +2.5262 %, not a factor 2,
    #      because regularisers entering through a monotone function of a
    #      quadrature sum compose MAX-LIKE (both/max = 1.025262).  The
    #      factor 2 was the RATIO OF THE TWO ALTERNATIVES, which prices
    #      choosing wrongly between them rather than applying both.
    #
    #  (b) the ratio 1.1509 is measured against a floor spelled with HALF
    #      this thesis's own degeneracy convention.  graphene's
    #      D(E) = 2|E|/(pi hbar^2 v_F^2) carries g = 4 and gives
    #      C_Q(n) = 2e^2 sqrt(n)/(sqrt(pi) hbar v_F) = 2e sqrt(n)/A_F, of
    #      which cq_floor above is EXACTLY half (0.000 ulp from 1/2 with one
    #      consistent hbar).  On this thesis's own convention the ratio is
    #      2.3018 and 'the same size' weakens to 'the same order'.  1.1509 is
    #      kept as printed because which spelling is right depends on what n0
    #      denotes in the source, which is a question for the literature.
    #
    #  (c) and cq_floor's hand-entered 1.054571817e-34 is the SAME literal
    #      A_FERMI carries, so the comparison above was internally
    #      consistent; mixing scipy's hbar in instead moves it by +6.13e-10
    #      relative and degrades the exact factor 2 in (b) to ~3e-10.

    print()
    print('SECTION 8.  FIGURE')
    print('-' * 78)
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))

    dVs = np.linspace(-2.0, 2.6, 600)
    ax[0].plot(dVs, np.abs(n_arg_dV(dVs)) * 1e-4, lw=2,
               label=r'committed: $C_q(\Delta V)$')
    ax[0].plot(dVs, np.abs(n_arg_EF_thermal(dVs)) * 1e-4, lw=2, ls='--',
               label=r'corrected: $C_q(E_F/e)$')
    ax[0].axhline(gfet.n_puddle * 1e-4, color='0.5', lw=1, ls=':',
                  label=r'$n_{\rm puddle}$')
    ax[0].set_xlabel(r'$\Delta V = V_g - V_{\rm Dirac}$  [V]')
    ax[0].set_ylabel(r'$|n|$  [cm$^{-2}$]')
    ax[0].set_yscale('log')
    ax[0].set_title('(a) raw induced charge')
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)

    Vs = np.logspace(-4, 1.2, 400)
    ax[1].loglog(Vs, cq_thermal(Vs), lw=2,
                 label=r'thermal, $\ln(2(1+\cosh\eta))$')
    ax[1].loglog(Vs, CQ_PREFACTOR * Vs, lw=1.5, ls='--',
                 label=r'$T=0$ limit, PREF$\cdot V$')
    ax[1].axhline(float(cq_thermal(0.0)), color='0.5', lw=1, ls=':',
                  label=r'$C_q(0)={\rm PREF}\cdot(kT/e)\ln 4$')
    ax[1].axvline(KB * gfet.T / E_CHARGE, color='C3', lw=1, ls='-.',
                  label=r'$kT/e$')
    ax[1].set_xlabel('argument  [V]')
    ax[1].set_ylabel(r'$C_q$  [F m$^{-2}$]')
    ax[1].set_title('(b) the thermal form and its $T=0$ limit')
    ax[1].legend(fontsize=8, loc='lower right')
    ax[1].grid(alpha=0.3, which='both')

    xx = np.arange(3)
    ax[2].bar(xx - 0.19, [abs(ma0), abs(mb0), abs(mboth0)], 0.38,
              label='vs $T=0$ (4.28)  [2026-10-10]')
    ax[2].bar(xx + 0.19, [abs(ma), abs(mb), abs(mboth)], 0.38,
              label='vs 300 K reference  [today]')
    ax[2].set_xticks(xx)
    ax[2].set_xticklabels(['argument\nonly', 'quadrature\nonly', 'both'],
                          fontsize=9)
    ax[2].set_ylabel(r'$|$move in raw charge$|$  [%]')
    ax[2].set_title('(c) the decomposition depends on the reference')
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=0.3, axis='y')

    fig.suptitle("Section 7.9 item 16: fixing $C_q$'s argument, and the "
                 'temperature the obvious fix deletes', fontsize=11)
    fig.tight_layout()
    fig.savefig('quantum_capacitance_argument.png', dpi=150)
    print('  wrote quantum_capacitance_argument.png')

    print()
    print('=' * 78)
    print('CHECKS: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return failed


if __name__ == '__main__':
    raise SystemExit(main())
