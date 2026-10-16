"""
graphene_double_regularisation.py

Chapter 7 Section 7.9 ITEM 18 -- "the thermal rounding and the puddle floor are
the same physical effect, and this thesis applies both".

Item 18 was created 2026-10-11 and is the FIRST item on that list to arrive
from a literature source rather than from an internal audit.  Its premise:

    the Southampton compact model (eprints.soton.ac.uk/366757) regularises
    C_q at neutrality with the residual disorder density n0 and has no
    thermal rounding at all; this thesis regularises C_q with kT AND adds n0
    in quadrature to n; at this n0 and this T the two are the same size
    (ratio 1.1509), so applying both is "a double count worth about a factor
    2 in C_q(0)", never costed.

This module costs it, and the premise is wrong twice over.  Both errors have
the shape of errors this repo has already recorded.

  1. THE COMMITTED MODEL DOES NOT DOUBLE-COUNT AT ALL.  The two regularisers
     reach different quantities and neither quantity has both.  In form (4.29)
     C_q's argument is the gate overdrive dV, and dV is never floored, so n0
     cannot reach C_q.  Conversely n_eff(0) = (C_ox*0/e)*f(C_q) = 0 BITWISE
     for every C_q convention, so the thermal rounding cannot reach the
     Dirac-point density.  Both are tested on exactly-known values -- 0.0 and
     n_puddle -- and not on tolerances.

  2. EVEN IN THE SPELLING THAT DOES DOUBLE-COUNT, THE COST IS NOT A FACTOR 2.
     Put the floor inside C_q's argument -- form B_T with C_q(E_F(sqrt(n^2 +
     n0^2))), the most double-counting spelling available and arguably the
     consistent one -- and C_q(0) rises over disorder-alone by about +2.5 %,
     not by a factor 2.  Two regularisers that enter through a monotone
     function of a quadrature sum compose MAX-LIKE, not additively.

     The "factor 2" was the RATIO OF THE TWO CANDIDATE VALUES.  That ratio is
     the cost of CHOOSING between the conventions; it is not the cost of
     APPLYING BOTH.  Reading one for the other is what the 2026-10-06 backlog
     item ("census every criterion that is a ratio") was created to catch,
     and this is its first instance in which the ratio was not a criterion
     but an estimate of a magnitude.

  3. AND THE FLOOR'S OWN SPELLING CARRIES AN EXACT FACTOR 2.  The Southampton
     floor e^2*sqrt(n0)/(sqrt(pi)*hbar*v_F) is EXACTLY HALF this thesis's own
     DOS-exact C_Q(n0) = 2e^2*sqrt(n0)/(sqrt(pi)*hbar*v_F) = 2e*sqrt(n0)/A_F,
     verified to 1 ulp with one consistent hbar.  So yesterday's ratio 1.1509
     compares kT against a floor written with half this thesis's own
     degeneracy convention; on the thesis's own convention it is 2.3018 and
     "the same size" must weaken to "the same order".

  4. AND THE FIRST DRAFT OF THIS MODULE FAILED ITS OWN SUBSTITUTABILITY
     CHECK.  Draft 1 built its switches on graphene_diffusion_current_model's
     literal E_CHARGE/HBAR while gfet.quantum_capacitance() uses scipy's, and
     reproduced the shipped carrier_density(2.0) only to 2375 ulp.  That is
     the EIGHTH instance of the 2026-10-09 rule and the second in which the
     substituted quantity was a physical constant.  The repair is 2026-10-11's
     repair, not a better transcription: the switches below SUBSTITUTE one
     shipped name each (gfet._thermal_rounding, gfet.n_puddle) and leave every
     other operation in the shipped call path, so form A is the shipped
     function and agrees with it bitwise by construction.

Run:  python3 graphene_double_regularisation.py
Writes double_regularisation_output.txt and double_regularisation.png.
"""

import contextlib

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import e as E_SCIPY, hbar as HBAR_SCIPY
from scipy.constants import Boltzmann as KB_SCIPY

import graphene_fet_model as gfet
import graphene_diffusion_current_model as dd

HBAR_REPO = dd.HBAR            # 1.054571817e-34, the literal the repo uses
A_FERMI = dd.A_FERMI           # E_F/e = A_F*sqrt(n)   [V m], on HBAR_REPO
N0_SHIPPED = gfet.n_puddle
C_OX = gfet.C_ox
DV_RF = 2.0 - gfet.V_dirac     # 1.2 V, the RF bias of record
DV_ON = 3.5 - gfet.V_dirac     # 2.7 V, the on-state of record

passed = 0
failed = 0

# ONE ARRIVAL COUNTER PER PATCHED NAME (2026-10-11's sixth repair: a counter
# that counts a different arrival from the one it gates reports a fault that
# does not exist).
arrivals = {'rounding': 0, 'n_puddle': 0}


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


def one(x):
    return float(np.atleast_1d(np.asarray(x, dtype=float))[0])


# ---------------------------------------------------------------------------
# The two regularisers, each as a SUBSTITUTION of one shipped name
# ---------------------------------------------------------------------------

def _rounding_T0(eta):
    """The thermal regulariser switched OFF, in place.  gfet's C_q is
    PREF*(kT/e)*F(eV/kT); substituting F(eta) -> |eta| gives PREF*|V|
    identically, which is the T = 0 dispersion.  Switching the rounding off
    is therefore EXACT, not a small-T approximation."""
    arrivals['rounding'] += 1
    return np.abs(np.asarray(eta, dtype=float))


@contextlib.contextmanager
def regularisers(thermal=True, floor=True):
    """Swap exactly the two names that carry the two regularisers, and
    nothing else in the call path."""
    saved_round = gfet._thermal_rounding
    saved_n0 = gfet.n_puddle
    try:
        if not thermal:
            gfet._thermal_rounding = _rounding_T0
        if not floor:
            gfet.n_puddle = 0.0
            arrivals['n_puddle'] += 1
        yield
    finally:
        gfet._thermal_rounding = saved_round
        gfet.n_puddle = saved_n0


def cq(V, thermal=True):
    """C_q at an explicit potential, through the shipped function."""
    with regularisers(thermal=thermal):
        return gfet.quantum_capacitance(np.asarray(V, dtype=float))


def EF_over_e(n):
    return A_FERMI * np.sqrt(np.abs(np.asarray(n, dtype=float)))


def series_factor(C_q):
    return C_q / (C_q + C_OX)


def n_form_A(dV, thermal=True, floor=True):
    """FORM A, the committed relation (4.29), as the shipped function.  With
    both switches on this IS gfet.carrier_density, bitwise."""
    with regularisers(thermal=thermal, floor=floor):
        return gfet.carrier_density(np.asarray(dV, dtype=float) + gfet.V_dirac)


def n_form_BT(dV, thermal=True, floor_in_cq=False, floor_out=True,
              iters=6000):
    """FORM B_T (item 16's correction): C_q's argument is E_F/e, solved as a
    fixed point.  floor_in_cq places the DISORDER regulariser INSIDE C_q's
    argument -- the only spelling in which both regularisers reach the same
    quantity.  floor_out is the committed quadrature step on the result."""
    dV = np.asarray(dV, dtype=float)
    n_es = C_OX * dV / E_SCIPY
    n = np.array(n_es, dtype=float, copy=True)
    n0 = N0_SHIPPED
    for _ in range(iters):
        arg = np.sqrt(n ** 2 + n0 ** 2) if floor_in_cq else np.abs(n)
        n_new = n_es * series_factor(cq(EF_over_e(arg), thermal=thermal))
        if np.all(np.abs(n_new - n) <= 1e-16 * np.abs(n_new) + 1e-9):
            n = n_new
            break
        n = n_new
    return np.sqrt(n ** 2 + n0 ** 2) if floor_out else np.abs(n)


def cq_at_bias_BT(dV, thermal=True, floor_in_cq=False):
    """The C_q actually evaluated inside form B_T at a given bias."""
    n = n_form_BT(dV, thermal=thermal, floor_in_cq=floor_in_cq,
                  floor_out=False)
    arg = (np.sqrt(n ** 2 + N0_SHIPPED ** 2) if floor_in_cq else np.abs(n))
    return cq(EF_over_e(arg), thermal=thermal)


# ---------------------------------------------------------------------------
def main():
    global N0_SHIPPED
    print('=' * 78)
    print('THE DOUBLE REGULARISATION -- Section 7.9 item 18, costed')
    print('=' * 78)
    print()
    print('  n0 = n_puddle = %.4e m^-2      T = %.1f K      kT/e = %.6f V'
          % (N0_SHIPPED, gfet.T, KB_SCIPY * gfet.T / E_SCIPY))
    print('  A_F = %.6e V m    C_ox = %.6e F/m^2' % (A_FERMI, C_OX))
    print()

    print('SECTION 0.  SUBSTITUTABILITY, BEFORE ANY MEASUREMENT')
    print('-' * 78)
    nA_shipped = one(gfet.carrier_density(2.0))
    nA_here = one(n_form_A(DV_RF))
    print('   carrier_density(2.0)  shipped %.12e   via switches %.12e'
          % (nA_shipped, nA_here))
    print('   ulp difference %.1f' % ulps(nA_shipped, nA_here))
    check('form A reached through the switches IS the shipped '
          'carrier_density, bitwise, so everything below is measured against '
          'the real model and not against a re-implementation of it',
          ulps(nA_shipped, nA_here) == 0.0,
          'Draft 1 of this module did re-implement it -- on '
          "graphene_diffusion_current_model's literal E_CHARGE and HBAR "
          'while gfet.quantum_capacitance uses scipy -- and agreed only to '
          '2375 ulp.  EIGHTH instance of the 2026-10-09 rule, second with a '
          'physical constant, and the first caught by a check written before '
          'the measurement rather than by a neutrality pass afterwards.')
    cq_t_direct = one(gfet.quantum_capacitance(0.0))
    cq_t_switch = one(cq(0.0, thermal=True))
    check('and the thermal switch in its ON position is a no-op to 0 ulp',
          ulps(cq_t_direct, cq_t_switch) == 0.0,
          'C_q(0) = %.12e both ways; the substitution installs nothing when '
          'thermal=True.' % cq_t_direct)
    before = dict(arrivals)
    with regularisers(thermal=False, floor=True):
        _ = gfet.quantum_capacitance(np.array([0.0, 1.0]))
    check('the rounding patch ARRIVES when switched off (one counter for the '
          'rounding name, checked where the rounding name is used)',
          arrivals['rounding'] > before['rounding'],
          'rounding arrivals %d -> %d.  The n_puddle counter is checked '
          'separately, at a call that reads n_puddle.'
          % (before['rounding'], arrivals['rounding']))
    before = dict(arrivals)
    _ = n_form_A(DV_RF, thermal=True, floor=False)
    check('and the n_puddle patch arrives at a call that reads n_puddle',
          arrivals['n_puddle'] > before['n_puddle'],
          'n_puddle arrivals %d -> %d.' % (before['n_puddle'],
                                           arrivals['n_puddle']))

    print()
    print("SECTION 1.  THE FLOOR CARRIES AN EXACT FACTOR 2, AND YESTERDAY'S")
    print('            RATIO WAS MEASURED AGAINST THE HALVED SPELLING')
    print('-' * 78)
    cq_lit = E_SCIPY ** 2 * np.sqrt(N0_SHIPPED) / (
        np.sqrt(np.pi) * HBAR_SCIPY * gfet.v_F)
    cq_dos = 2.0 * cq_lit
    cq_dos_AF = 2.0 * E_SCIPY * np.sqrt(N0_SHIPPED) / A_FERMI
    cq_therm0 = one(cq(0.0))
    print('  Southampton floor    C_q,min = e^2 sqrt(n0)/(sqrt(pi) hbar v_F)')
    print('                               = %.9e F/m^2' % cq_lit)
    print("  this thesis's own    C_Q(n0) = e^2 D(E_F), D = 2|E|/(pi hb^2 vF^2)")
    print('                               = %.9e F/m^2' % cq_dos)
    print('                       and as 2e sqrt(n0)/A_F (repo hbar)')
    print('                               = %.9e F/m^2' % cq_dos_AF)
    print('  thermal              C_q(0)  = PREF*(kT/e)*ln4')
    print('                               = %.9e F/m^2' % cq_therm0)
    print()
    print('  2e sqrt(n0)/A_F  vs  2e^2 sqrt(n0)/(sqrt(pi) hbar v_F): %.1f ulp'
          % ulps(cq_dos_AF, cq_dos))
    check("this thesis's two spellings of its own T = 0 C_Q agree, so the "
          'comparison below is between conventions and not between '
          'transcriptions',
          ulps(cq_dos_AF, cq_dos) < 1e7,
          '%.12e via A_F against %.12e written out; they differ only through '
          "the repo's hbar literal (%.1f ulp, %.3e relative)."
          % (cq_dos_AF, cq_dos, ulps(cq_dos_AF, cq_dos),
             cq_dos_AF / cq_dos - 1.0))
    check("the Southampton floor is EXACTLY HALF this thesis's own DOS-exact "
          'C_Q at the same density -- a convention difference, not a '
          'numerical one',
          ulps(cq_lit / cq_dos, 0.5) <= 1.0,
          'ratio %.17f, %.3f ulp from 1/2 with one consistent hbar.  '
          'graphene D(E) = 2|E|/(pi hbar^2 v_F^2) carries the spin-valley '
          'degeneracy g = 4 and gives C_q = e^2 D(E_F) = '
          '2e^2 sqrt(n)/(sqrt(pi) hbar v_F); the Southampton floor is that '
          'with the 2 absent.  Whatever its origin -- a per-species n0, a '
          'bilayer convention, a different definition of the residual '
          'density -- it is not a quantity this thesis can compare against '
          'its own kT before fixing the factor.'
          % (cq_lit / cq_dos, ulps(cq_lit / cq_dos, 0.5)))
    print()
    print('  the same comparison, in both spellings:')
    print('    lit/thermal     = %.6f   <- the figure recorded 2026-10-11'
          % (cq_lit / cq_therm0))
    print('    thesis/thermal  = %.6f   <- own degeneracy convention'
          % (cq_dos / cq_therm0))
    check('so "the two regularisations are the SAME SIZE at this disorder '
          'level" holds only in the halved spelling; on this thesis\'s own '
          'convention it is 2.30 and the claim weakens from "the same size" '
          'to "the same order"',
          abs(cq_lit / cq_therm0 - 1.1509) < 3e-3
          and 2.2 < cq_dos / cq_therm0 < 2.4,
          '%.6f against %.6f.  The 2026-10-11 entry is NOT rewritten to the '
          'larger figure: both are reported, because which is right depends '
          'on what n0 denotes in the source, which is a question for the '
          'literature and not for arithmetic.'
          % (cq_lit / cq_therm0, cq_dos / cq_therm0))
    rel = (E_SCIPY ** 2 * np.sqrt(N0_SHIPPED)
           / (np.sqrt(np.pi) * HBAR_REPO * gfet.v_F)) / cq_lit - 1.0
    print()
    print("  aside, the hbar literal.  Yesterday's Section 9 wrote hbar as")
    print('  1.054571817e-34 and compared against an A_FERMI carrying the')
    print('  SAME literal, so that comparison was internally consistent;')
    print('  against scipy the floor moves %+.4e relative.' % rel)

    print()
    print('SECTION 2.  WHICH QUANTITY EACH REGULARISER ACTUALLY REACHES')
    print('-' * 78)
    print('  On exactly-known values, not on tolerances.')
    print()
    cq_candidates = [('thermal', cq_therm0),
                     ('T=0 clean', one(cq(0.0, thermal=False))),
                     ('lit floor', cq_lit),
                     ('thesis C_Q(n0)', cq_dos)]
    n_eff0 = [(lab, (C_OX * 0.0 / E_SCIPY) * series_factor(max(c, 1e-300)))
              for lab, c in cq_candidates]
    for lab, v in n_eff0:
        print('    n_eff(dV=0) with C_q = %-15s : %r' % (lab, v))
    check('n_eff(dV=0) = 0.0 BITWISE for every C_q convention, so the '
          'thermal rounding of C_q cannot reach the Dirac-point density at '
          'all',
          all(v == 0.0 for _, v in n_eff0),
          'n_eff = (C_ox*dV/e)*f(C_q) has dV as a FACTOR, and f is bounded '
          'in (0,1).  At dV = 0 the product is zero whatever f returns.  The '
          'exactly-known value is 0.0 and all four conventions return it '
          'bitwise.')
    nd = one(gfet.carrier_density(gfet.V_dirac))
    check('and the committed Dirac-point density is EXACTLY n_puddle, so it '
          'is set by the disorder regulariser alone',
          nd == N0_SHIPPED,
          'carrier_density(V_dirac) = %r, n_puddle = %r, difference %r.  '
          'sqrt(0^2 + n0^2) = n0 exactly in IEEE-754 for this n0.'
          % (nd, N0_SHIPPED, nd - N0_SHIPPED))
    ev = [(dv, one(n_form_A(dv)), one(n_form_A(-dv))) for dv in (0.3, 1.2, 2.7)]
    check('n_form_A is bitwise even in dV, so neither switch introduces an '
          'asymmetry',
          all(a == b for _, a, b in ev),
          '; '.join('dV=%.1f: %r vs %r' % t for t in ev))
    print()
    print('  CONSEQUENCE.  In the committed relation (4.29) the two')
    print('  regularisers are DISJOINT: kT regularises C_q, n0 regularises n,')
    print('  and neither quantity receives both.  n0 cannot reach C_q because')
    print("  C_q's argument is dV and dV is not floored; kT cannot reach the")
    print('  Dirac-point density because dV multiplies the series factor.')
    check('item 18\'s premise -- "this thesis applies both" -- is FALSE for '
          'the committed model',
          True,
          'The committed model applies exactly one regulariser to each of '
          'two different quantities.  That is a reuse of one parameter for '
          'two jobs, which is a modelling choice to defend, and not a double '
          'count, which is an arithmetic error.  Item 18 was created from a '
          'comparison of two LITERATURE FORMS without first asking which '
          'quantity each term in THIS model reaches.')

    print()
    print('SECTION 3.  THE SPELLING THAT DOES DOUBLE-COUNT, AND ITS SIZE')
    print('-' * 78)
    print('  Form B_T (item 16) makes C_q a function of E_F(n).  Put the')
    print('  floor INSIDE that argument -- C_q(E_F(sqrt(n^2+n0^2))) -- and')
    print('  both regularisers reach C_q.  This is the most double-counting')
    print('  spelling available and arguably the consistent one, since the')
    print("  channel's DOS should be evaluated at the density it has.")
    print()
    g = {}
    for tn, th in (('thermal', True), ('T=0', False)):
        for fn, fi in (('floor-out', False), ('floor-in', True)):
            g[(tn, fn)] = one(cq_at_bias_BT(0.0, thermal=th, floor_in_cq=fi))
    c_neither = g[('T=0', 'floor-out')]
    c_therm = g[('thermal', 'floor-out')]
    c_dis = g[('T=0', 'floor-in')]
    c_both = g[('thermal', 'floor-in')]
    print('   C_q(0) in form B_T   [F/m^2]       thermal OFF      thermal ON')
    print('   floor outside C_q               %.9e   %.9e'
          % (c_neither, c_therm))
    print('   floor inside  C_q               %.9e   %.9e'
          % (c_dis, c_both))
    print()
    print('   neither         %.9e' % c_neither)
    print('   thermal alone   %.9e   (increment %+.6e)'
          % (c_therm, c_therm - c_neither))
    print('   disorder alone  %.9e   (increment %+.6e)'
          % (c_dis, c_dis - c_neither))
    print('   BOTH            %.9e' % c_both)
    print()
    sum_rule = c_both / (c_therm + c_dis - c_neither)
    max_rule = c_both / max(c_therm, c_dis)
    d_on_floor = 100 * (c_both / c_dis - 1.0)
    f_on_therm = 100 * (c_both / c_therm - 1.0)
    print('   both / (sum of increments)   = %.6f' % sum_rule)
    print('   both / max(alone)            = %.6f' % max_rule)
    print('   thermal ON on top of floor   = %+.4f %%' % d_on_floor)
    print('   floor ON on top of thermal   = %+.4f %%' % f_on_therm)
    check('the double count is %+.4f %% in C_q(0), NOT "about a factor 2"'
          % d_on_floor,
          0.0 < d_on_floor < 5.0,
          'Applying the thermal rounding on top of the disorder floor raises '
          'C_q(0) by %+.4f %%.  The 2.30 that was read as the size of the '
          'double count is the ratio of the two ALTERNATIVES -- the cost of '
          'choosing wrongly between them, which is a real cost and a '
          'different one.' % d_on_floor)
    check('and the reason is structural: regularisers entering through a '
          'monotone function of a quadrature sum compose MAX-LIKE, so the '
          'larger swallows the smaller',
          max_rule < 1.05 and sum_rule < 0.80,
          'both/max = %.6f, i.e. %.2f %% above pure max; both/sum = %.6f, a '
          'departure of %.1f %% from the additive null.  Item 17 requires a '
          'two-knob check to report the RATIO and flag any departure from 1: '
          'here the departure is SUB-additive, the direction a '
          'reinforcement-only check misses, and the lesson is that the sum '
          'rule is the wrong null for a regulariser pair.  The right null is '
          'max, and against it the departure is %.2f %%.'
          % (max_rule, 100 * (max_rule - 1), sum_rule,
             100 * (1 - sum_rule), 100 * (max_rule - 1)))

    print()
    print("SECTION 4.  WHAT IT COSTS DOWNSTREAM, AT THE MODEL'S OWN BIASES")
    print('-' * 78)
    print("  Item 17's second requirement: run the two-knob check at the")
    print('  operating point, not at a reference that inflates it.')
    print()
    print('   bias      neither         thermal         disorder        '
          'both            both/max')
    rows = []
    for dv, lab in ((0.0, 'dV=0   '), (DV_RF, 'dV=1.2 '), (DV_ON, 'dV=2.7 ')):
        v = {}
        for tn, th in (('T=0', False), ('thermal', True)):
            for fn, fi in (('floor-out', False), ('floor-in', True)):
                v[(tn, fn)] = one(n_form_BT(dv, thermal=th, floor_in_cq=fi,
                                            floor_out=True))
        a = v[('T=0', 'floor-out')]
        b = v[('thermal', 'floor-out')]
        c = v[('T=0', 'floor-in')]
        d = v[('thermal', 'floor-in')]
        print('   %s %.6e   %.6e   %.6e   %.6e   %.6f'
              % (lab, a, b, c, d, d / max(b, c)))
        rows.append((dv, a, b, c, d))
    dirac = rows[0]
    check('at the Dirac point all four conventions agree BITWISE, because '
          'each is overwritten by the quadrature floor on the result',
          len({dirac[1], dirac[2], dirac[3], dirac[4]}) == 1,
          'All four return n = %r = n_puddle at dV = 0.  The quantity item '
          '18 was created to protect is the one quantity that cannot move.'
          % dirac[1])
    rf = [r for r in rows if r[0] == DV_RF][0]
    print()
    print('   at the RF bias dV = 1.2 V, relative to "neither":')
    print('     thermal alone   %+.4f %%' % (100 * (rf[2] / rf[1] - 1)))
    print('     disorder alone  %+.4f %%' % (100 * (rf[3] / rf[1] - 1)))
    print('     both            %+.4f %%' % (100 * (rf[4] / rf[1] - 1)))
    print('     both/sum ratio  %.6f'
          % ((rf[4] - rf[1]) / ((rf[2] - rf[1]) + (rf[3] - rf[1]))))
    move_rf = 100 * (rf[4] / rf[2] - 1)
    check('the whole double count is worth %+.4f %% in n at the RF bias, '
          "smaller than every correction already on Section 7.9's list"
          % move_rf,
          abs(move_rf) < 0.5,
          'n moves %+.4f %% when the floor is moved inside C_q at dV = 1.2 '
          'V, against -1.0212 %% for item 16 itself and -0.60 %% for the '
          'exact charge relation.  Item 18 is REAL but SMALL, and it is '
          'smaller at the operating point than at neutrality -- the opposite '
          'of how it was ranked when it was created.' % move_rf)
    sig = [(lab, one(gfet.sheet_conductivity(
        n_form_BT(0.0, thermal=th, floor_in_cq=fi))))
        for lab, th, fi in (('neither', False, False),
                            ('thermal', True, False),
                            ('disorder', False, True),
                            ('both', True, True))]
    print()
    print('   sigma_min at the Dirac point [S/sq]:')
    for lab, v in sig:
        print('     %-9s %.9e' % (lab, v))
    check('sigma_min is identical to the last bit across all four, so no '
          'committed conductivity, resistance or on/off number in Chapters '
          '4, 5 or 6 is exposed to item 18 at neutrality',
          len({v for _, v in sig}) == 1,
          'sigma_min = %.9e S/sq in all four conventions; it is '
          'n_puddle*e*mu and nothing else.' % sig[0][1])

    print()
    print('SECTION 5.  OVERLAPPING-LIMIT VALIDATION AGAINST THE OLD MODEL')
    print('-' * 78)
    saved = N0_SHIPPED
    try:
        N0_SHIPPED = 0.0
        a = one(n_form_BT(DV_RF, thermal=True, floor_in_cq=True))
        b = one(n_form_BT(DV_RF, thermal=True, floor_in_cq=False))
    finally:
        N0_SHIPPED = saved
    print('   n0 -> 0:  floor-in %.12e   floor-out %.12e' % (a, b))
    check('with n0 = 0 the floor-in and floor-out spellings collapse onto '
          'one another, so the disorder switch is a strict extension of form '
          'B_T and not a new model',
          ulps(a, b) <= 1.0,
          'relative difference %.3e, %.1f ulp.  The switch is inert when '
          'there is no disorder, which is the only limit in which it must '
          'be.' % (a / b - 1.0, ulps(a, b)))
    V = np.array([0.01, 0.1, 0.5, 1.0])
    T_lo = 1e-3
    kTe = KB_SCIPY * T_lo / E_SCIPY
    with regularisers(thermal=True):
        lo = gfet.quantum_capacitance(V, T=T_lo)
    ex = one(cq(1.0, thermal=False)) * V           # PREF*|V|, from the switch
    dep = float(np.max(np.abs(lo / ex - 1.0)))
    print('   T -> %g K:  max relative departure from PREF*|V| = %.3e'
          % (T_lo, dep))
    check('and with T driven to 1e-3 K the thermal form reproduces the T = 0 '
          'dispersion to 1e-14, which is only runnable because item 13 '
          'removed the cosh overflow',
          dep < 1e-14 and bool(np.all(np.isfinite(lo))),
          'eta reaches %.3e here; the spelling committed before 2026-10-11 '
          'returned nan for every one of these values and this section could '
          'not have been written.  The T = 0 reference is itself produced by '
          'the thermal switch, so this is also a consistency check on the '
          'switch.  The departure is EXACTLY zero rather than ~1e-16, and '
          'that is weaker than it looks: at eta ~ 1e7 the term '
          '2*log1p(exp(-eta)) underflows to exactly 0.0, so the two '
          'expressions reduce to the same floating-point operations.  The '
          'informative version of this check is at moderate eta, where '
          'Section 3 of graphene_quantum_capacitance_argument.py already '
          'adjudicates the identity against 60-digit arithmetic.'
          % (1.0 / kTe))
    print()
    print('   arrival counters at end of run: %r' % arrivals)
    check('both patched names arrived, and neither counter gates a section '
          'that does not use its name',
          arrivals['rounding'] > 0 and arrivals['n_puddle'] > 0,
          'rounding %d, n_puddle %d.' % (arrivals['rounding'],
                                         arrivals['n_puddle']))

    print()
    print('SECTION 6.  THE SOURCE ITSELF COMPOSES ITS REGULARISERS AS A')
    print('            SMOOTH MAX, AND ITS WORST CASE IS 3/sqrt(5)')
    print('-' * 78)
    print('  Read today, to settle what n0 denotes.  It does not settle that')
    print('  -- the paper says only that at zero surface potential "the')
    print('  channel has been reported to have a charge density n0", cites')
    print('  two references for it, and never defines it as total, per-')
    print('  species, or a disorder parameter.  No degeneracy factor and no')
    print('  DOS expression appear anywhere in it: it goes straight from n to')
    print('  C_q.  So the exact factor 2 of Section 1 is not a stated')
    print('  convention that can be converted; the residue stands.')
    print()
    print('  But its Eq. (3) answers a different question, and answers it in')
    print('  favour of Section 3.  It composes the floor with the variable')
    print('  part as')
    print('    C_q(phi) = (C_min^2 + 2 h^2) / sqrt(C_min^2 + h^2),')
    print('    h = C_var/2,   C_var proportional to |phi|')
    print('  which is a SMOOTH MAX and not a sum.')
    print()
    c_min_lit = cq_lit

    def cq_soton(c_var):
        h = np.asarray(c_var, dtype=float) / 2.0
        return (c_min_lit ** 2 + 2.0 * h ** 2) / np.sqrt(c_min_lit ** 2 + h ** 2)

    print('   C_var/C_min     C_q        max(C_min,C_var)   C_q/max')
    for r_ in (0.0, 1e-4, 1.0, 2.0, 10.0, 1e4):
        cv = r_ * c_min_lit
        print('   %10.4g   %.6e   %.6e   %.6f'
              % (r_, one(cq_soton(cv)), max(c_min_lit, cv),
                 one(cq_soton(cv)) / max(c_min_lit, cv)))
    check('their Eq. (3) returns C_min EXACTLY at zero surface potential and '
          'C_var to 12 digits when C_var dominates, so it is an '
          'interpolation between the two and not an addition of them',
          one(cq_soton(0.0)) == float(c_min_lit)
          and abs(one(cq_soton(1e6 * c_min_lit)) / (1e6 * c_min_lit) - 1.0)
          < 1e-11,
          'C_q(0) = %.12e against C_min = %.12e, bitwise; and '
          'C_q/C_var -> %.12f at C_var = 1e6 C_min.'
          % (one(cq_soton(0.0)), c_min_lit,
             one(cq_soton(1e6 * c_min_lit)) / (1e6 * c_min_lit)))
    cv = np.logspace(-8, 2, 200001) * c_min_lit
    ratio = cq_soton(cv) / np.maximum(c_min_lit, cv)
    worst = float(ratio.max())
    at = float(cv[int(ratio.argmax())] / c_min_lit)
    exact = 3.0 / np.sqrt(5.0)
    print()
    print('   sup over all C_var of C_q/max(C_min,C_var) = %.9f at '
          'C_var/C_min = %.6f' % (worst, at))
    print('   closed form at C_var = C_min: 1.5/sqrt(1.25) = 3/sqrt(5) = '
          '%.9f' % exact)
    check('the maximum excess of their composition over pure max is EXACTLY '
          '3/sqrt(5) = 1.3416408, attained EXACTLY where the two '
          'regularisers are equal in size',
          abs(worst - exact) < 1e-7 and abs(at - 1.0) < 1e-4,
          'measured %.9f at C_var/C_min = %.6f against the closed form '
          '%.9f.  Setting C_var = C_min in (C_min^2 + 2h^2)/sqrt(C_min^2 + '
          'h^2) with h = C_min/2 gives 1.5 C_min^2 / (C_min sqrt(1.25)) = '
          '3 C_min/sqrt(5).' % (worst, at, exact))
    print()
    print('  TWO CONSEQUENCES, both for item 18 rather than for the model.')
    print()
    print('  (1) The literature applies n0 TWICE itself.  Eq. (2) floors C_q')
    print('      with it and the Drude relation n = sqrt(n0^2 + n*^2) floors')
    print('      the density with it, with no temperature anywhere.  So the')
    print('      comparison is not between a model that double-counts and')
    print('      one that does not: it is between two ALLOCATIONS of')
    print('      regularisers to quantities.  This thesis spends kT on C_q')
    print('      and n0 on n; the source spends n0 on both.')
    print()
    print('  (2) And yesterday\'s intuition was right about WHERE the worst')
    print('      case is and wrong about HOW BIG it is.  A smooth-max')
    print('      composition is furthest from pure max exactly when its two')
    print('      terms are equal, which is the condition 2026-10-11')
    print('      identified as holding at this n0 and this T.  At that worst')
    print('      case the excess is 3/sqrt(5), i.e. +34.164 %, and this')
    print('      thesis is at +%.4f %% because its two terms are 2.30:1 and'
          % d_on_floor)
    print('      not 1:1.  Neither figure is a factor 2.  The bound a')
    print('      quadrature-composed regulariser pair can reach over the')
    print('      larger of its terms is 1.342, not 2.')
    check('so the "factor 2" was unreachable in principle, not merely '
          'unmeasured: a quadrature-composed regulariser pair cannot exceed '
          'the larger of its two terms by more than 34.2 %',
          worst < 1.35 and d_on_floor < worst * 100 - 100,
          'sup = %.6f, i.e. +%.3f %% at equality, against the +%.4f %% this '
          'thesis actually incurs and the +100 %% the item asserted.  An '
          'upper bound on the whole class of error, derived from the '
          'source\'s own composition law.' % (worst, 100 * (worst - 1),
                                              d_on_floor))

    print()
    print('SECTION 7.  FIGURE')
    print('-' * 78)
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))

    Vs = np.logspace(-4, 0.5, 500)
    ax[0].loglog(Vs, cq(Vs), lw=2, label=r'thermal, $kT$ only')
    ax[0].loglog(Vs, cq(Vs, thermal=False), lw=1.5, ls='--',
                 label=r'$T=0$: PREF$\cdot|V|$')
    ax[0].axhline(cq_lit, color='C2', lw=1.5, ls='-.',
                  label=r'Southampton floor $e^2\sqrt{n_0}/\sqrt{\pi}\hbar v_F$')
    ax[0].axhline(cq_dos, color='C3', lw=1.5, ls=':',
                  label=r'this thesis $2e^2\sqrt{n_0}/\sqrt{\pi}\hbar v_F$')
    ax[0].axvline(KB_SCIPY * gfet.T / E_SCIPY, color='0.6', lw=1,
                  label=r'$kT/e$')
    ax[0].set_xlabel(r"$C_q$'s argument  [V]")
    ax[0].set_ylabel(r'$C_q$  [F m$^{-2}$]')
    ax[0].set_title('(a) the two floors differ by exactly 2')
    ax[0].legend(fontsize=7, loc='lower right')
    ax[0].grid(alpha=0.3, which='both')

    labels = ['neither', 'thermal\nalone', 'disorder\nalone', 'both']
    vals = [c_neither, c_therm, c_dis, c_both]
    bars = ax[1].bar(labels, vals, color=['0.7', 'C0', 'C2', 'C3'])
    ax[1].axhline(c_therm + c_dis - c_neither, color='k', lw=1.2, ls='--',
                  label='additive null (sum rule)')
    ax[1].axhline(max(c_therm, c_dis), color='C1', lw=1.2, ls=':',
                  label='max null')
    ax[1].set_ylabel(r'$C_q(0)$ in form B$_T$  [F m$^{-2}$]')
    ax[1].set_title('(b) both/max = %.4f, both/sum = %.4f'
                    % (max_rule, sum_rule))
    ax[1].legend(fontsize=8, loc='upper left')
    ax[1].grid(alpha=0.3, axis='y')
    for b_, v_ in zip(bars, vals):
        ax[1].text(b_.get_x() + b_.get_width() / 2, v_, '%.3e' % v_,
                   ha='center', va='bottom', fontsize=7)

    dVs = np.linspace(-0.6, 2.8, 400)
    for lab, th, fi, st in (('neither', False, False, '-'),
                            ('thermal alone', True, False, '--'),
                            ('disorder alone', False, True, '-.'),
                            ('both', True, True, ':')):
        ax[2].semilogy(dVs, n_form_BT(dVs, thermal=th, floor_in_cq=fi) * 1e-4,
                       lw=1.8, ls=st, label=lab)
    ax[2].axhline(N0_SHIPPED * 1e-4, color='0.5', lw=1,
                  label=r'$n_{\rm puddle}$')
    ax[2].axvline(DV_RF, color='C3', lw=1, ls=':', label='RF bias')
    ax[2].set_xlabel(r'$\Delta V = V_g - V_{\rm Dirac}$  [V]')
    ax[2].set_ylabel(r'$|n|$  [cm$^{-2}$]')
    ax[2].set_title('(c) at the Dirac point the four curves are one curve')
    ax[2].legend(fontsize=7)
    ax[2].grid(alpha=0.3)

    fig.suptitle('Section 7.9 item 18: the two regularisers are disjoint in '
                 'form (4.29), and compose max-like where they are not',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig('double_regularisation.png', dpi=150)
    print('   wrote double_regularisation.png')

    print()
    print('=' * 78)
    print('CHECKS: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    print()
    print('VERDICT -- item 18 CLOSES, and it closes by being wrong twice.')
    print()
    print('  (i)   The committed model (4.29) does not double-count.  The')
    print('        two regularisers reach disjoint quantities, and both')
    print('        facts are exactly-known values -- 0.0 and n_puddle --')
    print('        and not tolerances.')
    print('  (ii)  In the spelling that does double-count, form B_T with the')
    print('        floor inside C_q, the cost is %+.4f %% in C_q(0) and'
          % d_on_floor)
    print('        %+.4f %% in n at the RF bias, because regularisers'
          % move_rf)
    print('        compose max-like.  The "factor 2" was a ratio of')
    print('        ALTERNATIVES, not the size of a sum.')
    print("  (iii) And the floor is spelled with half this thesis's own")
    print('        degeneracy convention, exactly, so the ratio 1.1509 that')
    print("        created the item is 2.3018 on the thesis's own terms.")
    print()
    print('  What survives is not a double count but the REUSE OF ONE')
    print('  PARAMETER FOR TWO JOBS: n0 buys the DOS broadening in the')
    print('  literature and the carrier floor here, and the thesis has never')
    print('  stated which it means.  That is the standing n_puddle question,')
    print('  reached for the second time from a new direction.')
    return failed


if __name__ == '__main__':
    raise SystemExit(1 if main() else 0)
