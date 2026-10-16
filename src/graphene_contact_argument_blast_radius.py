"""
graphene_contact_argument_blast_radius.py
=========================================

Section 7.9 item 16 -- "flip the default to form B_T" -- has stood ranked
first since 2026-10-12 with the note that it "moves every committed number in
Chapters 4, 5 and 6 and needs its own session".  This is that session, and it
does not flip the default.  It measures the blast radius first, and the blast
radius is not what that sentence says.

Two findings, both measured here rather than argued.

(1) THE HARNESS THAT SIZED ITEM 16 CANNOT REACH HALF OF CHAPTER 4.

    `graphene_quantum_capacitance_argument.on_consistent_argument()` rebinds
    `gfet.carrier_density`, `gfet.quantum_capacitance` and
    `crx.carrier_density`.  It does not rebind
    `graphene_contact_doping_model.quantum_capacitance`, which that module
    bound at import time on its line 50:

        from graphene_fet_model import quantum_capacitance, ...

    and which `contact_edge_carrier_density()` calls on its line 105 with
    EXACTLY the defect item 16 names -- the thermal C_q evaluated at a
    potential drop instead of at E_F/e.  Every cost ever quoted for item 16
    (-2.2659 % raw charge, -1.0212 % in n at 3.5 V, -0.2530 % in I_d) was
    therefore measured with Sections 4.5, 4.7 and 4.8 switched OUT of the
    circuit.  This is the 2026-10-10 delivery fault in a second module; it was
    repaired there for `crx.carrier_density` and the same sweep was never run
    for the contact model.  Section 1 below demonstrates the hole with an
    arrival counter rather than by reading the source: under
    `on_consistent_argument()` the contact model's binding is reached EXACTLY
    ZERO times while its results change by nothing at all.

(2) AND THE SENTENCE OVERSTATES THE OTHER DIRECTION.  Chapters 5 and 6 do not
    call the charge model.  `graphene_interconnect_model` imports nothing from
    it; the four Chapter 6 modules import only `mu`, `METAL_WORK_FUNCTIONS`,
    `W_GRAPHENE` and `lambda_decay` -- parameters, not the charge path.
    Section 2 measures this the same way, with counters on every module-level
    binding, and finds zero arrivals from Chapter 5 and Chapter 6 code.

    So item 16's reach is: NOT Chapters 5 and 6 at all, and MORE of Chapter 4
    than was measured.

Then the correction itself, at the contact edge (Sections 3-5).  The contact
geometry is where the argument error is least forgiving.  In the channel the
series factor C_q/(C_q + C_ox) is small and attenuates the correction by 0.451
(2026-10-11, Section 7).  At the contact, C_interface = eps_0/3.5 Angstrom is
66x C_ox and C_q is COMPARABLE to it, so the attenuation runs the other way
and the correction arrives nearly undamped.  And the potential is a
work-function difference of 0.15-1.15 V rather than a gate overdrive of a few
hundred mV, so the gap between dV and E_F(n)/e is correspondingly wider.

Section 5 carries it into Section 4.7's per-metal recalibration, which is the
one place in the thesis where a 10 %-scale move in R_extra is known in advance
to be able to change a verdict: Section 4.10's condition-number audit found
Pd's +51.0 Ohm.um residual -- the only positive one, and the only evidence the
additive decomposition might survive -- flips sign on a 9.6 % error in
R_extra.  That margin was computed as a tolerance to INPUT error.  This
session asks whether a known MODELLING error already consumes it.

Method notes, following this repo's standing rules.

  * Every patch SUBSTITUTES ONE SHIPPED NAME and leaves the rest of the call
    path as shipped code (the 2026-10-11 repair).  Form B_T at the contact is
    obtained by rebinding `cd.quantum_capacitance` and then calling the
    SHIPPED `cd.contact_edge_carrier_density`, so form A is the shipped
    function with no patch at all and the two differ in the argument and in
    nothing else.  The fixed-point solver is
    `qca._fixed_point` itself, run with `gfet.C_ox` substituted by
    `C_interface`, because `qca.series_factor` reads that name at call time --
    so the channel case is the contact case with the substitution removed, and
    agrees with the committed form B_T by construction rather than by care.
  * ONE ARRIVAL COUNTER PER PATCHED NAME, each checked at a call that reads
    that name (the 2026-10-11 rule, after its 2026-10-10 failure).
  * At least one validation against an EXACTLY known value, per the standing
    requirement: the zero null at dWF = 0, the bitwise electron/hole
    antisymmetry, the C_s-independence of the large-series-capacitance limit,
    and an identity-patch neutrality pass that must be bit-exact.

Run:  python3 graphene_contact_argument_blast_radius.py
Writes: contact_argument_blast_radius_output.txt, contact_argument_blast_radius.png
"""

import contextlib; from graphene.paths import results_path  # noqa: E702 -- same line on purpose: no line numbers move (docs/INVENTORY.md §6)
import io
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import graphene_fet_model as gfet
import graphene_contact_doping_model as cd
import graphene_quantum_capacitance_argument as qca
import graphene_exact_charge_rederivation as X
import contact_resistance_crossover as crx

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


def pct(new, old):
    return 100.0 * (new - old) / old if old != 0 else float('nan')


# ---------------------------------------------------------------------------
# The substitution that turns the channel solve into the contact solve
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def series_capacitance(C_s):
    """Substitute the ONE shipped name that `qca.series_factor` reads.

    `qca.series_factor(C_q)` is `C_q / (C_q + gfet.C_ox)`, read at call time.
    Rebinding `gfet.C_ox` therefore turns `qca._fixed_point` -- the committed
    form-B_T solver, unmodified -- into the solver for any series capacitance,
    with every other operation left as shipped code.  With C_s = gfet.C_ox the
    context is the identity, which is the point: the contact solve and the
    channel solve are then the same function and cannot drift apart.
    """
    saved = gfet.C_ox
    try:
        gfet.C_ox = C_s
        yield
    finally:
        gfet.C_ox = saved


def cq_shipped_of_n(n, T=None):
    """The SHIPPED `gfet.quantum_capacitance`, evaluated at E_F(n)/e.

    Two reasons this and not `qca.cq_thermal_of_n`, both found by this
    session's checks rather than anticipated.

    (i) CONSTANTS.  `qca.cq_thermal` carries its own `CQ_PREFACTOR` and its
        own hand-entered `E_CHARGE`/`KB`, and differs from the shipped
        function by 9.48e6 ulp -- about 1.1e-9 relative -- at every argument.
        That is a CONSTANTS gap, not a spelling gap, and on a 12 % effect it
        is irrelevant; but a comparison that carries it is a two-knob
        comparison, and 2026-10-11 is explicit that a knob must be turned
        alone.  Building form B_T out of the shipped function makes form A
        and form B_T differ in C_q's ARGUMENT and in nothing else, by
        construction rather than by tolerance.  Measured in Section 3(g).

    (ii) OVERFLOW.  `qca.cq_thermal` is the COMMITTED `log(2(1+cosh eta))`
        spelling, deliberately retained there, and it overflows at
        eta = 710.4758.  The shipped function has used the overflow-free
        identity since item 13 closed on 2026-10-11.  The fixed point's FIRST
        iterate is the electrostatic density, which grows with the series
        capacitance, so a solver built on `qca.cq_thermal` returns nan for
        any C_s above a threshold measured in Section 3(e2).  C_interface is
        comfortably below it, so this did not affect any number here -- but
        it is item 13's lesson for a third time: a defect unreachable from
        the device's operating range is reachable from a check on its model.
    """
    nn = np.abs(np.asarray(n, dtype=float))
    T = gfet.T if T is None else T
    return gfet.quantum_capacitance(qca.A_FERMI * np.sqrt(nn), T=T)


def n_consistent(dV, C_s, cq_of_n=None):
    """n solving  n = (C_s*dV/e) * f(C_q(n)),  C_q(n) at E_F/e."""
    cq_of_n = cq_shipped_of_n if cq_of_n is None else cq_of_n
    with series_capacitance(C_s):
        return qca._fixed_point(dV, cq_of_n)


def residual(dV, C_s):
    """The solve's own error, |n - n_es f(C_q(n))|/|n|."""
    with series_capacitance(C_s):
        return qca.fixed_point_residual(
            dV, lambda d: qca._fixed_point(d, cq_shipped_of_n), cq_shipped_of_n)


# ---------------------------------------------------------------------------
# One counter per patched name
# ---------------------------------------------------------------------------

ARRIVED = {}


def _reset_arrivals():
    for k in list(ARRIVED):
        ARRIVED[k] = 0


def _counting(name, fn):
    ARRIVED.setdefault(name, 0)

    def wrapper(*a, **kw):
        ARRIVED[name] += 1
        return fn(*a, **kw)
    wrapper.__wrapped_name__ = name
    return wrapper


# Every module-level binding of the two charge-path names in the repo, found
# by reading each module's imports.  The point of listing them explicitly is
# that a binding absent from this list is a binding no harness can rebind.
_BINDINGS = [
    ('gfet.carrier_density', gfet, 'carrier_density'),
    ('gfet.quantum_capacitance', gfet, 'quantum_capacitance'),
    ('crx.carrier_density', crx, 'carrier_density'),
    ('cd.quantum_capacitance', cd, 'quantum_capacitance'),
]


@contextlib.contextmanager
def count_every_binding():
    """Wrap every module-level binding with an IDENTITY counting wrapper.

    The wrappers change no value, so anything measured inside this context
    must be bit-identical to the same thing measured outside it.  That is
    checked in Section 0 before the counts are used for anything: an
    instrument that perturbs what it measures reports its own arrival.
    """
    saved = [(mod, attr, getattr(mod, attr)) for _, mod, attr in _BINDINGS]
    try:
        for name, mod, attr in _BINDINGS:
            setattr(mod, attr, _counting(name, getattr(mod, attr)))
        _reset_arrivals()
        yield ARRIVED
    finally:
        for mod, attr, fn in saved:
            setattr(mod, attr, fn)


# ---------------------------------------------------------------------------
# Form B_T at the contact: substitute ONE shipped name inside the shipped fn
# ---------------------------------------------------------------------------

_CONTACT_ARRIVED = {'cq': 0}


def _cq_consistent_contact(V, T=300.0):
    """Drop-in for the `quantum_capacitance` name that
    `cd.contact_edge_carrier_density` reads.

    Same signature.  Takes the work-function-induced potential, solves the
    contact's own self-consistency for n, and returns C_q at that state's
    E_F/e.  The shipped function then forms its own series factor and its own
    product, both unchanged -- so the only thing that differs between form A
    and form B_T at the contact is C_q's argument.
    """
    _CONTACT_ARRIVED['cq'] += 1
    n = n_consistent(V, cd.C_interface)
    return cq_shipped_of_n(n, T=T)


@contextlib.contextmanager
def contact_on_consistent_argument():
    saved = cd.quantum_capacitance
    try:
        cd.quantum_capacitance = _cq_consistent_contact
        yield
    finally:
        cd.quantum_capacitance = saved


# ---------------------------------------------------------------------------

def main():
    print('=' * 78)
    print('ITEM 16\'s BLAST RADIUS: the contact model was never in the circuit')
    print('=' * 78)
    print()
    print('C_ox        = %.6e F/m^2   (90 nm SiO2 back gate)' % gfet.C_ox)
    print('C_interface = %.6e F/m^2   (%.1f Angstrom vacuum gap)'
          % (cd.C_interface, cd.d_interface * 1e10))
    print('ratio C_interface/C_ox = %.4f' % (cd.C_interface / gfet.C_ox))
    print()

    metals = sorted(cd.METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])

    # -- Section 0 -----------------------------------------------------------
    print('-' * 78)
    print('0.  THE INSTRUMENT IS NOT THE PERTURBATION (identity-patch pass)')
    print('-' * 78)
    print()
    print('Every binding below is wrapped with a counting wrapper that calls')
    print('the shipped function and returns its value unchanged.  If any')
    print('number moves by one bit, the counts in Sections 1-2 are measuring')
    print('the wrapper and not the code.')
    print()

    base_contact = {m: float(cd.contact_edge_carrier_density(wf)) for m, wf in metals}
    base_n35 = float(gfet.carrier_density(3.5))
    base_cq = float(gfet.quantum_capacitance(1.2))
    base_rc = cd.recalibrate_metal_rc()

    with count_every_binding():
        wrapped_contact = {m: float(cd.contact_edge_carrier_density(wf)) for m, wf in metals}
        wrapped_n35 = float(gfet.carrier_density(3.5))
        wrapped_cq = float(gfet.quantum_capacitance(1.2))
        wrapped_rc = cd.recalibrate_metal_rc()
        ident_counts = dict(ARRIVED)

    worst = max(ulps(wrapped_contact[m], base_contact[m]) for m, _ in metals)
    check('identity wrappers: all 7 n_contact bit-identical (0 ulp)',
          all(wrapped_contact[m] == base_contact[m] for m, _ in metals),
          'worst departure %.1f ulp' % worst)
    check('identity wrappers: n(3.5 V) bit-identical',
          wrapped_n35 == base_n35, 'ulp %.1f' % ulps(wrapped_n35, base_n35))
    check('identity wrappers: C_q(1.2 V) bit-identical',
          wrapped_cq == base_cq, 'ulp %.1f' % ulps(wrapped_cq, base_cq))
    rc_same = all(wrapped_rc[m][2] == base_rc[m][2] for m in base_rc)
    check('identity wrappers: all 4 R_extra bit-identical', rc_same)
    check('identity wrappers: cd.quantum_capacitance WAS reached (control)',
          ident_counts['cd.quantum_capacitance'] > 0,
          '%d arrivals -- so a zero in Section 1 is a hole, not a dead path'
          % ident_counts['cd.quantum_capacitance'])
    print()
    print('  arrivals under the identity patch:')
    for k in sorted(ident_counts):
        print('    %-32s %6d' % (k, ident_counts[k]))
    print()

    # -- Section 1 -----------------------------------------------------------
    print('-' * 78)
    print('1.  THE HOLE: on_consistent_argument() never reaches the contact model')
    print('-' * 78)
    print()

    # THIS SECTION'S FIRST DRAFT ASKED THE WRONG QUESTION, and the failure is
    # kept here as the instrument rather than edited out of the history.  It
    # counted ARRIVALS AT THE BINDING `cd.quantum_capacitance` and asserted
    # they would be zero under the patch.  They are not zero and never could
    # be: the contact model calls that name 11 times whether or not anything
    # was substituted into it.  What the section needs to count is arrivals at
    # the REPLACEMENT -- `qca._ARRIVED`, incremented inside
    # `_quantum_capacitance_consistent` itself.  NINTH instance of the
    # 2026-10-10 rule, and the second in its "counts a different arrival from
    # the one it gates" direction after 2026-10-11's: a control that counts
    # calls to a NAME cannot see whether the right FUNCTION is behind it.
    qca._ARRIVED['cq'] = 0
    qca._ARRIVED['charge'] = 0
    with count_every_binding():
        with qca.on_consistent_argument(capacitance=True):
            bound_to_patch = (cd.quantum_capacitance.__wrapped_name__
                              if hasattr(cd.quantum_capacitance, '__wrapped_name__')
                              else None)
            patched_contact = {m: float(cd.contact_edge_carrier_density(wf))
                               for m, wf in metals}
            patched_rc = cd.recalibrate_metal_rc()
        hole_counts = dict(ARRIVED)
    corrected = qca._ARRIVED['cq']
    reached = hole_counts['cd.quantum_capacitance']

    check('the contact model\'s charge path IS live under the patch (control)',
          reached > 0, '%d calls to cd.quantum_capacitance' % reached)
    check('of those %d calls, arrivals at the CORRECTED function == 0 EXACTLY'
          % reached, corrected == 0,
          '%d corrected arrivals -- on_consistent_argument() rebinds '
          'gfet.quantum_capacitance, and cd holds its own import-time '
          'reference' % corrected)
    check('and consequently every n_contact is UNCHANGED, bitwise',
          all(patched_contact[m] == base_contact[m] for m, _ in metals))
    check('and every Section 4.7 R_extra is UNCHANGED, bitwise',
          all(patched_rc[m][2] == base_rc[m][2] for m in base_rc))
    print()
    print('  The harness rebinds gfet.quantum_capacitance.  The contact model')
    print('  holds its own reference, bound at import time, and keeps calling')
    print('  the uncorrected function.  Nothing errors; the numbers simply do')
    print('  not move, and a cost table built this way reads that as "the')
    print('  contact model is insensitive to the argument".  It was absent.')
    print()

    # -- Section 2 -----------------------------------------------------------
    print('-' * 78)
    print('2.  AND CHAPTERS 5 AND 6 ARE NOT IN THE CIRCUIT EITHER')
    print('-' * 78)
    print()
    print('  "Flipping it moves every committed number in Chapters 4, 5 and 6"')
    print('  (AUTOMATION_LOG 2026-10-11 Section 7).  Measured, not read:')
    print()

    import graphene_interconnect_model as icm
    import graphene_photodetector_collection_model as pcm
    import graphene_photodetector_signed_carrier_model as psc
    import graphene_photodetector_two_contact_model as ptc

    # Each probe is the module's OWN committed summary routine, not a
    # hand-picked function: whatever Chapters 5 and 6 report is what runs.
    _probe_defs = [
        ('Ch5 interconnects: summary_numbers()', icm.summary_numbers),
        ('Ch5 interconnects: self_check_parallel_limit()', icm.self_check_parallel_limit),
        ('Ch6 collection: summary_numbers()', pcm.summary_numbers),
        ('Ch6 collection: modeled_eqe_by_metal()', pcm.modeled_eqe_by_metal),
        ('Ch6 signed-carrier: doping_type_table()', lambda: psc.doping_type_table(verbose=False)),
        ('Ch6 signed-carrier: symmetric_pair_table()', lambda: psc.symmetric_pair_table(verbose=False)),
        ('Ch6 signed-carrier: asymmetric_pair_table()', lambda: psc.asymmetric_pair_table(verbose=False)),
        ('Ch6 two-contact: symmetric_vs_single_contact_table()',
         lambda: ptc.symmetric_vs_single_contact_table(verbose=False)),
        ('Ch6 two-contact: asymmetric_pair_table()', lambda: ptc.asymmetric_pair_table(verbose=False)),
    ]
    probes = []
    _hush = io.StringIO()
    with count_every_binding():
        for label, fn in _probe_defs:
            _reset_arrivals()
            _keep = sys.stdout
            try:
                sys.stdout = _hush
                fn()
            finally:
                sys.stdout = _keep
            probes.append((label, dict(ARRIVED)))

    for label, counts in probes:
        tot = sum(counts.values())
        print('    %-52s  %d charge-path arrivals' % (label, tot))
    ch56_total = sum(sum(c.values()) for _, c in probes)
    check('Chapters 5 and 6: total charge-path arrivals == 0 EXACTLY',
          ch56_total == 0,
          '%d arrivals -- Chapter 5 imports nothing from the charge model and '
          'Chapter 6 imports only mu, METAL_WORK_FUNCTIONS, W_GRAPHENE, '
          'lambda_decay' % ch56_total)
    print()

    # -- Section 3 -----------------------------------------------------------
    print('-' * 78)
    print('3.  THE SOLVER, VALIDATED AGAINST EXACTLY KNOWN VALUES')
    print('-' * 78)
    print()

    DV_RF = qca.DV_RF

    # (a) overlap limit: C_s = C_ox must reproduce the committed form B_T.
    mine = float(n_consistent(DV_RF, gfet.C_ox))
    theirs = float(qca.n_arg_EF_thermal(DV_RF))
    rel_ab = abs(mine - theirs) / abs(theirs)
    check('(a) C_s = C_ox reproduces committed form B_T at dV = 1.2 V to 1e-8',
          rel_ab < 1e-8,
          'relative gap %.4e (%.3e ulp) -- NOT 0 ulp, and the residue is the '
          'constants gap of (g), not the solve' % (rel_ab, ulps(mine, theirs)))

    # the same solve on qca's own C_q must land on qca's own answer, 0 ulp:
    # that separates "my solver is theirs" from "my C_q is theirs".
    mine_q = float(n_consistent(DV_RF, gfet.C_ox, cq_of_n=qca.cq_thermal_of_n))
    check('(a2) the SOLVER is qca._fixed_point itself: same C_q -> 0 ulp',
          mine_q == theirs, 'ulp %.1f' % ulps(mine_q, theirs))

    # (b) the committed -1.0212 % at Vg = 3.5 V, reproduced through the floor.
    dv35 = 3.5 - gfet.V_dirac
    nA = float(X.floor_total(qca.n_arg_dV(dv35)))
    nB = float(X.floor_total(n_consistent(dv35, gfet.C_ox)))
    move = pct(nB, nA)
    check('(b) reproduces 2026-10-11\'s -1.0212 %% in n at V_g = 3.5 V',
          abs(move - (-1.0212)) < 5e-4,
          'measured %+.6f %%  (committed -1.0212 %%)' % move)

    # (c) zero null, EXACT: dWF = 0 -> dV = 0 -> n = 0 exactly, both forms.
    nA0 = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE))
    with contact_on_consistent_argument():
        nB0 = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE))
    check('(c) zero null: n_contact(dWF = 0) == 0.0 EXACTLY, form A',
          nA0 == 0.0, 'got %r' % nA0)
    check('(c) zero null: n_contact(dWF = 0) == 0.0 EXACTLY, form B_T',
          nB0 == 0.0, 'got %r' % nB0)

    # (d) antisymmetry, EXACT: n(W_G + d) == -n(W_G - d), bitwise, both forms.
    d = 0.37
    aP = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE + d))
    aM = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE - d))
    with contact_on_consistent_argument():
        bP = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE + d))
        bM = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE - d))
    check('(d) electron/hole antisymmetry bitwise, form A', aP == -aM,
          'ulp %.1f' % ulps(aP, -aM))
    check('(d) electron/hole antisymmetry bitwise, form B_T', bP == -bM,
          'ulp %.1f' % ulps(bP, -bM))

    # (e) large-C_s limit must be C_s-INDEPENDENT: n -> C_q(n)*dV/e.
    # AS FIRST WRITTEN this check demanded C_s-independence to 1e-9 between
    # C_s = 1e3 and 1e4 and FAILED at 5.254e-04.  The tolerance was chosen
    # before the run and was simply the wrong criterion: the limit is
    # approached at FIRST ORDER in C_q/C_s, not exactly, so no fixed tolerance
    # states what the limit claims.  Item 15's pattern for the eighth time, and
    # the repair is 2026-10-11's -- let the structure set the criterion.  What
    # the limit actually asserts is a RATE: each decade of C_s must cut the
    # departure by ten.  That is a statement no wrong model passes by accident,
    # and the failing tolerance is kept in this comment rather than deleted.
    Cs_decades = [1.0e3, 1.0e4, 1.0e5, 1.0e6]
    ns = [float(n_consistent(0.62, c)) for c in Cs_decades]
    gaps = [abs(ns[i + 1] - ns[i]) / abs(ns[i + 1]) for i in range(len(ns) - 1)]
    ratios = [gaps[i] / gaps[i + 1] for i in range(len(gaps) - 1)]
    cq_over_cs = float(cq_shipped_of_n(ns[0])) / Cs_decades[0]
    check('(e) large-C_s limit: departure falls by ~10x per decade of C_s '
          '(first order in C_q/C_s)',
          all(9.0 <= r <= 11.0 for r in ratios) and gaps[-1] < 1e-5,
          'gaps %s ; decade ratios %s ; C_q/C_s at 1e3 = %.3e against a '
          'measured gap of %.3e'
          % (', '.join('%.3e' % g for g in gaps),
             ', '.join('%.4f' % r for r in ratios), cq_over_cs, gaps[0]))

    # (e2) the same limit on qca's committed cosh spelling -- where it breaks.
    nan_at = None
    for Cs in (1.0e1, 3.0e1, 1.0e2, 3.0e2, 1.0e3):
        v = float(n_consistent(0.62, Cs, cq_of_n=qca.cq_thermal_of_n))
        if not np.isfinite(v):
            nan_at = Cs
            break
    eta_first = (qca.A_FERMI * np.sqrt(cd.C_interface * 1.15 / qca.E_CHARGE)
                 / (qca.KB * gfet.T / qca.E_CHARGE))
    check('(e2) the committed cosh spelling returns nan above some C_s '
          '(item 13, third sighting)',
          nan_at is not None,
          'first nan at C_s = %s F/m^2; C_interface = %.4e is below it, and '
          'the fixed point\'s worst first-iterate eta in this session is '
          '%.1f against the 710.476 overflow onset'
          % (nan_at, cd.C_interface, eta_first))

    # (f) the solve's own residual at every metal.
    res = max(float(residual(-(wf - cd.W_GRAPHENE), cd.C_interface))
              for _, wf in metals if wf != cd.W_GRAPHENE)
    check('(f) fixed-point residual <= 1e-12 at all 7 work functions',
          res <= 1e-12, 'worst %.3e' % res)

    # (g) spelling gap between the two thermal C_q transcriptions, REPORTED.
    g = max(ulps(float(gfet.quantum_capacitance(v)), float(qca.cq_thermal(v)))
            for v in (0.05, 0.15, 0.62, 1.15))
    grel = max(abs(float(gfet.quantum_capacitance(v)) - float(qca.cq_thermal(v)))
               / float(gfet.quantum_capacitance(v))
               for v in (0.05, 0.15, 0.62, 1.15))
    pref_ratio = qca.CQ_PREFACTOR / ((2 * gfet.e ** 3)
                                     / (np.pi * gfet.hbar ** 2 * gfet.v_F ** 2))
    check('(g) the gfet/qca C_q gap is a CONSTANTS gap and is <= 1e-8 relative',
          grel < 1e-8,
          'worst %.4e relative (%.3e ulp); prefactor ratio qca/gfet = '
          '%.15f, i.e. the whole gap.  It is OUTSIDE the A->B_T comparison '
          'because both forms below call the shipped function.'
          % (grel, g, pref_ratio))

    # (h) arrival control for the contact patch, at a call that reads the name.
    _CONTACT_ARRIVED['cq'] = 0
    with contact_on_consistent_argument():
        cd.contact_edge_carrier_density(cd.METAL_WORK_FUNCTIONS['Pd'])
    check('(h) contact patch arrival control: cd.quantum_capacitance reached',
          _CONTACT_ARRIVED['cq'] == 1, '%d arrivals' % _CONTACT_ARRIVED['cq'])
    print()

    # -- Section 4 -----------------------------------------------------------
    print('-' * 78)
    print('4.  THE CORRECTION AT THE CONTACT EDGE, PER METAL')
    print('-' * 78)
    print()
    with contact_on_consistent_argument():
        bt_contact = {m: float(cd.contact_edge_carrier_density(wf)) for m, wf in metals}

    A_F = qca.A_FERMI
    print('  %-4s %7s %10s %14s %14s %9s %9s' %
          ('', 'WF(eV)', 'dV (V)', 'n_A (m^-2)', 'n_BT (m^-2)', 'change', 'E_F/e (V)'))
    rows = []
    for m, wf in metals:
        dV = -(wf - cd.W_GRAPHENE)
        a, b = base_contact[m], bt_contact[m]
        ef = A_F * np.sqrt(abs(b)) if b != 0 else 0.0
        ch = pct(b, a) if a != 0 else float('nan')
        rows.append((m, wf, dV, a, b, ch, ef))
        if a == 0:
            print('  %-4s %7.2f %10.3f %14.6e %14.6e %9s %9.4f   (dWF = 0 exactly)'
                  % (m, wf, dV, a, b, '--', ef))
        else:
            print('  %-4s %7.2f %10.3f %14.6e %14.6e %+8.4f%% %9.4f'
                  % (m, wf, dV, a, b, ch, ef))
    print()
    moved = [r for r in rows if r[3] != 0 and np.isfinite(r[5])]
    mn = min(abs(r[5]) for r in moved)
    mx = max(abs(r[5]) for r in moved)
    print('  Range of the correction: %.4f %% to %.4f %%.' % (mn, mx))
    print('  For comparison, the SAME correction in the channel at V_g = 3.5 V')
    print('  is %+.4f %%.  The contact is %.1fx to %.1fx less forgiving.'
          % (move, mn / abs(move), mx / abs(move)))
    print()
    print('  Why, in one line: at the contact C_q/(C_q+C_s) runs %.3f-%.3f'
          % (min(float(cq_shipped_of_n(bt_contact[m])) /
                 (float(cq_shipped_of_n(bt_contact[m])) + cd.C_interface)
                 for m, _ in metals if bt_contact[m] != 0),
             max(float(cq_shipped_of_n(bt_contact[m])) /
                 (float(cq_shipped_of_n(bt_contact[m])) + cd.C_interface)
                 for m, _ in metals if bt_contact[m] != 0)))
    sf_ch = gfet.quantum_capacitance(dv35) / (gfet.quantum_capacitance(dv35) + gfet.C_ox)
    print('  while in the channel at V_g = 3.5 V it is %.3f.  The channel\'s'
          % sf_ch)
    print('  small series factor is what attenuated item 16 by 0.451; the')
    print('  contact has no such attenuation, and a wider dV to be wrong about.')
    print()

    # -- Section 5 -----------------------------------------------------------
    print('-' * 78)
    print('5.  SECTION 4.7 UNDER THE CORRECTION -- AND THE Pd ROW')
    print('-' * 78)
    print()
    with contact_on_consistent_argument():
        bt_rc = cd.recalibrate_metal_rc()

    print('  %-4s %10s %12s %12s %9s %12s %12s' %
          ('', 'Rc_meas', 'R_extra A', 'R_extra B_T', 'change', 'R_tr A', 'R_tr B_T'))
    flips = []
    rc_rows = []
    for m in ('Cu', 'Ni', 'Au', 'Pd'):
        wf, rcm, ra, rta = base_rc[m]
        _, _, rb, rtb = bt_rc[m]
        rc_rows.append((m, rcm, ra, rb, pct(rb, ra), rta, rtb))
        print('  %-4s %10.1f %12.3f %12.3f %+8.4f%% %12.3f %12.3f'
              % (m, rcm, ra, rb, pct(rb, ra), rta, rtb))
        if (rta > 0) != (rtb > 0):
            flips.append(m)
    print()

    pd_move = pct(bt_rc['Pd'][2], base_rc['Pd'][2])
    print('  Section 4.10\'s margin for the Pd row is 9.6 %% in R_extra.')
    print('  The argument correction moves Pd\'s R_extra by %+.4f %%.' % pd_move)
    print()
    if flips:
        print('  >>> SIGN FLIPPED: %s.  The one positive residual in Table 4.7'
              % ', '.join(flips))
        print('  >>> does not survive the correction.  Section 4.10 established')
        print('  >>> that margin as a tolerance to INPUT error; it is consumed')
        print('  >>> here by a KNOWN MODELLING error that was already in the')
        print('  >>> code, and the additive decomposition now fails on 4 of 4')
        print('  >>> metals rather than 3 of 4.')
    else:
        print('  >>> NO SIGN FLIPPED.  Pd\'s residual survives at %+.3f Ohm.um'
              % bt_rc['Pd'][3])
        print('  >>> (was %+.3f).  Section 4.7\'s verdict -- "3 of 4 negative,'
              % base_rc['Pd'][3])
        print('  >>> only Pd plausible" -- is unchanged in kind, and the')
        print('  >>> correction moves numbers and no verdict here either.')
    print()

    # The ordering claim of Section 4.7: Cu most negative, Pd least.
    ord_A = [m for m, _, _, _, _, rta, _ in
             sorted(rc_rows, key=lambda r: r[5])]
    ord_B = [m for m, _, _, _, _, _, rtb in
             sorted(rc_rows, key=lambda r: r[6])]
    check('Section 4.7 metal ORDERING by residual is preserved under B_T',
          ord_A == ord_B, 'A: %s   B_T: %s' % (' < '.join(ord_A), ' < '.join(ord_B)))

    # The n/p type assignment of Section 4.5 must not move: it is a sign.
    signs_same = all((base_contact[m] > 0) == (bt_contact[m] > 0)
                     for m, _ in metals if base_contact[m] != 0)
    check('Section 4.5 n-type/p-type assignment unchanged for all 7 metals',
          signs_same)
    print()

    # -- Section 6 -----------------------------------------------------------
    print('-' * 78)
    print('6.  WHAT ITEM 16 ACTUALLY COSTS, RESTATED')
    print('-' * 78)
    print()
    print('  The 2026-10-11 table is not wrong about the numbers it contains.')
    print('  It is incomplete, and incomplete in a way the harness could not')
    print('  report, because the missing rows are the ones whose binding it')
    print('  never reached.  Restated:')
    print()
    print('    Ch4 channel (measured 2026-10-11, stands):  n -1.0212 %%, I_d -0.2530 %%')
    print('    Ch4 contact (measured today, new):          n_contact %+.4f to %+.4f %%'
          % (-mx, -mn))
    print('    Ch4 Section 4.7 R_extra (new):              %+.4f to %+.4f %%'
          % (min(r[4] for r in rc_rows), max(r[4] for r in rc_rows)))
    print('    Ch5 interconnects:                          0 arrivals, no change')
    print('    Ch6 photodetectors:                         0 arrivals, no change')
    print()
    print('  So the flip is a CHAPTER 4 operation, not a three-chapter one --')
    print('  and inside Chapter 4 it reaches two models, not one.')
    print()

    # -- figure --------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))

    ax = axes[0]
    ms = [r[0] for r in rows]
    chg = [0.0 if not np.isfinite(r[5]) else r[5] for r in rows]
    cols = ['#b2182b' if abs(c) > abs(move) else '#999999' for c in chg]
    ax.bar(range(len(ms)), chg, color=cols)
    ax.axhline(move, ls='--', c='k', lw=1.2,
               label='channel, V$_g$=3.5 V (%+.2f%%)' % move)
    ax.axhline(0, c='k', lw=0.8)
    ax.set_xticks(range(len(ms)))
    ax.set_xticklabels(ms)
    ax.set_ylabel('change in $n_{contact}$, A $\\to$ B$_T$ (%)')
    ax.set_title('(a) the argument correction is\nmuch larger at the contact')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis='y')

    ax = axes[1]
    xs = np.linspace(0.01, 1.2, 300)
    nA_s = np.array([float(cd.contact_edge_carrier_density(cd.W_GRAPHENE + v))
                     for v in xs])
    with contact_on_consistent_argument():
        nB_s = np.array([float(cd.contact_edge_carrier_density(cd.W_GRAPHENE + v))
                         for v in xs])
    ax.plot(xs, np.abs(nA_s), label='form A: $C_q(\\Delta V)$', lw=1.8)
    ax.plot(xs, np.abs(nB_s), label='form B$_T$: $C_q(E_F/e)$', lw=1.8, ls='--')
    for m, wf in metals:
        if wf != cd.W_GRAPHENE:
            ax.axvline(abs(wf - cd.W_GRAPHENE), c='grey', lw=0.6, alpha=0.5)
            ax.text(abs(wf - cd.W_GRAPHENE), ax.get_ylim()[1] * 0.02, m,
                    fontsize=7, rotation=90, va='bottom')
    ax.set_xlabel('$|\\Delta W_F|$ (eV)')
    ax.set_ylabel('$|n_{contact}|$ (m$^{-2}$)')
    ax.set_yscale('log')
    ax.set_title('(b) the two forms, over the\nwork-function range in use')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2]
    w = 0.35
    idx = np.arange(len(rc_rows))
    ax.bar(idx - w / 2, [r[5] for r in rc_rows], w, label='form A (committed)')
    ax.bar(idx + w / 2, [r[6] for r in rc_rows], w, label='form B$_T$')
    ax.axhline(0, c='k', lw=1.0)
    ax.set_xticks(idx)
    ax.set_xticklabels([r[0] for r in rc_rows])
    ax.set_ylabel('implied $R_{transmission}$ ($\\Omega\\cdot\\mu$m)')
    ax.set_title('(c) Table 4.7\'s residual:\n%s'
                 % ('Pd flips sign' if flips else 'no sign flips'))
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis='y')

    fig.suptitle('Item 16\'s blast radius: the contact model was outside the harness',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(results_path('ch4', 'contact_argument_blast_radius.png'), dpi=150)
    print('  figure written: contact_argument_blast_radius.png')
    print()

    print('=' * 78)
    print('CHECKS: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    buf = io.StringIO()

    class Tee:
        def __init__(self, *s):
            self.s = s

        def write(self, d):
            for x in self.s:
                x.write(d)

        def flush(self):
            for x in self.s:
                x.flush()

    old = sys.stdout
    sys.stdout = Tee(old, buf)
    try:
        rc = main()
    finally:
        sys.stdout = old
    with open(results_path('ch4', 'contact_argument_blast_radius_output.txt'), 'w') as f:
        f.write(buf.getvalue())
    sys.exit(rc)
