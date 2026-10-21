"""
graphene_default_flip.py
========================

Section 7.9 item 16 -- "flip the default to form B_T" -- has stood ranked
first since 2026-10-12.  Its preconditions were discharged on 2026-10-13
(the cost in the channel), 2026-10-14 (the identity class, so the set of
names to move is computed rather than listed) and 2026-10-17 (the contact's
own form B_T, and the measurement that every cost ever quoted for item 16 had
been taken with the contact model switched out of the circuit).  The log entry
of 2026-10-20 recorded that "the flip and nothing else is left".

That was wrong, and this module is the reason.

THE TWO HALVES OF ITEM 16 DO NOT COMPOSE.

`graphene_quantum_capacitance_argument.on_consistent_argument()` installs form
B_T in the channel.  `graphene_contact_argument_blast_radius.
contact_on_consistent_argument()` installs form B_T at the contact.  Each is
committed, each is validated, and each is correct on its own.  Running both at
once does NOT give form B_T everywhere: it moves `n_contact` by **-51.4 % to
-61.9 %** away from the committed contact-only form B_T, which is between
four and five times the entire effect item 16 exists to install
(-4.7696 % to -12.4964 %).

The cause is a single attribute read.  The contact replacement is deliberately
built on the SHIPPED `gfet.quantum_capacitance` -- see
`blast.cq_shipped_of_n`, whose docstring gives two good reasons for that
choice (a 9.48e6 ulp constants gap in the alternative, and an overflow in the
committed `log(2(1+cosh))` spelling).  It reads that name by attribute lookup
at call time:

    return gfet.quantum_capacitance(qca.A_FERMI * np.sqrt(nn), T=T)

and the CHANNEL half of the flip rebinds exactly that name.  So with both
halves installed, the contact's form-B_T solve calls the channel's form-B_T
solve as its inner C_q, and the contact's self-consistency is solved against a
C_q that has itself already been solved -- once, wrongly, against the wrong
series capacitance.

This is the 2026-10-14 delivery-fault family seen from a new side, and it is
worth stating in its general form, because neither of the two existing guards
sees it.  `bc.rebind()` guarantees that every binding of a name moves
together; `bc.rebind_by_operating_point()` guarantees that no member of the
class is silently left behind.  Both are statements about WHICH NAMES MOVE.
Neither says anything about a replacement that READS one of the moved names.
A replacement built on a name that the same flip rebinds is not stale and is
not unassigned -- it is correctly installed and reads a different object than
it was validated against.

Section 1 measures it.  Section 2 decomposes it into its two mechanisms and
finds that they PARTIALLY CANCEL, so repairing either one alone makes the
error larger.  Section 3 is the repair: bind the contact replacement to the
shipped object CAPTURED AT ENTRY rather than to the live name, which
reproduces the committed contact-only form B_T **bitwise**.  Sections 4-6 are
the flip itself, installed through the two managers the 2026-10-20 note named
-- `bc.rebind()` for the charge, `bc.rebind_by_operating_point()` for the
capacitance -- with the cost table regenerated with both halves in the circuit
for the first time.

WHAT THIS DOES NOT DO.  It does not edit `graphene_fet_model.py`, and the
chapter default is still form A.  `on_default_flip()` is an installable
default, validated and costed, and the one-line change that makes it the
chapter's default is now a decision rather than a build.  No committed number
in the thesis is rewritten here; Section 4.6.9 and Section 7.9 item 16 are
annotated in place with what composition costs.

Method notes, following this repository's standing rules.

  * Every patch substitutes ONE shipped name, and the contact patch now
    substitutes it for an object captured before the flip rather than for
    whatever the name holds at call time.  That is the whole repair.
  * ONE ARRIVAL COUNTER PER PATCHED NAME, each checked at a call that reads
    that name (the 2026-10-11 rule).
  * Validations against EXACTLY KNOWN values, not plausible ranges, after
    2026-09-17 and again 2026-10-19: the identity assignment must be bitwise
    neutral; the flip's channel half must reproduce `on_consistent_argument()`
    bitwise; the flip's contact half must reproduce
    `contact_on_consistent_argument()` bitwise; `n_contact(dWF = 0)` must be
    exactly 0.0; the electron/hole antisymmetry must be bitwise; and an
    incomplete assignment must raise and rebind NOTHING.
  * A POSITIVE CONTROL ON THE REPAIR.  A repair that is not load-bearing is
    indistinguishable from no repair, so V8 asserts that the unrepaired flip
    DOES differ -- if the two ever agree, the repair has become a no-op and
    this module's premise has silently expired.
"""

import contextlib
import io
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import graphene_fet_model as gfet
import graphene_contact_doping_model as cd
import graphene_quantum_capacitance_argument as qca
import graphene_contact_argument_blast_radius as blast
import graphene_binding_class as bc
import contact_resistance_crossover as crx

VG_ON = qca.VG_ON
VG_RF = qca.VG_RF
VDS_RF = qca.VDS_RF
DV_RF = qca.DV_RF

passed = 0
failed = 0


def check(label, ok, detail=''):
    global passed, failed
    if ok:
        passed += 1
        print('    PASS  %s' % label)
    else:
        failed += 1
        print('    FAIL  %s' % label)
    if detail:
        print('          %s' % detail)


def ulps(a, b):
    if a == b:
        return 0.0
    return abs(a - b) / max(np.spacing(abs(a)), np.spacing(abs(b)))


METALS = sorted(cd.METAL_WORK_FUNCTIONS.items(), key=lambda kv: kv[1])


# ---------------------------------------------------------------------------
# The contact replacement, bound to a CAPTURED object
# ---------------------------------------------------------------------------

_ARRIVED = {'channel_charge': 0, 'channel_cq': 0, 'contact_cq': 0}


def _reset_arrivals():
    for k in list(_ARRIVED):
        _ARRIVED[k] = 0


def make_contact_cq(shipped_cq, live=False):
    """Build the contact's form-B_T capacitance.

    `shipped_cq` is the object `gfet.quantum_capacitance` held BEFORE the flip.
    With live=True the replacement instead reads the live `gfet.
    quantum_capacitance` name, which is what `blast._cq_consistent_contact`
    does and is the defect Section 1 measures.  Everything else is identical
    between the two, so the only difference is where the inner C_q comes from.

    The solve itself is `qca._fixed_point` with `gfet.C_ox` substituted by
    `cd.C_interface`, exactly as `blast.n_consistent` does it -- the committed
    form-B_T solver, unmodified, so the contact solve and the channel solve
    cannot drift apart.

    KNOWN HAZARD, stated rather than hidden: that substitution mutates a
    module global for the duration of the solve, so this replacement is not
    reentrant with a concurrent channel call.  Nothing in this repository is
    concurrent, and V7 measures the interleaving that would expose it.
    """
    def cq_of_n(n, T=None):
        T = gfet.T if T is None else T
        nn = np.abs(np.asarray(n, dtype=float))
        arg = qca.A_FERMI * np.sqrt(nn)
        if live:
            return gfet.quantum_capacitance(arg, T=T)
        return shipped_cq(arg, T=T)

    def replacement(V, T=300.0):
        _ARRIVED['contact_cq'] += 1
        saved = gfet.C_ox
        try:
            gfet.C_ox = cd.C_interface
            n = qca._fixed_point(V, cq_of_n)
        finally:
            gfet.C_ox = saved
        return cq_of_n(n, T=T)

    return replacement


def _counting(key, fn):
    def wrapper(*a, **kw):
        _ARRIVED[key] += 1
        return fn(*a, **kw)
    return wrapper


# ---------------------------------------------------------------------------
# THE FLIP
# ---------------------------------------------------------------------------

CHANNEL_SERIES = ('channel: self-consistency against C_ox, the gated '
                  'geometry Section 4.6.9 derives form B_T for')
CONTACT_SERIES = ('contact: self-consistency against C_interface = %.4e '
                  'F/m^2, 65.93x C_ox, so the series factor is 0.555-0.806 '
                  'rather than 0.999' % cd.C_interface)


@contextlib.contextmanager
def on_default_flip(live_contact=False, report=None):
    """Install form B_T as the default EVERYWHERE the charge path reaches.

    Two rebindings, because the two names have different shapes, which is the
    2026-10-14 finding (ii) and the 2026-10-17 generalisation of it:

      * `carrier_density`'s identity class is {gfet, crx} and both members sit
        in the channel, so they take THE SAME replacement -- `bc.rebind()`.
      * `quantum_capacitance`'s identity class spans the channel and the
        contact, which are two operating points, so the class must be
        PARTITIONED -- `bc.rebind_by_operating_point()`, with a stated reason
        per member and an exception if any member is unassigned.

    This is the first caller to give `rebind_by_operating_point` a non-None
    contact assignment.  Until 2026-10-17 the contact had no replacement of
    its own and `bc.CONTACT_OPERATING_POINT` left it alone with a reason; that
    exclusion was correct exactly as long as that was true.

    live_contact=True installs the UNREPAIRED contact replacement (Section 1).
    """
    shipped_cq = gfet.quantum_capacitance
    contact_cq = make_contact_cq(shipped_cq, live=live_contact)

    assignment = {
        ('graphene_fet_model', 'quantum_capacitance'):
            (_counting('channel_cq', qca._quantum_capacitance_consistent),
             CHANNEL_SERIES),
        ('graphene_contact_doping_model', 'quantum_capacitance'):
            (contact_cq, CONTACT_SERIES),
    }
    with contextlib.ExitStack() as stack:
        stack.enter_context(bc.rebind(
            gfet, 'carrier_density',
            _counting('channel_charge', qca._carrier_density_consistent)))
        installed, left = stack.enter_context(bc.rebind_by_operating_point(
            gfet, 'quantum_capacitance', assignment, report=report))
        _reset_arrivals()
        yield installed, left


def census_contact():
    return {m: float(cd.contact_edge_carrier_density(wf)) for m, wf in METALS}


def census_recal():
    return cd.recalibrate_metal_rc()


# ---------------------------------------------------------------------------

def main():
    print('=' * 78)
    print("ITEM 16's FLIP: the two halves do not compose")
    print('=' * 78)
    print()
    print('C_ox        = %.6e F/m^2' % gfet.C_ox)
    print('C_interface = %.6e F/m^2   (ratio %.4f)'
          % (cd.C_interface, cd.C_interface / gfet.C_ox))
    print()

    base_contact = census_contact()
    base_chan = qca.census()
    base_recal = census_recal()

    # -- Section 0 ----------------------------------------------------------
    print('-' * 78)
    print('0.  THE IDENTITY ASSIGNMENT IS BITWISE NEUTRAL (the instrument is')
    print('    not the perturbation)')
    print('-' * 78)
    print()
    cls_q = sorted(bc.identity_class(gfet, 'quantum_capacitance'))
    cls_n = sorted(bc.identity_class(gfet, 'carrier_density'))
    print('    identity class of quantum_capacitance : %s'
          % ', '.join('%s.%s' % s for s in cls_q))
    print('    identity class of carrier_density     : %s'
          % ', '.join('%s.%s' % s for s in cls_n))
    print()
    print('    Giving every member of the capacitance class ITS OWN shipped')
    print('    object exercises the whole partition machinery and must move')
    print('    nothing.  Any departure here and Sections 1-6 are measuring')
    print('    the harness.')
    print()
    ident = {s: (getattr(sys.modules[s[0]], s[1]), 'identity: shipped object')
             for s in cls_q}
    with bc.rebind_by_operating_point(gfet, 'quantum_capacitance', ident):
        id_contact = census_contact()
        id_chan = qca.census()
        id_recal = census_recal()
    worst_c = max(ulps(id_contact[m], base_contact[m]) for m, _ in METALS)
    check('identity assignment: all 7 n_contact bitwise (0 ulp)',
          all(id_contact[m] == base_contact[m] for m, _ in METALS),
          'worst departure %.1f ulp' % worst_c)
    check('identity assignment: all %d channel census entries bitwise'
          % len(base_chan),
          all(id_chan[k] == base_chan[k] for k in base_chan),
          'keys: %s' % ', '.join(sorted(base_chan)))
    check('identity assignment: all %d R_extra bitwise' % len(base_recal),
          all(id_recal[m][2] == base_recal[m][2] for m in base_recal))

    # -- Section 1 ----------------------------------------------------------
    print()
    print('-' * 78)
    print('1.  THE DEFECT: the unrepaired flip moves n_contact by -51 % to')
    print('    -62 % away from the committed contact-only form B_T')
    print('-' * 78)
    print()
    with blast.contact_on_consistent_argument():
        bt_contact = census_contact()
    with qca.on_consistent_argument():
        chan_only_contact = census_contact()
    with on_default_flip(live_contact=True):
        naive_contact = census_contact()
    with on_default_flip(live_contact=False):
        flip_contact = census_contact()
        flip_arrivals = dict(_ARRIVED)

    print('    %-5s %14s %14s %14s %14s' % ('metal', 'form A', 'B_T contact',
                                            'flip unrepaired', 'flip repaired'))
    for m, _ in METALS:
        print('    %-5s %14.6e %14.6e %14.6e %14.6e'
              % (m, base_contact[m], bt_contact[m], naive_contact[m],
                 flip_contact[m]))
    print()
    print('    %-5s %12s %12s %12s %10s'
          % ('metal', 'B_T vs A %', 'unrep vs B_T %', 'rep vs B_T',
             '|unrep/A|'))
    devs, ratios = [], []
    for m, _ in METALS:
        if bt_contact[m] == 0.0:
            print('    %-5s %12s %12s %12s %10s'
                  % (m, 'n = 0', 'n = 0', 'exact 0', '-'))
            continue
        d_bt = 100 * (bt_contact[m] / base_contact[m] - 1.0)
        d_nv = 100 * (naive_contact[m] / bt_contact[m] - 1.0)
        devs.append(d_nv)
        # PER-METAL ratio.  The first draft of the summary beneath this table
        # divided the smallest-magnitude naive deviation by the LARGEST-
        # magnitude committed effect and the other way round, quoting
        # "4.1x to 13.0x" from two pairs of numbers belonging to different
        # metals.  That is this repository's recorded dominant failure mode
        # -- a claim answering a weaker question than the one it was written
        # to ask -- and it was caught by READING the transcript, not by any
        # check in this file.
        ratios.append((abs(d_nv) / abs(d_bt), m))
        print('    %-5s %+12.4f %+12.4f %12s %10.2f'
              % (m, d_bt, d_nv,
                 'bitwise' if flip_contact[m] == bt_contact[m]
                 else '%.1f ulp' % ulps(flip_contact[m], bt_contact[m]),
                 abs(d_nv) / abs(d_bt)))
    lo, lom = min(ratios)
    hi, him = max(ratios)
    print()
    print('    The effect item 16 installs at the contact is -4.7696 % to')
    print('    -12.4964 %%.  The unrepaired composition is %+.4f %% to %+.4f %%,'
          % (max(devs), min(devs)))
    print('    and METAL BY METAL it is %.2fx (%s) to %.2fx (%s) the effect it'
          % (lo, lom, hi, him))
    print('    is installing, in the SAME direction -- so it would not have')
    print('    looked like a sign error.')
    print()
    check('V1 the channel half ALONE moves no contact number (the 2026-10-17 '
          'hole, unchanged)',
          all(chan_only_contact[m] == base_contact[m] for m, _ in METALS),
          'zero arrivals at cd.quantum_capacitance, as measured 2026-10-17')
    check('V2 the repaired flip DOES reach cd.quantum_capacitance '
          '(the 2026-10-17 hole is closed)',
          flip_arrivals['contact_cq'] == len(METALS),
          'contact_cq arrivals %d for %d metals'
          % (flip_arrivals['contact_cq'], len(METALS)))
    # ONE COUNTER PER PATCHED NAME, EACH CHECKED AT A CALL THAT READS THAT
    # NAME (the 2026-10-11 rule).  The contact census above exercises only
    # the contact binding, so its channel counters are correctly zero, and
    # asserting on them THERE would report a fault that does not exist --
    # which is how that rule was arrived at in the first place.  The channel
    # counters are checked here, at a channel call.
    with on_default_flip():
        float(gfet.carrier_density(VG_ON))
        float(gfet.quantum_capacitance(DV_RF))
        chan_arrivals = dict(_ARRIVED)
    check('V2b and the two channel bindings are reached at a CHANNEL call',
          chan_arrivals['channel_charge'] > 0
          and chan_arrivals['channel_cq'] > 0,
          'channel_charge %d, channel_cq %d, contact_cq %d (correctly 0 here)'
          % (chan_arrivals['channel_charge'], chan_arrivals['channel_cq'],
             chan_arrivals['contact_cq']))

    # -- Section 2 ----------------------------------------------------------
    print()
    print('-' * 78)
    print('2.  TWO MECHANISMS, AND THEY PARTIALLY CANCEL')
    print('-' * 78)
    print()
    print('    (i)  THE INNER C_q.  The contact solve reads the live')
    print('         `gfet.quantum_capacitance`, which the channel half has')
    print('         rebound to a form-B_T solve.')
    print('    (ii) THE INNER SERIES CAPACITANCE.  That inner solve reaches')
    print('         `qca.series_factor`, which reads `gfet.C_ox` at call time')
    print('         -- and the contact solve has it substituted by')
    print('         C_interface for the duration.  So the inner channel solve')
    print('         runs against the contact\'s series capacitance too.')
    print()
    print('    (i) alone is isolated by restoring C_ox for the inner call.')
    print()

    def contact_cq_cause_i_only(shipped_cox):
        def cq_of_n(n, T=None):
            T = gfet.T if T is None else T
            nn = np.abs(np.asarray(n, dtype=float))
            arg = qca.A_FERMI * np.sqrt(nn)
            s = gfet.C_ox
            gfet.C_ox = shipped_cox
            try:
                return gfet.quantum_capacitance(arg, T=T)
            finally:
                gfet.C_ox = s

        def replacement(V, T=300.0):
            saved = gfet.C_ox
            try:
                gfet.C_ox = cd.C_interface
                n = qca._fixed_point(V, cq_of_n)
            finally:
                gfet.C_ox = saved
            return cq_of_n(n, T=T)
        return replacement

    shipped_cox = gfet.C_ox
    with qca.on_consistent_argument():
        saved = cd.quantum_capacitance
        cd.quantum_capacitance = contact_cq_cause_i_only(shipped_cox)
        try:
            cause_i = census_contact()
        finally:
            cd.quantum_capacitance = saved

    print('    %-5s %12s %12s %12s' % ('metal', '(i) only %', '(i)+(ii) %',
                                       '(ii) moves'))
    offsets = []
    for m, _ in METALS:
        if bt_contact[m] == 0.0:
            continue
        d_i = 100 * (cause_i[m] / bt_contact[m] - 1.0)
        d_b = 100 * (naive_contact[m] / bt_contact[m] - 1.0)
        offsets.append(d_b - d_i)
        print('    %-5s %+12.4f %+12.4f %+12.4f' % (m, d_i, d_b, d_b - d_i))
    print()
    check('V3 mechanism (ii) OFFSETS mechanism (i) for every metal',
          all(o > 0 for o in offsets),
          'offsets %+.4f to %+.4f percentage points -- repairing (ii) alone '
          'would make the error LARGER, which is why the composition failure '
          'has to be repaired at (i)' % (min(offsets), max(offsets)))

    # -- Section 3 ----------------------------------------------------------
    print()
    print('-' * 78)
    print('3.  THE REPAIR, VALIDATED BITWISE AGAINST BOTH COMMITTED HALVES')
    print('-' * 78)
    print()
    with on_default_flip():
        flip_chan = qca.census()
        flip_contact2 = census_contact()
        flip_recal = census_recal()
    with qca.on_consistent_argument():
        chan_only_chan = qca.census()

    check('V4 flip CONTACT half == committed contact-only form B_T, bitwise '
          'on all 7 metals',
          all(flip_contact2[m] == bt_contact[m] for m, _ in METALS),
          'worst %.1f ulp'
          % max(ulps(flip_contact2[m], bt_contact[m]) for m, _ in METALS))
    chan_keys = sorted(base_chan)
    chan_bit = all(flip_chan[k] == chan_only_chan[k] for k in chan_keys)
    check('V5 flip CHANNEL half == committed on_consistent_argument(), '
          'bitwise on all %d census entries' % len(chan_keys), chan_bit,
          'worst %.1f ulp at %s'
          % max(((ulps(flip_chan[k], chan_only_chan[k]), k)
                 for k in chan_keys)) if not chan_bit else 'all 0 ulp')
    nz = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE))
    with on_default_flip():
        nz_flip = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE))
    check('V6 zero null: n_contact(dWF = 0) is EXACTLY 0.0 under the flip',
          nz_flip == 0.0 and nz == 0.0,
          'shipped %r, flip %r' % (nz, nz_flip))
    dwf = 0.45
    with on_default_flip():
        nP = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE + dwf))
        nM = float(cd.contact_edge_carrier_density(cd.W_GRAPHENE - dwf))
    check('V7 electron/hole antisymmetry at dWF = %.2f eV is BITWISE under '
          'the flip' % dwf, nP == -nM,
          'n(+) %r   n(-) %r   sum %r' % (nP, nM, nP + nM))
    check('V8 POSITIVE CONTROL on the repair: the unrepaired flip differs on '
          'all 6 nonzero metals',
          all(naive_contact[m] != bt_contact[m]
              for m, _ in METALS if bt_contact[m] != 0.0),
          'if this ever passes trivially, the repair has become a no-op and '
          'this module\'s premise has expired')

    print()
    print('    INTERLEAVING (the stated hazard in make_contact_cq).  The')
    print('    contact solve mutates gfet.C_ox for its duration, so a channel')
    print('    call made DURING it would see the wrong series capacitance.')
    print('    Measured by interleaving the two under one flip:')
    print()
    with on_default_flip():
        inter = []
        for m, wf in METALS:
            inter.append(float(cd.contact_edge_carrier_density(wf)))
            inter.append(float(gfet.carrier_density(VG_ON)))
    il_contact = inter[0::2]
    il_chan = inter[1::2]
    check('V9 interleaved: every n_contact equals its value computed alone, '
          'bitwise',
          all(a == bt_contact[m] for a, (m, _) in zip(il_contact, METALS)))
    check('V10 interleaved: every n(3.5 V) equals its flip value, bitwise',
          all(a == flip_chan['n_3V5'] for a in il_chan),
          'n(3.5 V) under the flip = %.10e' % flip_chan['n_3V5'])

    # -- Section 4 ----------------------------------------------------------
    print()
    print('-' * 78)
    print('4.  THE PARTITION REFUSES AN INCOMPLETE MAP (negative control)')
    print('-' * 78)
    print()
    partial = {('graphene_fet_model', 'quantum_capacitance'):
               (qca._quantum_capacitance_consistent, CHANNEL_SERIES)}
    raised = None
    try:
        with bc.rebind_by_operating_point(gfet, 'quantum_capacitance',
                                          partial):
            pass
    except bc.UnassignedClassMember as exc:
        raised = exc
    check('V11 omitting the contact member raises UnassignedClassMember',
          raised is not None,
          str(raised).split('.  ')[0] if raised else 'NOTHING RAISED')
    after_contact = census_contact()
    after_chan = qca.census()
    check('V12 and NOTHING was rebound: every number bitwise equal to shipped',
          all(after_contact[m] == base_contact[m] for m, _ in METALS)
          and all(after_chan[k] == base_chan[k] for k in base_chan))
    noreason = {s: (getattr(sys.modules[s[0]], s[1]), '') for s in cls_q}
    raised2 = None
    try:
        with bc.rebind_by_operating_point(gfet, 'quantum_capacitance',
                                          noreason):
            pass
    except bc.UnassignedClassMember as exc:
        raised2 = exc
    check('V13 an assignment with no stated reason is refused too',
          raised2 is not None,
          str(raised2) if raised2 else 'NOTHING RAISED')

    # -- Section 5 ----------------------------------------------------------
    print()
    print('-' * 78)
    print('5.  THE FLIP\'S COST, WITH BOTH HALVES IN THE CIRCUIT AT ONCE FOR')
    print('    THE FIRST TIME')
    print('-' * 78)
    print()
    print('    Channel census (identical to on_consistent_argument(), V5):')
    print()
    print('    %-14s %18s %18s %10s' % ('quantity', 'form A', 'flip', 'change %'))
    for k in chan_keys:
        a, b = base_chan[k], flip_chan[k]
        ch = 100 * (b / a - 1.0) if a != 0 else float('nan')
        print('    %-14s %18.10e %18.10e %+10.4f' % (k, a, b, ch))
    print()
    print('    Contact census (NEW: never before measured with the channel')
    print('    half installed simultaneously):')
    print()
    print('    %-5s %18s %18s %10s' % ('metal', 'form A', 'flip', 'change %'))
    for m, _ in METALS:
        a, b = base_contact[m], flip_contact2[m]
        ch = 100 * (b / a - 1.0) if a != 0 else 0.0
        print('    %-5s %18.10e %18.10e %+10.4f' % (m, a, b, ch))
    print()
    print('    Section 4.7 recalibration.  R_extra is the junction')
    print('    contribution this model computes; R_tr is the implied residual')
    print('    Rc_measured - R_extra.')
    print()
    print('    %-5s %12s %12s %12s %10s %12s %12s'
          % ('metal', 'Rc_meas', 'R_extra A', 'R_extra flip', 'change %',
             'R_tr A', 'R_tr flip'))
    sign_flips = []
    for m in sorted(base_recal, key=lambda x: base_recal[x][0]):
        wf, rc, ra, rt = base_recal[m]
        _, _, rb, rtb = flip_recal[m]
        ch = 100 * (rb / ra - 1.0)
        print('    %-5s %12.4f %12.4f %12.4f %+10.4f %12.4f %12.4f'
              % (m, rc, ra, rb, ch, rt, rtb))
        if (rt > 0) != (rtb > 0):
            sign_flips.append((m, rt, rtb))
    print()
    if sign_flips:
        for m, rt, rtb in sign_flips:
            print('    %s\'s implied transmission residual changes sign: '
                  '%+.4f -> %+.4f Ohm.um.' % (m, rt, rtb))
        print('    Section 4.10 found Pd\'s residual -- the only positive one,')
        print('    and the only evidence the additive decomposition survives')
        print('    -- flips on a 9.6 % error in R_extra.  The flip consumes')
        print('    that margin, so under form B_T everywhere the additive')
        print('    decomposition fails for every calibrated metal.')
    else:
        print('    No implied-residual sign changes.')
    print()
    neg_A = sum(1 for m in base_recal if base_recal[m][3] < 0)
    neg_B = sum(1 for m in flip_recal if flip_recal[m][3] < 0)
    check('V14 additive decomposition: %d of %d residuals negative under form '
          'A, %d of %d under the flip'
          % (neg_A, len(base_recal), neg_B, len(flip_recal)),
          neg_B == len(flip_recal),
          'the 2026-10-17 measurement, reproduced with both halves installed')

    # -- Section 6 ----------------------------------------------------------
    print()
    print('-' * 78)
    print('6.  WHAT IS LEFT OF ITEM 16')
    print('-' * 78)
    print()
    print('    The flip is built, validated and costed, and it is installable')
    print('    with one `with on_default_flip():`.  What is NOT done is')
    print('    editing graphene_fet_model.py so that form B_T is what the')
    print('    chapter computes by default, because that is a decision about')
    print('    what the thesis claims and not a build.  The 2026-10-20 note')
    print('    that "the flip and nothing else is left" was wrong in one')
    print('    specific way, now repaired: composition was left as well.')
    print()
    print('    AND THE GENERAL FORM, which no guard in this repository')
    print('    currently covers.  bc.rebind() asks "do all bindings of this')
    print('    name move together".  bc.rebind_by_operating_point() asks "is')
    print('    any member of the class unassigned".  Neither asks "does any')
    print('    REPLACEMENT read a name this flip rebinds".  The audit that')
    print('    would have caught today\'s defect without running anything is')
    print('    a one-line question per replacement, and it is now the first')
    print('    item on the standing list.')
    print()

    make_figure(base_contact, bt_contact, naive_contact, base_recal,
                flip_recal, ratios)

    print('=' * 78)
    print('%d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return failed


def make_figure(base_contact, bt_contact, naive_contact, base_recal,
                flip_recal, ratios):
    """Two panels.  Rendered and INSPECTED, not merely produced -- the
    2026-10-20 session found two unreadable committed figures on its first
    look and made inspection a standing item.  Looking at this one found
    three more, and the third is the reason `ratios` is a PARAMETER:

      1. Panel (b)'s legend sat on top of Cu's two tallest bars and its
         R_tr annotation, which was grey-on-white behind the legend box and
         unreadable.  Fixed with headroom rather than by moving the legend
         into another collision.
      2. The same annotation collided with the panel's top edge.
      3. PANEL (a)'S TITLE WAS FALSE.  It read "by 4-5x", the figure of the
         superseded unpaired-ratio calculation this session corrected in
         Section 1; the per-metal answer is 4.86x to 10.77x.  The text was
         repaired and the TITLE was not, because a title is not in the
         transcript and nothing reads it.  It is now computed from the same
         `ratios` the table prints, so the two cannot disagree again.
         Neither this nor (1) was reachable from a transcript reading
         "18 passed, 0 failed"."""
    names = [m for m, _ in METALS]
    x = np.arange(len(names))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.2))

    w = 0.27
    ax1.bar(x - w, [base_contact[m] / 1e16 for m in names], w,
            label='form A (shipped)', color='#4878a8')
    ax1.bar(x, [bt_contact[m] / 1e16 for m in names], w,
            label=r'form $B_T$ (committed, and the repaired flip)',
            color='#5fa052')
    ax1.bar(x + w, [naive_contact[m] / 1e16 for m in names], w,
            label='unrepaired composition', color='#c44e52')
    ax1.axhline(0.0, color='k', lw=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(names)
    ax1.set_xlabel('contact metal (increasing work function)')
    ax1.set_ylabel(r'$n_{\rm contact}$  ($10^{16}\ \rm m^{-2}$)')
    lo = min(r for r, _ in ratios)
    hi = max(r for r, _ in ratios)
    ax1.set_title('(a)  Composing the two halves overshoots the\n'
                  'effect it installs by %.1fx to %.1fx' % (lo, hi),
                  fontsize=10)
    ax1.legend(fontsize=7.5, loc='lower left')
    ax1.grid(axis='y', alpha=0.3)

    rnames = sorted(base_recal, key=lambda m: base_recal[m][0])
    xr = np.arange(len(rnames))
    wr = 0.26
    ax2.bar(xr - wr, [base_recal[m][1] for m in rnames], wr,
            label=r'$R_c$ measured (literature)', color='#8c8c8c')
    ax2.bar(xr, [base_recal[m][2] for m in rnames], wr,
            label=r'$R_{\rm extra}$, form A', color='#4878a8')
    ax2.bar(xr + wr, [flip_recal[m][2] for m in rnames], wr,
            label=r'$R_{\rm extra}$, flip (form $B_T$ everywhere)',
            color='#5fa052')
    for i, m in enumerate(rnames):
        rt = flip_recal[m][3]
        ax2.annotate(r'$R_{\rm tr}=%+.0f$' % rt, (xr[i] + wr, flip_recal[m][2]),
                     textcoords='offset points', xytext=(0, 5),
                     ha='center', fontsize=7.5)
    # HEADROOM.  Without it the legend covers Cu's bars and its R_tr
    # annotation; see this function's docstring, defect (1).
    top = max(max(base_recal[m][1], base_recal[m][2], flip_recal[m][2])
              for m in rnames)
    ax2.set_ylim(0.0, 1.45 * top)
    ax2.set_xticks(xr)
    ax2.set_xticklabels(rnames)
    ax2.set_xlabel('calibrated metal')
    ax2.set_ylabel(r'resistance ($\Omega\cdot\mu$m)')
    ax2.set_title('(b)  Section 4.7: every implied transmission residual\n'
                  r'is negative under the flip', fontsize=10)
    ax2.legend(fontsize=7.5, loc='upper left')
    ax2.grid(axis='y', alpha=0.3)

    fig.tight_layout()
    fig.savefig('default_flip.png', dpi=150)
    plt.close(fig)
    print('    figure written: default_flip.png')


if __name__ == '__main__':
    buf = io.StringIO()

    class Tee:
        def write(self, s):
            buf.write(s)
            sys.__stdout__.write(s)

        def flush(self):
            sys.__stdout__.flush()

    sys.stdout = Tee()
    rc = main()
    sys.stdout = sys.__stdout__
    with open('default_flip_output.txt', 'w') as fh:
        fh.write(buf.getvalue())
    sys.exit(1 if rc else 0)
