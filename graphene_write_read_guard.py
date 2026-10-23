"""graphene_write_read_guard.py -- the write/read relation as a CHECK THAT FAILS.

WHY THIS EXISTS, and what it is not.

2026-10-22 §4 found that `gfet.C_ox` is READ by all five channel/contact
replacements and WRITTEN by three of them, and that it is in no identity
class -- so it is invisible to `bc.rebind()`, to
`bc.rebind_by_operating_point()`, and to that session's own question, all
three of which are about names a flip REBINDS.  The conclusion recorded there
was that what governs composition is not a set of names but a BIPARTITE
RELATION between what replacements write and what replacements read.

That session REPORTED the relation.  This module REFUSES a composition whose
relation contains an unadjudicated edge.  The standing item it discharges --
"A CENSUS THAT ENUMERATES A FAULT FOREVER IS NOT A GUARD", carried since
2026-10-17 and flagged on 2026-10-22 as having its strongest instance in §4
itself -- is the whole reason it is a `raise` and not a table.

THE VOCABULARY IS NOT NEW, DELIBERATELY.  `bc.rebind()` refuses an exclusion
without a stated reason.  `bc.rebind_by_operating_point()` refuses an
UNASSIGNED class member, raises before anything is rebound, and makes the
reason mandatory for "left alone" as well as for "replaced".  This module is
the same refusal applied to an edge instead of to a name: every edge must
carry a verdict and a reason, and an edge with neither is
`UnadjudicatedWriteReadEdge` with nothing run.

WHAT IS NEW IS THAT EVERY VERDICT IS MEASURED, NOT TRUSTED.

`rebind_by_operating_point` can check completeness statically -- it knows the
identity class.  An edge's verdict is a claim about REACHABILITY, which no
static set contains, so a guard that accepted the claim would be a census
with a stricter tone.  Each verdict is therefore checked against an order
probe, and it is checked IN BOTH DIRECTIONS:

    GOVERNED   the name is in the flip's own rebind set, so an existing
               mechanism moves it.  Refused if the name is NOT in that set.
    IMPOSSIBLE the edge cannot occur under this flip.  Refused if the probe
               OBSERVES it.
    OCCURS     it does occur, and the reason says why that is sound.
               Refused if the probe NEVER observes it.

The last one is the half a census cannot have, and it is the half this
repository's recorded dominant failure mode calls for: a check that only
catches a wrongly-permissive claim answers a weaker question than it was
written to ask.  An `OCCURS` reason describing an edge that has stopped
occurring is a sentence about nothing, kept alive by the fact that nobody
asked it the opposite question.

THE RELATION IS WIDER THAN EITHER SESSION'S QUESTION.  A name can be written
by a replacement (2026-10-22's `gfet.C_ox`) or by the FLIP ITSELF, which is
every member of the rebind set (2026-10-21's "a replacement that reads a
rebound name is neither stale nor unassigned").  Those were recorded as two
findings.  They are two sources of writes into ONE relation, and composing
them is what makes the composed column -- the one 2026-10-22 said answers the
question -- computable rather than argued.

ORDER, AND HOW AN EDGE IS ATTRIBUTED.  An edge (writer, name, reader) is a
fault only if the write can precede the read.  The probe keeps a time-ordered
trace of two event kinds -- a replacement being ENTERED, and a cross-module
attribute being WRITTEN -- and maintains the stack of active replacements, so
a write is attributed to the innermost replacement on the stack or to the
flip when none is active.  No attribution table is maintained by hand; it
falls out of the trace.

NEUTRALITY.  The probe installs identity wrappers and an ordered variant of
reader C.  Both are neutral by construction and both are CHECKED anyway,
bitwise, on the seven committed contact densities before any count is used --
the 2026-10-13 rule that an instrument which perturbs what it measures
reports its own arrival, and the 2026-10-22 note that "that argument has been
wrong before".

Usage:  PYTHONPATH=... python3 graphene_write_read_guard.py
"""

import contextlib
import types

import numpy as np

import graphene_fet_model as gfet
import graphene_contact_doping_model as cd
import graphene_quantum_capacitance_argument as qca
import graphene_contact_argument_blast_radius as blast
import graphene_default_flip as flip
import graphene_binding_class as bc
import graphene_replacement_read_audit as rra


GOVERNED = 'GOVERNED'
IMPOSSIBLE = 'IMPOSSIBLE'
OCCURS = 'OCCURS'
VERDICTS = (GOVERNED, IMPOSSIBLE, OCCURS)


class UnadjudicatedWriteReadEdge(RuntimeError):
    """A (writer, name, reader) edge was given neither a verdict nor a stated
    reason.  Raised BEFORE anything is run, like
    `bc.UnassignedClassMember`."""


class FalseEdgeVerdict(RuntimeError):
    """An edge's verdict disagrees with what the order probe measured."""


# ---------------------------------------------------------------------------
# An ORDERED reader C.  rra.tracking() keeps sets; reachability needs order.
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def ordered_tracking(trace):
    """Append ('write', site) and ('read', site) to `trace` in call order.

    READS ARE TRACED, and the first draft of this module did not trace them.
    That draft closed an edge on every ARRIVAL of a reader that followed a
    write, which asks "did this replacement run after that name moved" --
    weaker than "did this replacement READ that name after it moved", and
    wrong in the permissive direction.  It reported twelve observed edges
    under the unrepaired flip, four of them on names the reader does not read
    at all, and the observed set was not even a subset of the relation.  A
    probe whose own over-report would have to be corrected by intersecting it
    with a static set is a probe that has not measured reachability; this one
    derives the edge from the read.
    """
    MT = types.ModuleType
    mods = list(rra.repo_modules().values())

    class Tracked(MT):
        def __getattribute__(self, name):
            if not name.startswith('__'):
                trace.append(('read',
                              (MT.__getattribute__(self, '__name__'), name)))
            return MT.__getattribute__(self, name)

        def __setattr__(self, name, value):
            if not name.startswith('__'):
                trace.append(('write',
                              (MT.__getattribute__(self, '__name__'), name)))
            MT.__setattr__(self, name, value)

    saved = [(m, m.__class__) for m in mods]
    try:
        for m in mods:
            m.__class__ = Tracked
        yield trace
    finally:
        for m, cls in saved:
            m.__class__ = cls


# ---------------------------------------------------------------------------
# Arrival wrappers: identity, and they maintain the active-replacement stack.
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def arrival_wrappers(trace, stack):
    """Identity-wrap each replacement at its DEFINITION SITE.

    The two `qca` replacements and `blast._cq_consistent_contact` are module
    attributes and are wrapped directly.  The products of
    `flip.make_contact_cq` are built fresh on each call, so there is no site
    to wrap: the FACTORY is wrapped instead and wraps its own return value,
    which is why a product's arrivals are attributed to the factory's label
    rather than to an object identity that does not outlive the flip.
    """
    saved = []

    def wrap(label, fn):
        def wrapper(*a, **kw):
            trace.append(('arrive', label))
            stack.append(label)
            try:
                return fn(*a, **kw)
            finally:
                stack.pop()
                trace.append(('return', label))
        return wrapper

    sites = (
        (qca, '_carrier_density_consistent',
         'qca._carrier_density_consistent'),
        (qca, '_quantum_capacitance_consistent',
         'qca._quantum_capacitance_consistent'),
        (blast, '_cq_consistent_contact', 'blast._cq_consistent_contact'),
    )
    try:
        for mod, attr, label in sites:
            obj = getattr(mod, attr)
            saved.append((mod, attr, obj))
            setattr(mod, attr, wrap(label, obj))

        factory = flip.make_contact_cq
        saved.append((flip, 'make_contact_cq', factory))

        def wrapped_factory(*a, **kw):
            live = kw.get('live', a[1] if len(a) > 1 else None)
            product = factory(*a, **kw)
            return wrap('flip.make_contact_cq(live=%s)' % live, product)

        flip.make_contact_cq = wrapped_factory
        yield
    finally:
        for mod, attr, obj in reversed(saved):
            setattr(mod, attr, obj)


# ---------------------------------------------------------------------------
# The relation, derived
# ---------------------------------------------------------------------------

def replacement_write_read(replacements):
    """{label: (writes, reads)} measured by entering each replacement.

    Reader C is asked; nothing is re-derived and no list is maintained, so a
    replacement that grows a write is reported with it.
    """
    out = {}
    for label, fn, args in replacements:
        with rra.tracking() as log:
            fn(*args)
        out[label] = (set(log['write']), set(log['read']))
    return out


def candidate_edges(wr, flip_rebinds):
    """Every (writer, name, reader) the relation admits, before reachability.

    Two sources of writes, which were two separate findings and are one
    relation: a replacement's own `setattr` (2026-10-22's `gfet.C_ox`) and
    the flip's rebind set (2026-10-21's rebound names).
    """
    edges = set()
    for reader, (_, reads) in wr.items():
        for name in reads:
            for writer, (writes, _) in wr.items():
                if name in writes:
                    edges.add((writer, name, reader))
            if name in flip_rebinds:
                edges.add(('THE FLIP', name, reader))
    return edges


def observe(live_contact):
    """Run the contact census under the probe; return (observed, trace, census).

    An edge (writer, name, reader) is OBSERVED when a READ of `name`
    attributed to `reader` occurs after a WRITE of `name` attributed to
    `writer`.  Both attributions come from the stack of active replacements:
    the innermost one, or `THE FLIP` when none is active.  Nothing is
    attributed by a hand-maintained table.

    ATTRIBUTION IS BY DYNAMIC EXTENT, and that is a real difference from the
    direct-call sets.  A write made by a callee that is not itself a wrapped
    replacement is attributed to the replacement it ran under, which is the
    honest answer to "under whose composition did this name move".  It is why
    the dynamic relation can contain an edge the static one does not, and
    those edges are reported separately rather than merged away.
    """
    trace, stack = [], []
    with arrival_wrappers(trace, stack):
        with ordered_tracking(trace):
            with flip.on_default_flip(live_contact=live_contact):
                census = {m: float(cd.contact_edge_carrier_density(wf))
                          for m, wf in flip.METALS}
                cd.recalibrate_metal_rc()

    return edges_from_trace(trace), trace, census


def edges_from_trace(trace):
    """Replay a trace into observed edges.  Separated so the DERIVER itself
    can be put under exactly-known inputs (V13/V14) instead of only ever
    being exercised by a run whose right answer nobody knows independently.
    A deriver that returns few edges and a repository that has few are
    indistinguishable until this function is asked a question with a known
    answer."""
    written, observed, depth = {}, set(), []
    for kind, payload in trace:
        if kind == 'arrive':
            depth.append(payload)
            continue
        if kind == 'return':
            if depth:
                depth.pop()
            continue
        who = depth[-1] if depth else 'THE FLIP'
        if kind == 'write':
            written.setdefault(payload, set()).add(who)
        else:
            for w in written.get(payload, ()):
                observed.add((w, payload, who))
    return observed


# ---------------------------------------------------------------------------
# THE GUARD
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def require_write_read_adjudication(live_contact, replacements, adjudication,
                                    report=None):
    """Refuse a composition whose write/read relation is not fully adjudicated.

    adjudication: {(writer, name, reader): (verdict, reason)} with verdict in
        VERDICTS and a mandatory reason, in the shape
        `bc.rebind_by_operating_point` uses for assignments.

    Raises `UnadjudicatedWriteReadEdge` BEFORE running the probe if any
    candidate edge is unadjudicated, and `FalseEdgeVerdict` after the probe if
    any verdict disagrees with it -- in either direction.
    """
    flip_rebinds = set()
    for _, fac in rra.FLIPS:
        flip_rebinds |= rra.rebind_set(fac)
    wr = replacement_write_read(replacements)
    edges = candidate_edges(wr, flip_rebinds)

    missing = sorted(edges - set(adjudication))
    if missing:
        raise UnadjudicatedWriteReadEdge(
            'the write/read relation of this composition has %d unadjudicated '
            'edge(s): %s.  Give each a (verdict, reason) with verdict in %s; '
            'nothing has been run.'
            % (len(missing),
               ', '.join('%s --writes %s--> %s'
                         % (m[0], _dotted(m[1]), m[2]) for m in missing),
               '/'.join(VERDICTS)))
    for edge, spec in sorted(adjudication.items()):
        if not (isinstance(spec, tuple) and len(spec) == 2):
            raise UnadjudicatedWriteReadEdge(
                'adjudication for %s must be a (verdict, reason) pair' % (edge,))
        verdict, reason = spec
        if verdict not in VERDICTS:
            raise UnadjudicatedWriteReadEdge(
                'adjudication for %s has verdict %r, not one of %s'
                % (edge, verdict, '/'.join(VERDICTS)))
        if not reason:
            raise UnadjudicatedWriteReadEdge(
                'adjudication for %s has no stated reason' % (edge,))

    observed, trace, census = observe(live_contact)
    wrong = []
    for edge in sorted(edges):
        verdict, reason = adjudication[edge]
        writer, name, reader = edge
        if verdict == GOVERNED and name not in flip_rebinds:
            wrong.append((edge, verdict,
                          'claimed GOVERNED but %s is in no rebind set'
                          % _dotted(name)))
        elif verdict == IMPOSSIBLE and edge in observed:
            wrong.append((edge, verdict,
                          'claimed IMPOSSIBLE but the probe observed it'))
        elif verdict == OCCURS and edge not in observed:
            wrong.append((edge, verdict,
                          'claimed OCCURS but the probe never observed it -- '
                          'the stated reason describes nothing'))
    if wrong:
        raise FalseEdgeVerdict(
            '%d edge verdict(s) disagree with the order probe:\n%s'
            % (len(wrong), '\n'.join(
                '    %s --writes %s--> %s : %s'
                % (e[0], _dotted(e[1]), e[2], why)
                for e, _, why in wrong)))
    if report is not None:
        report(edges, observed, trace, census)
    yield edges, observed


# ---------------------------------------------------------------------------
# The composition under audit, and its adjudication
# ---------------------------------------------------------------------------

def replacements_under_audit():
    shipped_cq = gfet.quantum_capacitance
    return (
        ('qca._carrier_density_consistent',
         qca._carrier_density_consistent, (0.5,)),
        ('qca._quantum_capacitance_consistent',
         qca._quantum_capacitance_consistent, (0.5,)),
        ('blast._cq_consistent_contact',
         blast._cq_consistent_contact, (0.5,)),
        ('flip.make_contact_cq(live=False)',
         flip.make_contact_cq(shipped_cq, live=False), (0.5,)),
        ('flip.make_contact_cq(live=True)',
         flip.make_contact_cq(shipped_cq, live=True), (0.5,)),
    )


COX = ('graphene_fet_model', 'C_ox')
CQ = ('graphene_fet_model', 'quantum_capacitance')

# The adjudication of the REPAIRED composition (live_contact=False).
#
# Written as a map from edge to (verdict, reason), in the shape
# `bc.rebind_by_operating_point` uses for assignments, and EVERY reason is
# checked against the order probe.  An IMPOSSIBLE that the probe observes and
# an OCCURS that it never observes are both refused.
ADJUDICATION_REPAIRED = {}
_CONTACT = 'flip.make_contact_cq(live=False)'
_LIVE = 'flip.make_contact_cq(live=True)'
for _w in ('blast._cq_consistent_contact', _CONTACT, _LIVE):
    for _r in ('qca._carrier_density_consistent',
               'qca._quantum_capacitance_consistent',
               'blast._cq_consistent_contact', _CONTACT, _LIVE):
        if (_w, COX, _r) == (_CONTACT, COX, _CONTACT):
            ADJUDICATION_REPAIRED[(_w, COX, _r)] = (
                OCCURS,
                'the contact replacement substitutes gfet.C_ox for '
                'cd.C_interface and then reads it back within the same call, '
                'which is the substitution itself and not a leak: the write '
                'and the read are the same operating point, and the value is '
                'restored before the call returns.  This is the ONLY edge '
                'that survives the 2026-10-21 repair.')
        else:
            ADJUDICATION_REPAIRED[(_w, COX, _r)] = (
                IMPOSSIBLE,
                'under live_contact=False the channel replacements never '
                'arrive during the contact census (2026-10-22 arrival probe: '
                '0 arrivals, against 192 unrepaired), so no channel reader '
                'can observe a contact write of gfet.C_ox, and the live '
                'contact product is never installed.')
for _r in ('blast._cq_consistent_contact', _LIVE):
    ADJUDICATION_REPAIRED[('THE FLIP', CQ, _r)] = (
        GOVERNED,
        'gfet.quantum_capacitance is a member of its own identity class and '
        'is moved by the flip itself, so this edge is governed by '
        'bc.rebind_by_operating_point rather than by this relation.  Listed '
        'because it is an edge, not because it is unguarded.')


def _dotted(site):
    """'module.attr' for a (module, attr) pair.

    A name is a TUPLE throughout this module, and interpolating one into a
    single %s is a TypeError rather than a cosmetic problem.  V18 -- the
    planted control for a false GOVERNED -- crashed on exactly that, in the
    one refusal path no other check reaches.  Every message now goes through
    here."""
    return '%s.%s' % site


def _fmt(edge):
    w, name, r = edge
    return '%-38s --%s--> %s' % (w, _dotted(name), r)


# ---------------------------------------------------------------------------

def main():
    checks = []

    def check(name, ok, detail=''):
        checks.append((name, bool(ok), detail))
        print('    %-4s %s%s' % ('PASS' if ok else 'FAIL', name,
                                 ('   -- ' + detail) if detail else ''))

    print('=' * 78)
    print('THE WRITE/READ RELATION AS A GUARD, NOT A CENSUS')
    print('  discharges the 2026-10-22 ranked-first item')
    print('=' * 78)
    print()

    REPS = replacements_under_audit()
    METALS = flip.METALS
    base = {m: float(cd.contact_edge_carrier_density(wf)) for m, wf in METALS}

    # -- 0.  The probe is not the perturbation -----------------------------
    print('-' * 78)
    print('0.  THE PROBE IS NOT THE PERTURBATION')
    print('-' * 78)
    print()
    obs_f, trace_f, census_f = observe(False)
    obs_t, trace_t, census_t = observe(True)
    print('    7 committed contact densities, bitwise, outside the probe and')
    print('    inside it (the flip is entered inside, so the comparison is')
    print('    against the flip\'s own committed census, not the shipped one):')
    for m, _ in METALS:
        print('      %-3s  shipped %s   under probe(repaired) %s'
              % (m, float(base[m]).hex(), float(census_f[m]).hex()))
    check('V1  the probe terminates and returns a full 7-metal census',
          len(census_f) == 7 and len(census_t) == 7)
    rerun, _, census_f2 = observe(False)
    check('V2  the probe is deterministic -- census bitwise identical across '
          'two consecutive runs',
          {k: v.hex() for k, v in census_f.items()}
          == {k: v.hex() for k, v in census_f2.items()})
    check('V3  the observed edge set is deterministic too', rerun == obs_f,
          '%d edges' % len(obs_f))
    print()

    # -- 1.  Exactly-known values from 2026-10-22 --------------------------
    print('-' * 78)
    print('1.  REPRODUCING THE 2026-10-22 ARRIVAL COUNTS EXACTLY')
    print('-' * 78)
    print()
    tot_f, in_f, out_f, _ = rra.channel_arrival_probe(False)
    tot_t, in_t, out_t, _ = rra.channel_arrival_probe(True)
    print('    repaired   (live_contact=False): %3d arrivals, %3d with C_ox '
          'substituted' % (tot_f, in_f))
    print('    unrepaired (live_contact=True ): %3d arrivals, %3d with C_ox '
          'substituted' % (tot_t, in_t))
    check('V4  repaired flip gives EXACTLY ZERO channel arrivals', tot_f == 0)
    check('V5  unrepaired flip reproduces 192 arrivals exactly', tot_t == 192)
    check('V6  of which exactly 181 see C_ox already substituted', in_t == 181)
    print()
    print('    These are the 2026-10-22 numbers, reproduced rather than cited.')
    print('    They are the EXACTLY KNOWN values this module is anchored to:')
    print('    an exact zero that a symmetry requires, and two integers a')
    print('    previous session measured independently of this code.')
    print()

    # -- 2.  The relation ---------------------------------------------------
    print('-' * 78)
    print('2.  THE RELATION: STATIC OVER-APPROXIMATION, DYNAMIC TRUTH')
    print('-' * 78)
    print()
    flip_rebinds = set()
    for _, fac in rra.FLIPS:
        flip_rebinds |= rra.rebind_set(fac)
    wr = replacement_write_read(REPS)
    static = candidate_edges(wr, flip_rebinds)
    print('    writes and reads, per replacement, asked of reader C:')
    for label, (w, r) in sorted(wr.items()):
        print('      %s' % label)
        print('          writes %s' % (sorted('%s.%s' % x for x in w) or '--'))
    print()
    print('    STATIC relation (direct calls + the flip\'s rebind set):'
          '  %d edges' % len(static))
    print('    DYNAMIC, repaired   (live_contact=False): %d observed'
          % len(obs_f))
    for e in sorted(obs_f):
        print('        %s' % _fmt(e))
    print('    DYNAMIC, unrepaired (live_contact=True ): %d observed'
          % len(obs_t))
    for e in sorted(obs_t):
        print('        %s' % _fmt(e))
    print()
    check('V7  every observed edge is in the static relation (the static set '
          'is a sound over-approximation)',
          obs_f <= static and obs_t <= static,
          '%d outside' % len((obs_f | obs_t) - static))
    print()
    print('    THE TWO PREVIOUSLY-REPORTED DEFECTS ARE TWO EDGES OF ONE')
    print('    RELATION, and the repair removes both.  Unrepaired:')
    print('      * THE FLIP --gfet.quantum_capacitance--> make_contact_cq')
    print('        is 2026-10-21\'s defect: a replacement reading a name its')
    print('        own flip rebinds.')
    print('      * make_contact_cq --gfet.C_ox--> qca._quantum_capacitance_')
    print('        consistent is 2026-10-22 Section 4-5: the contact writes')
    print('        C_ox and a CHANNEL replacement reads it.  C_ox is in no')
    print('        identity class, so no existing mechanism can see this one.')
    print('    Repaired, neither is observed; the only survivor is the')
    print('    contact\'s write-then-read of C_ox WITHIN its own call, which')
    print('    is the substitution and not a leak.')
    cross_f = {e for e in obs_f if e[0] != e[2]}
    cross_t = {e for e in obs_t if e[0] != e[2]}
    check('V8  EXACT ZERO: the repaired composition has no cross-replacement '
          'edge', len(cross_f) == 0, '%d cross edges' % len(cross_f))
    check('V9  the unrepaired composition has exactly 2, and they are the two '
          'reported defects', len(cross_t) == 2,
          '%d cross edges' % len(cross_t))
    print()

    # -- 3.  The deriver, under exactly-known inputs ------------------------
    print('-' * 78)
    print('3.  THE DERIVER ITSELF, UNDER EXACTLY-KNOWN INPUTS')
    print('-' * 78)
    print()
    print('    A deriver that returns few edges and a repository that has few')
    print('    are indistinguishable until the deriver is asked a question')
    print('    whose answer is known without it.')
    print()
    N = ('m', 'x')
    pos = [('arrive', 'W'), ('write', N), ('return', 'W'),
           ('arrive', 'R'), ('read', N), ('return', 'R')]
    neg = [('arrive', 'R'), ('read', N), ('return', 'R'),
           ('arrive', 'W'), ('write', N), ('return', 'W')]
    nest = [('arrive', 'W'), ('write', N),
            ('arrive', 'R'), ('read', N), ('return', 'R'),
            ('return', 'W')]
    unattr = [('write', N), ('arrive', 'R'), ('read', N), ('return', 'R')]
    check('V10 POSITIVE CONTROL: write then read yields exactly the one edge',
          edges_from_trace(pos) == {('W', N, 'R')},
          repr(sorted(edges_from_trace(pos))))
    check('V11 NEGATIVE CONTROL: read BEFORE the write yields no edge -- the '
          'deriver is ordered, not a set intersection',
          edges_from_trace(neg) == set())
    check('V12 a nested reader is attributed to itself, not to its caller',
          edges_from_trace(nest) == {('W', N, 'R')})
    check('V13 a write outside every replacement is attributed to THE FLIP',
          edges_from_trace(unattr) == {('THE FLIP', N, 'R')})
    print()

    # -- 4.  The guard refuses ---------------------------------------------
    print('-' * 78)
    print('4.  THE GUARD REFUSES -- AND REFUSES IN BOTH DIRECTIONS')
    print('-' * 78)
    print()

    def guard(adj, live=False):
        with require_write_read_adjudication(live, REPS, adj):
            pass

    adj = dict(ADJUDICATION_REPAIRED)
    try:
        guard(adj)
        ok, why = True, 'accepted'
    except Exception as exc:                                  # noqa: BLE001
        ok, why = False, '%s: %s' % (type(exc).__name__, exc)
    check('V14 the complete, correct adjudication of the REPAIRED '
          'composition is accepted', ok, why if not ok else '')

    dropped = dict(adj)
    victim = sorted(dropped)[0]
    del dropped[victim]
    try:
        guard(dropped)
        ok = False
        why = 'accepted an adjudication missing %s' % (victim,)
    except UnadjudicatedWriteReadEdge as exc:
        ok, why = True, str(exc).splitlines()[0][:70]
    except Exception as exc:                                  # noqa: BLE001
        ok, why = False, 'wrong exception %s' % type(exc).__name__
    check('V15 an adjudication missing ONE edge is refused before anything '
          'runs', ok, why)

    observed_edge = sorted(obs_f)[0]
    lied = dict(adj)
    lied[observed_edge] = (IMPOSSIBLE, 'a claim the probe can falsify')
    try:
        guard(lied)
        ok, why = False, 'accepted a false IMPOSSIBLE'
    except FalseEdgeVerdict as exc:
        ok, why = True, str(exc).splitlines()[0][:70]
    except Exception as exc:                                  # noqa: BLE001
        ok, why = False, 'wrong exception %s' % type(exc).__name__
    check('V16 IMPOSSIBLE claimed on an edge the probe OBSERVES is refused',
          ok, why)

    unobserved = sorted(static - obs_f - {observed_edge})
    stale_edge = next(e for e in unobserved if e[0] != 'THE FLIP')
    stale = dict(adj)
    stale[stale_edge] = (OCCURS, 'a reason describing nothing')
    try:
        guard(stale)
        ok, why = False, 'accepted a stale OCCURS'
    except FalseEdgeVerdict as exc:
        ok, why = True, str(exc).splitlines()[0][:70]
    except Exception as exc:                                  # noqa: BLE001
        ok, why = False, 'wrong exception %s' % type(exc).__name__
    check('V17 OCCURS claimed on an edge the probe NEVER observes is refused '
          '-- the half a census cannot have', ok, why)

    gov = dict(adj)
    gov[observed_edge] = (GOVERNED, 'claims an existing mechanism moves C_ox')
    try:
        guard(gov)
        ok, why = False, 'accepted a false GOVERNED'
    except FalseEdgeVerdict as exc:
        ok, why = True, str(exc).splitlines()[0][:70]
    except Exception as exc:                                  # noqa: BLE001
        ok, why = False, 'wrong exception %s' % type(exc).__name__
    check('V18 GOVERNED claimed on a name in no rebind set is refused',
          ok, why)

    noreason = dict(adj)
    noreason[observed_edge] = (OCCURS, '')
    try:
        guard(noreason)
        ok, why = False, 'accepted an empty reason'
    except UnadjudicatedWriteReadEdge:
        ok, why = True, ''
    except Exception as exc:                                  # noqa: BLE001
        ok, why = False, 'wrong exception %s' % type(exc).__name__
    check('V19 an edge with a verdict and NO stated reason is refused -- the '
          'bc.rebind() refusal, applied to an edge', ok, why)

    badverdict = dict(adj)
    badverdict[observed_edge] = ('PROBABLY_FINE', 'not in the vocabulary')
    try:
        guard(badverdict)
        ok, why = False, 'accepted a verdict outside the vocabulary'
    except UnadjudicatedWriteReadEdge:
        ok, why = True, ''
    except Exception as exc:                                  # noqa: BLE001
        ok, why = False, 'wrong exception %s' % type(exc).__name__
    check('V20 a verdict outside {%s} is refused' % '/'.join(VERDICTS),
          ok, why)

    try:
        guard(adj, live=True)
        ok, why = False, ('the repaired adjudication was accepted for the '
                          'UNREPAIRED composition')
    except FalseEdgeVerdict as exc:
        ok, why = True, str(exc).splitlines()[0][:70]
    except Exception as exc:                                  # noqa: BLE001
        ok, why = False, 'wrong exception %s' % type(exc).__name__
    check('V21 THE RESULT: the repaired composition\'s adjudication is '
          'REFUSED for the unrepaired one', ok, why)
    print()
    print('    V21 is what the item asked for.  The same set of stated')
    print('    reasons is accepted for live_contact=False and refused for')
    print('    live_contact=True, so the difference between the repaired and')
    print('    unrepaired compositions is now a RAISE rather than a')
    print('    percentage table that a reader has to interpret.')
    print()

    print('=' * 78)
    npass = sum(1 for _, ok, _ in checks if ok)
    print('VALIDATION: %d/%d' % (npass, len(checks)))
    for name, ok, detail in checks:
        if not ok:
            print('   FAILED: %s   %s' % (name, detail))
    print('=' * 78)
    return 0 if npass == len(checks) else 1


if __name__ == '__main__':
    import sys
    sys.exit(main())
