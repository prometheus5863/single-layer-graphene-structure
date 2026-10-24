"""graphene_capture_surface_guard.py -- a capture is adjudicated by a
COUNTERFACTUAL, not by a read set.

WHY THIS EXISTS.

2026-10-22 Section 6 found that a capture protects the NAME IT HOLDS and not
the names the HELD OBJECT reads at call time, and recorded that `4.6.17`'s
repair is therefore "correct and contingent".  2026-10-23 built a guard over
the write/read relation and ranked this item FIRST on the ground that a
captured object's own live reads are reached by neither that guard -- which
only sees names some replacement WRITES -- nor by `bc.*`, which only sees
names a flip REBINDS.  A name that is merely READ, through a captured object,
by a replacement whose own source never mentions it, is invisible to both.

THE FIRST THING MEASURED HERE IS THAT THE OBVIOUS INSTRUMENT CANNOT DO IT.

The cheapest entry the standing item suggested was `rra.static_reads` on each
captured object.  Run on the two contact replacements, that instrument
reports:

    static live-read set of flip.make_contact_cq(live=False)  : 8 names
    static live-read set of flip.make_contact_cq(live=True )  : the SAME 8

bitwise the same set, because `cq_of_n` contains BOTH branches in its
bytecode and a bytecode walk cannot know that the `live=False` branch never
takes the one that reads `gfet.quantum_capacitance`.  The repaired and
unrepaired contact replacements -- the two objects the entire 2026-10-21
repair is the difference between -- are INDISTINGUISHABLE to the static
reader.  This is not a bug in `static_reads`; it is sound and
over-approximate by construction, and `4.6.17`'s own correctness argument
is written in its vocabulary.  It means the argument was never checkable in
that vocabulary.

SO THE VERDICT IS A COUNTERFACTUAL.  For a capture and a name, rebind the
name to a perturbed object, call the holder, and compare BITWISE:

    PROTECTED  the rebind leaves the holder's value bitwise unchanged.
               Refused if any rebind MOVES it.
    EXPOSED    the rebind moves it, and the reason says why that is sound.
               Refused if no rebind moves it -- a reason describing nothing,
               the half a census cannot have (2026-10-23).
    UNREACHED  the rebind leaves it bitwise unchanged for a DIFFERENT
               reason: the holder never looks the name up at all, so the
               exact zero is not the capture's doing.  Refused if the name
               IS on the read surface.
    GOVERNED   the name is in some flip's rebind set, so an existing
               mechanism moves it deliberately.  Refused if the name is in
               no rebind set.

WHY `UNREACHED` IS A SEPARATE WORD, and the one job the static reader keeps.
Walking the holder's closure cells by OBJECT IDENTITY recovers not one held
name but the held object's whole IDENTITY CLASS -- here
{gfet.quantum_capacitance, cd.quantum_capacitance}, reproducing
`bc.identity_class(gfet, 'quantum_capacitance')` exactly and independently
(C7).  So the capture neutralises a rebind of EITHER member, which is a
STRONGER protection than `4.6.17` claimed.  But only one member is a name
the holder ever looks up; a rebind of the other is bitwise-neutral whether
or not anything was captured.  Reporting both as PROTECTED would credit the
capture with work it did not do -- the same conflation as `4.6.17`, with the
sign reversed.  Separating them needs a judgement the counterfactual cannot
make, and it is exactly the one judgement `static_reads` CAN make soundly:
an over-approximate read set has reliable NEGATIVES, so "not in the static
read surface" is trustworthy even though "in it" is not.  Section 0's
instrument is kept for that and for nothing else.

AND ONE REFUSAL THAT IS THE RECORDED ERROR ITSELF.  `4.6.17` concluded
"captured, therefore protected".  Here, a `PROTECTED` verdict on a HELD name
obliges every name the object held there reads live to carry its own verdict,
and any of those that the counterfactual MOVES must be adjudicated `EXPOSED`
with a stated reason.  A `PROTECTED` on the held name sitting next to a
missing or `PROTECTED` verdict on a read that moves is exactly the 2026-10-22
finding, and it is now a `raise` rather than a sentence in a log.

WHAT THE COUNTERFACTUAL FINDS, stated here because it is the result.

For the REPAIRED contact replacement the capture is worth exactly one name
and that name has an exact twin it does not cover:

    gfet.quantum_capacitance   PROTECTED  (exact 0; the 2026-10-21 repair)
    gfet._thermal_rounding     EXPOSED    (+58.34 %)

and the perturbed value is BITWISE THE SAME NUMBER in both cases --
`0x1.a92b527d6852cp-4` -- because `gfet._thermal_rounding` is the single live
read of the captured object itself.  The repair closed one route to that
value and left an equally powerful one open THROUGH THE OBJECT IT CAPTURED.
`gfet._thermal_rounding` is in no rebind set and is written by no
replacement, so it is reached by `bc.rebind_by_operating_point`, by the
2026-10-23 write/read guard, and by `4.6.17`'s argument equally: not at all.

NEUTRALITY AND EXACTNESS.  Every verdict rests on a bitwise comparison, so
the probe must not perturb what it measures: an IDENTITY rebind -- `setattr`
of the object already there -- is run for every name on every surface and
must move the value by EXACTLY ZERO (C1).  That is the symmetry this module
is anchored to, in the 2026-10-17 sense, and it is checked before any verdict
is used.

Usage:  PYTHONPATH=... python3 graphene_capture_surface_guard.py
"""

import contextlib
import types

import numpy as np

import graphene_fet_model as gfet
import graphene_binding_class as bc
import graphene_contact_doping_model as cd
import graphene_quantum_capacitance_argument as qca
import graphene_default_flip as flip
import graphene_replacement_read_audit as rra
import graphene_write_read_guard as wrg


PROTECTED = 'PROTECTED'
UNREACHED = 'UNREACHED'
EXPOSED = 'EXPOSED'
GOVERNED = 'GOVERNED'
VERDICTS = (PROTECTED, UNREACHED, EXPOSED, GOVERNED)

PROBE_V = 0.5
PERTURBATION = 1.5


class UnadjudicatedCaptureSurface(RuntimeError):
    """A (capture, name) pair on the surface was given neither a verdict nor a
    stated reason, or a `PROTECTED` held name has an unadjudicated live read.
    Raised BEFORE anything is measured, like `bc.UnassignedClassMember`."""


class FalseCaptureVerdict(RuntimeError):
    """A verdict disagrees with the counterfactual."""


def _dotted(site):
    """'module.attr' for a (module, attr) pair.

    A name is a TUPLE throughout this module.  2026-10-23's V18 found a
    `TypeError` in the one refusal branch nothing else reached, caused by
    interpolating such a pair into a single %s; every message here goes
    through this function for that reason."""
    return '%s.%s' % site


# ---------------------------------------------------------------------------
# What a holder HOLDS, and what the held objects READ
# ---------------------------------------------------------------------------

def held_objects(fn):
    """Every callable reachable from `fn`'s closure cells, transitively.

    Derived from the cells, not from a list: `flip.make_contact_cq`'s product
    holds `cq_of_n`, which in turn holds `shipped_cq`, and only the second is
    the captured module attribute.  A factory that grew a third cell would be
    reported with three without this module being edited -- the property the
    hand-maintained lists of 2026-10 did not have.
    """
    out, stack, seen = [], [fn], {id(fn)}
    while stack:
        cur = stack.pop()
        for val in rra._closure_map(cur).values():
            if not callable(val) or id(val) in seen:
                continue
            seen.add(id(val))
            out.append(val)
            if isinstance(val, types.FunctionType):
                stack.append(val)
    return out


def held_names(fn):
    """{(module, attr): object} for held objects that ARE a module attribute.

    This is the set a capture genuinely freezes: the identity `fn` holds is
    the object that name pointed at when the cell was filled, so a later
    rebind of the name cannot reach it.  A held object that is NOT a module
    attribute -- `cq_of_n`, built fresh per call -- freezes nothing nameable
    and is therefore not a held NAME, only a held object.
    """
    mods = rra.repo_modules()
    out = {}
    for obj in held_objects(fn):
        for mname, mod in mods.items():
            for attr, val in vars(mod).items():
                if val is obj and not attr.startswith('__'):
                    out[(mname, attr)] = obj
    return out


def capture_surface(fn):
    """(surface, read_surface, held, reads_of_held) for a holder `fn`.

    read_surface:   names `fn` or anything it holds LOOKS UP -- its own live
                    reads together with the live reads of every held object.
                    Over-approximate by construction and deliberately so: its
                    negatives are sound, which is all `UNREACHED` needs, and
                    its positives decide nothing (Section 0).
    surface:        read_surface together with the HELD NAMES.  A held name
                    belongs on the surface even when it is never looked up,
                    because the capture makes a claim about it.
    held:           {(module, attr): object} the capture freezes by identity;
                    the held object's whole identity class, not one name.
    reads_of_held:  {object id: ({held names}, {live reads})} -- the surface
                    the capture does NOT cover, keyed by the held OBJECT,
                    because the reads belong to the object and not to
                    whichever alias it was found under.
    """
    live, _ = rra.static_reads(fn)
    read_surface = set(live)
    held = held_names(fn)
    by_obj = {}
    for site, obj in held.items():
        names, reads = by_obj.setdefault(id(obj), (set(), set()))
        names.add(site)
        if not reads:
            reads |= set(rra.static_reads(obj)[0])
        read_surface |= reads
    return read_surface | set(held), read_surface, held, by_obj


# ---------------------------------------------------------------------------
# The counterfactual
# ---------------------------------------------------------------------------

def perturb(obj, k=PERTURBATION):
    """A deliberately wrong stand-in for `obj`, same shape.

    k=1.0 on a callable still returns a WRAPPER rather than the object, so the
    identity control (C1) must use `rebind_identity` instead: a wrapper that
    scales by 1.0 would hide a probe that perturbs by wrapping.
    """
    if isinstance(obj, bool):
        raise TypeError('no perturbation defined for a bool')
    if isinstance(obj, (int, float)):
        return type(obj)(obj * k) if obj else type(obj)(k)
    if callable(obj):
        def wrapper(*a, **kw):
            return k * np.asarray(obj(*a, **kw), dtype=float)
        return wrapper
    raise TypeError('no perturbation defined for %r' % type(obj).__name__)


@contextlib.contextmanager
def rebound(site, factory):
    """Temporarily put `factory(old)` at `site`, then restore it."""
    mod = rra.repo_modules()[site[0]]
    old = getattr(mod, site[1])
    setattr(mod, site[1], factory(old))
    try:
        yield old
    finally:
        setattr(mod, site[1], old)


def leverage(holder, site, k=PERTURBATION, V=PROBE_V):
    """(base, moved, exact) -- the holder's value before and under the rebind.

    `exact` is True when the two are BITWISE equal.  Nothing here is compared
    to a tolerance: the question "does a rebind of this name reach this
    object" has an exact answer, and asking it with a tolerance is the
    2026-10-09 tolerance item in a new place.
    """
    base = float(holder(V))
    with rebound(site, lambda old: perturb(old, k)):
        moved = float(holder(V))
    return base, moved, base.hex() == moved.hex()


def identity_leverage(holder, site, V=PROBE_V):
    """The same measurement with the SAME OBJECT put back -- must be exact."""
    base = float(holder(V))
    with rebound(site, lambda old: old):
        again = float(holder(V))
    return base, again, base.hex() == again.hex()


# ---------------------------------------------------------------------------
# THE GUARD
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def require_capture_adjudication(label, holder, adjudication, report=None):
    """Refuse a capture whose surface is not fully adjudicated.

    adjudication: {(module, attr): (verdict, reason)}, verdict in VERDICTS and
        the reason mandatory, in the shape `bc.rebind_by_operating_point` uses
        for assignments and `wrg.require_write_read_adjudication` uses for
        edges.

    Raises `UnadjudicatedCaptureSurface` BEFORE measuring anything if a
    surface name is unadjudicated, if a spec is malformed, or if a PROTECTED
    held name has an unadjudicated live read; and `FalseCaptureVerdict` after
    measuring if a verdict disagrees with the counterfactual -- in either
    direction.
    """
    surface, read_surface, held, reads_of_held = capture_surface(holder)
    flip_rebinds = set()
    for _, fac in rra.FLIPS:
        flip_rebinds |= rra.rebind_set(fac)

    for site in sorted(set(surface) & set(adjudication)):
        spec = adjudication[site]
        if not (isinstance(spec, tuple) and len(spec) == 2):
            raise UnadjudicatedCaptureSurface(
                '%s: adjudication for %s must be a (verdict, reason) pair'
                % (label, _dotted(site)))
        verdict, reason = spec
        if verdict not in VERDICTS:
            raise UnadjudicatedCaptureSurface(
                '%s: adjudication for %s has verdict %r, not one of %s'
                % (label, _dotted(site), verdict, '/'.join(VERDICTS)))
        if not reason:
            raise UnadjudicatedCaptureSurface(
                '%s: adjudication for %s has no stated reason'
                % (label, _dotted(site)))

    # THE 2026-10-22 FINDING, AS A REFUSAL.  "Captured, therefore protected"
    # is only an argument about the held NAME.  Every live read of the object
    # held there is a name the capture does not cover, and must be answered.
    # This runs BEFORE the generic missing-name check so that the message a
    # reader gets names the mechanism rather than just the gap.
    for names, reads in sorted(reads_of_held.values()):
        claimed = sorted(n for n in names
                         if adjudication.get(n, (None,))[0] == PROTECTED)
        if not claimed:
            continue
        unanswered = sorted(r for r in reads if r not in adjudication)
        if unanswered:
            raise UnadjudicatedCaptureSurface(
                '%s: %s is adjudicated PROTECTED, but the object held there '
                'reads %s live and %s carries no verdict.  A capture protects '
                'the name it holds, not the names the held object reads at '
                'call time (2026-10-22 Section 6).'
                % (label, ', '.join(_dotted(c) for c in claimed),
                   ', '.join(_dotted(r) for r in sorted(reads)),
                   ', '.join(_dotted(r) for r in unanswered)))

    missing = sorted(set(surface) - set(adjudication))
    if missing:
        raise UnadjudicatedCaptureSurface(
            '%s: %d name(s) on the capture surface are unadjudicated: %s.  '
            'Give each a (verdict, reason) with verdict in %s; nothing has '
            'been measured.'
            % (label, len(missing), ', '.join(_dotted(m) for m in missing),
               '/'.join(VERDICTS)))

    # -- measure -----------------------------------------------------------
    rows, wrong = [], []
    for site in sorted(surface):
        verdict, reason = adjudication[site]
        base, ident, ident_exact = identity_leverage(holder, site)
        if not ident_exact:
            wrong.append((site, verdict,
                          'the IDENTITY rebind of %s moved the value (%s -> '
                          '%s): the probe perturbs what it measures and no '
                          'verdict on this surface can be trusted'
                          % (_dotted(site), base.hex(), ident.hex())))
            continue
        base, moved, exact = leverage(holder, site)
        rows.append((site, verdict, base, moved, exact))
        if verdict == PROTECTED and not exact:
            wrong.append((site, verdict,
                          'claimed PROTECTED but a rebind of %s moved the '
                          'value by %+.4f %%'
                          % (_dotted(site), 100.0 * (moved - base) / base)))
        elif verdict == PROTECTED and site not in read_surface:
            wrong.append((site, verdict,
                          'claimed PROTECTED but %s is not on the read '
                          'surface at all, so the exact zero is not the '
                          "capture's doing -- the verdict is UNREACHED"
                          % _dotted(site)))
        elif verdict == UNREACHED and not exact:
            wrong.append((site, verdict,
                          'claimed UNREACHED but a rebind of %s moved the '
                          'value by %+.4f %%'
                          % (_dotted(site), 100.0 * (moved - base) / base)))
        elif verdict == UNREACHED and site in read_surface:
            wrong.append((site, verdict,
                          'claimed UNREACHED but %s IS on the read surface'
                          % _dotted(site)))
        elif verdict == EXPOSED and exact:
            wrong.append((site, verdict,
                          'claimed EXPOSED but a rebind of %s left the value '
                          'bitwise unchanged -- the stated reason describes '
                          'nothing' % _dotted(site)))
        elif verdict == GOVERNED and site not in flip_rebinds:
            wrong.append((site, verdict,
                          'claimed GOVERNED but %s is in no rebind set'
                          % _dotted(site)))

    # A PROTECTED held name whose held object's reads MOVE and are themselves
    # called PROTECTED is the recorded error with both halves present.
    for names, reads in sorted(reads_of_held.values()):
        claimed = sorted(n for n in names
                         if adjudication.get(n, (None,))[0] == PROTECTED)
        if not claimed:
            continue
        for r in sorted(reads):
            v = adjudication[r][0]
            hit = [row for row in rows if row[0] == r]
            if v in (PROTECTED, UNREACHED) and hit and not hit[0][4]:
                wrong.append((r, v,
                              'is a live read of the object held at %s, is '
                              'called %s, and a rebind of it moves the '
                              'value: the capture covers the held name and '
                              'not this' % (_dotted(claimed[0]), v)))
    if wrong:
        raise FalseCaptureVerdict(
            '%s: %d verdict(s) disagree with the counterfactual:\n%s'
            % (label, len(wrong), '\n'.join(
                '    %-46s %-10s %s' % (_dotted(s), v, why)
                for s, v, why in wrong)))
    if report is not None:
        report(label, surface, held, reads_of_held, rows)
    yield surface, rows


# ---------------------------------------------------------------------------
# Planted controls -- one per verdict AND one per refusal branch (2026-10-23)
# ---------------------------------------------------------------------------

def control_exposed():
    """Holds nothing; reads `gfet.T` live.  Must measure EXPOSED on gfet.T."""
    def holder(V):
        return float(gfet.T) * float(V)
    return holder


def control_protected():
    """Reads no module attribute at all.  Must measure PROTECTED on every
    name.  The negative half that makes the positive control mean something:
    a probe that reported movement unconditionally would pass the other."""
    def holder(V):
        return 2.0 * float(V)
    return holder


def control_capture_hides_a_read():
    """The 2026-10-22 shape, synthesised.

    Holds `gfet.quantum_capacitance` by identity in a cell -- so that name is
    genuinely PROTECTED -- while the held object's own live read of
    `gfet._thermal_rounding` reaches straight through.  A guard that accepts
    PROTECTED on the held name without answering for that read has reproduced
    the recorded error on a case whose right answer is known in advance.
    """
    captured = gfet.quantum_capacitance

    def holder(V):
        return float(captured(float(V), T=300.0))
    return holder


# ---------------------------------------------------------------------------
# The captures under audit, and their adjudication
# ---------------------------------------------------------------------------

CI = ('graphene_contact_doping_model', 'C_interface')
COX = ('graphene_fet_model', 'C_ox')
T = ('graphene_fet_model', 'T')
TR = ('graphene_fet_model', '_thermal_rounding')
QC = ('graphene_fet_model', 'quantum_capacitance')
CDQC = ('graphene_contact_doping_model', 'quantum_capacitance')
AF = ('graphene_quantum_capacitance_argument', 'A_FERMI')
FP = ('graphene_quantum_capacitance_argument', '_fixed_point')
SF = ('graphene_quantum_capacitance_argument', 'series_factor')

_PHYSICAL = ('a model ingredient of the contact self-consistency itself, not '
             'a name any flip rebinds or any replacement writes; the contact '
             'solve is MEANT to depend on it, and the dependence is recorded '
             'here so that a future flip over this name is a refusal rather '
             'than a silent change of Chapter 4 numbers.')

ADJUDICATION_REPAIRED = {
    CDQC: (UNREACHED,
           'the OTHER member of the held object\'s identity class, found by '
           'walking the closure cells on object identity and agreeing with '
           'bc.identity_class(gfet, "quantum_capacitance") exactly.  The '
           'capture neutralises a rebind of this member too -- which is a '
           'STRONGER protection than 4.6.17 claimed -- but the holder never '
           'looks this name up, so the bitwise zero here is not the '
           "capture's doing and calling it PROTECTED would credit the "
           'capture with work it did not do.  It is in '
           'blast.contact_on_consistent_argument()\'s rebind set, so GOVERNED '
           'would also be a true sentence about it; UNREACHED is the stronger '
           'one, because it is the reason the rebind is harmless HERE.'),
    QC: (PROTECTED,
         'the 2026-10-21 repair, re-measured as a counterfactual rather than '
         'cited: make_contact_cq(live=False) holds the object this name '
         'pointed at before the flip, so a rebind of the name cannot reach '
         'it and the value is BITWISE unchanged.  This is the only name on '
         'which the repaired and unrepaired replacements differ, and the '
         'static read set does not distinguish them on it.'),
    TR: (EXPOSED,
         'THE FINDING.  gfet._thermal_rounding is the single live read of '
         'the CAPTURED object, so the capture above does not cover it: a '
         'rebind of it moves the repaired replacement to 0x1.a92b527d6852cp-4, '
         'which is BITWISE the value a rebind of gfet.quantum_capacitance '
         'moves the UNREPAIRED one to.  The repair bought protection against '
         'one name and left an exactly equally powerful one open through the '
         'object it captured.  In no rebind set, written by no replacement, '
         'so reached by neither bc.rebind_by_operating_point nor the '
         '2026-10-23 write/read guard.'),
    COX: (PROTECTED,
          'EXACT ZERO, and it is the quantitative form of the one OCCURS '
          'edge the 2026-10-23 guard left standing.  That edge was '
          'adjudicated "the substitution itself and not a leak" on the '
          'argument that the write and the read are the same operating point; '
          'here the external value of gfet.C_ox is shown to reach the '
          "replacement's output not at all, because the replacement "
          'substitutes cd.C_interface for the whole solve and restores the '
          'old value before returning.  An argument yesterday, a bitwise '
          'zero today.'),
    CI: (EXPOSED, 'C_interface is what the contact solve is referenced TO; '
                  + _PHYSICAL),
    T: (EXPOSED, 'the lattice temperature enters the inner C_q through '
                 'cq_of_n(n, T=None); ' + _PHYSICAL),
    AF: (EXPOSED, 'the Fermi-energy prefactor E_F = A_FERMI*sqrt(n); '
                  + _PHYSICAL),
    FP: (EXPOSED, 'the committed form-B_T fixed-point solver, shared with the '
                  'channel so the two cannot drift apart; ' + _PHYSICAL),
    SF: (EXPOSED, 'the series factor C_q/(C_q+C) the contact solve closes on; '
                  + _PHYSICAL),
}


def _fmt_row(site, verdict, base, moved, exact):
    pct = 0.0 if exact else 100.0 * (moved - base) / base
    return ('    %-48s %-10s %s %s %s'
            % (_dotted(site), verdict, base.hex(), moved.hex(),
               'exact 0' if exact else '%+8.4f %%' % pct))


def main():
    checks = []

    def check(name, ok, detail=''):
        detail = detail if isinstance(detail, str) else repr(detail)
        checks.append((name, bool(ok), detail))
        print('    %-4s %s%s' % ('PASS' if ok else 'FAIL', name,
                                 ('   -- ' + detail) if detail else ''))

    print('=' * 78)
    print('THE CAPTURE SURFACE, ADJUDICATED BY COUNTERFACTUAL')
    print('  discharges the 2026-10-23 ranked-first item')
    print('=' * 78)
    print()

    shipped = gfet.quantum_capacitance
    rep = flip.make_contact_cq(shipped, live=False)
    unrep = flip.make_contact_cq(shipped, live=True)

    # -- 0.  The static reader cannot do this ------------------------------
    print('-' * 78)
    print('0.  THE OBVIOUS INSTRUMENT HAS ZERO DISCRIMINATING POWER')
    print('-' * 78)
    print()
    sr_r, _ = rra.static_reads(rep)
    sr_u, _ = rra.static_reads(unrep)
    print('    static live reads of the REPAIRED   replacement: %d' % len(sr_r))
    print('    static live reads of the UNREPAIRED replacement: %d' % len(sr_u))
    for s in sorted(sr_r | sr_u):
        print('        %-48s repaired=%s unrepaired=%s'
              % (_dotted(s), s in sr_r, s in sr_u))
    check('C0  the static read sets of the repaired and unrepaired contact '
          'replacements are IDENTICAL -- the instrument the standing item '
          'suggested cannot tell them apart', sr_r == sr_u,
          '%d names, symmetric difference %d'
          % (len(sr_r), len(sr_r ^ sr_u)))
    print()
    print('    Both contain gfet.quantum_capacitance as a LIVE read, because')
    print('    cq_of_n holds both branches in its bytecode.  Sound and')
    print('    over-approximate; useless for the question.  Every verdict')
    print('    below is a counterfactual instead.')
    print()

    # -- 1.  The probe is not the perturbation -----------------------------
    print('-' * 78)
    print('1.  THE PROBE IS NOT THE PERTURBATION (exact, before any verdict)')
    print('-' * 78)
    print()
    surface, read_surface, held, reads_of_held = capture_surface(rep)
    base_hex = float(rep(PROBE_V)).hex()
    check('C1  EXACT ZERO: an IDENTITY rebind of every one of the %d surface '
          'names moves the value by exactly 0' % len(surface),
          all(identity_leverage(rep, s)[2] for s in surface),
          'base %s' % base_hex)
    check('C2  the holder is deterministic -- bitwise identical on two '
          'consecutive calls',
          float(rep(PROBE_V)).hex() == float(rep(PROBE_V)).hex())
    check('C3  the repaired and unrepaired replacements agree bitwise when '
          'nothing is rebound (so every difference below is the rebind)',
          float(rep(PROBE_V)).hex() == float(unrep(PROBE_V)).hex(), base_hex)
    print()
    print('    held by identity (the capture): %s'
          % (sorted(_dotted(h) for h in held) or '--'))
    for names, reads in sorted(reads_of_held.values()):
        print('    live reads OF the object held at %s: %s'
              % (', '.join(_dotted(n) for n in sorted(names)),
                 sorted(_dotted(r) for r in reads)))
    print('    read surface %d names, full surface %d (the extra one is a '
          'held name nothing looks up)' % (len(read_surface), len(surface)))
    print()

    # -- 2.  The counterfactual --------------------------------------------
    print('-' * 78)
    print('2.  THE COUNTERFACTUAL, PER NAME (perturbation x%.1f)' % PERTURBATION)
    print('-' * 78)
    print()
    print('    %-48s %-10s %-22s %-22s %s'
          % ('name', 'verdict', 'base', 'rebound', 'move'))
    rows_r, rows_u = [], []
    for s in sorted(surface):
        b, m, e = leverage(rep, s)
        rows_r.append((s, ADJUDICATION_REPAIRED[s][0], b, m, e))
        print(_fmt_row(s, ADJUDICATION_REPAIRED[s][0], b, m, e))
        rows_u.append((s,) + leverage(unrep, s))
    print()
    diff = [s for (s, b, m, e), (_, vb, vm, ve)
            in zip(rows_u, [(r[0], r[2], r[3], r[4]) for r in rows_r])
            if e != ve]
    check('C4  the repaired and unrepaired replacements differ on EXACTLY ONE '
          'surface name, and it is gfet.quantum_capacitance',
          diff == [QC], '%s' % [_dotted(d) for d in diff])
    prot = {s for s, _, _, _, e in rows_r if e}
    check('C5  exactly 3 of the %d surface names are bitwise unmoved, and '
          'they are the two the capture earns plus the alias nothing looks '
          'up' % len(surface), prot == {QC, COX, CDQC},
          sorted(_dotted(p) for p in prot))
    check('C5b VALIDATION AGAINST THE EXISTING IMPLEMENTATION: walking the '
          'closure cells on object identity reproduces bc.identity_class '
          'exactly -- two independent derivations of the same set',
          set(held) == set(bc.identity_class(gfet, 'quantum_capacitance')),
          sorted(_dotted(h) for h in held))
    lev_tr = leverage(rep, TR)[1].hex()
    lev_qc_unrep = leverage(unrep, QC)[1].hex()
    check('C6  THE FINDING, BITWISE: a rebind of gfet._thermal_rounding moves '
          'the REPAIRED replacement to exactly the value a rebind of '
          'gfet.quantum_capacitance moves the UNREPAIRED one to',
          lev_tr == lev_qc_unrep, '%s' % lev_tr)
    held_reads = set()
    for _, reads in reads_of_held.values():
        held_reads |= reads
    check('C7  gfet._thermal_rounding is the SOLE live read of the captured '
          'object -- so the route C6 measures is the capture itself',
          held_reads == {TR}, sorted(_dotted(r) for r in held_reads))
    flip_rebinds = set()
    for _, fac in rra.FLIPS:
        flip_rebinds |= rra.rebind_set(fac)
    check('C8  and it is unguarded on BOTH existing axes: in no flip rebind '
          'set, and not a written name in the 2026-10-23 relation',
          TR not in flip_rebinds
          and TR not in {e[1] for e in wrg.observe(False)[0]},
          'rebind sets hold %d names' % len(flip_rebinds))
    check('C9  gfet.C_ox is EXACTLY protected for BOTH replacements -- the '
          'quantitative form of the one OCCURS edge the write/read guard '
          'left standing', leverage(rep, COX)[2] and leverage(unrep, COX)[2])
    print()

    # -- 3.  Planted controls, one per verdict -----------------------------
    print('-' * 78)
    print('3.  PLANTED CONTROLS: ONE PER VERDICT, ON KNOWN ANSWERS')
    print('-' * 78)
    print()
    ce, cp = control_exposed(), control_protected()
    check('C10 POSITIVE: a holder that reads gfet.T live MOVES under a rebind '
          'of it', not leverage(ce, T)[2])
    check('C11 NEGATIVE: a holder that reads no module attribute is bitwise '
          'unchanged under a rebind of the same name',
          leverage(cp, T)[2])
    check('C12 a name the holder does not read at all is bitwise unchanged -- '
          'the probe does not report movement it cannot cause',
          leverage(ce, AF)[2])
    ch = control_capture_hides_a_read()
    check('C13 the synthetic 2026-10-22 shape: the HELD name is protected and '
          'the held object\'s own read is not',
          leverage(ch, QC)[2] and not leverage(ch, TR)[2])
    check('C13b the alias member is bitwise unmoved for the synthetic shape '
          'too, and is likewise not on its read surface -- UNREACHED is a '
          'property of the holder, not of this one capture',
          leverage(ch, CDQC)[2]
          and CDQC not in capture_surface(ch)[1])
    print()

    # -- 4.  The guard refuses, one control per refusal branch -------------
    print('-' * 78)
    print('4.  THE GUARD REFUSES -- A PLANTED CONTROL FOR EVERY BRANCH')
    print('-' * 78)
    print()

    def attempt(adj, holder=rep, label='repaired contact capture'):
        try:
            with require_capture_adjudication(label, holder, adj):
                pass
            return None
        except Exception as exc:                              # noqa: BLE001
            return exc

    def expect(name, adj, exc_type, holder=rep):
        exc = attempt(adj, holder)
        if exc_type is None:
            check(name, exc is None,
                  '' if exc is None else '%s: %s' % (type(exc).__name__, exc))
            return
        ok = isinstance(exc, exc_type)
        check(name, ok,
              'no exception' if exc is None
              else (str(exc).splitlines()[0][:72] if ok
                    else 'wrong exception %s' % type(exc).__name__))

    adj = dict(ADJUDICATION_REPAIRED)
    expect('C14 the complete, correct adjudication of the REPAIRED capture is '
           'accepted', adj, None)

    dropped = dict(adj)
    del dropped[sorted(dropped)[0]]
    expect('C15 an adjudication missing ONE surface name is refused before '
           'anything is measured', dropped, UnadjudicatedCaptureSurface)

    lie_p = dict(adj)
    lie_p[TR] = (PROTECTED, 'claims the capture covers the held object\'s '
                            'own read')
    expect('C16 PROTECTED claimed on a name a rebind MOVES is refused -- and '
           'this is the 2026-10-22 error verbatim', lie_p,
           FalseCaptureVerdict)

    lie_e = dict(adj)
    lie_e[COX] = (EXPOSED, 'a reason describing nothing')
    expect('C17 EXPOSED claimed on a name a rebind leaves bitwise unchanged '
           'is refused -- the half a census cannot have', lie_e,
           FalseCaptureVerdict)

    lie_g = dict(adj)
    lie_g[TR] = (GOVERNED, 'claims an existing mechanism moves it')
    expect('C18 GOVERNED claimed on a name in no rebind set is refused',
           lie_g, FalseCaptureVerdict)

    ok_g = dict(adj)
    ok_g[QC] = (GOVERNED, 'gfet.quantum_capacitance IS in the flip\'s rebind '
                          'set, so GOVERNED is a true statement about it')
    expect('C19 GOVERNED claimed on a name that IS in a rebind set is '
           'ACCEPTED -- the verdict is reachable, not decorative', ok_g, None)

    noreason = dict(adj)
    noreason[TR] = (EXPOSED, '')
    expect('C20 a verdict with NO stated reason is refused', noreason,
           UnadjudicatedCaptureSurface)

    badv = dict(adj)
    badv[TR] = ('PROBABLY_FINE', 'not in the vocabulary')
    expect('C21 a verdict outside {%s} is refused' % '/'.join(VERDICTS),
           badv, UnadjudicatedCaptureSurface)

    notpair = dict(adj)
    notpair[TR] = EXPOSED
    expect('C22 an adjudication that is not a (verdict, reason) pair is '
           'refused', notpair, UnadjudicatedCaptureSurface)

    hidden = dict(adj)
    del hidden[TR]
    hidden_held = {k: v for k, v in hidden.items()}
    exc = attempt(hidden_held)
    ok = isinstance(exc, UnadjudicatedCaptureSurface) \
        and 'held there reads' in str(exc)
    check('C23 THE RECORDED ERROR AS A REFUSAL: PROTECTED on the held name '
          'with the held object\'s live read unanswered is refused, and the '
          'message says why', ok,
          str(exc).splitlines()[0][:72] if exc else 'no exception')

    undeserved = dict(adj)
    undeserved[CDQC] = (PROTECTED, 'credits the capture with a zero it did '
                                   'not cause')
    expect('C24 PROTECTED claimed on a name NOT on the read surface is '
           'refused -- the conflation of 4.6.17, sign reversed', undeserved,
           FalseCaptureVerdict)

    wrong_u = dict(adj)
    wrong_u[QC] = (UNREACHED, 'claims the holder never looks the name up')
    expect('C25 UNREACHED claimed on a name that IS on the read surface is '
           'refused', wrong_u, FalseCaptureVerdict)

    expect('C26 THE RESULT: the repaired capture\'s adjudication is REFUSED '
           'for the UNREPAIRED replacement', adj, FalseCaptureVerdict, unrep)
    print()
    print('    C26 is what the item asked for.  The same set of stated')
    print('    reasons is accepted for the repaired contact replacement and')
    print('    refused for the unrepaired one -- which Section 0 showed the')
    print('    static reader cannot distinguish at all.')
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
