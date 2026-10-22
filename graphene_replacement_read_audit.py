"""
graphene_replacement_read_audit.py
==================================

THE QUESTION, AND WHY IT IS RANKED FIRST

2026-10-21 found that the two halves of item 16 do not compose: the contact's
form-B_T replacement reads `gfet.quantum_capacitance` by attribute lookup at
call time, and the channel half of the same flip rebinds exactly that name.
`n_contact` landed 4.86x (Ni) to 10.77x (Cu) past the effect item 16 exists to
install, in the right direction and in range, off a green transcript.

It was found by accident.  The two guards this repository owns both answer
WHICH NAMES MOVE -- `bc.rebind()` guarantees no binding of a rebound name is
stale, `bc.rebind_by_operating_point()` guarantees no class member is left at
the shipped object -- and a replacement that READS a moved name is neither
stale nor unassigned.  So the log entry of 2026-10-21 created this item and
ranked it first, with the work stated as one question per replacement: which
module-level names does this function read at call time, and does any flip
that installs it rebind one of them?

This module answers that question with THREE independent readers, and the
answer to the item as asked is in Section 3.  The instrument built to answer
it then found two things the item did not ask about: a second exposure
strictly larger than the one asked about (Sections 4-5), and a limit on the
2026-10-21 repair itself (Section 7).

THE THREE READERS

  A  STATIC.  A bytecode walk of the transitive closure of a replacement's
     live name lookups: every (module, attr) reached by LOAD_GLOBAL of a repo
     module followed by LOAD_ATTR, plus every repo function reached through
     such a name or through a global, followed into recursively.  Closure
     cells are resolved from the real function object, so a CAPTURED object is
     distinguished from a LIVE name read -- which is the whole of the
     2026-10-21 repair.  Branch-INSENSITIVE, hence an over-approximation: it
     sees a read on a branch that does not execute.  That is the direction a
     guard should err in, and Section 6 V12 turns it into a result.

  B  SENTINELS.  An identity counting wrapper is installed at every member of
     the flip's rebind set, the replacement is called, and the counters that
     fired are the names it read.  Object substitution, so it sees a read
     whether the reader spells it `mod.attr` or holds it as a module global
     from `from mod import attr` -- the spelling that caused all three of the
     2026-10 delivery faults.  This is `blast.count_every_binding()`'s
     technique, pointed at a different question.

  C  MODULE-TYPE INTERCEPTION.  Each repo module's `__class__` is swapped for
     a ModuleType subclass whose `__getattribute__` and `__setattr__` record
     the access and then delegate.  Exactly neutral by construction -- it
     returns the real object -- and it sees two things B cannot: reads of
     NON-CALLABLES (`gfet.C_ox` is a float; a counting wrapper cannot wrap it)
     and WRITES.  It is blind where B is strong: a read spelled as a module
     global inside the owning module is a dict lookup and never reaches the
     type, so B and C are blind in opposite directions and agreeing is
     evidence rather than restatement.

B and C are independent mechanisms and agree exactly on all five replacements
(V3).  A is their over-approximation (V8).

WHAT THE ITEM ASKED (Section 3)

Of the four replacements the item names, EXACTLY ONE reads a name that is
rebound while it is installed: `blast._cq_consistent_contact`, the known
positive.  So the 2026-10-21 diagnosis was complete as well as correct, and
this audit converts that from a hope into a measurement.

But the item's own WORDING -- "does any flip THAT INSTALLS IT rebind one of
them" -- answers NO for that very replacement, and the audit reports it so
(V4b).  `blast.contact_on_consistent_argument()` rebinds
`cd.quantum_capacitance`; the replacement reads `gfet.quantum_capacitance`,
which a DIFFERENT flip rebinds.  The defect is CROSS-FLIP, so the relation is
replacement x every flip that can be installed CONCURRENTLY, not the flip
that installs it.  The item as worded would have cleared the one replacement
whose defect created the item.

WHAT IT DID NOT ASK (Sections 4-5)

Reader C reports that all five replacements READ `gfet.C_ox` at call time and
that THREE of them WRITE it.  `gfet.C_ox` is in nobody's rebind set, so it is
invisible to both existing guards AND to the item's own question, which is
about names a flip rebinds.  The relation that matters is therefore not a set
of names but a BIPARTITE RELATION between what replacements write and what
replacements read, and the 2026-10-21 mechanism (ii) -- "`series_factor`
reaches `gfet.C_ox` at call time, and the contact solve has it substituted for
its duration", costed there at +1.0204 to +2.2965 percentage points -- is one
edge of it rather than a quirk of one function.

Section 5 measures that edge's reachability directly, and the measurement
sharpens 2026-10-21 rather than repeating it.  Under the UNREPAIRED flip the
channel replacement is reached 192 times during the contact census, and 181 of
those 192 arrivals happen with `gfet.C_ox` already substituted by
`cd.C_interface`.  Under the REPAIRED flip it is reached 0 times.  So the two
mechanisms 2026-10-21 decomposed and found to partially cancel are not two
defects: they are ONE defect -- the live read -- with two consequences, and
181/192 of its arrivals carry both at once.  That is why capturing the shipped
object repairs both in one move, and why repairing the conspicuous half alone
was measured there as a regression.

TWO LIMITS ON THE 2026-10-21 REPAIR

From reader A (V12): `make_contact_cq`'s `live=` flag means the live read is
still IN the function, on a branch a caller can take.  A is branch-insensitive
and flags `live=False` exactly as it flags `live=True`; B and C, which
execute, separate them.  The repair is a default, not a removal, and
`make_contact_cq` is public.

From the two captured controls (Section 7, V16b-d): CAPTURE IS ONE LEVEL DEEP.
A cell protects the NAME it holds and not the names the held OBJECT reads at
call time.  The repair is safe because the shipped `gfet.quantum_capacitance`
reads nothing in the rebind set -- a contingent fact about that one function,
not a property of capturing.  Its neighbour `gfet.carrier_density` reads
`gfet.quantum_capacitance` at call time, so the identical construction applied
to it is NOT immune, and the planted captured-leaky control demonstrates that
rather than arguing it.  Had the contact needed the shipped charge relation
instead of the shipped capacitance, the same repair would have been written,
would have looked right, and would have left the defect in place.
"""

import contextlib
import dis
import os
import sys
import types

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import graphene_fet_model as gfet
import graphene_contact_doping_model as cd
import graphene_quantum_capacitance_argument as qca
import graphene_contact_argument_blast_radius as blast
import graphene_default_flip as flip
import graphene_binding_class as bc

REPO = os.path.dirname(os.path.abspath(__file__))


def repo_modules():
    """Every imported module whose file lives in this repository.

    Same criterion as bc.repo_module_objects() and for the same reason: the
    class that matters is the one in the process doing the patching, and
    importing a module to look at it would change what is measured.
    """
    out = {}
    for name, mod in list(sys.modules.items()):
        if mod is None or name.startswith('__'):
            continue
        f = getattr(mod, '__file__', None)
        if not f:
            continue
        try:
            if os.path.dirname(os.path.abspath(f)) == REPO:
                out[name] = mod
        except (TypeError, ValueError):
            continue
    return out


# ---------------------------------------------------------------------------
# READER A -- static, branch-insensitive, over-approximating
# ---------------------------------------------------------------------------

_MISSING = object()
_LOAD_OPS = ('LOAD_GLOBAL', 'LOAD_DEREF', 'LOAD_CLASSDEREF', 'LOAD_NAME')
_ATTR_OPS = ('LOAD_ATTR', 'LOAD_METHOD')


def _closure_map(fn):
    """{freevar_name: cell contents} for a real function object.

    Resolving closures from the OBJECT rather than from the code is what lets
    reader A tell a captured object from a live name: `shipped_cq` in
    `flip.make_contact_cq` is a cell holding the shipped function, and a cell
    is immune to a rebind of the name the object also happens to sit at.
    """
    co = fn.__code__
    if not co.co_freevars or not fn.__closure__:
        return {}
    out = {}
    for name, cell in zip(co.co_freevars, fn.__closure__):
        try:
            out[name] = cell.cell_contents
        except ValueError:          # cell not yet filled
            pass
    return out


def static_reads(fn, maxdepth=8):
    """(live, captured) for `fn`, transitively.

    live:     {(module_name, attr)} reached by a LIVE lookup -- a name that a
              rebind of that name would change under this function's feet.
    captured: {(module_name, qualname)} for repo functions reached through a
              CLOSURE CELL -- immune to any rebind, and reported so that a
              repair by capture is visible rather than inferred.

    Recursion is bounded by `maxdepth` and by a visited set on code identity,
    so a self-recursive or mutually recursive replacement terminates (V12).
    """
    mods = repo_modules()
    mod_by_id = {id(m): m.__name__ for m in mods.values()}
    glob_ids = {id(m.__dict__) for m in mods.values()}
    live, captured, seen = set(), set(), set()

    def walk(code, g, freemap, depth):
        if depth > maxdepth or id(code) in seen:
            return
        seen.add(id(code))
        ins = list(dis.get_instructions(code))
        for i, x in enumerate(ins):
            if x.opname not in _LOAD_OPS:
                continue
            name = x.argval
            if x.opname in ('LOAD_DEREF', 'LOAD_CLASSDEREF'):
                base, origin = freemap.get(name, _MISSING), 'closure'
            else:
                base, origin = g.get(name, _MISSING), 'global'
            if base is _MISSING:
                continue
            nxt = ins[i + 1] if i + 1 < len(ins) else None
            if id(base) in mod_by_id and nxt is not None \
                    and nxt.opname in _ATTR_OPS:
                mname, attr = mod_by_id[id(base)], nxt.argval
                live.add((mname, attr))
                tgt = getattr(mods[mname], attr, _MISSING)
                if isinstance(tgt, types.FunctionType) \
                        and id(tgt.__globals__) in glob_ids:
                    walk(tgt.__code__, tgt.__globals__,
                         _closure_map(tgt), depth + 1)
            elif isinstance(base, types.FunctionType) \
                    and id(base.__globals__) in glob_ids:
                owner = base.__globals__.get('__name__')
                if origin == 'global':
                    live.add((owner, name))
                else:
                    captured.add((owner, getattr(base, '__name__', name)))
                walk(base.__code__, base.__globals__,
                     _closure_map(base), depth + 1)
        for const in code.co_consts:
            if isinstance(const, types.CodeType):
                walk(const, g, freemap, depth)

    walk(fn.__code__, fn.__globals__, _closure_map(fn), 0)
    return live, captured


# ---------------------------------------------------------------------------
# READER B -- identity sentinels on a named set
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def sentinels(sites):
    """Install an identity counting wrapper at each (module, attr) in `sites`.

    The wrapper changes no value, so anything measured inside must be
    bit-identical to the same thing measured outside; V1 checks that before
    any count is used, on the 2026-10-13 rule that an instrument which
    perturbs what it measures reports its own arrival.

    A non-callable member cannot be wrapped and is reported as unwrappable
    rather than skipped silently -- reader C is the one that sees those.
    """
    mods = repo_modules()
    hits, unwrappable, saved = {}, [], []

    def wrap(key, obj):
        def wrapper(*a, **kw):
            hits[key] = hits.get(key, 0) + 1
            return obj(*a, **kw)
        return wrapper

    try:
        for site in sorted(sites):
            mname, attr = site
            obj = getattr(mods[mname], attr)
            if not callable(obj):
                unwrappable.append(site)
                continue
            saved.append((mname, attr, obj))
            setattr(mods[mname], attr, wrap(site, obj))
        yield hits, unwrappable
    finally:
        for mname, attr, obj in reversed(saved):
            setattr(mods[mname], attr, obj)


# ---------------------------------------------------------------------------
# READER C -- module-type attribute interception (reads AND writes)
# ---------------------------------------------------------------------------

def _tracking_type(log):
    MT = types.ModuleType

    class Tracked(MT):
        def __getattribute__(self, name):
            if not name.startswith('__'):
                log['read'].add((MT.__getattribute__(self, '__name__'), name))
            return MT.__getattribute__(self, name)

        def __setattr__(self, name, value):
            if not name.startswith('__'):
                log['write'].add((MT.__getattribute__(self, '__name__'), name))
            MT.__setattr__(self, name, value)

    return Tracked


@contextlib.contextmanager
def tracking():
    """Record every cross-module attribute read and write on repo modules.

    Exactly neutral: `__getattribute__` delegates to ModuleType's and returns
    the real object, so unlike a wrapper there is nothing to be bit-unequal
    about.  V2 checks it anyway, because that argument has been wrong before.

    BLIND SPOT, stated rather than discovered later: a name read as a module
    GLOBAL inside its own module (`from gfet import quantum_capacitance`, then
    a bare call) is a dict lookup on the frame's globals and never reaches the
    module type.  Reader B sees exactly those.  The two readers are blind in
    opposite directions, which is why V3 asking them to agree is a test and
    not a tautology.
    """
    mods = list(repo_modules().values())
    log = {'read': set(), 'write': set()}
    Tracked = _tracking_type(log)
    saved = [(m, m.__class__) for m in mods]
    try:
        for m in mods:
            m.__class__ = Tracked
        yield log
    finally:
        for m, cls in saved:
            m.__class__ = cls


# ---------------------------------------------------------------------------
# What a flip rebinds, DERIVED FROM THE FLIP
# ---------------------------------------------------------------------------

def rebind_set(manager_factory):
    """The set of (module, attr) a manager rebinds, measured at its entry.

    Not a re-implementation of `identity_class()` and not a list: a rebind IS
    a `setattr` on a module, reader C records writes, so entering the manager
    under tracking and collecting the writes asks the manager itself what it
    moves.  A manager that grew a third name would be reported as having
    three without this module being edited -- which is the property the
    hand-maintained lists of 2026-10 did not have.
    """
    with tracking() as log:
        with manager_factory():
            moved = set(log['write'])
    return moved


FLIPS = (
    ('qca.on_consistent_argument()',
     lambda: qca.on_consistent_argument()),
    ('blast.contact_on_consistent_argument()',
     lambda: blast.contact_on_consistent_argument()),
    ('flip.on_default_flip(live_contact=False)',
     lambda: flip.on_default_flip(live_contact=False)),
)


# ---------------------------------------------------------------------------
# Planted controls
# ---------------------------------------------------------------------------

def _control_positive(V):
    """Reads a rebound name through a live attribute lookup.  Must be flagged
    by all three readers; if it is not, the readers are not working."""
    return gfet.carrier_density(np.asarray(V, dtype=float))


def _control_negative(V):
    """Reads nothing in any rebind set.  Must be flagged by none."""
    return 2.0 * np.asarray(V, dtype=float)


def _control_recursive(k, V=1.0):
    """Self-recursive AND reads a rebound name: reader A must terminate and
    still flag it (V12).  Recursion depth is deliberately tiny at runtime so
    B and C see the same read."""
    if k <= 0:
        return gfet.carrier_density(np.asarray(V, dtype=float))
    return _control_recursive(k - 1, V)


def make_captured_control(captured_fn):
    """The 2026-10-21 repair as a control: the same read, through a CELL.

    Reader A must report it as captured and NOT as a live read; B and C, which
    execute, must not flag it at all.  This is the negative half of the pair
    that makes V5 meaningful -- a reader that cannot tell a cell from a name
    would pass V5 for the wrong reason.
    """
    def control(V):
        return captured_fn(np.asarray(V, dtype=float))
    return control


# ---------------------------------------------------------------------------
# Reachability of the write/read edge (Section 6)
# ---------------------------------------------------------------------------

def channel_arrival_probe(live_contact):
    """Count channel-replacement arrivals during the contact census, and how
    many of them happen with `gfet.C_ox` already substituted.

    The probe wraps the two channel replacements at their DEFINITION site in
    `qca`, before the flip reads them, and records `gfet.C_ox` as seen from
    inside the call.  It is an identity wrapper, so the censuses it runs under
    are the committed ones; what it adds is the value of one global at each
    arrival.

    Returns (total, inside, outside, values_seen).
    """
    seen = []
    saved = (qca._quantum_capacitance_consistent,
             qca._carrier_density_consistent)

    def probe(fn):
        def wrapper(*a, **kw):
            seen.append(gfet.C_ox)
            return fn(*a, **kw)
        return wrapper

    qca._quantum_capacitance_consistent = probe(saved[0])
    qca._carrier_density_consistent = probe(saved[1])
    try:
        with flip.on_default_flip(live_contact=live_contact):
            census = {m: float(cd.contact_edge_carrier_density(wf))
                      for m, wf in flip.METALS}
            cd.recalibrate_metal_rc()
    finally:
        (qca._quantum_capacitance_consistent,
         qca._carrier_density_consistent) = saved
    inside = sum(1 for c in seen if c == cd.C_interface)
    return len(seen), inside, len(seen) - inside, census


# ---------------------------------------------------------------------------

def _hexes(d):
    return {k: float(v).hex() for k, v in sorted(d.items())}


def main():
    checks = []

    def check(name, ok, detail=''):
        checks.append((name, bool(ok), detail))
        print('    %-4s %s%s' % ('PASS' if ok else 'FAIL', name,
                                 ('   -- ' + detail) if detail else ''))

    print('=' * 78)
    print('DOES ANY REPLACEMENT READ A NAME ITS OWN FLIP REBINDS?')
    print('=' * 78)
    print()
    print('C_ox        = %.6e F/m^2' % gfet.C_ox)
    print('C_interface = %.6e F/m^2   (ratio %.4f)'
          % (cd.C_interface, cd.C_interface / gfet.C_ox))
    print()

    shipped_cq = gfet.quantum_capacitance
    METALS = flip.METALS
    base_contact = {m: float(cd.contact_edge_carrier_density(wf))
                    for m, wf in METALS}

    # -- Section 0 ----------------------------------------------------------
    print('-' * 78)
    print('0.  NEITHER INSTRUMENT IS THE PERTURBATION')
    print('-' * 78)
    print()
    all_sites = set()
    for _, fac in FLIPS:
        all_sites |= rebind_set(fac)
    with sentinels(all_sites) as (hits_b, unwrappable_b):
        under_b = {m: float(cd.contact_edge_carrier_density(wf))
                   for m, wf in METALS}
    with tracking():
        under_c = {m: float(cd.contact_edge_carrier_density(wf))
                   for m, wf in METALS}
    print('    7 committed contact densities, bitwise, under each reader:')
    for m, _ in METALS:
        print('      %-3s  %s  B %s  C %s' % (
            m, float(base_contact[m]).hex(),
            'same' if base_contact[m].hex() == under_b[m].hex() else 'DIFFERS',
            'same' if base_contact[m].hex() == under_c[m].hex() else 'DIFFERS'))
    print()
    check('V1  reader B (sentinels) is bitwise neutral on all 7',
          _hexes(under_b) == _hexes(base_contact))
    check('V2  reader C (module type) is bitwise neutral on all 7',
          _hexes(under_c) == _hexes(base_contact))
    print()

    # -- Section 1 ----------------------------------------------------------
    print('-' * 78)
    print('1.  WHAT EACH FLIP REBINDS, ASKED OF THE FLIP')
    print('-' * 78)
    print()
    rebinds = {}
    for label, fac in FLIPS:
        moved = rebind_set(fac)
        rebinds[label] = moved
        print('    %s' % label)
        for s in sorted(moved):
            print('        %s.%s' % s)
    composed = set().union(*rebinds.values())
    rebinds['THE COMPOSITION (all three at once)'] = composed
    print()
    print('    THE COMPOSITION (all three at once)')
    for s in sorted(composed):
        print('        %s.%s' % s)
    print()
    print('    Derived from each manager\'s own setattr calls, not from a')
    print('    re-run of identity_class() -- so a manager that grows a name')
    print('    is reported with it and this module is not edited.')
    ic = set(bc.identity_class(gfet, 'carrier_density')) | \
        set(bc.identity_class(gfet, 'quantum_capacitance'))
    check('V3a composed rebind set == the two identity classes exactly',
          composed == ic,
          '%d sites' % len(composed))
    print()

    # -- Section 2-4: the readers on every replacement ----------------------
    print('-' * 78)
    print('2.  THE REPLACEMENTS, UNDER ALL THREE READERS')
    print('-' * 78)
    print()
    captured_clean = make_captured_control(gfet.quantum_capacitance)
    captured_leaky = make_captured_control(gfet.carrier_density)
    REPS = (
        ('qca._carrier_density_consistent',
         qca._carrier_density_consistent, (0.5,), True),
        ('qca._quantum_capacitance_consistent',
         qca._quantum_capacitance_consistent, (0.5,), True),
        ('blast._cq_consistent_contact',
         blast._cq_consistent_contact, (0.5,), True),
        ('flip.make_contact_cq(live=False)  [the 10-21 repair]',
         flip.make_contact_cq(shipped_cq, live=False), (0.5,), True),
        ('flip.make_contact_cq(live=True)   [the 10-21 defect]',
         flip.make_contact_cq(shipped_cq, live=True), (0.5,), True),
        ('CONTROL positive  (live read of a rebound name)',
         _control_positive, (0.5,), False),
        ('CONTROL negative  (reads nothing rebound)',
         _control_negative, (0.5,), False),
        ('CONTROL recursive (self-recursive + live read)',
         _control_recursive, (3, 0.5), False),
        ('CONTROL captured-clean (cell holding gfet.quantum_capacitance)',
         captured_clean, (0.5,), False),
        ('CONTROL captured-leaky (cell holding gfet.carrier_density)',
         captured_leaky, (0.5,), False),
    )
    results = {}
    for label, fn, args, is_item in REPS:
        a_live, a_cap = static_reads(fn)
        with sentinels(composed) as (hits, unwrappable):
            fn(*args)
        b_set = set(hits)
        with tracking() as log:
            fn(*args)
        c_read, c_write = set(log['read']), set(log['write'])
        results[label] = dict(a_live=a_live, a_cap=a_cap, b=b_set,
                              c_read=c_read, c_write=c_write,
                              is_item=is_item,
                              unwrappable=set(unwrappable))
        print('    %s' % label)
        print('        A static  live reads n=%-3d  in rebind set: %s'
              % (len(a_live), sorted(a_live & composed) or '--'))
        print('                  captured (immune): %s'
              % (sorted(a_cap) or '--'))
        print('        B sentinel fired          : %s'
              % (sorted(b_set) or '--'))
        print('        C tracked reads in set    : %s'
              % (sorted(c_read & composed) or '--'))
        print('        C tracked WRITES          : %s'
              % (sorted(c_write) or '--'))
        print()
    print('    Reader B could not wrap %d member(s) of the rebind set: %s'
          % (len(results[REPS[0][0]]['unwrappable']),
             sorted(results[REPS[0][0]]['unwrappable']) or '(all callable)'))
    print()
    # -- Section 3 ----------------------------------------------------------
    print('-' * 78)
    print('3.  THE CROSS TABLE -- AND THE ITEM AS WORDED IS TOO NARROW')
    print('-' * 78)
    print()
    INSTALLED_BY = {
        'qca._carrier_density_consistent':
            ('qca.on_consistent_argument()',
             'flip.on_default_flip(live_contact=False)'),
        'qca._quantum_capacitance_consistent':
            ('qca.on_consistent_argument()',
             'flip.on_default_flip(live_contact=False)'),
        'blast._cq_consistent_contact':
            ('blast.contact_on_consistent_argument()',),
        'flip.make_contact_cq(live=False)  [the 10-21 repair]':
            ('flip.on_default_flip(live_contact=False)',),
        'flip.make_contact_cq(live=True)   [the 10-21 defect]':
            ('flip.on_default_flip(live_contact=False)',),
    }
    flip_labels = [lab for lab, _ in FLIPS]
    hdr = '    %-38s' % 'replacement \\ flip'
    for lab in flip_labels:
        hdr += ' %-9s' % lab.split('.')[0][:8]
    hdr += ' %-9s' % 'COMPOSED'
    print(hdr)
    print('    %s' % ('-' * 74))
    narrow_hits, broad_hits = {}, {}
    for label, fn, args, is_item in REPS:
        if not is_item:
            continue
        r = results[label]
        row = '    %-38s' % label[:38]
        own = set(INSTALLED_BY[label])
        n_narrow = n_broad = 0
        for lab in flip_labels:
            conflict = r['b'] & rebinds[lab]
            row += ' %-9s' % (('YES(%d)' % len(conflict)) if conflict else '.')
            if conflict:
                n_broad += 1
                if lab in own:
                    n_narrow += 1
        comp = r['b'] & composed
        row += ' %-9s' % (('YES(%d)' % len(comp)) if comp else '.')
        narrow_hits[label], broad_hits[label] = n_narrow, n_broad
        print(row)
    print()
    print('    THE ITEM ASKED: "does any flip THAT INSTALLS IT rebind one of')
    print('    the names it reads".  Answered literally, for the one')
    print('    replacement that carries the 2026-10-21 defect:')
    key = 'blast._cq_consistent_contact'
    own_sets = [rebinds[l] for l in INSTALLED_BY[key]]
    own_conflict = results[key]['b'] & set().union(*own_sets)
    print('        %s' % key)
    print('        installed by : %s' % ', '.join(INSTALLED_BY[key]))
    print('        that flip rebinds : %s'
          % sorted(set().union(*own_sets)))
    print('        it reads          : %s' % sorted(results[key]['b']))
    print('        intersection      : %s' % (sorted(own_conflict) or 'EMPTY'))
    print()
    print('    EMPTY.  The item\'s own wording would have cleared the one')
    print('    replacement whose defect created the item.  The defect is')
    print('    CROSS-FLIP: `blast.contact_on_consistent_argument()` rebinds')
    print('    `cd.quantum_capacitance` and the replacement reads')
    print('    `gfet.quantum_capacitance`, which a DIFFERENT flip rebinds.')
    print('    So the relation is replacement x every flip that can be')
    print('    installed CONCURRENTLY, not the flip that installs it, and')
    print('    the composed column is the one that answers the question.')
    print()
    check('V4  the known positive is detected under the composition',
          results[key]['b'] & composed ==
          {('graphene_fet_model', 'quantum_capacitance')},
          'exactly {gfet.quantum_capacitance}')
    check('V4b and is NOT detected by the item\'s own narrow wording',
          own_conflict == set(),
          'the wording is the finding')
    item_four = ['qca._carrier_density_consistent',
                 'qca._quantum_capacitance_consistent',
                 'blast._cq_consistent_contact',
                 'flip.make_contact_cq(live=False)  [the 10-21 repair]']
    nonempty = [l for l in item_four if results[l]['b'] & composed]
    check('V5  the repair: make_contact_cq(live=False) conflict is EMPTY',
          results['flip.make_contact_cq(live=False)  [the 10-21 repair]']['b']
          & composed == set())
    check('V6  positive control on V5: live=True conflict == V4\'s set',
          results['flip.make_contact_cq(live=True)   [the 10-21 defect]']['b']
          & composed == {('graphene_fet_model', 'quantum_capacitance')},
          'if both were empty the instrument would be broken')
    check('V7  of the four replacements the item names, exactly 1 conflicts',
          len(nonempty) == 1 and nonempty[0] == key,
          'the accident of 10-21 found the whole exposure, as framed')
    print()

    # -- Section 4 ----------------------------------------------------------
    print('-' * 78)
    print('4.  WHAT THE ITEM DID NOT ASK: THE WRITE/READ RELATION')
    print('-' * 78)
    print()
    writers, readers = {}, {}
    for label, fn, args, is_item in REPS:
        if not is_item:
            continue
        for site in results[label]['c_write']:
            writers.setdefault(site, []).append(label)
        for site in results[label]['c_read']:
            readers.setdefault(site, []).append(label)
    edges = {s: (writers[s], readers.get(s, []))
             for s in sorted(writers) if readers.get(s)}
    print('    A name WRITTEN by one replacement and READ by another is')
    print('    outside both existing guards (it is in nobody\'s identity')
    print('    class) and outside the item\'s question (which is about names')
    print('    a flip rebinds).  Reader C sees it because it intercepts')
    print('    writes and because it can see a float.')
    print()
    for site, (ws, rs) in edges.items():
        print('    %s.%s' % site)
        print('        in any rebind set? %s'
              % ('YES' if site in composed else 'NO -- invisible to both guards'))
        print('        WRITTEN by %d/%d replacements:' % (len(ws), 5))
        for w in ws:
            print('            %s' % w)
        print('        READ by %d/%d replacements:' % (len(rs), 5))
        for r in rs:
            print('            %s' % r)
    print()
    check('V8a the write/read relation is non-empty',
          len(edges) >= 1, '%d edge name(s)' % len(edges))
    check('V8b and none of its names is in any rebind set',
          all(s not in composed for s in edges),
          'so neither existing guard can see it')
    print()

    # -- Section 5 ----------------------------------------------------------
    print('-' * 78)
    print("5.  REACHABILITY OF THE EDGE, AND 10-21's TWO MECHANISMS ARE ONE")
    print('-' * 78)
    print()
    tot_r, in_r, out_r, cen_r = channel_arrival_probe(live_contact=False)
    tot_d, in_d, out_d, cen_d = channel_arrival_probe(live_contact=True)
    print('    channel-replacement arrivals during the contact census:')
    print('        REPAIRED   flip (live_contact=False): %3d  '
          '(inside the C_interface window: %d)' % (tot_r, in_r))
    print('        UNREPAIRED flip (live_contact=True) : %3d  '
          '(inside the C_interface window: %d, outside: %d)'
          % (tot_d, in_d, out_d))
    print()
    print('    2026-10-21 decomposed the unrepaired composition into (i) the')
    print('    channel solve being reached at all and (ii) that solve seeing')
    print('    C_ox substituted, and measured them as partially cancelling')
    print('    (+1.0204 to +2.2965 pp).  These counts say they are not two')
    print('    defects: %d of the %d arrivals carry BOTH at once, so (ii)'
          % (in_d, tot_d))
    print('    cannot occur without (i).  One cause, two consequences --')
    print('    which is why capturing the shipped object repairs both in one')
    print('    move, and why repairing the conspicuous half alone was')
    print('    measured there as a regression.')
    print()
    check('V9  the repaired flip reaches the channel replacement 0 times',
          tot_r == 0, 'exactly 0')
    check('V9b the unrepaired flip reaches it, and mostly inside the window',
          tot_d > 0 and in_d > 0 and in_d <= tot_d,
          '%d arrivals, %d inside' % (tot_d, in_d))
    check('V9c inside + outside == total (the probe partitions its own count)',
          in_d + out_d == tot_d)
    check('V9d recorded expectation: 0 / 192 / 181',
          (tot_r, tot_d, in_d) == (0, 192, 181),
          'got %d / %d / %d -- if this moves, the composition moved'
          % (tot_r, tot_d, in_d))
    with blast.contact_on_consistent_argument():
        contact_only = {m: float(cd.contact_edge_carrier_density(wf))
                        for m, wf in METALS}
    check('V9e the repaired flip reproduces the committed contact-only '
          'form B_T BITWISE on all 7 metals',
          _hexes(cen_r) == _hexes(contact_only),
          'the 10-21 claim, re-measured by a different instrument')
    paired = [(m, cen_d[m], contact_only[m]) for m, _ in METALS
              if contact_only[m] != 0.0]
    zeros = [m for m, _ in METALS if contact_only[m] == 0.0]
    worst = max(abs(d / c - 1.0) for _, d, c in paired)
    best = min(abs(d / c - 1.0) for _, d, c in paired)
    check('V9f and the unrepaired flip does not (positive control on V9e)',
          _hexes(cen_d) != _hexes(contact_only),
          'deviation %.4f %% to %.4f %% over the %d metals with a non-zero '
          'denominator' % (100.0 * best, 100.0 * worst, len(paired)))
    check('V9g the metal(s) with n_contact exactly 0.0 are 0.0 under BOTH '
          'the repaired and the unrepaired flip',
          bool(zeros) and all(cen_d[m] == 0.0 and cen_r[m] == 0.0
                              for m in zeros),
          'exact zero at dW_F = 0: %s -- a defect in the composition cannot '
          'move it, so it is not evidence either way and is excluded from '
          'the ratio above rather than quietly divided by'
          % ', '.join(zeros))
    print()
    # -- Section 6 ----------------------------------------------------------
    print('-' * 78)
    print('6.  THE READERS AGAINST EACH OTHER, AND THE PLANTED CONTROLS')
    print('-' * 78)
    print()
    real_five = [l for l, _, _, i in REPS if i]
    all_subj = [l for l, _, _, _ in REPS]
    disagree = [l for l in all_subj
                if results[l]['b'] != (results[l]['c_read'] & composed)]
    print('    B (object substitution) vs C (module type), on the rebind set:')
    for label, fn, args, is_item in REPS:
        r = results[label]
        print('        %-52s B %-34s C %s'
              % (label[:52], sorted(r['b']) or '--',
                 sorted(r['c_read'] & composed) or '--'))
    print()
    check('V10 B and C agree exactly on all %d real replacements'
          % len(real_five),
          not [l for l in real_five if l in disagree],
          'two mechanisms blind in opposite directions, same answer')
    check('V10b and DISAGREE on %d planted control(s) -- C\'s stated blind '
          'spot is real and the controls reach it' % len(disagree),
          len(disagree) >= 1,
          'on %s; each is a read of a rebound name made from INSIDE its own '
          'module, a globals-dict lookup that never reaches the module type.'
          '  B sees them, C cannot.  The five real replacements happen not '
          'to exercise it, so without the controls this would have been '
          'invisible and C would have looked interchangeable with B.'
          % ', '.join(l.split('(')[0].strip() for l in disagree))
    contain = [l for l in results
               if not (results[l]['b'] <= results[l]['a_live'])]
    check('V11 A over-approximates B for every subject (A superset of B)',
          not contain)
    sl_f = static_reads(flip.make_contact_cq(shipped_cq, live=False))[0]
    sl_t = static_reads(flip.make_contact_cq(shipped_cq, live=True))[0]
    key_f = 'flip.make_contact_cq(live=False)  [the 10-21 repair]'
    key_t = 'flip.make_contact_cq(live=True)   [the 10-21 defect]'
    check('V12 A is branch-insensitive: static(live=False) == '
          'static(live=True)', sl_f == sl_t)
    check('V12b so A and B disagree on exactly the repaired replacement',
          (sl_f & composed) != results[key_f]['b']
          and (sl_t & composed) == results[key_t]['b'],
          'the live read is still IN the function, on a branch a caller '
          'can take')
    print()
    print('    That disagreement is the honest limit of the 10-21 repair:')
    print('    `make_contact_cq` is public and `live=True` re-arms the')
    print('    defect.  A default is not a removal.')
    print()
    cp = 'CONTROL positive  (live read of a rebound name)'
    cn = 'CONTROL negative  (reads nothing rebound)'
    cr = 'CONTROL recursive (self-recursive + live read)'
    cc = 'CONTROL captured-clean (cell holding gfet.quantum_capacitance)'
    cl = 'CONTROL captured-leaky (cell holding gfet.carrier_density)'
    check('V13 planted positive is flagged by all three readers',
          results[cp]['b'] and (results[cp]['c_read'] & composed)
          and (results[cp]['a_live'] & composed))
    check('V14 planted negative is flagged by none',
          not results[cn]['b'] and not (results[cn]['c_read'] & composed)
          and not (results[cn]['a_live'] & composed))
    check('V15 reader A terminates on a self-recursive subject and flags it',
          bool(results[cr]['a_live'] & composed)
          and bool(results[cr]['b']))
    check('V16 captured-CLEAN control: A reports the cell as captured, '
          'and no reader flags a live read',
          ('graphene_fet_model', 'quantum_capacitance')
          in results[cc]['a_cap']
          and not (results[cc]['a_live'] & composed)
          and not results[cc]['b']
          and not (results[cc]['c_read'] & composed),
          'a reader that cannot tell a cell from a name would pass V5 for '
          'the wrong reason')
    check('V16b captured-LEAKY control: the cell is captured AND B still '
          'flags a conflict',
          ('graphene_fet_model', 'carrier_density') in results[cl]['a_cap']
          and results[cl]['b'] ==
          {('graphene_fet_model', 'quantum_capacitance')},
          'CAPTURE IS ONE LEVEL DEEP -- see Section 7')
    shipped_qc_live = static_reads(shipped_cq)[0]
    shipped_cd_live = static_reads(gfet.carrier_density)[0]
    check('V16c the 10-21 repair is safe only CONTINGENTLY: the object it '
          'captures reads nothing in the rebind set',
          not (shipped_qc_live & composed),
          'gfet.quantum_capacitance\'s own live reads: %s'
          % (sorted(shipped_qc_live) or '--'))
    check('V16d and the neighbouring shipped object would NOT have been '
          'safe to capture (positive control on V16c)',
          bool(shipped_cd_live & composed),
          'gfet.carrier_density reads %s at call time'
          % sorted(shipped_cd_live & composed))
    with sentinels(composed):
        with tracking() as log2:
            blast._cq_consistent_contact(0.5)
    check('V17 C reports the same read set with B\'s sentinels installed',
          (set(log2['read']) & composed) == results[key]['b'],
          'neither reader is creating the read it reports')
    print()

    # -- Section 7 ----------------------------------------------------------
    print('-' * 78)
    print('7.  CAPTURE IS ONE LEVEL DEEP')
    print('-' * 78)
    print()
    print('    The 2026-10-21 repair binds the contact replacement to the')
    print('    shipped `gfet.quantum_capacitance` CAPTURED AT MANAGER ENTRY,')
    print('    and V5 confirms the conflict is gone.  The two captured')
    print('    controls say why, and the reason is weaker than it looks.')
    print()
    print('    A cell protects the NAME it holds.  It does not protect the')
    print('    names the held OBJECT reads at call time:')
    print()
    print('        gfet.quantum_capacitance  live reads: %s'
          % (sorted(shipped_qc_live) or '--'))
    print('            in the rebind set: %s'
          % (sorted(shipped_qc_live & composed) or 'NONE -- safe to capture'))
    print('        gfet.carrier_density      live reads: %s'
          % sorted(shipped_cd_live))
    print('            in the rebind set: %s  <-- NOT safe to capture'
          % sorted(shipped_cd_live & composed))
    print()
    print('    Both captures look identical at the call site:')
    print('    `make_contact_cq(shipped, live=False)` and')
    print('    `make_captured_control(gfet.carrier_density)` are the same')
    print('    construction.  One is immune and one is not, and which it is')
    print('    depends entirely on what the captured object reads.  So the')
    print('    repair of 2026-10-21 is CORRECT and its correctness is')
    print('    CONTINGENT: had the contact needed the shipped charge')
    print('    relation rather than the shipped capacitance, the same repair')
    print('    would have been written, would have looked right, and would')
    print('    have left the defect in place.')
    print()
    print('    This is the third member of the family, stated in the form')
    print('    the other two take:')
    print('        2026-10-14  a repair to one binding is not a repair to')
    print('                    the class.')
    print('        2026-10-17  a class can span two operating points, so it')
    print('                    must be partitioned.')
    print('        2026-10-21  a replacement that READS a rebound name is')
    print('                    neither stale nor unassigned.')
    print('        TODAY       a capture that fixes one name does not fix')
    print('                    the names the captured object reads.  The')
    print('                    safe criterion is the TRANSITIVE live-read')
    print('                    set of the captured object, which is what')
    print('                    reader A computes.')
    print()

    # -- figure -------------------------------------------------------------
    item_labels = [l for l, _, _, i in REPS if i]
    short = ['qca._carrier\n_density', 'qca._quantum\n_capacitance',
             'blast._cq_\nconsistent_contact', 'flip.make_contact_cq\n'
             '(live=False)', 'flip.make_contact_cq\n(live=True)']
    n_live = [len(results[l]['a_live']) for l in item_labels]
    n_conf = [len(results[l]['b']) for l in item_labels]
    n_wr = [len(results[l]['c_write']) for l in item_labels]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 6.2))
    x = np.arange(len(item_labels))
    ax1.bar(x - 0.26, n_live, 0.26, label='live name reads (reader A, total)',
            color='#4c72b0')
    ax1.bar(x, n_conf, 0.26,
            label='of those, in the composed rebind set (B)', color='#c44e52')
    ax1.bar(x + 0.26, n_wr, 0.26, label='module globals WRITTEN (reader C)',
            color='#dd8452')
    ax1.set_xticks(x)
    ax1.set_xticklabels(short, fontsize=7.5)
    ax1.set_ylabel('count of (module, attr) names')
    n_conflicting = sum(1 for v in n_conf if v)
    n_writing = sum(1 for v in n_wr if v)
    # Headroom BEFORE the legend is placed, so the legend cannot sit on the
    # tallest bar -- which is what it did on the first render of this figure.
    ax1.set_ylim(0, max(n_live) * 1.42)
    ax1.set_title('(a) %d of %d replacements read a rebound name; %d of %d '
                  'write a module global'
                  % (n_conflicting, len(item_labels), n_writing,
                     len(item_labels)), fontsize=9.5)
    ax1.legend(fontsize=7.5, loc='upper left', framealpha=0.95)
    ax1.grid(axis='y', alpha=0.3)
    for xi, c in enumerate(n_conf):
        if c:
            ax1.annotate('CONFLICT', (xi, c), textcoords='offset points',
                         xytext=(0, 5), ha='center', fontsize=7.5,
                         color='#c44e52', fontweight='bold')

    # Panel (b).  Both legend entries are built from the SAME colours the
    # bars are drawn with, by reusing the artists, rather than from a
    # hand-written swatch: the first render legended the overlay green while
    # the overlay it labelled was drawn dark red.
    cats = ['repaired\n(live_contact=False)', 'unrepaired\n(live_contact=True)']
    tot, ins = [tot_r, tot_d], [in_r, in_d]
    h_tot = ax2.bar(cats, tot, color='#c44e52', width=0.55,
                    label='channel-replacement arrivals')
    h_in = ax2.bar(cats, ins, color='#6d1f23', width=0.30,
                   label='of those, with gfet.C_ox already substituted '
                         'by cd.C_interface')
    ax2.set_ylim(0, max(tot) * 1.30 if max(tot) else 1.0)
    for xi, (t, i) in enumerate(zip(tot, ins)):
        ax2.annotate('%d arrival%s' % (t, '' if t == 1 else 's'),
                     (xi, t), textcoords='offset points', xytext=(0, 5),
                     ha='center', fontsize=9,
                     fontweight='bold' if t == 0 else 'normal')
        if i:
            ax2.annotate('%d inside' % i, (xi, i),
                         textcoords='offset points', xytext=(0, -15),
                         ha='center', fontsize=8, color='white')
    # An exact zero drawn as a zero-height bar is indistinguishable from no
    # data, so it is marked rather than left to the axis.
    if tot[0] == 0:
        ax2.axhline(0, color='#55a868', lw=3, xmin=0.08, xmax=0.42)
        ax2.annotate('exactly 0 -- the repair removes BOTH mechanisms',
                     (0, 0), textcoords='offset points', xytext=(0, 26),
                     ha='center', fontsize=8.5, color='#2e7d4f')
    ax2.set_ylabel('channel-replacement arrivals during the contact census')
    ax2.set_title('(b) %d of %d arrivals carry both 10-21 mechanisms at '
                  'once; the repair gives %d'
                  % (in_d, tot_d, tot_r), fontsize=9.5)
    ax2.legend(handles=[h_tot, h_in], fontsize=8, loc='upper left',
               framealpha=0.95)
    ax2.grid(axis='y', alpha=0.3)

    n_item_narrow = sum(
        1 for l in item_four
        if results[l]['b'] & set().union(*[rebinds[x]
                                           for x in INSTALLED_BY[l]]))
    n_item_composed = sum(1 for l in item_four if results[l]['b'] & composed)
    fig.suptitle('Does any replacement read a name a flip rebinds?   '
                 '%d of the %d replacements the item names, under the '
                 'composition -- but %d of %d under the item\'s own wording '
                 '("the flip that installs it")'
                 % (n_item_composed, len(item_four), n_item_narrow,
                    len(item_four)), fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig('replacement_read_audit.png', dpi=150)
    plt.close(fig)
    print('    figure -> replacement_read_audit.png')
    print()

    print('=' * 78)
    npass = sum(1 for _, ok, _ in checks if ok)
    print('VALIDATIONS: %d/%d passed' % (npass, len(checks)))
    for name, ok, _ in checks:
        if not ok:
            print('    FAILED: %s' % name)
    print('=' * 78)
    return 0 if npass == len(checks) else 1


if __name__ == '__main__':
    sys.exit(main())
