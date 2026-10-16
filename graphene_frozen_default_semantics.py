#!/usr/bin/env python3
"""
graphene_frozen_default_semantics.py -- 2026-10-16

THREE STANDING CENSUSES COUNT THE SAME PATTERN AND DISAGREE, AND THE REASON
IS THAT THE PATTERN IS TWO CLASSES WITH OPPOSITE BEHAVIOUR.

Background.  Since 2026-09-28 this repository has tracked a fault class it
calls the INERT KNOB: a module-level constant captured as a function
parameter default,

    K = 1.0
    def f(x=K): ...

whose default is bound when the `def` STATEMENT executes.  A later rebinding
of `K` -- which is how every mutation audit and every patch manager in this
repository delivers a change -- therefore never reaches `f`.  A
mutation-based audit that rebinds `K` and watches `f` measures nothing and
reports a pass.  That is why three separate instruments census the pattern:

    graphene_figure_provenance_audit.census()                   Section C
    graphene_mutation_arrival_probe.repo_wide_frozen_census()   Section 7b
    graphene_cross_module_delivery_audit                        cross-module

ALL THREE AGREE IN THEIR COMMITTED TRANSCRIPTS -- 0, 0 (an empty census
behind a [PASS]) and 8 -- AND ALL THREE TRANSCRIPTS WERE LAST WRITTEN
2026-10-04.  Run today they report 21, 20 and 16.  Twelve days of sessions
regrew a fault class that 2026-10-04 had driven to zero, in front of three
censuses, none of which was being read -- and the fourth instrument whose
entire purpose is to notice exactly that (graphene_pristine_transcript_audit,
"A COMMITTED TRANSCRIPT IS A CLAIM ABOUT THE CODE BESIDE IT") was last run on
the same day.  It reports today, unprompted, that four of its subjects have
gone STALE or now ERROR.

THE DISCREPANCY IS THE INTERESTING PART.  21 against 20 is not a rounding
difference and not two definitions of "the repository": run over one path
set, both censuses name the same five files and four of the five agree
exactly.  They differ by ONE occurrence, in
graphene_velocity_saturation_model.py: line 589,
`_R_harmonic(V_g, V_hi, n_quad=N_QUAD)`, defined inside the body of another
function.

    A `def` IN THE MODULE BODY, OR IN A CLASS BODY, EXECUTES ONCE AT IMPORT,
    SO ITS CAPTURED DEFAULT IS FROZEN.  A `def` IN A FUNCTION BODY EXECUTES
    EVERY TIME THAT FUNCTION IS CALLED, SO ITS CAPTURED DEFAULT IS
    RE-EVALUATED EVERY TIME AND A REBINDING ARRIVES NORMALLY.  THE SECOND
    KIND IS NOT AN INERT KNOB.  IT IS LIVE.

Section 1 establishes that by construction and measurement rather than by
assertion, in all three configurations, against exactly known values.

So the two censuses are not both wrong, and neither is simply right:

  * figure_provenance's 21 is a correct count of the SYNTAX and an over-count
    of the FAULT by exactly the function-nested captures.  Its own historical
    figures -- "49 occurrences across 14 files, 14 of them inside the audit
    modules ... 37 remain, and the 13 in audit modules are the next run's top
    item" -- are syntactic counts of a class that is not uniformly the fault,
    and the 49/37/13 split has never been separated into the two.
  * mutation_arrival_probe's 20 is a correct count of the FAULT, reached for
    the wrong reason.  Its `classify()` confirms each statically-found site at
    run time with `getattr(module, fname)` -- added, correctly, to drop a
    function that a later redefinition or a decorator has replaced -- and a
    function-nested function is not a module attribute, so `getattr` returns
    None and the site is dropped SILENTLY.  No row, no class, no reason.  The
    right verdict by an accident of resolution, which is not the same thing
    as the right verdict: the mechanism cannot distinguish "nested, therefore
    live" from "no such site".  A capture inside a CLASS body would be
    dropped by the same line and would be frozen, so the accident does not
    generalise -- it is wrong in the undercounting direction for the next
    case rather than this one.

WHY UNDERCOUNTING IS THE DANGEROUS DIRECTION HERE.  The signature of this
fault class is that a mutation silently fails to arrive and the audit reports
a pass.  An instrument that silently fails to report a site has the identical
signature one level up.  Dropping a true site is the same shape as the fault
being censused.

A DEFECT IN THE FIRST DRAFT OF THIS MODULE, KEPT ON THE RECORD.  The first
draft classified by NESTING DEPTH: depth 0 frozen, depth >= 1 live.  Its own
control planted a method in a class body, called it depth 1, and asserted in
prose that it was therefore "NOT counted as an inert knob".  That is wrong: a
class body executes once at import, so a method's captured default is as
frozen as a module-level function's.  Measured in Section 1 (V3), which the
first draft did not measure.  The repository today contains no class-body
capture, so the first draft's count would have been right by luck -- the same
luck `classify()` is running on, which is what made the draft worth
recording rather than quietly fixing.  THE CRITERION IS NOT DEPTH.  IT IS
WHETHER ANY ENCLOSING SCOPE IS A FUNCTION BODY.

WHAT THIS MODULE DOES NOT DO.  It does not convert the 20 genuine inert
knobs, and it does not patch `classify()`.  It establishes the semantics,
reconciles the censuses to zero residue, and leaves the repository with a
count of the fault that is correct for a stated reason for the first time.
Converting 21 captures as one class would have rewritten a live capture into
a None-sentinel for no reason and recorded it as a fix.

Run:  python3 graphene_frozen_default_semantics.py
"""

import ast
import glob
import importlib
import os
import sys
import tempfile
import types

passed = 0
failed = 0

# Files excluded from mutation_arrival_probe's census, reproduced here so the
# two censuses can be compared over the SAME path set.  Read by name at the
# one place it is used and never as a parameter default -- see Section 5.
PROBE_EXCLUDES = ('initialize_repo', 'main', 'generate_plots',
                  'graphene_mutation_arrival_probe')


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
    return ok


# ---------------------------------------------------------------------------
# SECTION 1 -- the semantics, measured in every configuration
# ---------------------------------------------------------------------------

_PROBE_SRC = '''\
K = 1.0


def at_module_level(x=K):
    """`def` runs once, at import: the default is bound to 1.0 forever."""
    return x


class InClassBody:
    def as_a_method(self, x=K):
        """A class body also runs once, at import.  Frozen, like the above --
        and reachable by getattr only through the class, not the module."""
        return x


def enclosing():
    def in_a_function_body(x=K):
        """This `def` runs on every call of enclosing(), so the default is
        re-evaluated from whatever K names at that moment."""
        return x
    return in_a_function_body()


def literal_default(x=7.0):
    """A literal default: not a capture at all, and no census should report
    it."""
    return x
'''


def _probe_module():
    """Import a throwaway module carrying one capture of each kind.

    Built on disk and imported rather than exec()'d, because the thing being
    measured is WHEN a `def` statement runs, and a module imported from a
    file is the configuration every module in this repository is in.
    """
    d = tempfile.mkdtemp(prefix='frozen_default_semantics_')
    path = os.path.join(d, 'fds_probe.py')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(_PROBE_SRC)
    sys.path.insert(0, d)
    sys.modules.pop('fds_probe', None)
    return importlib.import_module('fds_probe')


def section1_semantics():
    print('-' * 74)
    print('SECTION 1 -- the classes, established against exactly known values')
    print('-' * 74)
    m = _probe_module()
    obj = m.InClassBody()

    before = (m.at_module_level(), obj.as_a_method(), m.enclosing(),
              m.literal_default())
    print('  before rebinding, K = %r' % (m.K,))
    print('    at_module_level()            = %r' % (before[0],))
    print('    InClassBody().as_a_method()  = %r   (def in a CLASS body)' % (before[1],))
    print('    enclosing()                  = %r   (def in a FUNCTION body)' % (before[2],))
    print('    literal_default()            = %r' % (before[3],))
    check('all four arms agree before the rebinding',
          before[0] == before[1] == before[2] == 1.0 and before[3] == 7.0,
          'a probe whose arms already disagree cannot attribute a later '
          'difference to the rebinding')

    m.K = 99.0
    after = (m.at_module_level(), obj.as_a_method(), m.enclosing(),
             m.literal_default())
    print('  after  rebinding, K = %r' % (m.K,))
    print('    at_module_level()            = %r   <-- did NOT arrive' % (after[0],))
    print('    InClassBody().as_a_method()  = %r   <-- did NOT arrive' % (after[1],))
    print('    enclosing()                  = %r  <-- DID arrive' % (after[2],))
    print('    literal_default()            = %r' % (after[3],))

    check('V1  a capture in the MODULE body is FROZEN: exactly the old value',
          after[0] == 1.0,
          'at_module_level() == 1.0 exactly after K was set to 99.0 -- '
          'equality, not a tolerance: it is the same object')
    check('V2  a capture in a FUNCTION body is LIVE: exactly the new value',
          after[2] == 99.0,
          'enclosing() == 99.0 exactly -- the inner `def` re-executed and '
          're-evaluated K')
    check('V3  a capture in a CLASS body is FROZEN, not live',
          after[1] == 1.0,
          'as_a_method() == 1.0 exactly.  THIS IS WHAT THE FIRST DRAFT OF '
          'THIS MODULE GOT WRONG by classifying on nesting depth: a class '
          'body is nested and executes once anyway')
    check('V4  so the criterion is "enclosed by a FUNCTION body", not depth',
          after[0] == after[1] != after[2],
          'the module-body and class-body arms agree with each other and '
          'disagree with the function-body arm; depth would have split them '
          'the other way, 1 against 2')
    check('V5  control: a literal default is untouched by the rebinding',
          after[3] == 7.0 == before[3],
          'so a census reporting it would be reporting a non-capture')

    # The getattr confirmation in classify() exists to drop a function that
    # has been replaced since import.  Any fix to the nested case must keep
    # that, so show it does something real.
    original = m.at_module_level
    m.at_module_level = lambda x=m.K: x
    check('V6  control: the getattr confirmation has a real job',
          m.at_module_level() == 99.0 and original() == 1.0,
          'a rebound module attribute no longer carries the captured '
          'default, so confirming static sites at run time is correct and '
          'must be preserved by any fix to the nested case')
    m.at_module_level = original

    # And the class-body capture is the case where the accident does not hold.
    check('V7  a CLASS-body capture is also invisible to getattr(module, ...)',
          not isinstance(getattr(m, 'as_a_method', None), types.FunctionType)
          and after[1] == 1.0,
          'frozen in fact, and dropped by the same line that drops the live '
          'nested one -- so classify() would UNDERCOUNT a real inert knob of '
          'this shape.  The repository has none today; that is luck, not a '
          'property of the mechanism')
    print()


# ---------------------------------------------------------------------------
# SECTION 2 -- the arbiter: a census that records the enclosing SCOPE KIND
# ---------------------------------------------------------------------------

def census_with_scope(paths):
    """Every `def f(..., p=GLOBAL)` in `paths`, with the scope that carries it.

    Yields (path, lineno, qualname, param, const, in_function_body).
    `in_function_body` is True iff some enclosing scope between the `def` and
    the module is a function body -- the condition under which the `def`
    statement re-executes and the capture is LIVE.  A `def` reached through
    class bodies only is FROZEN, because a class body, like the module body,
    runs once at import.

    Descends with explicit recursion rather than ast.walk, because the
    enclosing scope chain is exactly what ast.walk discards and it is the
    whole question.
    """
    rows = []
    for path in paths:
        try:
            src = open(path, encoding='utf-8', errors='replace').read()
            tree = ast.parse(src)
        except Exception:
            continue
        module_globals = set()
        for node in tree.body:
            if isinstance(node, ast.Assign):
                for tgt in node.targets:
                    if isinstance(tgt, ast.Name):
                        module_globals.add(tgt.id)
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                module_globals.add(node.target.id)

        def descend(node, trail, in_func):
            for child in ast.iter_child_nodes(node):
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    a = child.args
                    pairs = list(zip(
                        [x.arg for x in a.args[len(a.args) - len(a.defaults):]],
                        a.defaults))
                    pairs += [(x.arg, d) for x, d
                              in zip(a.kwonlyargs, a.kw_defaults) if d]
                    qual = '.'.join(trail + [child.name])
                    for pname, dflt in pairs:
                        for sub in ast.walk(dflt):
                            if isinstance(sub, ast.Name) and sub.id in module_globals:
                                rows.append((path, child.lineno, qual, pname,
                                             sub.id, in_func))
                                break
                    descend(child, trail + [child.name], True)
                elif isinstance(child, ast.ClassDef):
                    descend(child, trail + [child.name], in_func)
                else:
                    descend(child, trail, in_func)

        descend(tree, [], False)
    return rows


_CONTROL_SRC = '''\
A = 1.0
B = 2.0


def one(x=A): return x
def two(x=A, y=B): return x, y


class Holder:
    def method(self, x=A): return x

    class Inner:
        def deeper(self, y=B): return y


def outer():
    def inner_a(x=A): return x

    def inner_b(y=B): return y

    class LocalClass:
        def method_in_local_class(self, x=A): return x
    return inner_a(), inner_b(), LocalClass().method_in_local_class()


def literal_only(x=3.0, y="s"): return x, y
'''


def section2_arbiter_controls():
    print('-' * 74)
    print('SECTION 2 -- positive and negative controls on the arbiter itself')
    print("09-30's rule: a census that cannot fail certifies nothing.")
    print('-' * 74)
    d = tempfile.mkdtemp(prefix='frozen_default_control_')
    path = os.path.join(d, 'control.py')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(_CONTROL_SRC)
    rows = census_with_scope([path])
    frozen = [r for r in rows if not r[5]]
    live = [r for r in rows if r[5]]
    for r in sorted(rows, key=lambda r: r[1]):
        print('    %-7s %-38s %-4s <- %s'
              % ('LIVE' if r[5] else 'FROZEN', r[2], r[3], r[4]))

    # Planted by construction.  FROZEN (reached through module and class
    # bodies only): one.x, two.x, two.y, Holder.method.x, Holder.Inner.deeper.y
    # = 5.  LIVE (some enclosing scope is a function body): outer.inner_a.x,
    # outer.inner_b.y, and outer.LocalClass.method_in_local_class.x -- the
    # last one is a class body INSIDE a function body, so its `class`
    # statement re-executes on every call of outer() and the capture is live.
    check('C1 the arbiter finds exactly the planted FROZEN captures',
          len(frozen) == 5,
          'expected 5 (one.x, two.x, two.y, Holder.method.x, '
          'Holder.Inner.deeper.y), got %d' % len(frozen))
    check('C2 the arbiter finds exactly the planted LIVE captures',
          len(live) == 3,
          'expected 3 (inner_a.x, inner_b.y, and the method of a class '
          'defined INSIDE outer()), got %d' % len(live))
    check('C3 a method in a top-level class body is FROZEN',
          any(r[2] == 'Holder.method' and not r[5] for r in rows),
          'the case the first draft of this module classified as live')
    check('C4 a method nested two class bodies deep is still FROZEN',
          any(r[2] == 'Holder.Inner.deeper' and not r[5] for r in rows),
          'depth 2 and frozen: depth carries no information here')
    check('C5 a method of a class defined inside a function is LIVE',
          any(r[2] == 'outer.LocalClass.method_in_local_class' and r[5]
              for r in rows),
          'the `class` statement itself re-executes, so the capture does '
          'too -- a depth rule and a class rule both get this wrong')
    check('C6 the arbiter reports NO literal default',
          all(not r[2].startswith('literal_only') for r in rows),
          'literal_only(x=3.0, y="s") captures no module-level name')
    print()
    return rows


# ---------------------------------------------------------------------------
# SECTION 3 -- reconcile the standing censuses, to zero residue
# ---------------------------------------------------------------------------

def section3_reconciliation():
    print('-' * 74)
    print('SECTION 3 -- the two standing censuses, reconciled over one path set')
    print('-' * 74)
    import graphene_figure_provenance_audit as fpa
    import graphene_mutation_arrival_probe as probe

    all_paths = sorted(glob.glob('*.py'))
    shared = [p for p in all_paths
              if os.path.splitext(os.path.basename(p))[0] not in PROBE_EXCLUDES]
    print('  python files in repository                : %d' % len(all_paths))
    print("  minus mutation_arrival_probe's exclusions : %d" % len(shared))
    print('    excluded: %s' % ', '.join(PROBE_EXCLUDES))

    mine = census_with_scope(shared)
    frozen_rows = [r for r in mine if not r[5]]
    live_rows = [r for r in mine if r[5]]
    fpa_rows = fpa.census(shared)

    probe_frozen = 0
    probe_by_name = {}
    for path, name, verdict, detail in probe.repo_wide_frozen_census(shared):
        if verdict == 'IMPORT_FAILED':
            continue
        n_frozen = int(str(detail).split()[0])
        probe_frozen += n_frozen
        probe_by_name[(path, name)] = n_frozen

    print()
    print('  %-52s %s' % ('instrument', 'count'))
    print('  %-52s %s' % ('-' * 52, '-----'))
    print('  %-52s %5d   syntactic: every `def f(p=GLOBAL)`'
          % ('graphene_figure_provenance_audit.census', len(fpa_rows)))
    print('  %-52s %5d   semantic: sites confirmed at run time'
          % ('mutation_arrival_probe confirmed frozen sites', probe_frozen))
    print('  %-52s %5d   this module, not enclosed by a function'
          % ('arbiter: FROZEN captures (the inert knobs)', len(frozen_rows)))
    print('  %-52s %5d   this module, enclosed by a function'
          % ('arbiter: LIVE captures (not the fault)', len(live_rows)))

    check('R1 the arbiter reproduces the SYNTACTIC census exactly',
          len(mine) == len(fpa_rows),
          'arbiter total %d vs figure_provenance %d -- two independent AST '
          'passes over the same files' % (len(mine), len(fpa_rows)))
    check('R2 the arbiter reproduces the SEMANTIC census exactly',
          len(frozen_rows) == probe_frozen,
          'arbiter FROZEN %d vs mutation_arrival_probe confirmed %d'
          % (len(frozen_rows), probe_frozen))
    check('R3 the discrepancy is EXACTLY the live captures -- zero residue',
          len(fpa_rows) - probe_frozen == len(live_rows),
          '%d - %d = %d, and the arbiter finds %d live capture(s).  A '
          'reconciliation with a remainder would mean a third cause is '
          'present and unaccounted for.'
          % (len(fpa_rows), probe_frozen, len(fpa_rows) - probe_frozen,
             len(live_rows)))
    # R2 compares TOTALS, and two totals can agree while the rows behind
    # them disagree in compensating directions.  R4 is the same comparison
    # per (file, constant), which is the form that can fail on a single
    # misclassified row.
    arb_by_name = {}
    for r in frozen_rows:
        arb_by_name[(r[0], r[4])] = arb_by_name.get((r[0], r[4]), 0) + 1
    mismatches = []
    for key in sorted(set(arb_by_name) | set(probe_by_name)):
        a = arb_by_name.get(key, 0)
        b = probe_by_name.get(key, 0)
        if a != b:
            mismatches.append((key, a, b))
    check('R4 and they agree PER (file, constant), not only in total',
          not mismatches,
          '%d name(s) compared, %d disagreement(s)%s'
          % (len(set(arb_by_name) | set(probe_by_name)), len(mismatches),
             '' if not mismatches else ': ' + '; '.join(
                 '%s.%s arbiter %d vs probe %d' % (os.path.basename(k[0]),
                                                   k[1], a, b)
                 for k, a, b in mismatches)))

    n_class_body = sum(1 for r in mine if not r[5] and '.' in r[2])
    print()
    print('  class-body captures in the repository today : %d' % n_class_body)
    print('    reported, not asserted.  This is the shape that is FROZEN and')
    print('    invisible to getattr(module, name) at the same time, so it is')
    print('    the shape on which classify()\'s accident fails.  There are')
    print('    none today; that is luck, not a property of the mechanism.')

    print()
    print('  the LIVE capture(s) -- counted by figure_provenance, dropped')
    print('  silently by mutation_arrival_probe, and live in fact:')
    if not live_rows:
        print('    (none)')
    for path, lineno, qual, pname, const, _ in live_rows:
        print('    %-46s line %d' % (path, lineno))
        print('        %s(%s=%s)' % (qual, pname, const))

    print()
    print('  per-file split of the %d syntactic captures:' % len(mine))
    print('    %-46s %7s %7s %7s' % ('file', 'total', 'FROZEN', 'LIVE'))
    by_file = {}
    for r in mine:
        by_file.setdefault(r[0], []).append(r)
    for path in sorted(by_file):
        rs = by_file[path]
        print('    %-46s %7d %7d %7d'
              % (path, len(rs), sum(1 for r in rs if not r[5]),
                 sum(1 for r in rs if r[5])))
    return frozen_rows, live_rows


# ---------------------------------------------------------------------------
# SECTION 4 -- the drop, located in classify() rather than inferred
# ---------------------------------------------------------------------------

def section4_locate_the_drop(live_rows):
    print()
    print('-' * 74)
    print('SECTION 4 -- the silent drop, demonstrated on the real site')
    print('-' * 74)
    import graphene_mutation_arrival_probe as probe

    if not live_rows:
        print('  no live capture in the repository today; nothing to show.')
        return
    path, lineno, qual, pname, const, _ = live_rows[0]
    modname = os.path.splitext(os.path.basename(path))[0]
    mod = importlib.import_module(modname)
    leaf = qual.split('.')[-1]

    static = probe._frozen_sites_static(mod, const)
    info = probe.classify(mod, const)
    static_names = sorted({s[0] for s in static})
    confirmed_names = sorted({c[0] for c in info['frozen_sites']})
    dropped = set(static_names) - set(confirmed_names)

    print('  module   : %s' % modname)
    print('  constant : %s' % const)
    print('  site     : %s(%s=%s) at line %d' % (qual, pname, const, lineno))
    print('  _frozen_sites_static finds %d site(s): %s'
          % (len(static), ', '.join(static_names)))
    print('  classify() confirms        %d site(s): %s'
          % (len(info['frozen_sites']), ', '.join(confirmed_names)))
    print('  dropped                   : %s' % ', '.join(sorted(dropped)))
    print('  verdict reported          : %s' % info['verdict'])

    check('D1 the static pass sees the site',
          leaf in static_names,
          '%s is found in the source by _frozen_sites_static' % leaf)
    check('D2 classify() drops exactly that site, and only it',
          dropped == {leaf},
          'dropped = %s, expected {%r}' % (sorted(dropped), leaf))
    check('D3 the drop is by module-attribute resolution, as diagnosed',
          not isinstance(getattr(mod, leaf, None), types.FunctionType),
          'getattr(%s, %r) is a %s -- not a FunctionType, so the '
          'confirmation loop `continue`s with no row and no reason recorded'
          % (modname, leaf, type(getattr(mod, leaf, None)).__name__))
    print()
    print('  Consequence.  This dropped site happens to deserve dropping: it')
    print('  is live, as Section 1 V2 measured.  But it was dropped for being')
    print('  unresolvable, not for being live, and nothing is reported either')
    print('  way.  Section 1 V7 shows the same line drops a CLASS-body')
    print('  capture, which is frozen -- so the next case of this shape is')
    print('  undercounted, and an instrument that silently omits a true site')
    print('  has the same signature as the fault it censuses.')


# ---------------------------------------------------------------------------
# SECTION 5 -- this module censuses itself
# ---------------------------------------------------------------------------

def section5_self_check():
    print()
    print('-' * 74)
    print('SECTION 5 -- the arbiter censuses itself')
    print('-' * 74)
    me = os.path.basename(__file__)
    rows = census_with_scope([me])
    print('  captures found in %s: %d' % (me, len(rows)))
    for r in rows:
        print('    %-7s %s(%s=%s)'
              % ('LIVE' if r[5] else 'FROZEN', r[2], r[3], r[4]))
    check('S1 this module introduces no capture of either class',
          not rows,
          'PROBE_EXCLUDES is read by name at the one place it is used and '
          'never as a parameter default -- a module that censuses a pattern '
          'and then contains it would be the 2026-09-28 band-structure fault '
          'in a new costume')


def main():
    print('=' * 74)
    print('FROZEN-DEFAULT SEMANTICS -- 2026-10-16')
    print('A captured default is frozen unless a FUNCTION body encloses its def.')
    print('Reconciling three censuses that have disagreed since 2026-10-04.')
    print('=' * 74)
    print()
    section1_semantics()
    section2_arbiter_controls()
    frozen_rows, live_rows = section3_reconciliation()
    section4_locate_the_drop(live_rows)
    section5_self_check()

    print()
    print('=' * 74)
    print('SUMMARY')
    print('  syntactic captures in the repository today : %d'
          % (len(frozen_rows) + len(live_rows)))
    print('    FROZEN -- the genuine inert knobs          : %d' % len(frozen_rows))
    print('    LIVE   -- re-evaluated, not the fault      : %d' % len(live_rows))
    print('  On 2026-10-04 all three censuses read zero.  The count of the')
    print('  FAULT is %d, not %d, and this is the first time the repository has'
          % (len(frozen_rows), len(frozen_rows) + len(live_rows)))
    print('  had that number for a stated reason rather than by coincidence.')
    print('  NOT done here: converting the %d, and patching classify() so it'
          % len(frozen_rows))
    print('  reports a scope-based class instead of dropping what it cannot')
    print('  resolve.  Both are named in AUTOMATION_LOG as the next items.')
    print('  Validations: %d passed, %d failed' % (passed, failed))
    print('=' * 74)
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
