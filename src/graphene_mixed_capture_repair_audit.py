#!/usr/bin/env python3
"""
graphene_mixed_capture_repair_audit.py -- 2026-10-27

MEASURES the repair that this session's first two commits made, against the
PRE-REPAIR SOURCE, and refuses to take either side of it on trust.

WHAT WAS REPAIRED
-----------------
The 2026-10-26 session regenerated every committed transcript and found that
five module-level constants were MIXED -- some readers live, some captured in
a function's `__defaults__`:

    graphene_chord_decay_length.py     D1_CHORD   3 frozen / 9 live
                                       D2_CHORD   3 frozen / 4 live
                                       D_SPAN     1 frozen / 3 live
    graphene_capture_surface_guard.py  PERTURBATION  2 frozen / 1 live
                                       PROBE_V       2 frozen / 1 live

`graphene_mutation_arrival_probe.py` calls MIXED "the worst class, not a
middle one": a rebind of a MIXED name moves the live readers and freezes the
captured ones, so where the live readers are an expected-value expression and
the captured ones are the measurement, the mutation manufactures the
disagreement the audit then reports.  That is not a hypothetical here.
`PERTURBATION` and `PROBE_V` are the two names EVERY counterfactual in this
repository is taken with, so an oracle/measurement split on them is a defect
in the instrument rather than in any model the instrument measures.

WHY A SEPARATE SUITE, AND WHY IT READS THE OLD BLOBS
----------------------------------------------------
A repair that is only asserted is a comment.  The repository's precedent is
Section 6 of the arrival probe, which measures the T_HOP fix on the before
and after states rather than claiming it; this module does the same for the
2026-10-27 repair, and for the same reason: the number that proves a repair
worked is the number the BROKEN code produced.

The pre-repair state is reached by BLOB HASH, not by a revision expression --
`HEAD~1` names a position, not a state, and Section 8 of the arrival probe
makes that a repository-wide guard.  The two blobs below are the exact
pre-repair sources, read out of the object store and imported under distinct
module names so that both states are live in one process at once.

WHAT IT CHECKS
--------------
  B1-B5  the five names were MIXED before, with the exact frozen/live counts
         the 2026-10-26 transcript recorded.  A repair measured against a
         remembered "it was broken" is measured against nothing.
  A1-A5  the five names are LIVE after.
  M1-M6  what the capture actually did to Eq. (D), MEASURED on D1_CHORD --
         and NOT what this suite was written predicting.  The prediction was
         the textbook MIXED signature: the live term moves, the captured term
         freezes.  M2 failed, twice, and both failures are kept in the
         section's own docstring because they are the result.  The captured
         endpoint sat in a function whose BODY re-centres the polynomial at a
         LIVE read of the same constant, so under a rebind the before-state
         averaged L' over the OLD interval while evaluating a polynomial
         centred on the NEW d1: a number belonging to NO consistent choice of
         d1.  Incoherent, not frozen -- and strictly worse, because a frozen
         number is at least one the model once produced.  The checks are
         therefore about COHERENCE: under a rebind to v, the default path must
         agree BITWISE with the explicit call at d1=v.  Before, it agrees with
         neither v (off by 1.617099e-02) nor the pre-rebind value (off by
         3.838368e-02); after, it agrees with v to the last bit.
  E1-E3  the repair is EXACTLY behaviour-preserving: bitwise-equal output
         between the two states on an unmutated grid, and the default path
         bitwise equal to the explicit-argument path in the repaired state.
         E3 is the sentinel's own negative control.
  C1-C2  must-change and identity controls on the repaired default path: a
         rebind must arrive (C1) and a rebind to the same value must leave
         the result bitwise unchanged (C2).  Without C1 this suite would
         pass against a sentinel that silently ignored the module attribute.
  P1-P2  planted controls.  A classifier that can no longer SAY MIXED would
         pass A1-A5 on an empty tree, so MIXED is planted and required, and
         the before-state import is required to be a different object from
         the repaired module rather than the same file read twice.
  R1-R2  the residual: exactly ONE capture site remains repo-wide, and it is
         the one `graphene_frozen_default_semantics.py` names as its own only
         real demonstration site.  See the closing note: the arrival probe's
         "no module-level constant is captured as a default anywhere" cannot
         go green while that suite requires a real live-by-scope site to
         exist, and the two standing requirements are in direct conflict.

         RESOLVED 2026-10-28, by option (i) of the two this file's closing
         note names, and the conflict above is kept as written because it is
         the diagnosis the repair was built from.  The arrival probe's 7b is
         now three checks: a SEMANTIC one with no exemptions (green), a
         SYNTACTIC one governed by a named registry (green, one entry), and a
         HYGIENE one that fails if the registry stops describing the tree.
         R1-R2 are unaffected -- the residual is still one site and still the
         same site; what changed is that the probe no longer reports it as a
         failure without saying which site or why.
"""

import ast
import glob
import importlib
import importlib.util
import io
import contextlib
import os
import subprocess
import sys

import graphene_mutation_arrival_probe as map_

# The pre-repair blobs, by content-addressed id.  A blob hash rather than a
# revision expression, deliberately: see the module docstring and Section 8
# of the arrival probe.
BLOB_CHORD_BEFORE = '4ef1e3c837c91c5d268fdc7605fcfb040d066fa7'
BLOB_CSG_BEFORE = '8f74a6493fb5510435d424a1de2d009efcae7d05'

# The census the 2026-10-26 transcript recorded, as (module, name, frozen, live).
# Module-level and read at call time (2026-09-30 mutation-arrival rule).
BEFORE_CENSUS = [
    ('graphene_chord_decay_length', 'D1_CHORD', 3, 9),
    ('graphene_chord_decay_length', 'D2_CHORD', 3, 4),
    ('graphene_chord_decay_length', 'D_SPAN', 1, 3),
    ('graphene_capture_surface_guard', 'PERTURBATION', 2, 1),
    ('graphene_capture_surface_guard', 'PROBE_V', 2, 1),
]

# The one site deliberately NOT converted, and the suite that requires it.
RESIDUAL_SITE = ('graphene_velocity_saturation_model.py', '_R_harmonic', 'N_QUAD')
RESIDUAL_OWNER = 'graphene_frozen_default_semantics.py'

# The pinned blobs are imported under this prefix, never under the real module
# name: see main().  `_live_readers` compares `fn.__module__` against
# `module.__name__`, so a prefixed name is still classified correctly.
BEFORE_PREFIX = 'pinned_before_'

# The grid E1/E2 compare the two states on.  Inside D_SPAN and inside the
# admissible region for both shapes, so neither state refuses.
GRID_D = (2.5, 3.0, 3.3, 4.0, 4.6, 5.0)
GRID_SHAPES = ((-0.30, 0.02), (-0.55, 0.10), (0.20, -0.01))

# The mutation M1-M4 and C1 apply.  Not 3.3 * anything: a factor would make
# "the measurement moved" ambiguous if a reader happened to be scale-free.
D1_MUTATED = 3.9

_PASS, _FAIL = [], []


def check(label, ok, detail=''):
    (_PASS if ok else _FAIL).append(label)
    print('  [%s] %s' % ('PASS' if ok else 'FAIL', label))
    if detail:
        for line in str(detail).rstrip('\n').split('\n'):
            print('        ' + line)
    return bool(ok)


def hexof(x):
    return float(x).hex()


# ---------------------------------------------------------------------------
# Materialising a pinned blob as an importable module
# ---------------------------------------------------------------------------
def blob_source(blob):
    p = subprocess.run(['git', 'cat-file', 'blob', blob],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError('git cat-file %s failed: %s' % (blob, p.stderr.strip()))
    return p.stdout


def import_blob(blob, modname, destdir):
    """Import a pinned blob under `modname`, with its own __file__ on disk.

    The file must exist on disk because `classify()` reads the source through
    `inspect.getsource(module)`; an exec of a string would classify as having
    no captured defaults at all, which is the answer this suite exists to
    distinguish from a real one.
    """
    path = os.path.join(destdir, modname + '.py')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(blob_source(blob))
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[modname] = mod
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# The repo-wide syntactic capture census, independent of the arrival probe
# ---------------------------------------------------------------------------
def syntactic_capture_sites(paths=None):
    """(path, function, name) for every module-level constant captured as a
    default.  A second, independent AST pass: R1 compares its count against
    the arrival probe's own census, and two passes that agree is evidence
    where one pass is an assertion."""
    out = []
    for path in sorted(paths or glob.glob('*.py')):
        try:
            tree = ast.parse(open(path, encoding='utf-8', errors='replace').read())
        except Exception:
            continue
        gnames = {t.id for n in tree.body if isinstance(n, ast.Assign)
                  for t in n.targets if isinstance(t, ast.Name)}
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            dflts = list(node.args.defaults) + [d for d in node.args.kw_defaults if d]
            for d in dflts:
                for sub in ast.walk(d):
                    if isinstance(sub, ast.Name) and sub.id in gnames:
                        out.append((path, node.name, sub.id))
    return out


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------
def section_before(before_mods):
    print()
    print('SECTION B -- the five names, classified in the PRE-REPAIR SOURCE')
    print('  %-34s %-14s %-8s %s' % ('module (pinned blob)', 'name', 'verdict',
                                     'frozen / live'))
    for i, (modname, name, nfrozen, nlive) in enumerate(BEFORE_CENSUS, 1):
        mod = before_mods[modname]
        info = map_.classify(mod, name)
        # frozen_sites + frozen_unresolved is what the verdict is taken on
        # (classify's `frozen_total`); counting only the first would
        # undercount exactly the site class the 2026-10-16 amendment named.
        got_f = len(info['frozen_sites']) + len(info['frozen_unresolved'])
        got_l = len(info['live_readers'])
        print('  %-34s %-14s %-8s %d / %d'
              % (modname, name, info['verdict'], got_f, got_l))
        check('B%d %s.%s was MIXED before the repair' % (i, modname, name),
              info['verdict'] == map_.CLASS_MIXED,
              'verdict=%s' % info['verdict'])
        check('B%d\'  ... with the 2026-10-26 counts (%d frozen / %d live)'
              % (i, nfrozen, nlive),
              (got_f, got_l) == (nfrozen, nlive),
              'measured %d frozen / %d live' % (got_f, got_l))


def section_after():
    print()
    print('SECTION A -- the same five names in the REPAIRED tree')
    for i, (modname, name, _, _) in enumerate(BEFORE_CENSUS, 1):
        mod = importlib.import_module(modname)
        info = map_.classify(mod, name)
        print('  %-34s %-14s %-8s %d live reader(s)'
              % (modname, name, info['verdict'], len(info['live_readers'])))
        check('A%d %s.%s is LIVE after the repair' % (i, modname, name),
              info['verdict'] == map_.CLASS_LIVE,
              'verdict=%s' % info['verdict'])


def section_manufactured_disagreement(before_mods):
    """Eq. (D) = dlogP(d) - mean_dlogP(), before and after, under a rebind.

    THE FIRST DRAFT OF THIS SECTION WAS WRONG, AND THE RUN SAID SO.  It
    predicted that mutating D1_CHORD would move `dlogP(d)` and FREEZE
    `mean_dlogP()`, the textbook MIXED signature of a frozen measurement
    against a moved oracle.  M2 failed: the before-state `mean_dlogP()` moved
    too.  The reason is worth more than the prediction was.

    `mean_dlogP` captured d1 as a default, but its BODY calls `poly_P`, and
    `poly_P` reads D1_CHORD LIVE -- it re-centres the polynomial at
    `d - D1_CHORD`.  So under a rebind the before-state averaged L' over the
    OLD interval [d1_old, d2] while evaluating a polynomial re-centred about
    the NEW d1.  That is not a frozen number and not a stale number: it is a
    value that corresponds to NO consistent choice of d1 at all, and it is
    strictly worse than the frozen case this section was written expecting,
    because a frozen number is at least a number the model once produced.

    So the check is restated as the exact statement that distinguishes the
    two states, which is COHERENCE and not movement: under a rebind of
    D1_CHORD to v, the default path must agree bitwise with the explicit call
    that passes d1=v.  Before the repair it agrees with neither v nor the old
    d1; after, it agrees with v exactly.
    """
    print()
    print('SECTION M -- what the capture actually did to Eq. (D), MEASURED')
    before = before_mods['graphene_chord_decay_length']
    after = importlib.import_module('graphene_chord_decay_length')
    d, (p1, p2) = 4.0, GRID_SHAPES[0]
    out = {}
    for tag, mod in (('before', before), ('after', after)):
        old = mod.D1_CHORD
        base_live = float(mod.dlogP(d, p1, p2))
        base_mean = float(mod.mean_dlogP(p1, p2))
        base_eqd = float(mod.discrepancy(d, p1, p2))
        try:
            setattr(mod, 'D1_CHORD', D1_MUTATED)
            mut_live = float(mod.dlogP(d, p1, p2))
            mut_mean = float(mod.mean_dlogP(p1, p2))
            mut_eqd = float(mod.discrepancy(d, p1, p2))
            # The COHERENT value under the rebind: the attribute and the
            # endpoint agree, both at the new d1.  Its partner -- the coherent
            # value at the OLD d1 -- is `base_mean`, taken before the rebind
            # when attribute and endpoint also agreed.
            #
            # A SECOND WRONG DRAFT, recorded because the measurement is what
            # corrected it: this pair was first computed as two explicit calls
            # taken WITH THE REBIND IN PLACE, which made the `old` one read
            # `mean_dlogP(p1, p2, old, D2)` against a polynomial already
            # re-centred at the new D1_CHORD -- the same incoherent quantity
            # the check was trying to distinguish itself from.  It therefore
            # agreed with the before-state default path to 0.000000e+00 and
            # M2 failed.  A reference value has to be taken in the state it
            # claims to describe.
            coherent_new = float(mod.mean_dlogP(p1, p2, D1_MUTATED, mod.D2_CHORD))
        finally:
            setattr(mod, 'D1_CHORD', old)
        out[tag] = dict(base_live=base_live, base_mean=base_mean,
                        base_eqd=base_eqd, mut_live=mut_live, mut_mean=mut_mean,
                        mut_eqd=mut_eqd, coherent_new=coherent_new,
                        coherent_old=base_mean)
        print('  %-6s D1_CHORD %.4f -> %.4f' % (tag, old, D1_MUTATED))
        print('         dlogP(d)       [live reader]   %+.12e -> %+.12e'
              % (base_live, mut_live))
        print('         mean_dlogP()   [default path]  %+.12e -> %+.12e'
              % (base_mean, mut_mean))
        print('         discrepancy()  [Eq. (D)]       %+.12e -> %+.12e'
              % (base_eqd, mut_eqd))
        print('         the two values that are COHERENT -- attribute and '
              'endpoint agreeing:')
        print('           d1 = %.4f everywhere : %+.12e   (after the rebind)'
              % (D1_MUTATED, coherent_new))
        print('           d1 = %.4f everywhere : %+.12e   (before it)'
              % (old, base_mean))

    b, a = out['before'], out['after']
    check('M1 BEFORE: the rebind moved the live reader dlogP(d)',
          hexof(b['base_live']) != hexof(b['mut_live']))
    check('M2 BEFORE: the default path matches NEITHER coherent value -- it is '
          'incoherent, not merely frozen',
          hexof(b['mut_mean']) != hexof(b['coherent_new'])
          and hexof(b['mut_mean']) != hexof(b['coherent_old']),
          'default %+.12e\nnew      %+.12e  (off by %.6e)\nold      %+.12e  '
          '(off by %.6e)'
          % (b['mut_mean'], b['coherent_new'],
             abs(b['mut_mean'] - b['coherent_new']), b['coherent_old'],
             abs(b['mut_mean'] - b['coherent_old'])))
    check('M3 AFTER: the default path equals the NEW-interval value BITWISE',
          hexof(a['mut_mean']) == hexof(a['coherent_new']),
          'default and explicit agree to the last bit under the rebind')
    check('M4 AFTER: and it is not the old-interval value either -- M3 is not '
          'the rebind being ignored',
          hexof(a['mut_mean']) != hexof(a['coherent_old']))
    check('M5 the two states respond to the SAME rebind differently in '
          'Eq. (D)',
          hexof(abs(b['mut_eqd'] - b['base_eqd']))
          != hexof(abs(a['mut_eqd'] - a['base_eqd'])),
          'before %.12e vs after %.12e'
          % (abs(b['mut_eqd'] - b['base_eqd']), abs(a['mut_eqd'] - a['base_eqd'])))
    check('M6 and the unmutated Eq. (D) is identical in both states',
          hexof(b['base_eqd']) == hexof(a['base_eqd']),
          'the repair moves no published number: Section E is the repo-wide form')


def section_exact(before_mods):
    print()
    print('SECTION E -- the repair is EXACTLY behaviour-preserving')
    before = before_mods['graphene_chord_decay_length']
    after = importlib.import_module('graphene_chord_decay_length')

    worst, n = None, 0
    for (p1, p2) in GRID_SHAPES:
        for d in GRID_D:
            n += 1
            b, a = float(before.discrepancy(d, p1, p2)), float(after.discrepancy(d, p1, p2))
            if hexof(b) != hexof(a):
                worst = (d, p1, p2, b, a)
    check('E1 discrepancy() is BITWISE equal in both states over %d (d, shape) '
          'points' % n, worst is None,
          'first mismatch at %r' % (worst,) if worst else
          'no mutation applied: a repair that changes a published number is a '
          'different change')

    worst2, n2 = None, 0
    for (p1, p2) in GRID_SHAPES:
        for d in GRID_D:
            n2 += 1
            dflt = float(after.discrepancy(d, p1, p2))
            expl = float(after.discrepancy(d, p1, p2, after.D1_CHORD, after.D2_CHORD))
            if hexof(dflt) != hexof(expl):
                worst2 = (d, p1, p2, dflt, expl)
    check('E2 the default path equals the EXPLICIT-argument path bitwise, %d '
          'points' % n2, worst2 is None,
          'first mismatch at %r' % (worst2,) if worst2 else
          'the sentinel substitutes the module attribute and nothing else')

    # E3 is E2's negative control: it must be POSSIBLE for the two paths to
    # disagree, or E2 certifies nothing.  Pass an endpoint that is not the
    # module's and require a difference.
    e3 = (hexof(float(after.discrepancy(4.0, *GRID_SHAPES[0])))
          != hexof(float(after.discrepancy(4.0, *GRID_SHAPES[0],
                                           d1=after.D1_CHORD + 0.4))))
    check('E3 and the explicit path CAN disagree when given a different '
          'endpoint -- E2 is not vacuous', e3)
    check('E4 admissible() default span equals the explicit D_SPAN path',
          after.admissible(*GRID_SHAPES[0])
          == after.admissible(*GRID_SHAPES[0], span=after.D_SPAN))


def section_controls():
    print()
    print('SECTION C -- must-change and identity controls on the repaired path')
    after = importlib.import_module('graphene_chord_decay_length')
    d, (p1, p2) = 4.0, GRID_SHAPES[1]
    old = after.D2_CHORD
    base = float(after.discrepancy(d, p1, p2))
    try:
        after.D2_CHORD = old + 0.5
        moved = float(after.discrepancy(d, p1, p2))
    finally:
        after.D2_CHORD = old
    check('C1 MUST_CHANGE: a rebind of D2_CHORD reaches the default path',
          hexof(base) != hexof(moved),
          'base %.12e -> %.12e' % (base, moved))
    try:
        after.D2_CHORD = old
        again = float(after.discrepancy(d, p1, p2))
    finally:
        after.D2_CHORD = old
    check('C2 IDENTITY: a rebind to the SAME value leaves it bitwise '
          'unchanged', hexof(base) == hexof(again))


def section_planted(before_mods, destdir):
    print()
    print('SECTION P -- planted controls (a classifier that cannot say MIXED '
          'certifies nothing)')
    path = os.path.join(destdir, 'planted_mixed_module.py')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('PLANTED = 2.0\n'
                 'def frozen_reader(x=PLANTED):\n'
                 '    return x * 10.0\n'
                 'def live_reader():\n'
                 '    return PLANTED * 10.0\n')
    spec = importlib.util.spec_from_file_location('planted_mixed_module', path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules['planted_mixed_module'] = mod
    spec.loader.exec_module(mod)
    info = map_.classify(mod, 'PLANTED')
    check('P1 the classifier still reports MIXED on a planted mixed name',
          info['verdict'] == map_.CLASS_MIXED, 'verdict=%s' % info['verdict'])
    live = importlib.import_module('graphene_chord_decay_length')
    check('P2 the pinned before-state is a DIFFERENT object from the '
          'repaired module',
          before_mods['graphene_chord_decay_length'] is not live
          and before_mods['graphene_chord_decay_length'].__file__ != live.__file__
          and before_mods['graphene_chord_decay_length'].__name__ != live.__name__,
          'before %s  (%s)\nafter  %s  (%s)'
          % (before_mods['graphene_chord_decay_length'].__name__,
             os.path.basename(before_mods['graphene_chord_decay_length'].__file__),
             live.__name__, os.path.basename(live.__file__)))
    check('P3 and sys.modules still maps the real names to the REPAIRED '
          'modules', all(sys.modules[n].__file__ == os.path.abspath(n + '.py')
                         for n in ('graphene_chord_decay_length',
                                   'graphene_capture_surface_guard')),
          'the shadowing that made the first run of this suite report 9 '
          'failures is now a check, not a memory')


def section_residual():
    print()
    print('SECTION R -- the residual capture site')
    sites = [s for s in syntactic_capture_sites()
             if os.path.basename(s[0]) not in
             ('graphene_mutation_arrival_probe.py',
              'graphene_mixed_capture_repair_audit.py')]
    for path, fn, name in sites:
        print('  %-44s %-24s %s' % (path, fn, name))
    check('R1 exactly ONE capture site remains repo-wide (11 of 12 converted '
          'today)', len(sites) == 1, '%d site(s): %r' % (len(sites), sites))
    check('R2 and it is the site %s names as its own demonstration site'
          % RESIDUAL_OWNER,
          bool(sites) and sites[0] == RESIDUAL_SITE,
          'found %r, expected %r' % (sites[0] if sites else None, RESIDUAL_SITE))
    owner = open(RESIDUAL_OWNER, encoding='utf-8').read()
    check('R3 and %s still names it, so R2 is not a coincidence of spelling'
          % RESIDUAL_OWNER,
          RESIDUAL_SITE[1] in owner and RESIDUAL_SITE[2] in owner)


def main():
    print('=' * 78)
    print('MIXED-CAPTURE REPAIR AUDIT -- 2026-10-27')
    print('the 2026-10-26 ranked-first item, measured against the pinned '
          'pre-repair blobs')
    print('=' * 78)

    destdir = os.path.join(os.environ.get('TMPDIR', '/tmp'),
                           'mixed_capture_before_%d' % os.getpid())
    os.makedirs(destdir, exist_ok=True)
    # The repaired modules are imported FIRST and the blobs are registered
    # under PREFIXED names.  The first run of this suite registered each blob
    # under its real module name, which shadowed the repaired module in
    # sys.modules -- so every "after" lookup returned the "before" module and
    # A1-A5, M2, M4, C1 and P2 all failed.  P2 is the control that catches
    # exactly this, and it is the reason the two states are required to be
    # distinct objects with distinct __file__ rather than merely asked for.
    for real in ('graphene_chord_decay_length', 'graphene_capture_surface_guard'):
        importlib.import_module(real)
    before_mods = {
        'graphene_chord_decay_length':
            import_blob(BLOB_CHORD_BEFORE,
                        BEFORE_PREFIX + 'graphene_chord_decay_length', destdir),
        'graphene_capture_surface_guard':
            import_blob(BLOB_CSG_BEFORE,
                        BEFORE_PREFIX + 'graphene_capture_surface_guard', destdir),
    }
    print('  pinned blobs materialised:')
    print('    %s  graphene_chord_decay_length.py' % BLOB_CHORD_BEFORE)
    print('    %s  graphene_capture_surface_guard.py' % BLOB_CSG_BEFORE)

    section_before(before_mods)
    section_after()
    section_manufactured_disagreement(before_mods)
    section_exact(before_mods)
    section_controls()
    section_planted(before_mods, destdir)
    section_residual()

    print()
    print('=' * 78)
    print('TOTAL: %d passed, %d failed' % (len(_PASS), len(_FAIL)))
    print('=' * 78)
    print()
    print('THE RESIDUAL, AND A CONFLICT BETWEEN TWO STANDING REQUIREMENTS')
    print('  11 of the 12 syntactic capture sites are converted.  The 12th,')
    print('  %s in %s,' % (RESIDUAL_SITE[2], RESIDUAL_SITE[0]))
    print('  is a `def` inside a function BODY, so its default is re-evaluated')
    print('  on every call of the enclosing function and a rebind DOES arrive:')
    print('  it is live_by_scope, not frozen, and %s'
          % RESIDUAL_OWNER)
    print('  measures it as live in its own Section 1 V2.  It is also that')
    print('  suite\'s ONLY REAL demonstration site: its Section 4 D1-D5 exist')
    print('  to show the silent-drop mechanism on a site that actually exists.')
    print()
    print('  So the arrival probe\'s Section 7b check "no module-level constant')
    print('  is captured as a default anywhere" CANNOT go green while that')
    print('  suite keeps its site, and its sibling "and in particular none is')
    print('  MIXED" now can.  The two standing requirements conflict, and the')
    print('  conflict is in the VOCABULARY, not in the tree: 7b counts')
    print('  syntactic captures and names the class it cares about only in the')
    print('  second check.  Not resolved here.  Resolving it means either')
    print('  splitting 7b into FROZEN-or-MIXED (achievable) and syntactic')
    print('  (not), or moving the demonstration site into a planted fixture --')
    print('  and the second trades a real site for a synthetic one, which is')
    print('  the trade that suite was written to argue against.')
    print()
    print('  RESOLVED 2026-10-28 -- option (i), and the paragraph above is')
    print('  kept verbatim because it is the diagnosis the repair was built')
    print('  from.  Section 7b is now three checks over the same census:')
    print('    7b-i   SEMANTIC  no site FREEZES a read (FROZEN / MIXED /')
    print('                     IMPORT_FAILED).  No exemptions.  GREEN.')
    print('    7b-ii  SYNTACTIC no site is unregistered.  One registered')
    print('                     entry, carrying its class and its reason.')
    print('    7b-iii HYGIENE   the registry still describes the tree:')
    print('                     VANISHED and DRIFTED reported separately.')
    print('  Eight planted controls (C7c) in the probe, of which arm 7 is')
    print('  the load-bearing one: the REGISTERED site, planted as FROZEN,')
    print('  must still fire 7b-i.  An exemption that could silence the')
    print('  semantic check would have made the split a weakening, so')
    print('  census_semantic_faults() does not take the registry at all.')
    print('  The residual measured here is unchanged: one site, same site.')
    return 0 if not _FAIL else 1


if __name__ == '__main__':
    raise SystemExit(main())
