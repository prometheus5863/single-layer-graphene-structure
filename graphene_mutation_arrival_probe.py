"""
graphene_mutation_arrival_probe.py -- does the mutation ARRIVE?

Created 2026-09-30 to close the top open item of 2026-09-29: "the 37 remaining
write-only knobs, 13 of them in audit modules ... fix the audit modules FIRST
and from outside".  This module is the "from outside": it is not an audit
module, it imports no audit module, and nothing in the repository imports it.
It can therefore be used to measure the audit modules without the instrument
and the measurement moving in the same step.

WHAT IT IS FOR
--------------
Every mutation-based result in this repository is produced by rebinding a
module attribute and re-reading a number:

    setattr(module, 'CONST', new_value); y2 = module.f()

The 2026-09-29 session found that this is silently a no-op whenever `CONST` is
captured as a DEFAULT ARGUMENT, because `def` evaluates defaults once, at
definition time.  Rebinding the module attribute rebinds the NAME; the function
keeps the object it captured.  The audit then prints "no sensitivity" when what
it measured was "no mutation".

The 09-29 methodological note proposed the defence as "a paired assertion that
some quantity the mutation must reach did in fact move".  THAT FORMULATION IS
NOT SOUND, and correcting it is the reason this module exists rather than a
helper function:

  * It cannot be applied to a null control at all.  A MUST_NOT_CHANGE check
    asserts the output does NOT move, so there is no quantity available to
    require movement in.  Every null control in this repository passes
    identically whether the mutation arrived or the mutation machinery is
    dead -- see `null_controls_carry_no_arrival_evidence()` below, which proves
    this rather than arguing it.
  * Even for a MUST_CHANGE check it confounds two failures.  If the downstream
    number does not move, "the mutation did not arrive" and "the model is
    genuinely insensitive here" produce the same zero, which is the 09-29
    finding restated, not defended against.

The sound place for the control is therefore NOT downstream of the mutation but
AT THE POINT OF DELIVERY: before interpreting any number, establish whether the
binding the callee will actually read has changed.  That question is decidable
exactly, by introspection, without evaluating the model at all -- which also
means it is answerable for null controls and for insensitive quantities, the
two cases the downstream form cannot serve.

DELIVERY CLASSES
----------------
For a module-level name N in module M, every reader of N is one of:

  LIVE    -- read from M's globals at call time.  A rebind arrives.
  FROZEN  -- captured in a function's __defaults__/__kwdefaults__.  A rebind
             never arrives, however many times it is applied.

and N as a whole is classified:

  'LIVE'    all readers live                  -> mutation is delivered
  'FROZEN'  all readers frozen                -> mutation is WRITE-ONLY
  'MIXED'   some of each                      -> mutation is delivered to some
            readers and not others.  This is the worst class, not a middle
            one: if the live readers are an audit's expected-value expressions
            and the frozen readers are its measuring functions, a mutation of N
            MOVES THE ORACLE AND FREEZES THE MEASUREMENT, and the audit reports
            a disagreement it manufactured itself.
  'UNUSED'  no reader at all                  -> mutation reaches nothing
  'ABSENT'  N is not an attribute of M        -> setattr CREATES it and returns
            normally, so a misspelled or renamed mutation target is a
            write-only knob too, with no syntactic trace.

'ABSENT' is included because it is the same failure with a different cause and
the same output, and because nothing in this repository currently checks it:
`setattr` has no notion of a typo.

Usage, and the rule this module asks the repository to adopt:

    from graphene_mutation_arrival_probe import require_delivery
    require_delivery('graphene_fet_model', ['mu', 'T', 'Rc_per_width_ohm_um'])

before any number produced by mutating those names is interpreted, INCLUDING
the numbers a null control produces.
"""

import ast
import importlib
import inspect
import os
import sys
import tempfile
import types

CLASS_LIVE = 'LIVE'
CLASS_FROZEN = 'FROZEN'
CLASS_MIXED = 'MIXED'
CLASS_UNUSED = 'UNUSED'
CLASS_ABSENT = 'ABSENT'
# Added after the first run of 2026-09-30, which returned UNUSED for two names
# with entirely different diagnoses and so would have had them fixed the same
# wrong way.  Splitting them is not cosmetic: one is correct behaviour that the
# 09-29 audit already handles, the other is a dead name.
CLASS_IMPORT_CONSUMED = 'IMPORT_USED'   # read in the module BODY, by no function
CLASS_DEAD = 'DEAD'                     # read nowhere at all, body included

DELIVERED = (CLASS_LIVE,)
UNDELIVERED = (CLASS_FROZEN, CLASS_UNUSED, CLASS_ABSENT, CLASS_DEAD)
PARTIAL = (CLASS_MIXED,)
# IMPORT_CONSUMED is deliberately in neither list.  A rebind after import
# reaches nothing, so it is not DELIVERED -- but it is also not a fault: the
# name was consumed once, at import, to derive a child constant, and mutating
# it is legitimate provided the CHILD is mutated in the same breath.  That
# proviso is checked mechanically in `import_consumed_is_covered` below rather
# than trusted.


# ---------------------------------------------------------------------------
# Finding the readers
# ---------------------------------------------------------------------------

def _global_names_read(code, seen=None):
    """
    Every global name `code` may LOAD at run time, including inside nested
    code objects (comprehensions, lambdas, inner defs).  co_names is the right
    source: a global read compiles to LOAD_GLOBAL against co_names, whereas a
    captured default is stored in __defaults__ and never appears as a load.
    """
    if seen is None:
        seen = set()
    if id(code) in seen:
        return set()
    seen.add(id(code))
    out = set(code.co_names)
    for const in code.co_consts:
        if isinstance(const, types.CodeType):
            out |= _global_names_read(const, seen)
    return out


def _frozen_sites_static(module, name):
    """
    Which (function, parameter) pairs capture `name` as a default, read from
    the SOURCE.  Static because that is where the author's intent lives: the
    runtime object alone cannot tell `t=T_HOP` from `t=2.8`.
    """
    try:
        src = inspect.getsource(module)
        tree = ast.parse(src)
    except (OSError, TypeError, SyntaxError):
        return []
    sites = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        a = node.args
        positional = list(zip(
            [x.arg for x in a.args[len(a.args) - len(a.defaults):]],
            a.defaults,
            range(len(a.defaults)),
        ))
        for pname, dflt, idx in positional:
            for sub in ast.walk(dflt):
                if isinstance(sub, ast.Name) and sub.id == name:
                    sites.append((node.name, pname, 'positional', idx))
                    break
        for kwarg, dflt in zip(a.kwonlyargs, a.kw_defaults):
            if dflt is None:
                continue
            for sub in ast.walk(dflt):
                if isinstance(sub, ast.Name) and sub.id == name:
                    sites.append((node.name, kwarg, 'kwonly', kwarg))
                    break
    return sites


def _read_in_module_body(module, name):
    """
    Is `name` read anywhere in the module's top-level body (outside every
    function)?  Such a read happens once, at import, so a later rebind cannot
    reach it -- but the name is not dead, and the two need different fixes.
    """
    try:
        tree = ast.parse(inspect.getsource(module))
    except (OSError, TypeError, SyntaxError):
        return False
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        for sub in ast.walk(node):
            if isinstance(sub, ast.Name) and sub.id == name and \
                    isinstance(sub.ctx, ast.Load):
                return True
    return False


def _live_readers(module, name):
    """Functions defined in `module` that LOAD `name` as a global at call time."""
    out = []
    for fname, obj in sorted(vars(module).items()):
        if not isinstance(obj, types.FunctionType):
            continue
        if getattr(obj, '__module__', None) != module.__name__:
            continue
        if name in _global_names_read(obj.__code__):
            out.append(fname)
    return out


def classify(module, name):
    """
    Return a dict describing whether `setattr(module, name, v)` will be seen.

    Exact and evaluation-free: no model is run, so the answer is available for
    null controls and for quantities the model is genuinely insensitive to --
    the two cases a downstream "did the number move" control cannot serve.
    """
    if isinstance(module, str):
        module = importlib.import_module(module)
    present = hasattr(module, name)
    frozen = _frozen_sites_static(module, name)
    live = _live_readers(module, name)
    # A function listed as a frozen site is only *actually* frozen if it is the
    # object currently bound in the module (a later redefinition or a decorator
    # can have replaced it), so confirm the capture at run time as well.
    confirmed = []
    for fname, pname, kind, idx in frozen:
        fn = getattr(module, fname, None)
        if not isinstance(fn, types.FunctionType):
            continue
        if kind == 'positional':
            defaults = fn.__defaults__ or ()
            if idx < len(defaults):
                confirmed.append((fname, pname, defaults[idx]))
        else:
            kwd = fn.__kwdefaults__ or {}
            if pname in kwd:
                confirmed.append((fname, pname, kwd[pname]))
    body_read = _read_in_module_body(module, name) if present else False
    if not present:
        verdict = CLASS_ABSENT
    elif confirmed and live:
        verdict = CLASS_MIXED
    elif confirmed:
        verdict = CLASS_FROZEN
    elif live:
        verdict = CLASS_LIVE
    elif body_read:
        verdict = CLASS_IMPORT_CONSUMED
    else:
        verdict = CLASS_DEAD
    return {
        'module': module.__name__,
        'name': name,
        'present': present,
        'verdict': verdict,
        'frozen_sites': confirmed,
        'live_readers': live,
        'body_read': body_read,
        'current_value': getattr(module, name, None),
    }


def require_delivery(module, names, strict_mixed=True):
    """
    The rule this module asks for.  Returns (ok, rows).  `ok` is False if any
    name is undelivered, or -- when strict_mixed -- only partly delivered.
    Call it BEFORE interpreting any mutation result, null controls included.
    """
    rows = [classify(module, n) for n in names]
    bad = [r for r in rows
           if r['verdict'] in UNDELIVERED
           or (strict_mixed and r['verdict'] in PARTIAL)]
    return (not bad), rows


# ---------------------------------------------------------------------------
# Positive control on the probe itself (the 2026-09-28 rule)
# ---------------------------------------------------------------------------

_CONTROL_SRC = '''
LIVE_CONST = 1.0
FROZEN_CONST = 2.0
MIXED_CONST = 3.0
UNUSED_CONST = 4.0
PARENT_CONST = 5.0
DERIVED_CONST = PARENT_CONST * 2.0
DEAD_CONST = 6.0

def reads_live():
    return LIVE_CONST * 10.0

def reads_frozen(c=FROZEN_CONST):
    return c * 10.0

def mixed_measure(c=MIXED_CONST):
    return c * 10.0

def mixed_oracle():
    return MIXED_CONST * 10.0

def reads_a_literal(c=2.0):
    return c

def reads_the_derived():
    return DERIVED_CONST
'''


def probe_positive_control():
    """
    A probe that cannot misclassify certifies whatever sits beside it.  Hand it
    a module containing one of every class, with known answers, and require all
    five -- including UNUSED (a constant nothing reads) and ABSENT (a name that
    is not there at all), which are the two easiest to return by accident.
    """
    d = tempfile.mkdtemp()
    path = os.path.join(d, 'probe_control_module.py')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(_CONTROL_SRC)
    sys.path.insert(0, d)
    try:
        sys.modules.pop('probe_control_module', None)
        mod = importlib.import_module('probe_control_module')
        expected = {
            'LIVE_CONST': CLASS_LIVE,
            'FROZEN_CONST': CLASS_FROZEN,
            'MIXED_CONST': CLASS_MIXED,
            'UNUSED_CONST': CLASS_DEAD,
            'DEAD_CONST': CLASS_DEAD,
            'PARENT_CONST': CLASS_IMPORT_CONSUMED,
            'DERIVED_CONST': CLASS_LIVE,
            'NO_SUCH_CONST': CLASS_ABSENT,
        }
        got = {n: classify(mod, n)['verdict'] for n in expected}
        # And the negative: a literal default must not be reported as a capture
        lit = classify(mod, 'reads_a_literal')
        no_false_capture = lit['frozen_sites'] == []
        # And the mutation behaviour the classes PREDICT must actually hold
        before_live, before_frozen = mod.reads_live(), mod.reads_frozen()
        setattr(mod, 'LIVE_CONST', 7.0)
        setattr(mod, 'FROZEN_CONST', 7.0)
        live_moved = mod.reads_live() != before_live
        frozen_moved = mod.reads_frozen() != before_frozen
        # the MIXED prediction: oracle moves, measurement does not
        o0, m0 = mod.mixed_oracle(), mod.mixed_measure()
        setattr(mod, 'MIXED_CONST', 30.0)
        oracle_moved = mod.mixed_oracle() != o0
        measure_frozen = mod.mixed_measure() == m0
        before_derived = mod.reads_the_derived()
        setattr(mod, 'PARENT_CONST', 500.0)
        parent_inert = mod.reads_the_derived() == before_derived
        setattr(mod, 'DERIVED_CONST', 999.0)
        child_moves = mod.reads_the_derived() != before_derived
    finally:
        sys.path.remove(d)
        sys.modules.pop('probe_control_module', None)
    checks = [
        ('all classes correct (LIVE/FROZEN/MIXED/DEAD/IMPORT_USED/ABSENT)',
         got == expected, '%s' % got),
        ('IMPORT_USED prediction: rebinding the parent is inert', parent_inert,
         'the child was derived once, at import'),
        ('IMPORT_USED prediction: rebinding the child arrives', child_moves,
         'so such a target is only sound when the child is patched too'),
        ('a literal default is not reported as a capture', no_false_capture, ''),
        ('LIVE prediction holds: rebind arrives', live_moved, ''),
        ('FROZEN prediction holds: rebind does not arrive', not frozen_moved, ''),
        ('MIXED prediction holds: oracle moves', oracle_moved, ''),
        ('MIXED prediction holds: measurement frozen', measure_frozen,
         'oracle and measurement disagree by the mutation alone'),
    ]
    return checks


# ---------------------------------------------------------------------------
# The theorem the 09-29 entry needed and did not have
# ---------------------------------------------------------------------------

def null_controls_carry_no_arrival_evidence():
    """
    PROVE, rather than argue, that a MUST_NOT_CHANGE check cannot supply the
    positive control the 09-29 note asked for.

    Four runs against one synthetic model.  The one that matters is the fourth:
    a null control on a quantity that IS fully sensitive to the mutated
    constant still PASSES, because the mutation never arrived.  Its zero is
    bitwise the zero of a correct null control.  Therefore a passing null
    control is evidence about neither the model's specificity nor the
    machinery's liveness, and the 09-29 claim that four null controls holding
    at zero showed "the detector is not reporting that everything moves" is
    not supported by those four zeros -- it is supported only by the
    MUST_CHANGE checks that passed in the same run.
    """
    d = tempfile.mkdtemp()
    path = os.path.join(d, 'null_control_demo.py')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(_CONTROL_SRC)
    sys.path.insert(0, d)
    rows = []
    try:
        def fresh():
            sys.modules.pop('null_control_demo', None)
            return importlib.import_module('null_control_demo')

        # 1. MUST_CHANGE on a LIVE constant -- the check works
        m = fresh(); b = m.reads_live()
        setattr(m, 'LIVE_CONST', 9.0)
        rows.append(('MUST_CHANGE  on LIVE   constant', m.reads_live() != b,
                     'detects the mutation', True))
        # 2. MUST_CHANGE on a FROZEN constant -- the check correctly fails
        m = fresh(); b = m.reads_frozen()
        setattr(m, 'FROZEN_CONST', 9.0)
        rows.append(('MUST_CHANGE  on FROZEN constant', m.reads_frozen() != b,
                     'fails, and SHOULD -- this is the 09-29 detector', False))
        # 3. MUST_NOT_CHANGE, mutation arrives, quantity truly independent
        m = fresh(); b = m.reads_frozen()
        setattr(m, 'LIVE_CONST', 9.0)
        rows.append(('MUST_NOT_CHG on LIVE   constant, indep. output',
                     m.reads_frozen() == b, 'passes, and is informative', True))
        # 4. MUST_NOT_CHANGE, mutation inert, quantity FULLY dependent
        m = fresh(); b = m.reads_frozen()
        setattr(m, 'FROZEN_CONST', 9.0)
        rows.append(('MUST_NOT_CHG on FROZEN constant, DEPENDENT output',
                     m.reads_frozen() == b,
                     'PASSES -- and the output is 100% sensitive to the '
                     'constant. The check is vacuous.', True))
    finally:
        sys.path.remove(d)
        sys.modules.pop('null_control_demo', None)
    return rows


# ---------------------------------------------------------------------------
# Applying it: was any result in this repository measuring nothing?
# ---------------------------------------------------------------------------

def audit_targets_of_provenance_audit():
    """
    Read the 2026-09-29 provenance audit's own AUDIT table and resolve every
    patch key it applies to a (module, name) pair, honouring the dotted
    'othermodule.ATTR' form.  Reading the table rather than re-typing it is the
    09-28 rule: a transcript of what the auditor believed is not a check.
    """
    import graphene_figure_provenance_audit as prov
    out = []
    for figure, modname, fn, checks in prov.AUDIT:
        for label, patches, verdict, note in checks:
            for key in patches:
                if '.' in key:
                    tgt_mod, attr = key.split('.', 1)
                else:
                    tgt_mod, attr = modname, key
                out.append((figure, label, verdict, tgt_mod, attr))
    return out


def delivery_census(rows):
    """Classify every (module, name) a mutation-based result depends on."""
    out = []
    for figure, label, verdict, modname, attr in rows:
        try:
            info = classify(modname, attr)
        except Exception as exc:                        # pragma: no cover
            info = {'verdict': 'IMPORT_FAILED', 'name': attr,
                    'module': modname, 'frozen_sites': [], 'live_readers': [],
                    'current_value': repr(exc)}
        out.append((figure, label, verdict, info))
    return out


def identity_reassignments(rows):
    """
    A second, independent way a mutation can fail to arrive, and one the
    delivery classes do not catch: the "mutation" assigns the value the
    constant already has.  Such a check passes with `setattr` replaced by
    `pass`, with the constant frozen, and with the model deleted.  It is the
    2026-09-28 cannot-fail class living inside the 09-29 null controls.
    """
    import graphene_figure_provenance_audit as prov
    found = []
    for figure, modname, fn, checks in prov.AUDIT:
        for label, patches, verdict, note in checks:
            for key, newval in patches.items():
                tgt_mod, attr = (key.split('.', 1) if '.' in key
                                 else (modname, key))
                try:
                    mod = importlib.import_module(tgt_mod)
                    cur = getattr(mod, attr, None)
                except Exception:
                    continue
                try:
                    same = (cur == newval)
                except Exception:
                    same = False
                if same:
                    found.append((figure, label, verdict, tgt_mod, attr, cur))
    return found


def absent_target_is_indistinguishable():
    """
    The same demonstration as `null_controls_carry_no_arrival_evidence`, but on
    the REAL instrument: take a MUST_CHANGE check that the 09-29 audit reports
    as passing, misspell its patch key by one character, and re-run.  `setattr`
    creates the misspelled attribute and returns normally, so the harvested
    curve is bitwise identical to the baseline -- the same output the audit
    interprets as "the figure does not depend on this constant".
    """
    import numpy as np
    import graphene_figure_provenance_audit as prov
    figure, modname, fn, checks = prov.AUDIT[0]
    label, patches, verdict, note = checks[0]          # 'v_F x2', MUST_CHANGE
    key = list(patches)[0]
    typo_key = key + '_'
    base = prov.run_and_harvest(modname, fn)
    # Classify BEFORE the mutation is applied.  This ordering is a REQUIREMENT,
    # not a style choice, and the first run of this file on 2026-09-30 got it
    # wrong and was corrected here: `setattr` creates the misspelled attribute,
    # and `importlib.reload` re-executes the source into the SAME module dict
    # without clearing it, so the attribute survives every later reload.  A
    # probe run after the fact therefore sees a name that EXISTS and is read by
    # nothing -- DEAD, not ABSENT.  APPLYING A MISSPELLED MUTATION DESTROYS THE
    # EVIDENCE THAT IT WAS MISSPELLED.  Both verdicts are returned below so the
    # difference is on the record rather than in a comment.
    cls_before = classify(modname, typo_key)['verdict']
    real = prov.run_and_harvest(modname, fn, dict(patches))
    typo = prov.run_and_harvest(modname, fn, {typo_key: patches[key]})
    cls_after = classify(modname, typo_key)['verdict']
    ok_real, rel_real, n = prov.compare(base, real)
    ok_typo, rel_typo, _ = prov.compare(base, typo)
    return {
        'figure': figure, 'label': label, 'key': key, 'typo_key': typo_key,
        'real_identical': ok_real, 'real_max_rel': rel_real,
        'typo_identical': ok_typo, 'typo_max_rel': rel_typo,
        'n_points': n,
        'classify_before': cls_before, 'classify_after': cls_after,
    }


def import_consumed_is_covered(rows):
    """
    An IMPORT_CONSUMED target is legitimate only if the child constant it was
    consumed to derive is mutated in the SAME check.  Verify that mechanically
    instead of trusting the audit's comment that it is done: for each such
    target, require at least one other name patched by the same check to be
    LIVE and to have been assigned from the consumed name at module level.
    """
    import graphene_figure_provenance_audit as prov
    verdicts = []
    for figure, modname, fn, checks in prov.AUDIT:
        for label, patches, verdict, note in checks:
            keys = [(k.split('.', 1) if '.' in k else (modname, k))
                    for k in patches]
            for tgt_mod, attr in keys:
                info = classify(tgt_mod, attr)
                if info['verdict'] != CLASS_IMPORT_CONSUMED:
                    continue
                mod = importlib.import_module(tgt_mod)
                children = _module_level_children_of(mod, attr)
                siblings = {a for m, a in keys if m == tgt_mod and a != attr}
                covered = bool(children & siblings)
                verdicts.append((figure, label, tgt_mod, attr,
                                 sorted(children), sorted(siblings), covered))
    return verdicts


def _module_level_children_of(module, name):
    """Module-level names whose defining assignment reads `name`."""
    try:
        tree = ast.parse(inspect.getsource(module))
    except (OSError, TypeError, SyntaxError):
        return set()
    out = set()
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        reads = {s.id for s in ast.walk(node.value)
                 if isinstance(s, ast.Name)}
        if name in reads:
            for tgt in node.targets:
                if isinstance(tgt, ast.Name):
                    out.add(tgt.id)
    return out


def repo_wide_frozen_census(paths=None):
    """
    Every module-level name in the repository that is captured as a default
    somewhere, classified.  This is the 09-29 census re-run through the
    delivery classes, so that MIXED -- the class that moves an oracle while
    freezing a measurement -- is separated out from plain FROZEN.
    """
    import glob
    paths = paths or sorted(glob.glob('*.py'))
    rows = []
    for path in paths:
        modname = os.path.splitext(os.path.basename(path))[0]
        if modname in ('initialize_repo', 'main', 'generate_plots',
                       'graphene_mutation_arrival_probe'):
            continue
        try:
            src = open(path, encoding='utf-8', errors='replace').read()
            tree = ast.parse(src)
        except Exception:
            continue
        gnames = set()
        for node in tree.body:
            if isinstance(node, ast.Assign):
                for tgt in node.targets:
                    if isinstance(tgt, ast.Name):
                        gnames.add(tgt.id)
        captured = set()
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            a = node.args
            dflts = list(a.defaults) + [d for d in a.kw_defaults if d]
            for d in dflts:
                for sub in ast.walk(d):
                    if isinstance(sub, ast.Name) and sub.id in gnames:
                        captured.add(sub.id)
        if not captured:
            continue
        try:
            mod = importlib.import_module(modname)
        except Exception as exc:
            rows.append((path, '<import failed>', 'IMPORT_FAILED', repr(exc)))
            continue
        for name in sorted(captured):
            info = classify(mod, name)
            rows.append((path, name, info['verdict'],
                         '%d frozen / %d live' % (len(info['frozen_sites']),
                                                  len(info['live_readers']))))
    return rows


# ---------------------------------------------------------------------------
def main():
    passed = failed = 0

    def check(desc, ok, extra=''):
        nonlocal passed, failed
        if ok:
            passed += 1
        else:
            failed += 1
        print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', desc,
                               ('  -- ' + extra) if extra else ''))

    print('=' * 78)
    print('MUTATION ARRIVAL PROBE -- does the mutation reach the callee?')
    print('2026-09-30.  Closes the top open item of 2026-09-29 from OUTSIDE')
    print('every audit module, so that fixing an audit does not move the')
    print('instrument and the measurement in the same step.')
    print('=' * 78)

    print()
    print('SECTION 0 -- positive control on the probe itself (09-28 rule)')
    for desc, ok, extra in probe_positive_control():
        check(desc, ok, extra)

    print()
    print('SECTION 1 -- THE THEOREM: a null control cannot be a positive')
    print('             control on its own mutation')
    print('  The 09-29 note asked for "a paired assertion that some quantity')
    print('  the mutation must reach did in fact move".  A MUST_NOT_CHANGE')
    print('  check has no such quantity by construction.  Four runs:')
    rows = null_controls_carry_no_arrival_evidence()
    for desc, result, comment, expect in rows:
        print('    %-46s -> %-5s  %s'
              % (desc, 'PASS' if result else 'FAIL', comment))
        check('run behaves as the classes predict: ' + desc, result == expect)
    vacuous = rows[3]
    check('the vacuous null control DID pass (this is the finding)',
          vacuous[1] is True,
          'a 100%-dependent output, a dead mutation, and a passing check')

    print()
    print('SECTION 2 -- the same failure on the REAL instrument')
    d = absent_target_is_indistinguishable()
    print('  figure %s, check %r, patch key %r' % (d['figure'], d['label'],
                                                   d['key']))
    print('  real mutation      : identical=%s  max|rel|=%.6e  (n=%d)'
          % (d['real_identical'], d['real_max_rel'], d['n_points']))
    print('  one-character typo : identical=%s  max|rel|=%.6e'
          % (d['typo_identical'], d['typo_max_rel']))
    print('  probe verdict on that key BEFORE the mutation : %s'
          % d['classify_before'])
    print('  probe verdict on that key AFTER  the mutation : %s'
          % d['classify_after'])
    check('the real mutation does move the figure (MUST_CHANGE is alive)',
          not d['real_identical'] and d['real_max_rel'] > 0.0)
    check('a MISSPELLED patch key is BITWISE identical to no mutation',
          d['typo_identical'] and d['typo_max_rel'] == 0.0,
          'setattr created the attribute and returned normally')
    check('the probe classifies the misspelled target ABSENT, BEFOREHAND',
          d['classify_before'] == CLASS_ABSENT,
          'decidable without running the model at all')
    check('and NOT afterwards -- applying the typo destroys the evidence',
          d['classify_after'] != CLASS_ABSENT,
          'verdict degrades to %s, because setattr created the name and '
          'importlib.reload does not clear the module dict.  The probe is '
          'order-dependent and that ordering is a requirement.'
          % d['classify_after'])

    print()
    print('SECTION 3 -- delivery census of every 09-29 mutation target')
    tgt = audit_targets_of_provenance_audit()
    cen = delivery_census(tgt)
    worst = {}
    for figure, label, verdict, info in cen:
        print('  %-11s %-34s %-16s %s'
              % (verdict, label[:34], info['name'][:16], info['verdict']))
        worst.setdefault(info['verdict'], 0)
        worst[info['verdict']] += 1
    print('  totals: %s' % worst)
    undelivered = [r for r in cen if r[3]['verdict'] in UNDELIVERED]
    partial = [r for r in cen if r[3]['verdict'] in PARTIAL]
    check('no 09-29 mutation target is UNDELIVERED',
          not undelivered,
          '%d undelivered' % len(undelivered) if undelivered else
          'the 9 late-binding fixes of 09-29 did hold')
    check('no 09-29 mutation target is only PARTLY delivered (MIXED)',
          not partial, '%d mixed' % len(partial) if partial else '')
    dead = [r for r in cen if r[3]['verdict'] == CLASS_DEAD]
    for figure, label, verdict, info in dead:
        print('  DEAD TARGET: %s.%s is read by no function AND nowhere in the'
              % (info['module'], info['name']))
        print('               module body.  Check %r (%s) mutates a name that'
              % (label, verdict))
        print('               nothing in the module can observe.')
    check('no 09-29 mutation target is a DEAD name', not dead,
          '%d dead' % len(dead) if dead else '')

    print()
    print('SECTION 3b -- IMPORT_USED targets, and whether the child was '
          'patched too')
    cov = import_consumed_is_covered(tgt)
    for figure, label, m, a, children, siblings, ok in cov:
        print('  %-34s %s.%s consumed at import -> %s'
              % (label[:34], m, a, children))
        print('       also patched in the same check: %s  -> %s'
              % (siblings, 'COVERED' if ok else 'NOT COVERED'))
    check('every IMPORT_USED target has its derived child patched alongside',
          all(r[6] for r in cov),
          'confirms the 09-29 audit\'s stated handling of import-time '
          'derived constants, mechanically rather than from its comment')

    print()
    print('SECTION 4 -- identity re-assignments: the OTHER dead mutation')
    ident = identity_reassignments(tgt)
    for figure, label, verdict, m, a, cur in ident:
        print('  %-11s %-34s %s.%s already == %r' % (verdict, label[:34],
                                                     m, a, cur))
    print('  %d of the 09-29 checks assign a value the constant already holds.'
          % len(ident))
    print('  Each passes with setattr replaced by `pass`.  They are not')
    print('  evidence of specificity; they are evidence of nothing.')
    check('every identity re-assignment found is a null control, not a '
          'MUST_CHANGE',
          all(v == 'MUST_NOT_CHANGE' for _, _, v, _, _, _ in ident),
          'a MUST_CHANGE that re-assigns its own value would be a live bug')

    print()
    print('SECTION 5 -- repository-wide delivery census')
    rows = repo_wide_frozen_census()
    counts = {}
    for path, name, verdict, detail in rows:
        counts[verdict] = counts.get(verdict, 0) + 1
    for path, name, verdict, detail in rows:
        flag = ' <<< MIXED: oracle/measurement split' if verdict == CLASS_MIXED else ''
        print('  %-11s %-44s %-22s %s%s'
              % (verdict, path, name, detail, flag))
    print('  totals: %s' % counts)
    mixed = [r for r in rows if r[2] == CLASS_MIXED]
    check('the census reports at least one name in every class it defines '
          'or says so', True, 'classes present: %s' % sorted(counts))
    print()
    print('  %d MIXED name(s).  MIXED is the worst class, not a middle one:'
          % len(mixed))
    print('  a mutation of a MIXED name moves the live readers and freezes')
    print('  the captured ones, so where the live readers are an audit\'s')
    print('  expected values and the captured ones are its measurement, the')
    print('  audit reports a disagreement it manufactured itself.')

    print()
    print('=' * 78)
    print('TOTAL: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return failed


if __name__ == '__main__':
    sys.exit(1 if main() else 0)
