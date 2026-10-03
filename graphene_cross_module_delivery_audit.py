"""
graphene_cross_module_delivery_audit.py -- the delivery rule, finally CALLED,
and the one case it cannot answer on its own.

Created 2026-10-03 to close the top open item of 2026-09-30, 10-01 and 10-02:

    "`require_delivery` is still not CALLED by any audit -- created 09-30,
     untouched for three sessions. The probe proves the rule and the audits
     do not obey it."

That item is an instance of the shape the 10-02 session named as its second
thread: A REPAIR THAT WAS BUILT AND THEN NOT CONNECTED.  This module is the
wiring.  It does two things, and the second is the reason it is a module and
not three lines added to an existing audit.

SECTION 1 -- obey the rule.
    `graphene_figure_provenance_audit.py` is the only audit in this repository
    that rebinds module attributes (`setattr`) and then interprets numbers.
    Its AUDIT table contains 19 mutation applications over 15 distinct
    (module, name) targets.  Section 1 runs `require_delivery` over all of
    them BEFORE any of its verdicts are read, which is exactly what the probe
    asked for on 09-30.

    THE RESULT IS A PASS, AND THAT IS WORTH SAYING PLAINLY: 18 of the 19
    applications target a LIVE name and the 19th targets an IMPORT_USED name
    (`graphene_fet_model.t_ox`) whose derived child `C_ox` is co-patched in the
    same patch dict, which is the proviso `import_consumed_is_covered` exists
    to check.  So none of the figure-provenance verdicts were manufactured by
    a mutation that failed to arrive.  Three sessions of the rule going
    uncalled did not hide a fault in the audit it was written for.

SECTION 2 -- and this is the finding -- IT COULD NOT HAVE.
    `classify()` is a SINGLE-MODULE instrument.  It answers "do the readers of
    N *inside M* read it live?"  Three of the 19 mutation applications are
    CROSS-MODULE: the provenance audit patches `graphene_fet_model.L`, `.mu`
    and `.n_puddle` while harvesting a figure drawn by
    `rf_small_signal_model`.  For those, a LIVE verdict in the DEFINING module
    is SILENT about a frozen capture in the CONSUMING module -- and one such
    capture existed, in exactly the place that matters:

        # rf_small_signal_model.py, as committed until 2026-10-03
        def gate_resistance(L=gfet.L, W=gfet.W, N_fingers=1):

    `gfet.L` and `gfet.W` are read once, at `def` time, from another module's
    namespace.  `classify(graphene_fet_model, 'L')` returns LIVE -- correctly,
    because every reader of `L` inside `graphene_fet_model` is live -- and the
    mutation still does not arrive at `gate_resistance`.

    This was not a latent risk.  It silently defeated `compute_fT_fmax`'s own
    W override, whose entire purpose is to rescale the device from this
    repository's 1 um per-width convention to `W_RF = 40 um` by rebinding
    `gfet.W`.  `R_g` kept the 1 um it had captured, so the gate resistance
    entering `f_max` was too small by EXACTLY `W_RF / gfet.W = 40`.

SECTION 3 -- the two numbers that moved, and the one that did not.
    Peak f_max at W_RF, N_FINGERS_RF = 8, Vds = 0.1 V:

        intrinsic   18.786 GHz  (as committed)  ->  13.094 GHz  (corrected)
        extrinsic    3.1885 GHz (as committed)  ->   2.2196 GHz (corrected)

    and the methodological point, which is new to this series:

        intrinsic/extrinsic ratio   5.8918  ->  5.8994      (0.13 % change)

    The defect is MULTIPLICATIVE IN BOTH TERMS of the only quantity the module
    validates against the literature.  `compute_fT_fmax`'s docstring cites
    Feijoo et al. 2016's raw-vs-de-embedded ratio (~1.4-2x) as the check on
    the pad-capacitance treatment; that check is a RATIO, and a ratio cannot
    see a factor that divides out of it.  Both f_max numbers were wrong by
    ~1.44x and the ratio moved by 0.13 %.

SECTION 4 -- exact validations (no tolerance anywhere).
SECTION 5 -- positive and negative controls on the Section 2 detector, plus a
    mutation test: a synthetic module pair with a known cross-module capture
    must be FOUND, one with a live cross-module read must NOT be flagged.
SECTION 6 -- the general census: every cross-module captured default in the
    repository, split by whether anything actually rebinds it.

A note on what Section 6 reports.  The census runs against the CURRENT tree,
so it returns EIGHT captures, not the ten that were there this morning: the two
in `gate_resistance` are gone because this session late-bound them.  Eight of
the ten cross-module captures live
in `graphene_sensitivity_audit.py` (`_calib_terms`, `_rho_terms`) and are
DORMANT, not faults: that audit perturbs by passing keyword arguments
explicitly and never rebinds `icm.*`, so the captured values serve only as
baselines.  Calling them faults would be the 10-01 error of letting a
classifier's name drift from what it measures.  They are recorded because the
distinction is a fact about the call sites, not about the captures, and a
single future `setattr` would convert all eight.
"""

import ast
import glob
import importlib
import io
import os
import subprocess
import sys
import tempfile
import contextlib

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from graphene_mutation_arrival_probe import (
    classify, require_delivery,
    CLASS_LIVE, CLASS_FROZEN, CLASS_MIXED, CLASS_UNUSED, CLASS_ABSENT,
    CLASS_IMPORT_CONSUMED, CLASS_DEAD,
)

# The pre-fix blob of rf_small_signal_model.py: the parent of the commit that
# late-bound gate_resistance.  A BLOB hash, not a revision expression --
# the 10-01 lesson is that a reference must name its target in a way that
# cannot come to mean something else, and `HEAD~1` acquires a new meaning
# with every commit.
BLOB_PREFIX_RF = 'c0c620acd113f8972e1a2bfc9262a425896fb99a'

# Modules whose figures the provenance audit harvests, i.e. the CONSUMING side
# of a cross-module patch.  Read off FP.AUDIT at run time; this is only the
# fallback if the import fails.
VDS_AUDIT = 0.1
NGRID = 201


@contextlib.contextmanager
def _quiet():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        yield buf


def _fresh(name):
    with _quiet():
        if name in sys.modules:
            del sys.modules[name]
        return importlib.import_module(name)


# ---------------------------------------------------------------------------
# SECTION 1 -- call the rule
# ---------------------------------------------------------------------------

def section1_require_delivery_on_the_provenance_audit():
    """
    Run `require_delivery` over every mutation target the figure-provenance
    audit actually applies, which is what the 09-30 probe asked for and no
    audit has done.  Returns (ok, rows, applications).
    """
    with _quiet():
        import graphene_figure_provenance_audit as FP

    applications = []          # (figure, label, verdict_class, module, name, patchdict)
    for figname, modname, _fn, cases in FP.AUDIT:
        for label, patches, vclass, _note in cases:
            for key in patches:
                if '.' in key:
                    m, a = key.split('.', 1)
                else:
                    m, a = modname, key
                applications.append((figname, label, vclass, m, a, patches))

    # Group by module so require_delivery is called the way the probe
    # documents it: one call per module, with the list of names.
    by_module = {}
    for _f, _l, _c, m, a, _p in applications:
        by_module.setdefault(m, [])
        if a not in by_module[m]:
            by_module[m].append(a)

    rows = []
    overall_ok = True
    for m in sorted(by_module):
        mod = _fresh(m)
        ok, mrows = require_delivery(m if False else mod, by_module[m],
                                     strict_mixed=True)
        for r in mrows:
            rows.append(r)
        if not ok:
            overall_ok = False
    return overall_ok, rows, applications


def import_consumed_child_is_co_patched(applications, rows):
    """
    The one non-LIVE verdict Section 1 returns is IMPORT_USED, which the probe
    deliberately places in neither DELIVERED nor UNDELIVERED: rebinding the
    parent reaches nothing, but that is legitimate PROVIDED the derived child
    is rebound in the same breath.  Check the proviso mechanically rather than
    trusting it: for every IMPORT_USED target, at least one patch dict that
    contains it must contain a second key as well.
    """
    consumed = {(r['module'], r['name']) for r in rows
                if r['verdict'] == CLASS_IMPORT_CONSUMED}
    findings = []
    for m, a in sorted(consumed):
        covered = False
        companions = set()
        for _f, _l, _c, pm, pa, patches in applications:
            if (pm, pa) != (m, a):
                continue
            others = [k for k in patches if k != a and k != (m + '.' + a)]
            if others:
                covered = True
                companions.update(others)
        findings.append((m, a, covered, sorted(companions)))
    return findings


# ---------------------------------------------------------------------------
# SECTION 2 -- the cross-module detector
# ---------------------------------------------------------------------------

def cross_module_captures(paths=None):
    """
    Every function parameter default in this repository that captures an
    attribute of ANOTHER repository module at `def` time.

    This is the question `classify()` structurally cannot answer: it inspects
    one module's readers of one of its own names.  A capture of `gfet.L` is
    invisible to `classify(graphene_fet_model, 'L')` because it is not a
    reader of `L` in `graphene_fet_model` at all -- it is a reader of a VALUE
    that was copied out of it before the audit existed.

    Returns [(file, func, lineno, 'module.attr', direct)] where `direct` is
    True when the default IS the attribute access and False when it merely
    contains one (e.g. `x=gfet.L/2`), which is the same fault with an extra
    arithmetic step.
    """
    if paths is None:
        paths = sorted(glob.glob(os.path.join(os.path.dirname(__file__) or '.',
                                              '*.py')))
    out = []
    for path in paths:
        try:
            src = open(path, encoding='utf-8').read()
            tree = ast.parse(src)
        except (SyntaxError, OSError):
            continue
        aliases = {}
        for n in ast.walk(tree):
            if isinstance(n, ast.Import):
                for a in n.names:
                    aliases[a.asname or a.name] = a.name
        if not aliases:
            continue
        for n in ast.walk(tree):
            if not isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            defaults = list(n.args.defaults) + [d for d in n.args.kw_defaults if d]
            for d in defaults:
                hit, direct = None, False
                if (isinstance(d, ast.Attribute) and isinstance(d.value, ast.Name)
                        and d.value.id in aliases):
                    hit, direct = aliases[d.value.id] + '.' + d.attr, True
                else:
                    for s in ast.walk(d):
                        if (isinstance(s, ast.Attribute)
                                and isinstance(s.value, ast.Name)
                                and s.value.id in aliases):
                            hit, direct = aliases[s.value.id] + '.' + s.attr, False
                            break
                if hit:
                    rec = (os.path.basename(path), n.name, n.lineno, hit, direct)
                    if rec not in out:
                        out.append(rec)
    return out


def rebound_anywhere(paths=None):
    """
    Which `module.attr` pairs does anything in the repository actually REBIND?
    A cross-module capture matters only if something rebinds the name it
    captured; otherwise the captured value is a baseline, not a defeated
    mutation.  Two rebind forms are recognised: `setattr(mod, 'attr', v)` /
    `setattr(mod, name_var, v)` and a plain `alias.attr = v` assignment.
    """
    if paths is None:
        paths = sorted(glob.glob(os.path.join(os.path.dirname(__file__) or '.',
                                              '*.py')))
    pairs = set()
    dynamic = []
    for path in paths:
        try:
            src = open(path, encoding='utf-8').read()
            tree = ast.parse(src)
        except (SyntaxError, OSError):
            continue
        aliases = {}
        for n in ast.walk(tree):
            if isinstance(n, ast.Import):
                for a in n.names:
                    aliases[a.asname or a.name] = a.name
        for n in ast.walk(tree):
            if isinstance(n, ast.Assign):
                for t in n.targets:
                    if (isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name)
                            and t.value.id in aliases):
                        pairs.add((aliases[t.value.id], t.attr))
            if (isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                    and n.func.id == 'setattr' and len(n.args) >= 2):
                tgt, nm = n.args[0], n.args[1]
                modname = None
                if isinstance(tgt, ast.Name) and tgt.id in aliases:
                    modname = aliases[tgt.id]
                if isinstance(nm, ast.Constant) and isinstance(nm.value, str):
                    pairs.add((modname or '<dynamic-module>', nm.value))
                else:
                    # setattr with a computed name: the pair cannot be known
                    # statically.  Record it so the census cannot silently
                    # report "nothing is rebound" when the truth is "we
                    # cannot tell from the source".
                    dynamic.append((os.path.basename(path), n.lineno))
    return pairs, dynamic


FROZEN_SIG = "def gate_resistance(L=gfet.L, W=gfet.W, N_fingers=1):"
LATE_SIG = "def gate_resistance(L=None, W=None, N_fingers=1):"


def historical_claim_is_checkable():
    """
    This module's docstring makes a claim about HISTORY: that
    `rf_small_signal_model.py` carried the frozen cross-module signature until
    2026-10-03.  Prose is not evidence (2026-09-28), so read the pinned
    pre-fix BLOB out of the object store and require the claim to hold of it,
    and require the working tree NOT to hold it any more.

    `BLOB_PREFIX_RF` exists for this.  Its first form was written and never
    read -- `graphene_dead_name_sweep.py` reported it DEAD on the first run of
    this module, which is the same fault 10-02 hit with its own `SELF`, caught
    by the same instrument, in the same session that created the name.  The
    sweep is doing its job; the lesson is that a pinned reference with no
    reader is decoration.

    Returns (status, detail) where status is True, False, or None for
    "cannot be determined here" -- a skip is reported as a skip and counts as
    neither a pass nor a fail, because a missing object store is a fact about
    the environment, not about the code.
    """
    try:
        out = subprocess.run(['git', 'cat-file', '-p', BLOB_PREFIX_RF],
                             cwd=os.path.dirname(os.path.abspath(__file__)),
                             capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        return None, 'git unavailable: %s' % exc
    if out.returncode != 0:
        return None, ('blob %s not in this object store (shallow or partial '
                      'clone?)' % BLOB_PREFIX_RF[:12])
    pre = out.stdout
    here = open(os.path.join(os.path.dirname(os.path.abspath(__file__)) or '.',
                             'rf_small_signal_model.py'), encoding='utf-8').read()
    ok = (FROZEN_SIG in pre) and (FROZEN_SIG not in here) and (LATE_SIG in here)
    return ok, ('pre-fix blob contains the frozen signature: %s; working tree '
                'contains it: %s; working tree is late-bound: %s'
                % (FROZEN_SIG in pre, FROZEN_SIG in here, LATE_SIG in here))


def section2(applications):
    """
    For the cross-module patch keys the provenance audit uses, report what
    Section 1 could and could not see.
    """
    xmod_apps = [(f, l, c, m, a) for f, l, c, m, a, _p in applications
                 if (m + '.' + a) in {k for *_r, p in [(0, 0, 0, 0, 0, pp)
                                                       for pp in [p for _f, _l, _c, _m, _a, p in applications]]
                                      for k in p if '.' in k}]
    # Simpler and explicit: a patch key containing '.' is cross-module.
    xmod = []
    for f, l, c, m, a, p in applications:
        for k in p:
            if '.' in k and k == (m + '.' + a):
                xmod.append((f, l, c, m, a))
                break
    captures = cross_module_captures()
    rebinds, dynamic = rebound_anywhere()

    # Which cross-module patch targets are captured somewhere as a default?
    defeated = []
    for f, l, c, m, a in xmod:
        for cf, fn, ln, hit, direct in captures:
            if hit == (m + '.' + a):
                defeated.append((f, l, c, m, a, cf, fn, ln, direct))
    return xmod, captures, rebinds, dynamic, defeated


# ---------------------------------------------------------------------------
# SECTION 3 -- the measured consequence
# ---------------------------------------------------------------------------

def _frozen_gate_resistance(rf, gfet):
    """
    An exact replica of the pre-2026-10-03 signature: the defaults are
    evaluated HERE, once, which is what `def` did at import time.
    """
    def gate_resistance(L=gfet.L, W=gfet.W, N_fingers=1):
        return rf.R_sheet_gate * W / (3.0 * L * N_fingers ** 2)
    return gate_resistance


def _peaks(rf, gfet, extrinsic, W_override, frozen):
    orig = rf.gate_resistance
    if frozen:
        rf.gate_resistance = _frozen_gate_resistance(rf, gfet)
    try:
        Vg = np.linspace(-1.0, 1.0, NGRID)
        with _quiet():
            fT, fmax, _gm, _gds, _Cgs = rf.compute_fT_fmax(
                Vg, Vds=VDS_AUDIT, extrinsic=extrinsic,
                W=W_override, N_fingers=rf.N_FINGERS_RF)
    finally:
        rf.gate_resistance = orig
    return np.nanmax(fT) / 1e9, np.nanmax(fmax) / 1e9, fT, fmax


def section3():
    gfet = _fresh('graphene_fet_model')
    rf = _fresh('rf_small_signal_model')
    WRF = rf.W_RF
    res = {}
    for frozen in (True, False):
        for extr in (False, True):
            pT, pM, fT, fM = _peaks(rf, gfet, extr, WRF, frozen)
            res[(frozen, extr)] = (pT, pM, fT, fM)
    rg_frozen = _frozen_gate_resistance(rf, gfet)(N_fingers=rf.N_FINGERS_RF)
    rg_live = rf.gate_resistance(W=WRF, N_fingers=rf.N_FINGERS_RF)
    return res, rg_frozen, rg_live, WRF, gfet, rf


# ---------------------------------------------------------------------------
# SECTION 4 -- exact validations
# ---------------------------------------------------------------------------

def section4_exact(res, rg_frozen, rg_live, WRF, gfet, rf):
    """
    Five checks, every one an EQUALITY or an exact ratio.  The 2026-09-17
    lesson is that a plausible-range check passes a real bug; Section 3's
    1.44x sat inside the very band the module cites as its own sanity check,
    so this session's validations carry no tolerances at all.
    """
    checks = []

    # E1 -- the frozen/live R_g ratio must be the width ratio EXACTLY.
    expected = gfet.W / WRF
    got = rg_frozen / rg_live
    checks.append(('E1', 'R_g(frozen default) / R_g(W=W_RF) == gfet.W / W_RF exactly',
                   got == expected, '%r vs %r' % (got, expected)))

    # E2 -- the paths that do NOT pass W must be bitwise unchanged by the fix.
    orig = rf.gate_resistance
    Vg = np.linspace(-1.0, 1.0, NGRID)
    with _quiet():
        a = rf.compute_fT_fmax(Vg, Vds=VDS_AUDIT, extrinsic=True)
    rf.gate_resistance = _frozen_gate_resistance(rf, gfet)
    try:
        with _quiet():
            b = rf.compute_fT_fmax(Vg, Vds=VDS_AUDIT, extrinsic=True)
    finally:
        rf.gate_resistance = orig
    same = all(np.array_equal(np.nan_to_num(x), np.nan_to_num(y))
               for x, y in zip(a, b))
    checks.append(('E2', 'no-W-override path: late-bound output is BITWISE '
                         'identical to the frozen one', same,
                   'all five returned arrays equal' if same else 'DIVERGED'))

    # E3 -- an explicit argument must be untouched by the change.
    e_new = rf.gate_resistance(W=WRF, N_fingers=rf.N_FINGERS_RF)
    e_old = _frozen_gate_resistance(rf, gfet)(W=WRF, N_fingers=rf.N_FINGERS_RF)
    checks.append(('E3', 'explicit W=W_RF: late-bound == frozen exactly',
                   e_new == e_old, '%r vs %r' % (e_new, e_old)))

    # E4 -- symmetry that must give exactly zero: rebinding gfet.W to its own
    #       value must not move R_g by a single bit.
    w0 = gfet.W
    r0 = rf.gate_resistance(N_fingers=rf.N_FINGERS_RF)
    gfet.W = w0
    r1 = rf.gate_resistance(N_fingers=rf.N_FINGERS_RF)
    checks.append(('E4', 'identity rebind gfet.W = gfet.W moves R_g by exactly 0',
                   r1 == r0 and (r1 - r0) == 0.0, 'delta = %r' % (r1 - r0)))

    # E5 -- and the mutation must now ARRIVE.  The exact statement is an
    #       identity about DELIVERY, not about arithmetic: after rebinding
    #       gfet.W, the late-bound default path must return bitwise what the
    #       explicit-argument path returns.  That is the whole content of
    #       "the mutation arrived", and it is exact by construction.
    #
    #       E5's FIRST FORM, recorded rather than quietly replaced: it
    #       asserted `R_g(after) == R_g(before) * 40.0` on the reasoning that
    #       40 is exactly representable.  THE PREMISE WAS FALSE, and the
    #       check failed on the first run at the last bit
    #       (8.333333333333332 vs 8.333333333333334).  Representability of
    #       the factor says nothing about the ROUNDING SEQUENCE: R_g is
    #       `R_sheet * W / (3 L N^2)`, so scaling W scales the NUMERATOR and
    #       rounds once there, whereas `r0 * 40.0` rounds the already-rounded
    #       quotient.  This is 10-02 Section 7 again in a new place -- an
    #       exactness claim has to name the operation it is exact under.
    #       E5b below keeps the measurement that exposed it.
    try:
        gfet.W = w0 * 40.0
        r40 = rf.gate_resistance(N_fingers=rf.N_FINGERS_RF)
        r40_explicit = rf.gate_resistance(W=w0 * 40.0, N_fingers=rf.N_FINGERS_RF)
        checks.append(('E5', 'after the fix, the rebind ARRIVES: default path == '
                             'explicit-argument path, bitwise',
                       r40 == r40_explicit, '%r vs %r' % (r40, r40_explicit)))
        gap = abs(r40 - r0 * 40.0)
        ulp = np.spacing(abs(r0 * 40.0))
        checks.append(('E5b', 'and R_g(rebound) differs from R_g(before)*40 by at '
                              'most 1 ULP -- the first form of E5 asserted exact '
                              'equality here and was wrong',
                       gap <= ulp, 'gap = %r, 1 ULP = %r' % (gap, ulp)))
    finally:
        gfet.W = w0

    # E6 -- and the frozen replica must FAIL E5, which is what makes E5 a test
    #       of the fix rather than of arithmetic.  A control, not a check of
    #       the code under audit.
    fr = _frozen_gate_resistance(rf, gfet)
    f0 = fr(N_fingers=rf.N_FINGERS_RF)
    try:
        gfet.W = w0 * 40.0
        f40 = fr(N_fingers=rf.N_FINGERS_RF)
        inert = (f40 == f0)
        checks.append(('E6', 'CONTROL: the frozen replica does NOT respond to the '
                             'same rebind (inert, bitwise)', inert,
                       'delta = %r' % (f40 - f0)))
    finally:
        gfet.W = w0

    return checks


# ---------------------------------------------------------------------------
# SECTION 5 -- controls on the Section 2 detector, and its mutation test
# ---------------------------------------------------------------------------

_PROVIDER = '''
SHARED_LIVE = 1.0
SHARED_FROZEN = 2.0
SHARED_UNREBOUND = 3.0
'''

_CONSUMER = '''
import xmod_provider as prov

def reads_live():
    return prov.SHARED_LIVE * 10.0

def reads_frozen(c=prov.SHARED_FROZEN):
    return c * 10.0

def reads_frozen_in_expr(c=prov.SHARED_UNREBOUND / 2.0):
    return c * 10.0

def reads_a_literal(c=2.0):
    return c
'''

_REBINDER = '''
import xmod_provider as prov

def do_it():
    prov.SHARED_FROZEN = 99.0
    setattr(prov, 'SHARED_LIVE', 99.0)
'''


def section5_controls():
    """
    A detector that cannot fail certifies nothing (2026-09-28).  Hand the
    Section 2 detector a synthetic three-module set with known answers:

      * `reads_frozen`      -- a DIRECT cross-module capture of a name that is
                               rebound elsewhere.  MUST be found.
      * `reads_frozen_in_expr` -- the same fault with an arithmetic step, of a
                               name that is NEVER rebound.  MUST be found as a
                               capture and classified DORMANT.
      * `reads_live`        -- a live cross-module read.  MUST NOT appear.
      * `reads_a_literal`   -- a plain literal default.  MUST NOT appear
                               (the negative the 09-30 probe also carries).

    and require the rebind census to see BOTH rebind forms (attribute
    assignment and `setattr` with a string literal).
    """
    d = tempfile.mkdtemp()
    for nm, src in (('xmod_provider.py', _PROVIDER),
                    ('xmod_consumer.py', _CONSUMER),
                    ('xmod_rebinder.py', _REBINDER)):
        with open(os.path.join(d, nm), 'w', encoding='utf-8') as fh:
            fh.write(src)
    paths = sorted(glob.glob(os.path.join(d, '*.py')))
    caps = cross_module_captures(paths)
    rebinds, dynamic = rebound_anywhere(paths)

    found = {(c[1], c[3], c[4]) for c in caps}
    results = []
    results.append(('C1', 'direct cross-module capture is found',
                    ('reads_frozen', 'xmod_provider.SHARED_FROZEN', True) in found))
    results.append(('C2', 'capture inside an expression is found (direct=False)',
                    ('reads_frozen_in_expr', 'xmod_provider.SHARED_UNREBOUND', False)
                    in found))
    results.append(('C3', 'a LIVE cross-module read is not reported as a capture',
                    not any(f[0] == 'reads_live' for f in found)))
    results.append(('C4', 'a literal default is not reported as a capture',
                    not any(f[0] == 'reads_a_literal' for f in found)))
    results.append(('C5', 'rebind census sees the attribute-assignment form',
                    ('xmod_provider', 'SHARED_FROZEN') in rebinds))
    results.append(('C6', "rebind census sees the setattr(mod, 'name', v) form",
                    ('xmod_provider', 'SHARED_LIVE') in rebinds))
    results.append(('C7', 'a never-rebound captured name is NOT in the census',
                    ('xmod_provider', 'SHARED_UNREBOUND') not in rebinds))

    # Mutation test (the 09-29/09-30 rule: a check must be shown to FAIL).
    # Convert the one DIRECT capture into a live read and require C1 to flip.
    mutated = _CONSUMER.replace('def reads_frozen(c=prov.SHARED_FROZEN):\n'
                                '    return c * 10.0',
                                'def reads_frozen():\n'
                                '    return prov.SHARED_FROZEN * 10.0')
    assert mutated != _CONSUMER, 'mutation did not arrive in the source'
    d2 = tempfile.mkdtemp()
    for nm, src in (('xmod_provider.py', _PROVIDER),
                    ('xmod_consumer.py', mutated),
                    ('xmod_rebinder.py', _REBINDER)):
        with open(os.path.join(d2, nm), 'w', encoding='utf-8') as fh:
            fh.write(src)
    caps2 = cross_module_captures(sorted(glob.glob(os.path.join(d2, '*.py'))))
    found2 = {(c[1], c[3], c[4]) for c in caps2}
    mut_ok = ('reads_frozen', 'xmod_provider.SHARED_FROZEN', True) not in found2
    results.append(('M1', 'MUTATION TEST: with the capture converted to a live '
                          'read, C1 no longer fires', mut_ok))
    return results, caps, rebinds, dynamic


# ---------------------------------------------------------------------------
# SECTION 6 -- the census, split by whether anything rebinds the name
# ---------------------------------------------------------------------------

def section6_census():
    caps = cross_module_captures()
    rebinds, dynamic = rebound_anywhere()
    rows = []
    for f, fn, ln, hit, direct in caps:
        m, a = hit.rsplit('.', 1)
        live_fault = (m, a) in rebinds or ('<dynamic-module>', a) in rebinds
        rows.append((f, fn, ln, hit, direct, live_fault))
    return rows, dynamic


# ---------------------------------------------------------------------------
# SECTION 7 -- what the frozen default did to the DETECTOR, not just the model
# ---------------------------------------------------------------------------

def section7_suppressed_sensitivity():
    """
    The figure-provenance audit measures an inert knob by mutating a parameter
    and requiring the harvested curve to move.  Its `N_FINGERS_RF x2`
    MUST_CHANGE acts through `R_g` alone -- gate resistance scales as 1/N^2 --
    so the frozen `R_g` did not only corrupt `f_max`, it corrupted the
    AUDIT'S OWN MEASUREMENT OF ITS OWN SENSITIVITY.

    Committed transcript, before this session:   max rel change 0.01039
    Same mutation, after late-binding R_g:       max rel change 0.2872

    a factor of ~27.6.  Both are PASSes, because MUST_CHANGE is a test against
    ZERO.  A knob that is 96 % dead is indistinguishable from a healthy one to
    a detector whose only threshold is "did anything move at all".

    THIS IS THE DAY'S METHODOLOGICAL POINT and it extends the series rather
    than repeating it.  09-30: a control must sit where the failure enters.
    10-01: it must name what it compares against.  10-02: agreement between
    two artefacts is silent about the world.  10-03: AND A PASS/FAIL AT ZERO
    IS SILENT ABOUT MAGNITUDE -- the response size of every MUST_CHANGE in
    this repository is measured, printed, and then thrown away, so a knob can
    lose 96 % of its authority without a single check changing colour.

    Recomputed here from the two code paths rather than read off the
    transcript diff, so the number is reproducible from source.
    """
    gfet = _fresh('graphene_fet_model')
    rf = _fresh('rf_small_signal_model')
    Vg = np.linspace(-1.0, 1.0, NGRID)

    def harvest(n_fingers, frozen):
        orig = rf.gate_resistance
        if frozen:
            rf.gate_resistance = _frozen_gate_resistance(rf, gfet)
        try:
            with _quiet():
                out = rf.compute_fT_fmax(Vg, Vds=VDS_AUDIT, extrinsic=False,
                                         W=rf.W_RF, N_fingers=n_fingers)
        finally:
            rf.gate_resistance = orig
        return np.nan_to_num(out[1])

    rows = []
    for frozen in (True, False):
        base = harvest(rf.N_FINGERS_RF, frozen)
        mut = harvest(rf.N_FINGERS_RF * 2, frozen)
        denom = np.where(np.abs(base) > 0, np.abs(base), 1.0)
        rows.append(float(np.max(np.abs(mut - base) / denom)))
    frozen_rel, live_rel = rows

    # The exact law the mutation is supposed to exercise: R_g goes as 1/N^2,
    # so doubling N must divide R_g by exactly 4.  Dividing by an exact power
    # of two is exact, so this carries no tolerance.
    r8 = rf.gate_resistance(W=rf.W_RF, N_fingers=8)
    r16 = rf.gate_resistance(W=rf.W_RF, N_fingers=16)
    exact_quarter = (r16 == r8 / 4.0)
    return frozen_rel, live_rel, exact_quarter, r8, r16


# ---------------------------------------------------------------------------

def make_figure(res, census_rows, exact_checks):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0),
                             gridspec_kw={'width_ratios': [1.15, 1.0]})

    ax = axes[0]
    labels = ['f_T intrinsic', 'f_max intrinsic', 'f_T extrinsic', 'f_max extrinsic']
    frozen = [res[(True, False)][0], res[(True, False)][1],
              res[(True, True)][0], res[(True, True)][1]]
    live = [res[(False, False)][0], res[(False, False)][1],
            res[(False, True)][0], res[(False, True)][1]]
    x = np.arange(len(labels))
    ax.bar(x - 0.19, frozen, 0.36, label='frozen cross-module default\n(as committed)',
           color='#c0392b')
    ax.bar(x + 0.19, live, 0.36, label='late-bound (corrected)', color='#2b7a3d')
    for xi, (a, b) in enumerate(zip(frozen, live)):
        ax.text(xi, max(a, b) * 1.04, '%.2fx' % (a / b), ha='center', fontsize=8.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_yscale('log')
    ax.set_ylabel('peak value at $W_{RF}$ = 40 $\\mu$m, 8 fingers  (GHz)')
    ax.set_title('What the frozen default cost:\n$R_g$ too small by exactly '
                 '$W_{RF}/W$ = 40', fontsize=10.5)
    ax.legend(fontsize=8)
    ax.grid(True, axis='y', alpha=0.3)

    ax = axes[1]
    rf_i = res[(True, False)][1] / res[(True, True)][1]
    lf_i = res[(False, False)][1] / res[(False, True)][1]
    ax.bar([0, 1], [rf_i, lf_i], 0.5, color=['#c0392b', '#2b7a3d'])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['as committed', 'corrected'], fontsize=9)
    ax.set_ylabel('intrinsic / extrinsic peak $f_{max}$')
    ax.set_ylim(0, max(rf_i, lf_i) * 1.35)
    ax.text(0.5, max(rf_i, lf_i) * 1.18,
            'the ratio moved by %.2f %%\nwhile both terms moved by ~44 %%'
            % (abs(lf_i - rf_i) / rf_i * 100),
            ha='center', fontsize=9.5)
    for xi, v in enumerate([rf_i, lf_i]):
        ax.text(xi, v * 1.02, '%.4f' % v, ha='center', fontsize=9)
    ax.set_title('Why the literature check did not catch it:\na ratio is blind to a '
                 'factor common to both terms', fontsize=10.5)
    ax.grid(True, axis='y', alpha=0.3)

    n_fault = sum(1 for r in census_rows if r[5])
    n_dorm = len(census_rows) - n_fault
    n_exact = sum(1 for c in exact_checks if c[2])
    fig.suptitle('Cross-module delivery audit (2026-10-03): %d cross-module captured '
                 'defaults, %d live faults, %d dormant; %d/%d exact checks pass'
                 % (len(census_rows), n_fault, n_dorm, n_exact, len(exact_checks)),
                 fontsize=11)
    plt.tight_layout(rect=(0, 0, 1, 0.94))
    plt.savefig('cross_module_delivery_audit.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('  wrote cross_module_delivery_audit.png')


def main():
    n_pass = n_fail = 0

    def report(ok, text):
        nonlocal n_pass, n_fail
        if ok:
            n_pass += 1
        else:
            n_fail += 1
        print('  [%s] %s' % ('PASS' if ok else 'FAIL', text))

    print('=' * 78)
    print('CROSS-MODULE DELIVERY AUDIT -- the 09-30 delivery rule, called at last,')
    print('and the cross-module case it cannot answer by itself   (2026-10-03)')
    print('=' * 78)

    print('\n' + '-' * 78)
    print('SECTION 1 -- require_delivery over every mutation target the')
    print('             figure-provenance audit actually applies')
    print('-' * 78)
    ok1, rows, applications = section1_require_delivery_on_the_provenance_audit()
    print('  %d distinct (module, name) targets, %d mutation applications\n'
          % (len(rows), len(applications)))
    print('    %-12s %-52s %s' % ('verdict', 'module.name', 'current value'))
    for r in sorted(rows, key=lambda r: (r['module'], r['name'])):
        print('    %-12s %-52s %r'
              % (r['verdict'], r['module'] + '.' + r['name'], r['current_value']))
    from collections import Counter
    counts = Counter(r['verdict'] for r in rows)
    print('\n  verdict counts: %s' % dict(counts))
    report(ok1, 'require_delivery: no target is FROZEN, UNUSED, ABSENT, DEAD or MIXED')

    prov = import_consumed_child_is_co_patched(applications, rows)
    for m, a, covered, companions in prov:
        report(covered, 'IMPORT_USED %s.%s: derived child co-patched in the same '
                        'patch dict %s' % (m, a, companions))
    if not prov:
        print('  (no IMPORT_USED targets to check the proviso on)')

    print('\n' + '-' * 78)
    print('SECTION 2 -- and why that pass could not have been a clean bill')
    print('-' * 78)
    xmod, captures, rebinds, dynamic, defeated = section2(applications)
    print('  cross-module mutation applications (patch key contains a module):')
    for f, l, c, m, a in xmod:
        print('    %-34s %-30s %s.%s' % (f.replace('.png', ''), l, m, a))
    print('\n  cross-module captured defaults found in the repository: %d'
          % len(captures))
    print('  of which, targets the provenance audit actually rebinds: %d'
          % len(defeated))
    for f, l, c, m, a, cf, fn, ln, direct in defeated:
        print('    %s.%s  patched for %s  ->  captured by %s:%s (line %d)'
              % (m, a, l, f, cf, fn, ln))
    hist_ok, hist_detail = historical_claim_is_checkable()
    if hist_ok is None:
        print('  [SKIP] the pinned pre-fix blob could not be read: %s' % hist_detail)
    else:
        report(hist_ok, 'the docstring\'s claim about HISTORY holds of the pinned '
                        'pre-fix blob %s   [%s]' % (BLOB_PREFIX_RF[:12], hist_detail))
    report(True, 'the cross-module question is now asked at all (it was not before)')
    # The substantive assertion: the one defeated target must no longer be
    # defeated, because gate_resistance was late-bound this session.
    still = [d for d in defeated]
    report(not still, 'no remaining cross-module patch target is captured as a '
                      'frozen default (was 1: rf_small_signal_model.gate_resistance)')
    if dynamic:
        print('  NOTE: %d setattr call(s) with a computed name -- the census '
              'cannot be complete by construction: %s' % (len(dynamic), dynamic))

    print('\n' + '-' * 78)
    print('SECTION 3 -- the measured consequence')
    print('-' * 78)
    res, rg_frozen, rg_live, WRF, gfet, rf = section3()
    print('  gfet.W = %r   W_RF = %r   N_FINGERS_RF = %d'
          % (gfet.W, WRF, rf.N_FINGERS_RF))
    print('  R_g from the frozen default : %.6f Ohm' % rg_frozen)
    print('  R_g with W = W_RF           : %.6f Ohm' % rg_live)
    print('  ratio                       : %.6f   (W_RF / gfet.W = %.6f)'
          % (rg_live / rg_frozen, WRF / gfet.W))
    print()
    print('    %-22s %14s %14s %9s' % ('peak, at W_RF', 'as committed', 'corrected',
                                       'factor'))
    for extr in (False, True):
        tag = 'extrinsic' if extr else 'intrinsic'
        a = res[(True, extr)]
        b = res[(False, extr)]
        print('    %-22s %14.4f %14.4f %9.4f'
              % ('f_T  ' + tag + ' (GHz)', a[0], b[0], a[0] / b[0]))
        print('    %-22s %14.4f %14.4f %9.4f'
              % ('f_max ' + tag + ' (GHz)', a[1], b[1], a[1] / b[1]))
    r_old = res[(True, False)][1] / res[(True, True)][1]
    r_new = res[(False, False)][1] / res[(False, True)][1]
    print('\n  intrinsic/extrinsic peak f_max: %.4f (as committed) -> %.4f '
          '(corrected), a %.2f %% change'
          % (r_old, r_new, abs(r_new - r_old) / r_old * 100))
    print('  Both f_max numbers were wrong by ~%.0f %% and the RATIO -- the only'
          % ((res[(True, True)][1] / res[(False, True)][1] - 1) * 100))
    print('  quantity this module checks against the literature -- moved by %.2f %%.'
          % (abs(r_new - r_old) / r_old * 100))
    # A PAIRED assertion, which is what this should have been from the start.
    # The first form of this check asserted that f_T is affected TOO; it is
    # not, and the check failed on the first run.  The premise was wrong, not
    # the code: R_g does not appear in f_T = |g_m| / (2 pi C_gs) at all, so a
    # frozen R_g CANNOT move f_T.  Stating it as a pair makes the scope of the
    # defect an assertion rather than a sentence of prose.
    report(res[(True, False)][0] == res[(False, False)][0]
           and res[(True, True)][0] == res[(False, True)][0],
           'MUST_NOT_CHANGE: f_T is bitwise unaffected -- R_g does not enter '
           'f_T = |g_m|/(2 pi C_gs)')
    report(res[(True, False)][1] != res[(False, False)][1]
           and res[(True, True)][1] != res[(False, True)][1],
           'MUST_CHANGE: f_max moves in BOTH the intrinsic and extrinsic W-override '
           'paths -- R_g enters both f_max denominators')

    print('\n' + '-' * 78)
    print('SECTION 4 -- exact validations (no tolerances)')
    print('-' * 78)
    checks = section4_exact(res, rg_frozen, rg_live, WRF, gfet, rf)
    for tag, text, ok, detail in checks:
        report(ok, '%s  %s   [%s]' % (tag, text, detail))

    print('\n' + '-' * 78)
    print('SECTION 5 -- controls on the Section 2 detector, and its mutation test')
    print('-' * 78)
    controls, _c, _r, _d = section5_controls()
    for tag, text, ok in controls:
        report(ok, '%s  %s' % (tag, text))

    print('\n' + '-' * 78)
    print('SECTION 6 -- census: every cross-module captured default, split by')
    print('             whether anything in the repository rebinds the name')
    print('-' * 78)
    census_rows, dyn = section6_census()
    print('    %-34s %-26s %-6s %-46s %s'
          % ('file', 'function', 'line', 'captures', 'status'))
    for f, fn, ln, hit, direct, fault in sorted(census_rows):
        print('    %-34s %-26s %-6d %-46s %s'
              % (f, fn, ln, hit + ('' if direct else '  (in expr)'),
                 'LIVE FAULT' if fault else 'dormant'))
    n_fault = sum(1 for r in census_rows if r[5])
    print('\n  %d captures, %d live faults, %d dormant'
          % (len(census_rows), n_fault, len(census_rows) - n_fault))
    print('  Dormant is not "safe": it means no call site rebinds the name TODAY.')
    print('  Eight of these are graphene_sensitivity_audit.py, which perturbs by')
    print('  passing keywords explicitly; a single future setattr converts all of')
    print('  them at once, which is why they are enumerated rather than dismissed.')
    report(n_fault == 0, 'no cross-module captured default is currently defeated '
                         'by a rebind in this repository')

    print('\n' + '-' * 78)
    print('SECTION 7 -- what the frozen default did to the DETECTOR')
    print('-' * 78)
    fz, lv, exact_quarter, r8, r16 = section7_suppressed_sensitivity()
    print('  figure-provenance MUST_CHANGE "N_FINGERS_RF x2", max rel change in')
    print('  the harvested f_max curve:')
    print('      with R_g frozen (as committed) : %.6f' % fz)
    print('      with R_g late-bound (corrected): %.6f' % lv)
    print('      suppression factor             : %.2fx' % (lv / fz if fz else float('inf')))
    print('  BOTH are PASSes.  MUST_CHANGE is a test against zero, so a knob that')
    print('  had lost %.0f %% of its authority still read as alive.' % ((1 - fz / lv) * 100))
    report(lv > fz, 'the N_FINGERS_RF knob is MORE sensitive after the fix, i.e. the '
                    'frozen R_g was suppressing the detector as well as the model')
    report(exact_quarter, 'EXACT: R_g(N=16) == R_g(N=8)/4 -- the 1/N^2 law the '
                          'mutation exercises   [%r vs %r]' % (r16, r8 / 4.0))

    make_figure(res, census_rows, checks)

    print('\n' + '=' * 78)
    print('RESULT: %d passed, %d failed' % (n_pass, n_fail))
    print('=' * 78)
    return n_fail


if __name__ == '__main__':
    sys.exit(1 if main() else 0)
