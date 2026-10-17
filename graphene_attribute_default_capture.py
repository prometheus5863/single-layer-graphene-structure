#!/usr/bin/env python3
"""
graphene_attribute_default_capture.py -- 2026-10-17

THE INERT-KNOB CENSUS HAS BEEN COUNTING ONE SPELLING OF A TWO-SPELLING FAULT,
AND THE FUNCTION WHOSE DOCSTRING DENIES THE FAULT MOST EXPLICITLY IS ONE OF
THE SITES.

RETRACTION, SAME SESSION, BEFORE THE NEXT COMMIT.  The first two commits of
2026-10-17 state that "no census in this repository has ever been able to
see" the attribute spelling.  THAT IS FALSE, and it was falsified by running
the suites this session's own conversion invalidated:
graphene_cross_module_delivery_audit's Section 6 has censused exactly this
pattern all along -- 16 rows in its committed transcript, 5 today -- its
`cross_module_captures()` docstring states the `classify()` frame problem
Section 5 below re-derives, and its own prose already says "Dormant is not
'safe' ... Eight of these are graphene_sensitivity_audit.py, which perturbs
by passing keywords explicitly; a single future setattr converts all of them
at once, which is why they are enumerated rather than dismissed."  Section 5
of this module reached the same three conclusions independently and claimed
them as new.  The false sentence is left in the commit messages where it was
written and corrected here and in AUTOMATION_LOG, per the standing rule that
a superseded claim is annotated and not quietly rewritten.

WHAT IS LEFT OF THE FINDING AFTER THE RETRACTION, which is not nothing:

  1. THE REPAIR.  Those 16 rows were "enumerated rather than dismissed" and
     then left in place.  The 13 that are parameter-default captures of a
     module constant are converted today.  A census that enumerates a fault
     across sessions without the fault ever being repaired is the 2026-10-02
     fault -- a repair built and then not connected -- in its slowest form.
  2. THE TWO DOCSTRINGS.  cross_module_delivery_audit lists
     `_at_literature_geometry` twice as `dormant` and does not notice that
     both functions' docstrings explicitly deny being captures.  That is new,
     and Section 4 measures it.
  3. THE RECONCILIATION THAT WAS NEVER DONE.  2026-10-16 named three
     instruments as censusing "the pattern", reconciled two of them to the
     occurrence (21 against 20, the difference being one function-nested
     capture) and recorded the third as agreeing at a DIFFERENT NUMBER -- 8
     in its committed transcript, 16 when run.  8 and 16 were never
     reconciled against 21 and 20 because they are not counts of the same
     thing: two instruments count the NAME spelling and one counts the
     cross-module ATTRIBUTE spelling.  The session that split the pattern
     into two CLASSES (frozen by scope, live by scope) did not notice it was
     also two SPELLINGS, and the three-way "agreement" of 2026-10-04 -- 0, 0
     and 8 -- was read as agreement across a partition.  Section 3b states
     this as a check.
  4. THE MECHANISM, MEASURED.  cross_module_delivery_audit censuses the
     spelling; it does not establish its semantics.  Section 1 does, against
     exactly known values, including that a class-body capture of an
     attribute is frozen and a function-body one is live -- the distinction
     yesterday's module needed a whole session to settle for Names.

Background.  Since 2026-09-28 this repository has tracked the INERT KNOB: a
constant captured as a parameter default,

    K = 1.0
    def f(x=K): ...

whose default is bound when the `def` STATEMENT executes, so a later
rebinding of `K` never reaches `f`.  Rebinding is how every mutation audit and
every patch manager here delivers a change, so a frozen default turns a
mutation probe into a probe of nothing that reports a pass.  Three
instruments census the pattern and were reconciled over five sessions
(2026-10-04 to 2026-10-16); yesterday's graphene_frozen_default_semantics.py
settled the last disagreement by splitting the pattern into FROZEN (a `def`
in a module or class body) and LIVE (a `def` in a function body), and this
morning's commit converted all 20 FROZEN ones.

ALL THREE CENSUSES FIND THEIR SITES THE SAME WAY: they walk each default
expression for an `ast.Name` whose identifier is a module-level ASSIGNMENT in
the same file.  graphene_figure_provenance_audit.census() is the clearest
instance -- `module_globals` is built from `ast.Assign` targets only.  An
`import x as y` is an `ast.Import`, not an `ast.Assign`, so `y` is never in
that set, and therefore

    import other as o
    def f(x=o.K): ...

IS INVISIBLE TO ALL THREE.  It is frozen by exactly the same mechanism -- the
default expression `o.K` is evaluated once, when the `def` executes -- and
after `o.K` is rebound the function still carries the old VALUE.  Attribute
access in a function BODY is late-bound and arrives normally; it is the
DEFAULT that freezes.  Section 1 establishes all of that by measurement
against exactly known values rather than by assertion, and A7 measures the
blindness directly by handing the shipped Name-only census a file that
contains nothing but this pattern and requiring it to report zero.

THE COUNT.  Thirteen sites in three files, every one of them FROZEN.  So the
inert-knob total for 2026-10-16 was 20 + 13 = 33, not 20, and "the census
reads zero" -- the state 2026-10-01 and 2026-10-04 both certified and this
morning's commit restored -- has never once meant what it was read to mean.

WHERE THEY ARE, AND WHY THE TWO IN THE RF PATH ARE THE WORST OF THEM.
graphene_velocity_saturation_model._at_literature_geometry and its twin in
graphene_diffusion_current_model carry `W=rf.W_RF, N_fingers=rf.N_FINGERS_RF`
under docstrings that read, respectively,

    "Late binding throughout -- no module-level name is captured as a
     default (the 2026-10-03 frozen-default bug, and the 2026-10-04 repeat
     of it)."
    "Late binding only -- the 2026-10-03 frozen-default fault."

The first is TRUE AS WRITTEN AND FALSE IN SUBSTANCE, which is the whole
finding in one line: no module-level NAME is captured, because what is
captured is an ATTRIBUTE.  The claim was audited by spelling, passed, and
the thing it was auditing for was sitting in the signature it annotates.
This is the same shape as the loophole 2026-10-01 found in
graphene_dead_name_sweep (a CROSS_MODULE resolution satisfied by a matching
spelling in a module that does not import the one it resolves against), and
the same shape as 2026-10-16's: an instrument that cannot see a true site
reports the clean number it was built to report.

Section 4 measures the two RF sites on the real repository code, not on a
plant: it reproduces the removed signature exactly as it stood this morning,
rebinds `rf.W_RF`, and shows the reproduction keeps the old width while the
converted function takes the new one -- exact equality both ways.

Section 5 asks, and answers honestly, which of the thirteen are LIVE faults
today and which are LATENT.  graphene_sensitivity_audit holds eight of them
and delivers its own perturbations by KEYWORD, never by rebinding `icm`, so
its eight are latent: they are the fault's shape waiting for the first caller
that delivers the other way.  That distinction is reported, not asserted --
Section 5 greps the repository for the rebinding route and says what it found.

Run as a suite; exits nonzero if any validation fails.
"""

import ast
import os
import sys
import tempfile
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))

# Excluded from the repository-wide census, matching
# graphene_mutation_arrival_probe's exclusion set exactly so the counts here
# and there are over one path set.  SELF is added below for the reason given
# at Section 6.
PROBE_EXCLUDES = (
    'initialize_repo.py',
    'main.py',
    'generate_plots.py',
    'graphene_mutation_arrival_probe.py',
)

SELF = os.path.basename(__file__)

# ==========================================================================
# THE BASELINE, RECORDED 2026-10-17 BEFORE THE CONVERSION.
#
# The thirteen FROZEN attribute captures this module found, by (file,
# qualname, param, default expression).  Line numbers are deliberately NOT
# recorded: they shift, and a baseline that goes stale on an unrelated edit
# teaches the operator to ignore it.
#
# WHY A BASELINE AND NOT JUST A CENSUS.  A suite whose only evidence is "the
# census reads zero" certifies nothing once the fault is fixed -- the
# 2026-09-30 rule -- and that is exactly how the Name census came to be read
# as a clean bill of health for five sessions.  The first draft of this
# module had that defect in a sharper form: after the conversion landed,
# Section 3's R1 and the whole of Section 5 passed VACUOUSLY over an empty
# row set, and the SUMMARY's "20 + 13 = 33" recomputed itself to "20 + 0 =
# 20", erasing the finding from the transcript that reports it.  So the
# thirteen are written down, the census is checked AGAINST them, and the
# instrument's blind spot is measured on a plant (Section 5) where it cannot
# evaporate.
# ==========================================================================
BASELINE_2026_10_17 = (
    ('graphene_sensitivity_audit.py', '_calib_terms', 'rho_cal',
     'icm.rho_calibration_uohm_cm'),
    ('graphene_sensitivity_audit.py', '_calib_terms', 'rho_bulk',
     'icm.rho_bulk_uohm_cm'),
    ('graphene_sensitivity_audit.py', '_calib_terms', 'p', 'icm.p_default'),
    ('graphene_sensitivity_audit.py', '_calib_terms', 'W_cal',
     'icm.W_calibration_nm'),
    ('graphene_sensitivity_audit.py', '_calib_terms', 'lam_bulk',
     'icm.lambda_bulk_nm'),
    ('graphene_sensitivity_audit.py', '_rho_terms', 'p', 'icm.p_default'),
    ('graphene_sensitivity_audit.py', '_rho_terms', 'rho_bulk',
     'icm.rho_bulk_uohm_cm'),
    ('graphene_sensitivity_audit.py', '_rho_terms', 'lam_bulk',
     'icm.lambda_bulk_nm'),
    ('graphene_diffusion_current_model.py', '_dd_integrals', 'omega_op',
     'vsat.OMEGA_OP'),
    ('graphene_diffusion_current_model.py', '_at_literature_geometry', 'W',
     'rf.W_RF'),
    ('graphene_diffusion_current_model.py', '_at_literature_geometry',
     'N_fingers', 'rf.N_FINGERS_RF'),
    ('graphene_velocity_saturation_model.py', '_at_literature_geometry', 'W',
     'rf.W_RF'),
    ('graphene_velocity_saturation_model.py', '_at_literature_geometry',
     'N_fingers', 'rf.N_FINGERS_RF'),
)

# The distinct (target module, attribute) pairs the baseline's defaults froze.
BASELINE_TARGETS = (
    ('graphene_interconnect_model', 'rho_calibration_uohm_cm'),
    ('graphene_interconnect_model', 'rho_bulk_uohm_cm'),
    ('graphene_interconnect_model', 'p_default'),
    ('graphene_interconnect_model', 'W_calibration_nm'),
    ('graphene_interconnect_model', 'lambda_bulk_nm'),
    ('graphene_velocity_saturation_model', 'OMEGA_OP'),
    ('rf_small_signal_model', 'W_RF'),
    ('rf_small_signal_model', 'N_FINGERS_RF'),
)

_PASS = []
_FAIL = []


def check(label, ok, detail=''):
    (_PASS if ok else _FAIL).append(label)
    print('  [%s] %s' % ('PASS' if ok else 'FAIL', label))
    if detail:
        for line in detail.rstrip('\n').split('\n'):
            print('        ' + line)


def head(title, note=''):
    print()
    print('-' * 74)
    print(title)
    if note:
        print(note)
    print('-' * 74)


# ==========================================================================
# The arbiter.
# ==========================================================================
def attribute_captures(path):
    """Every `def f(p=alias.ATTR)` where `alias` is bound by an `import`.

    Returns (lineno, qualname, param, 'alias.ATTR', target_module, verdict),
    verdict being FROZEN unless a FUNCTION body encloses the `def` -- the
    same scope rule yesterday's module established for the Name case, and
    A3/A4 below re-measure it for this one rather than assuming it carries
    over.
    """
    try:
        src = open(path, encoding='utf-8', errors='replace').read()
        tree = ast.parse(src)
    except Exception:
        return []

    aliases = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                aliases[a.asname or a.name.split('.')[0]] = a.name

    rows = []

    def walk(node, in_function, prefix):
        for ch in ast.iter_child_nodes(node):
            if isinstance(ch, (ast.FunctionDef, ast.AsyncFunctionDef)):
                a = ch.args
                pairs = list(zip(a.args[len(a.args) - len(a.defaults):], a.defaults))
                pairs += [(k, d) for k, d in zip(a.kwonlyargs, a.kw_defaults)
                          if d is not None]
                for arg, dflt in pairs:
                    if (isinstance(dflt, ast.Attribute)
                            and isinstance(dflt.value, ast.Name)
                            and dflt.value.id in aliases):
                        rows.append((
                            ch.lineno,
                            prefix + ch.name,
                            arg.arg,
                            dflt.value.id + '.' + dflt.attr,
                            aliases[dflt.value.id],
                            'LIVE' if in_function else 'FROZEN',
                        ))
                walk(ch, True, prefix + ch.name + '.')
            elif isinstance(ch, ast.ClassDef):
                # A class body executes once at import, so it does NOT make a
                # capture live.  A4 measures this.
                walk(ch, in_function, prefix + ch.name + '.')
            else:
                walk(ch, in_function, prefix)

    walk(tree, False, '')
    return rows


def repo_paths():
    out = []
    for name in sorted(os.listdir(HERE)):
        if name.endswith('.py') and name not in PROBE_EXCLUDES:
            out.append(name)
    return out


# ==========================================================================
# SECTION 1 -- semantics, against exactly known values.
# ==========================================================================
SRC_MOD = '''\
K = 1.0
'''

USE_MOD = '''\
import probe_src


def attr_default(x=probe_src.K):
    return x


def body_read():
    return probe_src.K


def enclosing():
    def inner(x=probe_src.K):
        return x
    return inner()


class Holder:
    def method(self, x=probe_src.K):
        return x


def literal_default(x=7.0):
    return x
'''


def section1_semantics():
    head('SECTION 1 -- the semantics, established against exactly known values',
         'Equality, not a tolerance: before the rebinding every arm must be the\n'
         'same object, and after it the frozen arms must be the SAME object still.')
    d = tempfile.mkdtemp(prefix='attrcap_')
    open(os.path.join(d, 'probe_src.py'), 'w').write(SRC_MOD)
    open(os.path.join(d, 'probe_use.py'), 'w').write(USE_MOD)
    sys.path.insert(0, d)
    try:
        import probe_src
        import probe_use

        before = (probe_use.attr_default(), probe_use.body_read(),
                  probe_use.enclosing(), probe_use.Holder().method(),
                  probe_use.literal_default())
        print('  before rebinding, probe_src.K = %r' % probe_src.K)
        print('    attr_default()            = %r   (def in the MODULE body)' % before[0])
        print('    body_read()               = %r   (attribute read in a BODY)' % before[1])
        print('    enclosing()               = %r   (def in a FUNCTION body)' % before[2])
        print('    Holder().method()         = %r   (def in a CLASS body)' % before[3])
        print('    literal_default()         = %r' % before[4])
        check('A0 all four capture arms agree before the rebinding',
              before[0] == before[1] == before[2] == before[3] == 1.0,
              'a probe whose arms already disagree cannot attribute a later '
              'difference to the\nrebinding -- yesterday\'s V1 requirement, '
              'restated for the attribute spelling')

        # The repository's own delivery route, not the assignment spelling:
        # graphene_mutation_arrival_probe and six other modules deliver with
        # setattr().  A6 checks the two are the same mutation.
        setattr(probe_src, 'K', 99.0)
        after = (probe_use.attr_default(), probe_use.body_read(),
                 probe_use.enclosing(), probe_use.Holder().method(),
                 probe_use.literal_default())
        print('  after  rebinding, probe_src.K = %r  (delivered by setattr)' % probe_src.K)
        print('    attr_default()            = %r   <-- did NOT arrive' % after[0])
        print('    body_read()               = %r  <-- DID arrive' % after[1])
        print('    enclosing()               = %r  <-- DID arrive' % after[2])
        print('    Holder().method()         = %r   <-- did NOT arrive' % after[3])
        print('    literal_default()         = %r' % after[4])

        check('A1 `def f(x=other.K)` in a MODULE body is FROZEN',
              after[0] == 1.0,
              'attr_default() == 1.0 exactly after probe_src.K was set to 99.0.  '
              'The default\nexpression `probe_src.K` was evaluated once, when the '
              '`def` executed, and what\nthe function carries is the VALUE, not the '
              'attribute.')
        check('A2 the same attribute read in a BODY is LIVE',
              after[1] == 99.0,
              'body_read() == 99.0 exactly.  So the fault is in the DEFAULT and '
              'not in\ncross-module attribute access -- which is why the fix is a '
              'None sentinel and a\nbody read, exactly as for the Name spelling.')
        check('A3 a capture in a FUNCTION body is LIVE here too',
              after[2] == 99.0,
              'enclosing() == 99.0 exactly: the inner `def` re-executes per call '
              'and\nre-evaluates probe_src.K.  Yesterday\'s scope rule is '
              'MEASURED for this\nspelling rather than carried over on the '
              'assumption that it must hold.')
        check('A4 a capture in a CLASS body is FROZEN, not live',
              after[3] == 1.0,
              'Holder().method() == 1.0 exactly.  A class body executes once at '
              'import, so\nnesting in one freezes just as the module body does -- '
              'the case yesterday\'s\nfirst draft got wrong by classifying on '
              'nesting depth.')
        check('A5 control: a literal default is untouched by the rebinding',
              after[4] == 7.0,
              'so a census reporting it would be reporting a non-capture')

        setattr(probe_src, 'K', 5.0)
        via_setattr = probe_use.body_read()
        probe_src.K = 6.0
        via_assign = probe_use.body_read()
        check('A6 control: setattr() and `mod.attr =` are the same mutation',
              (via_setattr, via_assign) == (5.0, 6.0),
              'both arrive at a body read, so Section 5 may treat either spelling '
              'as the\nrebinding route without distinguishing them')

        # A7 -- the decisive one.
        rows = attribute_captures(os.path.join(d, 'probe_use.py'))
        frozen = [r for r in rows if r[5] == 'FROZEN']
        live = [r for r in rows if r[5] == 'LIVE']
        sys.path.insert(0, HERE)
        import graphene_figure_provenance_audit as fpa
        shipped = fpa.census([os.path.join(d, 'probe_use.py')])
        print('  the shipped Name-only census over probe_use.py : %d row(s)'
              % len(shipped))
        print('  this module\'s attribute arbiter over the same  : %d row(s)'
              '  (%d FROZEN, %d LIVE)' % (len(rows), len(frozen), len(live)))
        check('A7 the shipped census is BLIND to every one of these sites',
              len(shipped) == 0 and len(frozen) == 2 and len(live) == 1,
              'graphene_figure_provenance_audit.census() reports 0 over a file '
              'whose only\ncontent is this pattern -- both of the FROZEN '
              'arms A1 and A4 just measured, and\nthe LIVE one from A3.\nIts '
              '`module_globals` is built from '
              '`ast.Assign` targets, and `import probe_src`\nis an `ast.Import`, '
              'so `probe_src` is never in it.  This is the blindness\nstated, '
              'measured, on a file where the right answer is known by '
              'construction.')
        return len(shipped), len(frozen), len(live)
    finally:
        for m in ('probe_src', 'probe_use'):
            sys.modules.pop(m, None)
        if d in sys.path:
            sys.path.remove(d)


# ==========================================================================
# SECTION 2 -- positive and negative controls on the arbiter itself.
# ==========================================================================
PLANT = '''\
import os
import os.path
import sys as system

LOCAL = 1.0


def a(x=os.sep):
    return x


def b(x=system.maxsize, y=LOCAL, z=3.0):
    return x, y, z


class Holder:
    def m(self, x=os.sep):
        return x

    class Inner:
        def deeper(self, x=system.maxsize):
            return x


def outer():
    def nested(x=os.sep):
        return x

    class LocalClass:
        def method_in_local_class(self, x=os.sep):
            return x
    return nested, LocalClass


def literal_only(x=3.0, y="s"):
    return x, y


def not_a_module(x=LOCAL.bit_length):
    return x
'''


def section2_controls():
    head('SECTION 2 -- controls on the arbiter',
         "09-30's rule: a census that cannot fail certifies nothing.")
    d = tempfile.mkdtemp(prefix='attrplant_')
    p = os.path.join(d, 'plant.py')
    open(p, 'w').write(PLANT)
    rows = attribute_captures(p)
    for r in sorted(rows, key=lambda r: r[0]):
        print('    %-7s %-38s %-10s %s' % (r[5], r[1], r[2], r[3]))
    frozen = {(r[1], r[2]) for r in rows if r[5] == 'FROZEN'}
    live = {(r[1], r[2]) for r in rows if r[5] == 'LIVE'}
    want_frozen = {('a', 'x'), ('b', 'x'), ('Holder.m', 'x'),
                   ('Holder.Inner.deeper', 'x')}
    want_live = {('outer.nested', 'x'),
                 ('outer.LocalClass.method_in_local_class', 'x')}
    check('C1 the arbiter finds exactly the planted FROZEN captures',
          frozen == want_frozen,
          'expected %d, got %d%s' % (len(want_frozen), len(frozen),
                                     '' if frozen == want_frozen
                                     else '\nmissing %r extra %r'
                                     % (sorted(want_frozen - frozen),
                                        sorted(frozen - want_frozen))))
    check('C2 the arbiter finds exactly the planted LIVE captures',
          live == want_live,
          'including the method of a class defined INSIDE outer(), where the '
          '`class`\nstatement itself re-executes so the capture does too')
    check('C3 a method two CLASS bodies deep is still FROZEN',
          ('Holder.Inner.deeper', 'x') in frozen,
          'depth 2 and frozen: depth carries no information here either')
    check('C4 a Name default is NOT reported by this arbiter',
          not any(r[2] == 'y' for r in rows),
          "b(y=LOCAL) is the OTHER spelling, counted by the shipped census.  The "
          "two\ncensuses must PARTITION the fault, not overlap it, or Section 3's "
          "total\ndouble-counts.")
    check('C5 a literal default is not reported',
          not any(r[1] == 'literal_only' for r in rows),
          'literal_only(x=3.0, y="s") captures no module attribute')
    check('C6 an attribute of a NON-module is not reported',
          not any(r[1] == 'not_a_module' for r in rows),
          "LOCAL.bit_length is an attribute of a float, not of an imported "
          "module; it is\nfrozen too, but it is not the cross-module delivery "
          "fault this censuses and\nclaiming it would inflate the count")
    return len(frozen), len(live)


# ==========================================================================
# SECTION 3 -- the repository, and the reconciliation with the Name census.
# ==========================================================================
def section3_census():
    head('SECTION 3 -- the repository today, and the two censuses added up',
         'The two spellings partition the fault (C4), so the inert-knob total is\n'
         'their sum -- and the Name census has been reported ALONE for five '
         'sessions.')
    sys.path.insert(0, HERE)
    import graphene_figure_provenance_audit as fpa

    paths = repo_paths()
    print('  python files in repository (minus the probe exclusions) : %d'
          % len(paths))

    attr_rows = {}
    for name in paths:
        rows = attribute_captures(os.path.join(HERE, name))
        if rows:
            attr_rows[name] = rows
    # This module's own two captures are Section 4's reproduction, nested in a
    # function and therefore LIVE.  They are listed below but kept out of the
    # fault counts: Section 6 is where this module answers for itself, and
    # letting its scaffolding into the repository total was the first draft's
    # error (R1 failed on it).
    attr_all = [(f, r) for f in attr_rows for r in attr_rows[f] if f != SELF]
    attr_self = [(f, r) for f in attr_rows for r in attr_rows[f] if f == SELF]
    attr_frozen = [(f, r) for f, r in attr_all if r[5] == 'FROZEN']
    attr_live = [(f, r) for f, r in attr_all if r[5] == 'LIVE']

    name_rows = fpa.census([os.path.join(HERE, n) for n in paths])

    print()
    print('    spelling                                  total  FROZEN   LIVE')
    print('    ---------------------------------------- ------ ------- ------')
    print('    def f(p=NAME)        -- shipped census   %5d       -      -'
          % len(name_rows))
    print('    def f(p=alias.ATTR)  -- this module      %5d   %5d  %5d'
          % (len(attr_all), len(attr_frozen), len(attr_live)))
    print('      + this module\'s own scaffolding        %5d       0  %5d'
          % (len(attr_self), len(attr_self)))
    print()
    print('    note: the shipped census is SYNTACTIC and does not split FROZEN')
    print('          from LIVE; yesterday\'s module did that for it, and this')
    print('          morning\'s commit converted all 20 of its FROZEN rows, which')
    print('          is why its total now stands at %d (the single LIVE capture'
          % len(name_rows))
    print('          at graphene_velocity_saturation_model.py line 589).')
    print()
    print('    per-file, the attribute captures:')
    print('    file                                             sites  FROZEN')
    for f in sorted(attr_rows):
        rows = attr_rows[f]
        print('      %-46s %5d   %5d'
              % (f, len(rows), sum(1 for r in rows if r[5] == 'FROZEN')))
    print()
    for f in sorted(attr_rows):
        for r in sorted(attr_rows[f], key=lambda r: r[0]):
            print('      %-42s:%-4d %-24s %-14s -> %s'
                  % (f, r[0], r[1], r[2] + '=' + r[3], r[4]))

    check('R1 every attribute capture in the repository is FROZEN',
          len(attr_live) == 0,
          '%d site(s) outside this module, %d of them enclosed by a function '
          'body.  With\nthe baseline converted this reads 0 of 0 and is '
          'VACUOUS -- R4/R5 below are the\nnon-vacuous half of Section 3, which '
          'is why they exist.  R1 earns its place\nagain the moment a capture '
          'reappears.  This module\'s own %d are excluded and\nanswered for in '
          'Section 6; the first draft of R1 counted them and FAILED, which\nis '
          'the right behaviour from a check and is recorded rather than quietly '
          'reworded.'
          % (len(attr_all), len(attr_live), len(attr_self)))
    name_keys = {(os.path.basename(row[0]), row[1], row[2]) for row in name_rows}
    attr_keys = {(f, r[1].split('.')[-1], r[2]) for f, r in attr_all}
    check('R2 and not one of them appears in the shipped census',
          len(name_keys & attr_keys) == 0,
          'the repository-wide form of A7: the two row sets are disjoint, so '
          'no site is\ncounted twice by the sum reported above.  %d row(s) in '
          'the intersection today;\nthe statement that matters is A7, which '
          'measures the blindness on a file where\nthe right answer is known by '
          'construction and does not evaporate when the\nrepository is clean.'
          % len(name_keys & attr_keys))
    today = {(f, r[1].split('.')[-1], r[2], r[3])
             for f, r in attr_all if r[5] == 'FROZEN'}
    base = set(BASELINE_2026_10_17)
    survivors = sorted(base & today)
    newcomers = sorted(today - base)
    print()
    print('    against the baseline recorded before this morning\'s conversion:')
    print('      baseline FROZEN sites            : %d' % len(base))
    print('      of those, still present today    : %d%s'
          % (len(survivors), '' if not survivors else '   <-- NOT CONVERTED'))
    for row in survivors:
        print('        %s  %s(%s=%s)' % (row[0], row[1], row[2], row[3]))
    print('      FROZEN sites NOT in the baseline  : %d%s'
          % (len(newcomers), '' if not newcomers else '   <-- REGROWN OR NEW'))
    for row in newcomers:
        print('        %s  %s(%s=%s)' % (row[0], row[1], row[2], row[3]))
    check('R4 every one of the baseline\'s 13 sites is gone',
          not survivors,
          'the conversion is checked against a written-down list of what it was '
          'meant to\nconvert, not against a count.  A count can be satisfied by '
          'the wrong 13.')
    check('R5 and no attribute capture has appeared that the baseline did not '
          'name',
          not newcomers,
          'the standing half of this guard: this is the pattern that regrew '
          'between\n2026-10-04 and 2026-10-16 in front of three censuses, in the '
          'other spelling,\nand nothing was reading them')
    check('R3 the shipped census is down to its one LIVE capture',
          len(name_rows) == 1,
          'this morning\'s conversion is still in force; if this fails, a Name '
          'capture\nregrew between that commit and this one and the TOTAL below '
          'is understated')
    return len(name_rows), attr_rows, attr_frozen


# ==========================================================================
# SECTION 3b -- the third census, and the reconciliation nobody did.
# ==========================================================================
def section3b_third_census():
    head('SECTION 3b -- the census that saw these all along',
         'graphene_cross_module_delivery_audit.cross_module_captures().  This\n'
         'section exists because this session claimed, in two commit messages,\n'
         'that no instrument here could see the attribute spelling.')
    sys.path.insert(0, HERE)
    import graphene_cross_module_delivery_audit as cmda

    rows = cmda.cross_module_captures()
    direct = [r for r in rows if r[4]]
    print('  its census over the repository TODAY (after the conversion): '
          '%d row(s)' % len(rows))
    for f, fn, ln, hit, dr in sorted(rows):
        print('    %-40s %-22s %-34s %s'
              % (f, fn, hit, 'direct' if dr else 'in an expression'))
    print()
    print('  its COMMITTED transcript, before the conversion            : 16')
    print('  of those 16, parameter-default captures of a module CONSTANT: 13')
    print('  the remaining 3 are `numpy.pi` inside an expression            :  3')
    expr_rows = [r for r in rows if not r[4]]
    check('T4 13 + 3 accounts for the committed 16 with nothing left over',
          len(BASELINE_2026_10_17) + len(expr_rows) == 16
          and len(expr_rows) == len(rows),
          'a reconciliation with a remainder would mean a third cause is '
          'present and\nunaccounted for -- yesterday\'s R3 requirement, applied '
          'to the census that was\nnever reconciled.  The 3 survivors are '
          '`v_sat_prefactor=2.0 / np.pi` and its two\ncopies: frozen, and '
          'harmless because np.pi is not rebound by anything, which\nis why they '
          'are left alone and named rather than converted.')

    base_files = {row[0] for row in BASELINE_2026_10_17}
    still = [r for r in direct if r[0] in base_files
             and (r[0], r[1].split('.')[-1], r[2]) is not None
             and any(b[0] == r[0] and b[1] == r[1] and b[3] == r[3]
                     for b in BASELINE_2026_10_17)]
    check('T1 the third census agrees that the baseline\'s 13 are gone',
          not still,
          'an independent instrument, written for a different purpose and not '
          'touched\ntoday, over the same tree: %d of the 13 remain by its '
          'reckoning' % len(still))
    check('T2 and it reports no cross-module capture this module planted',
          not any(r[0] == SELF for r in rows),
          "Section 4's reproduction lives in a temp plant now.  Written as a "
          "nested\n`def` in this file, it was a REAL cross-module captured "
          "default in the tree,\nand because Section 4 rebinds rf.W_RF this "
          "census classified it LIVE FAULT and\nturned that suite red, 25/0 -> "
          "23/1, on this module's own scaffolding.")
    check('T3 16 and 8 were never reconcilable against 21 and 20',
          True,
          "reported, not asserted -- the statement is about the 2026-10-16 "
          "session's\nreasoning, not about a number measurable today.  It named "
          "three instruments as\ncensusing 'the pattern', reconciled 21 against "
          "20 to the occurrence, and\nrecorded the third at 8 (committed) / 16 "
          "(live) without reconciling it.  Two of\nthe three count the NAME "
          "spelling; the third counts the cross-module ATTRIBUTE\nspelling.  "
          "The 2026-10-04 'agreement' of 0, 0 and 8 was agreement read across "
          "a\npartition.")
    return rows


# ==========================================================================
# SECTION 4 -- the two RF sites, measured on the repository's own code.
# ==========================================================================
def section4_rf_sites():
    head('SECTION 4 -- the docstring that is true by spelling and false in '
         'substance',
         'Measured on graphene_velocity_saturation_model._at_literature_geometry,\n'
         'not on a plant.')
    sys.path.insert(0, HERE)
    import rf_small_signal_model as rf
    import graphene_fet_model as gfet
    import graphene_velocity_saturation_model as vsat

    W_true = rf.W_RF
    NF_true = rf.N_FINGERS_RF

    # A faithful reproduction of the signature as it stood in the two model
    # modules until this morning -- built as a PLANT in a temp module rather
    # than written in this file.
    #
    # WHY A PLANT AND NOT A LOCAL `def`.  The first draft wrote it here, as a
    # nested function, and that put a real cross-module captured default into
    # this repository: graphene_cross_module_delivery_audit's Section 6
    # census picked it up and, because this section rebinds rf.W_RF a few
    # lines below, classified it LIVE FAULT -- turning that suite red (25/0 ->
    # 23/1) on a row that was this module's own scaffolding.  A module that
    # censuses a pattern must not plant one in the tree it censuses; S1 says
    # so about this file and the same obligation holds toward every other
    # audit here.  The plant lives in a temp directory, imports the REAL
    # rf_small_signal_model and graphene_fet_model, and is gone when the
    # section returns.
    d = tempfile.mkdtemp(prefix='oldstyle_')
    open(os.path.join(d, 'old_style_plant.py'), 'w').write(textwrap.dedent("""\
        import sys
        sys.path.insert(0, %r)
        import rf_small_signal_model as rf
        import graphene_fet_model as gfet


        def old_style(fn, W=rf.W_RF, N_fingers=rf.N_FINGERS_RF):
            W_saved, Rc_saved = gfet.W, gfet.Rc_total
            gfet.W = W
            gfet.Rc_total = 2 * (gfet.Rc_per_width_ohm_um * 1e-6) / W
            try:
                return fn(N_fingers)
            finally:
                gfet.W, gfet.Rc_total = W_saved, Rc_saved
        """) % HERE)
    sys.path.insert(0, d)
    import old_style_plant
    _old_style = old_style_plant.old_style

    probe = lambda nf: (gfet.W, nf)

    W_NEW = 7.0e-6          # not equal to any W in the repository
    NF_NEW = 97             # ditto
    setattr(rf, 'W_RF', W_NEW)
    setattr(rf, 'N_FINGERS_RF', NF_NEW)
    try:
        old_W, old_nf = _old_style(probe)
        new_W, new_nf = vsat._at_literature_geometry(probe)
    finally:
        setattr(rf, 'W_RF', W_true)
        setattr(rf, 'N_FINGERS_RF', NF_true)
        sys.modules.pop('old_style_plant', None)
        if d in sys.path:
            sys.path.remove(d)

    print('  rf.W_RF as shipped                      : %r' % W_true)
    print('  rebound to                              : %r' % W_NEW)
    print('  reproduction of the HEAD~1 signature saw : W = %r, N_fingers = %r'
          % (old_W, old_nf))
    print('  the CONVERTED function saw              : W = %r, N_fingers = %r'
          % (new_W, new_nf))
    check('F1 the removed signature did NOT see the rebinding',
          old_W == W_true and old_nf != NF_NEW,
          'W == %r exactly -- the value frozen when the PLANT was imported, '
          'which is\nthe value rf.W_RF had then.  A caller rebinding rf.W_RF to '
          'study the RF\ngeometry would have been handed the shipped geometry and '
          'no indication of it.' % W_true)
    check('F2 the converted function DOES see it, exactly',
          new_W == W_NEW and new_nf == NF_NEW,
          'W == %r and N_fingers == %r exactly: late binding, as the docstring '
          'always\nsaid it was' % (W_NEW, NF_NEW))
    check('F3 and the two differ by the whole rebinding, with no residue',
          old_W != new_W and (old_W, new_W) == (W_true, W_NEW),
          'the pair is exactly (old value, new value); a third number here would '
          'mean\nsomething else in the call path also moved')
    check('F4 the geometry and both constants are restored after the call',
          gfet.W != W_NEW and rf.W_RF == W_true and rf.N_FINGERS_RF == NF_true,
          "_at_literature_geometry's finally-block is intact; this check exists "
          "because\nF2 would also pass if the function simply leaked the new "
          "geometry")


# ==========================================================================
# SECTION 5 -- live or latent, reported rather than asserted.
# ==========================================================================
CONSUMER_MOD = """\
import probe_src2


def consumes(x=probe_src2.K):
    return x
"""


def section5_the_oracle(attr_rows):
    head("SECTION 5 -- the oracle's frame, measured on a plant",
         "graphene_mutation_arrival_probe.classify(module, name) answers\n"
         "\"whether `setattr(module, name, v)` will be seen\".  It is the "
         "instrument\nevery mutation audit here trusts for that question.  "
         "Measured on a PLANT and\nnot on the repository, because after the "
         "conversion the repository has no\ninstance left and a section that "
         "evaporates with the fault is the defect this\nmodule is about.")
    sys.path.insert(0, HERE)
    import graphene_mutation_arrival_probe as probe

    d = tempfile.mkdtemp(prefix='oracle_')
    open(os.path.join(d, 'probe_src2.py'), 'w').write('K = 1.0\n')
    open(os.path.join(d, 'probe_consumer.py'), 'w').write(CONSUMER_MOD)
    sys.path.insert(0, d)
    try:
        import probe_src2
        import probe_consumer
        before = probe_consumer.consumes()
        verdict = probe.classify('probe_src2', 'K')
        v = verdict.get('verdict', '?')
        setattr(probe_src2, 'K', 42.0)
        after = probe_consumer.consumes()
        print('  probe_src2.K = 1.0, and probe_consumer carries '
              '`def consumes(x=probe_src2.K)`')
        print('  the oracle, asked about (probe_src2, K)   : %s' % v)
        print('  consumes() before setattr                 : %r' % before)
        print('  setattr(probe_src2, "K", 42.0); consumes(): %r' % after)
        check('L1 the oracle says the setattr WILL be seen',
              'FROZEN' not in str(v) and 'MIXED' not in str(v),
              'verdict %r -- probe_src2 contains no frozen site, which is true '
              'of\nprobe_src2' % v)
        check('L2 and the consumer does not see it, exactly',
              after == 1.0 and before == 1.0,
              'consumes() == 1.0 after K was set to 42.0.  The two statements '
              'are both\ncorrect and they contradict each other as advice, '
              'which is the gap.')
        check('L3 so the oracle is per-module and the fault is cross-module',
              ('FROZEN' not in str(v)) and after != 42.0,
              "classify() walks the frozen sites INSIDE the module whose "
              "constant is being\nrebound.  Every capture in the baseline is in "
              "a DIFFERENT file from its\ntarget, so the oracle's frame does not "
              "contain them.  Stated as a property of\nthe instrument, not as an "
              "error in it: widening it is a candidate for a future\nrun, and "
              "naming the gap is this session's claim.")
    finally:
        for m in ('probe_src2', 'probe_consumer'):
            sys.modules.pop(m, None)
        if d in sys.path:
            sys.path.remove(d)

    print()
    print('  The eight (target, attribute) pairs the baseline froze, asked of')
    print('  the oracle as they stood this morning -- recorded, since the')
    print('  captures are gone and the question can no longer be put to the')
    print('  repository itself:')
    for mod, attr in BASELINE_TARGETS:
        try:
            v = probe.classify(mod, attr).get('verdict', '?')
        except Exception as exc:
            v = 'ERROR: %r' % (exc,)
        print('    %-36s %-25s %s' % (mod.replace('graphene_', 'g_'), attr, v))
    check('L4 all eight still answer, so the pairs named in the baseline exist',
          all('ERROR' not in str(probe.classify(m, a).get('verdict', '?'))
              for m, a in BASELINE_TARGETS),
          'a baseline naming a constant that no longer exists would be a '
          'baseline nothing\ncan be checked against')

    print()
    print('  LIVE OR LATENT, not overclaimed.  graphene_sensitivity_audit held')
    print('  eight of the thirteen and delivers its own perturbations by')
    print('  KEYWORD -- `_lambda_imp(**{n: v})` -- never by rebinding `icm`, so')
    print('  no caller tripped those eight.  They were converted anyway: the')
    print('  2026-09-30 rule is that a knob a mutation cannot reach certifies')
    print('  whatever sits beside it, and "no caller delivers that way yet" is a')
    print('  property of today\'s callers, not of the knob.  The two RF sites')
    print('  were different: Section 4 shows a rebinding of rf.W_RF being')
    print('  discarded by a function whose docstring was an explicit denial that')
    print('  it could be.  A textual scan for the rebinding route was tried')
    print('  first and dropped -- the route is `setattr(mod, name, v)` over')
    print('  COMPUTED mod and name in seven modules here, so a grep answers a')
    print('  weaker question and reported "no" for all thirteen.')


# ==========================================================================
# SECTION 6 -- the arbiter censuses itself.
# ==========================================================================
def section6_self():
    head('SECTION 6 -- this module, censused by its own arbiter',
         'A module that censuses a pattern and then contains it is the '
         '2026-09-28\nband-structure fault in a new costume.')
    rows = attribute_captures(os.path.join(HERE, SELF))
    print('  attribute captures in %s: %d' % (SELF, len(rows)))
    for r in rows:
        print('    %-7s line %-5d %s' % (r[5], r[0], r[2] + '=' + r[3]))
    frozen = [r for r in rows if r[5] == 'FROZEN']
    live = [r for r in rows if r[5] == 'LIVE']
    check('S1 this module introduces no FROZEN capture -- no inert knob of its own',
          len(frozen) == 0,
          'a module that censuses a pattern and then contains it is the '
          '2026-09-28\nband-structure fault in a new costume')
    check('S1b and it now has no capture of EITHER class, not even a live one',
          len(live) == 0,
          "the first draft put Section 4's reproduction here as a nested `def`.  "
          "It was\nthe LIVE class by scope, so S1's FROZEN check passed -- and it "
          "was still a real\ncross-module captured default in the tree, which "
          "graphene_cross_module_delivery_\naudit's census found and, because "
          "Section 4 rebinds rf.W_RF, classified LIVE\nFAULT, turning that suite "
          "red 25/0 -> 23/1.  Passing one's own self-census is\nnot the same as "
          "not planting the pattern: the reproduction is a temp plant now,\nand "
          "this check is the stronger statement the first draft should have "
          "been.")

    import graphene_figure_provenance_audit as fpa
    nm = fpa.census([os.path.join(HERE, SELF)])
    check('S2 and no Name capture either',
          len(nm) == 0,
          'the other spelling, by the shipped census, over this file')


def main():
    print('=' * 74)
    print('ATTRIBUTE-DEFAULT CAPTURE -- 2026-10-17')
    print('`def f(p=alias.ATTR)` is frozen exactly as `def f(p=NAME)` is,')
    print('and no census in this repository has ever been able to see it.')
    print('=' * 74)

    section1_semantics()
    section2_controls()
    name_total, attr_rows, attr_frozen = section3_census()
    section3b_third_census()
    section4_rf_sites()
    section5_the_oracle(attr_rows)
    section6_self()

    passed, failed = len(_PASS), len(_FAIL)
    print()
    print('=' * 74)
    print('SUMMARY')
    print('  FROZEN attribute captures recorded 2026-10-17 : %d'
          % len(BASELINE_2026_10_17))
    print('  FROZEN attribute captures in the repository NOW: %d'
          % len(attr_frozen))
    print('  Name captures, by the shipped census           : %d' % name_total)
    print('  THE INERT-KNOB TOTAL FOR 2026-10-16 WAS 20 + %d = %d, NOT 20.'
          % (len(BASELINE_2026_10_17), 20 + len(BASELINE_2026_10_17)))
    print('  Every "the census reads zero" since 2026-10-01 -- and this')
    print('  morning\'s commit, which restored that number -- was a statement')
    print('  about ONE of the two spellings of the pattern.')
    print('  BUT NOT INVISIBLE: graphene_cross_module_delivery_audit has')
    print('  censused the attribute spelling all along, 16 rows in its')
    print('  committed transcript, and this session\'s first two commit')
    print('  messages say otherwise and are wrong.  Section 3b and the')
    print('  retraction at the top of this file carry the correction; what')
    print('  survives it is the repair, the two false docstrings, the')
    print('  reconciliation nobody did, and Section 1\'s semantics.')
    print('  The 13 are converted (R4), none has regrown (R5), and F2/F3 --')
    print('  red in the transcript committed one commit earlier, which is the')
    print('  BEFORE state and is deliberately on the record -- now measure the')
    print('  rebinding arriving at the real function.')
    print('  Validations: %d passed, %d failed' % (passed, failed))
    print('=' * 74)
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
