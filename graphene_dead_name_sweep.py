"""
graphene_dead_name_sweep.py

A module-level constant that nothing reads is a PROVENANCE CLAIM WITH NO CODE
BEHIND IT.

Created 2026-10-01 to answer the 2026-09-30 item "Are there OTHER dead names
in this repository?", which existed because `ALPHA_ABS` was found dead only
because a mutation happened to target it -- i.e. by accident, and only for a
name some audit was already pointing at.  Nothing was looking for the class.

WHY A DEAD CONSTANT IS WORSE THAN AN UNUSED VARIABLE
----------------------------------------------------
Every one of these names carries a citation.  `PASSI_RC_DIRAC_OHM_UM` is six
measured numbers from a named table in a named paper; `D_CHEM_PHYS` is a
literature value for a specific physical term.  Writing the name and the
citation is a claim that the repository's results depend on that number.  If
no code reads it, the claim is false, and it is false in the direction that
flatters the work: a reader (or a thesis examiner) counts the citation as a
dependency the model honours.  2026-09-28 established that prose is a
detector; 09-30's ALPHA_ABS was the INVERTED case, prose asserting a
dependency the code does not have.  This file looks for the rest of that
class systematically instead of waiting for another accident.

THE EXEMPTION MECHANISM, AND WHY IT IS IN THE SOURCE AND NOT IN HERE
--------------------------------------------------------------------
Some unread constants are deliberate.  `FANG_ANTENNA_TEST_WAVELENGTH_NM` in
graphene_plasmonic_photodetector_model.py is one: its own comment says the
Fang et al. point is "Descriptive-only ... Deliberately NOT given an
F_max/fit here", with a reason that is still true.  A sweep with no way to
express that reports a documented decision as a fault, and a checker that
cries wolf on its own repository's considered choices is a checker that gets
ignored -- which is the failure mode that matters, because it is silent.

So the exemption lives in the SOURCE, as a

    # dead-name-exempt: <reason>

marker on the assignment line or the line above it, and the reason is
REQUIRED: an exemption with no reason is reported as a fault in its own
right.  Keeping the list here instead would put the repository's reasons in
the instrument that judges it, where nobody reading the constant would see
them.

CLASSES
  LIVE          read by some function in its own module
  BODY_ONLY     read only in the module body (the 09-30 IMPORT_USED case):
                not dead, and deliberately not a fault
  CROSS_MODULE  read by another module in this repository
  EXEMPT        declared descriptive-only in the source, with a reason
  DEAD          read nowhere, by nobody, and not declared
"""

import ast
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

CLASS_LIVE = 'LIVE'
CLASS_BODY_ONLY = 'BODY_ONLY'
CLASS_CROSS_MODULE = 'CROSS_MODULE'
CLASS_EXEMPT = 'EXEMPT'
CLASS_DEAD = 'DEAD'

EXEMPT_RE = re.compile(r'#\s*dead-name-exempt\s*:\s*(\S.*?)\s*$')


# ---------------------------------------------------------------------------
# The sweep
# ---------------------------------------------------------------------------

def module_level_constants(tree):
    """{NAME: lineno} for module-level ALL-CAPS assignments."""
    out = {}
    for node in tree.body:
        targets = []
        if isinstance(node, ast.Assign):
            targets = [t for t in node.targets if isinstance(t, ast.Name)]
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            targets = [node.target]
        for t in targets:
            n = t.id
            if n.isupper() and not n.startswith('_') and len(n) > 1:
                out.setdefault(n, node.lineno)
    return out


def _loads_in(node):
    return {n.id for n in ast.walk(node)
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}


def names_read_by_functions(tree):
    """Names loaded anywhere inside a function or class body."""
    out = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            out |= _loads_in(node)
    return out


def names_read_in_body(tree):
    """Names loaded by module-level statements that are not defs."""
    out = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        out |= _loads_in(node)
    return out


def names_referenced_at_ast_level(tree):
    """
    Every name another module could actually be USING: `Name` loads (which
    covers `from X import NAME`) and attribute names (which covers
    `mod.NAME`).  Deliberately NOT a text search.

    THE FIRST FORM OF THIS WAS A TEXT SEARCH AND IT HAD A LOOPHOLE THIS
    FILE'S OWN DOCSTRING EXERCISED (2026-10-01).  `re.search` over another
    module's source cannot tell a USE from a MENTION, so naming a constant in
    prose was enough to make it look alive -- and the first run reported
    `FANG_ANTENNA_TEST_WAVELENGTH_NM` as CROSS_MODULE because the paragraph
    above explaining why it is exempt contains its name.  Worse, this
    module's own `CLASS_LIVE` came back CROSS_MODULE purely because
    graphene_mutation_arrival_probe.py happens to define a constant with the
    same spelling: it was DEAD in here, and the loophole hid it.  A detector
    for dead names whose own dead name is concealed by its own mechanism is
    the 2026-09-28 cannot-fail class, and it was one run away from being
    committed.
    """
    out = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load):
            out.add(n.id)
        elif isinstance(n, ast.Attribute):
            out.add(n.attr)
        elif isinstance(n, ast.ImportFrom):
            for a in n.names:
                out.add(a.name)
    return out


def exemptions(src):
    """{NAME: reason} from `# dead-name-exempt:` markers, plus the no-reason
    offenders as {NAME: ''}."""
    lines = src.splitlines()
    out = {}
    for i, line in enumerate(lines):
        m = EXEMPT_RE.search(line)
        if not m:
            continue
        reason = m.group(1).strip()
        # The marker applies to the next ASSIGNMENT at or below it, skipping
        # the continuation comment lines a multi-line reason occupies.  The
        # first version of this looked only at i, i+1 and i+2, which is fine
        # for an inline marker and fails for a reason long enough to wrap --
        # and the only real exemption in this repository needs four lines, so
        # the check was reporting a declared exemption as DEAD.  Found by
        # running it, not by reading it.
        target = None
        for j in range(i, min(i + 8, len(lines))):
            if j > i and EXEMPT_RE.search(lines[j]):
                break          # a second marker: this one's scope has ended
            am = re.match(r'\s*([A-Z_][A-Z0-9_]*)\s*=', lines[j])
            if am:
                target = am.group(1)
                break
            stripped = lines[j].strip()
            if j > i and stripped and not stripped.startswith('#'):
                break          # real code that is not an assignment
        if target:
            out[target] = reason
    return out


def sweep(paths=None):
    """One row per module-level constant that is not LIVE."""
    if paths is None:
        paths = sorted(glob.glob(os.path.join(HERE, '*.py')))
    srcs = {}
    trees = {}
    for p in paths:
        with open(p, encoding='utf-8') as fh:
            srcs[p] = fh.read()
        trees[p] = ast.parse(srcs[p])
    ast_refs = {q: names_referenced_at_ast_level(trees[q]) for q in paths}
    rows = []
    n_live = 0
    for p in paths:
        consts = module_level_constants(trees[p])
        in_funcs = names_read_by_functions(trees[p])
        in_body = names_read_in_body(trees[p])
        exempt = exemptions(srcs[p])
        for name, lineno in sorted(consts.items()):
            if name in in_funcs:
                n_live += 1
                continue                                   # LIVE
            if name in exempt:
                rows.append((p, name, lineno, CLASS_EXEMPT, exempt[name]))
                continue
            others = [os.path.basename(q) for q in paths
                      if q != p and name in ast_refs[q]]
            if others:
                rows.append((p, name, lineno, CLASS_CROSS_MODULE,
                             ','.join(others)))
            elif name in in_body:
                rows.append((p, name, lineno, CLASS_BODY_ONLY, ''))
            else:
                rows.append((p, name, lineno, CLASS_DEAD, ''))
    return rows, n_live


# ---------------------------------------------------------------------------
# Controls on the sweep itself
# ---------------------------------------------------------------------------

_SYNTH = '''
"""synthetic module for the sweep's own controls."""
LIVE_ONE = 1.0
BODY_ONE = 2.0
DEAD_ONE = 3.0
EXEMPT_ONE = 4.0   # dead-name-exempt: a stated reason
EXEMPT_NO_REASON = 5.0
DERIVED = BODY_ONE * 2.0

def reads_live():
    return LIVE_ONE * 10.0
'''


def positive_control_on_the_sweep(tmpdir):
    """
    A sweep that reports nothing is indistinguishable from a repository with
    nothing to report, so it has to be shown capable of reporting.  Four
    classes are planted and all four must come back correctly -- including
    BODY_ONLY and EXEMPT, which are the two the sweep must NOT call dead.
    """
    path = os.path.join(tmpdir, 'synth_sweep_control.py')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(_SYNTH)
    rows, _n_live = sweep([path])
    got = {name: cls for _p, name, _ln, cls, _d in rows}
    expect = {
        'BODY_ONE': CLASS_BODY_ONLY,
        'DEAD_ONE': CLASS_DEAD,
        'EXEMPT_ONE': CLASS_EXEMPT,
        'EXEMPT_NO_REASON': CLASS_DEAD,
        'DERIVED': CLASS_DEAD,
    }
    # LIVE_ONE must be absent: it is read by a function
    live_absent = 'LIVE_ONE' not in got
    reason_ok = all(d for _p, n, _l, c, d in rows if c == CLASS_EXEMPT)
    return got, expect, live_absent, reason_ok


# ---------------------------------------------------------------------------
# The inverted detector: constants the PROSE names and the code does not read
# ---------------------------------------------------------------------------

def prose_claims_without_code(paths=None):
    """
    For each module, the constants its own docstring or comments NAME while no
    function in it reads them.  This is 2026-09-28's "prose is a detector"
    inverted, and it is a REPORT rather than an assertion: a docstring may
    legitimately mention a constant it does not depend on.  What it must not
    do is mention one as an input.
    """
    if paths is None:
        paths = sorted(glob.glob(os.path.join(HERE, '*.py')))
    out = []
    for p in paths:
        with open(p, encoding='utf-8') as fh:
            src = fh.read()
        tree = ast.parse(src)
        consts = module_level_constants(tree)
        in_funcs = names_read_by_functions(tree)
        doc = ast.get_docstring(tree) or ''
        comments = '\n'.join(l for l in src.splitlines() if l.lstrip().startswith('#'))
        prose = doc + '\n' + comments
        for name in sorted(consts):
            if name in in_funcs:
                continue
            if re.search(r'\b' + re.escape(name) + r'\b', prose):
                out.append((os.path.basename(p), name, consts[name]))
    return out


# ---------------------------------------------------------------------------
# A finding with a number in it: the chemical term exists twice
# ---------------------------------------------------------------------------

def modules_imported_by(tree):
    """Basenames of repository modules this module actually imports.

    Added 2026-10-02.  The CROSS_MODULE disposition says a constant is alive
    because ANOTHER module in this repository references its NAME.  That is
    only a liveness claim if the other module can actually reach this one --
    i.e. if it imports it.  Without that, CROSS_MODULE resolves on SPELLING,
    which is the same-spelling half of the loophole 10-01 found in this very
    module and fixed only the mention-vs-use half of.
    """
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                out.add(a.name.split('.')[0] + '.py')
        elif isinstance(node, ast.ImportFrom):
            if node.module and node.level == 0:
                out.add(node.module.split('.')[0] + '.py')
    return out


def unsupported_cross_module_rows(rows, paths=None):
    """CROSS_MODULE rows whose resolving module does NOT import the definer.

    Each one is a hidden DEAD: the constant is read by nothing that can see it,
    and the sweep calls it alive because some unrelated module happens to use
    the same identifier -- including, as found on 2026-10-02, a FUNCTION-LOCAL
    variable in an unrelated audit.  Returns (row, resolver, reason) triples.
    """
    if paths is None:
        paths = sorted(glob.glob(os.path.join(HERE, '*.py')))
    imports = {}
    for q in paths:
        with open(q, encoding='utf-8') as fh:
            imports[os.path.basename(q)] = modules_imported_by(ast.parse(fh.read()))
    bad = []
    for row in rows:
        path, name, lineno, cls, detail = row
        if cls != CLASS_CROSS_MODULE:
            continue
        definer = os.path.basename(path)
        for resolver in [d for d in detail.split(',') if d]:
            if definer not in imports.get(resolver, set()):
                bad.append((row, resolver,
                            '%s does not import %s' % (resolver, definer)))
    return bad


def chemical_term_copies():
    """
    `D_CHEM_PHYS` in graphene_contact_doping_nonlinear_model.py and
    `DC_ANCHOR` in graphene_per_metal_crossover_model.py are the SAME
    physical quantity -- Khomyakov et al.'s short-range chemical term at the
    physisorbed separation -- written twice, and only one of them is read.

    The physics is not wrong: the chemical term IS applied, through
    `w_cross_for_metal` -> `delta_c` -> `DC_ANCHOR` in the per-metal module,
    so no committed number is affected.  What is wrong is that there are two
    independent definitions of one literature value, one of them unread, so
    a revised extraction would move one and nothing would notice.

    The two modules cannot import from each other -- the per-metal module
    already imports the nonlinear one -- so the agreement check lives here,
    in the only file that imports both.  This is also the first code in the
    repository to read `D_CHEM_PHYS` at all.
    """
    import graphene_contact_doping_nonlinear_model as nl
    import graphene_per_metal_crossover_model as pm
    return float(nl.D_CHEM_PHYS), float(pm.DC_ANCHOR)


# ---------------------------------------------------------------------------

def main():
    import tempfile
    passed = failed = 0

    def check(label, ok, detail=''):
        nonlocal passed, failed
        if ok:
            passed += 1
            print('  [PASS] %s%s' % (label, '  -- ' + detail if detail else ''))
        else:
            failed += 1
            print('  [FAIL] %s%s' % (label, '  -- ' + detail if detail else ''))

    print('=' * 78)
    print('DEAD-NAME SWEEP -- module-level constants no code reads')
    print('2026-10-01.  Answers the 2026-09-30 item; ALPHA_ABS was found by')
    print('accident, and nothing was looking for the class.')
    print('=' * 78)

    tmp = tempfile.mkdtemp()

    print()
    print('SECTION 1 -- positive control on the sweep itself')
    got, expect, live_absent, reason_ok = positive_control_on_the_sweep(tmp)
    for name, cls in sorted(expect.items()):
        print('  planted %-18s expect %-13s got %s'
              % (name, cls, got.get(name, '(not reported)')))
    check('all five planted non-live classes come back correctly',
          all(got.get(k) == v for k, v in expect.items()))
    check('a constant read by a function is NOT reported', live_absent,
          'LIVE is the common case and must not be noise')
    check('an exemption with no reason is reported as DEAD, not EXEMPT',
          got.get('EXEMPT_NO_REASON') == CLASS_DEAD,
          'otherwise the marker is a way to silence the check for free')
    check('every EXEMPT row carries a reason', reason_ok)

    print()
    print('SECTION 2 -- the sweep over this repository')
    rows, n_live = sweep()
    for p, name, lineno, cls, detail in rows:
        print('  %-13s %-46s %-32s line %-5d %s'
              % (cls, os.path.basename(p), name, lineno, detail[:40]))
    dead = [r for r in rows if r[3] == CLASS_DEAD]
    tot = {c: sum(1 for r in rows if r[3] == c)
           for c in sorted({r[3] for r in rows})}
    tot[CLASS_LIVE] = n_live
    print('  totals: %s' % tot)
    check('no module-level constant in this repository is DEAD', not dead,
          '%d dead' % len(dead) if dead else
          '4 found and dispositioned 2026-10-01: PASSI_RC_DIRAC_OHM_UM given '
          'a reader, DELIVERED given a reader, D_CHEM_PHYS given an '
          'agreement check, FANG_ANTENNA_TEST_WAVELENGTH_NM declared exempt')

    print()
    print('SECTION 2b -- is every CROSS_MODULE resolution supported by an IMPORT?')
    print("-" * 70)
    print('  CROSS_MODULE says a constant is alive because another module here')
    print('  references its NAME.  That is a liveness claim only if the other')
    print('  module can reach this one.  10-01 fixed the mention-vs-use half of')
    print('  this loophole (text search -> AST) and left the same-spelling half.')
    bad = unsupported_cross_module_rows(rows)
    for (path, name, lineno, cls, detail), resolver, reason in bad:
        print('    UNSUPPORTED  %-46s %-24s line %-5d  %s'
              % (os.path.basename(path), name, lineno, reason))
    n_cm = sum(1 for r in rows if r[3] == CLASS_CROSS_MODULE)
    print('  CROSS_MODULE rows: %d   unsupported: %d' % (n_cm, len(bad)))
    # positive control: a synthetic unsupported resolution MUST be reported,
    # because a check that has never fired has no evidence that it can.
    synth_rows = [('/x/definer_mod.py', 'ONLY_SPELLING', 7,
                   CLASS_CROSS_MODULE, 'resolver_mod.py')]
    synth_bad = unsupported_cross_module_rows(synth_rows)
    check('2b positive control: a resolution with no import is reported',
          len(synth_bad) == 1,
          'a synthetic CROSS_MODULE row naming a module that does not exist '
          'and\ntherefore imports nothing must be flagged; got %d'
          % len(synth_bad))
    check('every CROSS_MODULE resolution is supported by a real import',
          not bad,
          'an unsupported resolution is a HIDDEN DEAD name: nothing that can\n'
          'see the constant reads it, and the sweep calls it alive because an\n'
          'unrelated module spells something the same way.\n'
          + '\n'.join('  %s.%s <- %s' % (os.path.basename(r[0][0]), r[0][1], r[2])
                       for r in bad))
    print()
    print("-" * 70)
    print('SECTION 3 -- the inverted detector: prose names it, no code reads it')
    claims = prose_claims_without_code()
    for base, name, lineno in claims:
        print('  %-46s %-32s line %d' % (base, name, lineno))
    print('  %d constant(s) named in prose and read by no function.' % len(claims))
    print('  This is a REPORT, not an assertion: a docstring may mention a')
    print('  constant it does not depend on.  What it must not do is cite one')
    print('  as an input.  Every row here is read by SOMETHING (see Section 2)')
    print('  or declared exempt, or Section 2 would have failed first.')

    print()
    print('SECTION 4 -- one literature value, two definitions')
    d_chem, dc_anchor = chemical_term_copies()
    print('  graphene_contact_doping_nonlinear_model.D_CHEM_PHYS = %.4f eV'
          % d_chem)
    print('  graphene_per_metal_crossover_model.DC_ANCHOR        = %.4f eV'
          % dc_anchor)
    check('the two copies of the short-range chemical term agree',
          d_chem == dc_anchor,
          'bitwise equal; the check exists because a revised extraction '
          'would otherwise move one copy silently')
    print('  The physics is unaffected -- the term is applied through')
    print('  w_cross_for_metal -> delta_c -> DC_ANCHOR -- so no committed')
    print('  number depends on D_CHEM_PHYS.  That is exactly why it went')
    print('  unread for ten sessions, and exactly why the duplication needs')
    print('  a check rather than a comment.')

    print()
    print('=' * 78)
    print('TOTAL: %d passed, %d failed' % (passed, failed))
    print('=' * 78)
    return failed


if __name__ == '__main__':
    sys.exit(1 if main() else 0)
