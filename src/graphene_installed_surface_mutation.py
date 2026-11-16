"""graphene_installed_surface_mutation.py -- mutation test for
graphene_installed_surface_audit.

A suite that passes against deliberately broken code has measured nothing,
and this repository's recorded dominant failure mode is a check that passes
its plausible reading and is false under mutation.  This runner injects one
defect at a time into a COPY of the repository -- never into the working
tree -- and requires the suite to FAIL, naming which checks caught it.

Each mutation breaks one LOAD-BEARING property of 2026-10-26's result, not
merely a detectable one:

  N1  compute the surface WITHOUT entering the flip -- the whole finding, and
      the state every previous session was in
  N2  take the counterfactual WITHOUT entering the flip, leaving the surface
      installed: the half of N1 that a surface-only repair would miss
  N3  give `flip_invariant` a default of True -- the silent-default the guard
      refuses to have, because 2026-10-25 passed with the question unasked
  N4  drop the PROTECTED-must-be-exact-AS-INSTALLED refusal: the one refusal
      the 2026-10-24 and 2026-10-25 guards could not make
  N5  collapse the four-valued dependence axis to two values, so
      "moved" and "moved in a way that reverses a verdict" share a bucket
  N6  put `vrs.ADJ_CONTACT` back as the capture holder's adjudication -- the
      mistake this module's own first draft made, which measures CLEAN in
      isolation because the two contact holders are bitwise indistinguishable
      there.  A regression test for a defect the inner guard caught
  N7  relax the UNREACHED rule from the installed READ surface to the full
      installed surface, which would refuse a name held in a closure cell --
      the fix made while writing this module
  N8  compare the holder's own value with a relative tolerance instead of
      bitwise, so a 6.4x change in the contact replacement still reads as
      flip-invariant
  N9  measure the rebind set INSIDE the flip -- the nested-entry bug made
      while exploring this, which reports 10 (flip, site) pairs instead of 8

TWO DEFECTS IN THE SUITE UNDER TEST WERE FOUND BY THIS RUNNER AND REPAIRED
THERE, NOT HERE.  N4 (delete the PROTECTED-as-installed refusal) went
undetected because the first draft's G6 planted its false verdict on a holder
whose adjudication was incomplete in that configuration, so the missing-name
refusal fired first and the branch under test was never reached -- a control
that was green and could not have been anything else.  N8 went undetected
because A4 asserted the two hex values and never the invariance FLAG computed
from them.  Fifth consecutive session in which a check that passed its
plausible reading was false under mutation.

Usage:  PYTHONPATH=... python3 graphene_installed_surface_mutation.py
"""

import io
import os
import re
import shutil
import subprocess
import sys
import tempfile

SRC = os.path.dirname(os.path.abspath(__file__))
ISA = 'graphene_installed_surface_audit.py'


def _edit(root, fname, old, new, count=1):
    p = os.path.join(root, fname)
    s = io.open(p, encoding='utf-8').read()
    if s.count(old) != count:
        raise AssertionError('%s: expected %d occurrence(s) of %r, found %d'
                             % (fname, count, old[:60], s.count(old)))
    io.open(p, 'w', encoding='utf-8').write(s.replace(old, new))


MUTATIONS = (
    ('N1', 'compute the surface WITHOUT entering the flip -- the whole '
           'finding',
     lambda r: _edit(r, ISA,
                     '''def installed_surface(holder, factory):
    """`vrs.complete_read_surface(holder)` computed with the flip ENTERED."""
    with factory():
        return vrs.complete_read_surface(holder)''',
                     '''def installed_surface(holder, factory):
    """`vrs.complete_read_surface(holder)` computed with the flip ENTERED."""
    return vrs.complete_read_surface(holder)''')),

    ('N2', 'take the counterfactual WITHOUT entering the flip',
     lambda r: _edit(r, ISA,
                     '''def installed_leverage(holder, site, factory):
    with factory():
        return csg.leverage(holder, site)''',
                     '''def installed_leverage(holder, site, factory):
    return csg.leverage(holder, site)''')),

    ('N3', 'give flip_invariant a silent default of True',
     lambda r: _edit(r, ISA,
                     'def require_installed_adjudication(label, holder, factory, adjudication,\n'
                     '                                   flip_invariant, report=None):',
                     'def require_installed_adjudication(label, holder, factory, adjudication,\n'
                     '                                   flip_invariant=True, report=None):\n'
                     '    flip_invariant = True if flip_invariant is None else flip_invariant')),

    ('N4', 'drop the PROTECTED-must-be-exact-AS-INSTALLED refusal',
     lambda r: _edit(r, ISA,
                     "        elif verdict == PROTECTED and not exact:\n"
                     "            wrong.append((site, verdict,\n"
                     "                          'claimed PROTECTED but a rebind of %s moves the '",
                     "        elif verdict == PROTECTED and False:\n"
                     "            wrong.append((site, verdict,\n"
                     "                          'claimed PROTECTED but a rebind of %s moves the '")),

    ('N5', 'collapse the four-valued dependence axis to two values',
     lambda r: _edit(r, ISA,
                     '''    if iso_exact and not inst_exact:
        return LIVE_ONLY_AS_INSTALLED
    if inst_exact and not iso_exact:
        return DEAD_AS_INSTALLED''',
                     '''    if iso_exact != inst_exact:
        return BOTH_LIVE''')),

    ('N6', 'put vrs.ADJ_CONTACT back as the capture holder\'s adjudication',
     lambda r: _edit(r, ISA,
                     'ADJ_CAPTURE = dict(csg.ADJUDICATION_REPAIRED)',
                     'ADJ_CAPTURE = dict(vrs.ADJ_CONTACT)')),

    ('N7', 'relax the UNREACHED rule to the FULL installed surface',
     lambda r: _edit(r, ISA,
                     '        if verdict == UNREACHED and site in inst_read:',
                     '        if verdict == UNREACHED and site in sites:')),

    ('N8', "compare the holder's own value with a relative tolerance",
     lambda r: _edit(r, ISA,
                     "    return iso, inst, iso.hex() == inst.hex()",
                     "    return iso, inst, abs(iso - inst) <= 1e3 * abs(iso)")),

    ('N9', 'measure the rebind set INSIDE the flip (nested entry) -- the '
           'population bug made while exploring this, which reports 10 '
           'pairs instead of 8',
     lambda r: _edit(r, ISA,
                     "        sites = sorted(rra.rebind_set(factory))\n"
                     "        with factory():\n"
                     "            for site in sites:",
                     "        with factory():\n"
                     "            sites = sorted(rra.rebind_set(factory))\n"
                     "            for site in sites:")),
)

# The needle each mutation must be caught BY.  A bare non-zero exit would be
# satisfied by a mutant that aborts for an unrelated reason -- the fault the
# sibling repository's 2026-10-25 selftest found in its own first draft -- so
# each row names the check whose failure is the intended detection.
EXPECTED = {
    'N1': ('G1', 'G2', 'G6', 'A3', 'A5'),
    'N2': ('A1', 'A2', 'A3', 'A5', 'G1', 'G6', 'X2'),
    'N3': ('G5',),
    'N4': ('G6',),
    'N5': ('C3', 'A3'),
    'N6': ('G1',),
    'N7': ('G1',),
    'N8': ('A4',),
    'N9': ('P1', 'P2'),
}


def run(root):
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join(
        x for x in (env.get('PYTHONPATH', ''), root) if x)
    proc = subprocess.run([sys.executable, ISA], cwd=root, env=env,
                          capture_output=True, text=True, timeout=900)
    out = proc.stdout + proc.stderr
    failed = re.findall(r'^    FAIL (\S+)', out, re.M)
    crashed = 'Traceback' in out
    m = re.search(r'^VALIDATION: (\d+)/(\d+)', out, re.M)
    score = '%s/%s' % m.groups() if m else 'no score'
    return proc.returncode, failed, crashed, score, out


def main():
    print('=' * 78)
    print('MUTATION TEST: graphene_installed_surface_audit.py')
    print('=' * 78)
    print()
    base_rc, base_failed, base_crash, base_score, _ = run(SRC)
    print('    UNMUTATED  rc=%d  %s  failed=%s  crashed=%s'
          % (base_rc, base_score, base_failed or '[]', base_crash))
    ok_base = base_rc == 0 and not base_failed and not base_crash
    print('    %-4s N0 the unmutated suite passes (the negative control)'
          % ('PASS' if ok_base else 'FAIL'))
    print()

    results = []
    # 2026-10-27: a FIXED path under /tmp, which is shared between sessions on
    # this VM and survives them.  Another session's run had already created
    # /tmp/isa_mut with mode 700 under a different uid, so `shutil.rmtree(root,
    # ignore_errors=True)` silently did nothing, `os.makedirs` on the parent
    # succeeded because the parent existed, and `copytree` then died with
    # PermissionError on the first mutation -- AFTER the unmutated control had
    # printed PASS, so the transcript on disk began with a green control and
    # ended in a traceback.  A mutation suite that cannot create its own
    # working tree reports NO score, which is not the same as 0/9 and must
    # not be mistaken for it.  Per-run and under TMPDIR, so it cannot collide
    # with another session or with another run of itself.
    work = tempfile.mkdtemp(prefix='isa_mut_',
                            dir=os.environ.get('TMPDIR', '/tmp'))
    for tag, why, apply in MUTATIONS:
        root = os.path.join(work, tag)
        shutil.rmtree(root, ignore_errors=True)
        os.makedirs(os.path.dirname(root), exist_ok=True)
        shutil.copytree(SRC, root, ignore=shutil.ignore_patterns(
            '.git', '*.png', '__pycache__', 'notes', 'thesis_draft'))
        try:
            apply(root)
        except AssertionError as exc:
            results.append((tag, why, False, False,
                            'BADINJ could not inject: %s' % exc))
            shutil.rmtree(root, ignore_errors=True)
            continue
        rc, failed, crashed, score, out = run(root)
        detected = rc != 0
        want = EXPECTED[tag]
        on_terms = any(f in want for f in failed)
        detail = '%s  rc=%d  caught by %s   expected one of %s' % (
            score, rc, ','.join(failed) if failed
            else ('a crash' if crashed else 'nothing'), ','.join(want))
        results.append((tag, why, detected, on_terms, detail))
        shutil.rmtree(root, ignore_errors=True)

    print('-' * 78)
    for tag, why, detected, on_terms, detail in results:
        ok = detected and on_terms
        print('    %-4s %s  %s' % ('PASS' if ok else 'FAIL', tag, why))
        print('         %s%s' % (detail,
                                 '' if on_terms
                                 else '   <== NOT on its own terms'))
    print('-' * 78)
    ndet = sum(1 for _, _, d, t, _ in results if d and t)
    print()
    print('DETECTED ON THEIR OWN TERMS: %d/%d injected defects%s'
          % (ndet, len(results),
             '' if ok_base else '   (CONTROL FAILED -- score is meaningless)'))
    print('=' * 78)
    return 0 if (ndet == len(results) and ok_base) else 1


if __name__ == '__main__':
    sys.exit(main())
