"""graphene_capture_surface_mutation.py -- mutation test for graphene_capture_surface_guard.

A self-checking suite that passes against deliberately broken code has
measured nothing.  2026-09-17 in the sibling repository found a testbench
that scored 55/55 against broken RTL; 2026-10-23's own V18 found a refusal
branch nothing reached.  This runner injects one defect at a time into a COPY
of the repository -- never into the working tree -- and requires the suite to
FAIL, naming which checks caught it.

Each mutation is chosen to break one LOAD-BEARING property of the suite, not
to be merely detectable:

  M1  undo the 2026-10-21 repair itself (live=False reads the live name)
  M2  compare with a relative tolerance instead of bitwise
  M3  drop the "PROTECTED but not on the read surface" refusal
  M4  attribute a held object to its owner's name only, dropping aliases
  M5  make the perturbation the identity, so nothing can move
  M6  drop the pre-measurement held-object refusal
  M7  attribute held-object reads to the holder's own reads (hide the route)

Usage:  PYTHONPATH=... python3 graphene_capture_surface_mutation.py
"""

import io
import os
import re
import shutil
import subprocess
import sys

SRC = os.path.dirname(os.path.abspath(__file__))
GUARD = 'graphene_capture_surface_guard.py'
FLIP = 'graphene_default_flip.py'


def _edit(root, fname, old, new, count=1):
    p = os.path.join(root, fname)
    s = io.open(p, encoding='utf-8').read()
    if s.count(old) != count:
        raise AssertionError('%s: expected %d occurrence(s) of %r, found %d'
                             % (fname, count, old[:60], s.count(old)))
    io.open(p, 'w', encoding='utf-8').write(s.replace(old, new))


MUTATIONS = (
    ('M1', 'undo the 2026-10-21 repair: live=False reads the live name',
     lambda r: _edit(r, FLIP,
                     '        if live:\n'
                     '            return gfet.quantum_capacitance(arg, T=T)\n'
                     '        return shipped_cq(arg, T=T)',
                     '        if live:\n'
                     '            return gfet.quantum_capacitance(arg, T=T)\n'
                     '        return gfet.quantum_capacitance(arg, T=T)')),
    ('M2', 'compare with a 1e-3 relative tolerance instead of bitwise',
     lambda r: _edit(r, GUARD,
                     '    return base, moved, base.hex() == moved.hex()',
                     '    return base, moved, abs(moved - base) '
                     '<= 1e-3 * abs(base)')),
    ('M3', 'drop the PROTECTED-but-unreached refusal',
     lambda r: _edit(r, GUARD,
                     '        elif verdict == PROTECTED '
                     'and site not in read_surface:',
                     '        elif verdict == PROTECTED and False:')),
    ('M4', 'attribute a held object to one name, dropping identity aliases',
     lambda r: _edit(r, GUARD,
                     '                if val is obj and '
                     "not attr.startswith('__'):\n"
                     '                    out[(mname, attr)] = obj',
                     '                if val is obj and '
                     "not attr.startswith('__') \\\n"
                     "                        and mname == 'graphene_fet_"
                     "model':\n"
                     '                    out[(mname, attr)] = obj')),
    ('M5', 'make the perturbation the identity, so nothing can move',
     lambda r: _edit(r, GUARD,
                     '        def wrapper(*a, **kw):\n'
                     '            return k * np.asarray(obj(*a, **kw), '
                     'dtype=float)',
                     '        def wrapper(*a, **kw):\n'
                     '            return 1.0 * np.asarray(obj(*a, **kw), '
                     'dtype=float)')),
    ('M6', 'drop the pre-measurement held-object refusal',
     lambda r: _edit(r, GUARD,
                     '        if unanswered:\n'
                     '            raise UnadjudicatedCaptureSurface(',
                     '        if unanswered and False:\n'
                     '            raise UnadjudicatedCaptureSurface(')),
    ('M7', 'report no reads for held objects, hiding the route entirely',
     lambda r: _edit(r, GUARD,
                     '            reads |= set(rra.static_reads(obj)[0])',
                     '            reads |= set()')),
)


def run(root):
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join(
        x for x in (env.get('PYTHONPATH', ''), root) if x)
    proc = subprocess.run([sys.executable, GUARD], cwd=root, env=env,
                          capture_output=True, text=True, timeout=600)
    out = proc.stdout + proc.stderr
    failed = re.findall(r'^    FAIL (\S+)', out, re.M)
    crashed = 'Traceback' in out
    m = re.search(r'^VALIDATION: (\d+)/(\d+)', out, re.M)
    score = '%s/%s' % m.groups() if m else 'no score'
    return proc.returncode, failed, crashed, score, out


def main():
    print('=' * 78)
    print('MUTATION TEST: graphene_capture_surface_guard.py')
    print('=' * 78)
    print()
    base_rc, base_failed, base_crash, base_score, _ = run(SRC)
    print('    UNMUTATED  rc=%d  %s  failed=%s  crashed=%s'
          % (base_rc, base_score, base_failed or '[]', base_crash))
    ok_base = base_rc == 0 and not base_failed and not base_crash
    print('    %-4s the unmutated suite passes (the control)'
          % ('PASS' if ok_base else 'FAIL'))
    print()

    results = []
    work = os.path.join('/tmp', 'csg_mut')
    for tag, why, apply in MUTATIONS:
        root = os.path.join(work, tag)
        shutil.rmtree(root, ignore_errors=True)
        os.makedirs(os.path.dirname(root), exist_ok=True)
        shutil.copytree(SRC, root, ignore=shutil.ignore_patterns(
            '.git', '*.png', '__pycache__', 'notes', 'thesis_draft'))
        try:
            apply(root)
        except AssertionError as exc:
            results.append((tag, why, False, 'could not inject: %s' % exc))
            shutil.rmtree(root, ignore_errors=True)
            continue
        rc, failed, crashed, score, out = run(root)
        detected = rc != 0
        detail = '%s  rc=%d  caught by %s' % (
            score, rc, ','.join(failed) if failed
            else ('a crash' if crashed else 'nothing'))
        results.append((tag, why, detected, detail))
        shutil.rmtree(root, ignore_errors=True)

    print('-' * 78)
    for tag, why, detected, detail in results:
        print('    %-4s %s  %s' % ('PASS' if detected else 'FAIL', tag, why))
        print('         %s' % detail)
    print('-' * 78)
    ndet = sum(1 for _, _, d, _ in results if d)
    print()
    print('DETECTED: %d/%d injected defects%s'
          % (ndet, len(results),
             '' if ok_base else '   (CONTROL FAILED -- score is meaningless)'))
    print('=' * 78)
    return 0 if (ndet == len(results) and ok_base) else 1


if __name__ == '__main__':
    sys.exit(main())
