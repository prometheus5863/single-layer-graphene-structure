"""graphene_value_read_mutation.py -- mutation test for
graphene_value_read_surface.

A suite that passes against deliberately broken code has measured nothing.
This runner injects one defect at a time into a COPY of the repository --
never into the working tree -- and requires the suite to FAIL, naming which
checks caught it.

Each mutation breaks one LOAD-BEARING property, not merely a detectable one:

  M1  delete the bare-scalar branch -- the whole finding
  M2  report EVERY module scalar, read or not (the over-reporting failure the
      V5b/V5c controls exist for)
  M3  compare with a 1e-3 relative tolerance instead of bitwise
  M4  key the module table on `m.__name__` again -- the bug fixed while
      writing this module, so the fix is regression-tested
  M5  drop the pre-measurement refusal for unadjudicated bare-scalar reads
  M6  drop the UNREACHED-is-unavailable refusal
  M7  widen the apparatus exclusion to also drop the two findings of the day
  M8  count every scalar as blind, destroying the census's precision
  M9  walk only the holder itself, not the functions it calls

Usage:  PYTHONPATH=... python3 graphene_value_read_mutation.py
"""

import io
import os
import re
import shutil
import subprocess
import sys

SRC = os.path.dirname(os.path.abspath(__file__))
VRS = 'graphene_value_read_surface.py'
GUARD = 'graphene_capture_surface_guard.py'


def _edit(root, fname, old, new, count=1):
    p = os.path.join(root, fname)
    s = io.open(p, encoding='utf-8').read()
    if s.count(old) != count:
        raise AssertionError('%s: expected %d occurrence(s) of %r, found %d'
                             % (fname, count, old[:60], s.count(old)))
    io.open(p, 'w', encoding='utf-8').write(s.replace(old, new))


MUTATIONS = (
    ('M1', 'delete the bare-scalar branch -- the whole finding',
     lambda r: _edit(r, VRS,
                     "            elif owner is not None and isinstance(base, (int, float)) \\\n"
                     "                    and not isinstance(base, bool) \\\n"
                     "                    and getattr(mods[owner], name, rra._MISSING) is base:\n"
                     "                vals.add((owner, name))",
                     "            elif False:\n"
                     "                vals.add((owner, name))")),
    ('M2', 'report EVERY module scalar, read or not',
     lambda r: _edit(r, VRS,
                     '    walk(fn.__code__, fn.__globals__, rra._closure_map(fn), 0)\n'
                     '    return vals',
                     '    walk(fn.__code__, fn.__globals__, rra._closure_map(fn), 0)\n'
                     '    return vals | blind_class()[0]')),
    ('M3', 'compare with a 1e-3 relative tolerance instead of bitwise',
     lambda r: _edit(r, GUARD,
                     '    return base, moved, base.hex() == moved.hex()',
                     '    return base, moved, abs(moved - base) '
                     '<= 1e-3 * abs(base)')),
    ('M4', "key the module table on `m.__name__` again (the fixed bug)",
     lambda r: _edit(r, VRS,
                     '    mod_by_id = {id(m): k for k, m in mods.items()}\n'
                     '    glob_owner = {id(m.__dict__): k for k, m in mods.items()}',
                     '    mod_by_id = {id(m): m.__name__ for m in mods.values()}\n'
                     '    glob_owner = {id(m.__dict__): m.__name__\n'
                     '                  for m in mods.values()}')),
    ('M5', 'drop the pre-measurement unadjudicated-bare-read refusal',
     lambda r: _edit(r, VRS,
                     '    if missing:\n'
                     '        raise UnadjudicatedValueRead(',
                     '    if missing and False:\n'
                     '        raise UnadjudicatedValueRead(')),
    ('M6', 'drop the UNREACHED-is-unavailable refusal',
     lambda r: _edit(r, VRS,
                     '        if verdict == UNREACHED:\n'
                     '            wrong.append((site, verdict,',
                     '        if verdict == UNREACHED and False:\n'
                     '            wrong.append((site, verdict,')),
    ('M7', "widen the apparatus exclusion to hide the day's two findings",
     lambda r: _edit(r, VRS,
                     "_PLANTED_APPARATUS = (('graphene_value_read_surface', "
                     "'_PLANTED_BARE'),)",
                     "_PLANTED_APPARATUS = (('graphene_value_read_surface', "
                     "'_PLANTED_BARE'),\n"
                     "                      ('graphene_quantum_capacitance_"
                     "argument', 'A_FERMI'),\n"
                     "                      ('graphene_quantum_capacitance_"
                     "argument', 'CQ_PREFACTOR'))")),
    ('M8', "count every scalar as blind, destroying the census's precision",
     lambda r: _edit(r, VRS,
                     '                    if (mname, x.argval) in total:\n'
                     '                        blind.add((mname, x.argval))',
                     '                    if (mname, x.argval) in total:\n'
                     '                        blind |= total')),
    ('M9', 'walk only the holder itself, not the functions it calls',
     lambda r: _edit(r, VRS,
                     '            elif isinstance(base, types.FunctionType) \\\n'
                     '                    and id(base.__globals__) in glob_owner:\n'
                     '                walk(base.__code__, base.__globals__,\n'
                     '                     rra._closure_map(base), depth + 1)',
                     '            elif isinstance(base, types.FunctionType) \\\n'
                     '                    and id(base.__globals__) in glob_owner:\n'
                     '                pass')),
)


def run(root):
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join(
        x for x in (env.get('PYTHONPATH', ''), root) if x)
    proc = subprocess.run([sys.executable, VRS], cwd=root, env=env,
                          capture_output=True, text=True, timeout=900)
    out = proc.stdout + proc.stderr
    failed = re.findall(r'^    FAIL (\S+)', out, re.M)
    crashed = 'Traceback' in out
    m = re.search(r'^VALIDATION: (\d+)/(\d+)', out, re.M)
    score = '%s/%s' % m.groups() if m else 'no score'
    return proc.returncode, failed, crashed, score, out


def main():
    print('=' * 78)
    print('MUTATION TEST: graphene_value_read_surface.py')
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
    work = os.path.join('/tmp', 'vrs_mut')
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
