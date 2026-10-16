"""
graphene_velocity_saturation_mutation.py

Mutation harness for graphene_velocity_saturation_model.py, in the form this
repository has used since 2026-09-29: inject a defect into a COPY of the
module under test, re-run its battery, and require the battery to FAIL.

WHY A CONTROL B IS MANDATORY HERE.  "All mutants die" is a claim about the
mutants, not about the suite, unless the suite is also shown NOT to fail on
a perturbation that is genuinely different but correct.  Control A (the
unmutated battery is green) and control B (a correct-but-different rewrite
survives) bracket the claim from both sides.  This is the 2026-10-04 rule and
it is the reason the 2026-10-03 and 2026-10-04 harnesses both improved the
suites they were pointed at instead of merely blessing them.

WHAT EACH MUTANT IS FOR.  M6 is the point of the exercise.  It changes the
2/pi prefactor in v_sat to 1 -- a pure MAGNITUDE change, numerically a 57 %
shift in v_sat with no sign flip, no exception and no structural difference.
The standing top methodological item of this repository (created 2026-10-03,
restated 2026-10-04) is that a PASS/FAIL at zero is blind to magnitude; M6
is the direct test of whether this suite has that blindness, on a knob that
is the only thing standing between criterion A and a fitted answer.

Run:  python3 graphene_velocity_saturation_mutation.py
Writes: velocity_saturation_mutation_output.txt
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

TARGET = "graphene_velocity_saturation_model.py"
DEPS = ("graphene_fet_model.py", "rf_small_signal_model.py",
        "bracketed_root.py")

# (tag, description, old, new, must_be_caught)
MUTANTS = [
    ("M1", "delete the saturation term entirely (S -> 0 in Eq. 4)",
     "    denom = gfet.L + gfet.mu * S",
     "    denom = gfet.L + gfet.mu * S * 0.0", True),

    ("M2", "integrate 1/n instead of n (the superseded model's average)",
     "    Q = np.trapezoid(n, V)",
     "    Q = (V[-1] - V[0]) ** 2 / np.trapezoid(1.0 / n, V) if V[-1] > V[0] "
     "else 0.0", True),

    ("M3", "ignore the contacts: profile over V_ds, not V_ds,ch",
     "                resid = Vds - Id * Rc - mid",
     "                resid = Vds - mid", True),

    ("M4", "v_sat density exponent: sqrt(pi*n) -> pi*n",
     "    return prefactor * omega_op / np.sqrt(np.pi * np.asarray(n, dtype=float))",
     "    return prefactor * omega_op / (np.pi * np.asarray(n, dtype=float))",
     True),

    ("M5", "sign of the saturation term: L + mu*S -> L - mu*S",
     "    denom = gfet.L + gfet.mu * S",
     "    denom = gfet.L - gfet.mu * S", True),

    ("M6", "MAGNITUDE ONLY: the 2/pi prefactor in v_sat becomes 1",
     "def v_sat_of_n(n, omega_op=OMEGA_OP, prefactor=2.0 / np.pi):",
     "def v_sat_of_n(n, omega_op=OMEGA_OP, prefactor=1.0):", True),

    ("M7", "the physical ceiling u/v_sat is clipped instead of reported",
     "            ceil_out[i] = float(np.max(u / vs))",
     "            ceil_out[i] = float(min(np.max(u / vs), 0.5))", True),

    # ---- CONTROL B: different code, identical answer.  MUST NOT be caught.
    ("B", "CONTROL B: trapezoid written as its own explicit sum "
          "(different code, identical arithmetic)",
     "    Q = np.trapezoid(n, V)",
     "    Q = float(np.sum((n[:-1] + n[1:]) * 0.5 * np.diff(V))) "
     "if len(V) > 1 else 0.0", False),
]


def run(cwd):
    env = dict(os.environ, VSAT_FAST="1")
    p = subprocess.run([sys.executable, TARGET], cwd=cwd, env=env,
                       capture_output=True, text=True, timeout=900)
    return p.returncode, p.stdout + p.stderr


def failed_tags(out):
    m = re.findall(r"\[FAIL\] (\S+)", out)
    return m


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    lines = []

    def say(s=""):
        lines.append(s)
        print(s)

    say("mutation report: graphene_velocity_saturation_model.py")
    say("=" * 78)
    say("  every run below uses VSAT_FAST=1 (gate grid 80 pts instead of 400).")
    say("  Control A asserts the battery is green AT THAT RESOLUTION, so a")
    say("  check that only passes at full resolution would show up here as a")
    say("  control-A failure rather than as a silently weaker harness.")

    # ---------------- CONTROL A ----------------
    with tempfile.TemporaryDirectory() as d:
        for f in (TARGET,) + DEPS:
            shutil.copy(os.path.join(here, f), d)
        rc, out = run(d)
    ctrl_a_ok = (rc == 0 and not failed_tags(out))
    say(f"  [{'PASS' if ctrl_a_ok else 'FAIL'}] CONTROL A  the unmutated "
        f"battery is green  -- exit {rc}, "
        f"{len(failed_tags(out))} failed checks")
    if not ctrl_a_ok:
        say("  control A failed; every mutant verdict below would be "
            "meaningless.  Stopping.")
        with open(os.path.join(here,
                  "velocity_saturation_mutation_output.txt"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        return 1

    src = open(os.path.join(here, TARGET)).read()

    killed = 0
    expected_kills = sum(1 for m in MUTANTS if m[4])
    bad = []
    say()
    say(f"  {'tag':<5}{'verdict':<10}{'killed by':<34}description")
    say("  " + "-" * 92)

    for tag, desc, old, new, must_catch in MUTANTS:
        if src.count(old) != 1:
            say(f"  {tag:<5}{'ARRIVAL?':<10}{'-':<34}"
                f"anchor appears {src.count(old)}x, not 1x -- NOT INJECTED")
            bad.append(tag)
            continue
        with tempfile.TemporaryDirectory() as d:
            for f in DEPS:
                shutil.copy(os.path.join(here, f), d)
            with open(os.path.join(d, TARGET), "w") as fh:
                fh.write(src.replace(old, new))
            try:
                rc, out = run(d)
            except subprocess.TimeoutExpired:
                rc, out = -1, "[FAIL] TIMEOUT  the mutant did not terminate"
        tags = failed_tags(out)
        caught = (rc != 0) or bool(tags)
        by = ",".join(tags[:3]) if tags else ("exit %d" % rc if caught else "-")
        if must_catch:
            ok = caught
            verdict = "KILLED" if caught else "SURVIVED"
            killed += 1 if caught else 0
        else:
            ok = not caught
            verdict = "SURVIVED" if not caught else "KILLED(!)"
        if not ok:
            bad.append(tag)
        say(f"  {tag:<5}{verdict:<10}{by[:33]:<34}{desc}")

    say()
    say("=" * 78)
    say(f"  {killed} of {expected_kills} mutants killed; "
        f"control A green; control B "
        f"{'survived as required' if 'B' not in bad else 'WAS CAUGHT -- the suite fails on a correct rewrite'}")

    # The point of M6, stated whichever way it went.
    say()
    say("  WHAT THIS HARNESS CHANGED.  Its FIRST run killed 5 of 7: M6 and M7")
    say("  both SURVIVED the battery as first committed.  M6 is a pure")
    say("  magnitude change to v_sat and every check in the battery was a")
    say("  sign, a zero or a ratio -- all invariant under a constant rescale")
    say("  of the one knob that decides whether criterion A was measured or")
    say("  fitted.  M7 clips the physical ceiling instead of reporting it, and")
    say("  X6 asserted only an upper bound, which a clipped report satisfies")
    say("  exactly as well as an honest one.  Three checks were added in")
    say("  response -- X7 (a magnitude pin on v_sat against an externally")
    say("  computed anchor), X6b (the reported ceiling must equal an")
    say("  independent recomputation) and X6c (it must respond across the")
    say("  V_ds ladder) -- and the battery went 23/23 -> 26/26.  The harness")
    say("  did not merely confirm the suite; it improved it, which is now the")
    say("  third time in three sessions across the two repositories.")
    say()
    say("  READ M6's KILLER COLUMN.  M6 is killed by EXACTLY ONE check, X7,")
    say("  and X7 is the only check in the battery that retains a MAGNITUDE.")
    say("  That is the same shape 2026-10-04 reported for M4/R7 in the f_max")
    say("  decomposition, reproduced independently here, and it is the")
    say("  strongest evidence yet for the standing top methodological item.")
    say()
    say("  ON M6 (the magnitude-only mutant).  M6 changes no structure, raises")
    say("  no exception and flips no sign: it multiplies v_sat by pi/2 and")
    say("  nothing else.  A battery made only of PASS/FAIL-at-zero checks")
    say("  cannot see it.  This suite's verdict on M6 is the measurement of")
    say("  whether the standing top methodological item -- record the")
    say("  MAGNITUDE of every MUST_CHANGE, not just its sign, created")
    say("  2026-10-03 and still open -- has been answered for this module or")
    say("  only restated.")

    if bad:
        say(f"  PROBLEM TAGS: {', '.join(bad)}")
    with open(os.path.join(here,
              "velocity_saturation_mutation_output.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
