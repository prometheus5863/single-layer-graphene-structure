"""
graphene_self_heating_mutation.py

Mutation test for graphene_self_heating_model.py.

A self-checking suite that has never been shown to FAIL is a suite of
unknown value.  On 2026-09-17 this repository committed a testbench that
passed 55/55 against deliberately broken RTL; the lesson was to require the
suite to fail on a defect before believing it when it passes.

Each mutant below injects ONE defect into a COPY of the model (never the
shipped file), runs the full suite, and requires it to go red.  A mutant that
SURVIVES is reported as a gap in the suite, not quietly dropped -- and this
run has survivors, which are listed at the end and promoted to the standing
list rather than hidden.
"""

import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TARGET = "graphene_self_heating_model.py"

# (label, old_text, new_text, which checks SHOULD catch it)
MUTANTS = [
    ("M-A 2*pi -> pi in the Eq. (7) branch cut",
     "return (w_OP / V_F) ** 2 / (2.0 * np.pi)",
     "return (w_OP / V_F) ** 2 / (1.0 * np.pi)",
     "B1/B2"),
    ("M-B kappa_ox offset 1.687 -> 1.600",
     "return 0.52 * np.log(T_ox) - 1.687",
     "return 0.52 * np.log(T_ox) - 1.600",
     "A1 and the A3/A4/A5 shares"),
    ("M-C R_Si loses its factor of 2 (spreading prefactor)",
     "R_Si = 1.0 / (2.0 * kappa_Si(T0) * np.sqrt(A))",
     "R_Si = 1.0 / (1.0 * kappa_Si(T0) * np.sqrt(A))",
     "A3/A4/A5 shares"),
    ("M-D R_ox uses t_ox^2 (a plausible-looking typo)",
     "R_ox = t_ox / (k_ox * A)",
     "R_ox = t_ox ** 2 / (k_ox * A)",
     "A2 total and the shares"),
    ("M-E v_sat loses the 2/pi prefactor",
     "return (2.0 / np.pi) * (w_OP / np.sqrt(np.pi * n)) * root / (N_OP + 1.0)",
     "return (w_OP / np.sqrt(np.pi * n)) * root / (N_OP + 1.0)",
     "V2 against the unoccupied form"),
    ("M-F the Bose occupation factor is dropped entirely",
     "root / (N_OP + 1.0)",
     "root / 1.0",
     "V3 monotonicity in T"),
    ("M-G self-heating gets the WRONG SIGN (channel cools)",
     "T_new = T0 + P_total * R_th",
     "T_new = T0 - P_total * R_th",
     "the dT >= 0 assertion in sweep_vds()"),
    ("M-H Eq. (4) density exponent 2.2 -> 2.0",
     "ALPHA = 2.2",
     "ALPHA = 2.0",
     "M4"),
    ("M-I Eq. (4) temperature exponent 3 -> 4",
     "BETA = 3.0",
     "BETA = 4.0",
     "M5"),
    ("M-J the Eq. (7) branch-cut GUARD is removed (clamping instead of "
     "refusing)",
     "    if np.any(n <= n_cut):",
     "    if False and np.any(n <= n_cut):",
     "V1, and C4's two-sidedness"),
]


def run_in(dirpath):
    env = dict(os.environ)
    env["MPLCONFIGDIR"] = "/tmp/mpl"
    env.setdefault("PYTHONPATH", "")
    proc = subprocess.run([sys.executable, TARGET, "--checks-only"],
                          cwd=dirpath, env=env,
                          capture_output=True, text=True, timeout=600)
    return proc.returncode, proc.stdout + proc.stderr


def main():
    print("=" * 74)
    print("MUTATION TEST -- graphene_self_heating_model.py")
    print("=" * 74)

    # baseline: the shipped file must be GREEN, or nothing below means
    # anything
    with tempfile.TemporaryDirectory() as d:
        for f in (TARGET, "graphene_fet_model.py"):
            shutil.copy(os.path.join(HERE, f), d)
        rc, out = run_in(d)
    base_green = (rc == 0)
    print("\nBASELINE (unmutated): %s  [%s]"
          % ("exit %d" % rc, "GREEN" if base_green else "RED"))
    for line in out.splitlines():
        if "passed," in line:
            print("  %s" % line.strip())
    if not base_green:
        print("\nBASELINE IS RED -- mutation results would be meaningless.")
        return 1

    detected, survived = [], []
    print("\n%-62s %s" % ("mutant", "verdict"))
    print("-" * 74)
    for label, old, new, expect in MUTANTS:
        with tempfile.TemporaryDirectory() as d:
            for f in (TARGET, "graphene_fet_model.py"):
                shutil.copy(os.path.join(HERE, f), d)
            path = os.path.join(d, TARGET)
            src = open(path, encoding="utf-8").read()
            if old not in src:
                print("%-62s %s" % (label[:62], "ANCHOR MISSING"))
                survived.append((label, expect, "mutation anchor not found"))
                continue
            if src.count(old) != 1:
                print("%-62s %s" % (label[:62], "ANCHOR AMBIGUOUS"))
                survived.append((label, expect,
                                 "anchor occurs %d times" % src.count(old)))
                continue
            open(path, "w", encoding="utf-8").write(src.replace(old, new))
            try:
                rc, out = run_in(d)
            except subprocess.TimeoutExpired:
                rc, out = -1, "TIMEOUT"
        nfail = sum(1 for ln in out.splitlines() if ln.startswith("  [FAIL]"))
        crashed = ("Traceback" in out) or ("AssertionError" in out)
        caught = (rc != 0) or nfail > 0 or crashed
        print("%-62s %s" % (label[:62],
                            "DETECTED (%d check failures%s)"
                            % (nfail, ", hard failure" if crashed else "")
                            if caught else "*** SURVIVED ***"))
        (detected if caught else survived).append((label, expect, out[:0]))

    print("-" * 74)
    print("%d/%d mutants detected, %d survived"
          % (len(detected), len(MUTANTS), len(survived)))
    if survived:
        print("\nSURVIVORS -- these are GAPS IN THE SUITE and are reported, "
              "not dropped:")
        for label, expect, why in survived:
            print("  - %s" % label)
            print("      expected to be caught by: %s" % expect)
            if why:
                print("      note: %s" % why)
    return 0 if not survived else 2


if __name__ == "__main__":
    raise SystemExit(main())
