"""
graphene_perfect_contact_mutation.py

Mutation harness for graphene_perfect_contact_counterfactual.py, in the form
this repository has used since 2026-09-29: inject a defect into a COPY of the
module under test, re-run its battery, and require the battery to FAIL.
Control A (the unmutated battery is green) and control B (a correct-but-
different rewrite survives) bracket the claim from both sides -- the
2026-10-04 rule, and the reason the last four harnesses each improved the
suite they were pointed at instead of blessing it.

WHAT THE MUTANTS ARE FOR, and why this target needs a different set than
2026-10-05's.  That harness attacked the saturation PHYSICS.  This module
contains almost no physics: it is a counterfactual, and a counterfactual is
a claim about what was switched off.  So the mutants here attack ARRIVAL and
COHERENCE rather than arithmetic:

  M1  the counterfactual never arrives (the rc0 flag is dropped).  This is
      the 2026-09-30 item -- a mutation that does not arrive is
      indistinguishable from a system that does not respond -- applied to
      the module's own switch.
  M2  the OTHER half of the counterfactual never arrives (the rs0 flag is
      dropped).  M1 and M2 together are the test of whether the module's
      central finding -- that R_c enters in two places and removing one
      inverts the sign -- is actually measured or merely asserted.
  M3  "EXACTLY" becomes "very nearly": R_c = 0 is replaced by R_c = 1e-12,
      which restores the bisection branch.  Numerically this changes almost
      nothing; structurally it destroys the one claim in the module title.
      If the battery cannot tell 0 from 1e-12 it has no business saying
      "exactly".  On this harness's first two runs M3 SURVIVED, because X1a
      was written against a LITERAL 0.0 while the counterfactual used
      RC_ZERO: the exactness claim was being made about a value the module
      under test did not use.  X1a now reads RC_ZERO.
  M4  MAGNITUDE ONLY: the 2026-10-05 prediction 1.699 becomes 1.80, a 6 %
      change with no sign flip, no exception and no structural difference.
      This is the direct test of the standing top methodological item
      (created 2026-10-03, restated 2026-10-04 and 2026-10-05), and D1 is
      the only check in the battery that pins a magnitude.
  M5  the g_m stencil becomes one-sided, so f_T moves in the fifth digit.
      On this harness's FIRST run M5 SURVIVED: f_max/f_T is nearly blind to
      f_T here (term B is 1.11 % of the f_max denominator), so the four
      originally-pinned numbers could not see a 1e-5 shift in f_T.  Check
      X2b -- a Richardson check on the ORDER of the stencil -- was added in
      response.  A pin against the committed f_T would NOT have worked: the
      committed numbers come from np.gradient on the 400-point gate grid and
      this module differences at dVg = 1e-4, so the two differ by ~1.4e-5 as
      ESTIMATORS and a pin across them would fail on the unmutated module.
      The order of the stencil is the thing M5 actually changes.
  M6  the residual-requirement formula drops term B.  On this harness's
      FIRST run M6 SURVIVED too, because R1 asserted only that the residual
      was finite and greater than one and nothing read its VALUE.  Check R1b
      -- the same factor recomputed by bisection on f_max/f_T, sharing no
      algebra with the closed form -- was added in response.
  M7  `at_literature_geometry` stops restoring gfet's globals.  K2 exists
      for exactly this and nothing else does.

Run:    python3 graphene_perfect_contact_mutation.py
Writes: perfect_contact_mutation_output.txt
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

TARGET = "graphene_perfect_contact_counterfactual.py"
DEPS = ("graphene_fet_model.py", "rf_small_signal_model.py",
        "graphene_velocity_saturation_model.py", "bracketed_root.py")

# (tag, description, old, new, must_be_caught)
MUTANTS = [
    ("M1", "ARRIVAL: the R_c = 0 flag is dropped (counterfactual never fires)",
     "    if rc0:\n        kw[\"Rc_total\"] = RC_ZERO\n"
     "    Vb = np.array([V_g - dVg, V_g, V_g + dVg])",
     "    if False:\n        kw[\"Rc_total\"] = RC_ZERO\n"
     "    Vb = np.array([V_g - dVg, V_g, V_g + dVg])", True),

    ("M2", "ARRIVAL: the R_s = 0 flag is dropped (half the counterfactual)",
     "    Rs = 0.0 if rs0 else rf.source_access_resistance()",
     "    Rs = rf.source_access_resistance()", True),

    ("M3", "'EXACTLY' becomes 1e-12, restoring the bisection branch",
     "RC_ZERO = 0.0",
     "RC_ZERO = 1e-12", True),

    ("M4", "MAGNITUDE ONLY: the 2026-10-05 prediction 1.699 becomes 1.80",
     "PREDICTED_UNDILUTED_FACTOR = 1.699",
     "PREDICTED_UNDILUTED_FACTOR = 1.80", True),

    ("M5", "the g_m stencil becomes one-sided (f_T moves in the 5th digit)",
     "    gm = (Id[2] - Id[0]) / (2.0 * dVg)",
     "    gm = (Id[2] - Id[1]) / dVg", True),

    ("M6", "the residual-requirement formula drops term B",
     "        residual_gds_factor = (c[\"termA\"] / (need_denom - c[\"termB\"])\n"
     "                               if need_denom > c[\"termB\"] else np.inf)",
     "        residual_gds_factor = c[\"termA\"] / need_denom", True),

    ("M7", "at_literature_geometry stops restoring gfet's globals",
     "    finally:\n        gfet.W, gfet.Rc_total = W_saved, Rc_saved",
     "    finally:\n        pass", True),

    # ---- CONTROL B: different code, identical answer.  MUST NOT be caught.
    ("B", "CONTROL B: f_max/f_T formed directly as 1/(2 sqrt(denom)) instead "
          "of f_max/f_T (algebraically identical, different code)",
     "                Rg=Rg, Rs=Rs, fT=fT, fmax=fmax, ratio=fmax / fT,",
     "                Rg=Rg, Rs=Rs, fT=fT, fmax=fmax,\n"
     "                ratio=(1.0 / (2 * np.sqrt(denom))\n"
     "                       if denom > 0 else np.inf),", False),
]


def run(cwd):
    env = dict(os.environ, PCC_FAST="1")
    p = subprocess.run([sys.executable, TARGET], cwd=cwd, env=env,
                       capture_output=True, text=True, timeout=900)
    return p.returncode, p.stdout + p.stderr


def failed_tags(out):
    return re.findall(r"\[FAIL\] (\S+)", out)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    lines = []

    def say(s=""):
        lines.append(s)
        print(s)

    def write():
        with open(os.path.join(here,
                  "perfect_contact_mutation_output.txt"), "w") as fh:
            fh.write("\n".join(lines) + "\n")

    say("mutation report: graphene_perfect_contact_counterfactual.py")
    say("=" * 78)
    say("  every run below uses PCC_FAST=1 (gate grid 80 pts instead of 400).")
    say("  Only Section 1's sweeps use that grid; every number the module's")
    say("  findings rest on is at a FIXED bias and is resolution-independent,")
    say("  so control A at fast resolution is the same claim as at full.")

    # ---------------- CONTROL A ----------------
    with tempfile.TemporaryDirectory() as d:
        for f in (TARGET,) + DEPS:
            shutil.copy(os.path.join(here, f), d)
        rc, out = run(d)
    ctrl_a_ok = (rc == 0 and not failed_tags(out))
    say("  [%s] CONTROL A  the unmutated battery is green  -- exit %d, "
        "%d failed checks"
        % ("PASS" if ctrl_a_ok else "FAIL", rc, len(failed_tags(out))))
    if not ctrl_a_ok:
        say("  control A failed; every mutant verdict below would be "
            "meaningless.  Stopping.")
        say("  --- control A output tail ---")
        for l in out.strip().splitlines()[-25:]:
            say("  " + l)
        write()
        return 1

    src = open(os.path.join(here, TARGET)).read()

    killed = 0
    expected_kills = sum(1 for m in MUTANTS if m[4])
    bad = []
    survivors = []
    say()
    say("  %-5s%-11s%-34s%s" % ("tag", "verdict", "killed by", "description"))
    say("  " + "-" * 100)

    for tag, desc, old, new, must_catch in MUTANTS:
        if src.count(old) != 1:
            say("  %-5s%-11s%-34sanchor appears %dx, not 1x -- NOT INJECTED"
                % (tag, "ARRIVAL?", "-", src.count(old)))
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
        by = ",".join(tags[:4]) if tags else ("exit %d" % rc if caught else "-")
        if must_catch:
            ok = caught
            verdict = "KILLED" if caught else "SURVIVED"
            killed += 1 if caught else 0
            if not caught:
                survivors.append(tag)
        else:
            ok = not caught
            verdict = "SURVIVED" if not caught else "KILLED(!)"
        if not ok:
            bad.append(tag)
        say("  %-5s%-11s%-34s%s" % (tag, verdict, by[:33], desc))

    say()
    say("=" * 78)
    say("  %d of %d mutants killed; control A green; control B %s"
        % (killed, expected_kills,
           "survived as required" if "B" not in bad
           else "WAS CAUGHT -- the suite fails on a correct rewrite"))
    if survivors:
        say("  SURVIVORS: %s -- each one is a hole in the battery, named "
            "here rather than left implicit." % ", ".join(survivors))

    say()
    say("  WHAT TO READ IN THE KILLER COLUMN.  Three mutants are the point of")
    say("  this harness and each should be killed by a DIFFERENT check:")
    say("    M1 by the checks that compare counterfactual against baseline")
    say("       (C1-C4, D1, D2) -- if M1 is killed only by X2 then the module")
    say("       measures reproduction and not response;")
    say("    M3 by X1a ALONE, because X1a is the only check that distinguishes")
    say("       R_c = 0 from R_c very small -- a module whose title says")
    say("       EXACTLY and whose battery cannot tell 0 from 1e-12 is making")
    say("       a claim it does not test;")
    say("    M5 by X2b ALONE, because X2b is the only check that reads the")
    say("       ORDER of the g_m stencil rather than a value computed with it;")
    say("    M4 by D1 ALONE, because D1 is the only check in the battery that")
    say("       retains a MAGNITUDE.  That is another independent instance of")
    say("       the 2026-10-03 shape -- the only detector that catches a")
    say("       pure-magnitude change is the only check that keeps a")
    say("       magnitude -- and it is why the standing item asks for a")
    say("       repository-wide census rather than one more worked example.")
    say()
    say("  WHAT THIS HARNESS CHANGED, which is the point of running it.  Its")
    say("  first run killed 3 of 7 and its second 6 of 7.  The three")
    say("  survivors across those runs -- M5, M6 and M3 -- were ONE fault in")
    say("  three places: a quantity or a claim that nothing in the battery")
    say("  actually read.")
    say("    M5 moved f_T in the fifth digit and nothing noticed, because")
    say("      f_max/f_T is almost independent of f_T at this geometry --")
    say("      term B is 1.11 % of the denominator, so a 1e-5 shift in f_T")
    say("      moves the ratio by ~5e-8.  A pin against the committed f_T")
    say("      would not have worked either, because the committed value is a")
    say("      np.gradient on a 0.015 V grid and this module differences at")
    say("      1e-4: a different ESTIMATOR, differing by 1.4e-5 on the")
    say("      unmutated module.  X2b therefore reads the ORDER of the")
    say("      stencil by a Richardson ratio (4 for central, 2 for one-sided)")
    say("      and X2c reports the estimator gap instead of asserting it")
    say("      away.")
    say("    M6 moved the residual g_ds requirement and nothing noticed,")
    say("      because R1 asserted only that it was finite and > 1.  R1b")
    say("      recomputes the same factor by bisection on f_max/f_T, sharing")
    say("      no algebra with the closed form.")
    say("    M3 made R_c = 1e-12 and nothing noticed, because X1a tested a")
    say("      LITERAL 0.0 -- an exactness check pointed at a value the")
    say("      counterfactual did not use.  X1a now reads RC_ZERO.")
    say("  The battery went 17/17 -> 20/20, and Section 6b later took it to")
    say("  23/23.  This is the FIFTH consecutive")
    say("  harness in this repository to improve the suite it was pointed at")
    say("  rather than bless it.  Two of the three holes are the standing top")
    say("  methodological item (a report with no magnitude pin on it) and the")
    say("  third is the 2026-09-30 arrival item, inside an exactness check.")
    say()
    if bad:
        say("  PROBLEM TAGS: %s" % ", ".join(bad))
    write()
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
