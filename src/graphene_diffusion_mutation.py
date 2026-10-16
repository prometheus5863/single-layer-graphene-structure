"""
graphene_diffusion_mutation.py

Mutation harness for graphene_diffusion_current_model.py, in the form this
repository has used since 2026-09-29: inject a defect into a COPY of the
module under test, re-run its battery, and require the battery to FAIL.
Control A (the unmutated battery is green) and control B (a correct-but-
different rewrite SURVIVES) bracket the claim from both sides.

WHAT THIS TARGET NEEDS THAT THE LAST TWO DID NOT.  2026-10-05's harness
attacked saturation PHYSICS; 2026-10-06's attacked ARRIVAL and COHERENCE in
a counterfactual.  This module's central claim is STRUCTURAL -- that the
diffusion contribution to the channel integral is a BOUNDARY term -- so the
mutants here attack the structure of the claim as well as its arithmetic:

  M1  ARRIVAL: lambda never reaches the integrals.
  M2  HALF-ARRIVAL: lambda reaches Q_D but not S_D.  This is 2026-10-06's
      coherence trap in a new location -- a quantity that enters the model
      in TWO places with a mutation that arrives at only one.
  M3  SIGN: kappa = -lambda*dV_F/dV becomes +lambda*dV_F/dV.
  M4  MAGNITUDE ONLY: the boundary term's A_F/3 becomes A_F/2, a 50 %
      change with no sign flip and no structural difference.  The direct
      test of the standing top methodological item.
  M5  THE CONTROL OF THE CONTROL: X6's deliberately lopsided grid becomes
      the uniform one, so X6a passes for the wrong reason.  X6b exists for
      exactly this.
  M6  THE STRUCTURAL CLAIM: the boundary term is replaced by the quadrature
      route to the same number.  Numerically this is a <0.1 % change; it
      destroys the only claim in the module title.  X6a is the only check
      that reads it.
  M7  the stencil CENTRE becomes its EDGE, where np.gradient falls back to
      a one-sided difference.  This is 2026-10-06's M5 fault directly, and
      G5d -- added in response to this harness's first run -- is what reads
      it.

CONTROL C, AND WHY IT IS A CONTROL AND NOT A MUTANT.  This harness's FIRST
run carried the stencil spacing 6/399 -> 6/400 as a mutant, and it
SURVIVED.  The reflex is to file that as a hole and add a check.  Measured
first instead: the change moves f_max/f_T by -6.560e-11 relative and f_T by
5.275e-08, against a committed value stated to 6 figures.  A central
difference is accurate to O(h^2) precisely so that its answer does not
depend on h, so this mutant is SEMANTICALLY INERT at the precision the claim
is made at, and a check sensitive to it would be a check on an artefact.
It is therefore kept as control C -- it MUST survive -- with its measured
size in its own description.

That is the verification repository's 2026-10-01 rule arriving in this one
for the first time: *a constant that can be mutated with no observable
effect is not a suite weakness.* The two repositories' methodological
series have been separate until now; this is the first rule to cross.

Run:    python3 graphene_diffusion_mutation.py
Writes: diffusion_mutation_output.txt
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

TARGET = "graphene_diffusion_current_model.py"
DEPS = ("graphene_fet_model.py", "rf_small_signal_model.py",
        "graphene_velocity_saturation_model.py", "bracketed_root.py")

# (tag, description, old, new, must_be_caught)
MUTANTS = [
    ("M1", "ARRIVAL: lambda never reaches the channel integrals",
     '    d = _dd_integrals(V_g, Vds_ch, lam=lam, **kw)',
     '    d = _dd_integrals(V_g, Vds_ch, lam=0.0, **kw)', True),

    ("M2", "HALF-ARRIVAL: lambda reaches Q_D but NOT S_D (the 10-06 "
           "coherence trap)",
     '    S = d["S"] + d["Sd"]',
     '    S = d["S"]', True),

    ("M3", "SIGN: kappa = -lambda*dV_F/dV becomes +lambda*dV_F/dV",
     '        kappa = -lam * np.gradient(V_F, V)',
     '        kappa = +lam * np.gradient(V_F, V)', True),

    ("M4", "MAGNITUDE ONLY: the boundary term's A_F/3 becomes A_F/2 (+50 %)",
     '    Qd = lam * (A_FERMI / 3.0) * (n[0] ** 1.5 - n[-1] ** 1.5)',
     '    Qd = lam * (A_FERMI / 2.0) * (n[0] ** 1.5 - n[-1] ** 1.5)', True),

    ("M5", "CONTROL-OF-THE-CONTROL: X6's lopsided grid becomes the uniform "
           "one",
     '    Vn = VDS_REF * (np.linspace(0.0, 1.0, 201) ** 2.5)',
     '    Vn = VDS_REF * (np.linspace(0.0, 1.0, 201) ** 1.0)', True),

    ("M6", "THE STRUCTURAL CLAIM: the boundary term is replaced by the "
           "quadrature route to the same number",
     '    return dict(V=V, n=n, V_F=V_F, vs=vs, kappa=kappa,\n'
     '                Q=Q, Qd=Qd, Qd_quad=Qd_quad, S=S, Sd=Sd)',
     '    return dict(V=V, n=n, V_F=V_F, vs=vs, kappa=kappa,\n'
     '                Q=Q, Qd=Qd_quad, Qd_quad=Qd_quad, S=S, Sd=Sd)', True),

    ("M7", "the stencil CENTRE becomes its EDGE, where np.gradient is "
           "one-sided (the 10-06 M5 fault directly)",
     '    k = 2\n',
     '    k = 0\n', True),

    # ---- CONTROL C: a mutation that is semantically INERT at the precision
    # the committed number is stated to.  MUST NOT be counted as a hole.
    ("C", "CONTROL C: the stencil spacing 6/399 -> 6/400.  MEASURED effect on "
          "f_max/f_T: -6.560e-11 relative, eleven orders below the committed "
          "6 figures -- INERT, not a suite weakness",
     'H_VG_COMMITTED = 6.0 / 399.0',
     'H_VG_COMMITTED = 6.0 / 400.0', False),

    # ---- CONTROL B: different code, identical answer.  MUST NOT be caught.
    ("B", "CONTROL B: n**1.5 written as n*sqrt(n) (algebraically identical, "
          "different code)",
     '    Qd = lam * (A_FERMI / 3.0) * (n[0] ** 1.5 - n[-1] ** 1.5)',
     '    Qd = lam * (A_FERMI / 3.0) * (n[0] * np.sqrt(n[0])\n'
     '                                  - n[-1] * np.sqrt(n[-1]))', False),
]


def run(cwd):
    env = dict(os.environ, DIFF_FAST="1")
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
        with open(os.path.join(here, "diffusion_mutation_output.txt"),
                  "w") as fh:
            fh.write("\n".join(lines) + "\n")

    say("mutation report: graphene_diffusion_current_model.py")
    say("=" * 78)
    say("  every run below uses DIFF_FAST=1 (gate grid 60 pts instead of 240).")
    say("  Only Section 3's sweep uses that grid; every number the module's")
    say("  findings rest on is at a FIXED bias, so control A at fast")
    say("  resolution is the same claim as at full resolution.")

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
        mutated = src.replace(old, new)
        # CONTROL B of the DIFF, per the verification repository's 2026-10-01
        # rule: the replaced text must be ABSENT from the mutant, or the
        # defect was only half-injected and the verdict is about a half-mutant.
        if old in mutated:
            say("  %-5s%-11s%-34sreplaced text STILL PRESENT -- half-injected"
                % (tag, "ARRIVAL?", "-"))
            bad.append(tag)
            continue
        with tempfile.TemporaryDirectory() as d:
            for f in DEPS:
                shutil.copy(os.path.join(here, f), d)
            with open(os.path.join(d, TARGET), "w") as fh:
                fh.write(mutated)
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
    say("  WHAT TO READ IN THE KILLER COLUMN.  Each mutant should be killed")
    say("  by a DIFFERENT check, and three of them by exactly one check:")
    say("    M4 by X5 ALONE, because X5 is the only check that compares the")
    say("       boundary term's VALUE against an independent computation.")
    say("       G2's band is a factor of 3 wide and a 50 % change fits inside")
    say("       it, which is the 2026-10-03 magnitude shape again: a check")
    say("       stated as a band is blind to changes smaller than the band.")
    say("    M5 by X6b ALONE, because X6b is the only check that tests")
    say("       whether X6a's regridding regrids anything.")
    say("    M6 by X6a ALONE, because X6a is the only check that reads the")
    say("       STRUCTURE of the term rather than its value -- the")
    say("       quadrature route agrees with the closed form to <0.1 %, so")
    show = ("       every magnitude check in the battery passes on M6.")
    say(show)
    say("    M7 by G5d ALONE, because G5a compares lambda = 0 against vsat")
    say("       on the SAME stencil and both move together.  G5d, which pins")
    say("       the stencil against the COMMITTED 0.683186, was added after")
    say("       this harness's first run and is the only check that can see a")
    say("       stencil change at all.")
    say()
    say("  CONTROL C IS THE OTHER HALF OF THAT STORY.  The first run carried")
    say("  the stencil SPACING 6/399 -> 6/400 as a mutant and it survived.")
    say("  Measured before being filed as a hole: -6.560e-11 relative in")
    say("  f_max/f_T, against a value committed to 6 figures.  A central")
    say("  difference is O(h^2) accurate precisely so that it does not depend")
    say("  on h, so the mutant is INERT and a check sensitive to it would be")
    say("  a check on an artefact.  It is control C now, and it is the")
    say("  verification repository's 2026-10-01 rule -- a constant that can be")
    say("  mutated with no observable effect is not a suite weakness --")
    say("  crossing into this repository for the first time.")
    say()
    say("=" * 78)
    say("  RESULT: %d of %d killed, control A green, control B %s, "
        "%d anomalies"
        % (killed, expected_kills,
           "survived" if "B" not in bad else "CAUGHT", len(bad)))
    say("=" * 78)
    write()
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
