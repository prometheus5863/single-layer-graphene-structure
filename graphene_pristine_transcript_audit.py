#!/usr/bin/env python3
"""
graphene_pristine_transcript_audit.py -- 2026-10-02

THE 2026-10-01 TOP ITEM, EXECUTED AND MADE STANDING.

10-01 found that the 09-30 before/after measurement had been failing since the
commit that introduced it: it passed in the working tree of the run that wrote
it and in no clone afterwards.  The entry closed by naming the cheapest test
for that class of fault -- run the suite from a pristine clone of the commit --
and recorded that this repository had never run it.  This module is that test,
generalised from one suite to all of them, and from "does it still pass" to the
stronger question:

    A COMMITTED TRANSCRIPT IS A CLAIM ABOUT THE CODE BESIDE IT.
    DOES THE COMMITTED CODE STILL PRODUCE THE COMMITTED TRANSCRIPT?

Those are different questions.  A suite can exit 0, report all its checks
green, and still have a transcript on disk that describes a repository that no
longer exists -- because the transcripts here are captured stdout, written by
hand with a shell redirect, and a session that changes a module is free to
forget the redirect.  Nothing in this repository was checking that, and
`git status` cannot see it: the transcript is unmodified, it is the CODE that
moved out from under it.  That is the same shape as 10-01's "a citation is a
claim, and nothing was checking it", one level out -- a transcript is a
citation of an entire run.

WHY THIS READS .git/HEAD AS A FILE INSTEAD OF NAMING THE REVISION
-----------------------------------------------------------------
Section 8 of graphene_mutation_arrival_probe.py is a standing guard that no
source file here may use a relative git revision as an argument, because
10-01's fault was exactly that: a position in history came to name a different
file.  This module needs "whatever is committed right now", which is the one
legitimate meaning of that symbol -- but a guard that has to reason about
intent is not a guard.  So this module never writes the revision name: it reads
the pointer file, resolves it to a 40-hex commit id, and from then on speaks
only in content-addressed ids, which is the discipline 10-01 introduced.  The
resolved id is printed, so a reader can check what was actually run.

CLASSIFICATION
--------------
  IDENTICAL      pristine stdout == committed transcript, byte for byte
  NUMERIC_DRIFT  differs ONLY in the digits of numeric tokens, and every
                 differing pair agrees to NUMERIC_RTOL.  Reported separately
                 and NOT as a pass, because it means the transcript is not
                 the exactly-reproducible artefact it is presented as
  STALE          differs in structure, text, or by more than NUMERIC_RTOL:
                 the transcript describes code that is no longer here
  ERROR          the suite did not complete in the pristine clone

Run:  python3 graphene_pristine_transcript_audit.py
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

# --------------------------------------------------------------------------
# Tolerances and limits.  Module-level and read at call time (2026-09-30
# mutation-arrival rule): a mutation applied to these must reach the code.
# --------------------------------------------------------------------------
NUMERIC_RTOL = 1e-9       # relative agreement required to call a diff numeric
NUMERIC_ATOL = 1e-30      # absolute floor, for tokens that are legitimately 0
SUITE_TIMEOUT_S = 900     # per-suite wall clock in the pristine clone
GIT_POINTER_RELPATH = ".git/HEAD"   # a FILE path, not a revision -- see above

# (module, committed transcript).  A suite with no committed transcript is
# listed with None: it is still RUN, so a crash is caught, but there is no
# claim on disk to compare against.
SUITES = [
    ("graphene_mutation_arrival_probe.py",    "mutation_arrival_probe_output.txt"),
    ("graphene_dead_name_sweep.py",           "dead_name_sweep_output.txt"),
    ("graphene_band_structure_audit.py",      "band_structure_audit_output.txt"),
    ("graphene_default_scale_audit.py",       "default_scale_audit_output.txt"),
    ("graphene_figure_provenance_audit.py",   "figure_provenance_audit_output.txt"),
    ("graphene_gds_quadrature_audit.py",      "gds_quadrature_audit_output.txt"),
    ("graphene_log_sensitivity_step_audit.py", "log_sensitivity_step_audit_output.txt"),
    ("graphene_rootfinder_audit.py",          "rootfinder_audit_output.txt"),
    ("graphene_transport_optical_audit.py",   "transport_optical_audit_output.txt"),
    # REPAIRED 2026-10-28.  This row carried None -- "no committed
    # transcript to compare" -- from the day the table was written, and that
    # was false: default_scale_probe_audit_output.txt IS this module's
    # transcript, byte-identical to its stdout as of today.  The filename
    # does not contain the module's name, which is why no one matched them.
    # A None row is worse than an absent row: an absent one shows up as a
    # gap in the coverage census added below, while this one reported RAN
    # every session and looked like coverage.
    ("graphene_covariance_probe_audit.py",
     "default_scale_probe_audit_output.txt"),
    ("graphene_differential_crossover_model.py", "differential_crossover_output.txt"),
    # Added 2026-10-03, in the same session that created the suite and the
    # transcript.  Registering a new (code, transcript) pair in the SAME
    # commit series that creates it is the one step the 10-02 session named
    # as the repository's recurring failure -- a repair built and then not
    # connected.  An unregistered transcript is precisely a claim on disk
    # that nothing checks.
    ("graphene_cross_module_delivery_audit.py",
     "cross_module_delivery_audit_output.txt"),
    # rf_small_signal_model.py is not an audit, but its transcript carries
    # the f_T / f_max numbers Chapter 4 quotes, and 2026-10-03 changed them.
    ("rf_small_signal_model.py",               "rf_small_signal_output.txt"),
    # Added 2026-10-04, in the same session that created the suite and the
    # transcript, per the 2026-10-03 wiring rule.  This one matters more
    # than most: its transcript is the EVIDENCE for a claim that
    # contradicts Chapter 4's prose (that R_g is not the binding constraint
    # on f_max), and a transcript carrying a contradiction that nothing
    # checks is the worst case of the 10-02 rule, not the mildest.
    ("graphene_fmax_shortfall_decomposition.py",
     "fmax_shortfall_decomposition_output.txt"),
    # Added 2026-10-17, in the same commit that creates the suite, per the
    # 2026-10-03 wiring rule.  This one is registered while it is RED: F2 and
    # F3 fail because the 13 attribute captures it censuses are not converted
    # until the next commit.  Registering it red is the point -- the
    # transcript on disk is then the BEFORE state, and the next commit's
    # regeneration is the after, which is stronger evidence than a single
    # green transcript asserting the repair worked.  2026-10-16's amendment
    # to run_suite is what makes this possible at all: before it, a suite
    # exiting 1 was reported ERROR with its transcript never compared.
    ("graphene_attribute_default_capture.py",
     "attribute_default_capture_output.txt"),
    # Added 2026-10-18, in the same commit that creates the suite, per the
    # 2026-10-03 wiring rule, and RED for the same reason as the row above:
    # R2 fails because the two standing items it checks are still on the list
    # in the form the tree contradicts.  This suite is the one that most needs
    # registering rather than running by hand, because what it audits is the
    # artefact every session reads to choose its work -- and the fault it
    # found (a retired item revived from outside the repository) is invisible
    # to a session that trusts that artefact, which is every session.
    ("graphene_standing_list_audit.py",
     "standing_list_audit_output.txt"),
    # Added 2026-10-26, in the same commit that regenerates the six stale
    # transcripts this table's own STALE verdict exposed.  The 2026-10-25
    # standing item asked for every transcript to be regenerated each session
    # rather than only the ones the session touched; doing that found five
    # MIXED names and twelve frozen/live numeric literals that had been
    # invisible for six sessions, because the transcripts predate the modules
    # containing them.  A suite not in this table is a suite whose transcript
    # can go stale silently, so today's suite joins the table on the day it is
    # written rather than on the day someone notices.
    ("graphene_installed_surface_audit.py",
     "installed_surface_audit_output.txt"),
    # Added 2026-10-27, in the same commit that creates the suite, per the
    # 2026-10-03 wiring rule.  This one measures a REPAIR -- the 11 of 12
    # capture sites converted today -- against the pinned pre-repair blobs,
    # and a suite whose whole content is "the before state was broken and the
    # after state is not" is exactly the kind whose transcript must not be
    # allowed to go stale: once the blobs are no longer reachable or the
    # sentinel form is edited, the claim becomes unfalsifiable without this
    # table noticing.
    ("graphene_mixed_capture_repair_audit.py",
     "mixed_capture_repair_output.txt"),
    # ---------------------------------------------------------------------
    # Added 2026-10-28: the 2026-10-27 ranked-second item, worked.
    #
    # 10-27 estimated "18 rows against roughly 40 committed transcripts".
    # The real figure is 49 committed transcripts against 17 compared ones
    # -- the estimate understated the queue by nine, which is itself the
    # argument for counting it mechanically rather than by eye.  The
    # coverage census below now does that every run.
    #
    # Every row here was measured before it was added: run from a clean tree
    # at HEAD and compared byte-for-byte against its committed transcript.
    # The time in each comment is that measurement, and it is also the
    # budget argument -- a shell call on this VM is capped at 180 s and
    # Section 2 already spent about 70 s before today.
    ("graphene_potential_declaration.py",    "potential_declaration_output.txt"),      # 0.7s
    ("graphene_capture_leverage_figure.py",  "capture_surface_leverage_output.txt"),   # 0.9s
    ("graphene_replacement_read_audit.py",   "replacement_read_audit_output.txt"),     # 0.9s
    ("graphene_capture_surface_guard.py",    "capture_surface_guard_output.txt"),      # 1.0s
    ("graphene_write_read_guard.py",         "write_read_guard_output.txt"),           # 1.0s
    ("graphene_second_anchor_bound.py",      "second_anchor_bound_output.txt"),        # 1.0s
    ("graphene_contact_consistent_argument.py",
     "contact_consistent_argument_output.txt"),                                        # 1.0s
    ("graphene_double_regularisation.py",    "double_regularisation_output.txt"),      # 1.2s
    ("graphene_chord_decay_length.py",       "chord_decay_length_output.txt"),         # 1.2s
    ("graphene_value_read_surface.py",       "value_read_surface_output.txt"),         # 1.2s
    ("graphene_quantum_capacitance_overflow.py",
     "quantum_capacitance_overflow_output.txt"),                                       # 1.3s
    ("graphene_frozen_default_semantics.py", "frozen_default_semantics_output.txt"),   # 1.4s
    ("graphene_fet_signed_carrier_model.py", "fet_signed_carrier_output.txt"),         # 1.6s
    ("graphene_exact_charge_rederivation.py",
     "exact_charge_rederivation_output.txt"),                                          # 3.6s
    ("graphene_quantum_capacitance_argument.py",
     "quantum_capacitance_argument_output.txt"),                                       # 8.3s
    # Registered RED, per the 2026-10-17 and 10-18 precedent: this suite
    # exits 1 on 47 unregistered potential-like names, which is open work
    # and not a broken suite.  run_suite compares its transcript anyway.
    ("graphene_potential_census_audit.py",   "potential_census_output.txt"),           # 2.0s
    # audit_output.txt was STALE BY A MONTH when it was measured today, and
    # stale in the direction that matters: it claims "exact zeros : 3/3"
    # where this code has printed "BITWISE zeros : 2/3" since 2026-09-27,
    # when the stepping was corrected to symmetric.  The transcript
    # OVERSTATES the tree, which is the same direction as the nineteen-day
    # false claim 10-27 found, and for the same reason -- nothing compared
    # it.  Regenerated in this commit and registered in it.
    ("graphene_sensitivity_audit.py",        "audit_output.txt"),                      # 1.0s
]

# ---------------------------------------------------------------------------
# COVERAGE: every committed transcript is REGISTERED or EXCLUDED FOR A REASON
# ---------------------------------------------------------------------------
# The 2026-10-27 session found four stale transcripts and all four were
# outside SUITES; one had asserted a repaired defect for nineteen days.
# Today's own first commit made a fifth stale, inside one session, for the
# same reason.  The lesson is not "register more" -- it is that nothing
# could SEE the gap, so the gap could not be worked down.
#
# Below, every committed *_output.txt that is not in SUITES must appear here
# with a reason, and the census in Section 5 fails if any transcript is in
# neither table.  Reason classes, and whether the reason is CHECKED or
# merely DECLARED:
#
#   SELF_WRITING   CHECKED.  The module writes this transcript into its own
#                  working tree, so running it inside the pristine clone
#                  overwrites the file the comparison is taken against and
#                  the verdict is circular -- 2026-10-27 section 5 item 3,
#                  which is how that session nearly missed the
#                  import-binding staleness.  Registering one of these
#                  needs run_suite to compare against the object store
#                  instead, which is a change to the audit and not a row.
#   SELF           CHECKED.  This module. There is already a runtime
#                  assertion in main() forbidding it from auditing itself.
#   NOT_A_STDOUT   DECLARED.  The file is not a capture of the module's
#                  stdout at all; it is a hand-written excerpt with its own
#                  header and date.  Byte comparison is not the right
#                  instrument and pretending otherwise would make the suite
#                  permanently red for no finding.
#   UNMEASURED     DECLARED, and the honest name for "next session's work".
#                  Not yet timed against the 180 s call cap.
TRANSCRIPT_EXCLUSIONS = {
    "contact_argument_blast_radius_output.txt":
        ("SELF_WRITING", "graphene_contact_argument_blast_radius.py"),
    "default_flip_output.txt":
        ("SELF_WRITING", "graphene_default_flip.py"),
    "diffusion_current_output.txt":
        ("SELF_WRITING", "graphene_diffusion_current_model.py"),
    "diffusion_mutation_output.txt":
        ("SELF_WRITING", "graphene_diffusion_mutation.py"),
    "import_binding_sweep_output.txt":
        ("SELF_WRITING", "graphene_import_binding_sweep.py"),
    "perfect_contact_counterfactual_output.txt":
        ("SELF_WRITING", "graphene_perfect_contact_counterfactual.py"),
    # The three below are ALSO over the 180 s call cap (2026-10-27 section 5
    # item 4, measured there), so they carry two independent reasons.
    "perfect_contact_mutation_output.txt":
        ("SELF_WRITING", "graphene_perfect_contact_mutation.py"),
    "velocity_saturation_output.txt":
        ("SELF_WRITING", "graphene_velocity_saturation_model.py"),
    "velocity_saturation_mutation_output.txt":
        ("SELF_WRITING", "graphene_velocity_saturation_mutation.py"),
    "pristine_transcript_audit_output.txt":
        ("SELF", "graphene_pristine_transcript_audit.py"),
    "sensitivity_converged_output.txt":
        ("NOT_A_STDOUT", "graphene_crossover_sensitivity_model.py"),
    "capture_surface_mutation_output.txt":
        ("UNMEASURED", "graphene_capture_surface_mutation.py"),
    "installed_surface_mutation_output.txt":
        ("UNMEASURED", "graphene_installed_surface_mutation.py"),
    "value_read_mutation_output.txt":
        ("UNMEASURED", "graphene_value_read_mutation.py"),
}


def writes_own_transcript(repo, module, transcript):
    """Does `module` open `transcript` for writing?  Source-level, no import."""
    path = os.path.join(repo, module)
    if not os.path.exists(path):
        return False
    try:
        src = open(path).read()
    except Exception:
        return False
    return bool(re.search(
        r"open\(\s*(?:os\.path\.join\([^)]*)?[\"']%s[\"']\s*\)?\s*,\s*[\"']w"
        % re.escape(transcript), src, re.S))


def committed_transcripts(repo):
    """Every *_output.txt tracked by git at the working tree's HEAD."""
    out = git(["ls-files", "*_output.txt"], repo)
    return sorted(l.strip() for l in out.split("\n") if l.strip())

# This module is deliberately absent from SUITES.  It would have to run itself
# inside its own pristine clone, and its transcript would then be a claim that
# only it can check -- the self-certification 09-30 ruled out.
#
# SELF IS ASSERTED ON, NOT MERELY DEFINED.  Its first form was written and
# never read, and graphene_dead_name_sweep.py did not report it DEAD -- it
# reported CROSS_MODULE, resolved against a FUNCTION-LOCAL variable of the same
# spelling in graphene_log_sensitivity_step_audit.py, a module that does not
# import this one.  That is the same-spelling half of the loophole 10-01 found
# in the sweep and fixed only the mention-vs-use half of, and the sweep's new
# Section 2b exists because this module walked straight into it.  The comment
# above is the REASON this module is excluded; the assertion below is what
# makes the exclusion checkable, which a comment is not.
SELF = os.path.basename(__file__)

_PASS = []
_FAIL = []


def check(label, ok, detail=""):
    (_PASS if ok else _FAIL).append(label)
    print("  [%s] %s" % ("PASS" if ok else "FAIL", label))
    if detail:
        for line in detail.rstrip("\n").split("\n"):
            print("        " + line)
    return ok


# --------------------------------------------------------------------------
# Section 1 -- resolve the committed state to a content-addressed id
# --------------------------------------------------------------------------
def git(args, cwd, check_rc=True):
    p = subprocess.run(["git"] + args, cwd=cwd, capture_output=True, text=True)
    if check_rc and p.returncode != 0:
        raise RuntimeError("git %s failed: %s" % (" ".join(args), p.stderr.strip()))
    return p.stdout.strip()


def resolve_committed_id(repo):
    """Read the pointer FILE and resolve it to a 40-hex commit id.

    Never names the revision: see the module docstring.  Handles both a
    symbolic pointer (normal checkout) and a detached 40-hex pointer.
    """
    pointer = os.path.join(repo, GIT_POINTER_RELPATH)
    with open(pointer) as fh:
        raw = fh.read().strip()
    if raw.startswith("ref: "):
        refname = raw[5:].strip()
        commit = git(["rev-parse", "--verify", refname], repo)
    else:
        commit = raw
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise RuntimeError("pointer did not resolve to a commit id: %r" % commit)
    return commit


def make_pristine_clone(repo, commit, dest):
    """A clone carries only COMMITTED objects, so the working tree cannot leak
    in.  `git archive` would be smaller but strips .git, and at least one suite
    reads a pinned blob out of the object store, so a real clone is required.
    """
    git(["clone", "--quiet", "--no-checkout", "file://" + os.path.abspath(repo), dest], ".")
    git(["checkout", "--quiet", "--detach", commit], dest)
    got = git(["rev-parse", "--verify", commit], dest)
    if got != commit:
        raise RuntimeError("clone is not at the requested commit")
    return dest


# --------------------------------------------------------------------------
# Section 2 -- the classifier
# --------------------------------------------------------------------------
_NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def _skeleton(line):
    """The line with every numeric token replaced by a placeholder."""
    return _NUM.sub("\x00", line)


def _numbers(line):
    return [float(m.group(0)) for m in _NUM.finditer(line)]


def close_enough(a, b):
    if a == b:
        return True
    denom = max(abs(a), abs(b))
    if denom <= NUMERIC_ATOL:
        return True
    return abs(a - b) / denom <= NUMERIC_RTOL


def classify(committed, produced):
    """IDENTICAL / NUMERIC_DRIFT / STALE, plus evidence."""
    if committed == produced:
        return "IDENTICAL", []
    cl = committed.split("\n")
    pl = produced.split("\n")
    if len(cl) != len(pl):
        return "STALE", ["line count %d -> %d" % (len(cl), len(pl))]
    worst = 0.0
    evidence = []
    for i, (a, b) in enumerate(zip(cl, pl)):
        if a == b:
            continue
        if _skeleton(a) != _skeleton(b):
            evidence.append("line %d is not a numeric difference:" % (i + 1))
            evidence.append("  committed: " + a.strip()[:96])
            evidence.append("  pristine : " + b.strip()[:96])
            return "STALE", evidence
        na, nb = _numbers(a), _numbers(b)
        if len(na) != len(nb):
            return "STALE", ["line %d: numeric token count changed" % (i + 1)]
        for x, y in zip(na, nb):
            if not close_enough(x, y):
                evidence.append("line %d: %r vs %r exceeds rtol %g"
                                % (i + 1, x, y, NUMERIC_RTOL))
                return "STALE", evidence
            d = max(abs(x), abs(y))
            if d > NUMERIC_ATOL:
                worst = max(worst, abs(x - y) / d)
    return "NUMERIC_DRIFT", ["worst relative difference %.3e" % worst]


# --------------------------------------------------------------------------
# Section 3 -- run a suite in the pristine clone
# --------------------------------------------------------------------------
def run_suite(clone, module):
    env = dict(os.environ)
    env["MPLBACKEND"] = "Agg"          # no display, no blocking show()
    env["PYTHONHASHSEED"] = "0"
    env["SOURCE_DATE_EPOCH"] = "0"
    import time
    t0 = time.time()
    try:
        p = subprocess.run([sys.executable, module], cwd=clone, env=env,
                           capture_output=True, text=True, timeout=SUITE_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return None, SUITE_TIMEOUT_S, "timeout after %ds" % SUITE_TIMEOUT_S, None
    dt = time.time() - t0
    # AMENDED 2026-10-16.  This used to discard p.stdout whenever the exit
    # code was nonzero and report ERROR, which conflated two different
    # things: a suite that CRASHED, and a suite that RAN to completion and
    # reported a failing guard.  Every suite here ends with
    # `return 0 if failed == 0 else 1`, so a red guard is a nonzero exit by
    # design -- and the transcript of exactly the suite most worth watching
    # was never compared at all.  That became load-bearing the same day this
    # repository began deliberately committing red transcripts (see
    # mutation_arrival_probe and dead_name_sweep).  A run is now ERROR only
    # if it produced NO usable stdout; otherwise the transcript is compared
    # as always and the exit code is carried alongside the verdict.
    if not p.stdout:
        return None, dt, "exit %d (no stdout): %s" % (
            p.returncode, p.stderr.strip()[-300:]), p.returncode
    return p.stdout, dt, None, p.returncode


def main():
    repo = os.path.dirname(os.path.abspath(__file__))
    # Runtime assertion, not a comment: this module must not audit itself.
    assert SELF not in [m for m, _ in SUITES], (
        "%s is listed in SUITES: it would run inside its own pristine clone "
        "and its transcript would be a claim only it can check" % SELF)
    print("=" * 74)
    print("PRISTINE TRANSCRIPT AUDIT -- 2026-10-02")
    print("Does the COMMITTED code still produce the COMMITTED transcripts?")
    print("=" * 74)

    print("\n" + "-" * 74)
    print("SECTION 1 -- the committed state, resolved to a content-addressed id")
    print("-" * 74)
    commit = resolve_committed_id(repo)
    print("  pointer file      : %s" % GIT_POINTER_RELPATH)
    print("  resolved commit   : %s" % commit)
    dirty = git(["status", "--porcelain", "--untracked-files=no"], repo)
    n_dirty = len([l for l in dirty.split("\n") if l.strip()])
    print("  tracked files modified in the WORKING TREE : %d" % n_dirty)
    if n_dirty:
        print("    (they are excluded by construction -- that is the point)")
        for line in dirty.split("\n")[:12]:
            if line.strip():
                print("      " + line)

    tmp = tempfile.mkdtemp(prefix="graphene_pristine_")
    clone = os.path.join(tmp, "pristine")
    try:
        make_pristine_clone(repo, commit, clone)
        print("  pristine clone    : built and detached at the id above")

        print("\n" + "-" * 74)
        print("SECTION 2 -- every suite, run from the pristine clone")
        print("-" * 74)
        results = {}
        cached = {}
        exits = {}
        for module, transcript in SUITES:
            if not os.path.exists(os.path.join(clone, module)):
                results[module] = ("MISSING", 0.0, ["not present at this commit"])
                print("  %-44s MISSING" % module)
                continue
            out, dt, err, rc = run_suite(clone, module)
            if err is not None:
                results[module] = ("ERROR", dt, [err])
                print("  %-44s ERROR    %6.1fs  %s" % (module, dt, err[:60]))
                continue
            cached[module] = out
            exits[module] = rc
            if transcript is None:
                results[module] = ("RAN", dt, ["no committed transcript to compare"])
                print("  %-44s RAN      %6.1fs  (no transcript)" % (module, dt))
                continue
            tpath = os.path.join(clone, transcript)
            if not os.path.exists(tpath):
                results[module] = ("NO_FILE", dt, [transcript + " absent"])
                print("  %-44s NO_FILE  %6.1fs" % (module, dt))
                continue
            with open(tpath) as fh:
                committed_text = fh.read()
            verdict, evidence = classify(committed_text, out)
            results[module] = (verdict, dt, evidence)
            tag = "" if rc == 0 else "  exit %d (guard red, transcript still checked)" % rc
            print("  %-44s %-14s %6.1fs%s" % (module, verdict, dt, tag))
            for e in evidence[:6]:
                print("        " + e)

        # ------------------------------------------------------------------
        print("\n" + "-" * 74)
        print("SECTION 3 -- positive controls on the detector itself")
        print("09-30's rule: a control that has never fired is a control with")
        print("no evidence that it can.  These run against REAL captured output.")
        print("-" * 74)
        donor = None
        for module, transcript in SUITES:
            if transcript and results.get(module, ("",))[0] == "IDENTICAL":
                donor = (module, transcript)
                break
        if donor is None:
            check("a donor suite classified IDENTICAL is available", False,
                  "no suite reproduced its transcript, so the controls below "
                  "cannot be run against real data")
        else:
            dmod, dtr = donor
            real = cached[dmod]
            print("  donor: %s" % dmod)

            # C1: unperturbed must stay IDENTICAL
            v, _ = classify(real, real)
            check("C1 unperturbed output classifies IDENTICAL", v == "IDENTICAL",
                  "got %s" % v)

            # C2: a WITHIN-TOLERANCE perturbation must be NUMERIC_DRIFT.
            #
            # THE FIRST FORM OF THIS CONTROL WAS WRONG AND IT FIRED ON ME, so
            # it is recorded rather than quietly re-tuned (09-30 rule).  It
            # bumped the LAST PRINTED DIGIT of a %.6f number and asserted the
            # verdict would be NUMERIC_DRIFT.  For 7.500000 -> 7.500001 that is
            # a relative change of 1.3e-7, which is a hundred times OUTSIDE
            # NUMERIC_RTOL, so STALE was the correct verdict and the control's
            # premise was false.  The lesson is that "the last digit" is a fact
            # about a FORMAT STRING and the tolerance is a fact about the
            # NUMBER, and a control must be built from the quantity it is
            # testing: this form perturbs by a known fraction of NUMERIC_RTOL,
            # so it tests the tolerance rather than the formatting.  The real
            # drift this detector was built to classify is ~2.6e-16, which is
            # seven orders of magnitude inside the tolerance -- the detector was
            # right about the real data and the synthetic control was wrong.
            m = None
            for line in real.split("\n"):
                mm = re.search(r"\d+\.\d{6,}", line)
                if mm and float(mm.group(0)) != 0.0:
                    m = (line, mm)
                    break
            if m is None:
                check("C2 a high-precision nonzero number exists to perturb", False)
            else:
                line, mm = m
                tok = mm.group(0)
                val = float(tok)
                nudged = val * (1.0 + NUMERIC_RTOL / 100.0)
                bumped = repr(nudged)
                rel = abs(nudged - val) / abs(val)
                perturbed = real.replace(line, line.replace(tok, bumped, 1), 1)
                v, ev = classify(perturbed, real)
                check("C2 a within-tolerance change classifies NUMERIC_DRIFT",
                      v == "NUMERIC_DRIFT",
                      "relative perturbation %.3e against rtol %.3e -> %s"
                      % (rel, NUMERIC_RTOL, v))
                check("C2b that perturbation is genuinely inside the tolerance",
                      rel < NUMERIC_RTOL,
                      "a control that perturbs OUTSIDE the tolerance is not "
                      "testing the tolerance;\nthis is the fault the first "
                      "form of C2 had, asserted here so it cannot return")

            # C3: a changed WORD must be STALE, never NUMERIC_DRIFT
            wline = None
            for line in real.split("\n"):
                if re.search(r"[A-Za-z]{6,}", line):
                    wline = line
                    break
            if wline is None:
                check("C3 a word-bearing line exists to perturb", False)
            else:
                mm = re.search(r"[A-Za-z]{6,}", wline)
                w = mm.group(0)
                perturbed = real.replace(wline, wline.replace(w, "zzzzzz", 1), 1)
                v, ev = classify(perturbed, real)
                check("C3 a changed word classifies STALE", v == "STALE",
                      "got %s  (%s -> zzzzzz)" % (v, w))

            # C4: a deleted line must be STALE
            lines = real.split("\n")
            if len(lines) > 4:
                perturbed = "\n".join(lines[:2] + lines[3:])
                v, _ = classify(perturbed, real)
                check("C4 a deleted line classifies STALE", v == "STALE", "got %s" % v)

            # C5: a NUMBER changed beyond tolerance must be STALE, not drift.
            # This is the one that matters: it separates "the transcript is not
            # byte-reproducible" from "the transcript is wrong".
            if m is not None:
                line, mm = m
                tok = mm.group(0)
                big = "%.6f" % (float(tok) * 2.0 + 1.0)
                perturbed = real.replace(line, line.replace(tok, big, 1), 1)
                v, _ = classify(perturbed, real)
                check("C5 a number changed beyond rtol classifies STALE",
                      v == "STALE", "got %s  (%s -> %s)" % (v, tok, big))

            # C6: END TO END.  Inject a defect into the donor's own source in a
            # throwaway clone and require the verdict to stop being IDENTICAL.
            # The five controls above test the comparator on strings; this one
            # tests that a real code change actually reaches a real verdict.
            scratch = os.path.join(tmp, "mutated")
            make_pristine_clone(repo, commit, scratch)
            src_path = os.path.join(scratch, dmod)
            with open(src_path) as fh:
                src = fh.read()
            mutated = src.replace("print(", "print_SHIM(", 1)
            injected = mutated != src
            if injected:
                mutated = ("def print_SHIM(*a, **k):\n"
                           "    import builtins\n"
                           "    builtins.print(*a, **k)\n"
                           "    builtins.print('INJECTED-DEFECT')\n") + mutated
                with open(src_path, "w") as fh:
                    fh.write(mutated)
                out2, _, err2, _rc2 = run_suite(scratch, dmod)
                if err2 is not None:
                    check("C6 the mutated donor still runs", False, err2[:200])
                else:
                    with open(os.path.join(scratch, dtr)) as fh:
                        ctext = fh.read()
                    v, _ = classify(ctext, out2)
                    check("C6 an injected defect stops the verdict being IDENTICAL",
                          v != "IDENTICAL",
                          "mutated donor classifies %s (was IDENTICAL clean)" % v)
            else:
                check("C6 a mutation site was found in the donor", False)

        # C7 -- the path added 2026-10-16: a suite that exits NONZERO but
        # ran to completion must have its stdout kept and its transcript
        # compared, not be written off as ERROR.  Before the amendment this
        # control fails, because run_suite returned None for any nonzero
        # exit.  Synthetic, so it does not depend on which real suite
        # happens to be red today.
        c7dir = tempfile.mkdtemp(prefix="graphene_pristine_c7_")
        c7mod = "red_suite.py"
        c7text = "LINE ONE\nvalue 3.000000000000\nTOTAL: 1 passed, 1 failed\n"
        with open(os.path.join(c7dir, c7mod), "w") as fh:
            fh.write("import sys\n"
                     "sys.stdout.write(%r)\n"
                     "raise SystemExit(1)\n" % c7text)
        out7, _dt7, err7, rc7 = run_suite(c7dir, c7mod)
        check("C7 a suite that exits nonzero still yields its stdout",
              err7 is None and out7 == c7text and rc7 == 1,
              "exit %r, stdout kept: %s -- a red guard is a nonzero exit BY "
              "DESIGN in every suite here (`return 0 if failed == 0 else 1`), "
              "so discarding its output means never checking the transcript "
              "of the suite most worth watching"
              % (rc7, out7 is not None))
        if out7 is not None:
            v7, _ = classify(c7text, out7)
            check("C7b and its transcript is then compared like any other",
                  v7 == "IDENTICAL",
                  "a red suite can be byte-reproducible; red and STALE are "
                  "independent properties and this audit asks only about the "
                  "second.  got %s" % v7)

        # ------------------------------------------------------------------
        print("\n" + "-" * 74)
        print("SECTION 4 -- the standing requirement")
        print("-" * 74)
        stale = sorted(m for m, (v, _, _) in results.items() if v == "STALE")
        drift = sorted(m for m, (v, _, _) in results.items() if v == "NUMERIC_DRIFT")
        errored = sorted(m for m, (v, _, _) in results.items() if v in ("ERROR", "NO_FILE", "MISSING"))
        ident = sorted(m for m, (v, _, _) in results.items() if v == "IDENTICAL")

        red = sorted(m for m, rc in exits.items() if rc not in (0, None))
        print("  IDENTICAL      %2d   %s" % (len(ident), ", ".join(s[:34] for s in ident)))
        print("  NUMERIC_DRIFT  %2d   %s" % (len(drift), ", ".join(drift)))
        print("  STALE          %2d   %s" % (len(stale), ", ".join(stale)))
        print("  ERROR/MISSING  %2d   %s" % (len(errored), ", ".join(errored)))
        print("  exited nonzero %2d   %s" % (len(red), ", ".join(red)))
        print("    (a red guard is not a stale transcript -- reported, and")
        print("     checked for CONSISTENCY with the transcript below)")

        check("no suite fails to complete in a pristine clone", not errored,
              "\n".join("%s: %s" % (m, results[m][2][0]) for m in errored))

        # A committed transcript and its suite's exit code are two claims
        # about the same run, and they can contradict each other.  This is
        # the check that would have caught 2026-10-04 -> 2026-10-15 directly:
        # a transcript reading "0 failed" beside code that exits 1.
        contradictions = []
        for module, transcript in SUITES:
            rc = exits.get(module)
            if rc is None or transcript is None:
                continue
            tpath = os.path.join(clone, transcript)
            if not os.path.exists(tpath):
                continue
            with open(tpath) as fh:
                ttext = fh.read()
            m = re.search(r"TOTAL:?\s+(\d+)\s+passed,\s+(\d+)\s+failed", ttext)
            if not m:
                continue
            claimed_failed = int(m.group(2))
            if (claimed_failed > 0) != (rc != 0):
                contradictions.append(
                    "%s: transcript claims %d failed, suite exits %d"
                    % (module, claimed_failed, rc))
        check("every transcript's own pass/fail summary agrees with its "
              "suite's exit code", not contradictions,
              "a transcript saying '0 failed' beside code that exits nonzero "
              "is the exact state this repository shipped from 2026-10-04 to "
              "2026-10-15, and no check here could see it:\n"
              + "\n".join("  " + c for c in contradictions))
        check("no committed transcript is STALE", not stale,
              "a STALE transcript is a claim about code that is no longer here;\n"
              "it cannot be seen by `git status`, because the transcript is\n"
              "unmodified and it is the CODE that moved:\n"
              + "\n".join("  %s -> %s" % (m, results[m][1] and results[m][2][0]) for m in stale))
        check("every transcript is exactly reproducible (no numeric drift)",
              not drift,
              "NUMERIC_DRIFT means the numbers agree to %g but the bytes differ,\n"
              "so these transcripts are reproducible RESULTS and not reproducible\n"
              "ARTEFACTS.  Byte comparison cannot certify them across environments."
              % NUMERIC_RTOL)

        # ------------------------------------------------------------------
        print("\n" + "-" * 74)
        print("SECTION 5 -- coverage: is every committed transcript ACCOUNTED FOR?")
        print("Created 2026-10-28.  The 10-27 session found four stale")
        print("transcripts and all four were outside this table; one had")
        print("asserted a repaired defect for nineteen days.  What was")
        print("missing was not rows -- it was a count nobody had to make.")
        print("-" * 74)
        all_tr = committed_transcripts(repo)
        registered = {t for _, t in SUITES if t}
        excluded = set(TRANSCRIPT_EXCLUSIONS)
        unaccounted = sorted(set(all_tr) - registered - excluded)
        by_class = {}
        for t in sorted(excluded):
            by_class.setdefault(TRANSCRIPT_EXCLUSIONS[t][0], []).append(t)
        print("  committed *_output.txt   %3d" % len(all_tr))
        print("  registered in SUITES     %3d" % len(registered))
        print("  excluded with a reason   %3d" % len(excluded))
        for cls in sorted(by_class):
            print("      %-14s %2d   %s"
                  % (cls, len(by_class[cls]),
                     ", ".join(t[:-len("_output.txt")] for t in by_class[cls])))
        print("  UNACCOUNTED              %3d%s"
              % (len(unaccounted),
                 "   " + ", ".join(unaccounted) if unaccounted else ""))

        check("every committed transcript is either registered in SUITES or "
              "excluded with a reason", not unaccounted,
              "neither compared nor excluded, so nothing can tell whether "
              "they are stale:\n" + "\n".join("  " + t for t in unaccounted)
              if unaccounted else
              "%d transcripts: %d compared, %d excluded.  On 2026-10-27 the "
              "same tree was 49 / 17 / 0-with-a-reason, and the shortfall "
              "was invisible because nothing counted it."
              % (len(all_tr), len(registered), len(excluded)))

        # The exclusions whose reason is a FACT ABOUT THE SOURCE are checked
        # here rather than believed.  SELF_WRITING is the load-bearing one:
        # it is the reason that would silently become false if a module
        # stopped writing its own transcript, at which point the transcript
        # should be registered and nothing would say so.
        sw_bad, sw_ok = [], 0
        for t in sorted(excluded):
            cls, mod = TRANSCRIPT_EXCLUSIONS[t]
            if cls != "SELF_WRITING":
                continue
            if writes_own_transcript(repo, mod, t):
                sw_ok += 1
            else:
                sw_bad.append("%s: %s does not open it for writing" % (t, mod))
        check("every SELF_WRITING exclusion is TRUE of the source",
              not sw_bad,
              "\n".join("  " + b for b in sw_bad) if sw_bad else
              "%d modules, each one confirmed to write the transcript it is "
              "excused from comparing" % sw_ok)

        # Hygiene, the same shape as the 7b registry check added today: an
        # exclusion naming a transcript that no longer exists is a row that
        # describes nothing, and it hides the fact that the file is gone.
        dead = sorted(excluded - set(all_tr))
        dup = sorted(excluded & registered)
        check("no exclusion names a transcript that is not committed, and "
              "none is excluded AND registered",
              not dead and not dup,
              "dead: %s; both: %s" % (dead, dup) if (dead or dup) else
              "%d exclusions, each naming a committed file, none of them "
              "also in SUITES" % len(excluded))

        # Planted controls for the three checks above, on synthetic sets --
        # no file is written into the tree other audits read.
        print()
        print("  C5 -- the coverage predicates, on PLANTED inputs")
        planted_all = ["a_output.txt", "b_output.txt"]
        planted_reg = {"a_output.txt"}
        check("  C5.1 a transcript in neither table is reported UNACCOUNTED",
              sorted(set(planted_all) - planted_reg - set()) == ["b_output.txt"],
              "the set arithmetic Section 5 runs on, on two named files")
        check("  C5.2 SELF_WRITING is refuted for a module that does not "
              "write the file",
              not writes_own_transcript(
                  repo, "graphene_mutation_arrival_probe.py",
                  "mutation_arrival_probe_output.txt"),
              "the arrival probe's transcript is captured by redirection, "
              "not written by the module, so the predicate must say no")
        check("  C5.3 and CONFIRMED for one that does",
              writes_own_transcript(repo, "graphene_default_flip.py",
                                    "default_flip_output.txt"),
              "graphene_default_flip.py line ~745 opens it 'w'")

        print("\n" + "=" * 74)
        print("SUMMARY")
        print("=" * 74)
        print("  commit audited : %s" % commit)
        for lab in _PASS:
            print("  PASS  %s" % lab)
        for lab in _FAIL:
            print("  FAIL  %s" % lab)
        print("\nTOTAL: %d passed, %d failed" % (len(_PASS), len(_FAIL)))
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
