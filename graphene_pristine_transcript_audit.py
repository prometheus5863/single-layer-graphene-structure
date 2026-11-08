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

AMENDED 2026-11-04 -- THE COMPARISON NOW READS THE OBJECT STORE
---------------------------------------------------------------
Until today Section 2 read the committed transcript out of the pristine
clone's WORKING TREE.  That is sound for a suite whose transcript is captured
by shell redirection, and circular for one that writes its own transcript:
the run overwrote the file a moment before it was read, so the comparison was
`stdout == stdout` and could only say IDENTICAL.  Nine transcripts were
excluded for that reason and the exclusion was correct given the read.

`read_committed_blob` reads the blob at the resolved commit instead.  The
object store is immutable, so the circularity is removed at the root rather
than sequenced around.  Four of the nine rows are now registered -- not six,
as the standing list has claimed since 10-28; two of the six are over the
call-cap budget and were measured so today.

IT FOUND SOMETHING ON THE FIRST RUN.  `import_binding_sweep_output.txt` is
STALE: it was regenerated on 2026-10-27 by the session that regenerated every
transcript in the tree, and 28 commits later it describes a repository with
71 modules where this one has 76.  Under the old read that row could never
have been registered, and had it been, it would have reported IDENTICAL.

AMENDED 2026-11-06 -- A CRASHED SUITE WAS BEING REPORTED AS A STALE TRANSCRIPT
-----------------------------------------------------------------------------
The 2026-11-05 census reported

    graphene_mixed_capture_repair_audit.py   STALE   0.0s  exit 1
          line count 153 -> 5

and the standing requirement "no suite fails to complete in a pristine clone"
passed on that same run.  Both readings were wrong about the same event.  That
suite did not run and produce a five-line transcript: it printed its four
banner lines and died on

    ModuleNotFoundError: No module named 'scipy'

raised from graphene_contact_doping_model.py, an import three modules deep.
Reproduced 2026-11-06 in a pristine clone at 94b019c with scipy absent from
the path: stdout is exactly the four banner lines, FIVE by `split("\n")`,
which is the "-> 5" the census printed, and the exit code is 1.  The literal
bytes are kept below as CRASH_STDOUT_NO_SCIPY so the reproduction is a control
and not a memory.

THE CAUSE IS IN `run_suite`, AND IT IS THE 2026-10-16 AMENDMENT'S BLIND SIDE.
That amendment stopped discarding stdout on a nonzero exit, because every
suite here ends `return 0 if failed == 0 else 1` and a red guard is a nonzero
exit by design -- correct, and C7 is its control.  But the test it left behind
for "did this run work at all" is `if not p.stdout`, i.e. EMPTINESS.  A suite
that crashes before its first print is ERROR; a suite that crashes one line
after it has stdout, so it is compared, and whatever the comparison says
becomes a statement about the TRANSCRIPT.  The verdict then names the wrong
object entirely: the transcript was fine, the run never happened.

Worse, the two readings disagree in a way nothing could see.  STALE feeds the
"no committed transcript is STALE" check, which went red for nine rows and so
carried this one invisibly among eight real ones; the completion check reads
only the ERROR bucket, which was empty, and went green.  A crash was therefore
reported as a documentation defect and as a clean bill of health at once.

THE REPAIR IS A MECHANISM-GATED VERDICT, WHICH IS THE 2026-10-16 RULE ITSELF.
A suite here announces its own completion: it prints `TOTAL: n passed, m
failed` as its last act.  So "ran to completion" is checkable, and the
expectation is per-suite and self-calibrating -- it comes from the suite's OWN
COMMITTED TRANSCRIPT rather than from a hard-coded list of which modules print
a summary.  A run whose committed transcript carries a summary line and whose
stdout does not did not reach its own end: class INCOMPLETE, a fault, counted
by the completion check and not by the staleness one.  A suite whose
transcript never had a summary line makes no completion claim and is compared
exactly as before, so the new class cannot mis-fire on it.

The regex for that summary line already existed inline in Section 5's
exit-code consistency check.  It is now ONE module-level constant read by
both, because two spellings of one predicate is a fault this repository has
been bitten by twice (the A1 trigger, 2026-11-05; the v_F names, 2026-10-25).

WHAT THIS DOES NOT CLAIM.  With scipy on the path today no suite is
INCOMPLETE, so the new rule and the superseded one return the SAME verdict for
all 43 registered suites.  The overlap is reported in Section 4 rather than
asserted, and the divergence is shown by the controls instead.

AMENDED 2026-11-06 -- THIS MODULE WAS THE ONE SUITE ITS OWN CHECK EXEMPTED
--------------------------------------------------------------------------
Section 5 asks, of every registered suite, whether the transcript's own
pass/fail summary agrees with the suite's exit code, on the grounds that "a
transcript saying '0 failed' beside code that exits nonzero is the exact state
this repository shipped from 2026-10-04 to 2026-10-15".  `main()` ended with a
bare `return 0`.

So this module exited 0 on every run, whatever it found, and its committed
transcript has read

    TOTAL: 24 passed, 2 failed

beside an exit code of 0 since at least 2026-11-04 -- the mirror image of the
state the check exists to catch, in the file that does the checking.  Measured
2026-11-06, not inferred: the run printed two FAILs and `$?` was 0.

NOTHING COULD SEE IT, AND THE REASON IS A GOOD EXCLUSION OVERREACHING.  This
module is deliberately absent from SUITES, because auditing itself inside its
own pristine clone would make its transcript a claim only it can check -- that
is the 09-30 self-certification rule and it is still right.  But it was being
read as an exemption from EVERY property, and the exit-code property is not
like the byte comparison: it needs no pristine clone, no second run, and no
self-certification.  It is a fact about one function's last statement.  An
exclusion justified for one check had silently acquired the authority of all
of them.

So the repair is not to register the module.  It is a control, C9, that asks
the question directly of the source: every `return` in `main()` must be the
convention's conditional, and none may be a bare constant.  Checked by AST
rather than by grepping for the text, because a source-level regex for this
convention was tried first and reported 27 of 43 suites violating it, which is
false -- five suites were observed exiting nonzero on this very run.  A
measurement that disagrees with a direct observation is the measurement that
is wrong, and it is not reported here as a finding.  Whether every OTHER suite
really obeys the convention is therefore UNMEASURED, and says so on the
standing list.

Run:  python3 graphene_pristine_transcript_audit.py
"""

import base64
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

# --------------------------------------------------------------------------
# Tolerances and limits.  Module-level and read at call time (2026-09-30
# mutation-arrival rule): a mutation applied to these must reach the code.
# --------------------------------------------------------------------------
NUMERIC_RTOL = 1e-9       # relative agreement required to call a diff numeric
NUMERIC_ATOL = 1e-30      # absolute floor, for tokens that are legitimately 0
SUITE_TIMEOUT_S = 900     # per-suite wall clock in the pristine clone
GIT_POINTER_RELPATH = ".git/HEAD"   # a FILE path, not a revision -- see above

# THE SUITE SUMMARY LINE.  ONE CONSTANT, TWO READERS -- added 2026-11-06.
# Section 5's exit-code consistency check carried this regex inline, and the
# completion test added today needs the identical predicate.  Two spellings of
# one predicate have now cost this repository two sessions, so this is a
# constant from the start rather than after they drift.
SUMMARY_RE = re.compile(r"TOTAL:?\s+(\d+)\s+passed,\s+(\d+)\s+failed")

# Classes that mean THE RUN DID NOT HAPPEN, as distinct from the run happening
# and its transcript being wrong.  Named once so the census bucket and the
# standing check cannot disagree about what counts as a failure to complete --
# which is exactly how 2026-11-05's crash came to be green in one place and red
# in another.
COMPLETION_FAULTS = ("ERROR", "NO_FILE", "MISSING", "INCOMPLETE")

# VERBATIM stdout of graphene_mixed_capture_repair_audit.py in a pristine clone
# at 94b019c with scipy absent -- the 2026-11-05 census row, reproduced
# 2026-11-06.  Four banner lines; FIVE by split("\n"), which is the census's
# "line count 153 -> 5".  Kept as a literal so the reproduction cannot decay
# into a claim about a run nobody can repeat.
CRASH_STDOUT_NO_SCIPY = '==============================================================================\nMIXED-CAPTURE REPAIR AUDIT -- 2026-10-27\nthe 2026-10-26 ranked-first item, measured against the pinned pre-repair blobs\n==============================================================================\n'

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
    # ---------------------------------------------------------------------
    # Added 2026-11-04: the standing "run_suite SHOULD READ THE OBJECT STORE"
    # item, worked -- and these rows are the whole point of it.
    #
    # Every row below is SELF_WRITING: the module opens its own transcript
    # for writing, so running it inside the pristine clone overwrites the
    # file the comparison used to be taken against.  Under the old
    # working-tree read those four rows were UNREGISTRABLE, not merely
    # unregistered: registering one would have compared the suite's stdout
    # against a file that same stdout had just written, which is IDENTICAL by
    # construction and carries no information.  Section 2b measures that
    # directly rather than asserting it.
    #
    # The standing list has carried "one change moves SIX of nine
    # SELF_WRITING rows" since 2026-10-28.  THAT ESTIMATE IS WRONG, and
    # today's measurement is the correction: only FOUR of the nine can move.
    # `graphene_diffusion_mutation.py` ( > 80 s, timed out at 80 s) and
    # `graphene_perfect_contact_counterfactual.py` (119.7 s) are over the
    # call-cap budget as surely as the three rows 10-27 already measured, so
    # they acquire a SECOND reason today and stay excluded.  Six was an
    # estimate made by reading the exclusion table; four is a number from a
    # stopwatch, and the gap is the same shape as 10-27's "18 rows against
    # roughly 40" against a real 49 / 17.
    #
    # Times are measured at this commit in a pristine clone, which is the
    # budget argument: Section 2 cost 96.2 s before today.
    ("graphene_default_flip.py",             "default_flip_output.txt"),               # 1.3s
    ("graphene_contact_argument_blast_radius.py",
     "contact_argument_blast_radius_output.txt"),                                      # 3.1s
    ("graphene_import_binding_sweep.py",     "import_binding_sweep_output.txt"),       # 9.3s
    ("graphene_diffusion_current_model.py",  "diffusion_current_output.txt"),          # 18.6s
    # ---------------------------------------------------------------------
    # Added 2026-11-04: four of the five transcripts Section 5 has been
    # reporting UNACCOUNTED.  These are not SELF_WRITING and never were --
    # they were simply never registered, which is the "REGISTER THE NEW
    # TRANSCRIPTS" item that has been carried since 10-28 and has widened
    # every session since.  Each was measured at this commit and each
    # reproduces its committed transcript byte for byte, so each joins as a
    # green row rather than a red one.
    ("graphene_impurity_limited_crossover.py",
     "impurity_limited_crossover_output.txt"),                                         # 1.2s
    ("graphene_kappa_impurity_bound.py",     "kappa_impurity_bound_output.txt"),       # 1.2s
    ("graphene_branch_asymmetric_mobility.py",
     "branch_asymmetric_mobility_output.txt"),                                         # 1.3s
    ("graphene_density_dependent_mobility.py",
     "density_dependent_mobility_output.txt"),                                         # 1.5s
]

# ---------------------------------------------------------------------------
# THE CALL-CAP BUDGET, as a checked number rather than a remembered one
# ---------------------------------------------------------------------------
# A shell call on this VM is capped at 180 s.  Section 2 runs every suite in
# SUITES serially, so the table above is a spending decision and the cap is
# the budget.  Two facts were being carried in prose and are now checked in
# Section 5: that the registered table still fits, and that each module
# excluded FOR being too expensive really is too expensive.
CALL_CAP_S = 180.0

# ---------------------------------------------------------------------------
# THE BUDGET CHECK COULD NOT BE READ IN THE ONE CASE IT EXISTS FOR (2026-11-07)
# ---------------------------------------------------------------------------
# The end-of-module budget check added 2026-11-04 asks "does Section 2 still
# fit inside one shell call?".  Today it was asked for the first time in the
# negative and said nothing, for two independent reasons.
#
#   (1) IT IS PRINTED LAST.  It sits after Sections 2b, 3, 4 and 5, so it is
#       reached only by a run that FINISHED -- and a run finishes only if the
#       whole module fits in the call.  The single condition the check exists
#       to detect is therefore also the condition that guarantees the check is
#       never printed.  Measured today: the module was killed at 179 s, the
#       largest timeout this harness permits, inside Section 3; the budget
#       block is four sections further on and produced no output at all.
#
#   (2) THE QUANTITY IS A PART AND THE CAP IS THE WHOLE.  `spent` is the sum
#       of run_suite's per-suite timers only.  It excludes the clone build,
#       Sections 1, 2b, 3, 4 and 5, and the loop's own per-suite overhead.
#       Measured today: Section 2 spends 176.4 s over 43 suites (169.5 s for
#       the 40 the capped run reported, plus 2.2 + 2.2 + 2.5 s for
#       kappa_impurity_bound, branch_asymmetric_mobility and
#       density_dependent_mobility, which no capped run has ever reached).
#       Against CALL_CAP_S the old check reports 3.6 s of HEADROOM, and the
#       module does not complete in 179 s.  The reported headroom was not
#       headroom for anything.
#
# Both repairs are below and they are the same repair twice: measure the
# quantity the cap is about, and evaluate it where it can still be printed.
#
#   - `section2_crossing` is pure, so its answers can be planted exactly.
#   - The Section 2 loop evaluates it incrementally and prints a flushed
#     warning naming the suite that crosses the allowance.  A run that will
#     not finish now says so while it still can.
#   - The end-of-module verdict is gated on ELAPSED wall clock, not on
#     `spent`.  `spent` is kept and printed as a REPORT, following the
#     2026-11-06 repair of V3b and P9a: gate on an invariant, print the
#     contingent value.  The invariant is containment -- elapsed includes
#     spent, so `elapsed >= spent` must hold on every run, exactly.
#
# SECTION2_ALLOWANCE_S leaves room for the rest of the module rather than
# spending the whole call on Section 2.  The reserve is measured, not
# declared: with SUITES truncated to one suite the whole module runs in 0.7 s
# of which Section 2 is 0.2 s, so Sections 1/2b/3/4/5 and the clone build
# cost ~0.5 s on a one-suite result set.  That is NOT the reserve to use: the
# same sections cost more on a 43-entry result set, and the 179 s run died
# inside Section 3 with Section 2 complete, which bounds the tail BELOW at
# ~2 s rather than above.  The allowance is therefore set well back from the
# cap, and the honest statement is that the tail's cost as a function of
# result-set size has not been measured -- it is on the standing list.
SECTION2_ALLOWANCE_S = 150.0


def section2_crossing(durations, allowance):
    """Index of the first suite at which the CUMULATIVE spend exceeds
    `allowance`, or None if it never does.

    Pure, and deliberately strict (`>`), so that a cumulative total landing
    exactly on the allowance is NOT a crossing.  That boundary is the one
    worth planting controls on and it is exact in binary for the values the
    controls use.
    """
    total = 0.0
    for i, dt in enumerate(durations):
        total += dt
        if total > allowance:
            return i
    return None
# ---------------------------------------------------------------------------
# THE SPLIT (2026-11-08).  SECTION 2 NO LONGER HAS TO FIT IN ONE CALL
# ---------------------------------------------------------------------------
# 2026-11-07 measured Section 2 at 176.4 s over 43 suites against a 180 s cap
# and recorded the honest statement that the repair that day fixed the
# REPORTING and not the problem: the module cannot complete, and the next
# registered suite of any size makes it unrunnable rather than slow.  This is
# the problem.
#
# The shape is a two-phase run with the per-suite verdicts persisted between
# calls:
#
#     python3 graphene_pristine_transcript_audit.py --part 1/4   (and 2,3,4)
#     python3 graphene_pristine_transcript_audit.py --assemble
#
# WHAT IS DELIBERATELY NOT CHANGED.  A plain invocation with no arguments runs
# exactly as before, down to the bytes of its transcript.  The split adds two
# modes; it does not replace the single-call one.  That matters because the
# committed transcript is a claim about the single-call path and re-pointing it
# at a new path would withdraw that claim silently.
#
# WHY A PART WRITES THE POST-RUN TRANSCRIPT BYTES AND NOT JUST A VERDICT.
# Sections 2b and 3 compare, for every SELF_WRITING suite, the bytes its own
# run left in the clone against the committed blob.  Those bytes exist only
# inside the call that ran the suite -- the clone is a tempdir and is deleted
# with it.  An assemble phase that rebuilt a fresh clone would find the
# committed bytes there instead and E2/E3 would quietly become vacuous: every
# self-writing suite would read "rewrote it byte-identically" because nothing
# had rewritten anything.  So a part carries them, and --assemble restores
# them into its own clone before Section 2b runs.  THIS IS THE ONE PLACE THE
# SPLIT COULD HAVE SILENTLY WEAKENED A CHECK, and it is why the equivalence
# test beside this module compares the two paths' output rather than trusting
# that they agree.
#
# WHY THE PARTITION IS BY COUNT AND NOT BY COST.  Cost-packing needs a
# per-suite cost table, and the only one that exists (MEASURED_OVER_CAP_S) has
# two rows out of 43.  A partition built on an estimate of the thing it is
# supposed to be bounding is the shape of fault the 11-07 budget check was.
# `block_slice` is therefore arithmetic on the table's LENGTH, which is known
# exactly, and --assemble emits the measured per-suite table so that a session
# with the measurement in hand can pack by cost and say what it is packing
# against.  That is the next item, not this one.
CACHE_ENV = "GRAPHENE_CENSUS_CACHE"
ONLY_ENV = "GRAPHENE_CENSUS_ONLY"


def cache_dir():
    """Where part files live.  Outside the repository by default: a part file
    is a measurement of one call and committing it would make the census's
    own intermediate state part of the tree it censuses."""
    d = os.environ.get(CACHE_ENV) or os.path.join(
        tempfile.gettempdir(), "graphene_census_parts")
    return d


def suite_table():
    """SUITES, optionally filtered by GRAPHENE_CENSUS_ONLY.

    The filter exists for the equivalence test, which has to run the whole
    module twice and cannot afford the real table.  It is read here rather
    than in main() so that every mode -- single-call, part and assemble --
    sees the same table, which is what makes a two-path comparison on a
    truncated table a statement about the paths and not about the filter.
    """
    only = os.environ.get(ONLY_ENV, "").strip()
    if not only:
        return list(SUITES)
    want = [s.strip() for s in only.split(",") if s.strip()]
    return [(m, t) for m, t in SUITES if m in want]


def block_slice(total, n_parts, index):
    """The `index`-th of `n_parts` contiguous blocks of range(total).

    Pure, and arithmetic on `total` alone -- no cost estimate anywhere in it.
    The first `total % n_parts` blocks are one longer than the rest, which is
    the only distribution that is both contiguous and balanced to within one.

    THE PROPERTY WORTH PLANTING IS NOT THE SIZES.  It is that the blocks
    PARTITION range(total) exactly: concatenating them in order must give
    [0, 1, ..., total-1], with nothing dropped and nothing twice.  A census
    that silently drops a suite would report a smaller STALE count and read
    as an improvement, so the partition is the invariant and the sizes are a
    consequence.  See C11.
    """
    if n_parts < 1:
        raise ValueError("n_parts must be >= 1, got %r" % (n_parts,))
    if not 0 <= index < n_parts:
        raise ValueError("index %r out of range for %r parts" % (index, n_parts))
    base, extra = divmod(total, n_parts)
    start = index * base + min(index, extra)
    size = base + (1 if index < extra else 0)
    return list(range(start, start + size))


def partition_of(total, n_parts):
    """Every block, in order.  The object C11 asserts on."""
    return [block_slice(total, n_parts, i) for i in range(n_parts)]


def flatten(blocks):
    out = []
    for b in blocks:
        out.extend(b)
    return out


def part_path(n_parts, index):
    return os.path.join(cache_dir(),
                        "part-%d-of-%d.json" % (index + 1, n_parts))


def table_fingerprint(table):
    """A content-addressed id for the suite table a part was run against.

    Parts are only comparable if they were run against the SAME table at the
    SAME commit.  Both go in the part file and --assemble refuses a mismatch
    by exact string equality rather than by a warning -- an assembled census
    built from two different tables would report a coverage figure that
    describes neither.  See C12.
    """
    h = hashlib.sha256()
    for m, t in table:
        # FIXED BEFORE ITS FIRST USE, BY ITS OWN CONTROL (2026-11-08).  The
        # first form was `t or ""`, which maps BOTH None (a suite with no
        # committed transcript) and "" (a suite registered as having an empty
        # one) to the same byte string, so a registry change between those
        # two states left the fingerprint unmoved.  C12.4 was written to
        # assert they differ and it FAILED on the first run of this module --
        # which is the whole argument for planting a control on each claim a
        # comment makes rather than on the claims that seem doubtful.
        h.update(("%s\x00%s\x1f"
                  % (m, "\x01<unregistered>" if t is None else t)
                  ).encode("utf-8"))
    return h.hexdigest()[:16]


def write_part(n_parts, index, commit, table, results, cached, exits,
               post_bytes, spent, elapsed):
    os.makedirs(cache_dir(), exist_ok=True)
    payload = {
        "schema": 1,
        "n_parts": n_parts,
        "index": index,
        "commit": commit,
        "table_fingerprint": table_fingerprint(table),
        "modules": [m for m, _ in table],
        "spent_s": spent,
        "elapsed_s": elapsed,
        "results": {m: {"verdict": v, "dt": dt, "evidence": ev}
                    for m, (v, dt, ev) in results.items()},
        "cached": cached,
        "exits": exits,
        # base64 because these are bytes and a transcript may legitimately
        # contain anything; the census must not re-encode what it compares.
        "post_bytes": {t: base64.b64encode(b).decode("ascii")
                       for t, b in post_bytes.items() if b is not None},
    }
    with open(part_path(n_parts, index), "w", encoding="utf-8") as fh:
        json.dump(payload, fh)
    return part_path(n_parts, index)


def read_parts(n_parts, commit, table):
    """Every part file for an `n_parts` split, or an explanation of what is
    missing.  Returns (parts, problems)."""
    parts, problems = [], []
    want_fp = table_fingerprint(table)
    for i in range(n_parts):
        p = part_path(n_parts, i)
        if not os.path.exists(p):
            problems.append("part %d of %d has not been run (%s is absent)"
                            % (i + 1, n_parts, p))
            continue
        try:
            with open(p, encoding="utf-8") as fh:
                d = json.load(fh)
        except Exception as exc:                      # pragma: no cover
            problems.append("part %d of %d is unreadable: %s" % (i + 1, n_parts, exc))
            continue
        if d.get("commit") != commit:
            problems.append("part %d of %d was run at commit %s, not %s"
                            % (i + 1, n_parts, d.get("commit"), commit))
            continue
        if d.get("table_fingerprint") != want_fp:
            problems.append("part %d of %d was run against suite table %s, "
                            "not %s" % (i + 1, n_parts,
                                        d.get("table_fingerprint"), want_fp))
            continue
        parts.append(d)
    return parts, problems


def assemble(parts, table):
    """Merge part payloads into the dicts Section 2 would have built.

    Returns (results, cached, exits, post_bytes, spent, problems).  A module
    present in two parts is a problem, not a last-writer-wins merge: it would
    mean the partition is not one.
    """
    results, cached, exits, post_bytes, problems = {}, {}, {}, {}, []
    spent = 0.0
    for d in sorted(parts, key=lambda x: x["index"]):
        for m, r in d["results"].items():
            if m in results:
                problems.append("%s appears in more than one part" % m)
                continue
            results[m] = (r["verdict"], r["dt"], r["evidence"])
        for m, out in d["cached"].items():
            cached[m] = out
        for m, rc in d["exits"].items():
            exits[m] = rc
        for t, b64 in d["post_bytes"].items():
            post_bytes[t] = base64.b64decode(b64)
        spent += d["spent_s"]
    missing = [m for m, _ in table if m not in results]
    if missing:
        problems.append("no part ran: %s" % ", ".join(missing))
    return results, cached, exits, post_bytes, spent, problems


def parse_mode(argv):
    """(mode, n_parts, index).  mode is 'single', 'part' or 'assemble'.

    `--part i/n` is 1-based for the operator and 0-based internally, because
    an operator typing `--part 0/4` has made an error this should name rather
    than silently run part 1.
    """
    mode, n_parts, index = "single", None, None
    rest = [a for a in argv[1:] if a not in ("--",)]
    i = 0
    while i < len(rest):
        a = rest[i]
        if a == "--assemble":
            mode = "assemble"
        elif a == "--part":
            i += 1
            if i >= len(rest):
                raise SystemExit("--part needs an argument like 2/4")
            mode, n_parts, index = _parse_part_spec(rest[i])
        elif a.startswith("--part="):
            mode, n_parts, index = _parse_part_spec(a.split("=", 1)[1])
        elif a.startswith("--parts="):
            n_parts = int(a.split("=", 1)[1])
        else:
            raise SystemExit("unrecognised argument %r" % a)
        i += 1
    return mode, n_parts, index


def _parse_part_spec(spec):
    if "/" not in spec:
        raise SystemExit("--part wants i/n, got %r" % spec)
    a, b = spec.split("/", 1)
    one_based, n = int(a), int(b)
    if n < 1:
        raise SystemExit("--part i/n needs n >= 1, got %r" % spec)
    if not 1 <= one_based <= n:
        raise SystemExit("--part %s is outside 1..%d" % (spec, n))
    return "part", n, one_based - 1


def _raises_systemexit(argv):
    """True iff parse_mode(argv) refuses.  A named helper so C13.2 reads as a
    requirement rather than as a try/except buried in a check call."""
    try:
        parse_mode(argv)
    except SystemExit:
        return True
    return False

# Measured 2026-11-04 in a pristine clone at this commit.  Only modules whose
# cost was measured today appear here; the three rows 2026-10-27 measured over
# the cap are named in their exclusion comments and not re-timed.
MEASURED_OVER_CAP_S = {
    "graphene_perfect_contact_counterfactual.py": 119.7,
    "graphene_diffusion_mutation.py": 80.0,   # a LOWER bound: timed out at 80 s
}

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
#   DRIFTS_BY_MEASUREMENT
#                  DECLARED, added 2026-11-04.  The module's stdout contains
#                  a measurement of the repository that every session moves,
#                  so a byte comparison is red by construction and carries no
#                  finding.  Distinct from UNMEASURED: the remedy is a change
#                  to the MODULE, not a row in SUITES.
#   OVER_CAP       CHECKED, as of 2026-11-04.  The module is SELF_WRITING
#                  too, but that reason stopped being load-bearing today:
#                  the object-store read removes it.  What remains is cost.
#                  Section 5 checks the measurement against CALL_CAP_S
#                  rather than taking the word of this comment.
TRANSCRIPT_EXCLUSIONS = {
    # AMENDED 2026-11-04.  Four rows left this table for SUITES -- see the
    # block at the end of SUITES.  The two below stayed, and their reason
    # CHANGED: SELF_WRITING no longer excuses anything now that run_suite
    # reads the object store, so each was timed, and each is over budget.
    "diffusion_mutation_output.txt":
        ("OVER_CAP", "graphene_diffusion_mutation.py"),
    "perfect_contact_counterfactual_output.txt":
        ("OVER_CAP", "graphene_perfect_contact_counterfactual.py"),
    # The three below are ALSO over the 180 s call cap (2026-10-27 section 5
    # item 4, measured there), so they carry two independent reasons.  They
    # keep the SELF_WRITING class because that is the reason 10-27 recorded
    # and their cost was measured elsewhere; the OVER_CAP rows above are the
    # two measured today.
    "perfect_contact_mutation_output.txt":
        ("SELF_WRITING", "graphene_perfect_contact_mutation.py"),
    "velocity_saturation_output.txt":
        ("SELF_WRITING", "graphene_velocity_saturation_model.py"),
    "velocity_saturation_mutation_output.txt":
        ("SELF_WRITING", "graphene_velocity_saturation_mutation.py"),
    "pristine_transcript_audit_output.txt":
        ("SELF", "graphene_pristine_transcript_audit.py"),
    # Added 2026-11-08 with the split.  The equivalence test RUNS THIS MODULE
    # four times as a subprocess, so registering it in SUITES would have the
    # census execute, inside its own pristine clone, a module that executes
    # the census -- the self-certification 09-30 ruled out, nested one level
    # deeper and costing ~30 s a round.  SELF is the right class for exactly
    # that reason and not for cost: the exclusion would stand even if it were
    # free.
    "census_split_equivalence_output.txt":
        ("SELF", "graphene_census_split_equivalence.py"),
    "sensitivity_converged_output.txt":
        ("NOT_A_STDOUT", "graphene_crossover_sensitivity_model.py"),
    "capture_surface_mutation_output.txt":
        ("UNMEASURED", "graphene_capture_surface_mutation.py"),
    "installed_surface_mutation_output.txt":
        ("UNMEASURED", "graphene_installed_surface_mutation.py"),
    "value_read_mutation_output.txt":
        ("UNMEASURED", "graphene_value_read_mutation.py"),
    # Added 2026-11-04.  The fifth UNACCOUNTED transcript, and the only one
    # of the five that could not simply be registered.  It is not stale by
    # anybody's mistake: `graphene_chapter_status_audit.py` prints the PROSE
    # WORD COUNT of every chapter, so its transcript goes stale the moment
    # any session writes a sentence -- today it reads 37686 words for
    # Chapter 4 against a tree holding 38305.  Registering it would make
    # Section 2 permanently red for a drift that carries no finding, and a
    # permanently red check is one nobody reads.  DRIFTS_BY_MEASUREMENT is
    # the honest name for that, and unlike UNMEASURED it is not a promise to
    # do the work later: the work is to give the module a stable projection
    # (a section/threshold verdict without the raw counts), which is a change
    # to that module and not a row here.
    "chapter_status_audit_output.txt":
        ("DRIFTS_BY_MEASUREMENT", "graphene_chapter_status_audit.py"),
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

# Set by __main__ from the command line.  Module-level rather than threaded
# through main() because Section 2's loop, the part writer and the budget
# block all need it and main() is one function.
MODE = "single"
N_PARTS = 1
PART_INDEX = 0

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


def read_committed_blob(repo, commit, relpath):
    """The committed BYTES of `relpath` at `commit`, or None if absent.

    ADDED 2026-11-04, and this is the whole of the repair.  Until today the
    comparison in Section 2 read the transcript out of the pristine clone's
    WORKING TREE, which is sound only for a module that does not write its
    own transcript.  For one that does, the suite's own run overwrote the
    file a moment before it was read, so the comparison was
    `stdout == stdout` and reported IDENTICAL whatever the committed blob
    said.  Nine transcripts were excluded from Section 2 for exactly that
    reason, and the exclusion was correct given the read.

    The object store is immutable and content-addressed, so reading the blob
    removes the circularity at the root instead of sequencing around it.
    Deliberately NOT `git()`: that helper `.strip()`s its output, and a
    byte-for-byte comparison cannot be built on a stripped read.  Bytes, not
    text, for the same reason.
    """
    p = subprocess.run(["git", "show", "%s:%s" % (commit, relpath)],
                       cwd=repo, capture_output=True)
    if p.returncode != 0:
        return None
    return p.stdout


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


def has_summary(text):
    """Does `text` carry a suite's own end-of-run summary line?"""
    return bool(text) and bool(SUMMARY_RE.search(text))


def completion_verdict(stdout, committed_text):
    """COMPLETE / INCOMPLETE / ERROR -- did the RUN happen?

    Pure, so the controls in Section 3 can drive every branch without running
    a suite.  Deliberately separate from `classify`: that function answers
    "is the committed transcript still true", and conflating the two is the
    2026-11-05 defect this amendment repairs.

      ERROR       no stdout at all -- the historical test, unchanged.
      INCOMPLETE  the committed transcript announces a summary line and this
                  run's stdout has none.  The suite is known to print one, so
                  its absence means the run did not reach the end.
      COMPLETE    otherwise -- including a suite whose transcript never had a
                  summary line, which makes no completion claim and must be
                  compared exactly as before.
    """
    if not stdout:
        return "ERROR"
    if has_summary(committed_text) and not has_summary(stdout):
        return "INCOMPLETE"
    return "COMPLETE"


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
    # The real quantity the call cap is about.  Started before anything,
    # including the clone build, because the cap applies to all of it.
    t_start = time.monotonic()
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
        # Both reads of every registered transcript, taken BEFORE any suite
        # runs.  The working-tree read is kept only so that Section 2b can
        # measure the two against each other; nothing compares against it.
        pre_tree, pre_blob = {}, {}
        for _m, _t in SUITES:
            if not _t:
                continue
            _p = os.path.join(clone, _t)
            pre_tree[_t] = open(_p, "rb").read() if os.path.exists(_p) else None
            pre_blob[_t] = read_committed_blob(clone, commit, _t)
        # ------------------------------------------------------------------
        # THE SPLIT, at the one point in main() where it matters (2026-11-08).
        # In 'part' mode the loop below runs a contiguous block of the table
        # and nothing else; in 'assemble' mode it runs NOTHING and the dicts
        # come off disk, with the post-run transcript bytes restored into this
        # call's clone so that Sections 2b and 3 see what the running call saw.
        # ------------------------------------------------------------------
        _table = suite_table()
        _post_bytes = {}
        _assembled_parts = []
        if MODE == "part":
            _idx = block_slice(len(_table), N_PARTS, PART_INDEX)
            _s2_table = [_table[i] for i in _idx]
            print("  PART %d of %d -- suites %s of %d, by block_slice"
                  % (PART_INDEX + 1, N_PARTS,
                     ("%d..%d" % (_idx[0], _idx[-1])) if _idx else "(none)",
                     len(_table)))
        elif MODE == "assemble":
            _s2_table = []
            _assembled_parts, _probs = read_parts(N_PARTS, commit, _table)
            if _probs:
                print("  CANNOT ASSEMBLE:")
                for p in _probs:
                    print("    " + p)
                print("\n  Run the missing parts first:")
                for i in range(N_PARTS):
                    print("    python3 %s --part %d/%d" % (SELF, i + 1, N_PARTS))
                return 2
            (results, cached, exits, _post_bytes, _assembled_spent,
             _probs) = assemble(_assembled_parts, _table)
            if _probs:
                print("  CANNOT ASSEMBLE:")
                for p in _probs:
                    print("    " + p)
                return 2
            # Restore the bytes each part's own run left behind.  Without this
            # Sections 2b/3 would compare the committed blob against itself.
            for _t, _b in _post_bytes.items():
                _p = os.path.join(clone, _t)
                if os.path.exists(os.path.dirname(_p)):
                    with open(_p, "wb") as _fh:
                        _fh.write(_b)
            print("  ASSEMBLED from %d parts at commit %s"
                  % (len(_assembled_parts), commit))
            for _d in sorted(_assembled_parts, key=lambda x: x["index"]):
                print("    part %d of %d: %2d suites, suite time %6.1f s, "
                      "elapsed %6.1f s"
                      % (_d["index"] + 1, _d["n_parts"], len(_d["results"]),
                         _d["spent_s"], _d["elapsed_s"]))
            for module, _tr in _table:
                _v, _dt, _ev = results[module]
                print("  %-44s %-14s %6.1fs" % (module, _v, _dt))
                for _e in _ev[:6]:
                    print("        " + _e)
        else:
            _s2_table = list(_table)
        # The incremental budget, evaluated WHERE IT CAN STILL BE PRINTED.
        # `warned` makes the warning fire once, on the suite that crosses.
        s2_durations = []
        s2_warned = False
        for module, transcript in _s2_table:
            if not os.path.exists(os.path.join(clone, module)):
                results[module] = ("MISSING", 0.0, ["not present at this commit"])
                print("  %-44s MISSING" % module)
                continue
            out, dt, err, rc = run_suite(clone, module)
            s2_durations.append(dt)
            if not s2_warned and section2_crossing(s2_durations,
                                                   SECTION2_ALLOWANCE_S) is not None:
                s2_warned = True
                print("  !! SECTION 2 BUDGET WARNING -- crossed its %.1f s "
                      "allowance at %s" % (SECTION2_ALLOWANCE_S, module))
                print("     cumulative suite time %.1f s, elapsed wall clock "
                      "%.1f s of a %.1f s call cap."
                      % (sum(s2_durations), time.monotonic() - t_start,
                         CALL_CAP_S))
                print("     This run will probably not reach the end-of-module "
                      "budget block, which is why this line is printed here "
                      "and not there (2026-11-07).")
                sys.stdout.flush()
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
            # THE COMMITTED TRANSCRIPT, FROM THE OBJECT STORE.  Not from the
            # clone's working tree: see read_committed_blob.  The bytes were
            # taken before any suite ran, but that is belt-and-braces -- a
            # blob at a fixed commit cannot change under us, which is the
            # property the old read lacked.
            blob = pre_blob.get(transcript)
            if blob is None:
                results[module] = ("NO_FILE", dt,
                                   [transcript + " is not committed at this id"])
                print("  %-44s NO_FILE  %6.1fs" % (module, dt))
                continue
            committed_text = blob.decode("utf-8", "replace")
            # DID THE RUN HAPPEN, before asking whether the transcript is
            # true.  See the 2026-11-06 amendment: the answers are about
            # different objects and the old code could only give the second.
            if completion_verdict(out, committed_text) == "INCOMPLETE":
                ev = ("exit %d; the committed transcript ends in a summary "
                      "line and this run's stdout has none, so the suite did "
                      "not reach its own end -- this is a failed RUN, not a "
                      "stale transcript" % rc)
                results[module] = ("INCOMPLETE", dt, [ev])
                print("  %-44s %-14s %6.1fs  exit %d (did NOT reach its own "
                      "summary line)" % (module, "INCOMPLETE", dt, rc))
                print("        " + ev)
                continue
            verdict, evidence = classify(committed_text, out)
            results[module] = (verdict, dt, evidence)
            tag = "" if rc == 0 else "  exit %d (guard red, transcript still checked)" % rc
            print("  %-44s %-14s %6.1fs%s" % (module, verdict, dt, tag))
            for e in evidence[:6]:
                print("        " + e)


        # ------------------------------------------------------------------
        # A PART STOPS HERE (2026-11-08).  Sections 2b-5 need the WHOLE
        # result set -- E1 counts transcripts, Section 5 checks coverage of
        # the registry -- so running them on a block would not be a partial
        # answer but a wrong one.  The part writes what it measured and says
        # what to run next.
        # ------------------------------------------------------------------
        if MODE == "part":
            _post = {}
            for _m, _t in _s2_table:
                if not _t:
                    continue
                _p = os.path.join(clone, _t)
                _post[_t] = open(_p, "rb").read() if os.path.exists(_p) else None
            _spent = sum(dt for _, dt, _ in results.values())
            _elapsed = time.monotonic() - t_start
            _written = write_part(N_PARTS, PART_INDEX, commit, _table,
                                  results, cached, exits, _post,
                                  _spent, _elapsed)
            print()
            print("  PART %d of %d BUDGET" % (PART_INDEX + 1, N_PARTS))
            print("    suite time this part     %6.1f s  over %d suites"
                  % (_spent, len(_s2_table)))
            print("    ELAPSED wall clock       %6.1f s" % _elapsed)
            print("    one-shell-call cap       %6.1f s" % CALL_CAP_S)
            print("    headroom, by elapsed     %6.1f s" % (CALL_CAP_S - _elapsed))
            # The same containment invariant as the single-call path, and it
            # is still exact: this call's elapsed clock contains the suite
            # timers this call charged.  It is asserted PER CALL because that
            # is the only scope in which it is true -- the assembled total
            # spans several calls and cannot be compared against one cap,
            # which is the thing --assemble has to say differently.
            check("part %d/%d: elapsed wall clock contains its suite time"
                  % (PART_INDEX + 1, N_PARTS),
                  _elapsed >= _spent,
                  "elapsed %.1f s, suite time %.1f s, difference %.1f s"
                  % (_elapsed, _spent, _elapsed - _spent))
            check("part %d/%d fits inside one shell call"
                  % (PART_INDEX + 1, N_PARTS),
                  _elapsed < CALL_CAP_S,
                  "%.1f s elapsed of a %.1f s cap" % (_elapsed, CALL_CAP_S))
            print()
            print("  written: %s" % _written)
            _left = [i for i in range(N_PARTS)
                     if not os.path.exists(part_path(N_PARTS, i))]
            if _left:
                print("  still to run: %s"
                      % ", ".join("--part %d/%d" % (i + 1, N_PARTS) for i in _left))
            else:
                print("  every part is on disk; next: python3 %s --assemble "
                      "--parts=%d" % (SELF, N_PARTS))
            print("\n" + "=" * 74)
            print("PART SUMMARY")
            print("=" * 74)
            for lab in _PASS:
                print("  PASS  %s" % lab)
            for lab in _FAIL:
                print("  FAIL  %s" % lab)
            print("\nTOTAL: %d passed, %d failed" % (len(_PASS), len(_FAIL)))
            return 1 if _FAIL else 0

        # ------------------------------------------------------------------
        print("\n" + "-" * 74)
        print("SECTION 2b -- the two reads, measured against each other")
        print("Created 2026-11-04.  The object-store read replaces a")
        print("working-tree read, so it owes two things: that it agrees with")
        print("the old read everywhere the old read was sound, and that it")
        print("disagrees somewhere -- or it changed nothing.")
        print("-" * 74)

        # E1 -- THE OVERLAPPING LIMIT, AND IT IS AN EXACT ONE.
        # For a module that does not write its own transcript, nothing in the
        # pristine clone can have touched the file, so the working-tree bytes
        # and the object-store bytes must be EQUAL -- not close, equal.  This
        # is the check that the new read is a read of the same thing.
        e1_bad, e1_n = [], 0
        for module, transcript in _table:
            if not transcript or transcript not in pre_tree:
                continue
            if writes_own_transcript(clone, module, transcript):
                continue
            t, b = pre_tree[transcript], pre_blob[transcript]
            if t is None or b is None:
                e1_bad.append("%s: working tree %s, object store %s"
                              % (transcript,
                                 "absent" if t is None else "present",
                                 "absent" if b is None else "present"))
            elif t != b:
                e1_bad.append("%s: %d working-tree bytes vs %d committed bytes"
                              % (transcript, len(t), len(b)))
            else:
                e1_n += 1
        check("E1 for every suite that does NOT write its own transcript, the "
              "object-store read is byte-identical to the old working-tree read",
              not e1_bad,
              "\n".join("  " + b for b in e1_bad) if e1_bad else
              "%d transcripts, 0 differing bytes.  An EXACT agreement, not a "
              "tolerance: a pristine checkout of a commit and that commit's "
              "blobs are the same bytes, so any difference here would mean "
              "the clone is not pristine." % e1_n)

        # E2 -- A CONTROL, NOT A REQUIREMENT.  RECAST 2026-11-05.
        #
        # SUPERSEDED FRAMING, LEFT IN PLACE SO THE CHANGE IS LEGIBLE.  Created
        # 2026-11-04 as an unconditional requirement, on this reasoning: "for a
        # registered SELF_WRITING suite the file in the clone has been
        # overwritten by the suite's own run, so the old read and the new read
        # now differ.  If none of them differ, the repair is cosmetic and this
        # check says so."  It passed on 11-04 and was counted among that day's
        # 22 passes.  On 11-05 it FAILED, with the evidence line "every
        # registered self-writing suite rewrote its transcript
        # byte-identically" -- it failed because the tree had stopped drifting.
        # A check that can only pass while something is wrong is a control that
        # requires a defect to exist, which is the fault this repository has
        # now recorded in several other places.
        #
        # The circularity the object-store read removes is STRUCTURAL: the old
        # read compared a transcript against a file the suite's own run had
        # just written, whether or not those bytes happened to differ that day.
        # E1 is what shows the new read is a read of the same thing; E3 is what
        # shows that a disagreement, where there is one, is adjudicated in the
        # strict direction.  Neither needs the tree to be unhealthy.  So the
        # occurrence is now REPORTED and not scored (E2r, CONFIRMED /
        # NOT_EXERCISED), the possibility it guards against is established by
        # PLANTING both branches rather than by waiting for one to arrive
        # (E2c1 / E2c2), and the one genuine requirement left -- that the
        # predicate has something to run on at all -- is what E2 now checks.
        def _circular(post_bytes, committed_bytes):
            """True iff the old working-tree read would have compared a
            suite's output against bytes its own run had changed."""
            return post_bytes != committed_bytes

        # The two planted branches.  EXACT, both of them: one must be True and
        # the other must be False, with no tolerance and no dependence on the
        # state of the repository.  That independence is the whole point -- it
        # is what the 11-04 form of E2 lacked.
        check("E2c1 the circularity predicate is TRUE when the run rewrote the "
              "transcript differently (planted)",
              _circular(b"line one\nline two CHANGED\n",
                        b"line one\nline two\n"),
              "planted post-run bytes differ from the planted committed blob by "
              "one word, so a working-tree read would have been circular")
        check("E2c2 and FALSE when the run rewrote it byte-identically "
              "(planted)",
              not _circular(b"line one\nline two\n", b"line one\nline two\n"),
              "planted post-run bytes equal the planted committed blob, so the "
              "old read would have seen exactly what the new read sees.  This "
              "is the state the live tree was in on 2026-11-05, and it is not "
              "a fault")

        e2_same, e2_diff, e2_missing = [], [], []
        for module, transcript in _table:
            if not transcript or module not in cached:
                continue
            if not writes_own_transcript(clone, module, transcript):
                continue
            post = os.path.join(clone, transcript)
            if not os.path.exists(post):
                e2_missing.append(transcript)
                continue
            post_bytes = open(post, "rb").read()
            if _circular(post_bytes, pre_blob.get(transcript)):
                e2_diff.append(transcript)
            else:
                e2_same.append(transcript)

        # THE REQUIREMENT THAT SURVIVES.  A suite registered as self-writing
        # must actually leave the file behind, or the predicate has nothing to
        # run on and the exclusion table is wrong about that suite.  This can
        # fail for a real reason; the old E2 could not.
        check("E2 the circularity predicate is total over the registered "
              "self-writing suites", not e2_missing,
              "\n".join("  " + t + ": registered self-writing, but absent "
                         "from the clone after the run" for t in e2_missing)
              if e2_missing else
              "%d registered self-writing suite(s), every one left its "
              "transcript behind for the predicate to read"
              % (len(e2_diff) + len(e2_same)))

        # THE REPORT, WHICH IS NOT SCORED.  Either branch is a sound state of
        # the repository, so neither is a pass and neither is a failure.
        if e2_diff:
            print("  [CONFIRMED] E2r the circularity was exercised on this run")
            for t in e2_diff:
                print("        " + t + ": the old read would have compared "
                      "this against a file the run had just written")
            if e2_same:
                print("        (and %d were rewritten identically, which the "
                      "old read could not have distinguished)" % len(e2_same))
        else:
            print("  [NOT_EXERCISED] E2r no registered self-writing suite "
                  "rewrote its transcript differently on this run")
            print("        %d suite(s) rewrote byte-identically.  The "
                  "circularity is still real -- E2c1/E2c2 above plant both of "
                  "its branches -- it simply did not occur today, which is "
                  "what a healthy tree looks like." % len(e2_same))

        # E3 -- THE TWO READS, ADJUDICATED.  E2r reports whether the bytes
        # differed on this run (it said CONFIRMED on 11-04 and NOT_EXERCISED on
        # 11-05); this shows that where they do differ it CHANGES A VERDICT,
        # which is the only form of the claim that
        # matters.  Where they disagree the object-store verdict must be the
        # stricter one: a read that disagreed by being more lenient would be a
        # worse read, not a better one.
        _RANK = {"IDENTICAL": 0, "NUMERIC_DRIFT": 1, "STALE": 2}
        e3_flips, e3_wrong = [], []
        for module, transcript in _table:
            if not transcript or module not in cached:
                continue
            if not writes_own_transcript(clone, module, transcript):
                continue
            post = os.path.join(clone, transcript)
            if not os.path.exists(post):
                continue
            old_text = open(post, encoding="utf-8", errors="replace").read()
            new_text = pre_blob[transcript].decode("utf-8", "replace")
            v_old, _ = classify(old_text, cached[module])
            v_new, _ = classify(new_text, cached[module])
            if v_old == v_new:
                continue
            e3_flips.append("%s: working-tree read %s -> object-store read %s"
                            % (transcript, v_old, v_new))
            if _RANK.get(v_new, 9) <= _RANK.get(v_old, 9):
                e3_wrong.append("%s: the object-store read is the WEAKER "
                                "verdict (%s vs %s)" % (transcript, v_new, v_old))
        check("E3 where the two reads disagree, the object-store verdict is "
              "the stricter one", not e3_wrong,
              "\n".join("  " + b for b in e3_wrong) if e3_wrong else
              ("%d verdict(s) changed, every one in the strict direction:\n"
               % len(e3_flips)) + "\n".join("  " + f for f in e3_flips)
              if e3_flips else
              "no verdict changed: the self-writing suites all rewrote their "
              "transcripts to the committed bytes, so this run has no "
              "disagreement to adjudicate")

        # ------------------------------------------------------------------
        print("\n" + "-" * 74)
        print("SECTION 3 -- positive controls on the detector itself")
        print("09-30's rule: a control that has never fired is a control with")
        print("no evidence that it can.  These run against REAL captured output.")
        print("-" * 74)
        donor = None
        for module, transcript in _table:
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

        # C8 -- THE 2026-11-06 AMENDMENT, AGAINST EXACTLY KNOWN ANSWERS.
        # Every answer below is known before the code runs, and none of them
        # is read off the live tree: the crash stdout is a literal, and the
        # committed transcript it is compared against is a string written
        # here.  C8b is the load-bearing one -- it is the measured reason the
        # 11-05 census called a crash a stale transcript.
        # C8b's FIRST DRAFT FAILED, AND THE CONTROL WAS WRONG, NOT THE CODE.
        # It built this fixture with three body lines, which makes it FIVE by
        # split("\n") -- the same count as the crash stdout.  `classify` then
        # never reached its line-count branch and reported "line 1 is not a
        # numeric difference" instead of "line count N -> 5", so the control
        # went red against correct code.  The real 11-05 pair was 153 against
        # 5; a control for it has to DIFFER IN LENGTH, which is the property
        # the first draft forgot to give it.  The FAIL is recorded here rather
        # than quietly corrected, as 2026-11-05's P5 was.
        #
        # AND THE SECOND DRAFT MISCOUNTED TOO, which is why the line below is
        # a runtime assert and not a comment: ten lines each ending in "\n"
        # are ELEVEN by split("\n"), the trailing empty included.  The assert
        # caught it before the control could go red for a third reason that
        # had nothing to do with the code under test.
        c8_committed = ("HEADER\n" + "".join("body line %d\n" % i
                                             for i in range(1, 9))
                        + "TOTAL: 12 passed, 0 failed\n")
        assert len(c8_committed.split("\n")) == 11, "C8's committed fixture"
        check("C8a the reproduced crash stdout is 5 lines -- the census's "
              "'-> 5', from the literal and not from a memory",
              len(CRASH_STDOUT_NO_SCIPY.split("\n")) == 5,
              "got %d" % len(CRASH_STDOUT_NO_SCIPY.split("\n")))
        _v8, _e8 = classify(c8_committed, CRASH_STDOUT_NO_SCIPY)
        check("C8b the SUPERSEDED reading: a crash classifies STALE, with a "
              "line-count evidence line ending in '-> 5'",
              _v8 == "STALE" and bool(_e8)
              and _e8[0] == "line count 11 -> 5",
              "got %s %r -- this is what 2026-11-05 printed about a suite "
              "that had not run" % (_v8, _e8[:1]))
        check("C8c and the NEW reading calls the same pair INCOMPLETE",
              completion_verdict(CRASH_STDOUT_NO_SCIPY, c8_committed)
              == "INCOMPLETE",
              "got %s"
              % completion_verdict(CRASH_STDOUT_NO_SCIPY, c8_committed))
        check("C8d a run that DID reach its summary line is COMPLETE even "
              "when its output differs -- real staleness is not swallowed",
              completion_verdict("HEADER\nchanged\nTOTAL: 1 passed, 0 failed\n",
                                 c8_committed) == "COMPLETE",
              "a completion test that hid staleness would be a worse fault "
              "than the one it repairs")
        check("C8e a committed transcript with NO summary line makes no "
              "completion claim, so a truncated run stays COMPLETE",
              completion_verdict("HEADER\n", "HEADER\nline two\n")
              == "COMPLETE",
              "self-calibration: the expectation comes from the suite's own "
              "transcript, not from a hard-coded list of which modules print "
              "a summary")
        check("C8f empty stdout is still ERROR -- the historical branch is "
              "unchanged", completion_verdict("", c8_committed) == "ERROR",
              "got %s" % completion_verdict("", c8_committed))
        check("C8g INCOMPLETE is in COMPLETION_FAULTS, so a crash cannot be "
              "green in the completion check while red in the staleness one "
              "-- the 2026-11-05 state",
              "INCOMPLETE" in COMPLETION_FAULTS and "STALE" not in COMPLETION_FAULTS,
              "COMPLETION_FAULTS = %r" % (COMPLETION_FAULTS,))
        check("C8h the summary regex is ONE constant: Section 5 reads "
              "SUMMARY_RE and no second spelling survives in the source",
              len(re.findall(r"TOTAL:\?\\s\+", open(os.path.join(repo, SELF)).read())) == 1,
              "a second inline copy of the pattern is the drift this "
              "constant exists to prevent")

        # ------------------------------------------------------------------
        print("\n" + "-" * 74)
        print("SECTION 4 -- the standing requirement")
        print("-" * 74)
        stale = sorted(m for m, (v, _, _) in results.items() if v == "STALE")
        drift = sorted(m for m, (v, _, _) in results.items() if v == "NUMERIC_DRIFT")
        # COMPLETION_FAULTS, not a second spelling of the list.  INCOMPLETE
        # joined it 2026-11-06 and the check below picked it up with no edit,
        # which is the whole point of naming the set once.
        errored = sorted(m for m, (v, _, _) in results.items() if v in COMPLETION_FAULTS)
        incomplete = sorted(m for m, (v, _, _) in results.items() if v == "INCOMPLETE")
        ident = sorted(m for m, (v, _, _) in results.items() if v == "IDENTICAL")

        red = sorted(m for m, rc in exits.items() if rc not in (0, None))
        print("  IDENTICAL      %2d   %s" % (len(ident), ", ".join(s[:34] for s in ident)))
        print("  NUMERIC_DRIFT  %2d   %s" % (len(drift), ", ".join(drift)))
        print("  STALE          %2d   %s" % (len(stale), ", ".join(stale)))
        print("  DID NOT COMPLETE %d   %s" % (len(errored), ", ".join(errored)))
        print("    of which INCOMPLETE (crashed after printing) %d   %s"
              % (len(incomplete), ", ".join(incomplete)))
        print("  exited nonzero %2d   %s" % (len(red), ", ".join(red)))
        print("    (a red guard is not a stale transcript -- reported, and")
        print("     checked for CONSISTENCY with the transcript below)")

        check("no suite fails to complete in a pristine clone", not errored,
              "\n".join("%s: %s" % (m, results[m][2][0]) for m in errored))

        # THE OVERLAPPING LIMIT, REPORTED AND NOT ASSERTED.  The superseded
        # rule called a run completed iff its stdout was non-empty.  Where
        # nothing is INCOMPLETE the two rules agree by construction, so the
        # number below is the measurement of how much today's change moved --
        # and on a day with every dependency present the honest answer is
        # nothing.  A day on which it moves is a day a suite crashed.
        _old_rule_completed = [m for m, (v, _, _) in results.items()
                               if v not in ("ERROR", "NO_FILE", "MISSING")]
        _new_rule_completed = [m for m, (v, _, _) in results.items()
                               if v not in COMPLETION_FAULTS]
        print()
        print("  the 2026-11-06 rule against the one it replaces:")
        print("    suites the EMPTINESS rule called completed   : %d"
              % len(_old_rule_completed))
        print("    suites the SUMMARY-LINE rule calls completed : %d"
              % len(_new_rule_completed))
        print("    verdicts moved                               : %d   %s"
              % (len(incomplete), ", ".join(incomplete) or "none"))

        # A committed transcript and its suite's exit code are two claims
        # about the same run, and they can contradict each other.  This is
        # the check that would have caught 2026-10-04 -> 2026-10-15 directly:
        # a transcript reading "0 failed" beside code that exits 1.
        contradictions = []
        for module, transcript in _table:
            rc = exits.get(module)
            if rc is None or transcript is None:
                continue
            tpath = os.path.join(clone, transcript)
            if not os.path.exists(tpath):
                continue
            with open(tpath) as fh:
                ttext = fh.read()
            m = SUMMARY_RE.search(ttext)   # the constant, not a second spelling
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

        # ------------------------------------------------------------------
        # THE BUDGET, CHECKED (2026-11-04).  Two sentences that have been
        # carried in prose since 10-27 become numbers here: the registered
        # table still fits in one shell call, and every module excluded FOR
        # its cost really does not.  An exclusion whose reason is a number
        # nobody recomputes is a declaration, and the 10-28 table already
        # says declarations are the weak class.
        spent = sum(dt for _, dt, _ in results.values())
        elapsed = time.monotonic() - t_start
        print()
        if MODE == "assemble":
            # THE BUDGET MEANS SOMETHING ELSE ONCE THE CENSUS IS SPLIT
            # (2026-11-08).  `spent` here is the sum over every part, and
            # `elapsed` is this call only, so the single-call containment
            # `elapsed >= spent` is NOT a statement about anything -- it is
            # false by construction on any real split, and asserting it would
            # be the 11-07 fault inverted: comparing a whole against a part.
            #
            # What the split has to deliver is that EVERY CALL FITS.  That is
            # three facts, each checked where it is true: each part's own
            # elapsed clock contains its own suite time (asserted inside the
            # part, and re-read here from the record so an assemble on stale
            # parts cannot hide it); each part fits the cap; and this
            # assemble call fits the cap.  The total is printed as the cost of
            # the census and is deliberately NOT compared against one cap.
            _p_elapsed = [(d["index"], d["elapsed_s"], d["spent_s"])
                          for d in sorted(_assembled_parts,
                                          key=lambda x: x["index"])]
            _worst = max((e for _, e, _ in _p_elapsed), default=0.0)
            print("  CENSUS BUDGET, SPLIT ACROSS %d PARTS + THIS ASSEMBLE CALL"
                  % N_PARTS)
            print("    suite time, all parts    %6.1f s  over %d suites"
                  % (spent, len([1 for _, t in _table if t])))
            print("    largest part, elapsed    %6.1f s  <- the quantity the "
                  "cap is about" % _worst)
            print("    this assemble call       %6.1f s" % elapsed)
            print("    one-shell-call cap       %6.1f s" % CALL_CAP_S)
            print("    whole census, elapsed    %6.1f s  (REPORT -- spans %d "
                  "calls and is not bounded by one cap)"
                  % (sum(e for _, e, _ in _p_elapsed) + elapsed, N_PARTS + 1))
            _bad_cont = ["part %d: elapsed %.1f s < suite time %.1f s"
                         % (i + 1, e, s) for i, e, s in _p_elapsed if e < s]
            check("every part's recorded elapsed clock contains its own "
                  "suite time",
                  not _bad_cont,
                  "\n".join(_bad_cont) if _bad_cont else
                  "%d parts, each with elapsed >= spent in its own call.  "
                  "Re-read from the part files rather than recomputed, so a "
                  "part written by a call that mis-timed itself is caught "
                  "here and not only there." % len(_p_elapsed))
            _over = ["part %d: %.1f s" % (i + 1, e)
                     for i, e, _ in _p_elapsed if e >= CALL_CAP_S]
            check("every part fitted inside one shell call",
                  not _over,
                  ", ".join(_over) if _over else
                  "largest part %.1f s of a %.1f s cap" % (_worst, CALL_CAP_S))
            check("the assemble call fits inside one shell call",
                  elapsed < CALL_CAP_S,
                  "%.1f s elapsed of a %.1f s cap.  Sections 2b-5 run here "
                  "with no suite execution at all, so this is the tail's cost "
                  "on a %d-entry result set -- the quantity 11-07 recorded as "
                  "unmeasured." % (elapsed, CALL_CAP_S, len(results)))
            print()
            print("  MEASURED PER-SUITE COST, from the parts (2026-11-08).")
            print("  The table a cost-packed partition would need.  Emitted")
            print("  rather than used: see block_slice's docstring.")
            for _m, _dt, _ in sorted(
                    ((m, r[1], r[0]) for m, r in results.items()),
                    key=lambda x: -x[1]):
                print("    %-48s %6.1f s" % (_m, _dt))
        else:
            print("  SECTION 2 BUDGET")
            print("    suite time this run      %6.1f s  over %d suites"
                  % (spent, len([1 for _, t in _table if t])))
            print("    ELAPSED wall clock       %6.1f s  <- the quantity the cap "
                  "is about" % elapsed)
            print("    one-shell-call cap       %6.1f s" % CALL_CAP_S)
            print("    headroom, by elapsed     %6.1f s" % (CALL_CAP_S - elapsed))
            print("    headroom, by suite time  %6.1f s  (REPORT -- the number "
                  "the 11-04 check used, and it overstates)" % (CALL_CAP_S - spent))
            print("    unaccounted by suites    %6.1f s  (clone build, the loop's "
                  "own overhead, Sections 1/2b/3/4/5)" % (elapsed - spent))
            check("elapsed wall clock contains Section 2's suite time",
                  elapsed >= spent,
                  "elapsed %.1f s, suite time %.1f s, difference %.1f s.  A "
                  "negative difference would mean the suite timers charge time "
                  "this call did not spend -- the instrument broken, rather than "
                  "the budget exceeded." % (elapsed, spent, elapsed - spent))
            check("the whole module still fits inside one shell call",
                  elapsed < CALL_CAP_S,
                  "%.1f s elapsed of a %.1f s cap.  Measured 2026-11-07: Section "
                  "2 alone spends 176.4 s over 43 suites, and the module did not "
                  "complete inside the 179 s maximum this harness allows -- so a "
                  "run that prints this line at all is already the better case.  "
                  "2026-11-08: --part/--assemble is the path that does not need "
                  "this to be true." % (elapsed, CALL_CAP_S))
        oc_bad = []
        _excluded_modules = set(m for _c, m in TRANSCRIPT_EXCLUSIONS.values())
        for mod, secs in sorted(MEASURED_OVER_CAP_S.items()):
            if mod not in _excluded_modules:
                oc_bad.append("%s is timed as over-cap but is not excluded" % mod)
                continue
            # AMENDED 2026-11-07 with the budget verdict above: this
            # re-admission test used `spent`, so it measured a candidate
            # against a budget that had already overstated what was free.
            # AMENDED 2026-11-08.  Under --assemble the question is no
            # longer "does it fit beside everything else in one call" but
            # "does it fit in a part of its own", so the comparison is
            # against the per-part allowance and the verdict is reported
            # rather than scored: re-admitting a module is a decision about
            # the registry and this is not the session that makes it.
            _fits = (secs <= SECTION2_ALLOWANCE_S if MODE == "assemble"
                     else elapsed + secs <= CALL_CAP_S)
            if MODE == "assemble" and _fits:
                print("    NOTE %s: %.1f s fits inside the %.1f s per-part "
                      "allowance, so COST is no longer a reason to exclude "
                      "it once the census is split.  Not re-admitted here."
                      % (mod, secs, SECTION2_ALLOWANCE_S))
                continue
            if _fits:
                oc_bad.append("%s: %.1f s + %.1f s elapsed = %.1f s, which is "
                              "INSIDE the %.1f s cap, so cost is no longer a "
                              "reason to exclude it"
                              % (mod, secs, elapsed, elapsed + secs, CALL_CAP_S))
        _oc_names = ", ".join(
            "%s %.1fs" % (m.replace("graphene_", "").replace(".py", ""), s)
            for m, s in sorted(MEASURED_OVER_CAP_S.items()))
        # The detail string has to say which question was asked, because
        # 2026-11-08 the answer stopped being the same question.  Under
        # --assemble the NOTE lines above record that both excluded modules
        # WOULD fit in a part of their own, so a detail claiming they exceed
        # the cap would contradict the lines directly above it.
        if MODE == "assemble":
            _oc_detail = ("%d modules, each measured: %s.  Under the split "
                          "the cost exclusion is no longer sound for either "
                          "of them -- see the NOTE lines above.  This check "
                          "passes because nothing is MIS-registered, not "
                          "because the exclusions are still justified, and "
                          "re-admitting them is a registry decision left on "
                          "the standing list."
                          % (len(MEASURED_OVER_CAP_S), _oc_names))
        else:
            _oc_detail = ("%d modules, each measured today: %s.  Adding any "
                          "one to the %.1f s already elapsed exceeds the "
                          "%.1f s cap."
                          % (len(MEASURED_OVER_CAP_S), _oc_names,
                             elapsed, CALL_CAP_S))
        check("every module excluded for cost is measured over the remaining "
              "budget", not oc_bad,
              "\n".join("  " + b for b in oc_bad) if oc_bad else _oc_detail)

        # ------------------------------------------------------------------
        # C10 -- PLANTED CONTROLS ON THE BUDGET CROSSING, WITH EXACTLY KNOWN
        # ANSWERS (2026-11-07).  `section2_crossing` is the only new piece of
        # logic the budget repair adds, and the repair exists because the old
        # budget check had never been exercised in the one state it reports.
        # So the new one is planted on inputs whose answers are not a
        # plausible range but an exact integer -- including both sides of the
        # boundary, which is where a `>` / `>=` slip would live and which the
        # values below make exact in binary (1.5, 2.5, 4.0, 8.0 are all
        # representable, so the cumulative sums are not approximations).
        print()
        print("  C10 -- the budget crossing, on inputs with exact answers")
        check("  C10.1 a set that never reaches the allowance does not cross",
              section2_crossing([1.0, 2.0, 3.0], 10.0) is None,
              "cumulative 1, 3, 6 against 10 -- must be None, not index 2.  "
              "This is the must-give-nothing symmetry: a budget that is not "
              "exceeded has no crossing suite to name")
        check("  C10.2 the crossing index is the FIRST suite past the "
              "allowance, exactly",
              section2_crossing([4.0, 4.0, 4.0], 8.0) == 2,
              "cumulative 4, 8, 12 against 8 -- index 2, because 8 is not "
              "> 8.  An off-by-one or a >= would give 1 and this control "
              "would be red")
        check("  C10.3 landing EXACTLY on the allowance is not a crossing",
              section2_crossing([1.5, 2.5], 4.0) is None,
              "cumulative 1.5 then exactly 4.0 against 4.0.  All three "
              "values are exact in binary, so this is a statement about the "
              "comparison and not about float error")
        check("  C10.4 a single suite over the allowance crosses at index 0",
              section2_crossing([999.0], 150.0) == 0,
              "the degenerate case: the first suite alone exceeds, so the "
              "warning must name it rather than wait for a second")
        check("  C10.5 an empty suite table cannot cross",
              section2_crossing([], 0.0) is None,
              "no durations, zero allowance -- still None, because 0.0 is "
              "not > 0.0 and there is nothing to iterate")
        # And the relation the live check is gated on, exercised on a plant
        # rather than only on this run's own numbers: containment is what
        # makes `elapsed >= spent` safe to assert unconditionally.
        check("  C10.6 containment is what the live budget check asserts, "
              "and it holds on a plant where the tail is all the cost",
              (0.5 + 176.4) >= 176.4 and (0.0 + 176.4) >= 176.4,
              "elapsed = tail + spent with tail >= 0, so elapsed >= spent "
              "for every tail including zero.  The contingent quantity is "
              "the tail; the invariant is the inequality")
        check("  C10.7 the allowance leaves room below the cap, or the "
              "warning could not be printed before the kill",
              SECTION2_ALLOWANCE_S < CALL_CAP_S,
              "allowance %.1f s, cap %.1f s.  An allowance equal to the cap "
              "reproduces the 11-04 defect exactly: the crossing and the "
              "kill coincide and nothing is printed."
              % (SECTION2_ALLOWANCE_S, CALL_CAP_S))

        # ------------------------------------------------------------------
        # C11-C13 -- THE SPLIT, ON EXACTLY KNOWN ANSWERS (2026-11-08)
        #
        # The split's one job is that the suites a split run covers are the
        # suites a single-call run would have covered.  The partition is
        # therefore planted on exact integer answers, both at the boundary
        # cases (fewer suites than parts, zero suites) and as a set identity
        # over a grid, because a partition that dropped a suite would make
        # the census report a SMALLER stale count and read as an improvement.
        # ------------------------------------------------------------------
        print()
        print("  C11 -- the partition, on inputs with exact answers")
        check("  C11.1 the blocks of a 43-into-4 split have exactly the "
              "sizes 11, 11, 11, 10",
              [len(b) for b in partition_of(43, 4)] == [11, 11, 11, 10],
              "43 = 4*10 + 3, so the first three blocks take the remainder "
              "one each.  An exact integer answer, not a balance criterion")
        check("  C11.2 and they partition range(43) EXACTLY -- nothing "
              "dropped, nothing twice",
              flatten(partition_of(43, 4)) == list(range(43)),
              "the concatenation in order is [0..42].  This is the invariant "
              "the census depends on; the sizes above are a consequence of it "
              "plus contiguity")
        check("  C11.3 more parts than suites gives empty blocks rather than "
              "an error, and still partitions",
              ([len(b) for b in partition_of(3, 5)] == [1, 1, 1, 0, 0]
               and flatten(partition_of(3, 5)) == [0, 1, 2]),
              "3 suites into 5 parts: three blocks of one and two empty.  An "
              "empty part is a call that runs nothing and writes an empty "
              "part file, which --assemble must accept")
        check("  C11.4 zero suites gives every block empty",
              flatten(partition_of(0, 3)) == [] and
              [len(b) for b in partition_of(0, 3)] == [0, 0, 0],
              "the degenerate case, and the one a filter can produce")
        _grid_bad = []
        for _tot in range(0, 30):
            for _np in range(1, 8):
                if flatten(partition_of(_tot, _np)) != list(range(_tot)):
                    _grid_bad.append("(%d, %d)" % (_tot, _np))
        check("  C11.5 the partition identity holds over a 30x7 grid of "
              "(suites, parts)",
              not _grid_bad,
              "210 cases, 0 failures" if not _grid_bad else
              "failed at " + ", ".join(_grid_bad))
        _err = []
        for _bad in ((5, 0, 0), (5, 2, 2), (5, 2, -1)):
            try:
                block_slice(*_bad)
            except ValueError:
                continue
            _err.append(repr(_bad))
        check("  C11.6 an out-of-range part index or a zero part count "
              "RAISES rather than returning a silently empty block",
              not _err,
              "block_slice(5,0,0), (5,2,2) and (5,2,-1) all raise ValueError. "
              "A silently empty block is the failure mode that would make a "
              "census miss a third of its suites and still print a summary"
              if not _err else "did not raise for " + ", ".join(_err))

        print()
        print("  C12 -- the part fingerprint, which is what refuses a "
              "mismatched assembly")
        _t1 = [("a.py", "a_output.txt"), ("b.py", "b_output.txt")]
        _t2 = [("a.py", "a_output.txt"), ("c.py", "c_output.txt")]
        _t3 = [("b.py", "b_output.txt"), ("a.py", "a_output.txt")]
        check("  C12.1 the same table gives the same fingerprint",
              table_fingerprint(_t1) == table_fingerprint(list(_t1)),
              "a content-addressed id, so equality is the requirement")
        check("  C12.2 a different table gives a different one",
              table_fingerprint(_t1) != table_fingerprint(_t2),
              "one suite swapped out of two")
        check("  C12.3 and so does the SAME suites in a different ORDER",
              table_fingerprint(_t1) != table_fingerprint(_t3),
              "order matters because block_slice is positional: two parts "
              "run against differently-ordered tables would cover the same "
              "suites under the single-call reading and different ones under "
              "the split.  This is the case a set-based fingerprint would "
              "have missed")
        check("  C12.4 a transcript-less suite is distinguished from one "
              "whose transcript is the empty string",
              table_fingerprint([("a.py", None)])
              != table_fingerprint([("a.py", "")]),
              "None is encoded as a sentinel rather than as '', so a "
              "registry change from unregistered to registered-as-empty "
              "moves the fingerprint.  THIS CONTROL FAILED ON ITS FIRST RUN: "
              "the fingerprint's first form used `t or \"\"` and mapped both "
              "to the same bytes")

        print()
        print("  C13 -- --part argument parsing, on answers known by "
              "inspection")
        check("  C13.1 --part 2/4 is part index 1 of 4",
              parse_mode(["x", "--part", "2/4"]) == ("part", 4, 1),
              "1-based for the operator, 0-based inside.  C13.2 is why the "
              "translation is checked rather than assumed")
        check("  C13.2 --part 0/4 and --part 5/4 are REFUSED",
              all(_raises_systemexit(["x", "--part", s])
                  for s in ("0/4", "5/4", "2", "2/0")),
              "an operator typing 0/4 has made an error, and running part 1 "
              "instead would produce a part file that covers the wrong block "
              "under the right name")
        check("  C13.3 no arguments is the single-call mode, unchanged",
              parse_mode(["x"]) == ("single", None, None),
              "the committed transcript is a claim about this path, so the "
              "default must not move")
        check("  C13.4 --assemble --parts=4 parses to an assemble of 4",
              parse_mode(["x", "--assemble", "--parts=4"]) == ("assemble", 4, None),
              "and the order of the two flags does not matter: %r"
              % (parse_mode(["x", "--parts=4", "--assemble"]),))

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
        # NOTE 2026-11-04: this control's donor is now REGISTERED in SUITES
        # rather than excluded.  The control is unaffected, and deliberately
        # left pointing at it: writes_own_transcript is a predicate on the
        # SOURCE, not on the tables, and a module being registered no longer
        # has any bearing on whether it writes its own transcript.  That
        # independence is exactly what today's change bought.
        check("  C5.3 and CONFIRMED for one that does",
              writes_own_transcript(repo, "graphene_default_flip.py",
                                    "default_flip_output.txt"),
              "graphene_default_flip.py line ~745 opens it 'w', and is now a "
              "registered suite as well -- the predicate does not care")

        # C9 -- THIS MODULE'S OWN EXIT CODE, ASKED OF THE SOURCE BY AST.
        # Added 2026-11-06.  Section 5 checks this property for every suite
        # but this one, and the SELF exclusion is why.  The question needs no
        # clone and no second run: it is the shape of main()'s returns.  AST
        # and not a regex -- see the amendment; the regex attempt was wrong.
        import ast as _ast
        _self_src = open(os.path.join(repo, SELF)).read()
        _main = next(n for n in _ast.walk(_ast.parse(_self_src))
                     if isinstance(n, _ast.FunctionDef) and n.name == "main")
        _rets = [n for n in _ast.walk(_main) if isinstance(n, _ast.Return)
                 and n.value is not None]
        _consts = [r for r in _rets if isinstance(r.value, _ast.Constant)]
        _ifexps = [r for r in _rets if isinstance(r.value, _ast.IfExp)]
        # AMENDED 2026-11-08, AND THE SPLIT IS WHAT FOUND IT.
        #
        # SUPERSEDED FORM, kept so the change is legible:
        #
        #     not _consts and len(_ifexps) == 1
        #
        # It failed the moment --part/--assemble were added, because those
        # give main() three more exits: `return 2` twice on the
        # cannot-assemble path and a second failure-conditional return at the
        # end of a part.  Read literally the old form says "main() has
        # exactly one exit and it is conditional", which is a statement about
        # main()'s SHAPE.  The property it was protecting is narrower and
        # does not mention the number of exits at all:
        #
        #     NO PATH OUT OF main() REPORTS SUCCESS UNCONDITIONALLY.
        #
        # So the three requirements below replace the two above, and they are
        # strictly stronger in the direction that matters: the old form
        # permitted a second conditional gated on the wrong name (it only
        # ever inspected `_ifexps[0]`), and permitted no constant exit at all
        # even for a refusal to run, which is a different thing from a clean
        # run and must be distinguishable from one.
        _zero_consts = [r for r in _consts if r.value.value in (0, None)]
        _bad_ifexp = [r for r in _ifexps
                      if not any(isinstance(nn, _ast.Name) and nn.id == "_FAIL"
                                 for nn in _ast.walk(r))]
        check("C9 no path out of main() reports success unconditionally",
              not _zero_consts,
              "%d conditional exit(s) and %d constant exit(s), of which %d "
              "return 0 or None.  A bare `return 0` here is a transcript "
              "reading 'N failed' beside an exit code of 0, which is the "
              "2026-10-04..10-15 state this audit exists to catch"
              % (len(_ifexps), len(_consts), len(_zero_consts)))
        check("C9b every failure-conditional exit is gated on the FAIL list "
              "this run actually appends to, not on some other name",
              bool(_ifexps) and not _bad_ifexp,
              "%d conditional exit(s), all of them naming _FAIL.  The old "
              "form of this check inspected only the FIRST conditional "
              "return, so the part-mode exit added today would have been "
              "unchecked by it -- a conditional on the wrong name exits 0 "
              "while red and passes a control that only reads its sibling "
              "(2026-10-16)" % len(_ifexps)
              if not _bad_ifexp else
              "%d conditional exit(s) do not mention _FAIL" % len(_bad_ifexp))
        check("C9c and every constant exit is NONZERO -- a refusal to run is "
              "not a clean run",
              all(isinstance(r.value.value, int) and r.value.value != 0
                  for r in _consts),
              "%d constant exit(s): %s.  These are the cannot-assemble paths, "
              "which report that the census did not happen.  Returning 0 "
              "there would make a missing part file indistinguishable from a "
              "green census, which is exactly the conflation C9 exists for"
              % (len(_consts), [r.value.value for r in _consts]))

        print("\n" + "=" * 74)
        print("SUMMARY")
        print("=" * 74)
        print("  commit audited : %s" % commit)
        for lab in _PASS:
            print("  PASS  %s" % lab)
        for lab in _FAIL:
            print("  FAIL  %s" % lab)
        print("\nTOTAL: %d passed, %d failed" % (len(_PASS), len(_FAIL)))
        # FIXED 2026-11-06.  This was a bare `return 0`: see the amendment in
        # the module docstring.  The convention this module's own Section 5
        # checks for in every other suite now holds for this one too.
        return 1 if _FAIL else 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    MODE, _n, _i = parse_mode(sys.argv)
    if MODE == "part":
        N_PARTS, PART_INDEX = _n, _i
    elif MODE == "assemble":
        N_PARTS = _n or int(os.environ.get("GRAPHENE_CENSUS_PARTS", "4"))
        PART_INDEX = 0
    sys.exit(main())
