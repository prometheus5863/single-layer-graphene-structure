#!/usr/bin/env python3
"""
graphene_census_split_equivalence.py -- 2026-11-08

THE SPLIT CENSUS MUST SEE WHAT THE SINGLE-CALL CENSUS SAW.

graphene_pristine_transcript_audit.py gained a --part/--assemble mode today,
because 2026-11-07 measured its Section 2 at 176.4 s over 43 suites against a
180 s shell-call cap and recorded that the module can no longer complete in one
call.  The split is a change to the SHAPE of an instrument, and the hazard in
that is specific and worth naming before any code:

    A SPLIT THAT DROPPED A SUITE WOULD REPORT A SMALLER STALE COUNT,
    AND A SMALLER STALE COUNT READS AS AN IMPROVEMENT.

So does a split that silently emptied a check.  Sections 2b and 3 of the census
compare, for every suite that writes its own transcript, the bytes that suite's
own run left in the pristine clone against the committed blob.  Those bytes
exist only inside the call that ran the suite.  An --assemble phase that
rebuilt a fresh clone would find the COMMITTED bytes sitting there instead, and
every one of those comparisons would come back "rewrote it byte-identically" --
not because the tree is healthy but because nothing had run.  The checks would
be green, numerous, and about nothing.  That is the same fault class as the
11-07 budget check, which reported 3.6 s of headroom on a run that could not
finish.

This module is the test that the two paths agree, and it is a COMPARISON OF
OUTPUT rather than an argument about code.  It runs the census both ways on a
deliberately tiny suite table and requires the regions that are claims about
the repository to be byte-identical, after masking the only thing that cannot
be: elapsed times.

WHAT IS *REQUIRED* IDENTICAL, AND WHAT IS *DECLARED* DIFFERENT
---------------------------------------------------------------
An equivalence test that demanded the whole transcript match would be false,
and the honest form is to say exactly where the two paths are allowed to
diverge and to require identity everywhere else:

  REGION A  Sections 1-2          DIFFERENT by declaration.  --assemble prints
                                  its provenance (which parts, at which
                                  commit, each part's own clock) where a
                                  single call prints the suites as they run.
  REGION B  Sections 2b, 3, 4     REQUIRED IDENTICAL.  Every claim the census
            and Section 5's       makes about the repository lives here.  This
            coverage checks       is the region the post-run-bytes hazard
                                  above would have hollowed out.
  REGION C  the budget block      DIFFERENT by declaration, and NECESSARILY:
                                  `elapsed >= spent` is a statement about one
                                  call, and the assembled totals span several.
                                  See the amendment in the census module.
  REGION D  C1-C13, the controls  REQUIRED IDENTICAL.
  REGION E  the SUMMARY list      DIFFERENT by declaration -- it lists the
                                  budget checks by name, and those differ.

MUTATION TESTED.  A comparison that cannot fail is not evidence, and 2026-09-17
and 2026-11-07 both recorded suites that passed against deliberately broken
code.  Section 3 injects two defects into a COPY of the census and requires
this comparison to catch both:

  M1  block_slice drops the last index of every block -- the partition
      silently stops being one, which is the under-reporting hazard.
  M2  --assemble does not restore the post-run transcript bytes -- the
      hollowing-out hazard.

Exit code is failure-conditional, per the convention Section 5 of the census
checks for in every suite.
"""
import base64
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CENSUS = "graphene_pristine_transcript_audit.py"

# Three cheap suites, chosen for cost and for COVERAGE OF THE HAZARD: at least
# one of them must write its own transcript, or Region B's self-writing
# comparisons would be empty and M2 would be undetectable by construction.
# That requirement is checked, not assumed -- see S0 below.
# S0b FAILED ON THE FIRST RUN OF THIS MODULE, AND THAT IS WHY IT IS HERE.
# The first table was potential_declaration, capture_leverage_figure and
# write_read_guard -- three of the cheapest suites, chosen for cost alone.
# None of the four self-writing suites in the census's table is among them,
# so Region B's self-writing comparisons were empty, M2 changed nothing, and
# M2 came back UNCAUGHT.  The precondition is what distinguished "the split
# does not restore the bytes and nothing noticed" from "there were no bytes
# to restore"; without it, M2's result would have been read as a defect in
# the split.  default_flip is the cheapest of the four (~1.3 s).
ONLY = [
    "graphene_potential_declaration.py",
    "graphene_default_flip.py",            # SELF_WRITING -- see S0b
    "graphene_write_read_guard.py",
]
N_PARTS = 2

_PASS, _FAIL = [], []


def check(label, ok, detail=""):
    (_PASS if ok else _FAIL).append(label)
    print("  [%s] %s" % ("PASS" if ok else "FAIL", label))
    if detail:
        for line in detail.rstrip("\n").split("\n"):
            print("        " + line)
    return ok


# ---------------------------------------------------------------------------
# Masking: the one thing the two paths cannot agree on
# ---------------------------------------------------------------------------
# Deliberately narrow.  A mask that swallowed all numbers would also swallow
# the stale counts, the transcript byte counts and the drift magnitudes -- the
# quantities this test exists to compare.  Only a number immediately followed
# by a seconds unit is masked.
_T = [
    (re.compile(r"\d+\.\d+ ?s\b"), "<T>s"),
    (re.compile(r"\d+\.\d+e[-+]\d+"), "<E>"),
]


def mask(text):
    for rx, sub in _T:
        text = rx.sub(sub, text)
    return text


REGION_MARKS = [
    ("A", re.compile(r"^SECTION 1 ", re.M)),
    ("B", re.compile(r"^SECTION 2b ", re.M)),
    ("C", re.compile(r"^  (?:SECTION 2 BUDGET|CENSUS BUDGET)", re.M)),
    ("D", re.compile(r"^  C10 -- ", re.M)),
    ("E", re.compile(r"^SUMMARY$", re.M)),
]
REQUIRED_IDENTICAL = ("B", "D")


def regions(text):
    """{name: text} for the five regions, split on the markers above.

    A missing marker is returned as an empty region rather than raising: a
    run that died early must produce a comparison that FAILS, not a crash
    that looks like an error in this module.
    """
    pos = []
    for name, rx in REGION_MARKS:
        m = rx.search(text)
        pos.append((name, m.start() if m else None))
    out = {}
    for i, (name, start) in enumerate(pos):
        if start is None:
            out[name] = ""
            continue
        nxt = None
        for _, s2 in pos[i + 1:]:
            if s2 is not None:
                nxt = s2
                break
        out[name] = text[start:nxt if nxt is not None else len(text)]
    return out


def run(args, cache, census=None, extra_env=None):
    """One census invocation.

    `census` defaults at CALL time and not at def time: a default ARGUMENT
    would capture the module-level constant when this def executes, which is
    the pattern graphene_mixed_capture_repair_audit.py censuses repo-wide and
    which this module's first draft added two new instances of -- caught by
    that census on the day it was written.
    """
    census = CENSUS if census is None else census
    env = dict(os.environ)
    env["GRAPHENE_CENSUS_ONLY"] = ",".join(ONLY)
    env["GRAPHENE_CENSUS_CACHE"] = cache
    env["MPLBACKEND"] = "Agg"
    env.update(extra_env or {})
    t0 = time.time()
    p = subprocess.run([sys.executable, census] + args, cwd=HERE, env=env,
                       capture_output=True, text=True, timeout=600)
    return p.stdout, p.returncode, time.time() - t0


def both_paths(cache_root, census=None):
    """(single_output, assembled_output, per-call seconds).

    census defaults at call time -- see run().
    """
    census = CENSUS if census is None else census
    shutil.rmtree(cache_root, ignore_errors=True)
    os.makedirs(cache_root, exist_ok=True)
    single, _rc, t_single = run([], cache_root, census)
    times = [("single", t_single)]
    for i in range(N_PARTS):
        _o, _rc, _t = run(["--part", "%d/%d" % (i + 1, N_PARTS)], cache_root, census)
        times.append(("part %d" % (i + 1), _t))
    asm, _rc, t_asm = run(["--assemble", "--parts=%d" % N_PARTS], cache_root, census)
    times.append(("assemble", t_asm))
    return single, asm, times


def compare(single, asm):
    """[(region, verdict, detail)] over all five regions."""
    rs, ra = regions(mask(single)), regions(mask(asm))
    out = []
    for name, _rx in REGION_MARKS:
        same = rs[name] == ra[name]
        out.append((name, same, len(rs[name]), len(ra[name])))
    return out


def main():
    print("=" * 74)
    print("CENSUS SPLIT EQUIVALENCE -- 2026-11-08")
    print("Does --part/--assemble see what the single call saw?")
    print("=" * 74)

    census_src = open(os.path.join(HERE, CENSUS), encoding="utf-8").read()

    # -- Section 0: the test's own preconditions ---------------------------
    print("\n" + "-" * 74)
    print("SECTION 0 -- this test's preconditions, checked rather than assumed")
    print("-" * 74)
    sys.path.insert(0, HERE)
    import graphene_pristine_transcript_audit as cen

    registered = dict(cen.SUITES)
    missing = [m for m in ONLY if m not in registered]
    check("S0a every suite this test names is registered in the census's own "
          "SUITES table",
          not missing,
          "%d suites: %s" % (len(ONLY), ", ".join(ONLY)) if not missing else
          "not registered: " + ", ".join(missing))
    self_writing = [m for m in ONLY
                    if m not in missing
                    and cen.writes_own_transcript(HERE, m, registered[m])]
    check("S0b at least one of them WRITES ITS OWN TRANSCRIPT, or Region B's "
          "self-writing comparisons are empty and M2 is undetectable",
          bool(self_writing),
          "self-writing: %s.  This is the precondition that makes the "
          "mutation in Section 3 a real test rather than a run of a vacuous "
          "one -- M2 breaks exactly the path these suites exercise."
          % (", ".join(self_writing) or "NONE"))
    check("S0c the census still exposes the three pieces this test drives",
          all(hasattr(cen, n) for n in
              ("block_slice", "suite_table", "parse_mode")),
          "block_slice, suite_table, parse_mode")

    # -- Section 1: the two paths, on the real census -----------------------
    print("\n" + "-" * 74)
    print("SECTION 1 -- the baseline comparison")
    print("-" * 74)
    cache = tempfile.mkdtemp(prefix="census_eq_")
    try:
        single, asm, times = both_paths(cache)
        for name, t in times:
            print("  %-10s %6.1f s" % (name, t))
        rows = compare(single, asm)
        print()
        print("  region  verdict      single   assembled")
        for name, same, n1, n2 in rows:
            print("  %-7s %-12s %7d %11d"
                  % (name, "IDENTICAL" if same else "DIFFERENT", n1, n2))
        bad = [n for n, same, _, _ in rows if n in REQUIRED_IDENTICAL and not same]
        check("E1 Regions B and D are BYTE-IDENTICAL between the two paths, "
              "after masking elapsed times",
              not bad,
              "Region B is Sections 2b, 3, 4 and Section 5's coverage checks "
              "-- every claim the census makes about the repository.  Region "
              "D is the planted controls.  Masking is narrow by design: only "
              "a number followed by a seconds unit is replaced, so stale "
              "counts, byte counts and drift magnitudes are compared as "
              "written." if not bad else
              "differing regions: " + ", ".join(bad))
        declared = [n for n, same, _, _ in rows
                    if n not in REQUIRED_IDENTICAL and not same]
        check("E2 and Regions A, C and E DO differ -- the declaration is "
              "exercised, not vacuous",
              sorted(declared) == ["A", "C", "E"],
              "A is the provenance header, C is the budget block (which MUST "
              "differ: `elapsed >= spent` is a one-call statement and the "
              "assembled totals span four calls), E is the summary's list of "
              "budget-check names.  If any of the three came back identical, "
              "either the split is not doing what it claims or this test is "
              "comparing the wrong thing.  Differing: %s" % declared)
        check("E3 every part, and the assemble call, fits the one-call cap",
              all(t < cen.CALL_CAP_S for _n, t in times[1:]),
              "worst %.1f s against a %.1f s cap, on a %d-suite table.  This "
              "is not the claim that the REAL table fits -- that is measured "
              "by the census itself -- only that the mechanism's own overhead "
              "does not eat the budget it exists to free."
              % (max(t for _n, t in times[1:]), cen.CALL_CAP_S, len(ONLY)))

        # -- Section 2: what Region B actually contains --------------------
        # A region can be identical because it is empty.  2026-11-05 retired
        # a check for being satisfiable only by a defect; the mirror of that
        # is a check satisfiable by a vacuum, and this is where one would
        # hide.
        print("\n" + "-" * 74)
        print("SECTION 2 -- Region B is not identical by being empty")
        print("-" * 74)
        rb = regions(mask(asm))["B"]
        n_checks = len(re.findall(r"^  \[(?:PASS|FAIL)\]", rb, re.M))
        n_e1 = len(re.findall(r"object-store read is byte-identical", rb))
        print("  Region B: %d bytes, %d check lines" % (len(rb), n_checks))
        check("E4 Region B carries a non-trivial number of checks",
              n_checks >= 8 and len(rb) > 1500,
              "%d check lines over %d bytes.  A split that hollowed out "
              "Sections 2b-5 would leave this region identical AND tiny, so "
              "E1 is gated on content as well as on equality." % (n_checks, len(rb)))
        check("E5 and it contains the two-reads comparison by name",
              n_e1 >= 1,
              "the E1 check of Section 2b, which is the one that would read "
              "as green on a clone nothing had run in")

        # -- Section 3: mutation ------------------------------------------
        print("\n" + "-" * 74)
        print("SECTION 3 -- MUTATION.  Two defects, both of which must be caught")
        print("-" * 74)
        mutants = [
            ("M1 block_slice drops the last index of every block",
             "    return list(range(start, start + size))",
             "    return list(range(start, max(start, start + size - 1)))",
             "the partition silently stops being one.  A census built from "
             "such parts covers fewer suites, reports a smaller STALE count, "
             "and reads as an improvement"),
        ]
        for label, old, new, why in mutants:
            if census_src.count(old) != 1:
                check(label + " -- PLANTED", False,
                      "the anchor %r occurs %d times in the census, not once: "
                      "the mutation could not be planted, so this row is a "
                      "FAULT IN THE TEST and not a passing control"
                      % (old, census_src.count(old)))
                continue
            mdir = tempfile.mkdtemp(prefix="census_mut_", dir=HERE)
            mpath = os.path.join(mdir, CENSUS)
            try:
                # The mutant must sit in a directory that is NOT the repo, or
                # it would resolve `repo` to the real tree and audit it; it
                # must still be able to find the repo, which it does because
                # make_pristine_clone clones from `repo`.  A subdirectory of
                # the repo is the one place both hold.
                os.symlink(os.path.join(HERE, ".git"), os.path.join(mdir, ".git"))
                with open(mpath, "w", encoding="utf-8") as fh:
                    fh.write(census_src.replace(old, new, 1))
                mcache = tempfile.mkdtemp(prefix="census_eq_mut_")
                try:
                    ms, ma, _mt = both_paths(
                        mcache, census=os.path.join(os.path.basename(mdir), CENSUS))
                    mrows = compare(ms, ma)
                    mbad = [n for n, same, _, _ in mrows
                            if n in REQUIRED_IDENTICAL and not same]
                    empty = not ma.strip()
                    caught = bool(mbad) or empty
                    if mbad:
                        how = "Regions %s now differ" % ", ".join(mbad)
                    elif empty:
                        how = "the mutated assemble produced no output at all"
                    else:
                        how = ("NOT CAUGHT -- Regions B and D came back "
                               "identical with the defect in place, so this "
                               "comparison is blind to it")
                    check(label, caught, "%s.\n        %s" % (why, how))
                finally:
                    shutil.rmtree(mcache, ignore_errors=True)
            finally:
                shutil.rmtree(mdir, ignore_errors=True)

        # -- Section 4: the restoration, with the divergence PLANTED -------
        #
        # M2 BELONGS HERE AND NOT IN SECTION 3, AND FINDING THAT OUT IS THE
        # SECOND THING THIS MODULE MEASURED.  Run end-to-end, the mutation
        # "--assemble does not restore the post-run bytes" came back NOT
        # CAUGHT even with a self-writing suite in the table, and the reason
        # is not that the comparison is weak: every self-writing suite in
        # this tree currently rewrites its transcript BYTE-IDENTICALLY (the
        # state 2026-11-05 recorded when it retired the old E2 for being
        # satisfiable only by a defect).  If post-run bytes equal the
        # committed blob, restoring them and not restoring them are the same
        # operation, and no end-to-end comparison on a healthy tree can tell
        # the two apart.
        #
        # So the divergence is PLANTED rather than waited for -- the same
        # discipline E2c1/E2c2 apply to the circularity predicate, one level
        # out.  The part file's stored post_bytes are edited to differ from
        # the committed blob by one line, and the two assemble paths must
        # then disagree:
        #
        #     restoration intact   -> Section 2b reports [CONFIRMED]   E2r
        #     restoration removed  -> Section 2b reports [NOT_EXERCISED] E2r
        #
        # A control that requires the tree to be unhealthy is not a control.
        print("\n" + "-" * 74)
        print("SECTION 4 -- M2, with the divergence planted into the part file")
        print("-" * 74)
        pcache = tempfile.mkdtemp(prefix="census_eq_plant_")
        mdir = tempfile.mkdtemp(prefix="census_mut_", dir=HERE)
        try:
            for i in range(N_PARTS):
                run(["--part", "%d/%d" % (i + 1, N_PARTS)], pcache)
            target, patched = "default_flip_output.txt", None
            for i in range(N_PARTS):
                fp = cen.part_path(N_PARTS, i)
                fp = os.path.join(pcache, os.path.basename(fp))
                with open(fp, encoding="utf-8") as fh:
                    d = json.load(fh)
                if target in d.get("post_bytes", {}):
                    raw = base64.b64decode(d["post_bytes"][target])
                    d["post_bytes"][target] = base64.b64encode(
                        raw + b"\nPLANTED DIVERGENCE -- one line, not a number\n"
                    ).decode("ascii")
                    with open(fp, "w", encoding="utf-8") as fh:
                        json.dump(d, fh)
                    patched = i
                    break
            check("E6 the planted divergence could be written into a part "
                  "file, so Section 4 is running on what it says it is",
                  patched is not None,
                  "%s found in part %s of %d and extended by one line"
                  % (target, (patched + 1) if patched is not None else "?",
                     N_PARTS)
                  if patched is not None else
                  "%s appears in no part file's post_bytes: the suite is in "
                  "the table but its bytes were not carried, which would be "
                  "a defect in the part writer" % target)

            os.symlink(os.path.join(HERE, ".git"), os.path.join(mdir, ".git"))
            m2_old = '            for _t, _b in _post_bytes.items():'
            m2_new = '            for _t, _b in list(_post_bytes.items())[:0]:'
            assert census_src.count(m2_old) == 1
            with open(os.path.join(mdir, CENSUS), "w", encoding="utf-8") as fh:
                fh.write(census_src.replace(m2_old, m2_new, 1))
            rel_mut = os.path.join(os.path.basename(mdir), CENSUS)

            good, _rc, _t = run(["--assemble", "--parts=%d" % N_PARTS], pcache)
            mut, _rc, _t = run(["--assemble", "--parts=%d" % N_PARTS], pcache,
                               census=rel_mut)
            g_conf = "[CONFIRMED] E2r" in good
            m_conf = "[CONFIRMED] E2r" in mut
            check("E7 with the restoration intact, the planted divergence "
                  "REACHES Section 2b and is reported",
                  g_conf and target in good,
                  "Section 2b prints [CONFIRMED] E2r and names %s.  This is "
                  "the path the split had to preserve: the bytes a part's own "
                  "run left behind, carried across a call boundary and "
                  "compared against the committed blob." % target)
            check("M2 --assemble that does NOT restore the post-run bytes is "
                  "CAUGHT by the same plant",
                  g_conf and not m_conf,
                  "intact: [CONFIRMED] E2r.  mutated: %s.  The mutant "
                  "compares the committed blob against itself, so every "
                  "self-writing suite reads 'rewrote it byte-identically' -- "
                  "green, numerous, and about nothing.  Run end-to-end on "
                  "this tree the same mutation is INVISIBLE, because the "
                  "tree's self-writing suites already rewrite identically; "
                  "planting the divergence is what makes the control a "
                  "control rather than a wait."
                  % ("[CONFIRMED] E2r" if m_conf else
                     "[NOT_EXERCISED] E2r"
                     if "[NOT_EXERCISED] E2r" in mut else "no E2r line"))
        finally:
            shutil.rmtree(pcache, ignore_errors=True)
            shutil.rmtree(mdir, ignore_errors=True)

        print("\n" + "=" * 74)
        print("SUMMARY")
        print("=" * 74)
        for lab in _PASS:
            print("  PASS  %s" % lab)
        for lab in _FAIL:
            print("  FAIL  %s" % lab)
        print("\nTOTAL: %d passed, %d failed" % (len(_PASS), len(_FAIL)))
        return 1 if _FAIL else 0
    finally:
        shutil.rmtree(cache, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
