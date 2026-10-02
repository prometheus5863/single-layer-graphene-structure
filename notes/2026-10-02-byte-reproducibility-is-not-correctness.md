# A transcript is a claim about the code beside it — and agreeing is not being right

**2026-10-02.** Session note for `single-layer-graphene-structure`.
Instrument: `graphene_pristine_transcript_audit.py`.
Transcript: `pristine_transcript_audit_output.txt`.

---

## What this session was asked to do

The 2026-10-01 entry closed with a new top item, created by its own worst
finding. That finding was that the 09-30 before/after measurement had been
failing since the commit that introduced it: it passed in the working tree of
the run that wrote it and in no clone afterwards. The entry named the cheapest
test for that class of fault — *run the suite from a pristine clone of the
commit* — and recorded that this repository had never run it.

This session ran it, for every suite, and asked a question one step stronger
than "does it still pass":

> **Does the committed code still produce the committed transcript?**

Those are different questions, and the gap between them is where this
repository had been losing track of itself. The transcripts here
(`*_output.txt`) are captured stdout, written by hand with a shell redirect. A
session that changes a module is free to forget the redirect. Nothing checked
it, and — this is the part that matters — **`git status` cannot see the
result**, because the transcript is unmodified. It is the code that moved out
from under it.

Opening state at `2d51dd8`: **4 STALE, 1 NUMERIC_DRIFT**, and one of the new
instrument's own controls failing on a premise of its own. Closing state:
**0 STALE**.

---

## 1. The sharpest finding: the strongest verdict was wrong

`graphene_band_structure_audit.py` classified **IDENTICAL** — the audit's
highest verdict. Code and transcript agreed byte for byte. And the transcript
contained this, under a heading reading `RESULT 1 -- the shipped band structure
has NO DIRAC POINT`:

```
  minimum gap over the whole shipped path = 0.0000 eV
  gap at the index labelled K             = 0.0000 eV
  ... The shipped figure shows a 0.0 eV gap at its own K label.
```

**A 0.0 eV gap at K is the Dirac point.** The sentence offers the absence of
the bug as the evidence for the bug. `RESULT 1` reports the *superseded*
spectrum — that is its entire subject — and it was reading the *corrected*
path. The repository's own bug report has been quoting the fixed value as the
fault since the fix landed five weeks ago.

The lesson is not about band structure:

> **Byte-reproducibility is a claim that two artefacts agree. It is not a claim
> that either is right.** No amount of reproducibility checking reaches an
> error that has been faithfully reproduced into both halves. What reached this
> one was reading the numbers against the sentence printed beside them.

This is the first finding in this series that the day's own new instrument
could not have found, and it was found while looking at that instrument's
output.

## 2. V5 could only fail, and the oracle it needed was built the same day

`V5 bug reproduced exactly` is the validation whose job is to reproduce the
superseded spectrum. Transcript history, which is unambiguous:

| commit | date | V5 | tally |
|---|---|---|---|
| `64c6347` | 2026-09-28 | **PASS** | validations: 5/5 |
| `fef532d9` | 2026-09-28 | **FAIL** | validations: 4/5 |

…and FAIL in every clone since. Same mechanism as §1: V5 called
`calculate_band_structure()` with no argument — the corrected path — so from
the moment the fix landed it was asking the fixed code to reproduce the bug,
which it cannot do by construction.

The twist is where the fix was. `legacy=True`, which preserves the superseded
path exactly, **was added in `64c6347`, the same commit as the fix**, explicitly
under the 09-24 rule that *a deleted implementation destroys the only oracle
available for judging its replacement*. The oracle was preserved, and the one
validation whose purpose is to exercise it was never pointed at it.

**Why it survived five weeks:** the suite prints `FAIL` on a summary line and
exits 0. A standing FAIL in a tally is worse than a missing check, because it
teaches a reader that 4/5 is this suite's normal state — which is exactly the
condition under which a *real* regression here would be invisible. A check that
cannot fail is the 09-28 hazard; **a check that cannot pass is its mirror, and
it is louder and therefore easier to stop hearing.**

V5 now reads the legacy path; new **V5b** asserts the default path *has* a
Dirac point. Before: V5 could only fail. Now: V5 fails if the oracle is lost,
V5b fails if the fix is reverted. Validations 4/5 → **8/8**.

## 3. The third site: the 09-28 missing π is still in the repository

`plot_fermi_surface()` carried its own inline copy of the wrong zone corner:

```python
k_K = 4.0/(3*np.sqrt(3))      # = 0.76980 — π times too small
```

So *"Fermi Surface of Single-Layer Graphene (Around K Point)"* was a contour
plot centred where the gap is **14.40 eV**: a figure of the Dirac cone with no
Dirac cone in it.

| centre | value | gap there |
|---|---|---|
| superseded (inline literal) | 0.7698003589195 | **14.401937 eV** |
| corrected (`HIGH_SYMMETRY`) | 2.4183991523123 | 2.487e-15 eV |

**Why nothing caught it.** The function is called by no suite, and
`fermi_surface.png` is not committed — so the only artefact that could have
disagreed with the code does not exist. The 09-28 fix corrected the *constant*
and verified it against the k-path and the DOS, the two sites it knew about.

> **A literal is not fixed by fixing the constant it duplicates.**
> `graphene_dead_name_sweep.py` exists because an unread *name* is a hazard.
> This is the mirror case: a read *value* that is not a name at all, and
> therefore invisible to every name-based instrument here.

Now reads `fermi_surface_grid_centre()`, which exists so the choice is
assertable from outside a plotting function. **V6** requires the gap at that
centre to sit on the floor; **V6b** forbids an inline copy of the superseded
corner outside `HIGH_SYMMETRY_LEGACY`. V6b was mutation-tested: reinjecting the
literal takes the suite 8/8 → 7/8 and names the line.

## 4. A ranked table of 21 equal keys

The pristine audit's first run reported `graphene_differential_crossover_model.py`
STALE on a *structural* difference: the committed transcript named `Ti/Pd` in a
table and a pristine run of the *same commit* named `Cr/Cu`, **with every
printed number identical**.

The cause was printed four lines below the table by the function's own prose:
`asym` is zero for all 21 pairs to machine precision, because N is exactly odd
in s (Validation 6). So `sort(key=-asym)` sorted 21 equal keys, `rows[:6]` took
an arbitrary six, and *which* six was decided by the last bits of a quantity the
module had already proved is zero. Measured spread: **5.574e-15**.

No number and no verdict moves — P2 and P4 are still falsified by the identity.
What was wrong is that **a reader counts a ranked table as a claim that the top
row differs from the bottom one.** Now tie-broken lexicographically (so the
transcript is reproducible at all) and the degeneracy is printed, so the
ordering cannot be read as information. This answers the 09-23 "is there a
RANKING that is a step artefact" item in its *noise* form.

## 5. CROSS_MODULE resolved on spelling, and this session walked into it

The new module defined `SELF` and never read it. The dead-name sweep did **not**
report it DEAD — it reported `CROSS_MODULE`, resolved against a
**function-local variable of the same spelling** in
`graphene_log_sensitivity_step_audit.py`, a module that does not import the new
one and could not reach that constant if it tried.

10-01 recorded this exact fault in the sweep ("its own `CLASS_LIVE` came back
CROSS_MODULE purely because the probe defines a constant with the same
spelling") and **fixed half of it**: the mention-vs-use half, by moving from
`re.search` to an AST pass. The half left standing is that an AST reference in
another module is a liveness claim *only if that module can reach this one*.

New **Section 2b** checks every CROSS_MODULE resolution against a real import,
with a positive control. Measured: 3 rows, 1 unsupported. The two supported ones
are real and untouched (`per_metal` imports `CHEMISORBED`; the sweep itself
imports `nl.D_CHEM_PHYS`). `SELF` was given a *reader* rather than a deletion —
a runtime assertion that the module is not listed in `SUITES`, because the
exclusion had been carried by a comment and a comment is not checkable.

## 6. Two stale transcripts that predated today

Neither visible to `git status`:

- **`figure_provenance_audit_output.txt`** — stale since 10-01. Still reported
  **37 frozen default-capture sites across 11 files** (the sites 10-01
  converted; the census now reads 0 of 36) and still showed the `T_HOP`
  mutation as frozen, `8.379013 -> 8.379013`, when it now arrives,
  `-> 16.758026`. **The transcript was a session out of date on exactly the
  fault 10-01 fixed.**
- **`log_sensitivity_step_audit_output.txt`** — stale since 09-28, when the
  covariance probe audit was added: **11 call sites → 18**, and
  discarded-convergence-flag comparisons **1 → 3**. A reader of the committed
  transcript would have concluded the repository was cleaner than it is.

## 7. A control of mine fired, again

`C2` bumped the last printed digit of a `%.6f` number and asserted the verdict
would be `NUMERIC_DRIFT`. For `7.500000 -> 7.500001` that is a relative change
of **1.3e-7**, a hundred times *outside* `NUMERIC_RTOL = 1e-9`, so `STALE` was
correct and the control's premise was false.

> **"The last digit" is a fact about a format string. The tolerance is a fact
> about the number.** A control must be built from the quantity it is testing.

The classifier was right about the real data — the drift it was built to catch
is 2.6e-16, seven orders of magnitude inside the tolerance — and only the
synthetic control was wrong. C2 now perturbs by `NUMERIC_RTOL/100` of the value,
and **C2b** asserts that the perturbation really is inside the tolerance, so the
fault cannot return.

## 8. The one FAIL that is kept deliberately

`graphene_rootfinder_audit.py` classifies `NUMERIC_DRIFT` at a worst relative
difference of **2.586e-16**: the root agrees to a few ULP, but the last three
printed digits differ by environment
(`4.99733219460098000e-11` committed against `...97870e-11` here; residuals
identical). So the check *"every transcript is exactly reproducible"* is **FALSE
for this repository**, and it should keep saying so. These transcripts are
reproducible **results**, not reproducible **artefacts**, and byte comparison
cannot certify them across library versions. Lowering the bar to make the suite
green would delete the only statement of that limitation.

---

## Methodological note, continuing the series

- **09-25** — a procedure asked whether it has converged can answer yes and be 44% wrong.
- **09-26** — an anchored comparison is anchored in one variable.
- **09-27** — an instrument can be systematically smallest where the answer is worst.
- **09-28** — prose is a detector, and a check that cannot fail is worse than no check.
- **09-29** — a mutation that does not arrive is indistinguishable from a system that does not respond.
- **09-30** — a control has to sit where the failure enters, not where it shows.
- **10-01** — and it has to name what it compares against in a way that cannot come to mean something else.
- **10-02** — **AND WHEN IT AGREES, THAT IS A FACT ABOUT TWO ARTEFACTS, NOT ABOUT THE WORLD.**

Today's instrument asks whether the code and the transcript agree, and the day's
worst finding sat inside the one suite that agreed **perfectly**. `RESULT 1` was
reproducible, pinned, byte-identical, and false. 10-01 established that a
reference must name what it points at in a way that cannot drift; 10-02 adds
that **a reference that cannot drift can still point at the wrong thing from the
day it was written**, and agreement between two artefacts is silent about which.

There is a second thread, and it runs through §2, §3 and §5 together:
**every one of the three was a repair that was built and then not connected.**
09-28 preserved the legacy path and left the validation pointed elsewhere.
09-28 fixed the zone-corner constant and left a duplicate of its value three
functions away. 10-01 diagnosed the spelling loophole in the sweep and fixed one
of its two halves. In each case the session that found the fault understood it
correctly, wrote the right mechanism, and **stopped one wiring step short** —
and in each case what remained was invisible precisely because the fix's own
write-up read as complete. The lesson is not to be more careful. It is that
**the last step of a repair is the check that the repair is reachable from where
the fault was**, and none of these three had one.
