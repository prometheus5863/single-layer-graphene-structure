# Three unchecked claims, and a guard that is blind to its own module

2026-10-28. Session subject: the 2026-10-27 ranked-first and ranked-second
items, worked in order — and the second of them caught the first committing
the fault it guards against, on the same day.

---

## 1. The ranked-first item: a red check that could not go green

Section 7b of `graphene_mutation_arrival_probe.py` carried two checks over one
census:

| check | vocabulary | state on 2026-10-27 |
|---|---|---|
| `no module-level constant is captured as a default anywhere` | counts **syntactic** capture sites | RED, 1 site-group |
| `and in particular none is MIXED` | names a **delivery class** | GREEN, first time since 2026-10-01 |

The 10-27 entry recorded why the first cannot go green: the single remaining
site is `def _R_harmonic(V_g, V_hi, n_quad=N_QUAD)` in
`graphene_velocity_saturation_model.py`, a nested `def` inside a function
*body*, so its default is re-evaluated per call of the enclosing function and
a rebind arrives. `classify()` files it under `live_by_scope`, not
`frozen_sites`. It is also the only real `live_by_scope` site in the tree, and
`graphene_frozen_default_semantics.py` takes its Section 4 `D1`–`D5`
demonstration from `live_rows[0]` — so removing the site makes those five
checks **stop running** rather than start failing.

A red check that cannot go green has stopped being a guard. Two properties it
had lost:

1. It is indistinguishable on sight from a red check nobody has read in a
   month.
2. Its report — `1 site-groups remaining` — is indistinguishable from *a new
   capture having appeared while the known one was removed*. That is the
   standing "pin a red suite's failing membership, not its count" item, inside
   a single check.

### The split

Option (i) of the two the 10-27 entry named. Three checks over the same
census:

- **7b-i SEMANTIC** — no capture site anywhere *freezes a read*. `FROZEN`,
  `MIXED` and `IMPORT_FAILED` all fail it, because an unimportable module is
  UNKNOWN and unknown is not green. **Takes no exemptions**:
  `census_semantic_faults()` does not accept the registry as a parameter, so
  no caller can pass one. GREEN — this is the statement 7b was written to
  make, stated so that it *can* be true.
- **7b-ii SYNTACTIC** — no capture site in the tree is unregistered. One
  registered entry, carrying its expected class and its reason in full. A new
  capture anywhere still fails the suite, which is what the old check was for.
- **7b-iii HYGIENE** — every registry entry still matches a real census row
  carrying the class it claims. `VANISHED` (site gone — delete the entry) and
  `DRIFTED` (class moved under the entry) are reported separately because they
  call for opposite repairs.

What the split gives up: the claim that the syntactic count is zero, which was
never true of this tree and is not going to become true. What it gains: three
statements that can each be false, where there was one that could not become
true and one that could.

Eight planted controls, `C7c`. Seven require a predicate to fire on a planted
census; the eighth requires one *not* to. The load-bearing arm is the seventh:

> the **registered** site, planted as `FROZEN`, must **still** fire 7b-i.

Without it, the split would be a candidate weakening with no evidence against
the obvious failure — a registry that silences the semantic check. The eighth
is the positive control without which arm 4 would pass for the wrong reason
(every row unregistered, exemption path never reached).

Suite 33/1 → **43/0**. The transcript diff is confined to the 7b block and the
`TOTAL` line.

---

## 2. And then the split committed the fault it guards against

`graphene_figure_provenance_audit.py`'s Section C runs a repository-wide
census of the inert-knob pattern. It went from **1 occurrence in 1 file to 4
in 2**. The three new ones were mine, written an hour earlier:

```python
def census_semantic_faults(rows, faults=SEMANTIC_CAPTURE_FAULTS):
def census_unregistered_sites(rows, registry=SYNTACTIC_CAPTURE_EXEMPTIONS):
def registry_hygiene_faults(rows, registry=SYNTACTIC_CAPTURE_EXEMPTIONS):
```

Module-level constants captured as function defaults — in the module that
implements the guard against exactly that. Converted to the 2026-09-30
late-bound sentinel form; suite 43/0 and its transcript **bitwise unchanged**,
the standard every earlier conversion was held to.

**Section 7b did not catch it, and cannot.** `repo_wide_frozen_census()` skips
four modules by name, and one of them is
`graphene_mutation_arrival_probe` — the module that defines the census. The
guard is blind to its own source. The exclusion is defensible (the probe is
full of deliberately planted capture patterns for its own controls) but its
consequence had never been stated: *the only reason today's three were found
on the day they were written is that a different suite censuses the whole
tree, including this one.* New standing item.

It is also the second instance today of the §4.6.8 shape — knowing a failure
mode is not the same as not committing it — and this one is sharper than the
`fmax` instance it echoes, because the module committing the capture is the
one that defines the class.

---

## 3. The ranked-second item: the queue was nine larger than the estimate

10-27 wrote "18 rows against roughly 40 committed transcripts". Counted
mechanically today:

```
committed *_output.txt   49
compared in SUITES       17
```

The estimate understated the queue by nine, which is the argument for counting
it rather than eyeballing it. **SUITES 17 → 35 compared transcripts**, every
added row measured first (run from a clean tree at HEAD, compared
byte-for-byte) with its measured seconds in the comment beside it.

Three things fell out of the counting that a row-by-row pass would not have
produced.

### 3.1 A row that said `None` for a transcript that exists

`graphene_covariance_probe_audit.py` was registered with `transcript=None` —
"no committed transcript to compare" — from the day the table was written.
Its transcript is `default_scale_probe_audit_output.txt`, whose filename does
not contain the module's name, which is why nobody paired them. It is
byte-identical to the module's stdout today.

A `None` row is **worse than an absent row**. An absent row is a gap; this one
printed `RAN` every session and looked like coverage.

### 3.2 A fifth stale transcript, stale by a month, overstating the tree

`audit_output.txt` (`graphene_sensitivity_audit.py`) claims

```
  (c) exact zeros : 3/3  S(rho_bulk) = 0.0, S(p_default) = 0.0, S(lambda_bulk) = 0.0
```

The code has printed

```
  (c) BITWISE zeros : 2/3   zeros within the derived floor 5.7e-11 : 3/3
      S(rho_bulk) = -1.110e-11, ...
```

since **2026-09-27**, when the stepping was corrected to symmetric — and the
source's own annotation says one of the three exact zeros was floating-point
luck. Same direction as the nineteen-day false claim 10-27 found: the
transcript **overstates** the tree, which sends a reader to trust a bitwise
result that is not bitwise. Verified deterministic run-to-run before
registering, then regenerated.

### 3.3 Section 5: the gap is now a number the suite prints

```
  committed *_output.txt    49
  registered in SUITES      35
  excluded with a reason    14
      NOT_A_STDOUT    1
      SELF            1
      SELF_WRITING    9
      UNMEASURED      3
  UNACCOUNTED                0
```

The suite fails if any committed transcript is in neither table. Two of the
four exclusion classes are **checked rather than declared**:

- `SELF_WRITING` (9) — the module opens that transcript `'w'` in its own
  working tree, so the pristine comparison would read the file the run just
  wrote. Circular. `writes_own_transcript()` confirms all nine from the
  source, and `C5.2` is the negative control on a module that does not.
- `SELF` (1) — this module; `main()` already asserts it.
- `NOT_A_STDOUT` (1) — `sensitivity_converged_output.txt` is a hand-written
  excerpt with its own header and date. Byte comparison is the wrong
  instrument and pretending otherwise would make the suite permanently red for
  no finding.
- `UNMEASURED` (3) — the honest name for next session's work.

Plus the hygiene check today's 7b work argues for: no exclusion may name an
uncommitted transcript, and none may be excluded *and* registered.

Pristine audit 13/0 → **19/0**, 35 of 35 IDENTICAL, 0 STALE, 0 NUMERIC_DRIFT.
103 s wall against a 180 s shell cap, so 77 s of headroom — measured, not
estimated, and the three `UNMEASURED` exclusions are mutation suites that
could each be large.

---

## 4. The third unchecked claim, found by accident

Registering the two names today's first commit introduced raised the question
that produced the day's best finding.

The 7b exemption carries its reason in full, and the reason spells the
signature it is about: ``def _R_harmonic(V_g, V_hi, n_quad=N_QUAD)``.
`graphene_potential_census_audit.py`'s `_string_names_in()` collects
potential-like tokens out of string constants **by design** — the census is
"deliberately over-wide; every match must be adjudicated" — so two names in
**prose** entered it, and the census moved from 133 names across 35 modules to
135 across 36 *within one commit of the same session*. Nothing would have
caught that: `potential_census_output.txt` was not in `SUITES` either. Today's
own first commit is therefore the ranked-second item's argument, restated.

Both names are `DISPLAY_ONLY`. Which raised: **what checks a `DISPLAY_ONLY`
row?** Nothing, from 2026-10-07 until today. An adjudication nobody checks is
the same hole as an exemption nobody checks — the other half of this session
— so it got the same treatment.

**One false row in thirteen.**
`graphene_photodetector_collection_model.py::V_bias` read `DISPLAY_ONLY,
imported and printed only`. Line 72 is

```python
E_BIAS = V_bias / L_channel  # V/m, uniform bias-driven drift field
```

— the field its own docstring describes, and the `TRANSPORT_MARKER_RE` shape.
Re-adjudicated `TRANSPORT_ONLY`, the reading its sibling
`graphene_photodetector_model.py::V_bias` already carries for the same stated
reason (no charge relation exists anywhere in Chapter 6, so a transport-only
reading carries no channel ambiguity). The old reading is kept above the new
one.

**And the first form of the check was wrong, which is worth more than the
row.** Comparing the scopes the census already records reported **two**
violations. The second,
`graphene_photodetector_two_contact_model.py::V_bias`, is not one: its uses
are `{V_bias}` inside f-strings, and an f-string interpolation is a real
`ast.Name` inside the enclosing function, so the census files it under
`make_plots()` rather than under `<format string>`. **The scope vocabulary
cannot express "displayed inside a function"** — the same
two-vocabularies-over-one-census fault the 7b split repaired this morning, in
a third place on the same day. The shipped form is lexical: inside a
`JoinedStr` or a string constant is display, an `ast.alias` is provenance,
everything else is computable.

**And the planted controls' first draft failed too**, and that failure was the
control's: both planted modules opened with `V_probe = 1.0`, and an assignment
target is an `ast.Name` in `Store` context. Every real `DISPLAY_ONLY` row is
either a LaTeX fragment that is no Python name at all or a name imported from
elsewhere, so the controls now have that shape — and the behaviour that
tripped them is pinned as `P-DO3` rather than patched away: a module-level
binding is a value and not a label, so the predicate reports it. No row in the
tree has that shape today.

14/2 → **18/2**, the same two failures (47 unregistered names; one EXPECTED),
and `UNREGISTERED` is 47 before and after, because the two names this session
introduced are adjudicated by the session that introduced them.

---

## 5. What the day was actually about

Three claims, in three unrelated vocabularies, none of them checked:

| claim | lived unchecked | how it failed |
|---|---|---|
| 7b's syntactic count | 2026-10-01 → today | could not go green; its count hid membership |
| a `SUITES` row's `None` | since the table was written | printed RAN, looked like coverage |
| a `DISPLAY_ONLY` adjudication | 2026-10-07 → today | 1 of 13 was false |

In each case the repair is the same shape: *a registry entry carries its
reason, and a check confirms the reason is still true of the source.* That is
now implemented three times — `7b-iii`, `writes_own_transcript()`,
`display_only_violations()` — and all three were written today, which is the
first time this repository has had the pattern rather than instances of it.

And the counter-lesson, from §2: **a guard excluded from its own census is
blind in exactly the place where it is easiest to be wrong**, because the
module that defines a fault class is the module most likely to be editing that
class's shape.

## 6. Environment

`/sessions` is **100 % full** — 9.3 G of 9.8 G, owned by other sessions,
unreadable and undeletable — which is the fourth consecutive session to record
it. A clone into `$HOME` fails with `could not create work tree dir ... No
space left on device`, so the clone goes into `/tmp` on `/dev/sda1` (4.0 G
free), as the 10-25 through 10-27 entries ask. New this session: `pip` also
fails there, because its build and cache directories default under `$HOME`;
`TMPDIR=/tmp/pt PIP_CACHE_DIR=/tmp/pipcache --no-cache-dir` installs `scipy`
1.15.3 into `/tmp/pylibs`. `/tmp/pylibs` did **not** survive from 10-27, so
`/tmp` is per-VM and not per-machine, and the install is a cost every session
pays.

`numpy` 2.2.6, `matplotlib` 3.10.9, `MPLCONFIGDIR=/tmp/mplcache_1028` (a fresh
name, per the 10-27 note that the old one is owned by another session).

**Live web search not used.** Today's subject is this repository's own
instruments measured against themselves, and there is no external source for
it. PubMed/PMC not attempted, per the standing note.
