# The repair that refuted its own prediction, twice

**2026-10-27.** Session subject: the 2026-10-26 ranked-first item — fix the
five MIXED names and the twelve frozen/live capture sites its full-transcript
regeneration exposed. Live web search was **not used**: the subject is this
repository's own instruments measured against themselves, and there is no
external source for it. PubMed/PMC not attempted, per the standing note.

## 1. What was wrong

A module-level name is `LIVE` if every reader loads it from the module's
globals at call time, `FROZEN` if every reader captured it in a function's
`__defaults__`, and `MIXED` if some of each. The arrival probe's own
transcript calls `MIXED` *"the worst class, not a middle one"*: a rebind
reaches the live readers and not the captured ones, so where the live readers
are an expected-value expression and the captured ones are the measurement,
the mutation manufactures the disagreement the audit then reports.

Five names, twelve sites, three modules:

| module | name | before |
|---|---|---|
| `graphene_chord_decay_length.py` | `D1_CHORD` | 3 frozen / 9 live |
| | `D2_CHORD` | 3 / 4 |
| | `D_SPAN` | 1 / 3 |
| `graphene_capture_surface_guard.py` | `PERTURBATION` | 2 / 1 |
| | `PROBE_V` | 2 / 1 |
| `graphene_velocity_saturation_model.py` | `N_QUAD` | 0 / 2, 1 by scope |

`PERTURBATION` and `PROBE_V` are the two names **every counterfactual in this
repository is taken with**, so the split was a defect in the instrument, not
in any model the instrument measures.

Eleven sites converted to the 2026-09-30 sentinel form. Both affected
transcripts **bitwise unchanged**, which is the standard the 09-30 and 10-01
conversions were held to.

## 2. The prediction, and why it was wrong

`graphene_mixed_capture_repair_audit.py` imports the two **pre-repair blobs**
out of the object store by hash (`4ef1e3c8`, `8f74a649` — a blob hash, not
`HEAD~1`; a position is not a state) under prefixed module names, so both
states are live in one process.

Section M predicted the textbook `MIXED` signature on `D1_CHORD`: `L'(d)`
moves under a rebind, `⟨L'⟩` freezes. **It does not freeze.** `mean_dlogP`
captured `d1` as a default, but its body calls `poly_P`, and `poly_P` reads
`D1_CHORD` *live* to re-centre the polynomial at `d − D1_CHORD`. So under a
rebind the pre-repair code averaged `L'` over the **old** interval while
evaluating a polynomial centred on the **new** `d₁`:

```
d1: 3.3 -> 3.9 A,  shape (p1, p2) = (-0.30, 0.02),  d = 4.0 A

               default path      new-interval value   pre-rebind value
before   -3.156428052855e-01   -3.318137959586e-01  -3.540264836195e-01
                                 off by 1.617099e-02   off by 3.838368e-02
after    -3.318137959586e-01   -3.318137959586e-01  (bitwise equal)
```

The pre-repair number lies **between** the two coherent values and equals
neither. It corresponds to no consistent choice of `d₁` at all. That is
strictly worse than frozen: a frozen number is at least one the model once
produced. And the quantity carrying it is Eq. (D), the γ-free discrepancy
§7.9 item 2's annotation quotes.

## 3. The second wrong draft, which is the methodological point

Restating the check as **coherence** rather than movement caught a second
error of the same shape. The "old interval" reference was first computed as an
explicit call taken **with the rebind still in place** — so the reference read
`mean_dlogP(p1, p2, d1_old, d2)` against a polynomial already re-centred at
the new `D1_CHORD`. That is the same incoherent quantity the check was trying
to distinguish itself from. It agreed with the pre-repair default path to
`0.000000e+00` and the check failed, correctly.

**A reference value has to be taken in the state it claims to describe.** Both
wrong drafts are kept in the section's docstring; they are the finding, not
noise around it.

A third error, caught by the suite's own control: the pinned blobs were first
registered in `sys.modules` under their **real** module names, which shadowed
the repaired modules, so every "after" lookup returned the "before" module and
nine checks failed at once. `P2` (the two states must be distinct objects with
distinct `__file__`) and the new `P3` (`sys.modules` must still map the real
names to the repaired modules) are why that is now a check rather than a
memory.

## 4. The residual, and a conflict between two standing requirements

The twelfth site is `n_quad=N_QUAD`, a `def` inside a function **body**: the
default is re-evaluated on every call of the enclosing function, so a rebind
arrives. It is `live_by_scope`, not frozen, and
`graphene_frozen_default_semantics.py` measures it as live in its own Section 1
V2. It is also that suite's **only real** demonstration site — its Section 4
`D1`–`D5` exist to show that a static pass can find a capture its confirmation
step then silently discards, on a site that actually exists.

So the arrival probe's Section 7b splits:

- *"and in particular none is `MIXED`"* — **green for the first time since it
  was written on 2026-10-01**;
- *"no module-level constant is captured as a default anywhere"* — still red
  at one site-group, and **cannot** go green while that demonstration keeps
  its site.

The conflict is in 7b's **vocabulary**, not in the tree: it counts syntactic
captures and names the class it cares about only in the second check. Not
resolved here. The two repairs are (i) split 7b into FROZEN-or-MIXED
(achievable) and syntactic (not), or (ii) move the demonstration onto a
planted fixture — which trades a real site for a synthetic one, the trade that
suite was written to argue against. (i) should be costed first.

## 5. What the regeneration found that nobody was looking for

The 2026-10-25 standing item asks for **every** committed transcript to be
regenerated each session. Doing it:

- **`potential_declaration_output.txt` asserted a defect that had been
  repaired three days after the transcript was written, and carried the claim
  for nineteen days.** Written 2026-10-08: *"quantum_capacitance() returns a
  finite value up to dV = 18.3493 V and NaN/inf above ... It is still a defect
  ... left as a new item rather than patched today."* The overflow-free
  spelling was installed 2026-10-11 (`7019f61`, "Item 13 CLOSED"). Measured
  today: finite at 18.3, 18.4, 25, 40, 80, 200 V. 10-26 found transcripts that
  **understated** the tree; this one **overstated** it. Annotated in place —
  18.3493 V and η ≈ 710 are the pre-repair state and are kept.
  The bisection now returns its own search **ceiling**, and a bisection
  reporting its ceiling is not an onset, so `SEARCH_LO`/`SEARCH_HI` are named
  rather than written twice and a second check requires the two states to be
  distinguishable.
- **`import_binding_sweep_output.txt`**: parsed 56 modules against 71;
  from-imports 158 → 179, module aliases 70 → 124.
- **`potential_census_output.txt`**: 86 names across 25 modules against 133
  across 35; `COLLECTED_NO_ROLE` 60 → 99, eleven new unadjudicated `V_dirac`
  rows.
- **`frozen_default_semantics_output.txt`**: reading 58 python files against
  71 — it predates thirteen modules.

## 6. Three operational facts for the next session's first command

1. **Transcripts are stdout only.** Capturing stderr as well puts four numpy
   `RuntimeWarning` lines at the head of every transcript whose suite imports
   `graphene_quantum_capacitance_argument.py`, making three transcripts read
   STALE on environment noise alone. Same signature as the 2026-10-26
   `MPLCONFIGDIR` fault: lines that are not suite output but are
   indistinguishable from it to a byte comparison.
2. **Three modules write their own transcript into the working tree** —
   `graphene_import_binding_sweep.py`, `graphene_default_flip.py`,
   `graphene_contact_argument_blast_radius.py`. Running one in place
   overwrites the file a byte comparison is taken against, so the verdict is
   **circular** unless taken against the object store. That is how the
   import-binding staleness was caught and how it was nearly missed.
3. **Three suites cannot be regenerated in this environment at all.**
   `graphene_perfect_contact_mutation.py`,
   `graphene_velocity_saturation_model.py` and
   `graphene_velocity_saturation_mutation.py` each exceed 170 s, and a single
   shell call here is capped at 180 s with **no surviving background
   process** — `nohup` and `setsid` are both killed when the call returns.
   None is in `SUITES`. Their transcripts are left untouched rather than
   half-written.

Plus, in the repository's own source rather than the environment: two mutation
suites built their mutated copies under a **fixed** `/tmp` path another
session owned, so `rmtree(ignore_errors=True)` did nothing, `makedirs` on the
existing parent succeeded, and `copytree` died with `PermissionError` —
**after** the unmutated control had printed PASS. A mutation suite that cannot
create its working tree reports **no** score, which is not `0/9`. Both now use
`mkdtemp` under `TMPDIR`; both transcripts bitwise unchanged at 9/9.

## 7. A new entry point for the oldest standing item

Four committed PNGs do not re-render byte-identically even where their
transcripts do: `differential_crossover.png` (−38724 B),
`log_sensitivity_step_audit.png` (−80248 B),
`perfect_contact_counterfactual.png` (−26968 B), `import_binding_sweep.png`
(+5715 B). The first three are **reverted** — their transcripts are bitwise
identical, so there is no content change to track and the delta is a
rendering-environment difference. The fourth is committed, because its
transcript genuinely changed and figure and transcript should stay in step.

"Render and inspect every committed figure" has been carried for six sessions
with ~63 PNGs behind it and no entry point. This is one: four named files, a
measured byte delta each, and a hypothesis (matplotlib rendering drift, not
content) that a single rendered comparison would settle.
