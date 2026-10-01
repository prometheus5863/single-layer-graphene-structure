# A citation is a claim, and all three faults today were in the naming

**2026-10-01.** The 2026-09-30 top item closes: the last 23 captured-default
sites in the seven model modules are converted, and the delivery census that
measured them reports empty. Closing it turned up three instruments that could
not fail, and they are the session's real content. None of the three had a
wrong comparison. All three had a comparison against the wrong *name*.

## 1. The top item, and two corrections to how it was stated

09-30 filed it as "the 10 remaining knobs in 6 MODEL modules, 7 of them
MIXED". Counting rather than re-reading it:

- It is **seven** modules. `graphene_differential_crossover_model`'s
  `U_SMALL` was in the census listing and not in the module count.
- It is **23 sites**, not 10. Ten was the number of *names*. `W_CROSS_CHEM`
  alone is captured at 8 sites in one module and `ELL_DEFAULT` at 4 in
  another — which is exactly why 09-30 called those two the worst, since the
  9-frozen-against-2-live and 4-against-1 ratios it quoted were *site*
  counts. The fix is per site.

All 23 converted to late binding under a bitwise no-change requirement; all
seven module transcripts compare identical (7376, 5909, 7153, 3434, 1206,
3728, 2728 bytes). The hazard was prospective rather than live — nothing
currently mutates these names, so no committed number was wrong. What they
would have done is start the *next* mutation-based result built on any of them
at measuring nothing, which is why the item outranked the physics.

## 2. `HEAD~2` is not a reference to a state

With the model modules fixed, the probe was re-run and reported **two
failures** in Section 6, the before/after measurement of `T_HOP`'s
manufactured disagreement — the centrepiece of 09-30's own fix. It read its
`before` source with

```
git show HEAD~2:graphene_band_structure_audit.py
```

`HEAD~2` names a **position** in history, not a state of the repository, and
positions move. On 09-30 it resolved to the pre-fix file, because the fix was
still in the working tree and two commits had landed since the 09-29 log. The
moment the fix was itself committed, `HEAD~2` came to name the **post-fix**
file — so `before` and `after` became the same source, the oracle/measurement
split collapsed to two equal numbers, and both checks have failed ever since.

**Verified, not argued:** a pristine clone checked out at `c99c48d` — the
09-30 state exactly as committed — runs the probe to **27 passed, 2 failed**.
The 09-30 entry's "Validations: 29/29" is true of a working tree that was
never committed and false of every clone of this repository since.

The fix is to pin the reference by **content**: blob
`21a7e2ad96dd07f79d56807ca1f42eaccb8ecc29`. A blob hash is content-addressed,
so it cannot come to mean a different file — the one property the reference
needed and did not have. Section 8 is the standing guard: no source file in
this repository may reference `HEAD`, `HEAD~n`, `HEAD^` or `ORIG_HEAD`, because
fixing one site does not stop the next being written the same way.

And a **positive control on the reference**, which is 09-30's rule applied one
level out: the fetched blob must actually contain the capture the fix removed,
and must differ from the file on disk. **Its first form fired, on me.** It
looked for `'T_HOP=T_HOP'`; the signature is `def bands(k, t=T_HOP)` — the
parameter is not named after the constant. Recorded rather than quietly
re-tuned, because a control that has never fired is a control with no evidence
that it can.

The test that matters is whether the fix survives its own commit. A pristine
clone of the new HEAD runs the probe to **32 passed, 0 failed**. That is the
test the 09-30 version failed.

## 3. A citation is a claim, and nothing was checking it

09-30 asked whether there are **other** dead names here, because `ALPHA_ABS`
was found dead only because a mutation happened to point at it. Nothing was
looking for the class. `graphene_dead_name_sweep.py` looks: 113 module-level
constants, **112 LIVE, 1 BODY_ONLY, 2 CROSS_MODULE, 1 EXEMPT, 0 DEAD** after
four dispositions.

A dead constant is worse than an unused variable because every one of these
carries a **citation**. `PASSI_RC_DIRAC_OHM_UM` is six measured numbers from a
named table in a named paper. Writing the name is a claim that a result depends
on the number — and if nothing reads it the claim is false in the direction
that flatters the work, because a reader counts the citation as a dependency
honoured.

**The four, and what each turned out to be:**

1. **`PASSI_RC_DIRAC_OHM_UM`** — the second column of Passi et al.'s TLM
   table, whose own comment says the table is "on-state (V_BG = −40 V) **and**
   Dirac-point contact resistance". Only the on-state column was ever read.
   This is the session's physics result and it has its own section below.
2. **`D_CHEM_PHYS` = 0.9 eV** is a **second definition** of `DC_ANCHOR` =
   0.9 eV — one literature value, two copies, one unread. The physics is
   unaffected: the chemical term is applied through `w_cross_for_metal` →
   `delta_c` → `DC_ANCHOR`. That is exactly why it sat unread for ten
   sessions, and exactly why a revised extraction would move one copy
   silently. The two modules cannot import each other, so the agreement check
   lives in the sweep — the only file importing both, and the first code to
   read `D_CHEM_PHYS` at all.
3. **`DELIVERED`** in the probe — defined for symmetry with `UNDELIVERED` and
   `PARTIAL`, read by nothing. Given a reader: the census's conclusion stated
   *positively*, so a verdict class added later cannot slip past all three
   checks.
4. **`FANG_ANTENNA_TEST_WAVELENGTH_NM`** is **deliberately** unread; its own
   comment says the Fang et al. point is descriptive-only. So the exemption
   mechanism is a source-level `# dead-name-exempt: <reason>` marker, with the
   reason **required** — a reasonless marker is reported as DEAD. The
   exemption lives in the source and not in the checker, because a checker
   that reports a documented decision as a fault is a checker that gets
   ignored, and that failure mode is silent.

## 4. Two faults in the sweep itself, both found by running it

**Its cross-module check was a text search, and its own docstring exercised
the loophole.** `re.search` over another module's source cannot tell a **use**
from a **mention**, so the paragraph above explaining why
`FANG_ANTENNA_TEST_WAVELENGTH_NM` is exempt was enough to make it look alive.
Worse: this module's own `CLASS_LIVE` came back CROSS_MODULE purely because
`graphene_mutation_arrival_probe.py` happens to define a constant with the
same spelling. It was **DEAD in here**, and the loophole hid it. A dead-name
detector whose own dead name is concealed by its own mechanism is the 09-28
cannot-fail class, and it was one run from being committed. It is now an AST
pass over `Name` loads, `Attribute` names and `ImportFrom` aliases.

**Its exemption scope was three lines**, and the only real exemption in this
repository needs a four-line reason — so it reported a declared exemption as
DEAD. Found by running it, not by reading it.

## 5. The physics: the ~11× was compared against the wrong column

Reading the Dirac-point column gives a quantity the repository did not have.
**Rc(Dirac)/Rc(on-state)** is the *gate-tunable multiple* of the contact
resistance: near 1 would mean contact resistance is set at the interface and
is indifferent to the carrier density beneath it; large means it is dominated
by graphene's own density of states in and near the contact, which the gate
controls.

| hole D (nm) | on-state | Dirac | ratio |
|---|---|---|---|
| 0 | 519 | 1372 | 2.64 |
| 50 | 212 | 620 | 2.92 |
| 100 | 352 | 732 | 2.08 |
| **200** | **45** | **456** | **10.13** |
| 500 | 410 | 1354 | 3.30 |
| 1000 | 560 | 1590 | 2.84 |

2.08–3.30 at five of six diameters; **10.13** at D = 200 nm — the one
diameter whose 11× Chapter 4 quotes, **3.7× the mean of the other five**.

So **the headline reduction is a gate-bias statement, not a geometry
statement**: D = 0 → D = 200 nm is **11.53×** in the on state and **3.01×** at
the Dirac point, same two devices, same etched geometry.

And 3.01× is the comparison the model's own physics selects. `R_extra` is
computed from the *contact-induced* carrier density alone and the gate does
not enter it; at the Dirac point the gate contributes no channel carriers, so
the measured access resistance is dominated by exactly the term the model
isolates. **Against that column the model's 2.59× is within 14% of measured,
where against the on-state column it was short by a factor of 4.5.** Chapter 4
Section 4.8 held a gate-independent model against an on-state measurement,
which is a category error; 4.8.2 annotates it in place and keeps every old
number. Caveats kept rather than absorbed: the 2.59× is a pure-mode ratio
against a device that is a mixture at an unreported fill fraction, and Section
4.7's absolute-magnitude discrepancy is untouched.

**And 4.8.1's open item is sharpened rather than closed.** Whatever produces
the D = 200 nm optimum is **3.5× more of an effect in the on state than at the
Dirac point** (4.71× against 1.36×, taking D = 50 nm as reference), so it
scales with **carrier density**, not interface area. A model built entirely
from interface geometry — which 4.8.1's is — cannot produce a minimum whose
depth depends on the gate. Its failure was the wrong *class* of model, not a
missing refinement, and the evidence sat in the column the chapter cited and
did not read.

## 6. The methodological note

- **09-25** a procedure asked whether it has converged can answer yes and be
  44% wrong. **09-26** an anchored comparison is anchored in one variable.
  **09-27** an instrument can be systematically smallest where the answer is
  worst. **09-28** prose is a detector, and a check that cannot fail is worse
  than no check. **09-29** a mutation that does not arrive is indistinguishable
  from a system that does not respond. **09-30** a control has to sit where the
  failure enters, not where it shows.

- **10-01: AND IT HAS TO NAME WHAT IT COMPARES AGAINST IN A WAY THAT CANNOT
  COME TO MEAN SOMETHING ELSE.** Three instruments failed today and not one of
  them had a wrong comparison. The before/after measurement compared correctly
  — against a revision expression that came to denote a different file. The
  cross-module check compared correctly — against a name resolved by text, so
  a mention counted as a use. The exemption check compared correctly — against
  a marker whose scope was three lines and whose reason was four. **Every one
  is a naming failure, and a naming failure is invisible at the comparison,
  because the comparison is doing exactly what it says.** 09-30 put the
  control where the failure enters; today's lesson is that *where* is only
  half of it, and the other half is that the thing being pointed at has to stay
  the thing being pointed at.

- **The uncomfortable corollary.** Of today's three, the one that mattered most
  was in the instrument 09-30 built *to enforce* this series' own rules, and it
  broke by being committed. It worked in the run that wrote it and in no run
  afterwards, which is the worst possible failure schedule: it passes while
  under scrutiny and fails only once nobody is looking. **A check validated
  only in the session that wrote it has not been validated**, and the cheapest
  possible test for it — run the suite from a pristine clone of the commit —
  is a test this repository has never run. It ran today, twice, and found the
  fault immediately both times.
