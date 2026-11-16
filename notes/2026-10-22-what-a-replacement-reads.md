# What a replacement reads: the question was too narrow, and a capture is one level deep

**2026-10-22.** Ranked-first item from the 2026-10-21 entry: *does any
replacement read a name its own flip rebinds?* It was ranked first because the
2026-10-21 defect was found **by accident** — the two existing guards both
answer *which names move*, and a replacement that reads a moved name is
neither stale nor unassigned, so nothing in the repository would have asked.

Module: `graphene_replacement_read_audit.py`. Transcript:
`replacement_read_audit_output.txt`. Figure: `replacement_read_audit.png`.
30/30 validations, all against exactly known values, transcript byte-identical
across two consecutive runs. No live web search was used; the subject is this
repository's own code measured against itself.

---

## 1. The answer to the item as asked

Of the four replacements the item names, **exactly one** reads a name that is
rebound while it is installed:

| Replacement | own flip | the composition |
|---|---|---|
| `qca._carrier_density_consistent` | — | — |
| `qca._quantum_capacitance_consistent` | — | — |
| `blast._cq_consistent_contact` | **—** | **`gfet.quantum_capacitance`** |
| `flip.make_contact_cq(live=False)` | — | — |
| `flip.make_contact_cq(live=True)` | — | `gfet.quantum_capacitance` |

So 2026-10-21's diagnosis was **complete as well as correct**: the accident
found the whole of the exposure as the item framed it. That is worth having as
a measurement rather than a hope, and it is the less interesting half of the
session.

## 2. THE ITEM'S OWN WORDING CLEARS THE DEFECT THAT CREATED IT

The item, and §4.6.17(c) which it came from, both say: *which module-level
names does this function read at call time, and does **the flip** rebind any
of them?* "The flip" means the flip that installs the replacement. Asked that
way, for `blast._cq_consistent_contact`:

```
installed by      : blast.contact_on_consistent_argument()
that flip rebinds : [('graphene_contact_doping_model', 'quantum_capacitance')]
it reads          : [('graphene_fet_model', 'quantum_capacitance')]
intersection      : EMPTY
```

**Empty.** A guard written to the item's own wording would have passed the
composed flip and reported nothing.

The defect is **cross-flip**. The name read is rebound by a *different*
manager, and the two managers are concurrent only because item 16 composes
them. So the relation is not replacement × its own flip but

> replacement × every flip that can be installed **concurrently**

and the composed column is the one that answers the question. This is the
same shape as the repository's recorded dominant failure mode — a check
answering a weaker question than the one it was written to ask — reached this
time through the *scope of a quantifier* rather than through a ratio.

I have left §4.6.17(c)'s sentence unrevised and annotated it beside, because
it is the record of what was asked.

## 3. The instruments, and why there are three

| | mechanism | sees | blind to |
|---|---|---|---|
| **A** static | bytecode closure over live lookups, cells resolved from the real function object | reads on branches that do not execute; captured-vs-live | nothing dynamic; over-approximates |
| **B** sentinels | identity counting wrapper at each rebind-set member | reads spelled `mod.attr` **and** reads made from inside the owning module | non-callables; writes |
| **C** module type | `__getattribute__`/`__setattr__` on a `ModuleType` subclass | non-callables (`gfet.C_ox` is a float); **writes** | a read made from inside a name's own module — a globals-dict lookup never reaches the type |

B and C are blind in **opposite** directions, which is why asking them to
agree is a test. They agree exactly on all five replacements.

**The planted controls earned their place.** Three of the four controls are
flagged by B and invisible to C — each is a read of a rebound name made from
inside its own module. The five real replacements happen not to exercise that
path, so without the controls, C would have looked interchangeable with B and
the blind spot would have been invisible. A reader that agrees with another
reader only because neither was asked a question it could fail is not
corroborated.

Rebind sets are **derived from each manager's own `setattr` calls** (entering
it under reader C), not re-derived by calling `identity_class()` again. A
manager that grows a third name is reported with three without this module
being edited — which is the property the hand-maintained lists of October
did not have.

## 4. NEW: a write/read relation, outside every existing question

Reader C reports:

```
graphene_fet_model.C_ox
    in any rebind set? NO -- invisible to both guards
    WRITTEN by 3/5 replacements
    READ    by 5/5 replacements
```

`gfet.C_ox` is in nobody's identity class. It is therefore invisible to
`bc.rebind()`, to `bc.rebind_by_operating_point()`, **and** to this item's own
question, all three of which are about names a flip *rebinds*. The thing that
governs composition is not a set of names but a **bipartite relation between
what replacements write and what replacements read**, and 2026-10-21's
mechanism (ii) is one edge of it rather than a quirk of one function.

## 5. 2026-10-21's two mechanisms are ONE

Counting channel-replacement arrivals during the contact census, and recording
`gfet.C_ox` at each arrival:

| flip | arrivals | with `C_ox` already substituted |
|---|---|---|
| repaired (`live_contact=False`) | **0** | 0 |
| unrepaired (`live_contact=True`) | **192** | **181** |

2026-10-21 decomposed the failure into (i) the channel solve being reached at
all and (ii) that solve seeing `C_ox` substituted, measured them as partially
cancelling (+1.0204 to +2.2965 pp), and concluded that repairing the
conspicuous half alone would be a regression. The counts sharpen that: **(ii)
cannot occur without (i)**, because (ii) is a property of the context (i)'s
arrivals land in. One cause, two consequences, 181 of 192 arrivals carrying
both. That is *why* a single change at (i) repairs both — which the percentage
table could not have told us, and which I had been reading as a coincidence.

## 6. NEW: CAPTURE IS ONE LEVEL DEEP

This is the finding I did not expect and the one that bears on the repair.

A cell protects the **name** it holds. It does not protect the names the held
**object** reads at call time.

| captured object | its own live reads | in the rebind set? |
|---|---|---|
| `gfet.quantum_capacitance` | `{gfet._thermal_rounding}` | **no** — safe |
| `gfet.carrier_density` | `{gfet._thermal_rounding, gfet.quantum_capacitance}` | **yes** — not safe |

The two captures are the *same construction* and are indistinguishable at the
call site. `make_contact_cq(shipped, live=False)` and
`make_captured_control(gfet.carrier_density)` are written identically; one is
immune and one is not, and which it is depends entirely on what the captured
object reads. The planted captured-leaky control is still flagged by reader B,
which demonstrates the limit rather than arguing it.

So **§4.6.17's repair is correct, and its correctness is contingent.** Had the
contact needed the shipped *charge relation* instead of the shipped
*capacitance*, the same repair would have been written, would have looked
right, and would have left the defect in place.

Stated in the form the family's other members take:

- **2026-10-14** — a repair to one binding is not a repair to the class.
- **2026-10-17** — a class can span two operating points, so it must be
  partitioned.
- **2026-10-21** — a replacement that *reads* a rebound name is neither stale
  nor unassigned.
- **2026-10-22** — a capture that fixes one name does not fix the names the
  captured object reads. The safe criterion is the **transitive live-read set
  of the captured object**, which is what reader A computes.

## 7. The second limit: a default is not a removal

Reader A is branch-insensitive and reports `make_contact_cq(live=False)`
exactly as it reports `live=True`. That is not a false positive. The live read
is still *in* the function, on a branch a caller selects with a keyword, and
`make_contact_cq` is public. B and C, which execute, separate them; A, which
does not, is right that the defect is one argument away.

Reporting all three readers separately rather than reducing them to one
verdict is what makes that visible.

## 8. Defects in this session's own work

Three, all in the figure, all found by rendering it and looking at it:

1. **Panel (a)'s legend sat on top of the tallest bar.** Repaired by setting
   the y-limit before placing the legend rather than after.
2. **A FALSE LEGEND in panel (b).** The legend swatch for the "inside the
   `C_interface` window" overlay was drawn green; the overlay it labelled was
   drawn dark red. This is the second-kind figure defect 2026-10-21 named —
   a figure asserting something the data does not say — and it is the second
   consecutive session in which one has appeared. Repaired by building the
   legend from the **bar artists themselves** rather than from a hand-written
   swatch, so the two cannot drift.
3. **An exact zero drawn as a zero-height bar**, visually indistinguishable
   from missing data, in the panel whose entire point is that the repaired
   flip gives exactly 0. Marked explicitly now.

Third consecutive session in which inspecting a rendered figure found a defect
on the first try. The 2026-10-21 prescription — *any figure carrying a number
in its title or annotations should derive it from the computation* — was
followed here and held; the defect that got through was a **colour**, not a
number, which the prescription does not cover. Its generalisation: **any
property of a mark that the legend restates should be read from the mark.**

## 9. What this does not claim

It audits five replacements. The eleven surviving hand-listed rebinding sites
are not replacements and are not covered. It does not flip the default;
`graphene_fet_model.py` is unedited and the chapter default is still (4.29).
And the write/read relation of §4 is *reported*, not guarded: its reachability
under the repaired flip is currently zero, but that zero is a property of the
capture, not of the write.
