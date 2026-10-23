# Two findings, one relation — and the relation now refuses

*2026-10-23. Discharges the 2026-10-22 ranked-first item: guard the
write/read relation, do not just report it.*

## What the item asked, and what it turned out to be

2026-10-22 §4 ended with a conclusion that was wider than the session that
produced it: what governs composition is **not a set of names** but a
bipartite relation between what replacements write and what replacements
read. `gfet.C_ox` was the instance — written by three of the five
replacements, read by all five, and in no identity class, so invisible to
`bc.rebind()`, to `bc.rebind_by_operating_point()`, and to that session's own
question, all three of which are about names a flip *rebinds*.

The item asked for that relation to become a check that fails. It is
`graphene_write_read_guard.py`, and building it turned the relation into
something the two sessions could not have seen separately.

## The vocabulary is deliberately not new

This repository already has the refusal twice:

- `bc.rebind()` refuses an exclusion **without a stated reason**, because an
  unstated exclusion is how all three of the 2026-10 delivery faults were
  introduced.
- `bc.rebind_by_operating_point()` refuses an **unassigned class member**,
  raises before anything is rebound, and makes the reason mandatory for
  "left alone" as well as for "replaced".

The guard is that refusal applied to an **edge** rather than to a name. Each
`(writer, name, reader)` edge must carry a verdict from
`{GOVERNED, IMPOSSIBLE, OCCURS}` and a mandatory reason; an unadjudicated
edge raises `UnadjudicatedWriteReadEdge` and nothing runs.

## What *is* new: every verdict is measured, and in both directions

Completeness is static — `rebind_by_operating_point` can check it because it
knows the identity class. A verdict is a claim about **reachability**, which
no static set contains. A guard that accepted the claim would be a census
with a stricter tone, which is exactly the standing item
(*a census that enumerates a fault forever is not a guard*, carried since
2026-10-17) that this work was supposed to discharge.

So each verdict is checked against an order probe:

| verdict | refused when |
|---|---|
| `GOVERNED` | the name is in **no** rebind set |
| `IMPOSSIBLE` | the probe **observes** the edge |
| `OCCURS` | the probe **never** observes the edge |

The third row is the half a census cannot have, and it is the half this
repository's recorded dominant failure mode calls for. A guard that only
catches a wrongly-*permissive* claim answers a weaker question than it was
written to ask. An `OCCURS` reason describing an edge that has stopped
occurring is a sentence about nothing, and it stays true-looking forever
because nobody ever asks it the opposite question.

## THE RESULT: the two findings are two edges of one relation

A name can be written by a **replacement** (10-22's `gfet.C_ox`) or by the
**flip itself**, which is every member of the rebind set (10-21's "a
replacement that reads a rebound name is neither stale nor unassigned").
Those were recorded as two findings in two sessions. They are two *sources of
writes* into one relation, and composing them is what makes the composed
column — the one 10-22 said answers the question — computable rather than
argued.

Measured, with reads traced and the edge derived from the read:

**Unrepaired (`live_contact=True`) — 3 observed edges**

```
THE FLIP          --gfet.quantum_capacitance--> make_contact_cq(live=True)
make_contact_cq   --gfet.C_ox--> qca._quantum_capacitance_consistent
make_contact_cq   --gfet.C_ox--> itself
```

**Repaired (`live_contact=False`) — 1 observed edge**

```
make_contact_cq   --gfet.C_ox--> itself
```

The first unrepaired edge **is** 2026-10-21's defect. The second **is**
2026-10-22 §§4–5. The survivor is the contact substituting `gfet.C_ox` for
`cd.C_interface` and reading it back inside its own call — the substitution
itself, same operating point, restored before return — so it is not a leak.

**One change removes both cross-replacement edges.** That is §5's "one cause,
two consequences" as a measurement rather than as a reading of a coincidence,
and it is why the percentage table of 10-21 was being misread: the 181/192
was never two phenomena.

**V21 is the item, answered.** The repaired composition's adjudication — the
same set of stated reasons, unchanged — is **accepted** for
`live_contact=False` and **refused** for `live_contact=True`, naming three
disagreeing verdicts. The difference between the two compositions is now a
raise instead of a table a reader has to interpret.

## Anchored to exactly-known values

- **V4–V6** reproduce the 2026-10-22 arrival counts rather than citing them:
  **0** arrivals repaired, **192** unrepaired, exactly **181** of those with
  `C_ox` already substituted. An exact zero a symmetry requires, plus two
  integers a previous session measured independently of this code.
- **V8** EXACT ZERO: no cross-replacement edge survives the repair.
- **V9** the unrepaired composition has exactly **2**, and they are the two
  previously reported defects.
- **V7** every observed edge lies inside the static relation, so the static
  set is a sound over-approximation and not a different answer.
- **V2–V3** the probe is deterministic: census and edge set bitwise identical
  across consecutive runs, checked **before** any count was used.

## The deriver is tested against known traces, not only against the repo

**V10–V13.** A deriver that returns few edges and a repository that *has* few
are indistinguishable until the deriver is asked a question whose answer is
known without it. So the replay is a separate function and is fed synthetic
traces: write-then-read yields exactly the one edge; read-**before**-write
yields none, which is what makes it an ordered probe rather than a set
intersection; a nested reader is attributed to itself and not to its caller;
a write outside every replacement is attributed to `THE FLIP`.

## Two defects in this session's own work

**1. The first probe did not trace reads.** It closed an edge on every
*arrival* of a reader that followed a write. That asks *did this replacement
run after that name moved* — weaker than *did it read that name after it
moved*, and wrong in the permissive direction. It reported **twelve** observed
edges under the unrepaired flip, four of them on names the reader does not
read at all, and the observed set was **not even a subset** of the relation.
A probe whose over-report has to be corrected by intersecting it with a
static set has not measured reachability. Reads are traced now and the edge
is derived *from the read*.

Same shape as the recorded dominant failure mode, reached this time through
**the choice of which events to instrument** — not through a ratio (10-22),
not through the scope of a quantifier (10-21), but through an instrument that
was never asked whether it could see the thing its verdict named.

**2. A planted control found a crash in a refusal path that had never run.**
V18 — a false `GOVERNED` — raised `TypeError` instead of `FalseEdgeVerdict`.
A name is a `(module, attr)` **tuple** and the message interpolated it into a
single `%s`. Nothing else in the module reaches that branch. Without the
control the guard would have shipped with one of its three refusals unable to
report itself, and it would have looked fine, because the two refusals that
*are* exercised by the real composition both work.

The generalisation, which is the fourth member of the 10-14/10-17/10-21/10-22
family: **a refusal that nothing reaches is not a refusal.** Every verdict in
a vocabulary needs a planted control, not just the verdicts the live subject
happens to exercise.

## No figure, stated rather than silent

The result is a 3-edge set and a 1-edge set. A table carries that exactly and
a plot would decorate it. The last three sessions each found a defect in a
committed figure, two of them figures asserting something the data did not
say; there is no reason to add a fourth opportunity for a result this small.

Machine output: `write_read_guard_output.txt` (21/21).
