# 2026-10-24 — The capture, and the name with the same leverage

Session subject: the 2026-10-23 ranked-first item, *audit the captures, not
just the reads*. Discharged, with a result that reverses the reading of
§4.6.18 (d) without changing any of its numbers.

Artefacts: `graphene_capture_surface_guard.py` (29/29),
`graphene_capture_surface_mutation.py` (7/7 injected defects detected),
`graphene_capture_leverage_figure.py` (15/15) and
`capture_surface_leverage.png`.

---

## 1. The instrument the item proposed cannot answer the question

The standing item named its own cheapest entry: `rra.static_reads` on each
captured object. Run on the two objects the whole 2026-10-21 repair is the
difference between:

| | static live-read set |
|---|---|
| `flip.make_contact_cq(shipped_cq, live=False)` — repaired | 8 names |
| `flip.make_contact_cq(shipped_cq, live=True)` — unrepaired | the same 8 |

Symmetric difference **zero**. Both sets contain `gfet.quantum_capacitance`
as a *live* read, because `cq_of_n` carries both branches in its bytecode and
a bytecode walk cannot know which one runs. `static_reads` is sound and
over-approximate by construction, so this is not a defect in it. What it
means is that §4.6.17's correctness argument — *bind to the captured object,
therefore the rebind cannot reach it* — is written in a vocabulary that
cannot distinguish the repaired construction from the unrepaired one. The
argument was never checkable where it was made.

So every verdict in the new module is a **counterfactual**: rebind the name
to a deliberately wrong stand-in, call the holder, compare **bitwise**.
Nothing is compared to a tolerance; the question *does a rebind of this name
reach this object* has an exact answer.

## 2. The finding

Perturbation ×1.5, contact replacement at `V = 0.5 V`, base
`0x1.0c835a6063492p-4` F/m²:

| name on the capture surface | verdict | repaired | unrepaired |
|---|---|---|---|
| `gfet._thermal_rounding` | EXPOSED | **+58.3421 %** | +58.3421 % |
| `qca.A_FERMI` | EXPOSED | +58.3417 % | +58.3417 % |
| `qca.series_factor` | EXPOSED | +26.1720 % | +26.1720 % |
| `qca._fixed_point` | EXPOSED | +22.4740 % | +22.4740 % |
| `cd.C_interface` | EXPOSED | +18.1154 % | +18.1154 % |
| `gfet.T` | EXPOSED | +0.0034 % | +0.0034 % |
| `cd.quantum_capacitance` | UNREACHED | exact 0 | exact 0 |
| `gfet.C_ox` | PROTECTED | **exact 0** | exact 0 |
| `gfet.quantum_capacitance` | PROTECTED | **exact 0** | **+58.3421 %** |

The two replacements differ on **exactly one** of nine names, and it is the
one the repair is about. That is the repair, measured rather than argued.

**And the number is not unique to it.** A rebind of `gfet._thermal_rounding`
moves the *repaired* replacement to `0x1.a92b527d6852cp-4`, which is
**bitwise** the value a rebind of `gfet.quantum_capacitance` moves the
*unrepaired* one to. `gfet._thermal_rounding` is the **sole** live read of
the captured object, so the route is the capture itself: the 2026-10-21
repair closed one path to that number and left an exactly equally powerful
one open *through the object it captured*.

`gfet._thermal_rounding` is in no flip's rebind set and is written by no
replacement, so `bc.rebind_by_operating_point`, the 2026-10-23 write/read
guard and §4.6.17's argument reach it equally: not at all.

### What this does to §4.6.18 (d)

That table's row — `gfet.quantum_capacitance`, own live reads
`{gfet._thermal_rounding}`, *in the rebind set?* **no — safe to capture** —
is **correct in every entry and wrong in its conclusion**. Not being in a
rebind set is not safety; it is the absence of a mechanism that would notice.
The column answers *is there a flip over this name today*, and the row was
read as answering *can this name move the result*. It can, by 58.3421 %, and
by exactly as much as the name the capture protects. The table's numbers
stand; §4.6.20 annotates the reading in place.

## 3. A second result, in the opposite direction

Walking the holder's closure cells on **object identity** recovers not one
held name but the held object's whole identity class —
`{gfet.quantum_capacitance, cd.quantum_capacitance}` — reproducing
`bc.identity_class(gfet, 'quantum_capacitance')` exactly and independently
(C5b, two derivations of one set, which is this session's validation against
the existing implementation in the overlapping limit).

So the capture neutralises a rebind of **both** members, which is a
**stronger** protection than §4.6.17 claimed; `cd.quantum_capacitance` is in
`blast.contact_on_consistent_argument()`'s rebind set and the capture
neutralises that flip too. But only one member is ever looked up, so
reporting both as PROTECTED would credit the capture with work it did not do
— the same conflation as §4.6.18 (d) with the sign reversed. `UNREACHED`
exists to keep the two apart, and separating them needs the one judgement
`static_reads` **can** make soundly: an over-approximate read set has
reliable *negatives*, so *not on the read surface* is trustworthy even though
*on it* is not. Section 0's instrument is kept for that and nothing else.

## 4. An argument from yesterday, a bitwise zero today

The 2026-10-23 guard left exactly one `OCCURS` edge standing — the contact
replacement's write-then-read of `gfet.C_ox` within its own call —
adjudicated *the substitution itself and not a leak* on the argument that the
write and the read are the same operating point. Here the **external** value
of `gfet.C_ox` is shown to reach the replacement's output **not at all**, for
both variants, because the replacement substitutes `cd.C_interface` for the
whole solve and restores the old value before returning. The same claim, now
a measured exact zero.

This is also a device-physics statement worth having in Chapter 4: the
contact self-consistency is referenced to `C_interface` and is **exactly**
independent of the oxide capacitance, so nothing in §4.6's contact numbers
inherits an oxide-thickness error through that path.

## 5. Defects in this session's own work

Four, all in the figure's provenance checks, and each found by the check
failing on its own source or by rendering the figure — none by inspection of
the code.

1. **The literal-number check was a regex over source text.** It reported
   three violations that were *code*: a regex has no notion of where a string
   literal ends, so the code *between* two literals was scanned as if it were
   one. **Resolve-on-spelling, running total fifteen** across both
   repositories, and the dominant failure mode at once — an instrument
   reporting on something it cannot see. Replaced with `ast`.
2. **The `ast` version walked every argument of the call**, and flagged
   `color='0.35'` inside `arrowprops`. A colour is not drawn as a glyph.
   `TEXT_ARGS` now names which argument of each call is the text.
3. **It then flagged the precision in `'%+.4f %%'`.** A digit inside a
   `%`-conversion specifier is a formatting instruction, not a number a
   reader sees.
4. **The figure was committed in its third render.** Renders 1 and 2 each
   passed every check they had:
   - render 1, linear axis: `gfet.T`'s +0.0034 % bar was indistinguishable
     from the exact zeros, so the *one* distinction the figure exists to make
     was not visible in it. `F3` had passed — it asked whether the zeros were
     a different **artist**, not whether a reader could see the difference.
     It also clipped the end of its own title and x-axis label.
   - render 2: the legend box sat on top of the `quantum_capacitance` row,
     the single row the figure is about. `F2` passed, because counting legend
     entries says nothing about where the box lands.

Each is now a check that names the render it came from (`F7` clipping, `F8`
small-move visibility, `F9` legend occlusion, `F10` label overrun), and
`F0a`–`F0e` are planted controls for all five boundaries of the
literal-number check. The pattern across all four: **every one of them is a
check that asked a question adjacent to the one it was written to ask.**

The family, in the form its members take:

- **10-14** a repair to one binding is not a repair to the class.
- **10-17** a class can span two operating points, so it must be partitioned.
- **10-21** a replacement that *reads* a rebound name is neither stale nor
  unassigned.
- **10-22** a capture that fixes one name does not fix the names the captured
  object reads.
- **10-23** a refusal that nothing reaches is not a refusal.
- **10-24** **a name outside every rebind set is not thereby harmless — it is
  only unwatched.** The capture surface is where that distinction is
  measurable, and *in no rebind set* was being read as *cannot move the
  result* in a committed thesis table.

## 6. Mutation test

`graphene_capture_surface_mutation.py` injects seven defects one at a time
into a **copy** of the repository and requires the suite to fail. 7/7
detected, control clean. Two worth naming: undoing the 2026-10-21 repair
itself (M1) drops the score to 23/29, and replacing the bitwise comparison
with a 1 × 10⁻³ relative tolerance (M2) is caught as well — so the suite
genuinely depends on exactness and not on a threshold.

## 7. Environment

`scipy` absent again and installed (1.15.3) — **ninth consecutive** session;
the 2026-10-15 note that it is present remains false. The device VM's
`/sessions` filesystem is **still 100 % full** (9.3 G of 9.8 G, other
sessions' directories, unreadable and undeletable from here), so the
2026-10-23 workaround was needed again and worked unchanged: clone into
`/tmp` (`/dev/sda1`, 4.4 G free) and `TMPDIR=/tmp/piptmp pip install
--no-cache-dir --target /tmp/pylibs scipy` with `PYTHONPATH=/tmp/pylibs`.
This is now a standing environment fact rather than an incident.

No live web search was used and none was needed: the subject is this
repository's own composition measured against itself.
