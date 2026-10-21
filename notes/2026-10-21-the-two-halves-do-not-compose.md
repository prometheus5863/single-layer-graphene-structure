# 2026-10-21 — The two halves of item 16 do not compose, and the third failure mode no mechanism here covers

**Session focus:** item 16, "flip the default to form `B_T`", ranked first
since 2026-10-12 and untouched for two sessions. The 2026-10-20 log recorded
that "the flip and nothing else is left". This session set out to land it and
found that sentence false in one specific way, and the way is worth more than
the flip.

---

## 1. The result

The two halves of item 16 — `qca.on_consistent_argument()` in the channel and
`blast.contact_on_consistent_argument()` at the contact — are each committed,
each validated and each correct **alone**. Installing both at once does not
give form `B_T` everywhere.

| Metal | committed `B_T` vs form A | **both halves, unrepaired, vs `B_T`** | ratio |
|---|---|---|---|
| Ti | −6.2578 % | **−52.7772 %** | 8.43× |
| Cu | −4.7696 % | **−51.3850 %** | 10.77× |
| Ni | −12.4233 % | **−60.3949 %** | 4.86× |
| Au | −12.4878 % | **−60.7340 %** | 4.86× |
| Pd | −12.4964 % | **−60.8302 %** | 4.87× |
| Pt | −11.8738 % | **−61.9143 %** | 5.21× |

The unrepaired composition overshoots the effect item 16 exists to install by
**4.86× to 10.77×, in the same direction**. It is not a sign error, it is not
a `nan`, and it is not out of range for a quantity measured in units of
10¹⁶ m⁻². A session that landed the flip and regenerated transcripts would
have committed it.

## 2. The cause: one attribute read

`blast.cq_shipped_of_n` builds the contact's form-`B_T` solve on the
**shipped** `gfet.quantum_capacitance`, and that choice is *right* — its
docstring gives two measured reasons (the alternative carries its own
constants and differs by 9.48e6 ulp; the committed `log(2(1+cosh η))`
spelling overflows where the fixed point's first iterate lands). It reads the
name by attribute lookup at call time:

```python
return gfet.quantum_capacitance(qca.A_FERMI * np.sqrt(nn), T=T)
```

and the channel half rebinds exactly that name.

## 3. The general form — and why neither existing guard sees it

This repository has built two mechanisms against delivery faults:

- **2026-10-14, `bc.rebind()`** — the set a patch must reach is the name's
  runtime **identity class**, computed at entry. Guarantees: *no binding of a
  rebound name is left stale.*
- **2026-10-17, `bc.rebind_by_operating_point()`** — a class can span two
  operating points, so it must be **partitioned**, and an unassigned member
  raises. Guarantees: *no class member is silently left at the shipped
  object.*

Both are statements about **which names move**. Neither says anything about a
**replacement that reads one of the moved names**.

> A replacement built on a name the same flip rebinds is **not stale** and is
> **not unassigned**. It is correctly installed, computed by the right
> mechanism, and reported as installed by that mechanism's own report — and
> it reads a different object than the one it was validated against.

That is a third failure mode, and it is invisible to both guards. The audit
that would have caught today's defect without running anything is one
question per replacement: *which module-level names does this function read at
call time, and does the flip rebind any of them?* It is now the first item on
the standing list.

## 4. Two mechanisms, and they partially cancel

| Metal | (i) inner `C_q` alone | (i)+(ii) | (ii) contributes |
|---|---|---|---|
| Ti | −53.8545 % | −52.7772 % | **+1.0772** |
| Cu | −52.4054 % | −51.3850 % | **+1.0204** |
| Ni | −62.1287 % | −60.3949 % | **+1.7338** |
| Au | −62.5403 % | −60.7340 % | **+1.8063** |
| Pd | −62.6595 % | −60.8302 % | **+1.8293** |
| Pt | −64.2107 % | −61.9143 % | **+2.2965** |

Mechanism (ii) is that the inner channel solve reaches `qca.series_factor`,
which reads `gfet.C_ox` at call time, while the contact solve has `C_ox`
substituted by `C_interface` for its duration. **(ii) offsets (i) for every
metal.** (ii) is the more obviously wrong-looking of the two — a module global
mutated across a nested call is visible on inspection — and repairing it alone
would have made the error *larger*. This is worth stating generally: when a
composite defect has two parts of opposite sign, fixing the conspicuous part
first is a *regression*, and the only thing that distinguishes that from
progress is having measured both parts separately.

## 5. The repair

Bind the contact replacement to the shipped object **captured at manager
entry** rather than to the live name. `graphene_default_flip.on_default_flip()`
then installs the whole flip through the two managers the 2026-10-20 note
named — `bc.rebind()` for the charge, `bc.rebind_by_operating_point()` for the
capacitance. **This is the first caller to give that manager a non-`None`
contact assignment**; until 2026-10-17 the contact had no replacement of its
own, and `bc.CONTACT_OPERATING_POINT` left it alone with a stated reason,
which was correct exactly as long as that was true.

**18 validations, 0 failures, all against exactly known values:** the identity
assignment is bitwise neutral on all nineteen committed numbers; the flip's
contact half reproduces the committed contact-only `B_T` **bitwise** on all
seven metals; its channel half reproduces `on_consistent_argument()`
**bitwise** on all eight census entries; `n_contact(ΔW_F = 0)` is exactly
`0.0`; the antisymmetry is bitwise; an incomplete assignment raises and
rebinds nothing; an assignment with no stated reason is refused; interleaving
channel and contact calls under one flip leaves both bitwise (the stated
reentrancy hazard, measured rather than assumed); and **V8 is a positive
control on the repair itself** — it asserts the unrepaired flip *does* differ,
so if the two ever agree the module's premise is known to have expired rather
than silently satisfied. Transcript byte-identical across two consecutive
runs.

**Every cost item 16 has ever quoted is reproduced exactly** by the composed
flip — `n` −1.0212 %, `I_d` −0.2530 %, `n_contact` −4.7696 % to −12.4964 %,
§4.7 `R_extra` −1.2583 % to +47.6925 %, Pd's residual +50.9489 → −203.2764
Ω·µm, additive decomposition 4-of-4 negative. That reproduction is the
evidence the defect lives in the composition and in neither half.

## 6. Three defects in this session's own work

- **An unpaired ratio.** The first draft of §1's summary divided the
  smallest-magnitude unrepaired deviation by the **largest**-magnitude
  committed effect and the other way round, quoting "4.1× to 13.0×" from two
  pairs of numbers belonging to different metals. The honest per-metal answer
  is 4.86× to 10.77×. Same shape as the dominant failure mode this repository
  has logged for weeks: *a claim answering a weaker question than the one it
  was written to ask.* Caught by **reading the transcript**, not by a check.
- **Two literal `%%` leaks** in `print` calls carrying no `%` operator —
  cosmetic, and caught the same way.
- **A FALSE FIGURE TITLE, and this is the one worth keeping.** Panel (a) read
  *"overshoots the effect it installs by 4-5x"* — the figure of the
  superseded unpaired calculation. The prose was repaired and the title was
  not, **because a title is not in the transcript and nothing reads it.** The
  title is now computed from the same `ratios` the table prints, so the two
  cannot disagree again. Beside it, the legend in panel (b) sat on top of Cu's
  two tallest bars and its `R_tr` annotation. Neither was reachable from a
  transcript reading `18 passed, 0 failed`.

**That is the second consecutive session in which rendering a figure and
looking at it found defects on the first try, and this time one of them was a
false claim rather than an unreadable one.** The standing item is upgraded
accordingly: the question is not only "is the figure legible" but "does the
figure assert anything the text no longer says".

## 7. Environment

`scipy` absent again and `pip install`ed (1.15.3) — **sixth consecutive
session**; the 2026-10-15 note that it is present remains false. `numpy`
2.2.6 and `matplotlib` 3.10.9 present as always. The bridge was stable
throughout, unlike 2026-10-20; every commit was pushed immediately after
commit regardless, which is now the standing pattern. No live web search was
used and none was needed — today's result is this repository's own code
measured against itself. A figure was copied into the connected folder to be
staged and inspected and then moved to `_to_delete/`, since deletion inside
connected folders is blocked and requesting permission was denied by the
safety classifier on 2026-09-17.
