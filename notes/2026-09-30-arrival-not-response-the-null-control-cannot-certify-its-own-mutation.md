# Arrival is not response: why a null control cannot certify its own mutation

**2026-09-30.** Study note accompanying `graphene_mutation_arrival_probe.py`.
Closes the top open item of 2026-09-29 (the 37 remaining write-only knobs, 13
of them in audit modules) and, in the process, corrects the methodological note
that item was created by.

---

## 1. What 2026-09-29 concluded, and the part of it that is wrong

The 09-29 session found that a module-level constant captured as a default
argument is a **write-only knob**: `def` evaluates defaults once, at definition
time, so `setattr(module, 'CONST', v)` rebinds the *name* and leaves the
captured object alone. In a repository that audits by mutation this is worse
than a latent bug, because the audit then prints "no sensitivity" when what it
measured was "no mutation". That finding stands, and today's work is built on
it.

The defence it proposed does not stand. In the 09-29 words:

> a mutation harness must carry a **positive control on the mutation itself**
> -- a paired assertion that some quantity the mutation must reach did in fact
> move -- before it is entitled to interpret a zero anywhere else.

This is unsound in two separate ways.

**(a) It cannot be applied to a null control at all.** A `MUST_NOT_CHANGE`
check asserts that the output does *not* move. There is, by construction, no
quantity in that check in which to require movement. So the rule exempts
exactly the checks whose entire evidential value is a zero -- the ones most
exposed to a zero that means something else.

**(b) Even for a `MUST_CHANGE` check it confounds two failures.** If the
downstream number does not move, "the mutation did not arrive" and "the model
is genuinely insensitive here" produce the same zero. That *is* the 09-29
finding, restated as its own remedy.

## 2. The theorem, proved rather than argued

`null_controls_carry_no_arrival_evidence()` in the probe runs four cases
against one synthetic model with a `LIVE` constant and a `FROZEN` one:

| check | constant | outcome |
|---|---|---|
| `MUST_CHANGE` | LIVE | passes -- detects the mutation |
| `MUST_CHANGE` | FROZEN | fails -- correctly; this is the 09-29 detector |
| `MUST_NOT_CHANGE` on an independent output | LIVE | passes, and is informative |
| `MUST_NOT_CHANGE` on a **fully dependent** output | FROZEN | **passes** |

The fourth row is the result. The output in that row is 100% sensitive to the
constant being mutated; the check still passes, because the mutation never
arrived. Its zero is *bitwise* the zero of a correct null control. There is no
post-hoc examination of the number that can tell the two apart.

**Consequence for the 09-29 entry.** It reported that "**four null controls
held at exactly zero** ... so the detector is not reporting that everything
moves", and offered that as evidence of the instrument's specificity. Those
four zeros do not carry that evidence. Whatever confidence the 09-29 audit
earned in its own specificity came entirely from the `MUST_CHANGE` checks that
passed in the same run — the null controls were riding on them.

## 3. Where the control actually belongs

Not downstream of the mutation. **At the point of delivery.**

Before any number is interpreted, ask whether the binding the callee will
actually read has changed. That question is decidable by introspection, without
evaluating the model at all -- which is precisely why it works where the
downstream form fails: it is available for null controls, and it is available
for quantities the model is genuinely insensitive to.

The probe answers it with six classes:

| class | meaning | does a rebind arrive? |
|---|---|---|
| `LIVE` | every reader loads it from globals at call time | yes |
| `FROZEN` | every reader captured it as a default | **no** |
| `MIXED` | some of each | **partly, which is worse** |
| `IMPORT_USED` | read only in the module body, at import | no, but not a fault |
| `DEAD` | read nowhere at all | no |
| `ABSENT` | not an attribute of the module | no -- `setattr` creates it |

`MIXED` is not a middle case. If the live readers are an audit's expected-value
expressions and the frozen readers are its measuring functions, a mutation
**moves the oracle and freezes the measurement**, and the audit reports a
disagreement it manufactured itself.

## 4. Four findings against this repository

### 4.1 A one-character typo is indistinguishable from insensitivity

Taking the first `MUST_CHANGE` check of the 09-29 provenance audit and
misspelling its patch key by one character gives output **bitwise identical**
to no mutation over 2004 harvested points. `setattr` has no notion of a typo:
it creates the misspelled attribute and returns normally. So a renamed or
mistyped mutation target is a write-only knob with no syntactic trace anywhere
in the source.

And a second-order finding on top of it: the probe's verdict on that key
degrades from `ABSENT` to `DEAD` *once the mutation has been applied*, because
`importlib.reload` re-executes the source into the same module dict without
clearing it, so the attribute survives every later reload. **Applying a
misspelled mutation destroys the evidence that it was misspelled.** Arrival has
to be classified *before* the mutation; that ordering is a requirement, not an
implementation detail. (This module's own first run got it wrong.)

### 4.2 `ALPHA_ABS` was a decorative name, and the check policing it was vacuous

`graphene_photodetector_model.ALPHA_ABS` — declared at the top of the module
with a comment attributing it to Chapters 2-3 — was read by **no function and
nowhere in the module body**. Two things follow.

The module docstring says the model "reuses the existing 2.3%
universal-absorption result". It did not. The only route absorption could have
taken is `ALPHA_ABS`, and nothing read it; `EQE_BARE = 0.0015` is an
independent literature literal. This is the 09-28 *prose is a detector* class
for the third time and **inverted**: 09-28 and 09-29 both found code
contradicting its prose, whereas here the prose asserts a dependency the code
does not have.

And the 09-29 null control on it — logged as "a real finding if it holds: the
figure must not double-count absorption" — could not have failed. With the name
unread, that check would have passed identically had the figure double-counted
absorption, counted it once, or not counted it at all. The evidence that it was
vacuous was already in the 09-29 log, one item away, unconnected.

**The fix yields a real number.** `internal_collection_efficiency()` now reads
both constants live and reports `EQE_BARE / pi*alpha` = **6.5430%**. Roughly
fifteen of every sixteen absorbed photons are assumed lost before collection.
That is the quantitative content of the phrase "collection bottleneck" the
existing comment used without a number, and it sits underneath every bare-device
responsivity in the thesis. No shipped number moved:
`(EQE_BARE/ALPHA_ABS)*ALPHA_ABS == EQE_BARE` bitwise, and every pre-existing
line of `summary_numbers()` is unchanged.

### 4.3 Three of the 09-29 null controls are identity re-assignments

`graphene_fet_model.T -> 300.0` when `T` is already `300.0`;
`n_puddle -> 5e15` when it already is; `AU_FERMI_SHIFT_SURFACE_EV -> 0.14` when
it already is. Such a check passes with `setattr` replaced by `pass`, with the
constant frozen, and with the model deleted. It is the 09-28 cannot-fail class
living *inside* the 09-29 null controls, which were themselves offered as the
defence against over-reporting.

### 4.4 `MIXED` is the majority class, not the exception

The 09-29 entry treated the remaining 37 occurrences as one fault with `T_HOP`
singled out as "the worst". At the level of **names**, 15 of 20 were in
`T_HOP`'s class — and `A_CC`, in the same file, was one of them and went
unmentioned. The oracle/measurement split is also present outside the audits,
in `D0_SEPARATION`, `H_DIFF`, `W_CROSS_CHEM`, `ELL_DEFAULT` (4 frozen readers
against 1 live) and `N_POINTS` in three separate photodetector modules.

## 5. The fix, measured

All 14 audit-module sites are converted to late-bound reads. The hazard in such
a fix is quietly moving a published number, so it was made subject to an exact
requirement: all four audits run to completion before and after, and all four
stdout transcripts (8232, 9463, 6388, 8421 bytes) compare **identical**.

`T_HOP` before and after, taken by running the pre-fix file out of git beside
the post-fix one rather than by describing it:

```
before  T_HOP MIXED: oracle 2*T_HOP  5.600 -> 11.200   (moved)
                     bands()[0]     -7.225687 -> -7.225687   (FROZEN)
after   T_HOP LIVE : oracle 2*T_HOP  5.600 -> 11.200   (moved)
                     bands()[0]     -7.225687 -> -14.451374  (moved)
```

Doubling the hopping integral used to double every expected value in the module
and leave the computed band energies bitwise unchanged.

A standing guard (`audit_modules_are_fully_delivered`) now fails if any
`*audit*.py` captures a module constant again, so the fix does not depend on
the next session having read this note.

## 6. The rule this asks the repository to adopt

```python
from graphene_mutation_arrival_probe import require_delivery
ok, rows = require_delivery('graphene_fet_model', ['mu', 'T', 'n_puddle'])
```

Called **before** any mutation is applied and before any resulting number is
interpreted, **null controls included**. Total probe: 29/29, of which 8 are
positive controls on the probe itself.

---

## Methodological note, continuing the series

- **09-25:** a procedure asked whether it has converged can answer yes and be 44% wrong.
- **09-26:** an anchored comparison is anchored in one variable.
- **09-27:** an instrument can be systematically smallest where the answer is worst.
- **09-28:** prose is a detector, and a check that cannot fail is worse than no check.
- **09-29:** a mutation that does not arrive is indistinguishable, in the output, from a system that does not respond.
- **09-30:** **a control has to sit where the failure enters, not where it shows.**
  09-29 correctly identified the failure and then placed its control downstream,
  in the output, where the two candidate causes are already indistinguishable —
  which is the very property it had just established. A defence built at the
  point of measurement can only ever compare zeros. Arrival is a property of the
  *binding*, so it has to be established at the binding, before the run; and
  because that check needs no model evaluation, it costs nothing to place there.
  The corollary for this repository is uncomfortable and worth stating plainly:
  **the check most likely to be vacuous is the one whose passing you find
  reassuring**, because a check that passes attracts no scrutiny, and a vacuous
  check always passes.
