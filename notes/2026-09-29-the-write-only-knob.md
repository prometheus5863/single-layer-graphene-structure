# The write-only knob: when a mutation audit mutates nothing

*2026-09-29. Study note for the first half of the session. Companion to
`graphene_figure_provenance_audit.py` and `figure_provenance_audit_output.txt`.*

## The question that was asked

The 2026-09-28 entry left one item at the top of the repository's list. Three of
three band-structure figures had been found to be typed-in known answers rather
than computed quantities, which is why a 14.40 eV gap at the K point survived
five weeks of daily audits: the figures could not show the bug, because the
figures never touched the code that had it. Chapters 4–6's figures had never
been asked the same question.

The detector was named in that entry and is cheap: **mutate the model the
figure claims to come from, and require the figure to change.**

## The answer

They are computed. Six device figures, mutated against the models they claim to
come from, with every physically-meaningful mutation responding:
`quantum_capacitance` to `v_F` and to temperature; the GFET transfer
characteristics to oxide capacitance, mobility and the puddle floor; the RF
figures of merit to channel length, mobility and finger count; the edge-vs-top
contact figure to the Au Fermi shift, by a factor 519.

Three properties of that result matter more than the result:

- **Four null controls held at exactly zero.** A detector that reports
  "everything moved" certifies nothing, and a harness that accidentally
  re-imports a clean module on every run would report exactly that. Mobility
  and contact resistance were mutated by ten times against the quantum
  capacitance figure and changed it by bitwise zero, because `C_q` is
  electrostatic and cannot contain them.
- **Four closed-form validations passed, three of them bitwise.** `R_bare` is
  exactly linear in EQE, so doubling EQE doubles the curve with a ratio of
  exactly 2.0 over 200 decades of trap time. `G·BW` is exactly constant.
  `C_q` scales as exactly `1/v_F²` and is exactly even in gate voltage. These
  do not ask whether the figure moved; they ask whether it moved *by the amount
  the algebra says*, which is the only question that separates a computed
  figure from one that happens to correlate with the parameter.
- **The band-structure fault is genuinely absent from the device chapters.**
  That is worth saying plainly, because the prior after 09-28 was that it would
  be everywhere.

## What the audit found instead

Three `MUST_CHANGE` checks failed bitwise. They were not the band-structure
fault. They were one fault, and it does not yet appear anywhere in this log:

```python
t_liner_nm_default = 3.0

def cu_resistivity_with_liner(W_nm, t_liner_nm=t_liner_nm_default):
    ...
```

`t_liner_nm_default` is evaluated **once**, when `def` executes. Rebinding the
module attribute afterwards rebinds the *name* and leaves the function's
captured default untouched. The constant at the top of the file is a
**write-only knob**: it reads as a parameter and behaves as a comment.

In most repositories that is a mild API wart. Here it is specifically
dangerous, because **this repository's entire audit methodology is
mutation-based**. A mutation applied the obvious way — `setattr(module, CONST,
value)`, which is what a sensitivity study, a convergence check or a future
provenance audit would all do — is silently a no-op, and the audit that applied
it reports *no sensitivity* when what it has actually measured is *no
mutation*. The two look identical in the output.

Census: **49 occurrences across 14 files**, with a positive control on the
census itself (hand it a snippet containing the pattern and require it to be
reported; hand it a literal default and require it not to be). Fourteen of the
49 are inside audit modules.

## Two cases that are worse than inert

**The figure that misstates its own parameter.**
`plot_resistivity_vs_linewidth` builds its legend from the *live* global and its
curve from the *frozen* default:

```python
rho_cu_liner = cu_resistivity_with_liner(W_nm)          # frozen at 3 nm
ax.plot(W_nm, rho_cu_liner,
        label=f'Cu, {t_liner_nm_default:.0f}nm TaN/Co liner ...')   # live
```

Set the constant to 4 nm and the figure draws the 3 nm curve — 4.8980 µΩ·cm at
W = 20 nm — under a legend reading "Cu, 4nm TaN/Co liner". The 4 nm curve it
names is 7.0000 µΩ·cm. The figure is wrong **in writing, on the figure**, by
42.9% in the plotted quantity, and a reader has no way to detect it from the
figure alone. An inert knob is a knob that does nothing; this is a knob that
changes the caption and not the data.

**The knob that moves the oracle and freezes the measurement.**
`graphene_band_structure_audit.T_HOP` is commented "the hopping used throughout
the repo". It is captured as a default by `bands`, `fermi_velocity_analytic`
and `dos_from_bands` — the functions that **measure** — and it is used live, as
a module global, in that same file's expected-value expressions:
`2.0 * T_HOP * abs(phi(K))`, `0.8 / (4.0 * T_HOP)`, the van Hove comparison.

Doubling it:

| | before | after |
|---|---|---|
| MEASURED `bands(k)` | 8.379013 eV | 8.379013 eV (bitwise) |
| MEASURED `v_F` | 9.0627 × 10⁵ m/s | 9.0627 × 10⁵ m/s (bitwise) |
| ORACLE `2·T_HOP` | 5.600 eV | 11.200 eV |

So the mutation moves the expected value and freezes the measured one. The
audit reports a disagreement **manufactured by its own mutation machinery**, and
an investigator following it would be debugging a Hamiltonian that never
changed. Yesterday's Hamiltonian mutation escaped this only because it replaced
a *function* rather than a constant — luck, not design.

## What was fixed, and what deliberately was not

The nine sites behind the three failures were converted to late-bound defaults.
The hazard in a fix like that is that it quietly moves a published number, so it
was made subject to a bitwise requirement rather than a plausibility check:
every plotted value of all three figures (86 400, 8 112 and 14 592 bytes of
harvested curve data) and all three `summary_numbers()` transcripts were
captured before and after and compared exactly. **None moved.** The fix changes
only what happens when someone turns a knob that previously did nothing.

Section B of the audit no longer narrates the fault; it asserts the two
properties that would break if it came back.

The 37 remaining occurrences, **13 of them inside audit modules**, are left. Not
from lack of time: fixing an audit module from inside the audit that found it
moves the instrument and the measurement in the same step, and this repository
has already recorded (09-28) what happens when an instrument's own correctness
is assumed while it is being used. They are the next run's top item.

## The methodological point, continuing the series

09-25: a procedure asked whether it has converged can answer yes and be 44%
wrong. 09-26: an anchored comparison is anchored in one variable. 09-27: an
instrument can be systematically smallest where the answer is worst. 09-28:
prose is a detector, and a check that cannot fail is worse than no check.

**09-29 adds: a mutation that does not arrive is indistinguishable, in the
output, from a system that does not respond.** "No sensitivity" and "no
stimulus" produce the same number, and every mutation-based instrument in this
repository reports that number the same way. The defence is not care. It is
that a mutation harness must carry a **positive control on the mutation
itself** — a paired assertion that *some* quantity the mutation must reach did
in fact move — before it is entitled to interpret a zero anywhere else. The
null controls in this audit were built to stop the detector over-reporting; the
missing control was the opposite one, and it is the one the three failures
needed.
