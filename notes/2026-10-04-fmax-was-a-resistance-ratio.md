# 2026-10-04 — f_max was a resistance ratio, and the literature says so too

Closes the open item created 2026-10-03: *quantify the g_ds vs. R_g·C_gd
split in the remaining f_max shortfall.* Module:
`graphene_fmax_shortfall_decomposition.py`. Figure:
`fmax_shortfall_decomposition.png`. Transcript:
`fmax_shortfall_decomposition_output.txt` (30/30, 6/6 mutants killed,
registered in the pristine-transcript SUITES the same session).

## 1. The question, and why it was answerable exactly

Chapter 4 §4.6's 2026-10-03 annotation ended:

> …the remaining shortfall sits in g_ds and in the R_g·C_gd feedback term,
> not in the gate resistance. Quantifying that split is left open.

The f_max denominator is a **sum**, so the split needs no perturbation
theory:

```
f_max/f_T = 1/(2√denom),   denom = g_ds(R_g+R_s) + 2π f_T C_gd R_g
                                   \____ A ____/   \____ B ____/
```

Zeroing each term is exact. Every counterfactual below is produced by the
same code path as the baseline with one input deleted — a counterfactual
computed by a second, parallel expression would be testing that
expression.

## 2. The split

40 µm / 8 fingers, L = 200 nm, V_ds = 0.1 V, at the peak-f_T bias:

| term | value | share |
|---|---|---|
| A = g_ds(R_g+R_s) | 6.033531×10⁻¹ | 99.7846 % |
| — g_ds·R_g | 3.175543×10⁻¹ | 52.5182 % |
| — g_ds·R_s | 2.857988×10⁻¹ | 47.2664 % |
| B = 2π f_T C_gd R_g | 1.302534×10⁻³ | **0.2154 %** |

Deleting C_gd entirely moves f_max/f_T from 0.643007 to 0.643701: **+0.108 %**
against a shortfall of a factor of 2.2. B's share rises with drain bias but
only to 1.9908 % at V_ds = 1 V. **Chapter 4 named two causes and one of
them is not a cause.** The sentence is struck through in place, numbers
kept.

## 3. The result that contradicts five weeks of this repository's own work

| counterfactual | f_max/f_T |
|---|---|
| as committed | 0.643007 |
| C_gd = 0 | 0.643701 |
| 16 fingers (R_g/4) | 0.827025 |
| 64 fingers (R_g/64) | 0.927229 |
| **R_g = 0 exactly** | **0.935276** |
| R_s = 0 exactly | 0.885467 |
| g_ds ÷ 4.8483 | 1.410000 |

A **perfect gate** — not a better one — lands 28.06 % below the bottom of
Feijoo *et al.*'s 1.3–1.4 band, and 64 fingers is within 1 % of that
ceiling, so the multi-finger mechanism is spent by about sixteen. The
2026-10-03 session found a real 40× delivery bug in R_g and was right to.
The belief that R_g was the **lever** did not survive measurement. Nothing
in the audit record could have caught this: every check asked whether R_g
was computed correctly; none asked whether R_g mattered.

## 4. Why: the DC model is a resistor

`transfer_characteristic()` computes `Id = Vds / (R_channel(Vg,Vds) +
Rc_total)`. There is no velocity saturation, no pinch-off, no drain-field
cutoff. So the output conductance *is* the channel conductance. Measured:

| V_ds (V) | 0.05 | 0.1 | 0.3 | 0.6 | 1.0 |
|---|---|---|---|---|---|
| g_ds·R_total | 1.0026 | 1.0051 | 1.0153 | 1.0301 | 1.0487 |
| g_m/g_ds | 0.00514 | 0.01025 | 0.03043 | 0.05965 | 0.09648 |
| peak f_T (GHz) | 10.140 | 20.279 | 60.779 | 121.168 | 200.440 |
| f_max/f_T | 0.6442 | 0.6430 | 0.6385 | 0.6312 | 0.6218 |

g_ds ≈ 1/R_total to 0.3–4.9 %, so term A is a pure resistance ratio and

```
f_max/f_T ≈ ½ √( R_total / (R_g + R_s) )      0.645352 vs exact 0.643007
```

to 0.36 %. Reaching 1.41 needs R_total/(R_g+R_s) ≥ 4·1.41² = 7.95; this
device has 26.38/15.83 = 1.67. **f_max is primarily a measurement of output
resistance, and this model has none to measure.** The g_ds that would put
it in the band at its own R_g, R_s, C_gd is 7.859726×10⁻³ S (0.1965 mS/µm)
against 3.810651×10⁻² S (0.9527 mS/µm) — a factor of **4.8483**.

## 5. The literature agrees, and that is worth recording separately

This was derived from the model before it was looked up, and the lookup
corroborates the mechanism rather than the number. From the Chalmers
doctoral thesis on graphene FETs for high-frequency applications
([publications.lib.chalmers.se/records/fulltext/252776](https://publications.lib.chalmers.se/records/fulltext/252776/252776.pdf),
and the summary record at
[research.chalmers.se/publication/513392](https://research.chalmers.se/publication/513392)):

> "f_max is mostly limited by the high drain conductance g_ds due to the
> lack of current saturation."

and on the mechanism:

> "current in GFETs does not saturate via pinch-off, since graphene has no
> bandgap. Rather than pinch-off, the charge carrier type will change
> within the channel from one type to the other when the in-plane electric
> field is strong enough."

with a second, non-intrinsic cause — "charge carriers emitted from
interface states at high fields are preventing the current to saturate
and, hence, restricting f_max." The same source gives the theoretically
achievable saturation velocity as 0.2–0.8 v_F and a measured value around
1.4×10⁷ cm/s (1.4×10⁵ m/s), and records intrinsic f_max = 200 GHz at
L = 60 nm and 105 GHz at L = 100 nm (bilayer on SiC), against intrinsic
f_T = 427 GHz at L = 67 nm. IBM's experiment/simulation/theory study of
current saturation in submicrometre graphene transistors with thin gate
dielectric
([research.ibm.com](https://research.ibm.com/publications/current-saturation-in-submicrometer-graphene-transistors-with-thin-gate-dielectric-experiment-simulation-and-theory))
is the standard reference for the thin-dielectric route to partial
saturation, and the compact large-signal GFET modelling literature
([arXiv:1605.08235](https://arxiv.org/pdf/1605.08235)) is where a
saturation term for `transfer_characteristic()` would come from.

**One honest correction to my own text.** The module's R5 check carries the
comment "real GFETs: order unity and above" for g_m/g_ds. I did not find a
sourced g_m/g_ds figure for a comparable device, so that parenthetical is
an assertion and is marked as one; what *is* sourced is the qualitative
statement above, that f_max in real GFETs is g_ds-limited for exactly the
reason this model exhibits in extreme form. The check's threshold (0.2) is
set from this model's own measured range, not from literature, and the
check is therefore a **drift detector**, not a comparison against the
world. 2026-10-02: agreement between two artefacts is silent about the
world — and so is disagreement with a threshold I chose.

## 6. A third result, about f_T, and a withdrawn inference

Peak f_T scales ×19.8 over a ×20 drain-bias range, because g_m =
d(V_ds/R)/dV_g ∝ V_ds in a resistor model. A saturating device's f_T does
not. So the ≈20 GHz that §4.6 reported as "consistent with the literature
range for non-exotic gate lengths" is a fact about the bias
`plot_fT_fmax()` passes; the same model matches the 200 GHz record devices
at V_ds = 1 V. **No number in §4.6 is retracted. What is retracted is
reading their proximity to the literature as corroboration.** A model with
one knob that moves the answer ×20, compared against a literature range
spanning ×20, will agree somewhere.

Free diagnostic from the same table: f_max/f_T *falls* with drain bias here
(0.6442 → 0.6218) where a real device's rises as it saturates. The sign is
wrong, not just the magnitude, so the first thing a saturation term must
fix is a trend.

## 7. Two faults of my own, both caught by this repository's own guards

**(a) The module written about frozen defaults contained four.** The first
committed form of `decompose()` read

```python
def decompose(Vg=VG_SWEEP, Vds=VDS_REF, W=rf.W_RF,
              N_FINGERS=rf.N_FINGERS_RF, ...)
```

— two of its own module constants and two of
`rf_small_signal_model`'s, captured as function-parameter defaults at `def`
time. For `rf.W_RF` and `rf.N_FINGERS_RF` that is the **exact cross-module
shape of the 2026-10-03 bug**, committed by the session that had just read
the entire notes file about it. It was caught by
`graphene_mutation_arrival_probe.py` §7b on the first run after the module
landed — the first time that census has had anything to report since the 23
sites in 7 modules were converted on 2026-10-01. **Knowing a failure mode
is not the same as not committing it; a standing guard is.** All four are
now late-bound, with the builtin `Ellipsis` as the "unsupplied" sentinel
rather than a module-level `_DEFAULT = object()`, because a module-level
sentinel is itself a captured module-level name and §7b reports it
correctly by its own rule. The rule was left as written and the code moved
instead (2026-10-01: a classifier's name must not drift from what it
measures).

**(b) The guard detected it and then crashed.** §7b's renderer read five
fields out of rows the census has only ever built with four, so **any
non-empty census raised `IndexError` on its first element** and the module
exited 1 before Sections 8 onward. It had never fired, because the census
has returned zero rows every run since 2026-10-01 — **a reporting path that
only runs when the guard fires had never run.** The operator saw a stack
trace where four named sites were the intended output. A crash is a
detection, but a strictly weaker one than a failed check: no diagnosis,
every later section lost, and it rests on an arithmetic accident downstream
rather than on an assertion. Fixed as `_print_census()`, and **C7b now
renders a synthetic three-row census (FROZEN, MIXED, IMPORT_FAILED) on
every run**, so the fired-guard path is covered by the suite rather than by
the day the guard fires. Without C7b the fix would itself be untested until
the next captured default, which on this week's evidence is three days to
never. Suite 33/33 → 34/34, exit 1 → 0.

## 8. Methodological note

09-28: prose is a detector. 09-29: a mutation that does not arrive is
indistinguishable from a system that does not respond. 09-30: a control has
to sit where the failure enters. 10-01: and name what it compares against
in a way that cannot drift. 10-02: when it agrees, that is a fact about two
artefacts, not about the world. 10-03: a PASS/FAIL at zero is silent about
magnitude.

**10-04: AND A QUANTITY CAN BE COMPUTED CORRECTLY UNDER THE WRONG NAME.**
Every check in this repository asked whether f_max was computed correctly
from g_ds, R_g, R_s and C_gd. It was — bitwise, reproducibly, with the
40× delivery bug found and fixed. What none of them asked is **what the
resulting number is a measurement of**, and the answer turned out to be the
device's access-to-channel resistance ratio, which is not f_max. There is
no tolerance that detects this, no mutation that kills it, and no
transcript that goes stale: the arithmetic was right the whole time. The
only thing that found it was asking the denominator which of its terms
mattered — a question about *attribution*, not about *correctness*.

The corollary joins 10-03's. Yesterday: an instrument's agreement with its
own specification is silent about whether the specification covers the
case. Today: **an instrument's correctness is silent about whether it
measures the quantity its name claims.** And the sharper, more
uncomfortable form, because it is about where effort went rather than where
error did: the R_g thread was *correct work aimed at a non-binding
constraint*, and five weeks of audits could not have told anyone, because
auditing asks whether a number is right and never whether it is worth
having.

**Sources:**
[Chalmers doctoral thesis, GFETs for high-frequency applications (full text)](https://publications.lib.chalmers.se/records/fulltext/252776/252776.pdf) ·
[the same, record page](https://research.chalmers.se/publication/513392) ·
[IBM, current saturation in submicrometre graphene transistors with thin gate dielectric](https://research.ibm.com/publications/current-saturation-in-submicrometer-graphene-transistors-with-thin-gate-dielectric-experiment-simulation-and-theory) ·
[arXiv:1605.08235, large-signal compact model of GFETs, Part I](https://arxiv.org/pdf/1605.08235) ·
Feijoo *et al.*, *Sci. Rep.* **6**, 35717 (2016), as cited in
`notes/2026-08-26-fmax-parasitics-and-fT-fmax-ratio.md`.
