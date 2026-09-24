# The bracket-bug audit: one real instance, three misclassified suspects, and a fault the fix itself introduced

**Date:** 2026-09-24
**Repo:** single-layer-graphene-structure
**Open item addressed:** "Audit the repo's remaining analysis-layer code for
the bracket bug class — `alpha_from_separation`, the `ell` bisection of
2026-09-21, `contact_resistance_crossover.py` and
`graphene_sensitivity_audit.py`. Today's new root-finder is guarded; the four
older ones are not known to be." — created 2026-09-23.
**Code:** `bracketed_root.py`, `graphene_rootfinder_audit.py`,
`rootfinder_audit_output.txt`; the fix is in
`graphene_per_metal_crossover_model.p2_threshold`.

---

## 1. What the audit was looking for

On 2026-09-22 a bisection in `graphene_crossover_sensitivity_model.py` ran
with no check that its interval contained a sign change. When none existed it
walked one bound onto the other and reported that bound as a root — twenty-one
times, every one of the twenty-one a physically plausible work function, while
all five of that module's exact validations passed and printed underneath it.
On 2026-09-23 the answer was a `bracketed_bisect` that raises on an
unbracketed interval, and the log left an item asking whether anything else in
the repo had the same hole.

The item named four suspects. This is the audit of those four.

## 2. RESULT — the item's own list was wrong about three of its four entries

Three of the named "root-finders" do not solve anything:

| candidate | what it actually is | bracket bug reachable? |
|---|---|---|
| `p2_threshold` (the "`ell` bisection of 2026-09-21") | 200 unchecked halvings | **yes — the one real instance** |
| `alpha_from_separation` | affine in `d` | no; cannot iterate |
| `crossover_length` | closed form, `L_x = R_c W σ_sheet` | no; one evaluation |
| Family-B calibration "solve" | a single division | no |

This is not a reading of the source. Validation 1 tests it by identity, and
each identity is false for anything iterative:

- `crossover_length(2R_c) == 2 · crossover_length(R_c)` **bitwise** — exact
  linearity in `R_c`, which no bisection to finite tolerance reproduces to the
  last bit;
- `alpha_from_separation` at the midpoint of `[2.8, 3.8] Å` equals the mean of
  its endpoint values to **1 ulp** — affine;
- `λ_imp == λ_bulk / Σ(terms)` **bitwise**.

The two entries that were misclassified do carry hazards, and naming them was
not wasted: `alpha_from_separation` returns a *negative* α for `d < d₀`, and
the Family-B division has a denominator that can vanish. Both are real, both
are a different class, and the first already has its guard in
`alpha_for_metal()`. Recording that distinction is the substance of the
result — **an audit list assembled by recalling which functions "solve for
something" picks up functions that merely compute something, and it took a
mechanical test rather than a re-reading to separate them.**

## 3. RESULT — the one real instance, and it never produced a wrong number

`p2_threshold(tol=0.10, lo=0.05e-10, hi=5.0e-10)` returns the shortest decay
length ℓ at which every usable metal's |ΔW| is still within `tol` of its
flat-crossover value. From 2026-09-21 it ran 200 halvings with no bracket
check at all.

At the shipped defaults the interval **is** bracketed:

```
f(lo) = worst(0.05 Å) − 0.10 = +1.370649
f(hi) = worst(5.00 Å) − 0.10 = −0.090361
root  = 0.499733219 Å,  worst(root) − tol = −9.6e-16
```

So Chapter 6's **0.50 Å is a genuine root** and nothing published moves. The
fault was **latent**: `tol`, `lo` and `hi` are keyword arguments with
defaults, and the first caller to move any of them out of the bracketed region
gets a plausible number and no warning.

Naming an unguarded solver is not the same as showing it can produce a wrong
number, so Validation 3 wrote the two wrong numbers down *before* running:

| interval | why root-free | predicted | measured |
|---|---|---|---|
| `tol = 2.0` (above `worst(lo) = 1.4706`) | `f < 0` everywhere | collapses toward `lo` | **1 ulp above `lo`**, not `lo` |
| `tol = 0.005` (below `worst(hi) = 0.009638`) | `f > 0` everywhere | `hi` exactly | `hi` **exactly** |

Both are perfectly reasonable-looking decay lengths. Neither is a root:
`worst(0.0500 Å) = 1.4706` against a target of 2.0, and
`worst(5.00 Å) = 0.009639` against a target of 0.005.

### 3.1 A correction to the 2026-09-22 wording

The 2026-09-22 note recorded that an unguarded bisection "returns the
**ENDPOINT**". That is true in one direction and false in the other, and the
asymmetry is a rounding fact rather than a detail:

- when the loop takes `lo = mid`, it converges onto `hi` and lands there exactly;
- when the loop takes `hi = mid`, it **stalls one ulp above `lo`**, because
  once `hi` is the next double after `lo`, `0.5·(lo + hi)` rounds back up to
  `hi`.

The practical consequence is what makes this worth a correction rather than a
footnote: `result == lo or result == hi` looks like a cheap test for this
fault, and **it is false half the time while the result is meaningless
anyway.** The earlier wording would have reassured exactly the caller who
went looking. The guard has to be at the entrance; there is no reliable check
at the exit.

## 4. RESULT — the fix introduced a second fault, of a class the guard cannot see

The first version of the fix failed its own Validation 5, and the cause was
not the guard but the guard's **default tolerance**.

`bracketed_bisect(..., tol=1e-14)` compares an **absolute** interval width.
That default was chosen in `graphene_differential_crossover_model` for τ, a
variable of order 1 eV, where 1e-14 is ~machine precision: 48 halvings,
relative precision 1e-14. `p2_threshold` bisects a decay length in **metres**,
of order 5e-11. The same 1e-14 is then *coarse*:

| scale | halvings before `(hi−lo) < 1e-14` | relative precision |
|---|---|---|
| O(1) eV — τ, the original caller | 48 | 1e-14 |
| 5e-11 m — ℓ, the new caller | **16** | **2e-4** |

Measured:

```
tol=1e-14 (absolute default) : 4.99749374389648371e-11 m   residual −3.364e-06
rtol=eps  (scale-free)       : 4.99733219460097935e-11 m   residual +2.359e-16
```

2.5e11 ulps apart. Wrong in the fifth significant figure, and **nothing
raised**. The call now passes `tol=0.0, rtol=eps`; `rtol` defaults to 0.0 in
`bracketed_root.py` so every pre-existing caller is bitwise unaffected.

**The lesson, which is why this got its own validation rather than a quiet
patch: a guard moved to a new call site carries its defaults with it, and a
default tolerance is a claim about the scale of the caller's variable.**
Reviewing the guard for correctness would never have found this — the guard
was correct. It was caught only because the pre-existing unguarded loop had
been returning the right answer for three days and could be used as an oracle.
Deleting the old body before validating against it would have destroyed the
only evidence that the new one was worse.

That is the second time in eight days that a *correct* piece of machinery
produced a wrong number at a new call site (the first: 2026-09-20's
unrepresentative sample), and the first time the wrongness was invisible to
every check except a comparison against a known-good value.

## 5. RESULT — the monotonicity the bisection rests on is provable, not merely observed

`p2_threshold`'s docstring asserted from 2026-09-21 that the |ΔW| shift is
monotone decreasing in ℓ. That assertion is load-bearing: a bisection on a
non-monotone function can bracket a root and return the wrong one of several,
so a checked bracket buys nothing without it. It was asserted, not shown.

It is provable term by term. For metal *m*,

```
worst_m(ℓ) = Δ_c(anchor) · |exp(−(d_m − d_anchor)/ℓ) − 1| / |f_m|
```

and `d_m − d_anchor` has a fixed sign per metal, so the exponential approaches
1 monotonically from one side as ℓ grows. Each absolute value is therefore
monotone non-increasing, and a pointwise max of monotone non-increasing
functions is monotone non-increasing. No grid required.

The sharper test is an exact consequence rather than the proof: **Pt's
`d_eq` is 3.30 Å, equal to `D_ANCHOR` exactly**, so Pt's term is `0.0`
bitwise at every ℓ and the max is really over Cu and Au alone. Any code
disagreeing with that is not evaluating the formula above. Measured: 801 grid
points, Pt exactly zero at all of them, Cu (`d − d_anchor = −0.04 Å`) and Au
(`+0.01 Å`) each non-increasing, aggregate non-increasing.

Continuing the 2026-09-23 distinction — a claim derived from a symmetry and a
claim checked on a wide sample are not the same kind of knowledge — this
converts a sampled claim into a derived one. The repo's second such
conversion.

## 6. Housekeeping done, and one thing deliberately not done

`bracketed_bisect` now lives in `bracketed_root.py` and
`graphene_differential_crossover_model` imports it. Two copies of a guard is
a guard that can drift, and the audit's finding that a *second* module needed
it is what forced the move.

Verified rather than assumed: all six validations and all four results of that
module were run before and after the move and their printed output is
**bitwise identical, 107 lines**. (The two `FAIL` lines in that output are
2026-09-23's predictions P2 and P4, falsified then and unchanged now.)

Not done: `__pycache__` is tracked in this repo, so every run of any script
dirties the working tree. It is repo hygiene rather than physics and it has
no bearing on any result, so it is logged as a candidate rather than folded
into an audit commit.

## 7. What this does and does not settle

Settled: the bracket-bug class is now closed across the analysis layer, with
the one real instance guarded and the three non-instances shown to be
non-instances by identity.

Not settled, and worse than before: the class of faults that *the guards
themselves* introduce now has one member with a measured cost, and there is
no systematic check for it. Every numeric default in this repo — tolerances,
step sizes, grid counts, `rel_step` in `log_sensitivity` — is a claim about
the scale of a variable, and those claims were made at the call site where
each function was born. `graphene_sensitivity_audit.log_sensitivity` uses
`rel_step=1e-5`, which is relative and therefore safe; the grid counts in the
parity study were checked by a convergence study on 2026-09-23. Nothing else
has been looked at.

The honest summary of the day is that an audit asked for one fault class,
found it in one place out of four, and left with a new fault class it had
created itself while fixing the old one.
