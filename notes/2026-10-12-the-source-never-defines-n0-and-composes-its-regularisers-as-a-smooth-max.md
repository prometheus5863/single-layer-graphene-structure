# 2026-10-12 — the source never defines `n₀`, and it composes its two regularisers as a smooth max

**Item costed:** Chapter 7 §7.9 **item 18** — "costing the double
regularisation of `C_q` at neutrality". Created 2026-10-11, ranked **first**
for this session, and the first item on that list to have arrived from a
literature source rather than from an internal audit.

**Outcome:** closed. The item was wrong twice over, and a third finding —
reached by reading the source rather than by computing — puts an *upper bound*
on the whole class of error it described.

**Module:** `graphene_double_regularisation.py` (22 checks, 22 passed).
**Output:** `double_regularisation_output.txt`. **Figure:**
`double_regularisation.png`. **Thesis:** §4.6.11, with §7.9 items 17 and 18
annotated in place.

---

## 1. Live search and fetch

**Available and used: one fetch, successful.** No PubMed/PMC attempted, per
the standing note (three logged reCAPTCHA occurrences).

Fetched:
[Southampton bilayer-GFET compact model, `GFetIEEEtran.pdf`](https://eprints.soton.ac.uk/366757/1/GFetIEEEtran.pdf)
— the same document 2026-10-11 used to create item 18. It was read yesterday
for *what it regularises with*; it was read today for *what `n₀` means in it*,
because that was the one residue of yesterday's finding that arithmetic could
not resolve.

## 2. The question the fetch was meant to settle, and did not

Yesterday's §9 compared this thesis's thermal `C_q(0)` against the source's
disorder floor

```
C_q,min = q²√n₀ / (√π · ħ · v_f)          [their Eq. (2), at φ_s1 = 0]
```

and found the ratio **1.1509**. Today's §1 found that this floor is **exactly
half** this thesis's own DOS-exact `C_Q(n₀) = 2e²√n₀/(√π ħ v_F)` — 0.000 ulp
from ½ with one consistent `ħ` — because graphene's
`D(E) = 2|E|/(π ħ² v_F²)` carries the `g = 4` spin–valley degeneracy. So
either their `n₀` means something other than a total carrier density, or the
factor is simply absent.

**The paper does not say.** At `φ_s1 = 0` it states only that "the channel has
been reported to have a charge density, `n₀`", cites two references for it,
and never defines it as total, per-species, or as a disorder parameter. In the
resistance section the same symbol appears as "the minimum charge density".
**No degeneracy factor appears anywhere in the paper, and it gives no density
of states at all** — it goes directly from `n` to `C_q` through Eqs. (1)–(3).

> **The residue stands, and it stands for a reason worth recording: a
> cross-model comparison can be blocked not by disagreement but by the source
> having left a symbol undefined and delegated it to a citation.** The
> correction is not available by converting a convention, because there is no
> stated convention to convert.

Their fitted values are in the `10¹⁶ m⁻²` range — **1.2e16 m⁻² for their Test
A**, which is **2.4×** this thesis's `n_puddle = 5e15 m⁻²`. That is now a
second literature anchor on the standing `n_puddle` question (§4.6.7 had the
first, from the opposite direction), and it points the same way: this thesis's
floor may be low rather than high.

## 3. What the fetch did settle, and it settles it in the measurement's favour

Their **Eq. (3)** composes the floor with the gate-driven part as

```
C_q(φ_s) = (C_q,min² + 2h²) / √(C_q,min² + h²),      h = C_q,var/2
```

This is a **smooth max**, not a sum, and it is verifiable in closed form:

| property | value |
|---|---|
| `C_q(0)` | `9.710792902255e-03` = `C_q,min` **bitwise** |
| `C_q/C_var` as `C_var → ∞` | `1.000000000000` |
| `sup C_q / max(C_q,min, C_var)` | **1.341640786** |
| attained at | `C_var/C_q,min` = **1.000000** |
| closed form | `1.5/√1.25` = **3/√5** = `1.341640786` |

Two consequences, both for the item rather than for the model.

**(1) The literature applies `n₀` twice itself.** Eq. (2) floors `C_q` with
it, and the Drude relation `n = √(n₀² + n*²)` floors the density with it —
with no temperature anywhere in the model. So the comparison item 18 rested on
is **not** between a model that double-counts and one that does not. It is
between two **allocations of regularisers to quantities**: this thesis spends
`kT` on `C_q` and `n₀` on `n`; the source spends `n₀` on both.

**(2) Yesterday's intuition was right about *where* the worst case is and
wrong about *how big* it is.** A smooth-max composition is furthest from pure
max exactly when its two terms are **equal** — which is precisely the
condition yesterday identified as holding at this `n₀` and this `T` (ratio
1.1509, crossover at `n₀* = 0.755 n_puddle` and `T* = 1.151 × 300 K`). That
reasoning was sound. But at that worst case the excess over the larger term is
`3/√5`, i.e. **+34.164 %** — and this thesis sits at **+2.5262 %** because its
two terms are 2.30:1 rather than 1:1.

> **A quadrature-composed regulariser pair cannot exceed the larger of its two
> terms by more than 34.2 %. The "factor 2" was unreachable in principle, not
> merely unmeasured.**

That is an upper bound on the entire class of error item 18 named, derived
from the source's own composition law rather than from this thesis's.

## 4. The two internal findings, in brief

Both are in §4.6.11 in full; recorded here only because they are the reason
the item closes rather than narrows.

**(a) In the committed relation the two regularisers never meet.** `C_q`'s
argument in (4.29) is `dV`, which is never floored, so `n_puddle` cannot reach
`C_q`. And `n_eff = (C_ox·dV/e)·f(C_q)` has `dV` as a **factor** with
`f ∈ (0,1)`, so `n_eff(0) = 0.0` **bitwise** under all four `C_q` conventions
tried, and `carrier_density(V_dirac)` returns **exactly** `n_puddle`. Checked
against exactly-known values — `0.0` and `n_puddle` — not against tolerances.
Applying one regulariser to each of two quantities is a *reuse of one
parameter for two jobs*, a modelling choice to defend; it is not a double
count, which is an arithmetic error.

**(b) Where they do meet, the cost is +2.5262 %.** Only form B_T with the
floor inside `C_q`'s argument puts both on one quantity:
`C_q(0) = 1.991221e-02` against `1.942159e-02` disorder-alone and
`8.437399e-03` thermal-alone — **both/max = 1.025262**, both/sum = 0.714750.
Downstream, **+0.1527 %** in `n` at the RF bias and **0.0000 %** at the Dirac
point, where all four conventions agree bitwise because the quadrature floor
on the *result* overwrites every one of them. `σ_min` is identical to the last
bit in all four.

**The methodological residue.** The "factor 2" was the **ratio of the two
candidate values**, which prices *choosing wrongly between* the conventions
and not *applying both*. That is the 2026-10-06 backlog item — *census every
criterion that is a ratio* — in its **first instance where the ratio was not a
criterion but an estimate of a magnitude**. Five sessions open as a question
about verdicts; it was a question about magnitudes too.

## 5. And the module's first draft failed its own check

Recorded because the ordering is new: the substitutability check was written
**before** the measurement, not during a neutrality pass after it.

Draft 1 built its two switches on `graphene_diffusion_current_model`'s literal
`E_CHARGE`/`HBAR` while `quantum_capacitance()` uses `scipy`'s, and reproduced
the shipped `carrier_density(2.0)` only to **2375 ulp**. **Eighth** instance of
the 2026-10-09 rule; **second** with a physical constant.

The repair is 2026-10-11's repair rather than a better transcription: each
switch **substitutes one shipped name** — `gfet._thermal_rounding`,
`gfet.n_puddle` — and leaves the rest of the call path as shipped code. Form A
then agrees with the shipped function at **0 ulp** by construction. One
arrival counter per patched name, each checked at a call that reads that name.

One caveat kept in the output rather than smoothed over: the `T → 0`
validation reports *exactly* zero departure, which is weaker than `1e-16`
would be, because at `η ~ 10⁷` the term `2·log1p(exp(−η))` underflows to
exactly `0.0` and the two expressions reduce to the same floating-point
operations.

## 6. What this note does not claim

1. It does **not** claim `n_puddle = 5e15 m⁻²` is wrong. The source's 1.2e16
   is a fitted parameter for a *bilayer* device on a different stack.
2. It does **not** claim the factor 2 is resolved. The source leaves `n₀`
   undefined; both 1.1509 and 2.3018 are reported in §4.6.11 and neither is
   installed.
3. It does **not** claim form B_T's floor belongs inside `C_q`'s argument.
   That spelling was built to *measure* the double count, and whether it is
   the right physics is the same unexamined question as what `n_puddle` means.
4. The `3/√5` bound is a property of **their** composition law. This thesis's
   nested spelling is a different function; its +2.5262 % is measured, and the
   bound is quoted as context for the item's "factor 2", not as a bound on this
   thesis's own arithmetic.
