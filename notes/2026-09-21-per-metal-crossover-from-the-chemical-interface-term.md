# A per-metal p/n crossover from the chemical interface term Δ_c(d)

**Date:** 2026-09-21
**Open item closed:** "A per-metal crossover `w_cross = W_G + D_c(d)`" — created
2026-09-20 and listed **top** of that session's "not yet covered" list.

> **This file is committed in two parts, deliberately.** Sections 1–4 —
> including the four numbered predictions in Section 4 — were committed
> *before* the model in `graphene_per_metal_crossover_model.py` was written or
> run. Section 5 was appended afterwards. The commit history is the evidence
> that the predictions were not written to match the answer. 2026-09-20's
> methodological note observed that this repo's claims have twice failed on
> widened samples rather than on physics; pre-registration is the cheap half
> of the fix, and enumerating the whole space (Section 4, P0) is the other.

---

## 1. What this repo currently assumes, and why it is an assumption

Since 2026-09-18 every signed contact-field calculation in this repo has used

```python
W_CROSS_CHEM = 5.4            # eV
signed_offset(W_metal) = W_metal - W_CROSS_CHEM
```

a **single crossover work function shared by all seven metals** in
`METAL_WORK_FUNCTIONS`. The 5.4 eV is real and well-sourced — Giovannetti,
Khomyakov, Brocks, Karpan, van den Brink & Kelly, *Phys. Rev. Lett.* **101**,
026803 (2008) — but it is a *derived* number, not a constant of nature.

Khomyakov, Giovannetti, Rusu, Brocks, van den Brink & Kelly, *Phys. Rev. B*
**79**, 195425 (2009) ([arXiv:0902.1203](https://arxiv.org/abs/0902.1203)) —
the companion paper this repo already leans on for the nonlinear
`dW → dE_F` relation (2026-09-20) — writes the interface potential step as the
sum of a charge-transfer term and a **short-range chemical term**,

```
    ΔV(d) = Δ_tr(d) + Δ_c(d)                                        (their §III)
    Δ_c(d) = e^(−γ d) (a_0 + a_1 d + a_2 d²)                        (their Eq. 4)
```

Charge neutrality requires the *total* offset `dW' = W_M − W_G − Δ_c(d)` to
vanish, so the crossover is

```
    W_0(d) = W_G + Δ_c(d)                                           (*)
```

**a function of the metal–graphene separation `d`, not a constant.** The 5.4 eV
figure is (*) evaluated at the physisorbed equilibrium separation
`d ≈ 3.3 Å`, where `Δ_c ≈ 0.9 eV`. Metals that sit at a different separation
have a different crossover, and therefore a different `dW` — which moves
*every* number downstream of `signed_offset()`.

That is the difference between this open item and 2026-09-20's. The nonlinear
`dE_F(dW')` relation **compressed** existing offsets by a few percent without
moving their zero. A per-metal crossover **moves the zero itself**, separately
for each metal, on the same body of evidence.

## 2. The parameters that are, and are not, available

`Δ_c` is the *only* new ingredient (*) needs. Extracting its four fitted
constants from the PRB was attempted and **failed**:

* `ar5iv.labs.arxiv.org/html/0902.1203` renders Eq. 4 correctly in symbolic
  form, and two separate fetches with differently-worded extraction prompts
  both returned the equation **without** numerical values for `γ`, `a_0`,
  `a_1`, `a_2`. The rendered text passes straight from "we parameterize it"
  to the physical discussion.
* This is the **second** extraction failure on this same paper. 2026-09-20
  recorded two passes over the PDF returning mutually contradictory Table I
  contents. The failure mode is different (missing here, contradictory there)
  but the lesson is the same, and the response is the same: **no fitted
  coefficient from this paper is used anywhere in the model.**

What *is* available, and was cross-confirmed by both of 2026-09-20's passes
and by the ar5iv text today:

| quantity | value | role |
|---|---|---|
| `Δ_c(3.3 Å)` | ≈ 0.9 eV | the anchor |
| `W_G` | 4.5 eV | graphene's own work function |
| `W_0(3.3 Å)` | ≈ 5.4 eV | consistency check on the anchor: 4.5 + 0.9 |
| `d_eq` per metal | 2.05–3.41 Å | already tabulated in `D_EQ` (2026-09-20) |
| `d_0` | 2.4 Å | gap-capacitance cutoff (2026-09-20, unrelated to Δ_c) |

## 3. The model actually used: one unknown, anchored and swept

Rather than guess `(γ, a_0, a_1, a_2)`, `Δ_c` is written as a **one-parameter
family pinned to the anchor**:

```
    Δ_c(d ; ℓ) = 0.9 eV · exp( −(d − 3.3 Å) / ℓ )                    (N1)
    w_cross(metal ; ℓ) = W_G + Δ_c(d_eq(metal) ; ℓ)                  (N2)
```

Three properties make this the honest choice:

1. **It reproduces the repo's current model exactly in a limit.** As
   `ℓ → ∞`, (N1) returns 0.9 eV for every `d`, so (N2) returns 5.4 eV for
   every metal and the whole per-metal apparatus collapses **bitwise** onto
   `W_CROSS_CHEM`. There is an exact reduction test, not a plausible-range
   one — the distinction that caught a real bug on 2026-09-17 and again on
   2026-09-19.
2. **It is exact at the one datum that is trusted.** `Δ_c(3.3 Å) = 0.9 eV`
   by construction, bitwise.
3. **The single unknown has a clear meaning and a bounded range.** `ℓ` is the
   decay length of an interface term produced by wavefunction overlap and
   Pauli repulsion. Such terms fall off on the scale of the orbital tails
   involved: **ℓ ∈ [0.3 Å, 1.5 Å]** brackets that generously at both ends,
   and every conclusion below is reported as a sweep over it.

(N1) is a *simplification* of Eq. 4, not Eq. 4. Dropping the polynomial
prefactor is a real approximation and Section 5 reports where it breaks —
which turns out to be the most interesting result of the session.

## 4. Pre-registered predictions

**P0 (method, not physics).** Following 2026-09-20's lesson — "enumerate the
whole space instead of tabulating two representative cases" — every claim
below is to be checked on **all** metals and **all** ordered pairs that the
model admits, not a chosen subset.

**P1 — no doping-sign changes.** No metal in `METAL_WORK_FUNCTIONS` changes
from n-doping to p-doping or back, for any `ℓ ∈ [0.3, 1.5] Å`.

**P2 — small individual shifts.** For the metals the model admits, `|dW|`
changes by **≤ 10%** relative to the flat 5.4 eV convention, across the whole
`ℓ` range.

**P3 — amplification on same-sign pairs.** For a pair whose two contacts sit
on the *same* side of the crossover, the fractional change in the pair's net
response **exceeds** the fractional change in either contact's own `dW`. This
is the mechanism 2026-09-20 identified — a same-sign pair is a near
cancellation, so its surviving magnitude is set by the *ratio* of the two
offsets, and perturbing a ratio is a larger relative perturbation than
perturbing either term. Straddling pairs should show no such amplification.

**P4 — the chemisorbed metals are out of regime, again.** Extrapolating (N1)
inward to `d_eq ≈ 2.05–2.3 Å` (Ti, Ni, Pd) should give crossover work
functions **above the work function of every elemental metal**, which would
make all three maximally n-doping irrespective of their actual work function.
If so, (N1) must be refused there — and for a reason entirely independent of
2026-09-20's, which refused the same three metals because `d_eq < d_0` put
them outside the *gap-capacitance* term. Two different parts of the same
framework failing on the same three metals, for different reasons, would be a
structural statement about Chapter 6 rather than a numerical one.

Sections 1–4 end here. Everything below was written after the model ran.
