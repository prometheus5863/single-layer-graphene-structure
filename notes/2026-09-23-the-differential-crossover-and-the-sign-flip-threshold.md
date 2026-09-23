# The differential crossover, and the exact threshold at which a response sign flips

2026-09-23. Written **before** `graphene_differential_crossover_model.py`
exists, and committed first, so that the predictions below are scored rather
than rationalised. This is the fourth session in a row to pre-register, and
the third to do it in a note committed ahead of the code.

---

## 1. Where this comes from

2026-09-22 proved something reassuring and something incomplete.

**Reassuring.** Write the two contact offsets as

    dW_A = m - s,   dW_B = m + s,        s = (W_B - W_A)/2,
                                          m = (W_A + W_B)/2 - w_cross.

A *scalar* offset `δ` on the crossover moves `m` and leaves `s` **exactly**
fixed. The field splits as

    E(x) = m [f(L-x) - f(x)] + s [f(x) + f(L-x)],   f(x) = λ⁻¹ (1 + x/λ)⁻²

into an antisymmetric part carrying `m`, which contributes exactly zero to
the net response `N`, and a symmetric part carrying `s`, which sets
`sign(N)`. A scan over `δ ∈ [-3, +3]` eV, 601 points × 21 pairs = 12621
evaluations, found **zero sign changes**. So the *direction* of a
two-terminal photoresponse — the thing a design rule actually asserts — is
immune to the entire `~5.4 eV` tilde and to a common work-function error six
times its size.

**Incomplete.** That statement is about the *common mode* of the crossover
uncertainty, and only about the common mode. Section 6.12 does not predict a
single crossover: it predicts a **different** one per metal,
`w_cross(d_eq)`, because `Δ_c` depends on the metal-graphene separation.
Two contacts therefore do not share a `w_cross` at all, and the difference
between their crossovers is exactly the part of the uncertainty that 09-22's
theorem does **not** cover — because a per-metal offset changes `s`.

The 09-22 log named this the top open item and put a number on it: *"a
crossover difference of 0.02 eV would flip Au/Pd, where a 3 eV scalar offset
cannot."* This session is that item.

---

## 2. The perturbation, stated precisely

Give contact A crossover `w_cross + c - τ/2` and contact B crossover
`w_cross + c + τ/2`. Then

    dW_A = (W_A - w_cross) - c + τ/2,
    dW_B = (W_B - w_cross) - c - τ/2,

so in the `(m, s)` variables

    m → m - c,          s → s - τ/2.

`c` is 2026-09-22's scalar offset, renamed. `τ` is the new variable: the
**differential crossover offset**, the amount by which the two contacts'
crossovers disagree. It is the complement of `c`, and together `(c, τ)` span
the whole two-metal crossover uncertainty. Nothing here is a new model —
it is the 2026-09-18 signed-carrier field entered at the `dW` level, exactly
as `graphene_per_metal_crossover_model.py` already does.

Why this is worth a session rather than a paragraph: `s` is what sets
`sign(N)`, `τ` is the only thing that moves `s`, and `τ` is a quantity
Section 6.12 can actually estimate for three metals.

---

## 3. Derivation committed in advance

**D1.** `N = 0` exactly when `s = 0`, i.e. at

    τ* = 2 s = W_B - W_A                                            (★)

independently of `c`, of `λ`, of `L`, of the applied bias at zero bias, and
of every transport parameter in `_transport`. At `s = 0` the two contacts
have *identical* offsets `dW_A = dW_B = m`, the field is antisymmetric about
mid-channel, and Validation 3 of 2026-09-22 already established that such a
field gives `N = 0` to machine precision at every offset.

(★) is the useful form of the whole question, because it says: **the
differential crossover offset needed to reverse a pair's photoresponse is
exactly that pair's work-function gap, and nothing else enters.** No model
evaluation is required to obtain it, and the scalar uncertainty of
2026-09-22 is *exactly* irrelevant to it.

The margin table follows from `METAL_WORK_FUNCTIONS` alone. The
uncomfortable entries are visible by inspection: Au/Pd is `0.02` eV, Ni/Au
`0.06`, Ni/Pd `0.08` — and **Ni and Pd are chemisorbed**, i.e. exactly the
metals for which Section 6.12 Result 2 showed the anchored exponential
diverges and refused to return a crossover.

---

## 4. Predictions

Marked **[D]** where the reasoning above already constrains the answer (so a
pass is worth little and only a *failure* is informative), and **[P]** where
the answer is genuinely open to me before running. Saying which is which is
the honest version of this practice; 2026-09-22's four predictions were not
labelled this way and two of them were closer to [D] than they looked.

- **P1 [D]** A bracketed bisection on `N(τ)` will locate the root at `(★)`
  to `≤ 1e-12` eV for all 21 pairs, and the located root will not move by
  more than `1e-12` eV as `c` is swept over `[-1, +1]` eV. *The tolerance is
  the risky part, not the location.*

- **P2 [P]** `|N|` is **not** symmetric about the flip: for at least one
  pair, `|N(τ* + u)|` and `|N(τ* - u)|` at `u = 0.05` eV differ by more than
  1%. Reason to expect it: `m ≠ 0` there and `_transport` is nonlinear, so
  there is no symmetry forcing them equal. Reason it might fail: near the
  root `N` may be locally linear in `s` to well within 1%.

- **P3 [P]** Among the three pairs Section 6.12 can actually reach (Cu/Au,
  Cu/Pt, Au/Pt), the worst-case margin ratio `τ*/|τ_model|` over
  `ℓ ∈ [0.3, 1.5]` Å is **smallest for Cu/Au** and lies in `[2, 4]`. A ratio
  in that range would mean the sign survives, but not comfortably.

- **P4 [P]** The asymmetry of P2 is larger for straddling pairs than for
  same-sign pairs.

- **P5 [P]** No pair among the 21 has `τ*` smaller than the *smallest*
  differential offset the model produces anywhere in `ℓ ∈ [0.3, 1.5]` Å.
  (If this fails, some pair's sign is reversed by the model's own most
  conservative estimate.)

---

## 5. What this session may not conclude

- It may **not** conclude that Au/Pd, Ni/Au or Ni/Pd are sign-reversed.
  Section 6.12 cannot produce `τ_model` for Ni or Pd at all, and the
  extrapolated value is a reductio, not an estimate. The honest output for
  those pairs is a **threshold and a statement that it is unreachable**,
  which is a different and weaker claim than a prediction of reversal.
- It may **not** treat `τ_model` from the one-parameter `Δ_c` family as a
  measurement. It is one anchored exponential with a swept decay length and
  no fitted prefactor, and Section 6.12.7 already says so.
- Where it contradicts an earlier section it must annotate that section in
  place and keep its numbers, per the standing rule.

## 6. Audit item carried in

2026-09-22 found an unbracketed bisection that returned 21 plausible fake
roots and opened *"audit the repo's other analysis-layer code for the same
bug class"*. This session's own root-finder is the first place that lesson
has to be applied, and a root-finder whose answer is **known exactly in
advance** by (★) is the best possible test of it: the bracket check and the
exact cross-check are independent detectors of the same fault.
