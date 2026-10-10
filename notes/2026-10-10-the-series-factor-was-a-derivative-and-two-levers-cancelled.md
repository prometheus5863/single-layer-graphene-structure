# 2026-10-10 — The series factor was a derivative, and the fifth lever was two levers that cancelled

Session notes for Chapter 7 §7.9 **item 12**, created 2026-10-08, blocked by
item 9, unblocked 2026-10-09, closed today. Module
`graphene_exact_charge_rederivation.py`, 21/21 checks, transcript
`exact_charge_rederivation_output.txt`, figure
`exact_charge_rederivation.png`. Thesis: §4.6.8, §7.8.1g, §7.9 items 12/15/16.

---

## 1. The item and what it actually was

The item read: *re-derive §§4.3–4.6 on (4.28) and re-report, keeping the (4.29)
numbers annotated in place.* (4.28) is the exact charge relation under
2026-10-08's declaration that `V_ch` is the quasi-Fermi potential,

```
V_g − V_dirac − V_ch = e·n/C_ox + E_F(n)/e,   E_F(n)/e = A_F·√n      (4.28)
n = (C_ox·dV/e)·C_q(dV)/(C_q(dV) + C_ox)                             (4.29)
```

and the chapter has described (4.29) as *the linearisation of* (4.28) in three
places since 2026-10-08. **That description is checkable in one line and it is
wrong.** A linearisation shares its function's leading order by construction.
These two do not share one:

| `dV` (V) | d ln n/d ln dV, (4.28) | d ln n/d ln dV, (4.29) |
|---|---|---|
| 10⁻¹⁰ | 1.99999854 | 1.00000000 |
| 10⁻⁴ | 1.943766 | 1.000000 |
| 10⁻² | 1.274493 | 1.002201 |

(4.28) has `a·n + b·√n = dV`; as `dV → 0` the `√n` term wins and `n → (dV/b)²`.
(4.29) has `C_q(0)` finite and so is linear. Quadratic against linear.

### What (4.29) is instead

Differentiate (4.28):

```
d(dV)/dn = e/C_ox + e/C_Q(n),   C_Q(n) := e²·dn/dE_F = 2e√n/A_F      (4.30)
dn/d(dV) = (C_ox/e)·C_Q(n)/(C_Q(n) + C_ox)                           (4.31)
```

(4.31) **is** (4.29) with `dn/d(dV)` where (4.29) writes `n/dV`. The series
factor is not an approximation *of* the charge relation; it is the charge
relation's **exact derivative**, and (4.29) applies it over the whole drive as
if it were algebraic — a **one-point rectangle rule** for its own integral.

That reframes the whole comparison. A truncated expansion gets better as the
argument shrinks; a rectangle rule gets worse where the integrand's slope
varies fastest, which here is near the Dirac point. It is why 2026-10-08
measured +68 % there and was right to decline to quote it.

Both identities are checked rather than asserted: `e·dn/d(dV)` from (4.28)'s
closed-form root equals the series combination of `C_ox` with `C_Q(n)` to
**1 ulp**, and integrating (4.31) numerically back to the closed-form root
converges at observed order **3.994** with the `h⁴` extrapolation landing
**1.5×10⁻¹³** relative from it.

---

## 2. Two errors, separated for the first time

(4.29) deviates from (4.28) in two unrelated ways:

1. **quadrature** — the rectangle rule above;
2. **argument** — `quantum_capacitance()` evaluates `C_q` at the gate
   overdrive `dV`, while (4.30) requires `E_F/e = A_F√n`, which is 4.04 % of
   the overdrive at the RF bias. Noticed in §4.6.5's 2026-10-07 annotation and
   never costed until today.

The 2×2 ({rectangle, integral} × {`dV`, `E_F/e`}), raw charge, floor off:

| `dV` (V) | total | quadrature alone | **argument alone** | sum of the two |
|---|---|---|---|---|
| 0.01 | +68.2193 % | +32.1722 % | **+68.0927 %** | +100.2649 % |
| 1.20 (RF) | +5.2064 % | +2.6396 % | **+4.7416 %** | +7.3812 % |
| 2.70 | +3.4729 % | +1.7523 % | **+3.2193 %** | +4.9716 % |

Three results. The **argument** error is the larger single contributor, which
nobody could have known while the two were reported as one number. Neither is
the explanation alone. And they are **not additive** — the sum overstates by
1.42× at the RF bias and 1.47× at `dV` = 0.01 V, because they compose
multiplicatively in the series factor. A decomposition reported as a sum would
have erred in the flattering direction.

The puddle floor then attenuates all of it by **4.3×** (+5.21 % raw → +1.22 %
in what `carrier_density()` returns). Worth saying plainly: that attenuation is
a property of `n_puddle`, which 2026-10-09 measured to be 2.0678× off Chapter
3's floor. The re-derivation moving little is a fact about the floor, not
reassurance about the charge model.

---

## 3. The re-derivation, and the item was mis-scoped

| quantity | (4.29) | (4.28) | change |
|---|---|---|---|
| `n` at 3.5 V (m⁻²) | 8.170466e15 | 8.000036e15 | −2.0859 % |
| `R_ch` at the Dirac point (Ω) | 624.1509 | 624.1509 | 0.0000 % |
| `L_x`, Pd / 300 / 500 (nm) | 115.20 / 314.17 / 523.62 | 112.79 / 307.62 / 512.70 | −2.0859 % |
| terminal on/off | 1.240756 | 1.230617 | −0.8171 % |
| intrinsic on/off | 1.614570 | 1.581261 | −2.0630 % |
| `I_d` saturated, RF (A) | 8.190217e-05 | 8.140885e-05 | **−0.6023 %** |

The last line reproduces §4.6.6's published figure to four decimal places from
an independently written patch path. By 2026-10-02's rule that is a fact about
two artefacts, not a validation — what it buys is the assurance that the
(4.29) column is the committed model.

0.0000 % at the Dirac point because the floor supplies the entire density
there. On/off again moves less at the terminals than intrinsically: the
**third** parameter in three sessions to do so, so it is now recorded as a
property of the device topology rather than of any one parameter.

### The mis-scoping

Item 12 was written as a charge item. Under (4.31) the small-signal gate
capacitance **is** the charge relation's derivative, so re-deriving the charge
forces `C_gs` to move as well — and `f_T = g_m/(2πC_gs)` sees both:

| | (4.29) | (4.28) charge only | (4.28) + `C_Q` |
|---|---|---|---|
| peak-`f_T` bias (V) | −1.1115 | −1.1992 | −1.1867 |
| `f_T` (GHz) | 20.2787 | 19.8734 | 20.2666 |
| change | — | **−1.9986 %** | **−0.0594 %** |

The combined −0.06 % is a **cancellation of two ~2 % terms**, not a small
effect. Reporting the charge alone overstates by 34×; reporting only the total
calls the item negligible. Every check in this repository asks about one
quantity at a time, and none of them would see this.

And the **location** of the peak moves 6.76 % against 0.0594 % in height, a
ratio of 114. Both peaks stay on the hole branch, which answers §4.6.7's
peak-`f_T` branch-crossing question NO for this pair.

---

## 4. Five checks failed as first written. None was a code defect.

Continuing 2026-10-09's series, which said: *an exactness check that fails is a
measurement of whatever it was written on top of.*

**(a) The quadratic-exponent check**, swept only to `dV` = 10⁻⁴ V, returned
1.9438 against ±0.01. The two terms of (4.28) are equal at
`dV_x = 2·A_F²·C_ox/e` = **6.519 mV**, so 10⁻⁴ V is 65× below the crossover,
where the correction is O(√(dV/dV_x)) ≈ 12 %. It measured the **width of the
crossover**.

**(b) The integral identity** failed at 1.1×10⁻⁷ with order 3.49 on a uniform
grid — and it measured the **same 6.519 mV**. The integrand of (4.31) turns
over inside a layer of that width at the lower limit, and 4001 uniform nodes
over [0, 1.2 V] put 22 of them in it. A grid graded about `dV_x` reaches order
3.994 with 12× fewer points.

Two exactness checks, in different sections, written for different purposes,
both failing by measuring the same number. That is the strongest evidence this
repository has that 10-09's rule is a rule.

**(c) "Every number moves by under 2 %"** failed at 2.0859 %. The 2 % was a
round number chosen before the run, so it was **deleted** rather than widened —
widening it to 2.5 % is exactly the tuned tolerance §7.9 item 15 exists to
find. The magnitudes are now printed with no threshold and the check tests the
qualitative claims. Separately, a guessed 10⁻¹⁰ bound on (b) was replaced by
the **observed convergence order and the absence of a plateau**, a criterion
the method supplies rather than the author.

**(d) and (e) were wrong in their prose and failed on it.** One asserted the
capacitance half was the *larger* half (0.970×, so no). One asserted `f_T` was
"two orders below the 100–300 GHz RF requirement" — about one order, and
§7.8.1a does not rest on a frequency threshold at all. 2026-09-28: prose is a
detector. Today it detected the author's own checks.

**Today's addition to the series.** 10-09 listed three repairs available at the
moment an exactness check fails — loosen the tolerance, fix the model, find out
what it measured. There is a fourth, and this session needed it: **delete the
criterion and check the claim instead.** (c) was not a measurement that wanted
a better threshold; it was a threshold with no claim underneath it.

---

## 5. A delivery fault found on the way

`contact_resistance_crossover.py` opens with `from graphene_fet_model import
carrier_density`, so it holds its own reference, bound at import. Rebinding
`graphene_fet_model.carrier_density` does not reach it, and a re-derivation
that patched only the owning module would have left §4.4 on (4.29) **while
reporting it as re-derived**. This is 2026-10-03's frozen-default fault in its
monkey-patch form. An arrival control (2026-09-30's rule — a control has to sit
where the failure enters) now precedes the §4.4 numbers and fails if the patch
did not arrive where they are computed.

---

## 6. Literature

Live web search was available: **one search and one fetch, both successful.**
No PubMed/PMC attempted, per the standing note (three reCAPTCHA occurrences
logged).

[Compact modeling technology for the simulation of integrated circuits based on
graphene field-effect transistors, **arXiv:2209.00388**](https://arxiv.org/pdf/2209.00388)
— the Pasadas/Jiménez-lineage compact-model review — independently supports
**both** of today's diagnoses, and neither as a novelty:

- its `C_q` is defined as `∂Q_net/∂V_c` and **evaluated at the local chemical-
  potential shift `V_c`**, never at the gate overdrive. So this repository's
  argument error is a deviation from the literature's own convention rather
  than a defensible modelling choice;
- its series-capacitance expression appears as
  `dV/dV_c = 1 + C_q(V_c)/(C_t + C_b)` — explicitly a **differential**
  relation, used for charge partitioning and *not* to solve the static
  gate-to-channel relation, which it treats as implicit and hands to the
  circuit simulator to solve numerically;
- it states the Fermi–Dirac form "is not convenient for a compact model" and
  approximates `C_q ≈ k·c₁·√(1 + (V_c/c₁)²)`, which integrates to an `asinh`
  form. **The `c₁ → 0` limit of that is exactly (4.30)**, and `c₁` is the
  thermal rounding (4.28) drops.

That last point is the honest placement of (4.28): it is the T = 0 limit of a
standard literature approximation, "exact" with respect to the declaration and
the Dirac dispersion and not with respect to 300 K physics. The regime where
`c₁` matters is `E_F ≲ kT`, far inside the puddle floor, so it is stated as a
limitation and not costed.

Also notable: the review's closed-form-root status differs from this thesis's.
It does **not** have one, because it keeps thermal broadening; (4.28) has one
*because* it drops it. The convenience and the limitation are the same choice.

---

## 7. What this leaves open

- **Item 16, new and already costed.** Fix `quantum_capacitance()`'s argument.
  +4.7416 % of the raw charge at the RF bias — larger than the quadrature
  error, larger than anything else open on Chapter 4, and **not latent**: it is
  active in every committed number in Chapters 4, 5 and 6. Ranked low only
  because that module owns committed transcripts and the change needs its own
  output-neutrality pass. It would be the first change in this thesis to
  **move** committed numbers rather than annotate them.
- **A check that can see a cancellation.** §7.8.1g's new finding is that a
  small number can be small because two larger numbers oppose, and that no
  check in this repository would notice, because they all ask about one
  quantity at a time. That is a gap in the method, not in a model.
- **The `d|n|/dV` = 0 half of item 11**, untouched again.
- **Item 15**, five sessions overdue, and now with two worked repairs attached.
- **A second anchor for `Δ_c`** — still the top *physics* item, open since
  2026-09-21, **untouched for twenty consecutive sessions**. §7.7's unequal
  scrutiny gained a fifth instance today and this remains the standing one.
