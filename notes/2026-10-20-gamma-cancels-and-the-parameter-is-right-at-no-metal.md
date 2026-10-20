# Gamma cancels, and the parameter is right at no metal at all

*2026-10-20 — costing the 2026-10-18 item "is `ell` even the right parameter
to bound?", and finding the charge both wrong and understated.*

## The item

`graphene_per_metal_crossover_model`'s (N1) is a single exponential

```
Delta_c(d ; ell) = 0.9 eV * exp( -(d - 3.3 A) / ell )          (N1)
```

while Khomyakov et al.'s Eq. (4), which it approximates, is an exponential
times a quadratic:

```
Delta_c(d) = exp(-gamma d) * P(d),   P(d) = a0 + a1 d + a2 d^2  (4)
```

(N1) has **one** decay length; Eq. (4) has a different one at every
separation. 2026-10-18 noticed that `ell` as the 2026-10-15 bound constrains
it is a **chord** length over `[3.3, 5.0] Å`, while `w_cross_for_metal`
evaluates `delta_c` at each metal's own `d_eq`. It created the item, did not
cost it, and named the entry point.

`gamma, a0, a1, a2` are reported nowhere in the paper's text — four extraction
passes, and `arxiv.org` is refused by the egress proxy in both shells, so
Fig. 4 cannot be digitised either. So the obvious objection is that the item
cannot be costed at all without the numbers.

## It can, because gamma cancels

Write `L = ln P`. Then

```
1/ell_chord(d1,d2) = gamma - <L'>_[d1,d2]        (C)
1/ell_local(d)     = gamma - L'(d)               (Lo)
```

and subtracting,

```
1/ell_chord - 1/ell_local(d) = L'(d) - <L'>      (D)
```

**`gamma` cancels exactly.** And since `L' = P'/P` is scale-invariant, `a0`
cancels with it. Two of the four unreported coefficients leave the question
entirely; the gap depends only on `P`'s two **shape** ratios, which can be
swept.

The relative error needs no `gamma` either. `ell_chord` **is** the quantity
(N1) calls `ell` and the repository already sweeps over `[0.3, 1.5] Å`, so

```
ell_local(d) / ell = 1 / (1 - ell * Delta),   Delta = L'(d) - <L'>   (R)
```

is a function of the swept parameter and the shape alone.

## The headline: an exact separation, not a bound

By the mean value theorem there is a separation `xi` inside `[d1, d2]` at
which `1/ell_chord == 1/ell_local(xi)` **exactly**. That is the one separation
where `ell` is the right parameter. Measured over 204 admissible quadratic
shapes:

```
xi in [3.8482, 4.4952] A
largest tabulated d_eq = 3.41 A  (Al)
margin = +0.4382 A
```

**`xi` lies above every tabulated metal, for every admissible shape.** So the
chord parameter is the local decay length at *no* metal in the table — not at
the chemisorbed ones, and not at the four physisorbed ones that sit inside the
chord interval either. Being inside the interval buys them nothing, because
`Delta` vanishes only at the single point `xi`, not throughout.

## The item was wrong in both directions

| the item said | measured |
|---|---|
| `d_eq` range is `2.05–3.31 Å` | `2.05–3.41 Å`; **Al sits at 3.41 Å** |
| metals are "outside the chord" | Al, Ag, Au, Pt are **inside** `[3.3, 5.0]` |
| the charge applies to those below the anchor | the charge is a **universal** over all nine |

So the wording was too generous to itself in one place and far too modest in
another. The correct statement is stronger and simpler than the one it made.

## The cost

`Delta` is `ell`-independent by construction; the ratio is that one number
read at three values of the repository's swept parameter.

| metal | d_eq (Å) | Delta range (1/Å) | ratio @ ell=0.6 Å | bad @ ell=1.5 Å |
|---|---|---|---|---|
| Al | 3.41 | −0.182 … 0.767 | 0.90–1.85 | 3 |
| Ag | 3.33 | −0.198 … 0.792 | 0.89–1.90 | 5 |
| Au | 3.31 | −0.202 … 0.797 | 0.89–1.92 | 5 |
| Pt | 3.30 | −0.204 … 0.800 | 0.89–1.92 | 5 |
| Cu | 3.26 | −0.212 … 0.811 | 0.89–1.95 | 5 |
| Pd | 2.30 | −0.389 … 0.988 | 0.81–2.45 | 14 |
| Ti | 2.10 | −0.425 … 1.011 | 0.80–2.54 | 17 |
| Ni | 2.05 | −0.433 … 1.056 | 0.79–2.73 | 19 |
| Co | 2.05 | −0.433 … 1.056 | 0.79–2.73 | 19 |

"bad" counts shapes for which `ell_local` at that metal is **infinite or
negative** — Eq. (4)'s step is not decaying there, which no single exponential
represents at any `ell`. At `ell = 0.3` and `0.6 Å` there are none; at the top
of the swept range, `ell = 1.5 Å`, there are **92**, and they reach the
physisorbed metals too. The divergence condition is `ell * Delta = 1`, so it
is reached by enlarging either factor.

This **compounds with** the `d0` regime failure rather than replacing it. Ti,
Ni, Co and Pd are outside (N1)'s regime because `d_eq < d0 = 2.40 Å` whatever
`ell` is, **and** are evaluated with a parameter read off an interval they are
not in. Two independent reasons, the same four metals.

## What this does not settle

It prices the mismatch; it does not decide whether (N1) should be replaced.
Deciding needs Fig. 4's point set, which needs the PDF, which needs
`arxiv.org` — closed by policy since 2026-10-18. A copy of the PDF placed in
the connected folder opens it immediately.

## Four defects in this file's own first draft

All four are the same shape and it is the one this repository keeps finding:
**a check or a claim that answered a weaker question than the one it was
written to ask.**

- **Validation 5 claimed `BITWISE` and got `2.13e-14 eV`.** The algebra is
  exact; the *arithmetic* is not the same arithmetic. `delta_c()` forms the
  exponent as `(d - d_anchor)/ell` with both in metres — one division — while
  this file formed it as `gamma * (d - d1)` with `gamma = 1/ell` in Å⁻¹ — a
  reciprocal and a multiply, in different units. Bitwise was never available
  and the *claim* was the error. Split: `5a` matches the arithmetic and
  requires bitwise; `5b` **measures** the units route at **12 ULP** and reports
  the number, because every module here that converts metres to Å pays it.
- **Validation 4 demanded the asymptotic order at a finite step.** It required
  every successive error ratio within `0.05` of 2 and got
  `2.286, 2.145, 2.073, 2.036` — a textbook first-order sequence approaching 2
  from above. The fix changed the assertion's **shape**, not its tightness:
  monotone convergence to 2, which a second-order or wrong-order error fails,
  where a loosened tolerance would have passed anything.
- **Validation 3 silently skipped an inadmissible planted shape**, so it
  reported on two cases while claiming three. Refusals are named and counted
  now.
- **The cost table contradicted its own prose.** It tabulated `ell = 0.6 Å`,
  found `diverges = 0` on all nine rows, and the paragraph beneath said the
  divergences "cluster at the SMALL separations" — a statement about rows that
  were all zero. The sweep now covers `ell` and the claim is whatever the
  table says.

And a fifth, which only looking could catch: **the figure was unreadable
twice.** Admissibility only requires `P > 0` on `[2.0, 5.2] Å`, so a shape may
approach a zero of `P` just inside the right edge where `L' = P'/P` diverges;
one such shape drove the left panel's autoscale to `−100 1/Å` and flattened
every other curve onto the axis. Fixed, then the metal labels overprinted the
x-axis label with Ni clipped to read "li". Both views are clipped now with the
clipping **stated on the axis**, and nothing is dropped from the computation —
only from the view. Neither defect was reachable from the transcript, which
said `15 passed, 0 failed` throughout.
