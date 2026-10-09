# The sign was computed and then thrown away, and the floor was never a floor

**2026-10-09.** Session notes for Chapter 7 Section 7.9 item 9 (a signed-carrier
FET model), which turned out to also close Section 4.9's floor item and Section
7.9 item 11. Module: `graphene_fet_signed_carrier_model.py`. Transcript:
`fet_signed_carrier_output.txt`. Figure: `fet_signed_carrier_model.png`.

---

## 1. Why item 9 had to come before item 12

The 10-08 session declared `V_ch` to be the quasi-Fermi potential and, in doing
so, wrote down the exact charge relation that Chapter 4's linearised form
approximates:

    V_g − V_dirac − V_ch = e·n/C_ox + E_F(n)/e                          (4.28)

It then measured that re-deriving on (4.28) would move `I_d` at the RF bias by
−0.6023 %, made that item 12, the top item — and ordered item 9 ahead of it
with the reason that re-deriving first "would rebuild the branch-sign flaw in a
new equation."

That reason is correct and it is sharper than 10-08 stated it. **Both terms on
the right of (4.28) are odd in `n`.** `e·n/C_ox` obviously; `E_F(n)/e =
A_F·sign(n)·√|n|` because the Fermi level of a hole gas is below the Dirac
point. A relation whose two sides are odd has a root that is odd — and a root
that is odd cannot be computed by a function whose first line is

    dV = np.abs(dV)

which is the first line of `graphene_potential_declaration.n_quasi_fermi()`.
That function returns the right *magnitude*. What it cannot return, and what
nothing in the repository could return before today, is **which branch the
channel is on at a given point**. Chapter 4 would have been re-derived on an
equation that is defined on two branches, using a solver that knows about one.

`n_net_exact()` is that solver with the branch restored. It is the same closed
form — the magnitude solves `a·x² + b·x − |dV| = 0` on *both* branches, which
is why the 10-08 version got the magnitude right — with `sign(dV)` carried
through. Two checks, neither a plausible range:

- `n(−dV) == −n(+dV)` **bitwise** over 801 points.
- substituting the root back into (4.28) leaves a residual of **7.690e-16**
  relative over `dV ∈ [−2, 2] V`. An algebraic identity, in the sense 10-08
  used for its own 3.8e-16.

**Item 12 is unblocked.** That is the deliverable; everything below is what
turned up while building it.

## 2. The decomposition, and the identity nobody had noticed

Near neutrality a disordered graphene sheet is not nearly empty, it is an
electron–hole puddle landscape carrying both species at once. Write the net
density as signed and split it:

    n_net = n_e − n_h                  (electrostatics, signed, odd in dV)
    S     = √(n_net² + n_puddle²)
    n_e   = (S + n_net)/2,   n_h = (S − n_net)/2

Three things follow, and none of them is fitted:

**(i) `n_e + n_h` IS the shipped magnitude.** `carrier_density()` returns `S`.
So the quantity Chapter 4 has been flooring since August is the **total
conducting density** — not the net density, and not a per-species density.
Section 4 below is entirely dependent on this identification.

**(ii) `n_e·n_h == n_puddle²/4`**, to 3.603e-16 relative. The quadrature floor
was chosen in August for smoothness — it rounds the conductivity minimum
without a kink — and it is *algebraically* a constant electron–hole product,
i.e. a mass-action law with no free parameter. That was not the reason it was
chosen, and it is worth recording that a convenience choice turned out to have
a closure's structure.

**(iii) `n_e(+dV) == n_h(−dV)` bitwise.** Graphene's ambipolar symmetry becomes
machine-checkable rather than something one verifies by looking at a figure.

### What the literature actually does, which is not quite this

[Wiedmann et al., *Coexistence of electron and hole transport in graphene*,
Phys. Rev. B **84**, 115314 (2011); arXiv:1107.3929](https://arxiv.org/pdf/1107.3929)
measures both species through the Hall coefficient of a two-carrier
(compensated-semiconductor) model, their Eq. (2) with equal mobilities:

    1/R_H = e·(n + p)² / (n − p)

and takes `n − p = q ≃ α·V_g` from electrostatics — **the same split this
module uses**: net from the gate, total from disorder. They then note that the
divergence of `1/R_H` at the neutrality point "implies that `n + p` must remain
finite," which is an independent, measurement-side argument that **the floored
quantity is the total**. That is direct support for (i).

**It is not support for (ii), and this should not be overstated.** Wiedmann et
al. impose no electron–hole product; `n` and `p` are *outputs* of their fit, not
a closure. The constant-product law above is therefore a **stronger assumption
than the literature makes**, inherited from a smoothness choice, and it is
logged as a limitation rather than a result. Fetch recorded honestly: the paper
was retrieved and read, and it says something adjacent to what was wanted
rather than what was wanted.

Their sample B gives `n(q=0) = p(q=0) ≈ 4.2 × 10¹⁴ m⁻²`, i.e. a **total**
puddle density of ≈ 8.4 × 10¹⁴ m⁻². Chapter 4's `n_puddle` = 5 × 10¹⁵ m⁻² is
**5.95×** that. Noted, not acted on: that is a different sample of unstated
substrate quality and an external number, where Section 4's comparison is
internal to this thesis and therefore the stronger finding.

## 3. Three checks failed as first written, and all three were findings

This is the part of the session worth keeping.

### 3.1 The ambipolar symmetry failed "bitwise", by 2.0 /m²

Written as `n_e(+dV) == n_h(−dV)` with the bias passed as `V_dirac ± x`, it
failed at a maximum absolute difference of 2.0 /m². The symmetry is **not
broken**: stated in `dV` directly it is bitwise, residual exactly 0.

What is not exact is the *argument*. `(V_dirac − x) − V_dirac` is not
`−((V_dirac + x) − V_dirac)` in binary floating point — worst case 2.220e-16 V,
at x = −2.0 V, where 0.8 − 2.0 − 0.8 returns −1.9999999999999998. The
`C_ox/e` lever is 2.395e15 /(V·m²), and 2.22e-16 × 2.4e15 ≈ 0.5, which through
the series factor becomes the 2.0 /m².

Two consequences:

- **Every sweep in this repository is a `linspace` over `V_g` followed by a
  subtraction of `V_dirac`.** Any BITWISE symmetry claim about any of them
  carries this artefact, and such a claim has to name its parameterisation.
  `n_net_of_dV()` exists so that symmetry claims can be made in the variable
  the model is a function of.
- The residual is **accounted for, not thresholded**. Feeding the reconstructed
  `dV` values into the `dV`-level model reproduces the `V_g`-route species
  bitwise, which proves all 2.0 /m² belongs to the argument and none to the
  model. A tolerance of "≤ 2 ulp" would have passed too, and would have been a
  number tuned until the test went green.

This is 09-30's lesson moved upstream: a control has to sit where the failure
enters, and here the failure entered through the *input construction*, not
through the mechanism under test.

### 3.2 "n_e − n_h recovers n_net" failed at 1.772e-14, and passes at 1.5e-16

Same residual, same data, **119× apart**, because the denominator differs. The
absolute error is at most 1.0 /m², which is exactly **one ulp of `S`**.
Measured against `S` that is 1.5e-16 — a pass. Measured against `n_net` it
reaches 1.772e-14, because `n_net → 0` at neutrality while `S` stays pinned at
the puddle floor.

So the choice of denominator decided PASS from FAIL. **That is Section 7.9's
10-06 item — "census every criterion and check that is a ratio" — occurring in
a check written today, by the same process that is supposed to be auditing for
it.** The item has been open for three sessions; this is the first instance of
it generated *by the auditor rather than found in the audited code*.

Design rule, now in the module: **never reconstruct the net density from the
species pair.** Call `n_net()` or `n_net_of_dV()`. The species pair carries the
net density only to the precision of the total, which near neutrality is no
precision at all.

### 3.3 The σ identity failed "bitwise" at 3 ulp, and was not repaired by widening

`e·μ·(n_e + n_h)` versus `e·μ·S` differ by at most 1.626e-19 S/sq = 3 ulp,
because `0.5(S+n) + 0.5(S−n)` reassociates the sum and the two multiplications
add roundings.

The repair was **not** to raise a ulp budget from 2 to 4. 09-29's α anchor
settled the policy: when an exactness check fails, state the identity being
checked and name the slack, rather than loosening a number until it passes. The
criterion is now 1e-12 relative, justified as six orders below the four
significant figures Chapter 4 quotes and eleven orders above float noise — a
statement about the chapter, not about the test. The measured slack, 3.609e-16,
is printed next to it.

10-02 established that byte identity is not correctness. Today adds the
converse: **absence of byte identity is not a defect**, and three ulp of
summation order is not a physics claim.

## 4. Section 4.9's floor item, open eleven sessions, answered

`n_puddle` is the total (Section 2(i)), so at `dV = 0` the split gives
`n_e = n_h = n_puddle/2` **exactly** and

    σ_min(Ch. 4) = e·μ·(n_e + n_h) = e·μ·n_puddle = 3.204353e-04 S/sq
                 → 3.1208 kΩ/sq

Chapter 3 Section 3.3's *measured* floor, recomputed here from CODATA rather
than re-typed from the chapter text so the comparison cannot drift:

    σ_min(Ch. 3) = 4q_e²/h = 1.549618e-04 S/sq → 6.4532 kΩ/sq

(Section 3.3 quotes 6.45 kΩ/sq; this module reproduces it to 0.0032 kΩ/sq.)

**Chapter 4 is 2.0678× more conductive at neutrality than Chapter 3.** The
`n_puddle` that would reconcile them is 2.417989e+15 m⁻² = 2.418e11 cm⁻²,
against the shipped 5e11 cm⁻².

Quoted as a factor first, per 10-03: a PASS/FAIL on "do the two floors agree"
returns FAIL and is silent about whether the gap is 2× or 200×.

**Why this needed item 9 to be decidable.** `n_puddle` enters as a quadrature
floor on a *magnitude*, and a magnitude cannot say whether the floor is on the
net density, the total, or each species. Those three readings differ by factors
of 1, 1 and 2 in `σ_min` — the same order as the disagreement itself. So before
Section 2 the gap was not quotable at all. And the convention cannot be the
explanation: under the per-species reading the total would be `2·n_puddle` and
Chapter 4 would be **4.1357×** too conductive, i.e. the discrepancy moves the
**wrong way**.

**What is not claimed.** `4q_e²/h` is a quantum/ballistic minimum conductivity
observed experimentally; `e·μ·n_puddle` is a diffusive disorder floor. They are
not required to be equal. What *is* a defect is that this thesis quotes both,
in chapters that feed each other, and had never compared them — and that the
diffusive one comes out **below** the measured one, i.e. Chapter 4's channel at
neutrality is less resistive than any measured graphene sheet.

### What the factor costs, intrinsically and at the terminals

At `V_ds` = 0.1 V, `V_g,on` = 3.5 V, off-state taken at `V_g = V_dirac`:

| | shipped `n_puddle` | Ch.3-matched | on/off multiplier |
|---|---|---|---|
| terminal (R_c = 600 Ω) | on/off 1.2408 | 1.7820 | **1.4362×** |
| intrinsic (R_c = 0) | on/off 1.6146 | 2.8024 | **1.7357×** |

`I_off` falls by 35.19 % at the terminals and 51.58 % intrinsically. The
2.0678× floor error reaches the on/off ratio as 1.7357× intrinsically and only
1.4362× at the terminals: **the contacts attenuate the error by a further
1.2085×.** This is 10-08's negative-feedback structure — series contact
resistance attenuating an intrinsic change on the way to a terminal reading —
replicated on an entirely different parameter, which makes it a property of the
device topology rather than of the charge relation.

Nothing is rewritten. Chapter 4 keeps its numbers; the comparison is stated and
Section 4.6 is annotated in place.

## 5. Item 11 answered, by an assertion that was expected to pass

Section 5 of the module was written with a guard that `I_on` is essentially
untouched by the floor, on the authority of `carrier_density()`'s own
docstring — the floor "regularizes `n` at the Dirac point", "regularizes the
conductivity minimum".

**It failed. Changing only the floor moved `I_on` by −6.92 %.**

The reason, measured:

| `V_g` (V) | \|n_net\| (m⁻²) | `n_puddle`/\|n_net\| | floor adds to total |
|---|---|---|---|
| 0.90 | 2.357e+14 | 21.22 | +2023.9 % |
| 1.50 | 1.672e+15 | 2.990 | +215.2 % |
| 2.00 | 2.870e+15 | 1.742 | +100.9 % |
| 2.50 | 4.067e+15 | 1.229 | +58.5 % |
| 3.00 | 5.265e+15 | 0.950 | +37.9 % |
| **3.50** | **6.462e+15** | **0.774** | **+26.4 %** |

At `V_g` = 3.5 V — the largest overdrive swept anywhere in this thesis — the
"residual" puddle density is **77 % of the gate-induced density** and still adds
**+26.4 %** to the total. The floor's contribution falls below 1 % only at
`dV` = **14.73 V**, which is **5.46×** anything this thesis sweeps.

So Section 7.9 item 11, "whether `n_puddle` removes other regimes", created
10-07: **it does not remove regimes, it adds a floor to all of them.** It is not
a near-Dirac regularizer, it is a parallel conduction channel carried at every
bias, and the regime in which the docstring's description is true is a regime
the thesis never enters.

**09-28 said prose is a detector. Today adds: so is a sanity check you expected
to pass.** The docstring was the claim, the assertion was built on it, and the
assertion is what measured it.

## 6. Item 13 produced a wrong number on its first use outside its own module

10-08 found that `quantum_capacitance()` returns `inf`/`NaN` above
`dV` = 18.3493 V, made it item 13, and **ranked it last on purpose**: latent, no
committed number affected.

The first draft of Section 5b's 1 %-crossover bisection read
`contrib = ... if nnv > 0 else 1e9`. Above the overflow, `C_q` is `inf`,
`C_q/(C_q + C_ox)` is `inf/inf = nan`, `nan > 0` is **False**, so the NaN branch
was scored as "the floor dominates", the bisection walked away from the root,
and it returned its own 200 V bracket ceiling as the answer. The correct answer
is 14.73 V.

The bracket is now capped at 0.95 × the overflow onset and NaN is a hard stop
rather than a score. **Latent and harmless are not the same property**: item 13
affected no committed number and produced a wrong one within minutes of first
being used from outside the module that owns it. Its ranking is unchanged — it
is still a numerics change inside a module with committed transcripts across
three chapters — but "no committed number affected" is no longer the whole of
its description.

## 7. Methodological note, continuing the series

09-28: prose is a detector. 09-29: a mutation that does not arrive is
indistinguishable from a system that does not respond. 09-30: a control has to
sit where the failure enters. 10-01: and name what it compares against in a way
that cannot drift. 10-02: when it agrees, that is a fact about two artefacts.
10-03: a PASS/FAIL at zero is silent about magnitude. 10-04: a quantity can be
computed correctly under the wrong name. 10-05: a criterion with a number in it
can be missed in a way that locates the real constraint. 10-06: and decided
correctly on the quantity it names while missing the finding, if that quantity
is a ratio. 10-07: an item can be well-posed, correctly motivated and built
exactly as written while its premise is an artefact of an unnamed variable.
10-08: and closing such an item can show the ambiguity was the smaller half.

**10-09: AN EXACTNESS CHECK THAT FAILS IS A MEASUREMENT OF WHATEVER IT WAS
WRITTEN ON TOP OF.** Three failed today. One measured the floating-point
behaviour of the repository's sweep parameterisation; one measured the
denominator of its own criterion; one measured a docstring. None measured the
thing it was pointed at, and none of the three was a defect in the new code.

The checkable rule: **when an exactness check fails, find out what it measured
before deciding what to change.** The three repairs available at the moment of
failure were: loosen the tolerance, fix the model, or find out. Only the third
produced the three findings above, and the first would have silently accepted
all three. Mechanisable form: new Section 7.9 item 15 — for every tolerance in
this repository, record what fixed it at that value, and flag any that was set
after a failure.

Second rule, from Section 4: **a parameter whose meaning is fixed by its sign
convention cannot be compared against anything until the convention is
declared.** The floor item sat eleven sessions not because the calculation was
hard — it is one line — but because the quantity in it had three readings
differing by the same factor as the answer. 10-01 found the fault was in the
naming; today's variant is that the fault was in the *signature*.
