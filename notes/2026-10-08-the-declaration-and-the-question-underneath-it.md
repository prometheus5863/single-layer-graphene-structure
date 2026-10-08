# The declaration, and the question that was underneath it

**2026-10-08.** Chapter 7 §7.9 items 7 and 10, both created or rewritten on
2026-10-07, closed together. Item 7 asked for a declaration; item 10 asked for
a census that would find out whether the defect item 7 names is widespread.
Doing them in that order was the right order, and the second one is the part
that produced a surprise.

---

## 1. The census: one defect in 86 names, and three predictions that were wrong

`graphene_potential_census_audit.py` collects every potential-like name in the
repository — 86 of them across 25 modules — and requires each to carry an
explicit adjudication in a registry inside the file. A name that appears in the
repository and is absent from the registry is reported as a fault in its own
right, so the census cannot go stale silently. That is the whole reason to
mechanise it rather than answer item 10 in prose.

The name pattern is deliberately too wide. It catches `mu_acoustic`, which is a
phonon-limited mobility and not a chemical potential at all. That is the correct
bias for a census: a pattern narrow enough to have no false positives has
already decided which names are potentials, which is the thing being audited.
So `MOBILITY_NOT_A_POTENTIAL` is one of the readings a name can be adjudicated
to, and 12 names carry it.

**The 10-07 entry named three places it expected siblings of the `V_ch`
defect. All three are checked and none of them is one, for three different
reasons.**

| 10-07's candidate | What the census found |
|---|---|
| Chapter 5's interconnect channel potential | **Does not exist.** `graphene_interconnect_model.py` contains no potential-like name anywhere. It is a resistivity-and-geometry model with no channel potential. |
| Chapter 6's two-contact photodetector channel potential | `V_bias` there is **imported and printed and enters no computation** in that module, so it has no second reading to be ambiguous between. |
| §4.4's local Dirac-point shift | `V_channel_shift` appears **only in the docstring derivation**; the closed form actually evaluated never computes it. |

The first of those is kept as a check that is *expected to fail*, and fails on
every run, because a pass there would mean a potential variable had appeared in
Chapter 5 — which is information, and the only way to keep it is to keep asking.

**Measured base rate: one defect in 86 names.** The 10-07 entry's guess that
the class was widespread is not supported. Seven names remain
`UNADJUDICATED` — recorded as open work, not as a pass — and none of them is
channel-swept, so none can carry the defect; what is undecided about them is
whether they are potentials at all.

## 2. The gate, and the thirteen false positives that produced it

The first run of the classifier reported **thirteen** dual-role names,
including `V_g`, `V_hi`, `mu_t` and `V_SAT_ANCHOR`. Every one was a false
positive, and instructively so.

`V_g` genuinely appears in a charge relation and in a current relation, in the
same function. That is **not** a defect. A terminal bias is the same physical
quantity in both expressions — the gate electrode's potential, fixed by a
supply — so there is no second reading available to it.

What made `V_ch` ambiguous was never that it appeared twice. It is that `V_ch`
is a **channel** variable: a potential swept from 0 to the drain bias *inside a
single current evaluation*. Only such a variable has two inequivalent
readings, electrostatic and quasi-Fermi, that differ by a diffusion term.

So "both roles" is necessary and not sufficient. The defect verdict now
requires both roles **and** channel-swept, detected mechanically by following
`profile = linspace(0, Vds, n)` through to `for V_ch in profile`. Both roles
without the sweep is reported as `DUAL_ROLE_TERMINAL` and is not a fault; ten
names sit there.

Mutant **M4b** exists so that gate is not decoration: it removes the drain-bias
name from the sweep in a copy of `graphene_fet_model.py` and requires `V_ch` to
*leave* the defect class, which it does. Mutant **M4a** is its complement — it
strips the declaration and requires `V_ch` to fall *back* to `UNDETERMINED`,
which proves the declaration is what cleared the defect rather than a change in
the classifier.

Two other first-version faults are worth recording because they are the
repository's own recurring classes:

- Role markers fired on **prose**: on a comment about `sheet_conductivity`,
  and on each function's own `def` line matching its own name. Comments,
  docstrings, string literals and signatures are now blanked before matching.
  2026-09-28 established that prose is a detector; it is not evidence of a data
  dependency.
- Three registry rows had **no code behind them**, because the scanner saw only
  `ast.Name` and missed every `gfet.V_dirac` — cross-module potential use
  through an `Attribute`. The dead-row check earned its keep on its first run.

## 3. The declaration

> **`V_ch` is the quasi-Fermi (electrochemical) potential of the channel
> carriers, in volts, measured from the source.**

Made in `graphene_fet_model.py`, where `V_ch` is defined, as an in-source
marker — following `graphene_dead_name_sweep.py`'s convention that the
repository's reasons belong in the repository and not in the instrument that
judges it.

Three reasons, none of them new work:

1. **The compact model Eq. (4) descends from says so.** Pasadas and Jiménez,
   *IEEE TED* 63(7) 2016 ([arXiv:1605.08235](https://arxiv.org/pdf/1605.08235)),
   write `v = μF` with `F = −dV/dx` and state that *"V(x) is the quasi-Fermi
   level along the graphene channel"*. This was already quoted on 2026-10-07;
   what was missing was binding this repository's variable to it.
2. **Classical long-channel theory uses the same variable and says what it
   buys.** H.-S. P. Wong's long-channel MOSFET notes
   ([EE 316 slides](https://dunham.ece.uw.edu/ee531/Long_Channel_MOSFET.pdf))
   define the channel variable `V(y)` as the *"Quasi-Fermi potential along the
   channel"* and head the current-density equation it feeds *"Current density
   equation (both drift and diffusion)"*. That is the entire content of the
   declaration: in the quasi-Fermi variable, **one term is the complete
   current**.
3. **The code already assumes it.** `carrier_density()` multiplies the
   electrostatic estimate by `C_q/(C_q + C_ox)`, which is the correction one
   applies when the channel variable is the quasi-Fermi level and the graphene
   drop `E_F/e` is carried separately. Under a strictly electrostatic reading
   that factor is itself a double count of the graphene drop.

### What the literature does *not* say, recorded because it matters

A review of GFET compact models (Lu, Wang, Li and Liu, *A review for compact
model of graphene field-effect transistors*,
[arXiv:1703.09759](https://arxiv.org/pdf/1703.09759)) was fetched specifically
to see whether the identification is standard. **It is not made there.** The
review's §2.5 calls the normalised imref splitting *"a.k.a. channel voltage"*,
but never connects that to the `V(x)` of its compact-model sections, which it
describes only as "the voltage along the channel"; its Eq. (12) carries a
diffusion term that its Eq. (15) then writes with the same expression as
drift, which is internally inconsistent and unexplained; and its §5 compact
model uses a drift form while saying only that "a drift-diffusion carrier
transport is assumed", never showing the diffusion contribution is absorbed.

So the ambiguity found in this repository on 2026-10-07 **is inherited from the
compact-model literature, not invented here.** That does not make it acceptable
in a thesis, and it does raise the value of declaring it: the declaration is
the step the review skips.

## 4. The consequence: 10-07's term is a double count

`graphene_diffusion_current_model.py`'s own Eq. (9) writes the current as
`d(V − λV_F)/dx`, and its own line 72 names `V − V_F` the quasi-Fermi
potential. Under the declaration, `V` already *is* that potential, so `λ = 0`
is the complete drift-diffusion current and `λ = 1` adds a term that is already
there.

**2026-10-07's diffusion term is therefore a double count, not a missing
0.50 %.** Measured three ways at `V_g` = 2.0 V, `V_ds` = 0.1 V:

| Level | `I_d` (λ=0) | `I_d` (λ=1) | difference |
|---|---|---|---|
| intrinsic (fixed channel drop) | 1.432116e-04 A | 1.442091e-04 A | **+0.6965 %** |
| terminal (contact feedback on) | 8.190217e-05 A | 8.225263e-05 A | **+0.4279 %** |
| 10-07, at its own saturated-peak bias | — | — | +0.5018 % |

The terminal number *must* be smaller than the intrinsic one, because the
series contact resistance is a negative feedback on `I_d`. That inequality is
checked rather than assumed, and it is the kind of check that is worth having:
its sign is fixed by the circuit and not by a fit.

Nothing is rewritten. The committed `λ = 0` currents stand, §4.6.5 keeps its
numbers, and `λ` keeps its one remaining use — it is now the instrument that
*measures* the difference between the two readings, which is the only thing a
knob between two models is good for.

## 5. The question underneath, which is larger than the one closed

Declaring `V_ch` quasi-Fermi makes a testable claim about `carrier_density()`.
Under the declaration the gate drive divides exactly:

```
V_g − V_dirac − V_ch  =  e·n/C_ox  +  E_F(n)/e ,    E_F(n)/e = A_F·√n    (D1)
```

a quadratic in `√n` with a **closed-form root** and no new parameter — `A_F` is
this repository's own `A_FERMI`. The repository instead computes the linearised
series factor

```
n = (C_ox·dV/e) · C_q/(C_q + C_ox)                                        (D2)
```

(D2) is the linearisation of (D1), so the declaration is self-consistent only
to the extent the two agree. **They do not agree as well as the term just
retired.** Running Eq. (4) on (D1) with nothing else changed:

| | `I_d` at the RF bias |
|---|---|
| on (D2), λ = 0 — committed | 8.190217e-05 A |
| on (D1), λ = 0 — declared exactly | 8.140885e-05 A |
| **the declaration costs** | **−0.6023 %** |
| the double count it removes | +0.4279 % |

So the honest headline is **not** "the 0.50 % is resolved". It is that a
question about the charge model, 1.4 times larger, was sitting underneath the
0.50 % and nothing in this repository had asked it. Chapter 4 keeps its (D2)
numbers today and gains a stated limitation with a measured size;
re-deriving on (D1) is a session of its own and is §7.9's new top item.

### A number that would have been overstated by a factor of 31

The raw (D1)-vs-(D2) charge disagreement reaches **+68 % at `dV` = 0.01 V**,
and quoting that would have been wrong: `n_puddle` = 5×10¹⁵ m⁻² swamps both
forms there, and no committed result lives in that regime. Compared in what
`carrier_density()` actually *returns* — floor included — the worst
disagreement over the drives Chapter 4 sweeps is **+2.22 %**, and in `I_d` it
is **−0.60 %**.

This is 2026-10-06's rule turned on today's own result: a comparison decided on
the quantity it *names* can mislead when that is not the quantity at stake. All
three numbers are reported in the transcript and the last one is the headline.

## 6. Validation against exactly known values

- (D1)'s closed-form root satisfies (D1) to **3.8×10⁻¹⁶** relative over
  `dV` ∈ [0.01, 3.5] V — an algebraic identity, not a plausible range.
- `dV = 0` gives `n = 0` **exactly** in both forms before the puddle floor.
- `n(+dV) == n(−dV)` **bitwise** in both forms, which both carry `n` as a
  magnitude.
- The large-drive limit is measured and monotone for both, not asserted
  (2026-10-03: a PASS at zero is silent about magnitude).
- `graphene_fet_model.py`'s own output is verified **unchanged** by today's
  edit: `transfer_characteristic` and `carrier_density` are bitwise identical
  to the previous commit's over a 41-point sweep. The edit is comments and one
  docstring annotation, and that is checked rather than assumed — 2026-10-02
  recorded that byte reproducibility is not correctness, but it is exactly the
  right instrument for a comment-only change.

## 7. A new defect found on the way, measured and left open

`quantum_capacitance()` returns a finite value up to `dV` = **18.3493 V** and
`inf`/`NaN` above it, because `log(2(1 + cosh η))` overflows once
`η = dV/(kT/e)` passes about 710 (`kT/e` = 25.85 mV here). The exact large-`η`
limit of that expression is `|η| + log 2`.

This is far outside the ~2.7 V maximum drive this thesis sweeps, so **no
committed number is affected**. It is still a defect, and it is left as a new
§7.9 item rather than patched today, because that module owns committed
transcripts and a numerics change there needs its own verification pass.

## 8. Two census rows adjudicated — one against the census's own guess

The census left `V_F` `UNADJUDICATED` in three modules on the guess that it was
a Fermi *velocity* written with a `V_` prefix. In
`graphene_diffusion_current_model.py` it is not: `fermi_voltage(n)` returns
`E_F/e` in **volts**, and it is precisely the quantity that converts between
the two readings of `V_ch`. In `graphene_contact_doping_nonlinear_model.py` the
same name **is** a velocity, 1.0×10⁶ m/s.

One name, two quantities, two modules — 2026-10-01's fault (a name that cannot
drift is the whole requirement) in a third place. Both rows now carry readings;
the third stays `UNADJUDICATED` because its module was not read today, and
saying so is cheaper than guessing again.

## 9. Methodological note, continuing the series

09-28: prose is a detector. 09-29: a mutation that does not arrive is
indistinguishable from a system that does not respond. 09-30: a control has to
sit where the failure enters. 10-01: and name what it compares against in a way
that cannot drift. 10-02: when it agrees, that is a fact about two artefacts.
10-03: a PASS/FAIL at zero is silent about magnitude. 10-04: a quantity can be
computed correctly under the wrong name. 10-05: a criterion with a number in it
can be missed in a way that locates the real constraint. 10-06: and can be
decided correctly on the quantity it names and still miss the finding, if that
quantity is a ratio. 10-07: an item can be well-posed, correctly motivated and
built exactly as written while its premise is an artefact of an unnamed
variable.

**10-08: AND CLOSING SUCH AN ITEM CAN REVEAL THAT THE AMBIGUITY WAS THE SMALLER
HALF OF THE QUESTION. Naming the potential settled ±0.50 % and exposed −0.60 %
that the naming question had been standing in front of.**

The checkable rule: **when a declaration resolves an ambiguity, measure what
the declaration itself costs in the quantity you quote.** A declaration is not
free — it promises that the rest of the model is the model that variable
belongs to — and the size of that promise is computable. Here it was computed
by monkey-patching one function and re-running Eq. (4), which took one line.
The mechanisable form is new §7.9 item 12.

And the second rule, from §1: **a predicted class is not a measured class.**
Three sibling defects were predicted on 10-07 from one instance. The census
found zero of them, and found the reasons to be three different ones. Predicting
where a defect class lives is cheap; the census that measures it cost one
session and is now a standing instrument.
