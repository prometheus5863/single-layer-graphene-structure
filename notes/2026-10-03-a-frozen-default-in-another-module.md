# A frozen default in another module, and a detector that only tests against zero

**2026-10-03.** Session type: audit + correction. Live web search was not
used; this is an internal-consistency session on code and numbers already in
the repository, and no claim below rests on a literature fetch.

Top open item closed: *"`require_delivery` is still not CALLED by any audit
— created 09-30, untouched for three sessions."*

---

## 1. The rule, obeyed, and the pass that is worth stating plainly

`graphene_mutation_arrival_probe.py` was written on 2026-09-30 and asks the
repository to adopt one rule:

```python
from graphene_mutation_arrival_probe import require_delivery
require_delivery('graphene_fet_model', ['mu', 'T', 'Rc_per_width_ohm_um'])
```

before interpreting any number produced by mutating those names. For three
sessions nothing called it. `graphene_cross_module_delivery_audit.py` now
does, over every mutation target `graphene_figure_provenance_audit.py`
actually applies — the only audit here that rebinds attributes and then reads
numbers. Nineteen mutation applications, fifteen distinct `(module, name)`
targets.

**It passes.** Eighteen targets are `LIVE`; the nineteenth,
`graphene_fet_model.t_ox`, is `IMPORT_USED` and its derived child `C_ox` is
co-patched in the same patch dict, which is the proviso the probe's
`import_consumed_is_covered` exists to check and which is now checked
mechanically rather than trusted. So none of the figure-provenance verdicts
were manufactured by a mutation that failed to arrive. Three sessions of the
rule going uncalled hid nothing.

That is a negative result and it is the honest headline of Section 1. The
finding is in what the pass *could not have covered*.

## 2. `classify()` asks about one module. Three of the mutations cross two.

`classify(M, N)` answers: *do the readers of `N` inside `M` read it at call
time?* That is the right question for a rebind of `M.N` **if every reader is
in `M`**. Three of the nineteen applications patch `graphene_fet_model.L`,
`.mu` and `.n_puddle` while harvesting a figure drawn by
`rf_small_signal_model`. For those, a `LIVE` verdict in the **defining**
module is silent about a frozen capture in the **consuming** one.

One such capture existed:

```python
# rf_small_signal_model.py, as committed until today
def gate_resistance(L=gfet.L, W=gfet.W, N_fingers=1):
    return R_sheet_gate * W / (3.0 * L * N_fingers ** 2)
```

`classify(graphene_fet_model, 'W')` returns `LIVE` — **correctly**. Every
reader of `W` inside `graphene_fet_model` is live. The frozen reader is a
reader of a *value that was copied out of that module before the audit
existed*, and it is not a reader of the name at all. No amount of care
applying the 09-30 rule as written would have found it.

## 3. It was not latent. It defeated the model's own override.

`compute_fT_fmax(..., W=W_RF)` exists to rescale the device from this
repository's 1 µm per-width convention up to a literature-representative
40 µm / 8-finger geometry, and its docstring says so. It does it by
**rebinding `gfet.W`**. `R_g` kept the 1 µm it had captured, so the gate
resistance entering both `f_max` denominators was too small by exactly
`W_RF / W = 40` — 0.208 Ω where it should have been 8.333 Ω.

| peak value, 40 µm / 8 fingers, V_ds = 0.1 V | as committed | corrected | factor |
|---|---|---|---|
| f_T intrinsic | 20.2365 GHz | 20.2365 GHz | bitwise identical |
| f_max intrinsic | 18.786 GHz | 13.094 GHz | 1.4347 |
| f_T extrinsic | 3.4350 GHz | 3.4350 GHz | bitwise identical |
| f_max extrinsic | 3.1885 GHz | 2.2196 GHz | 1.4365 |

`f_T` is untouched because `R_g` does not appear in
`f_T = |g_m| / (2π C_gs)`. The normalized 1 µm path — panels 1 and 2 of
`rf_figures_of_merit.png`, and the 20.279 / 18.731 GHz numbers Chapter 4 §4.6
quotes — is **bitwise** unchanged, verified against the pinned pre-fix blob
rather than asserted.

## 4. Three reasons it survived, and the third is the general one

**(a) The transcript printed the right R_g beside an f_max computed from the
wrong one.** `summary_numbers()` calls `gate_resistance(W=W_RF,
N_fingers=N_FINGERS_RF)` with `W` **explicit** and printed
`R_g (8-finger) = 8.33 Ohm`. `_compute_fT_fmax_core` called the same function
with the default and got 0.208 Ω. The printed diagnostic and the number it
was supposed to diagnose came from two different values of the same quantity,
two lines apart, for five weeks.

**(b) Every check in reach was a ratio.** The module cites Feijoo *et al.*'s
raw-vs-de-embedded gap as its sanity check on the pad-capacitance treatment.
That check is a ratio, and the defect was multiplicative in both terms:

| | as committed | corrected | change |
|---|---|---|---|
| intrinsic/extrinsic peak f_max | 5.8918 | 5.8994 | **0.13 %** |
| extrinsic/intrinsic f_max degradation | 0.1697 | 0.1695 | **0.1 %** |
| both f_max values individually | — | — | **≈44 %** |

A factor common to numerator and denominator divides out of the only quantity
being validated. **A ratio cannot see a defect it shares.** This is why the
model agreed with the literature on the one number anyone was checking while
both of its inputs were wrong by nearly half.

**(c) `MUST_CHANGE` is a test against zero.** The figure-provenance audit's
`N_FINGERS_RF x2` mutation acts through `R_g` alone — gate resistance goes as
1/N². With `R_g` frozen 40× too small, `R_g` barely mattered in the `f_max`
denominator, so the audit measured:

```
  max rel change, N_FINGERS_RF x2, as committed : 0.01033
  max rel change, same mutation, corrected      : 0.28620
  suppression                                   : 27.7x
```

**Both are PASSes.** The knob had lost 96 % of its authority and still read
as alive, because the only threshold was *did anything move at all*. The
frozen default did not merely corrupt the model — it corrupted the audit's
measurement of its own sensitivity, and the audit had no way to say so.

## 5. Two faults of my own, recorded rather than re-tuned

**E5's first form was mis-specified.** It asserted
`R_g(after rebind) == R_g(before) * 40.0` exactly, reasoning that 40 is
representable. It failed on the first run at the last bit
(`8.333333333333332` vs `8.333333333333334`). The premise was false:
representability of the *factor* says nothing about the *rounding sequence*.
`R_g = R_sheet·W / (3 L N²)`, so scaling `W` rounds once in the numerator,
whereas `r0 * 40.0` rounds an already-rounded quotient. **An exactness claim
has to name the operation it is exact under.** E5 is now an identity about
delivery — the late-bound default path must equal the explicit-argument path
*bitwise* — and E5b keeps the 1-ULP measurement that exposed the error.

**Section 3's first assertion was simply false.** It claimed `f_T` was
affected too. It is not, and the check failed on its own false premise. It is
now a *paired* `MUST_NOT_CHANGE` / `MUST_CHANGE`, which makes the scope of
the defect an assertion instead of a sentence of prose.

**And one caught by an existing instrument.** `BLOB_PREFIX_RF` in the new
module was written and never read. `graphene_dead_name_sweep.py` reported it
`DEAD` on the first run — the same fault the 10-02 session hit with its own
`SELF`, found by the same sweep, in the session that created the name. It now
has a real reader: `historical_claim_is_checkable()` reads the pinned pre-fix
blob out of the object store and requires this module's claim *about history*
to hold of it, reporting `SKIP` rather than `PASS` when the object store
cannot answer. A pinned reference with no reader is decoration.

## 6. The census, and what "dormant" does and does not mean

Ten cross-module captured defaults existed this morning; eight remain, all in
`graphene_sensitivity_audit.py` (`_calib_terms`, `_rho_terms`, capturing
`graphene_interconnect_model` constants). They are **dormant, not faults**:
that audit perturbs by passing keyword arguments explicitly and never rebinds
`icm.*`, so the captured values serve only as baselines. Calling them faults
would be the 10-01 error of letting a classifier's name drift from what it
measures.

Dormant means *no call site rebinds the name today*. One future `setattr`
converts all eight at once, which is why they are enumerated rather than
dismissed. The census also reports that two `setattr` calls in the
figure-provenance audit use a **computed** name, so it cannot be complete by
construction — said out loud, so that "0 faults" does not come to mean "none
possible".

## 7. Where this leaves Chapter 4

§4.6 is annotated in place; the superseded numbers are kept beside the
corrected ones. The uncomfortable part is recorded there too: with `R_g`
correct, this model gives **f_max/f_T = 0.647** at the literature-scale
device rather than 0.928, i.e. *further* from Feijoo *et al.*'s de-embedded
1.3–1.4, not closer. At 8 fingers the corrected `R_g = 8.33 Ω` is already in
the engineered-low range that paper describes, so the remaining shortfall
sits in `g_ds` and in the `R_g·C_gd` feedback term rather than in the gate
resistance. That split is not quantified and is now an open item.

---

**Methodological note, continuing the series.** 09-28: prose is a detector,
and a check that cannot fail is worse than no check. 09-29: a mutation that
does not arrive is indistinguishable from a system that does not respond.
09-30: a control has to sit where the failure enters. 10-01: and it has to
name what it compares against in a way that cannot drift. 10-02: and when it
agrees, that is a fact about two artefacts, not about the world.
**10-03: AND A PASS/FAIL AT ZERO IS SILENT ABOUT MAGNITUDE.** Every
`MUST_CHANGE` in this repository measures a response size, prints it, and
throws it away; the verdict retains only the sign. A knob can lose 96 % of
its authority — or, in the limit this session did not test, 99.9 % — without
a single check changing colour. The second strand is narrower and sharper: a
single-module instrument asked a cross-module question and returned the right
answer to the wrong question, and the only reason that was discoverable is
that somebody asked what the instrument's scope *was* rather than what its
verdict *said*.
