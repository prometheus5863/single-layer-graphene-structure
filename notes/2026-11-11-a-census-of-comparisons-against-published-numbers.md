# A census of comparisons against published numbers, and a citation that had nothing behind it

**2026-11-11.** Discharges the ranked-first item of 2026-11-10: *census every
comparison against a published number and check it uses a significant-figure
test, not a relative tolerance*. The census is
`graphene_tolerance_census_audit.py` (60 checks, 0 failed, six mutants
detected). What it found, and what reading the source paper then did to the
finding, are the subject of this note.

---

## 1. What a comparison against a published number actually claims

A paper that prints `1.1e11 cm^-2` is not claiming a value. It is claiming an
interval. Two significant figures assert only that the true value lies in

```
[1.05e11, 1.15e11)
```

and nothing finer. So a check written as

```python
abs(got - 1.1e11) / 1.1e11 < 0.01
```

tests a precision the source never claimed — it rejects most of the band the
paper vouches for. That is the defect found on 2026-11-10 in check B2 of
`graphene_self_heating_model`, which failed at 1.027 % on a value, 1.1113e11,
that rounds to 1.1e11 exactly as required.

The generalisation is arithmetic. For `q` quoted to `s` significant figures,

```
h_abs = 0.5 * 10**(floor(log10|q|) - s + 1)        h_rel = h_abs / |q|
```

and a check is sound exactly when its tolerance equals that half-width, in the
tolerance's own units.

### Where the census corrects the item that commissioned it

The item was phrased as *significant-figure test, **not** a relative
tolerance* — as though the operator were the defect. It is not. A relative
tolerance set to `h_rel` **is** the round-to-`s`-figures test, up to the
half-way tie-break. The defect is a tolerance whose *number* was not derived
from the quoted precision, and every hand-written tolerance in this tree
qualifies, because none of them was.

This matters practically. A session that mechanically rewrote every relative
tolerance as `_rounds_to()` would have fixed three instances by accident and
left the class untouched. Verdicts are therefore reported as the ratio
`licensed / tested`, not as a pass/fail on which function was called.

---

## 2. What the tree contains

102 tolerance comparisons across 80 modules. Seven of them are claims about a
source; the rest are internal-consistency checks and are counted but not
adjudicated.

| module:line | form | tolerance | quoted | s.f. | licensed | verdict |
|---|---|---|---|---|---|---|
| `graphene_fet_signed_carrier_model:648` | `abs(a-b)` | ±0.01 | 6.45 kΩ/sq | 3 | ±0.005 | **VACUOUS ×2** |
| `graphene_self_heating_model:A1` | `abs(a-b)` | ±5e-4 | 1.2794 W/m·K | 5 | ±5e-5 | **VACUOUS ×10** |
| `graphene_self_heating_model:A2` | `abs(a-b)/b` | ±5 % | ~1e4 K/W | 1 | ±50 % | **OVER_TIGHT ×10** |
| `graphene_self_heating_model:A3,A4,A5` | `round(x)==n` | — | 4 %, 84 %, 12 % | — | — | sound |
| `graphene_kappa_impurity_bound:D2` | `float('%.2g'%x)` | — | 1.5e12 cm^-2 | 2 | — | sound |
| `graphene_self_heating_model:V3` | *sign only* | — | ~20 % | amb. | — | **no numeric test** |

Three observations.

**The `R_th` row is over-tight by exactly ten.** Dorgan prints `R_th ~ 10^4
K/W` — one significant figure, with a tilde, i.e. an order-of-magnitude claim
licensing ±50 %. The check tests ±5 %. It passes today at 9915.6 K/W with room
to spare, so the defect is latent; a model change moving `R_th` to 1.1e4 would
fail a check against a paper that would still agree with it.

**`V3` names a published magnitude and tests only its sign.** Dorgan's
"~20 % / ~2 %" statement about the temperature dependence of `v_sat` appears
in the label; the check asserts
`saturation_velocity(1e17, 500) < saturation_velocity(1e17, 300)`. The sign is
right. The magnitude has never been compared to anything.

**The significant-figure test already existed here, in three spellings.**
`_rounds_to()`, `round(x) == n`, and `float('%.2g' % x)`. The third is in
`graphene_kappa_impurity_bound`, which predates `_rounds_to` by five weeks,
and the first version of this census scored it `LIT_NO_TEST` — the worst
category available — because it only knew the first two spellings.
Under-reporting a correct check as an absent one is the same defect as
over-reporting a tolerance.

---

## 3. The finding: a vacuous tolerance passing a false equality

Check A1 of `graphene_self_heating_model` read

```python
_check("A1 kappa_ox(300K) == accepted thermal-SiO2 value 1.2794 W/m.K",
       abs(kappa_ox(300.0) - 1.2794) < 5e-4, ...)
```

The formula gives `kappa_ox(300) = 0.52*ln(300) - 1.687 = 1.2789668868`, which
rounds to **1.2790** at five figures, not 1.2794. The residual, 4.3311e-4, sits
between the licensed half-width (5e-5) and the coded one (5e-4). So the check
was green while the equality in its own label was false, and a tolerance set to
what the quotation licenses would have failed it.

This confirms the 2026-11-10 prediction (1) in its strong form: the census
would find a tolerance *vacuous* — loose enough to admit a wrong value — rather
than merely over-tight like B2. It did, and the wrong value was already
through.

---

## 4. Then reading the paper retracted the question

The first commit of the day called *which of the two numbers is wrong* an open
question, and offered an intercept of 1.686567 or a coefficient of 0.52007593
as candidate roundings of "either published coefficient".

Both branches are wrong, because **1.2794 is not a published number.**

> Dorgan, Bae and Pop, *Mobility and Saturation Velocity in Graphene on
> SiO2*, Appl. Phys. Lett. **97**, 082112 (2010);
> [arXiv:1005.2711](https://arxiv.org/pdf/1005.2711),
> [publisher PDF with SOM](https://poplab.stanford.edu/pdfs/Dorgan-GrapheneVsat-apl10.pdf)
>
> - **supp. Sec. 6:** "These can be approximated as `κox = ln(Tox^0.52) −1.687`
>   and `κSi = 2.4×10^4/T0` by simple fitting to the experimental data in
>   Refs. 10 and 11 respectively." `0.52` and `1.687` are the only digits
>   given. No units note, no stated range of validity. `T_ox` is specified
>   elsewhere as `(T0+T)/2`.
> - **supp. Sec. 5:** "`κox ≈ 1.3 Wm-1K-1` is the oxide thermal conductivity
>   (at 300 K)". This is the **only** numerical conductivity the paper prints.
> - **main text:** "At 300 K for our geometry `R_th ≈ 10^4 K/W`, or
>   `~2.8×10^-7 m^2 K/W` per unit of device area", with `R_ox` ~84 %, `R_Si`
>   ~12 %, `R_B` ~4 %, and `W = 7 μm`, `L = 4 μm`, `t_ox = 300 nm`.

There is no 1.2794 in the paper at any precision. The candidate roundings were
roundings of coefficients no source quotes beyond the two digits the module
already uses, and the phrase "accepted thermal-SiO2 value" named an authority
that does not exist.

### Where it came from

`notes/2026-11-10-...md` §4, the paragraph that resolved the `ln(T^0.52)`
ambiguity. It states `0.52*ln(300) - 1.687 = 1.2794` — an arithmetic slip in
the fourth decimal — and then calls that "the accepted value for thermal
SiO2". The next day the slip was copied into A1's label, where the words "the
accepted value" turned it into a citation. **An arithmetic error in a note
became a citation in a test in one step, and nothing in the tree could have
caught it.** The note has been annotated in place; its numbers are not deleted,
because none of them is wrong by more than 4.33e-4 and none of its conclusions
changes.

The note also claimed 1.2794 was "the value needed to reproduce the paper's
84 % share". It was not: 1.2789668868 gives 84.49 %, which rounds to the same
84 %. The share was never evidence for the fourth decimal.

### What survives

The identity `ln(x^0.52) ≡ 0.52·ln(x)` was resolved correctly, and the
transcription was never in doubt. `kappa_ox(300) = 1.2789668868` **does** round
to the 1.3 the paper prints, at the two figures the paper prints it to. Two
significant figures license ±0.05 — a band a thousand times wider than the
±5e-5 the phantom five-figure value appeared to license. A1 is now that
comparison, and it passes.

What was wrong was never the physics. It was the warrant.

---

## 5. The limit of this instrument, stated plainly

The census decides whether a comparison is a literature claim, and to what
precision, **by reading the check's label**. It therefore inherits every false
citation in the tree and cannot detect one. `1.2794` reads as five significant
figures and `1.3` as two; both readings are correct, and only one of the two
values has a source. Catching that needed the paper, not the AST.

This is now check 9k of the audit rather than a sentence in a note, so that a
clean census is not mistaken for verified citations.

---

## 6. Five defects the census found in itself

Each is fixed with the measurement that caught it recorded at the fix.

1. **Every tolerance was compared against the *relative* half-width** — wrong
   for `abs(a-b) < tol`, which lives in the quantity's units. The 6.45 kΩ/sq
   row was reported vacuous by 12.9 when the true factor is 2. Over-stating a
   defect is the same error as understating it.
2. **Only inline residuals were matched**, so the commonest shape in this tree
   — `rel = abs(a-b)/b` bound to a name, then checked — was invisible. **44 of
   102 rows, 43 %,** are reached only through a local name binding.
3. **The cue list included "reported"**. On the real tree, 16 of 29 cued labels
   — **55 %** — were cued by that word alone, and in this repository "reported"
   means reported *by an instrument in this tree*. "accepted" needed the phrase
   "accepted … value" for the same reason: bare, it runs at 33 % precision
   here. Requiring the phrase is what brought A1 into the census at all.
4. **"footnote" preceded "footnote-18" in the cue alternation**, and a number
   was allowed to begin immediately after a cue, so on check B2 the census read
   the *ordinal of the footnote*, −18, as the published value — a negative
   carrier density quoted to two significant figures, about to be adjudicated.
5. **The mutation harness tested nothing.** It did `import
   graphene_tolerance_census_audit as self_mod`, which loads a *second copy* of
   the file when it is run as a script, so every monkeypatch landed on a module
   nothing was using. All five mutants were "detected" only because the planted
   expectation for P8 was itself wrong. A mutation suite reporting 5/5 while
   testing nothing is the 2026-09-17 failure again, inside the tool written to
   prevent it. `sys.modules[__name__]` is the fix; with it, the mutants visibly
   change the planted verdicts.

And one that is not a defect: **P8 of the planted control was planted as
VACUOUS and is MATCHED.** 2.5e3 written to two figures licenses exactly 2 %, so
the 2 % tolerance is right to the last bit. The planted expectation was wrong
and the instrument was right — the only way round that is any use. It is left
at 0.02 with the arithmetic recorded rather than nudged to make the plant look
prescient.

---

## 7. Two standing items picked up instances

Both were found by *running* the census, not by reading for them.

**A MODULE MUST NOT PLANT THE PATTERN IT CENSUSES** (created 2026-10-17). The
first full run reported three `OVER_TIGHT` tolerances and two `LIT_NO_TEST`
rows that exist only inside this module's own planted controls and check suite
— one of them quoting a "published" value of `0`, lifted out of the label of a
check asserting that a quoted token of 0 licenses nothing. The tree had
acquired three tolerance defects that live inside the instrument reporting
them. The module is now excluded by default, and checks 8f and 8g *assert* the
self-plant rather than describing it.

**GIVE EVERY VERDICT IN EVERY VOCABULARY A PLANTED CONTROL** (nineteenth
instance). It arrived from the vocabulary side rather than the verdict side:
the three spellings of the significant-figure test in §2, and the two senses of
"accepted" and "reported" in §6.3. Every one of the seven verdicts now has a
planted instance, and the plant is run against six mutants.

---

## 8. Scoring the 2026-11-10 predictions

1. *"The census finds at least one more relative tolerance compared against a
   two-or-three-significant-figure literature number, and at least one of those
   is vacuous rather than over-tight."* — **Right on the substance, wrong on
   the form.** Two further defective tolerances were found and one is vacuous
   and already admitting a wrong value, which is the strong form of the claim.
   But neither new row is a *relative* tolerance against a *two-or-three*-figure
   number: the vacuous one is an **absolute** tolerance against a **five**-figure
   value, and the over-tight one is relative against a **one**-figure value. The
   prediction guessed the shape of the defect from the single instance it had.
2. *"The `v_F` item's count is now four names."* — **Not exercised.** Untouched
   today; carried.
3. *"`_rounds_to` will be the thing imported, and the first module to import it
   will be `graphene_standing_list_audit.py`, not a physics module."* —
   **Half right.** `_rounds_to` and `_round_sig` were indeed the imports. The
   first importer is `graphene_tolerance_census_audit.py`, which did not exist
   when the prediction was made and is an audit module rather than a physics
   one — so the prediction's *reasoning* holds and its named target does not.

---

## 9. What this leaves open

- **`R_th` per unit area.** The paper prints `~2.8×10^-7 m^2·K/W` and the tree
  has never checked it. One significant figure, so it licenses ±1.8 %; it is an
  independent check on `R_th` and the area together, where the existing A2 tests
  their product only.
- **The `V3` magnitude.** Dorgan's ~20 % / ~2 % statement is testable and only
  its sign is tested.
- **Where 1.2794 is used elsewhere.** `grep` finds it only in the annotated note
  and in A1c. Nothing else in the tree reads it.
- **The tolerance census over the *other* repository.**
  `Design-Verification-Roadmap` has its own self-checking suites with
  hand-written tolerances, and `graphene_tolerance_census_audit.census_repo()`
  takes a `root`.

## Sources

- Dorgan, Bae & Pop, "Mobility and Saturation Velocity in Graphene on SiO2",
  *Appl. Phys. Lett.* **97**, 082112 (2010) —
  [arXiv:1005.2711](https://arxiv.org/pdf/1005.2711) ·
  [publisher PDF incl. SOM](https://poplab.stanford.edu/pdfs/Dorgan-GrapheneVsat-apl10.pdf)
