# The citation was a section number, and the only false citation was a credit to CODATA

**2026-11-13.** The item ranked first on 2026-11-12 closes, and so does the
item ranked second. Dorgan is re-read from the source rather than from this
repository's notes, which confirms every transcription in the newest module and
contradicts nothing. The citation census then finds three defects, and the
first of them is in the instrument that commissioned it.

---

## 1. The item

The 2026-11-12 standing list ranked first:

> **A CENSUS OF CITATIONS, NOT TOLERANCES** — carried from 11-11, **now ranked
> first**, and today strengthened its case twice. 9k already said no AST-level
> instrument can catch a label naming a value its source never printed; 9l now
> adds a second class the AST cannot see — a width whose warrant is a
> *different* printed number from its centre. Both are citation facts, not
> arithmetic facts. The tree has 13 labels claiming a source and the worklist
> is mechanical to emit.

The two censuses ask different questions and the difference is the whole
design:

| | question | evidence |
|---|---|---|
| `graphene_tolerance_census_audit` | is this check's **width** the one the quoted precision licenses? | the Python AST |
| `graphene_citation_census_audit` | did the named source **print** the named number, at the precision the label spells it to? | a declared registry of what each source prints |

The second question cannot be answered by parsing Python, so it is not
answered by parsing Python. `SOURCES` is the registry: one entry per source
this tree cites by name, each printed value recorded **as the source spells
it**, with the place in the document and the date the document was read. The
registry is the evidence. The module is only the bookkeeping that forces each
label to meet it.

## 2. The source, read again

**V. E. Dorgan, M.-H. Bae and E. Pop, "Mobility and saturation velocity in
graphene on SiO2", *Appl. Phys. Lett.* **97**, 082112 (2010)**, read today from
[arXiv:1005.2711](https://arxiv.org/pdf/1005.2711), which carries the
Supporting Online Materials in the same file. This is a second, independent
read: the 2026-11-10 session built `graphene_self_heating_model` from its own
read of the same preprint, and the point of repeating it was to test the
module's transcriptions against the text rather than against the note that
recorded them.

What is printed, and where:

| value, as printed | what | where |
|---|---|---|
| `R_th ≈ 10^4 K/W` | device thermal resistance at 300 K | main text after Eq. (5); again in supp. Sec. 5 |
| `~2.8×10^-7 m^2K/W` | the same, per unit of device area | **the same sentence** |
| `~84%`, `~12%`, `~4%` | shares of `R_ox`, `R_Si`, `R_B` | main text after Eq. (5) |
| `κ_ox ≈ 1.3 Wm^-1K^-1 ... (at 300 K)` | oxide thermal conductivity | supp. Sec. 5, and **nowhere in the main text** |
| `κ_ox = ln(T_ox^0.52) - 1.687`, `T_ox = (T_0+T)/2` | the fit, and its argument | supp. Sec. 6 |
| `9.4 × 10^11 cm^-2` at `ħω_OP = 160 meV`; `1.1 × 10^11 cm^-2` at `55 meV` | validity floor of Eq. (7) | footnote 18 |
| `W = 7 μm, L = 4 μm, t_ox = 300 nm` | geometry | Fig. 1(a) caption; supp. Sec. 1 |
| `~20%` and `~2%` *decrease* between `~280 K` and `~500 K` | `v_sat` temperature dependence | main text after Eq. (7) |

**Two transcriptions were at risk and both survived.** The fit is printed with
`ln(T_ox^0.52)` rather than `0.52·ln(T_ox)` — identical, and the module's
docstring already says so — and its argument is the **mean** oxide temperature
`(T_0+T)/2`, which `thermal_resistances` already evaluates it at. At the 300 K
benchmark the distinction is exactly inert, because the mean of 300 and 300 is
300 with no rounding at all, so no A-series number depends on it; at 500 K the
two readings differ by 2.1 %, which is why it was worth checking rather than
assuming. And `κ_ox ≈ 1.3` really is in supp. Sec. 5 and not in the main text,
exactly as check A1's label claims.

**The module with the most published numbers per line has zero false
citations.** That is a result and section 7 of the new module states it as
one. It also refutes this repository's own prediction, below.

New information the read brought back that the tree does not yet carry: the
`~20 %` / `~2 %` statement is bracketed between `~280 K` and `~500 K`, and the
`~20 %` belongs to the SiO2-phonon (55 meV) branch while the `~2 %` belongs to
the graphene-OP (160 meV) branch. Check V3 tests only the sign.

## 3. Finding 1 — the published value was a section number

The tolerance census extracts a quoted value by matching a publication cue and
then taking the next number within forty characters. On check A1 of
`graphene_self_heating_model`:

> `A1 kappa_ox(300K) rounds to the 1.3 W/m.K Dorgan supp. Sec. 5 prints for the
> oxide at 300 K, at the 2 s.f. that value is printed to`

the cue that wins is `Dorgan`, and the next number after it is the **5 of
"Sec. 5"**. Since 2026-11-11 that census has been reporting A1 as a comparison
against a published value of 5, quoted to one significant figure, and printing
it in its table of the tree.

This is the same class as the footnote-ordinal bug that census fixed on the day
it was written — where `footnote` matched before `footnote-18` and the census
read the ordinal `-18` as a carrier density. The fix made `footnote-\d+` a
longer cue so the ordinal could not be read as the value. That does nothing
here, because this ordinal **follows the cue instead of preceding the number**.
And the guard the earlier bug left behind cannot see it either: check 8e rejects
a **negative** quoted value, and a section number is positive.

The repair is a detector, not a longer regex: a number preceded by a
sectioning word (`Sec.`, `supp.`, `Fig.`, `Table`, `Eq.`, `Ref.`, …) is an
`ORDINAL` and is not adjudicated against anything. On this tree it fires on
exactly one label, which is asserted rather than hoped for.

The check A1 performs is sound — it tests `_rounds_to(kappa_ox(300), 1.3, 2)`
in the code. The misparse is in the census's reading of the label, and the
value A1 actually cites is the `1.3 W/m.K` three words earlier.

## 4. Finding 2 — the one false citation, and it credits CODATA

`graphene_fet_signed_carrier_model`:

> `this module reproduces Section 3.3's quoted 6.45 kOhm/sq from CODATA`

CODATA prints `e` and `h`. **6.45 kΩ/□ is `1/(4e²/h)` computed from them, by
this repository, in Chapter 3.** No publication printed it, and the warrant
for its three significant figures is a number this tree wrote itself.

This is the class check 9k said no AST instrument could reach, and it is the
only instance in the tree. Its consequence was live: the tolerance census
adjudicated that row **VACUOUS by ×0.5** against a three-figure publication
interval that no publication licenses. The arithmetic was right, the quantity
was right, and the interval it was judged against was imaginary.

The repair is a wording — cite Chapter 3 for the quantity and CODATA for `e`
and `h` — and it is deliberately **not applied today**, because the two
censuses are cross-validated against each other's populations this session and
editing the label in the same commit would make that validation vacuous.

## 5. Finding 3 — a real citation, spelled to a precision its source did not print

`graphene_areal_thermal_resistance` check I1 writes Dorgan's thermal resistance
as `1.0e4`. The paper prints `10^4`. Two figures claimed, one printed — and the
difference in licensed width is exactly a factor of ten, `1/2` against `1/20`,
asserted as exact rationals rather than compared as floats.

Nothing is wrong today, because I1 is a `LIT_NO_TEST` row: it tests an identity,
not a tolerance, so the inflated precision licenses nothing. It is the same
defect A2's literal had before 2026-11-12 — `1.0e4` in the code where the label
said `~1e4` — one layer further out, in the **label**, which is where the
tolerance census takes its evidence from.

## 6. Two things that are not defects, and one that is unscorable

**`DERIVED`.** "Dorgan's 2.8e-11 m²" appears nowhere in Dorgan. It is `L·W`
from the printed geometry, and the possessive is doing work the paper did not
do. That is a different warrant class, not a weaker one — and the distinction
has teeth: because the geometry is *stated* and not fitted, the derivation
inherits no uncertainty, which is exactly why 2026-11-12 could show the paper's
two numbers are one statement divided by an exact constant. Checked here to 0
ulp.

**`AMBIG_PREC`.** V3 quotes `~20 %`, whose spelling fixes no significant-figure
count. Dorgan does print `~20%`, so the citation is real and only its precision
is unscorable. The census says so and invents nothing — and a mutant that
guesses two figures instead of refusing is detected on that control.

**One row does not identify its own warrant.** A3's "the paper's 4 %" matches
two printed values of Dorgan at once: the `~4 %` share of `R_B` and the `4 μm`
channel length. The census matches on **value**, not on meaning, so it scores
that `CITED` and counts it. That count is a measurement of the ambiguity, not a
dismissal of it, and it is the second stated limit of the new instrument.

## 7. The second warrant: A2 from OVER_TIGHT ×28 to MATCHED

The item ranked second was the repair of the census's first false positive, and
it is the same sentence of Dorgan twice over. A2's centre is the one-figure
`~1e4`; its width is the two-figure `~2.8e-7` **from the same sentence**. The
census read only the centre, licensed ±50 %, saw ±1.7857 % and returned
`OVER_TIGHT` at ratio ×28.

The repair is the one the item predicted — a second optional token in the label
vocabulary, `[width-warrant: 2.8e-7]`, spelled machine-readably because reading
a precision claim out of English is the whole defect of the first version. A2
becomes `MATCHED`. **Nothing in the comparison changed**: not the tolerance, not
the value, not the outcome; only which printed number the census reads the
licensed width off.

It is honoured on relative tolerances only. A relative width is dimensionless
and transfers between quantities related by an exact constant; an absolute width
lives in the units of its own quantity, and transferring it would require the
label to assert the two quantities share units, which a label cannot be trusted
to do.

And the new verdict is a statement about the equality **band**, not about
equality:

    A2's coded tolerance   0.05/2.8      = 0.01785714285714285962
    the census's own width  0.5e-8/2.8e-7 = 0.01785714285714285615
    float(Fraction(1,56))                 = 0.01785714285714285615

Two spellings of `1/56`, one ulp apart. At `rel_eps = 0` the verdict is
`VACUOUS` **by that one ulp** rather than `MATCHED`. Third session in three
where a one-ulp fact decided something.

## 8. Scoring yesterday's predictions

**(1) Wrong, and wrong in the direction it argued against.** It said the
citation census's first finding would be in `graphene_self_heating_model`
rather than in an older module, on the reasoning that volume of citation and
not age predicts where a census bites. The only `UNWARRANTED` row is in
`graphene_fet_signed_carrier_model`; `self_heating` has none. The reasoning was
half right in a way the prediction did not claim: `self_heating` produced the
only `ORDINAL` and the only `RETRACTED` row — both of the citation anomalies
that are *not* false citations — which is what a high volume of citation per
line actually predicts.

**(2) First half right, second half wrong and unscorable as posed.** It said
the split-warrant repair would turn `MATCHED`, not `OK`, on A2, because A2's
width equals the licensed half-width exactly. It does, and §7 shows that
equality holds to one ulp and not better. It also said the repair would expose
at least one row currently adjudicated `MATCHED` as mis-adjudicated. The tree
contained **zero** `MATCHED` rows outside the instrument modules yesterday, so
the population it named was empty before the repair ran — the X4 class for the
third time — and A2 is the first `MATCHED` row this tree has ever had.

**(3) Not exercised.** The areal band at the repo's own `t_ox` was not built.

## 9. Two self-caught errors, and one correction in the older module

- **The mutation harness patched a copy of the module it was testing.** It did
  `import graphene_citation_census_audit as M` and patched attributes on `M`.
  Run as `python3 thisfile.py` the running module is `__main__`, so the import
  built a *second* module object from the same source; every mutation landed on
  the copy while the controls went on calling the original, and all five
  mutants reported NOT DETECTED. This is the 2026-09-30 item — *arrival, not
  response: the null control cannot certify its own mutation* — in a new
  disguise. Check 3e is now the arrival check for the harness itself, and it
  asserts that an import by name under `__main__` does yield a different
  object.
- **The new module planted the pattern it censuses.** Its own check labels
  quote registry tokens, because its subject *is* the quoting of registry
  tokens, so its arrival put ten of its planted controls into the tolerance
  census's table of **the tree** — including three `ORDINAL` and two
  `UNWARRANTED` rows. That is the 2026-10-17 standing item, one module over,
  and a `SELF`-only exclusion cannot see it: the excluded thing is not "this
  file" but "a module that plants the pattern it censuses", and there are now
  two. The exclusion is now the declared class `PLANTING_MODULES`, and checks
  8h/8i measure what it removes instead of asserting it is harmless.
- **`_check_labels`' docstring was describing a tolerance the code does not
  have.** It listed `_c` among the check-function names the tree uses
  interchangeably; the code requires a name ending in `check`, which `_c` is
  not. A planted control written with `_c` on the docstring's authority
  censused to zero rows with no error. Corrected in place, and the measurement
  recorded: no module in the tree uses `_c` for a check, so nothing was being
  missed.

## 10. Count

`graphene_citation_census_audit.py`: 51 checks, 51 passing, nine verdicts each
with a planted control, 5/5 mutants detected.
`graphene_tolerance_census_audit.py`: 69 → 82 checks, 0 failing; tree verdicts
now `VACUOUS` 1, `MATCHED` 1, `OVER_TIGHT` 0.
`graphene_self_heating_model.py`: 47/47, unchanged except for A2's label.

Two items closed, three defects found, one of them in the commissioning
instrument. One source fetched and every transcription it was fetched to test
confirmed.

## Sources

- [V. E. Dorgan, M.-H. Bae and E. Pop, *Mobility and saturation velocity in graphene on SiO2*, Appl. Phys. Lett. 97, 082112 (2010); arXiv:1005.2711](https://arxiv.org/pdf/1005.2711)
- [Y. Zion, A. Butenko, A. Sharoni, M. Kaveh and I. Shlimak, *Charge carrier transport asymmetry in monolayer graphene*, arXiv:1801.03565](https://arxiv.org/abs/1801.03565) — the registry's second source, read 2026-10-29/30 and not re-read today.
