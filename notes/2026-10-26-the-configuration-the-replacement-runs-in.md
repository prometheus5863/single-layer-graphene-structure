# The configuration the replacement runs in

**2026-10-26.** Discharging the 2026-10-25 ranked-first item — *complete the
read surface for every replacement in `rra.FLIPS`, not four holders* — turns
up a property of the **measurement** rather than a gap in coverage.

## 1. The question the item did not ask

A replacement is a function a flip manager *installs* at a site. It only ever
runs with that manager entered. But all three instruments this repository
measures with are relative to the namespace at the moment they are used:

| instrument | why it is configuration-relative |
|---|---|
| `csg.capture_surface` | a static walk; a bare global resolving to a repo function is followed, and *which* function a name resolves to is whatever is bound now |
| `vrs.value_reads` | the same walk, recording module-level scalars reached bare |
| `csg.leverage` | calls the holder in the current namespace |

Every committed leverage figure — 2026-10-22's, 2026-10-24's, 2026-10-25's —
was taken by calling the replacement **directly, with no flip entered**. That
is the one configuration in which the replacement never runs. Nothing said
whether it mattered.

## 2. The population, derived

`rra.FLIPS` × `rra.rebind_set` gives **8 (flip, site) pairs and 6 distinct
installed objects**:

| flip | site | installed object |
|---|---|---|
| `qca.on_consistent_argument()` | `crx.carrier_density` | `_carrier_density_consistent` |
| `qca.on_consistent_argument()` | `gfet.carrier_density` | `_carrier_density_consistent` |
| `qca.on_consistent_argument()` | `gfet.quantum_capacitance` | `_quantum_capacitance_consistent` |
| `blast.contact_on_consistent_argument()` | `cd.quantum_capacitance` | `_cq_consistent_contact` |
| `flip.on_default_flip(live_contact=False)` | `crx.carrier_density` | `_counting(_carrier_density_consistent)` |
| `flip.on_default_flip(live_contact=False)` | `gfet.carrier_density` | `_counting(_carrier_density_consistent)` |
| `flip.on_default_flip(live_contact=False)` | `gfet.quantum_capacitance` | `_counting(_quantum_capacitance_consistent)` |
| `flip.on_default_flip(live_contact=False)` | `cd.quantum_capacitance` | `make_contact_cq(live=False)` |

2026-10-25 audited four holders. **Two of the six — the arrival-counting
wrappers — had never had a read surface computed at all.**

## 3. The answer, and it is mostly good news

Measured over all 6 × 3 cells:

| object | own flip | other flips |
|---|---|---|
| `_carrier_density_consistent` | 14/14 `BITWISE_SAME` | invariant |
| `_quantum_capacitance_consistent` | 11/11 | invariant |
| `_counting(_carrier_density_consistent)` | 15/15 | invariant (surface grows 15 → 17 under `qca`, no leverage moves) |
| `_counting(_quantum_capacitance_consistent)` | 12/12 | invariant |
| `make_contact_cq(live=False)` | 19/19 | invariant |
| **`_cq_consistent_contact`** | **15/15** | **not invariant under either other flip** |

**All six objects are bitwise flip-invariant under their own flip** — every
name, every counterfactual, and the holder's own value. *That* is the
proposition every committed leverage number in this repository rests on, and
it had never been stated, let alone checked. It is now check `A1`.

A **surface** change is not the same as a **verdict** change, and the table
separates them: four of the six objects gain surface names in some
configuration while every counterfactual stays bitwise identical.

## 4. The one row that moves

`blast._cq_consistent_contact` under `flip.on_default_flip(live_contact=False)`:
value `0x1.0c835a6063492p-4` → `0x1.4ef081d146707p-7`, **a factor of 6.4134**,
and **11 of 15 names cross the exact-zero boundary**.

| name | surface | isolated | as installed | class |
|---|---|---|---|---|
| `blast.cq_shipped_of_n` | both | +58.3421 % | +51.8166 % | `BOTH_LIVE` |
| `blast.n_consistent` | both | +22.4740 % | +3.7310 % | `BOTH_LIVE` |
| `cd.C_interface` | both | +18.1154 % | +3.1971 % | `BOTH_LIVE` |
| **`gfet.C_ox`** | both | **exact 0** | **+7.8110 %** | `LIVE_ONLY_AS_INSTALLED` |
| `gfet.T` | both | +0.0034 % | +0.1725 % | `BOTH_LIVE` |
| `gfet._thermal_rounding` | iso | +58.3421 % | **exact 0** | `DEAD_AS_INSTALLED` |
| `gfet.e` | iso | **+278.4911 %** | **exact 0** | `DEAD_AS_INSTALLED` |
| `gfet.hbar` | iso | −63.5273 % | **exact 0** | `DEAD_AS_INSTALLED` |
| `gfet.kB` | iso | +0.0242 % | **exact 0** | `DEAD_AS_INSTALLED` |
| `gfet.quantum_capacitance` | both | +58.3421 % | +51.8166 % | `BOTH_LIVE` |
| `gfet.v_F` | iso | −63.5273 % | **exact 0** | `DEAD_AS_INSTALLED` |
| `qca.A_FERMI` | both | +58.3417 % | +41.7502 % | `BOTH_LIVE` |
| `qca.CQ_PREFACTOR` | inst | exact 0 | **+52.4103 %** | `LIVE_ONLY_AS_INSTALLED` |
| `qca.E_CHARGE` | both | −21.2518 % | −32.2411 % | `BOTH_LIVE` |
| `qca.KB` | inst | exact 0 | **+35.8863 %** | `LIVE_ONLY_AS_INSTALLED` |
| `qca._fixed_point` | both | +22.4740 % | +14.5902 % | `BOTH_LIVE` |
| `qca.cq_thermal` | inst | exact 0 | **+52.4103 %** | `LIVE_ONLY_AS_INSTALLED` |
| `qca.cq_thermal_of_n` | inst | exact 0 | **+52.4103 %** | `LIVE_ONLY_AS_INSTALLED` |
| `qca.n_arg_EF_thermal` | inst | exact 0 | **+9.0738 %** | `LIVE_ONLY_AS_INSTALLED` |
| `qca.series_factor` | both | +26.1720 % | +15.4460 % | `BOTH_LIVE` |

**The value change is not a new result.** It is 2026-10-21's *"the two halves
do not compose"*, which measured exactly this composition and costed it at
−51 % to −62 % in `n_contact`. What is new is that the same non-composition
appears as a change in the **read surface** — the object every guard here
demands adjudication over. The standing item *adjudicate the relation for the
flips 10-23 did not compose* is therefore partly discharged: the relation is
now measured over all 6 × 3 cells, and exactly one row is non-composing.

## 5. `gfet.C_ox`: a correct verdict with an over-claiming reason

2026-10-25 adjudicated `gfet.C_ox` **PROTECTED** on this holder, on the
strength of an exact bitwise zero, with the stated reason

> "the contact solve closes against `C_interface`, so the oxide capacitance
> reaches this value not at all … the zero is a property of the contact model
> and not of the 2026-10-21 repair."

The verdict is correct in the configuration it was measured in and is **not
withdrawn**. The reason is not: the first clause is right and the
generalisation is wrong, because the same holder gives `gfet.C_ox`
**+7.8110 %** under the default flip. The exact zero is a property of the
contact model *composed with nothing else*. The reason is corrected in place
and the superseded wording kept.

And then the sharper statement, `G9`, exhaustive over the vocabulary:

| verdict for `gfet.C_ox` | in isolation | as installed |
|---|---|---|
| `PROTECTED` | ok | refused (moves) |
| `UNREACHED` | refused (it is read) | refused (moves) |
| `EXPOSED` | refused (bitwise zero) | ok |
| `GOVERNED` | refused (no rebind set) | refused |

**Exactly one verdict is admissible in each configuration and they are
different ones.** So the four-verdict vocabulary cannot express the state of
this name, and non-composition is *not dischargeable by adjudication alone* —
which is why that standing item has survived three sessions of being noted.

## 6. `flip._counting` is a capture, and nobody had noticed

Each counting wrapper's surface contains exactly one name the replacement it
wraps does not: that replacement's **own name**. A rebind of it leaves the
wrapper's value **bitwise unchanged**, because `_counting` holds the function
object in a closure cell. So:

* this is a **third capture** in the repository, after the 2026-10-21 contact
  repair and the `bc.rebind` family;
* the channel replacement **name** is immune under the default flip;
* the verdict is `UNREACHED`, not `PROTECTED`, by `csg`'s own 2026-10-24
  argument for `CDQC` — the wrapper never *looks the name up* (closure cell,
  not a global read), so the bitwise zero is not a capture neutralising a
  read, and `PROTECTED` would credit the capture with work it did not do.
  **This module's first draft said `PROTECTED` and the 2026-10-24 guard
  refused it** — the fourth time a rule planted on an earlier day has caught
  a claim written on a later one.
* the wrapper returns a value **bitwise identical** to the replacement it
  wraps, and every shared name has **bitwise identical** leverage on both,
  which is what transfers 2026-10-25's verdicts onto the wrappers *with
  proof* rather than by assumption.

## 7. Exact anchors, not plausible ranges

Fifth consecutive session in which an exactly-known value, not a
plausible-range check, is what made the result trustworthy:

* `gfet.hbar` and `gfet.v_F` move the capture holder to the **bitwise
  identical** `0x1.87bc5bc0c26d0p-6` in all four configurations, because the
  `C_q` prefactor depends on them only through `(ħ v_F)²`. A probe reaching
  the value by any other route could not reproduce that.
* the counting wrapper's value is bitwise equal to the wrapped
  replacement's.
* **the repaired and unrepaired contact replacements are bitwise
  indistinguishable in isolation** (`0x1.0c835a6063492p-4` both). The entire
  content of the 2026-10-21 repair appears *only* as installed — so five
  sessions measured both holders in the one configuration in which they
  cannot be told apart.
* the identity control flip leaves the value and the surface bitwise
  unchanged: the instrument is not the perturbation.

## 8. And two of today's own checks were passing for the wrong reason

The mutation runner (9 defects, each required to be caught *by a named
check*) found both:

* **N4** — deleting the new `PROTECTED`-must-be-exact-as-installed refusal
  scored **26/26 undetected**. The first draft's `G6` planted its false
  verdict on a holder whose adjudication was incomplete *in that
  configuration*, so the missing-name refusal fired first and the branch
  under test was never entered. A control that was green and could not have
  been anything else. Sixth instance of *every verdict needs a planted
  control*, and the first where the unreachable control sat in the same file
  as the branch.
* **N8** — comparing the holder's own value with a relative tolerance instead
  of bitwise also went undetected, because `A4` asserted the two hex values
  and never the **invariance flag** computed from them. A factor of 6.4 could
  read as flip-invariant.

Both repaired in the audit module. Final: audit **28/28**, mutation **9/9 on
their own terms**, `N0` control clean.

## 9. What is NOT covered

* `live_contact=True`, and the *simultaneous* composition of two managers:
  three configurations are measured, not all of them.
* "invariant under its own flip" is a bitwise statement at one probe point
  (`csg.PROBE_V`) with one perturbation (`csg.PERTURBATION`), not a theorem.
* `vrs.value_reads` is still over-approximate and its negatives are still
  unsound (2026-10-25); nothing here repairs that.
* No live literature search was used this session — the subject is this
  repository's own instruments measured against themselves, and there is no
  external source for it.

## 10. Consequence for §7.9 item 22 (`v_F`)

Item 22 reads *"Chapter 4 quotes the contact `C_q` at a `v_F` its own
substrate does not have"*, on the −63.5273 % that a ×1.5 `v_F` carries. That
figure is measured on the **isolated** holders. It **survives** for the path
the thesis actually uses: on `flip.make_contact_cq(live=False)`, the object
the default flip installs, `gfet.v_F` carries −63.5273 % in *every*
configuration. It does **not** survive for `blast._cq_consistent_contact`
composed with either other flip, where `gfet.v_F` is bitwise dead. Item 22
(i) is unchanged in substance and gains a scope sentence; item 23 is created
for the configurations this session did not measure.
