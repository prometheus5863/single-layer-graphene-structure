# Self-heating: the thermal network, and two exact benchmarks it must reproduce

**Date:** 2026-11-10
**Standing-list item discharged (in part):** *A SELF-HEATING TERM*, created
2026-10-06 from the Feijoo *et al.* reading, carried unchanged through
**thirty-four** sessions. It was created with a number attached — "self-heating
alone degrading `f_max` from 65 to 40 GHz" — and that number is what made it
unignorable: it bounds every *absolute* frequency this repository prints.

## 0. Why this item sat for thirty-four sessions, and what unblocked it

The item was created without a thermal model, and the obstacle was never the
physics — it was that a thermal resistance has to come from somewhere. A
self-heating term with an invented `R_th` would have been the 10-01 fault in a
new costume: a number with no provenance, propagating into `f_T`/`f_max`
and looking like a measurement. What unblocked it today is that the primary
source for this repository's high-field transport parameters **already
contains its own thermal model**, with the geometry it was applied to and the
share each term contributed. That makes `R_th` checkable rather than asserted.

## 1. Primary source

**V. E. Dorgan, M.-H. Bae and E. Pop, "Mobility and saturation velocity in
graphene on SiO2", *Appl. Phys. Lett.* 97, 082112 (2010)**, read as
[arXiv:1005.2711](https://arxiv.org/pdf/1005.2711) (preprint + supplement).

Fetched successfully as extracted PDF text. **Two expressions arrived
possibly garbled and both were resolved by arithmetic rather than by
assumption** — see §4; one resolved to an identity, one is recorded as an
ambiguity that the model now carries as a *sensitivity*, not as a choice.

Secondary/corroborating, not re-fetched today:

- **Feijoo, Pasadas, Bonmann *et al.*, "Does carrier velocity saturation help
  to enhance f_max in graphene field-effect transistors?", *Nanoscale Adv.* 2
  (2020), DOI 10.1039/c9na00733d** — the source of the 65 -> 40 GHz figure
  that created this item on 10-06. Already cited in the 10-06 entry.
- **Chen, Jang, Xiao *et al.*, *Appl. Phys. Lett.* 95, 161910 (2009)** —
  Dorgan's cited source for the graphene/SiO2 boundary conductance
  `h ~ 1e8 W m^-2 K^-1`. **This is a secondary citation here**: the value is
  taken as Dorgan quotes it and the Chen paper was not fetched. Labelled as
  such, per the standing 10-01 rule.

## 2. The three equations, as the paper states them

**Thermal network (Eq. 5).** A series resistance from channel to ambient:

    dT = T - T0 = P * (R_B + R_ox + R_Si)

    R_B  = 1 / (h * A)                 boundary (Kapitza) resistance
    R_ox = t_ox / (kappa_ox * A)       1-D conduction through the oxide
    R_Si = 1 / (2 * kappa_Si * sqrt(A))  spreading into the substrate
    A    = L * W

with `P` the dissipated electrical power, `h ~ 1e8 W m^-2 K^-1`, and
temperature-dependent conductivities (supplement Sec. 6):

    kappa_ox = 0.52 * ln(T_ox) - 1.687     [W/m.K],  T_ox = (T0 + T)/2
    kappa_Si = 2.4e4 / T0                  [W/m.K]

The paper solves `dT` **iteratively**, because `kappa_ox` depends on the
answer. Stated validity: `W, L >> t_ox`, and sample much smaller than the
wafer. Data with `dT > ~200 K` are **excluded from the paper's own fits**,
because heating past that point causes irreproducible behaviour (a `V_0`
shift). That exclusion is a usable ceiling and is carried into the model as a
reported regime boundary, not as a clamp.

**Mobility (Eq. 4).**

    mu(n,T) = mu0 / (1 + (n/n_ref)^alpha) * 1 / (1 + (T/T_ref - 1)^beta)

    mu0 = 4650 cm^2/V.s,  n_ref = 1.1e13 cm^-2,  T_ref = 300 K,
    alpha = 2.2,  beta = 3

Measured 300-500 K. Peak ~4500 cm^2/V.s at 300 K. The fit is **not** forced
below ~2e12 cm^-2 (their footnote 14). Note the temperature factor is
identically 1 at `T = T_ref` — used below as an exact check.

**Saturation velocity (Eq. 7).**

    v_sat(n,T) = (2/pi) * (w_OP / sqrt(pi*n))
                 * sqrt(1 - w_OP^2 / (4*pi*n*v_F^2))
                 * 1 / (N_OP + 1)

    N_OP = 1 / (exp(hbar*w_OP / (k_B*T)) - 1)

with `hbar*w_OP = 55 meV` (SiO2 surface phonon) and `160 meV` (graphene OP)
bracketing the data, and **~82 meV suggested by the paper itself for device
simulation**. Low-density ceiling `v_max = (2/pi)*v_F`. The paper reports the
model giving a **~20 %** `v_sat` decrease from ~280 to ~500 K for the 55 meV
branch and **~2 %** for the 160 meV branch — i.e. the temperature sensitivity
of `v_sat` is almost entirely a statement about *which phonon* you believe,
which is why 82 meV is carried as a parameter and swept rather than fixed.

Validity (their footnote 18): the square root is real only for

    n > (w_OP / v_F)^2 / (2*pi)

## 3. Two exact benchmarks, computed by hand before any code was written

This is the part that makes the item dischargeable. Both are external,
published, and exact — not plausible-range checks. On 2026-09-17 a
plausible-range check passed a real bug; these are the kind that does not.

**Benchmark A — the paper's own thermal budget.** Dorgan's device is
`W = 7 um`, `L = 4 um`, `t_ox = 300 nm`, so `A = 2.8e-11 m^2`. With
`kappa_ox(300 K) = 0.52*ln(300) - 1.687 = 1.2794` and
`kappa_Si = 2.4e4/300 = 80`:

| term | value (K/W) | share | paper's stated share |
|---|---|---|---|
| `R_B`  | 357.1 | 3.6 % | **4 %** |
| `R_ox` | 8378  | 84.5 % | **84 %** |
| `R_Si` | 1181  | 11.9 % | **12 %** |
| total  | **9916** | | **~1e4 K/W** |

All three shares and the total reproduce the paper's own quoted numbers. Any
implementation that does not hit this is wrong, and the agreement is on three
independent ratios, not one scalar — a single fudged constant cannot buy it.

**Benchmark B — the two validity thresholds.** `n_min = (w_OP/v_F)^2/(2*pi)`
evaluated at the paper's two phonon energies:

| `hbar*w_OP` | computed `n_min` | paper's footnote 18 |
|---|---|---|
| 160 meV | 9.405e11 cm^-2 | **9.4e11 cm^-2** |
| 55 meV  | 1.111e11 cm^-2 | **1.1e11 cm^-2** |

Two more exact hits, and this one independently pins `v_F`, `hbar` and the
factor of `2*pi` together — a sign or factor error in the `v_sat` prefactor
survives Benchmark A and dies here.

At the paper's suggested **82 meV** the same expression gives
`n_min = 2.47e11 cm^-2`, which is *below* this repository's operating
densities everywhere, so the branch cut is not a live constraint for us —
but it is now computed rather than hoped.

## 4. The two extraction ambiguities, and how each was settled

**`kappa_ox` resolved to an identity, not a guess.** The extracted text read
`kappa_ox = ln(T_ox^0.52) - 1.687`. That looked garbled. It is not:
`ln(x^0.52) == 0.52*ln(x)` identically, so the two spellings are the same
function, and it evaluates to 1.2794 W/m.K at 300 K — the accepted value for
thermal SiO2, *and* the value needed to reproduce the paper's 84 % share in
Benchmark A. Settled, with the arithmetic as the evidence.

> **ANNOTATED 2026-11-11. The identity is right. The arithmetic and the
> phrase "the accepted value" are both wrong, and this paragraph is where
> both errors entered the repository.**
>
> `0.52*ln(300) - 1.687 = 1.2789668868`, not `1.2794`. The slip is in the
> fourth decimal, 4.33e-4, and it was copied from here into the label of
> `graphene_self_heating_model` check A1 as `"accepted thermal-SiO2 value
> 1.2794 W/m.K"` — at which point an arithmetic slip in a note became a
> citation in a test.
>
> There is no accepted value of 1.2794 to be the accepted value of. Dorgan
> et al. print no numerical conductivity other than supp. Sec. 5's
> `kappa_ox ~ 1.3 Wm-1K-1 ... (at 300 K)`, and supp. Sec. 6 gives the fit
> with `0.52` and `1.687` as its only digits. The computed 1.2789668868
> rounds to 1.3 at the two figures the paper prints, so **the identity this
> paragraph settled was settled correctly** and nothing downstream of it is
> wrong by more than the slip.
>
> Nor was 1.2794 "the value needed to reproduce the paper's 84 % share":
> 1.2789668868 gives the same 84.5 %, which rounds to the same 84 %. The
> share was never evidence for the fourth decimal.
>
> The numbers in this note are left as written. The `R_ox = 3.52e5 K/W` of
> Section 5 below is unaffected — recomputing with 1.2789668868 gives
> 3.5185e5, which is 3.52e5 at three significant figures — and so is every
> number in the table above. See
> `notes/2026-11-11-a-census-of-comparisons-against-published-numbers.md`
> and `graphene_tolerance_census_audit.py`.

**The mobility temperature factor is a genuine ambiguity and is carried as
one.** The extracted Eq. (4) has the temperature factor *multiplying*:

    (i)  mu = mu0 / (1 + (n/n_ref)^a) * 1/(1 + (T/T_ref - 1)^b)     [product]

The other reading consistent with how such fits are usually written is a
single denominator:

    (ii) mu = mu0 / (1 + (n/n_ref)^a + (T/T_ref - 1)^b)             [sum]

Both are exactly 1 in the temperature factor at `T = 300 K`, so **Benchmark A
cannot distinguish them and neither can any 300 K number.** They differ off
room temperature: at 500 K and low density, (i) gives `0.771*mu0` and (ii)
gives `0.772*mu0` — nearly identical — but at 500 K *and* `n = n_ref` they
differ by ~13 %, because (i) multiplies the two suppressions and (ii) adds
them.

**Not settled by picking one.** The model implements both and the figure
reports the spread. A session that chose silently would have produced a
number indistinguishable from a measurement, which is the failure mode this
repository has spent two months learning to recognise. The discriminating
observation is Dorgan's Fig. 2b (mobility vs T at fixed `n = 3.5e12 cm^-2`),
which is a figure read, not a number quoted — so it is on the standing list
rather than resolved here.

## 5. What this says about *this* repository's device, in advance of the code

The device modelled in `graphene_fet_model.py` is `W = 1 um`, `L = 200 nm`,
`t_ox = 90 nm`, so `A = 2e-13 m^2` — **140x smaller in area than Dorgan's**.
Scaling the network:

    R_B  = 1/(1e8 * 2e-13)              = 5.00e4 K/W
    R_ox = 90e-9/(1.2794 * 2e-13)       = 3.52e5 K/W
    R_Si = 1/(2*80*sqrt(2e-13))         = 1.40e4 K/W
    total                               ~ 4.16e5 K/W

**42x the thermal resistance of the paper's device**, and the term ordering
is unchanged (oxide dominant) but far more lopsided: 85 % oxide, 12 %
boundary, 3 % spreading. Two consequences, both of which the code must report
rather than assume:

1. At the small-signal bias the RF chapter actually uses (`V_ds = 50 mV`),
   `I_d` is tens of microamps and `P` is a few microwatts, so `dT` is of
   order **1 K**. If that survives the self-consistent solve, then the
   10-06 worry — that the absolute `56.687 GHz` is uninsulated — is
   **narrowed**: it is uninsulated in principle and negligible at this bias.
2. At the bias where real RF GFETs are driven (`V_ds ~ 1 V`), the same
   network predicts `dT` of order hundreds of kelvin, i.e. **past Dorgan's
   own ~200 K reproducibility ceiling**. The interesting output is therefore
   not a correction factor but a *boundary*: the `V_ds` beyond which this
   device's own thermal environment puts it outside the regime its transport
   parameters were measured in.

**And the 1-D oxide term is an over-estimate here, which must be said out
loud.** Eq. (5)'s stated validity is `W, L >> t_ox`. This device has
`L = 200 nm` against `t_ox = 90 nm` — a ratio of 2.2, not `>>`. Heat spreads
laterally in the oxide over a length comparable to `t_ox`, so the true
`R_ox` is *smaller* than `t_ox/(kappa_ox*A)`. The consequence is directional
and therefore usable: **the network as written is an upper bound on `R_th`,
hence an upper bound on `dT`, hence a conservative (worst-case) bound on the
self-heating degradation.** A bound in a known direction is worth more than a
centred estimate of unknown sign, and it is the honest thing this geometry
supports. Tightening it needs a spreading-corrected `R_ox` — on the standing
list, not invented here.

## 6. What was not attempted, recorded so no session re-spends it

- **PubMed/PMC: not attempted**, per the standing reCAPTCHA note (three
  prior occurrences logged).
- **`pubs.aip.org`: not attempted today.** It 403'd on 10-06 for the *AIP
  Advances* row and is the known-403 publisher for this repository; the APL
  version of Dorgan is behind it, which is why the arXiv preprint was used.
  The preprint carries the supplement, which is where the thermal model
  lives, so nothing was lost.
- **Chen *et al.* (2009) not fetched** — `h` is a secondary citation, flagged
  in §1 rather than laundered into a primary one.
- **Dorgan Fig. 2b and Fig. 3c not read as data.** Both are figures, and both
  are what would settle §4's ambiguity and pin `v_sat(n,T)` numerically. That
  is a figure-reading task, which this session cannot do from extracted text,
  and it is now the specific blocking step for the product-vs-sum question
  rather than a vague wish for "better data".
