# The caveat was in the docstring and the claim was in the transcript, and only the transcript travels

**2026-10-18.** A literature pass on the top physics item — a *lower* bound on
`ℓ`, the two-sided half of the second anchor — which did not obtain the lower
bound, and found instead that the committed *upper* bound's headline number is
attributed to the source in the one artefact that does not carry the caveat.

Source throughout: Khomyakov, Giovannetti, Rusu, Brocks, van den Brink and
Kelly, *First-principles study of the interaction and charge transfer between
graphene and metals*, Phys. Rev. B **79**, 195425 (2009),
[arXiv:0902.1203](https://arxiv.org/abs/0902.1203); read today via
[ar5iv](https://ar5iv.labs.arxiv.org/html/0902.1203).

## 1. What was being looked for, and why it is not there

`graphene_per_metal_crossover_model` carries

```
Δ_c(d ; ℓ) = 0.9 eV · exp( −(d − 3.3 Å) / ℓ )        (N1)
```

with `ℓ` swept over `[0.3, 1.5] Å`. 2026-10-15 found the paper's second
separation — a calculation at `d = 5.0 Å` of which it says the crossover
"decreases to its expected value, `W_M ∼ W_G`, for large separations" — and
turned it into `Δ_c(5.0 Å) ≤ ε`, hence `ℓ ≤ ℓ_max(ε)`. That is one-sided. A
*lower* bound needs the opposite inequality, `Δ_c(5.0 Å) ≥ δ > 0`: positive
evidence that the 5.0 Å crossover is *not exactly* `W_G`, which only a
resolvable residual offset in Fig. 4's upper point set can supply. Hence the
standing route, "digitise Fig. 4's upper (5.0 Å) point set".

**That route did not open today, for an infrastructure reason and not a
physics one.** `arxiv.org` is refused by the egress proxy in *both* shells —
the device VM (`curl: (56) Received HTTP code 403 from proxy after CONNECT`)
and the cloud sandbox (`connect_rejected`) — so the PDF cannot be downloaded
into either place, and a figure cannot be read without it. The ar5iv HTML was
reachable and truncated again, for the **fourth** recorded time, ending mid-
derivation of `Δ_tr(d)` at `Δ_tr(d) = αN(d)z` — before **Eq. (7)**, which is
the equation Fig. 4's fit lines follow. So Eq. (7) remains unread after four
attempts across three sessions.

Two things were nevertheless recovered verbatim and are worth having.

**Eq. (4), the paper's own form for the short-range term:**

```
Δ_c(d) = e^(−γd) (a₀ + a₁d + a₂d²)
```

`γ, a₀, a₁, a₂` are **not reported numerically anywhere in the text** — which
is the fact that justified building (N1) in the first place, and the reason
four sessions of hoping for a number were hoping for something absent. (N1)
drops the quadratic prefactor; `graphene_per_metal_crossover_model`'s own
docstring already records this as "a real approximation", so it is not a new
finding, and it is restated here only because it bears on §3 below.

**The sentence the committed `ε` comes from:**

> "We see that while Eq. (1) holds for a relatively large separation
> `d = 5.0 Å` between the graphene sheet and the metal surface, there is a
> small deviation of ∼0.08 eV at the equilibrium separation `d ≈ 3.3 Å` which
> can be traced to a perturbation [of the graphene electronic structure by
> physisorption that cannot be described as a rigid shift]"

## 2. The finding: four artefacts carry the caveat, and the fifth is the one that travels

`ε = 0.08 eV` is borrowed. Everyone involved knew:

| artefact | says `ε` is borrowed? |
|---|---|
| `graphene_second_anchor_bound.py` docstring | **yes** — "borrowed from a sentence about a DIFFERENT quantity", "no committed number is changed on the strength of one borrowed epsilon" |
| `notes/2026-10-15-…md` | **yes**, immediately under the table, with the separation quoted |
| `AUTOMATION_LOG.md`, 2026-10-15 entry | **yes**, in bold under its table |
| thesis §7.9 item 2 | hedged — "the only small energy scale the paper attaches to 'small' nearby"; does not say which separation |
| **`second_anchor_bound_output.txt`** | **no — and it made the strongest claim of the five** |

Checked by grep: the committed transcript contained **zero** occurrences of
`borrow`, `different quantity` or `non-rigid`, while printing

* `ell <= 0.7024 A at the paper's own 0.08 eV scale.`
* `66.5 % of the committed swept range is excluded by the source.`
* and a table row marked `<- paper`.

So the qualification lived in four places, none of which is the artefact that
`graphene_pristine_transcript_audit` compares, that a reader of results looks
at, or that gets quoted — and the standing item has been quoting the
unqualified form ("ℓ ≤ 0.7024 Å at the paper's own 0.08 eV scale") for three
sessions.

**This is the repository's standing fault class running the other way.**
2026-10-17 found two docstrings that *denied* what their code did; this is a
docstring that is scrupulous about what its *output* overstates. Both are the
same defect — a claim and its qualification in different files, with only one
of them read — and 2026-10-01's formula for the dead-name class inverts to fit
it exactly: not a citation with no code behind it, but **a qualification with
no output behind it.**

## 3. And the caveat was itself too weak

All four prose artefacts call `0.08 eV` "a different quantity". The verbatim
sentence supports something stronger, in two independent ways:

1. **It is attached to the other separation.** `∼0.08 eV` is the deviation at
   `d ≈ 3.3 Å`. It is being used as a tolerance on `Δ_c` at `d = 5.0 Å`, the
   separation about which that sentence says Eq. (1) *holds*.
2. **It measures the failure of the rigid-shift picture, and `Δ_c` *is* a
   rigid shift.** The paper traces the deviation to a perturbation "that
   cannot be described as a rigid shift". A number quantifying how badly rigid
   shift describes 3.3 Å is not a conservative bound on a rigid-shift term at
   5.0 Å; it is a different kind of quantity, and the two defects compound
   rather than cancel.

The consequence none of the artefacts drew: **there is no privileged row.**
The family over `ε ∈ [0.02, 0.20] eV` is the result, and the single headline
number the SUMMARY reported should never have existed. That is what changed
today.

## 4. What changed, and what did not

**No number moved.** `ℓ_max(0.08 eV) = 0.7024 Å` and the 66.5 % exclusion are
unchanged; the numeric multiset of the transcript is a strict superset of its
predecessor's, verified by diff. What changed:

* a `PROVENANCE` section is printed, giving every sourced number its section,
  its verbatim quote, **the separation it applies at**, and what quantity it
  measures — so the transcript carries its own attribution;
* the family table's `<- paper` marker becomes `<- BORROWED, see PROVENANCE`;
* the SUMMARY's two false attributions are replaced, with the superseded
  wording **quoted in the transcript** rather than deleted;
* four new validations, of which one is the rule made mechanical: *no
  `print()` in the module may assert the source attribution unqualified, and
  quoting it as superseded is allowed only within the same string literal.*

**The first draft of that check failed, and it failed for the third time
today on the same cause.** It forbade a phrase anywhere in a `print()`
argument, which cannot distinguish *asserting* "the paper's own 0.08 eV
scale" from *quoting it as retracted* — and the new provenance section does
the second deliberately. A check that cannot tell those apart forces a choice
between keeping the record and passing the check, i.e. it pushes toward
deleting the retracted wording, which is the one move the annotation rule
exists to refuse. Scoping the supersession marker to the individual string
literal fixes it, and has the right shape: the retraction has to travel in the
same breath as the retracted words.

Three for three today: this repository's instruments fail first by answering a
weaker question than the one asked.

## 5. Still open

* **The lower bound on `ℓ`.** Unobtainable from the text — the paper's
  statement at 5.0 Å is a tilde with no tolerance. It needs Fig. 4's upper
  point set, which needs the PDF, which needs `arxiv.org` to be reachable.
  **This is now blocked by egress policy rather than by effort**, and that is
  a change in the item's status worth recording: three sessions carried it as
  "cheapest route: digitise Fig. 4", and the route is closed from both shells
  until the allowlist changes. A copy of the PDF placed in the connected
  folder would open it immediately.
* **Eq. (7)**, unread after four ar5iv truncations.
* **Whether `ℓ` is even the right parameter to bound.** Under Eq. (4) the
  single-exponential `ℓ` is a *chord* decay length over `[3.3, 5.0] Å`, while
  `w_cross_for_metal` evaluates `delta_c` at each metal's `d_eq` — 2.05–3.31 Å,
  i.e. mostly *below* the anchor and outside the chord. Not costed here, and
  flagged for a session of its own.
