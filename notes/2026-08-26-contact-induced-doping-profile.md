# Research notes: metal-work-function-dependent doping profile under graphene contacts

*Date: 2026-08-26. Follow-on to
`notes/2026-08-21-contact-resistance-and-quantum-capacitance.md`, which
modeled metal-graphene contact resistance as a single lumped, literature-
calibrated series resistance (`Rc_total` in `graphene_fet_model.py`). That
lumped model does not capture *why* contact resistance depends on the
specific metal used, nor the fact that the contact also locally dopes the
graphene channel over a finite spatial extent, forming an in-plane
p-n/p-p'/n-n' junction. This is the gap flagged as "not yet covered" in
`AUTOMATION_LOG.md` since the 2026-08-21 run ("Metal contact work-function-
dependent doping profile ... currently only a lumped series resistance, not
a spatially resolved p-n junction model"). WebSearch was available and used
for this research.*

## 1. Physical origin: work-function mismatch drives charge transfer

When a metal is deposited on graphene, the two materials equilibrate to a
common Fermi level. Because their work functions generally differ, charge
flows across the interface until the resulting electrostatic dipole brings
the bands into alignment. Intrinsic (undoped, free-standing) graphene has a
work function of approximately 4.5 eV. If the contacting metal's work
function is *higher* than graphene's, electrons flow from graphene into the
metal, leaving the graphene under (and near) the contact hole-doped
(p-type); if the metal's work function is *lower*, electrons flow the other
way and the graphene is electron-doped (n-type) (Giovannetti et al.,
*Doping graphene with metal contacts*, Phys. Rev. Lett. 101, 026803 (2008),
https://www.researchgate.net/publication/23230469_Doping_Graphene_with_Metal_Contacts).

A first-principles/continuum treatment by Khomyakov, Giovannetti, Rusu, van
den Brink et al. (*Nonlinear screening of charges induced in graphene by
metal contacts*, Phys. Rev. B 82, 115437 (2010), arXiv:0911.2027,
https://arxiv.org/abs/0911.2027) found that:

- The crossover metal work function separating n-type from p-type
  contact-induced doping is **≈5.4 eV** — notably *higher* than graphene's
  own 4.5 eV work function, because weak (physisorptive/van-der-Waals-like)
  bonding at the interface shifts the effective crossover away from the
  naive "equal work functions" expectation.
- Consequently, **most common contact metals (Ti, Cr, Al, Cu, Ag) n-type
  dope graphene**, while only the highest-work-function metals commonly
  used in device fabrication (**Au, Pt**, and to a lesser/more variable
  extent **Pd** depending on interface chemistry) p-type dope it.
  Reported Fermi-level shifts for weakly-bonded metals (Al, Ag, Cu, Au,
  Pt) are on the order of **~0.5 eV**, well below the shift a naive
  full-charge-transfer estimate from the raw ΔWF would predict — a direct
  signature of graphene's own low density of states (quantum capacitance)
  limiting how much charge it can absorb before its own Fermi level moves
  to relieve the potential difference. This is the same C_q-limits-induced-
  charge physics already used in `graphene_fet_model.py` for the gated
  channel, now applied to the ungated contact interface.

## 2. Spatial extent: the contact-induced doping is not confined to "under the metal"

Because graphene's density of states vanishes at the Dirac point, its
screening is anomalously weak (`Reads and writes charge` review;
Khomyakov et al. 2010 above). The induced electrostatic potential (and
hence the induced doping) does not stay localized under the metal pad —
it leaks into the exposed channel region beyond the contact edge. Khomyakov
et al. give the asymptotic decay of the induced potential with distance x
from the contact edge as:

- **x^(-1/2)** for *undoped* graphene far from the contact, and
- **x^(-1)** for *doped* graphene (i.e., once a finite carrier density,
  and hence finite screening, is established) —

both much slower than the exponential (Debye-like) screening expected in a
normal 2D metal, and consistent with reports that contact-induced doping
"extends hundreds of nanometers" beyond the nominal contact edge (Khomyakov
et al. 2010; also see Nagashio & Toriumi's contact transfer-length
literature already cited in the 2026-08-21 notes, Section 2, which
independently observed that the *electrical* contact transfer length scale
for graphene is likewise several hundred nm — physically consistent with
this doping-profile picture, since both trace to the same weak-screening
mechanism).

## 3. Consequence: an additional, metal-specific resistance component

A metal contact therefore does not just add a fixed interface (tunneling/
transfer) resistance — it also creates an adjoining **in-plane junction**
(p-p', n-n', or p-n depending on the channel's own doping state) in the
first few hundred nm of graphene beyond the contact edge. Two consequences
relevant to this thesis's device models:

1. **Metal choice affects total device resistance beyond the tabulated
   Rc value.** The literature Rc values already used in this thesis
   (Pd ≈110 Ω·µm, Ni ≈470 Ω·µm, etc., from the 2026-08-21 notes) are
   *lumped, metal-specific measured quantities* that implicitly include
   some of this junction effect, but do not let us predict how the
   *channel-side* doping profile (and hence local sheet resistance) varies
   with gate voltage relative to the fixed contact doping — that requires
   an explicit spatial model.
2. **The junction is asymmetric and gate-tunable.** As the channel gate
   voltage sweeps the bulk channel doping through the Dirac point (as in
   `graphene_fet_model.py`'s transfer characteristic), the contact region
   stays pinned near its metal-set doping level (only weakly gate-tunable,
   since the metal itself screens the back-gate field near the contact
   edge — a further consequence of the same weak in-plane screening).
   This means the *type* of junction (p-p'/n-n' unipolar vs. p-n bipolar)
   changes as V_g sweeps, which is a first-order effect for photodetector
   operation too: graphene photodetectors deliberately exploit a contact-
   doped p-n junction near the metal edge to spatially separate
   photogenerated carriers (see `thesis_draft/06-graphene-photodetectors.md`,
   Section 6.5's flagged prerequisite, "spatially resolved collection
   model" — this is precisely the missing ingredient that today's model
   partially addresses).

## 4. Model plan for `graphene_contact_doping_model.py`

Following the compact-analytic-model style already established in this
repo (`graphene_fet_model.py`):

- Reuse the existing quantum-capacitance formula (`quantum_capacitance()`
  in `graphene_fet_model.py`) to compute the *contact-edge* Fermi level
  shift self-consistently, but with the oxide capacitance term replaced by
  an effective, much larger "interface capacitance" (using a sub-nm
  effective separation appropriate for a direct/weakly-bonded metal-
  graphene contact rather than a ~90 nm SiO2 back-gate dielectric). This
  reproduces the qualitative literature result that graphene's own quantum
  capacitance — not the interface geometry — is the dominant charge-
  limiting element right at the contact.
- Represent the spatial decay away from the contact edge with a
  saturating power-law profile, `f(x) = 1 / (1 + x/lambda)`, which
  reduces to unity at the contact edge and falls off as ~lambda/x for
  x >> lambda, matching the Khomyakov et al. doped-graphene asymptotic
  (x^-1) result; lambda is set to a representative "few hundred nm" scale
  consistent with the reported doping extent.
- Compute the resulting extra sheet-resistance contribution by integrating
  the (position-dependent) sheet conductivity along the junction region and
  comparing to the resistance the same length of channel would have at the
  bulk (gate-set) channel doping — i.e., isolate the resistance
  *specifically* attributable to the contact-doping profile, as distinct
  from the already-tabulated lumped Rc.

## References

1. Giovannetti, G. et al. "Doping Graphene with Metal Contacts." *Phys.
   Rev. Lett.* 101, 026803 (2008).
   https://www.researchgate.net/publication/23230469_Doping_Graphene_with_Metal_Contacts
2. Khomyakov, P. A., Giovannetti, G., Rusu, P. C., Brocks, G., van den
   Brink, J., Kelly, P. J. "Nonlinear screening of charges induced in
   graphene by metal contacts." *Phys. Rev. B* 82, 115437 (2010).
   arXiv:0911.2027. https://arxiv.org/abs/0911.2027
3. "Effect of Graphene Doping Level near the Metal Contact Region on
   Electrical and Photoresponse Characteristics of Graphene
   Photodetector." *Sensors* 20(17), 4661 (2020).
   https://www.mdpi.com/1424-8220/20/17/4661 (PMC copy:
   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7506932/) — confirms the
   device-relevance of contact-region doping level for both electrical
   *and* photoresponse characteristics, directly motivating the link
   drawn above to Chapter 6.
4. `notes/2026-08-21-contact-resistance-and-quantum-capacitance.md` (this
   repo) — lumped Rc literature values and quantum-capacitance formula
   reused here.
