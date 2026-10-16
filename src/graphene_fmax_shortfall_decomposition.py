"""
graphene_fmax_shortfall_decomposition.py

Closes the open item created 2026-10-03: *quantify the g_ds vs. R_g*C_gd
split in the remaining f_max shortfall.*

BACKGROUND.  Chapter 4 Section 4.6 compares this repository's small-signal
estimate against Feijoo et al., Sci. Rep. 6, 35717 (2016), which reports
f_max/f_T = 1.3-1.4 for de-embedded graphene FETs at every gate length
measured.  On 2026-10-03 a frozen cross-module default was found to have
held R_g at 1/40 of its correct value at the literature-scale geometry; the
correction moved this model's f_max/f_T from 0.928 to 0.647 -- FURTHER from
that literature, not closer.  Chapter 4 then stated, without measuring it,
that "the remaining shortfall sits in g_ds and in the R_g*C_gd feedback
term, not in the gate resistance", and left the split open.

The two terms are separable by construction, because the f_max denominator
in rf_small_signal_model._compute_fT_fmax_core is a SUM:

    denom = g_ds*(R_g + R_s)  +  2*pi*f_T*C_gd*R_g
            \_____ TERM A ___/    \______ TERM B ______/

    f_max/f_T = 1 / (2*sqrt(denom))

so zeroing each term in turn is exact, not a perturbation.  This module
measures A and B, splits A further into its g_ds*R_g and g_ds*R_s parts
(the only two places R_g and R_s enter at all), and reports what each
counterfactual does to f_max/f_T.

RESULT, stated up front because it contradicts the sentence it was written
to quantify.  At the literature-scale 40 um / 8-finger geometry and
V_ds = 0.1 V, term B carries **0.22 %** of the denominator and term A
carries **99.78 %**.  Deleting C_gd entirely -- not reducing it, deleting
it -- moves f_max/f_T from 0.6430 to 0.6437, a **0.11 %** improvement
against a shortfall of a factor of 2.2.  The R_g*C_gd feedback term is not
a contributing cause of the shortfall at any bias tested (its share runs
0.1 % at V_ds = 0.05 V to 2.0 % at V_ds = 1 V).  Chapter 4's sentence named
two causes; **one of them is not a cause**, and that is recorded in place
rather than quietly dropped.

THE SECOND, LARGER RESULT.  Term A is not merely dominant, it is nearly an
identity.  Because `transfer_characteristic()` computes

    Id = V_ds / (R_channel(V_g, V_ds) + R_c,total)

the model's drain current is a bias-dependent RESISTOR with no saturation
mechanism of any kind -- no velocity saturation, no pinch-off, no
drain-field cutoff.  Its output conductance is therefore the channel
conductance itself:

    g_ds  =  d(Id)/d(V_ds)  ~  1/R_total        (measured: within 0.3 %
                                                 at V_ds = 0.05 V, 4.9 %
                                                 at V_ds = 1 V)

and term A collapses to a pure resistance ratio

    TERM A  ~  (R_g + R_s) / R_total        =>
    f_max/f_T  ~  (1/2) * sqrt( R_total / (R_g + R_s) )

At W = 40 um the whole device is R_total = 26.38 Ohm while
R_g + R_s = 15.83 Ohm, so the ratio inside the square root is 1.67 and
f_max/f_T is pinned near 0.64.  Reaching 1.41 would need
R_total/(R_g+R_s) >= 4*1.41^2 = 7.95.

**The consequence is that R_g cannot close the gap even in the limit
R_g = 0.**  With a perfect zero-resistance gate, R_s = 7.5 Ohm survives and
the ratio becomes 0.9353 exactly (and 0.938 by the resistance-ratio identity
below) -- still 28 % short of the bottom of the Feijoo band.  This module
checks that limit exactly rather than extrapolating to it.  So the
entire multi-finger R_g thread in this repository, including the 40x
delivery bug that 2026-10-03 spent a session on, was never operating on the
binding constraint.  Finding the bug was correct; believing R_g was the
lever was not.

WHAT THE SHORTFALL ACTUALLY IS, then: the absence of current saturation.
The model's intrinsic voltage gain g_m/g_ds runs 0.005 to 0.096 across
V_ds = 0.05-1 V.  MARKED AS AN ASSERTION, 2026-10-04: an earlier draft of
this paragraph added "i.e. one to two orders of magnitude below the g_m/g_ds
of order unity and above that real GFETs exhibit", and no sourced g_m/g_ds
figure for a comparable device was found to support it.  What IS sourced is
the mechanism, not the number -- the Chalmers GFET high-frequency thesis
states that "f_max is mostly limited by the high drain conductance g_ds due
to the lack of current saturation", and that graphene cannot pinch off
because it has no bandgap, so the carrier type inverts along the channel
instead.  R5's 0.2 threshold below is therefore set from THIS MODEL's own
measured range, which makes R5 a drift detector and not a comparison against
the world.  Citations in
notes/2026-10-04-fmax-was-a-resistance-ratio.md, Section 5.  A device whose g_ds is its
channel conductance has no output resistance to speak of, and f_max is
primarily a measurement of output resistance.  This module reports the
g_ds required for f_max/f_T = 1.41 at the model's own R_g, R_s and C_gd,
and the factor between that and the model's g_ds, so the size of the
missing physics is a number rather than an adjective.

A THIRD OBSERVATION, which is a caution about f_T and not about f_max.
Because g_m = d(V_ds/R)/d(V_g) is almost exactly proportional to V_ds in a
resistor model, peak f_T here scales nearly linearly with drain bias:
10.1 GHz at V_ds = 0.05 V, 20.3 GHz at 0.1 V, 200.4 GHz at 1 V.  Peak f_T
in a real device saturates with V_ds.  So the agreement between this
model's ~20 GHz and the literature's tens of GHz is a fact about the
chosen bias point, not an independent corroboration -- the same model
"agrees" with record 200 GHz devices if asked at V_ds = 1 V.  This is
reported, not repaired; repairing it requires the saturation term above.

EXACTNESS.  Five checks here are tolerance-free.  Per the 2026-10-03
finding that an exactness claim must name the operation it is exact under:
E1 is exact under FLOATING-POINT SUMMATION ORDER (A and B are recombined in
the same order the model adds them, so the reconstruction is bitwise);
E2 is exact under MULTIPLICATION BY ZERO; E3 is exact under the closed-form
1/N^2 law in `gate_resistance` evaluated at powers of two, where N^2 is
representable and the division is by an exact power of 4; E5 is exact under
IDENTITY OF CODE PATH (the default-geometry path of this module calls the
same core as rf_small_signal_model and must agree bitwise, which is the
overlapping-limit validation against the superseded model).  E4 is a
round-trip and is checked to 1 ULP rather than claimed exact, because
sqrt is involved.

MAGNITUDES.  Per the top open item of 2026-10-03 -- a PASS/FAIL at zero is
silent about magnitude -- every MUST_CHANGE control in this module prints
the size of the response it measured, and the two controls whose whole point
is a magnitude (C4, C5) assert a BAND rather than a sign.  This does not
close that item, which asks for a committed repository-wide per-mutation
sensitivity baseline; it is one module written in the form the item asks
for.

Outputs: fmax_shortfall_decomposition.png,
         fmax_shortfall_decomposition_output.txt (transcript, registered in
         graphene_pristine_transcript_audit.py's SUITES in this same
         session per the 2026-10-03 wiring rule).
"""

import numpy as np
import matplotlib.pyplot as plt

import graphene_fet_model as gfet
import rf_small_signal_model as rf

# Literature target this repository compares f_max against: Feijoo et al.,
# Sci. Rep. 6, 35717 (2016), f_max/f_T = 1.3-1.4 de-embedded at every gate
# length measured.  Named once, here, rather than re-typed, per 2026-10-01
# (a citation is a claim) and 2026-09-24 (import, do not re-type).
FEIJOO_RATIO_LOW = 1.3
FEIJOO_RATIO_HIGH = 1.4
FEIJOO_RATIO_REF = 1.41    # the specific 35.4/50 GHz de-embedded pair

VG_SWEEP = np.linspace(-1.5, 3.5, 400)   # the sweep plot_fT_fmax uses
VDS_REF = 0.1                            # the bias Chapter 4 quotes

_PASS, _FAIL = [], []


def check(label, ok, detail=""):
    (_PASS if ok else _FAIL).append(label)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f"  -- {detail}" if detail else ""))
    return ok


class Geometry:
    """A (W, N_fingers) device, with gfet's module-level width-dependent
    names rebound for the duration.  This is the same rebinding
    compute_fT_fmax performs; it is done here explicitly so that R_g and
    R_s are read at CALL time inside the rebound scope, which is the
    2026-10-03 correction, rather than captured at def time, which was the
    2026-10-03 fault."""

    def __init__(self, W=None, N_fingers=1):
        self.W = W
        self.N = N_fingers

    def __enter__(self):
        self._saved = (gfet.W, gfet.Rc_total)
        if self.W is not None:
            gfet.W = self.W
            gfet.Rc_total = 2 * (gfet.Rc_per_width_ohm_um * 1e-6) / self.W
        return self

    def __exit__(self, *a):
        gfet.W, gfet.Rc_total = self._saved
        return False


def decompose(Vg=None, Vds=None, W=Ellipsis, N_fingers=None,
              Rg_override=None, Rs_override=None, Cgd_frac=None,
              gds_scale=1.0):
    """Term-by-term decomposition of the f_max denominator.

    Returns a dict of arrays over Vg plus the scalars R_g, R_s.  The
    optional overrides exist so that each counterfactual is produced by
    the SAME code path as the baseline -- a counterfactual computed by a
    second, parallel expression would be testing that expression.

    EVERY DEFAULT HERE IS LATE-BOUND, and that is not a style choice.
    The first committed version of this function read

        def decompose(Vg=VG_SWEEP, Vds=VDS_REF, W=rf.W_RF,
                      N_fingers=rf.N_FINGERS_RF, ...)

    i.e. it captured two of this module's own module-level constants and
    two of rf_small_signal_model's as function-parameter DEFAULTS, read
    once at def time.  That is precisely the fault this module exists to
    quantify the consequences of, and for rf.W_RF / rf.N_FINGERS_RF it is
    the exact cross-module shape of the 2026-10-03 R_g bug.  It was caught
    by this repository's own standing guard
    (graphene_mutation_arrival_probe.py, Section 7b) on the first run
    after the module was committed -- the guard doing exactly its job, on
    a module written by someone who had just spent a session reading about
    the fault.  Knowing a failure mode is not the same as not committing
    it; a standing guard is.

    W needs a distinguishable sentinel rather than None, because W=None
    is a MEANINGFUL request here (use graphene_fet_model's own width, the
    normalized 1 um device, which is what E5's overlapping-limit check
    needs) and must not be confused with 'argument not supplied'.  The
    sentinel is the BUILTIN Ellipsis rather than a module-level
    `_DEFAULT = object()`, because a module-level sentinel is itself a
    module-level name captured as a default, and the Section 7b guard --
    correctly, by its own stated rule -- reports it.  A sentinel that is
    never rebound is harmless in fact but indistinguishable from a
    harmful capture by AST, and 2026-10-01's lesson is that a
    classifier's name must not drift from what it measures.  Using a
    builtin keeps the guard's rule exactly as written.
    """
    Vg = VG_SWEEP if Vg is None else Vg
    Vds = VDS_REF if Vds is None else Vds
    W = rf.W_RF if W is Ellipsis else W
    N_fingers = rf.N_FINGERS_RF if N_fingers is None else N_fingers
    with Geometry(W, N_fingers):
        gm, Id = rf.transconductance(Vg, Vds=Vds)
        gds = rf.output_conductance(Vg, Vds=Vds) * gds_scale
        Cgs = rf.gate_capacitance(Vg)
        frac = rf.Cgd_over_Cgs if Cgd_frac is None else Cgd_frac
        Cgd = frac * Cgs
        Rg = rf.gate_resistance(N_fingers=N_fingers) if Rg_override is None \
            else Rg_override
        Rs = rf.source_access_resistance() if Rs_override is None \
            else Rs_override
        Rc_total = gfet.Rc_total

    fT = np.abs(gm) / (2 * np.pi * Cgs)
    term_gds_Rg = gds * Rg
    term_gds_Rs = gds * Rs
    term_A = gds * (Rg + Rs)             # summed as the model sums it
    term_B = 2 * np.pi * fT * Cgd * Rg
    denom_raw = term_A + term_B          # same association order as the core
    denom = np.clip(denom_raw, 1e-30, None)
    fmax = fT / (2 * np.sqrt(denom))

    return dict(Vg=Vg, Vds=Vds, gm=gm, Id=Id, gds=gds, Cgs=Cgs, Cgd=Cgd,
                fT=fT, fmax=fmax, ratio=fmax / fT,
                term_A=term_A, term_B=term_B,
                term_gds_Rg=term_gds_Rg, term_gds_Rs=term_gds_Rs,
                denom=denom, denom_raw=denom_raw,
                Rg=Rg, Rs=Rs, Rc_total=Rc_total,
                i_peak=int(np.argmax(fT)))


def denom_required(ratio):
    """The f_max denominator a target f_max/f_T implies.  Inverting
    r = 1/(2*sqrt(d)) gives d = 1/(4 r^2) -- algebra, with no model in it."""
    return 1.0 / (4.0 * ratio ** 2)


def required_gds(d, ratio_target=None):
    """The g_ds that would put this device at `ratio_target`, holding R_g,
    R_s and C_gd fixed.  f_T does not depend on g_ds, so term B is a
    constant of this inversion and the solve is linear:
        g_ds_req = (denom_req - term_B) / (R_g + R_s)
    `ratio_target` is late-bound to FEIJOO_RATIO_REF, per the note in
    decompose().
    Returns (g_ds_req, factor) at the peak-f_T bias, or (nan, nan) when
    term B alone already exceeds the required denominator (i.e. the target
    is unreachable at any g_ds, including zero)."""
    ratio_target = FEIJOO_RATIO_REF if ratio_target is None else ratio_target
    i = d['i_peak']
    need = denom_required(ratio_target)
    num = need - d['term_B'][i]
    if num <= 0:
        return float('nan'), float('nan')
    g_req = num / (d['Rg'] + d['Rs'])
    return g_req, d['gds'][i] / g_req


def report_split(d, title):
    i = d['i_peak']
    print(f"\n{title}")
    print(f"  bias / geometry : V_ds = {d['Vds']:.2f} V, W = {gfet_W_of(d):.0f} um, "
          f"R_g = {d['Rg']:.4f} Ohm, R_s = {d['Rs']:.4f} Ohm, "
          f"R_c,total = {d['Rc_total']:.4f} Ohm")
    print(f"  peak f_T        : {d['fT'][i]/1e9:9.4f} GHz at V_g = {d['Vg'][i]:+.4f} V")
    print(f"  f_max there     : {d['fmax'][i]/1e9:9.4f} GHz   "
          f"f_max/f_T = {d['ratio'][i]:.6f}")
    tot = d['denom'][i]
    for name, val in (("A = g_ds*(R_g+R_s)", d['term_A'][i]),
                      ("    of which g_ds*R_g", d['term_gds_Rg'][i]),
                      ("    of which g_ds*R_s", d['term_gds_Rs'][i]),
                      ("B = 2*pi*f_T*C_gd*R_g", d['term_B'][i])):
        print(f"  {name:<24s} = {val:.6e}   ({100*val/tot:7.4f} % of denominator)")
    print(f"  denominator       = {tot:.6e}")
    print(f"  required for f_max/f_T = {FEIJOO_RATIO_REF}: {denom_required(FEIJOO_RATIO_REF):.6e} "
          f"  (must shrink {tot/denom_required(FEIJOO_RATIO_REF):.4f}x)")


def gfet_W_of(d):
    # the width this decomposition was evaluated at, recovered from
    # R_c,total rather than carried alongside it, so the printed width and
    # the resistance used cannot disagree
    return 2 * gfet.Rc_per_width_ohm_um / d['Rc_total']


# ---------------------------------------------------------------------------
# The identity battery.
#
# These are the checks a mutation can be tested against: each is a function
# of the committed code only, returns a bool, and is recomputed from
# scratch, so re-running the battery under a monkeypatch re-derives every
# number rather than re-reading a cached one.  (2026-09-30: a mutation that
# does not arrive is indistinguishable from a system that does not
# respond.)  Anything that merely PRINTS a measurement -- the bias table,
# the split -- lives outside the battery, because it cannot fail.
# ---------------------------------------------------------------------------

def battery():
    out = {}
    d = decompose()
    i = d['i_peak']

    with Geometry(rf.W_RF, rf.N_FINGERS_RF):
        _, fmax_m, _, gds_m, Cgs_m = rf._compute_fT_fmax_core(
            VG_SWEEP, Vds=VDS_REF, extrinsic=False, N_fingers=rf.N_FINGERS_RF)
    out['E1'] = (np.array_equal(d['fmax'], fmax_m),
                 f"max |diff| = {np.max(np.abs(d['fmax'] - fmax_m)):.3e}")
    out['E1c'] = (np.array_equal(d['gds'], gds_m) and
                  np.array_equal(d['Cgs'], Cgs_m), "")

    d0 = decompose(Rg_override=0.0)
    out['E2'] = (bool(np.all(d0['term_B'] == 0.0) and
                      np.all(d0['term_gds_Rg'] == 0.0)), "")
    out['E2b'] = (np.array_equal(d0['denom_raw'], d0['term_gds_Rs']), "")
    dC = decompose(Cgd_frac=0.0)
    out['E2c'] = (bool(np.all(dC['term_B'] == 0.0)) and
                  np.array_equal(dC['term_A'], d['term_A']), "")
    dZ = decompose(Rg_override=0.0, Rs_override=0.0, Cgd_frac=0.0)
    out['E2d'] = (bool(np.all(dZ['denom_raw'] == 0.0)) and
                  bool(np.all(dZ['denom'] == 1e-30)),
                  "every parasitic deleted: denominator exactly 0.0, clipped "
                  "to 1e-30 by the model's own guard")

    d16 = decompose(N_fingers=16)
    out['E3'] = (np.array_equal(d16['term_gds_Rg'] * 4.0, d['term_gds_Rg']),
                 f"R_g: {d['Rg']:.6f} -> {d16['Rg']:.6f} Ohm")
    out['E3b'] = (np.array_equal(d16['term_B'] * 4.0, d['term_B']), "")
    out['E3c'] = (np.array_equal(d16['term_gds_Rs'], d['term_gds_Rs']), "")

    r = d['ratio'][i]
    ulps = abs(1.0 / (2.0 * np.sqrt(denom_required(r))) - r) / np.spacing(r)
    out['E4'] = (ulps <= 1.0, f"{ulps:.3f} ULP")

    d_def = decompose(W=None, N_fingers=1)
    fT_r, fmax_r, _, _, _ = rf.compute_fT_fmax(VG_SWEEP, Vds=VDS_REF)
    out['E5'] = (np.array_equal(d_def['fmax'], fmax_r), "")
    out['E5b'] = (np.array_equal(d_def['fT'], fT_r), "")

    fA = d['term_A'][i] / d['denom'][i]
    fB = d['term_B'][i] / d['denom'][i]
    out['R1'] = (fB < 0.01, f"term B share = {100*fB:.4f} %")
    out['R2'] = (fA > 0.99, f"term A share = {100*fA:.4f} %")

    rows = []
    for Vds in (0.05, 0.1, 0.3, 0.6, 1.0):
        db = decompose(Vds=Vds)
        j = db['i_peak']
        Rtot = Vds / db['Id'][j]
        rows.append((Vds, db['fT'][j] / 1e9, db['ratio'][j],
                     db['term_A'][j] / db['denom'][j],
                     db['term_B'][j] / db['denom'][j],
                     abs(db['gm'][j]) / db['gds'][j], db['gds'][j] * Rtot))
    out['R3'] = (max(x[4] for x in rows) < 0.05,
                 f"largest B share = {100*max(x[4] for x in rows):.4f} % at "
                 f"V_ds = {max(rows, key=lambda x: x[4])[0]} V")
    out['R4'] = (all(abs(x[6] - 1.0) < 0.05 for x in rows),
                 f"g_ds*R_total in [{min(x[6] for x in rows):.4f}, "
                 f"{max(x[6] for x in rows):.4f}]")
    out['R5'] = (all(x[5] < 0.2 for x in rows),
                 f"g_m/g_ds in [{min(x[5] for x in rows):.5f}, "
                 f"{max(x[5] for x in rows):.5f}]")

    base = d['ratio'][i]
    rC = dC['ratio'][i]
    rG = d0['ratio'][i]
    out['C4'] = (0.0 < (rC - base) / base < 0.005,
                 f"response = {100*(rC-base)/base:+.5f} % -- the C_gd knob is "
                 f"ALIVE but its authority is negligible")
    out['C5'] = (rG < FEIJOO_RATIO_LOW,
                 f"R_g = 0 limit gives {rG:.6f}, "
                 f"{100*(FEIJOO_RATIO_LOW-rG)/FEIJOO_RATIO_LOW:.2f} % below the "
                 f"bottom of the Feijoo band")
    out['C6'] = ((rG - base) / (rC - base) > 100.0,
                 f"R_g response / C_gd response = {(rG-base)/(rC-base):.2f}x")

    Rtot = VDS_REF / d['Id'][i]
    approx = 0.5 * np.sqrt(Rtot / (d['Rg'] + d['Rs']))
    out['R6'] = (abs(approx - base) / base < 0.01,
                 f"(1/2)sqrt(R_total/(R_g+R_s)) = {approx:.6f} vs exact "
                 f"{base:.6f}, rel. err {abs(approx-base)/base:.4e}")
    g_req, g_fac = required_gds(d)
    out['R7'] = (3.0 < g_fac < 10.0, f"required g_ds factor = {g_fac:.4f}")

    return out, d, rows


# ---------------------------------------------------------------------------
# Validations
# ---------------------------------------------------------------------------

_LABELS = {
    'E1':  "reconstructed f_max is bitwise the model's f_max (400 pts)",
    'E1c': "g_ds and C_gs bitwise identical to the core's",
    'E2':  "R_g = 0 zeroes term B and g_ds*R_g EXACTLY",
    'E2b': "R_g = 0 leaves the denominator exactly equal to g_ds*R_s",
    'E2c': "C_gd = 0 zeroes term B exactly, term A bitwise unchanged",
    'E2d': "all parasitics deleted: denominator exactly 0.0, then clipped",
    'E3':  "fingers 8->16 divides g_ds*R_g by exactly 4 (1/N^2 law)",
    'E3b': "fingers 8->16 divides term B by exactly 4",
    'E3c': "fingers 8->16 leaves g_ds*R_s exactly invariant",
    'E4':  "ratio -> required denominator -> ratio closes within 1 ULP",
    'E5':  "default geometry: f_max bitwise equals rf.compute_fT_fmax",
    'E5b': "default geometry: f_T bitwise equals rf.compute_fT_fmax",
    'R1':  "term B (R_g*C_gd feedback) is below 1 % of the denominator",
    'R2':  "term A (g_ds) is above 99 % of the denominator",
    'R3':  "term B stays below 5 % of the denominator at EVERY bias",
    'R4':  "g_ds is the channel conductance to within 5 % at every bias",
    'R5':  "g_m/g_ds below 0.2 at every bias (DRIFT detector: the 0.2 is "
           "this model's own range, not literature)",
    'C4':  "MAGNITUDE: deleting C_gd moves f_max/f_T by under 0.5 % (band)",
    'C5':  "MAGNITUDE: R_g = 0 entirely STILL falls short of f_max/f_T = 1.3",
    'C6':  "CONTROL: the C_gd and R_g responses are not aliased (>100x apart)",
    'R6':  "the resistance-ratio identity reproduces f_max/f_T within 1 %",
    'R7':  "the required g_ds reduction is a factor between 3 and 10",
}


def validate():
    print("=" * 78)
    print("SECTION 1.  The identity battery")
    print("=" * 78)
    results, d, rows = battery()
    for k, (ok, detail) in results.items():
        check(f"{k:<4s} {_LABELS[k]}", ok, detail)

    print()
    print("=" * 78)
    print("SECTION 2.  The split, measured")
    print("=" * 78)
    report_split(d, f"Literature-scale device, V_ds = {VDS_REF} V "
                    f"({rf.W_RF*1e6:.0f} um / {rf.N_FINGERS_RF} fingers)")
    print()
    print("  Bias dependence of the split (peak-f_T bias at each V_ds):")
    print(f"  {'V_ds':>6} {'peak f_T':>13} {'f_max/f_T':>10} {'A share':>9} "
          f"{'B share':>9} {'g_m/g_ds':>9} {'g_ds*R_total':>13}")
    for x in rows:
        print(f"  {x[0]:6.2f} {x[1]:9.3f} GHz {x[2]:10.4f} "
              f"{100*x[3]:8.4f}% {100*x[4]:8.4f}% {x[5]:9.5f} {x[6]:13.6f}")
    print("  NOTE: peak f_T scales nearly LINEARLY with V_ds (x19.8 over a x20")
    print("  bias range), which a saturating device would not do.  The ~20 GHz")
    print("  this chapter quotes is therefore a fact about the chosen bias.")

    print()
    print("=" * 78)
    print("SECTION 3.  Counterfactuals: can either named cause close the gap?")
    print("=" * 78)
    i = d['i_peak']
    base = d['ratio'][i]
    print(f"  baseline f_max/f_T                                   = {base:.6f}")
    for name, dc in [
            ("delete C_gd entirely (term B -> 0)", decompose(Cgd_frac=0.0)),
            ("16 fingers instead of 8 (R_g / 4)", decompose(N_fingers=16)),
            ("64 fingers instead of 8 (R_g / 64)", decompose(N_fingers=64)),
            ("delete R_g entirely (perfect gate)", decompose(Rg_override=0.0)),
            ("delete R_s entirely (perfect contacts)",
             decompose(Rs_override=0.0)),
            ]:
        rr = dc['ratio'][i]
        print(f"  {name:<50s} = {rr:.6f}   "
              f"({100*(rr-base)/base:+9.4f} % vs baseline)")
    # R_g = R_s = C_gd = 0 is reported as what it is -- the model's 1e-30
    # denominator clip -- rather than as a device.  Printing the resulting
    # 5e14 as a counterfactual f_max/f_T would be quoting a guard constant.
    dZ = decompose(Rg_override=0.0, Rs_override=0.0, Cgd_frac=0.0)
    print(f"  {'delete R_g AND R_s AND C_gd (not a device)':<50s} = "
          f"denominator exactly {dZ['denom_raw'][i]:.1f}, clipped to "
          f"{dZ['denom'][i]:.0e} by the model's guard")
    g_req, g_fac = required_gds(d)
    dS = decompose(gds_scale=1.0 / g_fac)
    print(f"  {'divide g_ds by %.4f (saturation added)' % g_fac:<50s} = "
          f"{dS['ratio'][i]:.6f}   ({100*(dS['ratio'][i]-base)/base:+9.4f} % "
          f"vs baseline)")
    print()
    print(f"  g_ds required for f_max/f_T = {FEIJOO_RATIO_REF}: {g_req:.6e} S "
          f"({g_req/(rf.W_RF*1e6)*1e3:.4f} mS/um) vs model's "
          f"{d['gds'][i]:.6e} S ({d['gds'][i]/(rf.W_RF*1e6)*1e3:.4f} mS/um)")
    print(f"  => the missing physics (current saturation) is worth a factor of "
          f"{g_fac:.4f} in g_ds")
    print()
    print("  VERDICT on the open item.  Chapter 4 named two residual causes,")
    print("  g_ds and the R_g*C_gd feedback term.  The split is 99.78 / 0.22.")
    print("  The feedback term is not a contributing cause at any bias tested,")
    print("  and R_g is not the binding constraint even in the R_g = 0 limit")
    print("  (0.935 < 1.3).  The shortfall is the ABSENCE OF CURRENT")
    print("  SATURATION in transfer_characteristic(), worth ~4.85x in g_ds.")
    return results, d, rows


# ---------------------------------------------------------------------------
# SECTION 4.  Mutation harness
#
# Every mutation is applied to a COPY of the live state (a monkeypatch that
# is reverted), the full battery is re-run, and the mutation is required to
# flip at least one check.  Control A is the unmutated battery (must be all
# green).  Control B is a mutation that SHOULD NOT be caught, so that
# "every mutation is caught" is a claim about the mutations and not a
# property of a battery that fails on anything.
# ---------------------------------------------------------------------------

class patch:
    """Temporarily rebind attributes on a module."""

    def __init__(self, mod, **kw):
        self.mod, self.kw = mod, kw

    def __enter__(self):
        self.old = {k: getattr(self.mod, k) for k in self.kw}
        for k, v in self.kw.items():
            setattr(self.mod, k, v)

    def __exit__(self, *a):
        for k, v in self.old.items():
            setattr(self.mod, k, v)
        return False


def _battery_failures():
    res, _, _ = battery()
    return sorted(k for k, (ok, _) in res.items() if not ok)


def mutation_report():
    print()
    print("=" * 78)
    print("SECTION 4.  Mutation harness")
    print("=" * 78)

    ctrlA = _battery_failures()
    check("MA  CONTROL A: the unmutated battery is fully green",
          ctrlA == [], f"failures = {ctrlA}")

    mutations = []

    # M1: break the summation ORDER that E1's exactness premise names.
    # Same algebra, different association: gds*Rg + gds*Rs + B instead of
    # gds*(Rg+Rs) + B.  This must break E1 and nothing about the physics.
    def _assoc(Vg=None, **kw):
        d = _orig_decompose(Vg=Vg, **kw)
        d['denom_raw'] = d['term_gds_Rg'] + d['term_gds_Rs'] + d['term_B']
        d['denom'] = np.clip(d['denom_raw'], 1e-30, None)
        d['fmax'] = d['fT'] / (2 * np.sqrt(d['denom']))
        d['ratio'] = d['fmax'] / d['fT']
        return d
    mutations.append(("M1 re-associate the denominator sum (same algebra)",
                      lambda: patch(_MOD, decompose=_assoc)))

    # M2: the pre-2026-08-26 model -- R_s missing from the f_max denominator.
    mutations.append(("M2 drop R_s from the denominator (pre-08-26 model)",
                      lambda: patch(rf, source_access_resistance=lambda: 0.0)))

    # M3: 1/N instead of 1/N^2 for the multi-finger gate.
    def _bad_fingers(L=None, W=None, N_fingers=1):
        L = gfet.L if L is None else L
        W = gfet.W if W is None else W
        return rf.R_sheet_gate * W / (3.0 * L * N_fingers)
    mutations.append(("M3 1/N instead of 1/N^2 finger law",
                      lambda: patch(rf, gate_resistance=_bad_fingers)))

    # M4: reinstate the 2026-10-03 frozen cross-module default -- R_g read
    # at W = 1 um regardless of the rebind, i.e. 40x too small.
    def _frozen(L=None, W=None, N_fingers=1):
        return rf.R_sheet_gate * 1e-6 / (3.0 * 200e-9 * N_fingers ** 2)
    mutations.append(("M4 reinstate the 2026-10-03 frozen R_g default (40x small)",
                      lambda: patch(rf, gate_resistance=_frozen)))

    # M5: give the DC model a fake saturation by dividing g_ds by 5 --
    # the physics this chapter says is missing.  R4 and R5 must notice.
    _orig_gds = rf.output_conductance
    mutations.append(("M5 divide g_ds by 5 (fake saturation)",
                      lambda: patch(rf, output_conductance=(
                          lambda *a, **k: _orig_gds(*a, **k) / 5.0))))

    # M6: C_gd fraction 0.4 -> 40, a 100x feedback capacitance.
    mutations.append(("M6 C_gd/C_gs 0.4 -> 40 (100x feedback capacitance)",
                      lambda: patch(rf, Cgd_over_Cgs=40.0)))

    caught = 0
    for name, mk in mutations:
        with mk():
            fails = _battery_failures()
        ok = len(fails) > 0
        caught += ok
        check(f"    {name}", ok,
              f"killed by {len(fails)} check(s): {','.join(fails) if fails else 'NONE'}")

    # Control B: a mutation that must NOT be caught.  Rebinding R_sheet_gate
    # to its own value is an identity; if the battery fails on this, it is
    # failing on the act of patching rather than on the patch.
    with patch(rf, R_sheet_gate=rf.R_sheet_gate):
        ctrlB = _battery_failures()
    check("MB  CONTROL B: an identity rebind of R_sheet_gate is NOT caught",
          ctrlB == [], f"failures = {ctrlB}")

    print(f"\n  Mutants killed: {caught} of {len(mutations)}")
    if caught < len(mutations):
        print("  NOTE: a surviving mutant is reported, not hidden.")
    print("  WORTH READING OFF: M4 -- the 2026-10-03 frozen-default bug, 40x")
    print("  in R_g -- is killed by exactly ONE check, R7, and R7 is a")
    print("  MAGNITUDE check (how big a g_ds reduction the gap implies), not a")
    print("  sign check.  Every PASS/FAIL-at-zero check in this battery is")
    print("  blind to it, because shrinking R_g moves this module's conclusion")
    print("  in the direction it already argues.  That is the 2026-10-03 top")
    print("  open item reproduced from the other side: the only detector that")
    print("  saw the historical bug is the one that retains a magnitude.")
    return caught, len(mutations)


_MOD = None
_orig_decompose = None


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def plot_decomposition(d, rows):
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    Vg, tot = d['Vg'], d['denom']
    axes[0].stackplot(
        Vg,
        100 * d['term_gds_Rs'] / tot,
        100 * d['term_gds_Rg'] / tot,
        100 * d['term_B'] / tot,
        labels=[r'$g_{ds}R_s$', r'$g_{ds}R_g$', r'$2\pi f_T C_{gd} R_g$'],
        colors=['#4C72B0', '#DD8452', '#C44E52'], alpha=0.9)
    axes[0].axvline(Vg[d['i_peak']], color='k', linestyle='--', lw=1.2,
                    label=r'peak $f_T$ bias')
    axes[0].set_xlabel(r'$V_g$ (V)')
    axes[0].set_ylabel(r'share of the $f_{max}$ denominator (%)')
    axes[0].set_ylim(0, 100)
    axes[0].set_xlim(Vg[0], Vg[-1])
    axes[0].set_title(f'Where the $f_{{max}}$ denominator comes from\n'
                      f'({rf.W_RF*1e6:.0f} $\\mu$m / {rf.N_FINGERS_RF} fingers, '
                      f'$V_{{ds}}$ = {VDS_REF} V)', fontsize=10)
    axes[0].legend(loc='center right', fontsize=9)
    axes[0].grid(True, alpha=0.25)

    i = d['i_peak']
    g_fac = required_gds(d)[1]
    bars = [('model\nas committed', decompose(), '#4C72B0'),
            ('$C_{gd}=0$', decompose(Cgd_frac=0.0), '#C44E52'),
            ('64 fingers\n($R_g$/64)', decompose(N_fingers=64), '#DD8452'),
            ('$R_g=0$\n(perfect gate)', decompose(Rg_override=0.0), '#DD8452'),
            ('$R_g=R_s=0$\n(no access $R$)',
             decompose(Rg_override=0.0, Rs_override=0.0), '#8172B3'),
            ('$g_{ds}/%.2f$\n(saturation)' % g_fac,
             decompose(gds_scale=1.0 / g_fac), '#55A868')]
    vals = [b[1]['ratio'][i] for b in bars]
    axes[1].bar(range(len(vals)), vals, color=[b[2] for b in bars])
    axes[1].axhspan(FEIJOO_RATIO_LOW, FEIJOO_RATIO_HIGH, color='green',
                    alpha=0.15,
                    label='Feijoo et al. 2016\n(de-embedded, 1.3-1.4)')
    axes[1].set_xticks(range(len(vals)))
    axes[1].set_xticklabels([b[0] for b in bars], fontsize=8)
    axes[1].set_ylabel(r'$f_{max}/f_T$')
    axes[1].set_title('No parasitic deletion reaches the literature band;\n'
                      'only $g_{ds}$ does (one input deleted per bar)',
                      fontsize=10)
    for k, v in enumerate(vals):
        axes[1].text(k, v + 0.03, f'{v:.3f}', ha='center', fontsize=8)
    axes[1].legend(fontsize=8, loc='upper left')
    axes[1].grid(True, alpha=0.25, axis='y')

    ax = axes[2]
    for Vds, style in ((0.1, '-'), (1.0, '--')):
        db = decompose(Vds=Vds)
        Rtot = Vds / db['Id']
        ax.plot(db['Vg'], db['gds'] * 1e3, style, color='C3', lw=2.2,
                label=f'$g_{{ds}}$, $V_{{ds}}$ = {Vds} V')
        ax.plot(db['Vg'], 1e3 / Rtot, style, color='C0', lw=1.2,
                label=f'$1/R_{{total}}$, $V_{{ds}}$ = {Vds} V')
    ax.axvline(gfet.V_dirac, color='gray', linestyle=':', label='Dirac point')
    ax.set_xlabel(r'$V_g$ (V)')
    ax.set_ylabel('conductance (mS)')
    ax.set_title('Why $g_{ds}$ is large: the DC model is a resistor\n'
                 r'$g_{ds}\approx 1/R_{total}$ (no saturation mechanism)',
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)

    plt.tight_layout()
    plt.savefig('fmax_shortfall_decomposition.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    import sys
    _MOD = sys.modules[__name__]
    _orig_decompose = decompose
    print("f_max shortfall decomposition: g_ds vs. R_g*C_gd")
    print("(closes the open item created 2026-10-03; see the module docstring)")
    print()
    results, d, rows = validate()
    killed, total = mutation_report()
    print()
    print("=" * 78)
    print(f"SUMMARY: {len(_PASS)} passed, {len(_FAIL)} failed; "
          f"mutants killed {killed}/{total}")
    for f in _FAIL:
        print(f"  FAILED: {f}")
    print("=" * 78)
    plot_decomposition(d, rows)
    print("\nSaved: fmax_shortfall_decomposition.png")
