"""
bracketed_root.py -- the repo's single guarded bisection.

WHY THIS FILE EXISTS
--------------------
On 2026-09-22 a bisection in `graphene_crossover_sensitivity_model.py` ran
with no check that its interval was bracketed.  When no sign change existed
it walked one bound onto the other and returned that bound as a root -- 21
times, every one of the 21 physically plausible, while all five of that
module's exact validations passed and printed underneath.  The fault lived
in the analysis layer, which is why model-level validation could not see it.

On 2026-09-23 `graphene_differential_crossover_model.py` answered that with
a `bracketed_bisect` that raises on an unbracketed interval, and the answer
was good.  But it answered it *locally*: the function sat in the module that
needed it, so the 2026-09-23 log left an open item to audit the repo's other
numeric solvers for the same class.  That audit ran on 2026-09-24
(`graphene_rootfinder_audit.py`) and found one further unguarded bisection --
`p2_threshold()` in `graphene_per_metal_crossover_model.py` -- which meant the
guard needed one home rather than one copy per module.  This is that home.

A SECOND-ORDER FINDING FROM THAT AUDIT, recorded here because it changes what
the guard is protecting against
--------------------------------------------------------------------------
The 2026-09-22 wording -- an unguarded bisection "returns the ENDPOINT" -- is
not quite true, and the way it fails is asymmetric:

  * when no root exists because f exceeds the target everywhere, the loop
    assigns `lo = mid` and lands on `hi` **exactly**;
  * when no root exists because f is below the target everywhere, the loop
    assigns `hi = mid` and stalls **one ulp above `lo`** -- never reaching it,
    because once `hi` is the next double after `lo`, `0.5 * (lo + hi)` rounds
    back up to `hi`.

Both were pre-registered and measured on 2026-09-24 (Validation 3 of the
audit).  The practical consequence: a caller cannot detect the fault by
testing `result == lo or result == hi`.  Half the time that test is false and
the result is still meaningless.  The guard has to be at the entrance.
"""

import numpy as np


def bracketed_bisect(f, lo, hi, tol=1e-14, max_iter=200, rtol=0.0):
    """
    Bisection that REFUSES an unbracketed interval.

    Raises ValueError when `f(lo)` and `f(hi)` share a sign, rather than
    returning a sentinel: a sentinel is something a caller can forget to
    check and an exception is not.

    The guard is one line and it is the first line.

    `tol` IS AN ABSOLUTE INTERVAL WIDTH AND IT IS NOT SCALE-FREE
    -----------------------------------------------------------
    This is the trap that sharing the function exposed, and it is a
    different fault class from the one the guard prevents.  The default
    1e-14 was chosen in `graphene_differential_crossover_model` for a
    variable of order 1 eV, where it sits near machine precision: 48
    halvings, relative precision 1e-14.  The first caller outside that
    module bisected a decay length in METRES, order 5e-11, and 1e-14 is
    then a *coarse* tolerance -- it stops after 16 halvings at a relative
    precision of 2e-4.  Measured on 2026-09-24: the result landed 2.5e11
    ulps from the converged root, with residual 3.4e-6 where the
    pre-existing unguarded loop had 1e-16.  Nothing raised; the number was
    simply wrong in its fifth significant figure.

    So pass `rtol` (relative interval width) for any variable whose scale is
    not order 1, or `tol=0.0` to force the full `max_iter` halvings.  `rtol`
    defaults to 0.0 so that every pre-existing caller is bitwise unaffected.

    The general lesson, recorded because it cost a wrong number to find:
    **a guard moved to a new call site brings its own defaults with it, and
    a default is a claim about scale.**
    """
    flo, fhi = f(lo), f(hi)
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if np.sign(flo) == np.sign(fhi):
        raise ValueError(
            f"root not bracketed on [{lo:.6g}, {hi:.6g}]: "
            f"f(lo)={flo:+.6e}, f(hi)={fhi:+.6e} -- refusing to bisect")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        width = max(tol, rtol * max(abs(lo), abs(hi)))
        if fmid == 0.0 or (hi - lo) < width:
            return mid
        if np.sign(fmid) == np.sign(flo):
            lo, flo = mid, fmid
        else:
            hi, fhi = mid, fmid
    return 0.5 * (lo + hi)
