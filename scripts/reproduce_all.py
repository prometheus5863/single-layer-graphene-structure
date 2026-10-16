#!/usr/bin/env python3
"""Regenerate every figure and transcript, and optionally check them against
the Checkpoint 0 baseline.

    python scripts/reproduce_all.py                  # regenerate everything into results/
    python scripts/reproduce_all.py --check          # ... then compare with docs/baseline_hashes.txt
    python scripts/reproduce_all.py --only graphene_fet_model.py rf_small_signal_model.py

Each module in src/ is run once, in sorted order, in its own process, with
src/ as the working directory (the audit layer enumerates and copies its
sibling modules from there -- docs/INVENTORY.md section 6). Modules that write
their own outputs do so through graphene.paths.results_path; for the modules
whose transcript is their stdout, that stdout is written to the transcript
path below.

Byte-identity with the baseline proves the reorganisation changed nothing. It
does not prove the numbers are right (note 2026-10-02).
"""
import argparse, hashlib, json, os, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
BASELINE = os.path.join(ROOT, "docs", "baseline_hashes.txt")

# module -> transcript it produces on stdout (the files once written with `| tee`)
STDOUT_TRANSCRIPT = {
    "graphene_band_structure_audit.py": "results/transcripts/ch2/band_structure_audit_output.txt",
    "graphene_contact_consistent_argument.py": "results/transcripts/ch4/contact_consistent_argument_output.txt",
    "graphene_covariance_probe_audit.py": "results/transcripts/audits/default_scale_probe_audit_output.txt",
    "graphene_cross_module_delivery_audit.py": "results/transcripts/ch4/cross_module_delivery_audit_output.txt",
    "graphene_dead_name_sweep.py": "results/transcripts/ch4/dead_name_sweep_output.txt",
    "graphene_default_scale_audit.py": "results/transcripts/audits/default_scale_audit_output.txt",
    "graphene_differential_crossover_model.py": "results/transcripts/ch6/differential_crossover_output.txt",
    "graphene_double_regularisation.py": "results/transcripts/ch4/double_regularisation_output.txt",
    "graphene_exact_charge_rederivation.py": "results/transcripts/ch4/exact_charge_rederivation_output.txt",
    "graphene_fet_signed_carrier_model.py": "results/transcripts/ch4/fet_signed_carrier_output.txt",
    "graphene_figure_provenance_audit.py": "results/transcripts/audits/figure_provenance_audit_output.txt",
    "graphene_fmax_shortfall_decomposition.py": "results/transcripts/ch4/fmax_shortfall_decomposition_output.txt",
    "graphene_gds_quadrature_audit.py": "results/transcripts/ch4/gds_quadrature_audit_output.txt",
    "graphene_log_sensitivity_step_audit.py": "results/transcripts/audits/log_sensitivity_step_audit_output.txt",
    "graphene_mutation_arrival_probe.py": "results/transcripts/audits/mutation_arrival_probe_output.txt",
    "graphene_potential_census_audit.py": "results/transcripts/ch4/potential_census_output.txt",
    "graphene_potential_declaration.py": "results/transcripts/ch4/potential_declaration_output.txt",
    "graphene_pristine_transcript_audit.py": "results/transcripts/audits/pristine_transcript_audit_output.txt",
    "graphene_quantum_capacitance_argument.py": "results/transcripts/ch4/quantum_capacitance_argument_output.txt",
    "graphene_quantum_capacitance_overflow.py": "results/transcripts/ch4/quantum_capacitance_overflow_output.txt",
    "graphene_rootfinder_audit.py": "results/transcripts/ch6/rootfinder_audit_output.txt",
    "graphene_second_anchor_bound.py": "results/transcripts/ch7/second_anchor_bound_output.txt",
    "graphene_transport_optical_audit.py": "results/transcripts/ch3/transport_optical_audit_output.txt",
    "rf_small_signal_model.py": "results/transcripts/ch4/rf_small_signal_output.txt",
}
# Superseded demo scripts. They stay in src/ only because several audits count
# every *.py in that directory and print the count (docs/INVENTORY.md section 6);
# they are not part of the results and are never run here.
LEGACY = {"generate_plots.py", "simple_graphene_plots.py", "summary.py",
          "initialize_repo.py", "main.py"}

# Committed transcripts that are NOT regenerated, because they are records,
# not outputs (docs/INVENTORY.md section 3).
NOT_REGENERATED = {"audit_output.txt", "sensitivity_converged_output.txt"}

# Modules whose stdout is not deterministic by construction (they print the
# current commit, working-tree state, or wall-clock timings). Their files are
# still compared; their stdout is reported, not failed.
NONDETERMINISTIC_STDOUT = {"graphene_pristine_transcript_audit.py"}


def sha(b):
    return hashlib.sha256(b).hexdigest()


def pixel_sha(path):
    import numpy as np
    from PIL import Image
    a = np.asarray(Image.open(path).convert("RGBA"))
    return sha(a.tobytes() + str(a.shape).encode())


def snapshot():
    s = {}
    for d, dirs, files in os.walk(ROOT):
        dirs[:] = [x for x in dirs if x not in (".git", "__pycache__")]
        for f in files:
            p = os.path.join(d, f)
            st = os.stat(p)
            s[os.path.relpath(p, ROOT)] = (st.st_mtime_ns, st.st_size)
    return s


def run(modules):
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONHASHSEED="0", SOURCE_DATE_EPOCH="0")
    report = []
    for m in modules:
        before = snapshot()
        t0 = time.time()
        p = subprocess.run([sys.executable, m], cwd=SRC, env=env, capture_output=True, text=True)
        dt = time.time() - t0
        if m in STDOUT_TRANSCRIPT:
            out = os.path.join(ROOT, STDOUT_TRANSCRIPT[m])
            os.makedirs(os.path.dirname(out), exist_ok=True)
            with open(out, "w") as fh:
                fh.write(p.stdout)
        after = snapshot()
        touched = sorted(k for k, v in after.items() if before.get(k) != v)
        report.append(dict(module=m, rc=p.returncode, seconds=round(dt, 1),
                           stdout_sha256=sha(p.stdout.encode()), touched=touched))
        print(f"{m:58s} exit {p.returncode}  {dt:7.1f}s  {len(touched)} file(s)", flush=True)
        if p.returncode not in (0, 1):
            print(p.stderr[-2000:], file=sys.stderr)
    return report


def load_baseline():
    rows = []
    for line in open(BASELINE):
        if line.startswith("#") or not line.strip():
            continue
        h, px, kind, name, producer = line.split()
        rows.append((h, px, kind, name, producer))
    return rows


def check(report):
    base = load_baseline()
    by_module = {r["module"]: r for r in report}
    # where each output now lives: basename -> new relative path
    produced = {}
    for r in report:
        for t in r["touched"]:
            if "__pycache__" in t:
                continue
            produced.setdefault(os.path.basename(t), t)
    ok, bad, skipped, notes = 0, [], [], []
    for h, px, kind, name, producer in base:
        if producer not in by_module:
            skipped.append(f"{name} ({producer} is legacy, not regenerated)")
            continue
        if h == "(overwritten)":
            continue
        if kind == "stdout":
            got = by_module[producer]["stdout_sha256"]
            if got == h:
                ok += 1
            elif producer in NONDETERMINISTIC_STDOUT:
                notes.append(f"{name}: differs (non-deterministic by construction)")
            else:
                bad.append(f"{name}: stdout differs")
            continue
        new = produced.get(os.path.basename(name))
        if new is None:
            bad.append(f"{name}: not produced")
            continue
        p = os.path.join(ROOT, new)
        if sha(open(p, "rb").read()) == h:
            ok += 1
        elif name.endswith(".png") and px != "-" and pixel_sha(p) == px:
            ok += 1
            notes.append(f"{name}: bytes differ, pixels identical (metadata only)")
        else:
            bad.append(f"{name} -> {new}: content differs")
    print(f"\nCHECK against {os.path.relpath(BASELINE, ROOT)}: {ok} identical, "
          f"{len(bad)} mismatched, {len(skipped)} legacy rows skipped")
    for n in notes:
        print("  note:", n)
    for b in bad:
        print("  MISMATCH:", b)
    return not bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--report", help="write a JSON run report here")
    a = ap.parse_args()
    modules = sorted(f for f in os.listdir(SRC) if f.endswith(".py") and f not in LEGACY)
    if a.only:
        modules = [m for m in modules if m in a.only]
    report = run(modules)
    if a.report:
        json.dump(report, open(a.report, "w"), indent=1)
    if a.check:
        sys.exit(0 if check(report) else 1)


if __name__ == "__main__":
    main()
