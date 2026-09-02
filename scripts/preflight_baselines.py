"""Offline preflight for the multi-baseline comparison run (#83).

Fail fast BEFORE spending time (and cloud money) on the metal-aware baseline rescores.
This checks every frozen input the run depends on, hashes the two things that must not
have drifted (the pre-registration + protocol.yaml), reports exactly which run steps are
still outstanding and which tools/receptors are missing, and prints the run plan.

Stdlib + openafr.protocol only (NO rdkit, NO docking binaries), so it runs on ANY host
-- including the osx-arm64 machine that cannot run the scorers -- as the go/no-go gate a
cloud operator runs first. It never scores anything; it only inspects the tree.

Exit codes:
  0  READY-TO-GRADE   -- poses + both baseline output dirs present; run the grader.
  2  RUN-STEPS-REMAIN -- something the run needs is missing (the normal state pre-run).
  1  DRIFT            -- a frozen hash (pre-registration or protocol) does NOT match;
                         the comparison is not credible until this is resolved.

Usage:
    python scripts/preflight_baselines.py [POSE_DIR]
        POSE_DIR   dir of Vina .pdbqt poses to rescore (default work/screen2)
"""
import hashlib
import os
import pathlib
import shutil
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from openafr import baselines, protocol  # noqa: E402

PREREG = "work/PREREGISTRATION_baselines_multi.md"
PREREG_SHA_FILE = "work/PREREG_baselines_multi.sha256"
GRADE_RECEPTOR = "work/receptor_A.pdb"          # heme Fe source for mode C + grading
# Per-baseline receptor a cloud operator must prepare (absent from git; see the RUNBOOK).
BASELINE_RECEPTOR = {
    "gnina": "work/receptor.pdbqt",             # prep_receptor.py
    "plants": "work/receptor_plants.mol2",      # SPORES --mode complete
}
# Tools each run step shells out to (checked for presence only, never invoked here).
BASELINE_TOOLS = {"gnina": ["gnina"], "plants": ["PLANTS", "SPORES", "obabel"]}

OK, MISSING, BAD = "  ok ", "MISS ", "DRIFT"


def _sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def _check_frozen_sha(sha_file):
    """Verify every 'sha  path' line in a .sha256 file. Returns (all_ok, [lines])."""
    lines, all_ok = [], True
    for raw in pathlib.Path(ROOT / sha_file).read_text().splitlines():
        raw = raw.strip()
        if not raw:
            continue
        want, _, rel = raw.partition("  ")
        target = ROOT / rel.strip()
        if not target.exists():
            lines.append(f"  [{MISSING}] {rel} (pinned file absent)")
            all_ok = False
            continue
        got = _sha256(target)
        ok = got == want
        all_ok = all_ok and ok
        lines.append(f"  [{OK if ok else BAD}] {rel}  {got[:16]}...")
    return all_ok, lines


def _pose_count(pose_dir):
    d = ROOT / pose_dir
    if not d.is_dir():
        return None
    return sum(1 for f in os.listdir(d) if f.endswith(".pdbqt"))


def main():
    pose_dir = sys.argv[1] if len(sys.argv) > 1 else "work/screen2"
    print("=" * 70)
    print("PREFLIGHT — multi-baseline comparison run (#83)")
    print("=" * 70)

    drift = False

    # 1. Frozen hashes: the rules + the protocol must not have moved.
    print("\n[1] Frozen-hash checks (must all be ok):")
    prereg_ok, prereg_lines = _check_frozen_sha(PREREG_SHA_FILE)
    for ln in prereg_lines:
        print(ln)
    proto_path = ROOT / "protocol.yaml"
    proto_got = _sha256(proto_path)
    proto_ok = proto_got == protocol.PROTOCOL_SHA
    print(f"  [{OK if proto_ok else BAD}] protocol.yaml  {proto_got[:16]}...")
    drift = drift or not prereg_ok or not proto_ok

    # 2. Poses: the identical run-2 poses every scorer rescores.
    print(f"\n[2] Held-out poses ({pose_dir}):")
    n = _pose_count(pose_dir)
    if n is None:
        print(f"  [{MISSING}] absent — regenerate first (see RUNBOOK step 1: run_screen2.sh)")
    elif n == 0:
        print(f"  [{MISSING}] present but empty — docking did not populate it")
    else:
        print(f"  [{OK}] {n} pose files")

    # 3. The grading receptor (mode C reads the heme Fe from it) — committed.
    print("\n[3] Grading receptor (committed):")
    gr = ROOT / GRADE_RECEPTOR
    print(f"  [{OK if gr.exists() else MISSING}] {GRADE_RECEPTOR}")

    # 4. Per-baseline readiness: run step outstanding vs already-scored, tools, receptor.
    print("\n[4] Baselines (each is a cloud/Linux run step):")
    outstanding = []
    for spec in baselines.REGISTRY:
        out_dir = ROOT / spec.subdir
        scored = out_dir.is_dir() and any(
            f.endswith(spec.ext) for f in os.listdir(out_dir))
        recep = BASELINE_RECEPTOR.get(spec.key)
        recep_ok = recep is None or (ROOT / recep).exists()
        tools = BASELINE_TOOLS.get(spec.key, [])
        have = {t: shutil.which(t) is not None for t in tools}
        state = "SCORED" if scored else "not-yet-run"
        print(f"  {spec.label} ({spec.key}) — {state}")
        print(f"      receptor {recep or '(n/a)'}: "
              f"{'ok' if recep_ok else 'MISSING (prepare it — see RUNBOOK)'}")
        print(f"      tools: " + ", ".join(
            f"{t}={'ok' if have[t] else 'absent'}" for t in tools))
        if not scored:
            outstanding.append(spec)

    # Verdict.
    print("\n" + "=" * 70)
    if drift:
        print("VERDICT: DRIFT — a frozen hash does not match. Do NOT run; resolve first.")
        return 1
    if n and not outstanding:
        print("VERDICT: READY-TO-GRADE — poses + all baseline outputs present.")
        print("  run: python scripts/rescore_baselines.py " + pose_dir + " " + GRADE_RECEPTOR)
        return 0
    print("VERDICT: RUN-STEPS-REMAIN (the normal pre-run state).")
    if not n:
        print("  - regenerate poses: scripts/run_screen2.sh   (RUNBOOK step 1)")
    for spec in outstanding:
        print(f"  - run {spec.key}: scripts/run_{spec.key}_rescore.sh  ->  {spec.subdir}/")
    print("  - or just run the orchestrator: scripts/run_baselines_all.sh")
    print("See work/RUNBOOK_baselines_multi.md for the full pinned recipe.")
    return 2


if __name__ == "__main__":
    sys.exit(main())
