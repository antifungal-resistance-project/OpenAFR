"""Cert grader for FKS1 structural contact label (moonshot #4, Rung 3).

Validates that core-tier panel POSITIONS cluster in `drug_contact` / `pocket_adjacent`
proximity tier in the 9WZU contact map — anti-circularly, since the panel was determined
from clinical MIC studies while the contact map was derived from cryo-EM drug coordinates.

Pre-registered pass bar (work/PREREGISTRATION_fks1_structural_verdict.md):
  PASS:         >= 80% of core panel positions in drug_contact OR pocket_adjacent
  FAIL:         <  80%
  UNDERPOWERED: fewer than 3 core panel positions with a contact-map entry

Writes work/RESULTS_fks1_structural_verdict.md on completion.
"""
import hashlib
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from openafr import fks1_caller, fks1_structural  # noqa: E402

PREREG = ROOT / "work" / "PREREGISTRATION_fks1_structural_verdict.md"
RESULTS = ROOT / "work" / "RESULTS_fks1_structural_verdict.md"
PASS_BAR = 0.80
UNDERPOWERED_THRESHOLD = 3   # fewer than this many mapped core positions → UNDERPOWERED

# Pre-committed predictions from the preregistration (frozen before this script ran).
# Each key is a Ca residue position; value is the expected proximity tier.
EXPECTED_TIERS = {
    635: fks1_structural.DRUG_CONTACT,
    639: fks1_structural.POCKET_ADJACENT,
    691: fks1_structural.DRUG_CONTACT,
    1354: fks1_structural.POCKET_ADJACENT,
}


def _core_panel_positions():
    """All residue positions in the core (high-PPV) resistance panel."""
    low_ppv = set(fks1_caller.LOW_PPV_PANEL.keys())
    positions = set()
    for w in fks1_caller.FKS1_WINDOWS.values():
        for pos in w.panel:
            if pos not in low_ppv:
                positions.add(pos)
    return sorted(positions)


def _lookup_tier(pos):
    """Proximity tier for a position via direct contact-map lookup (bypasses the
    'skip core panel markers' guard in structural_proximity, which is correct for
    alert generation but wrong here — we need the tier for the core positions
    themselves)."""
    m = fks1_structural._load_map()
    entry = m["residues"].get(str(pos))
    if entry is None:
        return fks1_structural.NO_STRUCTURE, None
    return entry["tier"], entry.get("min_dist_to_drug_A")


def _wt_of(pos):
    for w in fks1_caller.FKS1_WINDOWS.values():
        if pos in w.panel:
            return w.panel[pos][0]
    return "?"


def run():
    prereg_sha = hashlib.sha256(PREREG.read_bytes()).hexdigest()

    core_positions = _core_panel_positions()
    in_contact_zone = []
    rows = []

    for pos in core_positions:
        tier, dist = _lookup_tier(pos)
        wt = _wt_of(pos)
        in_zone = tier in (fks1_structural.DRUG_CONTACT, fks1_structural.POCKET_ADJACENT)
        if in_zone:
            in_contact_zone.append(pos)
        expected = EXPECTED_TIERS.get(pos)
        match = "✓" if tier == expected else "✗ MISMATCH"
        dist_str = f"{dist:.2f} Å" if dist is not None else "N/A"
        rows.append((pos, wt, tier, dist_str, "YES" if in_zone else "no", match))

    n_mapped = sum(1 for _, _, t, _, _, _ in rows if t != fks1_structural.NO_STRUCTURE)
    n_in_zone = len(in_contact_zone)
    n_total = len(core_positions)

    if n_mapped < UNDERPOWERED_THRESHOLD:
        verdict = "UNDERPOWERED"
    elif n_in_zone / n_mapped >= PASS_BAR:
        verdict = "PASS"
    else:
        verdict = "FAIL"

    pct = 100.0 * n_in_zone / n_mapped if n_mapped else 0.0

    # --- write results ---
    lines = [
        "# Results — FKS1 structural contact label cert (moonshot #4, Rung 3)\n\n",
        f"> **Pre-registered bar:** ≥ 80% of core panel positions in "
        f"`drug_contact` OR `pocket_adjacent`  \n"
        f"> **Pre-registration SHA-256:** `{prereg_sha}`  \n"
        f"> **Contact map SHA-256:** `{fks1_structural.CONTACT_MAP_SHA256}`\n\n",
        "## Core panel position tiers\n\n",
        "| Ca position | WT | Proximity tier | Distance | In contact zone | "
        "Matches prediction |\n",
        "|---|---|---|---|---|---|\n",
    ]
    for pos, wt, tier, dist_str, zone, match in rows:
        lines.append(f"| {pos} | {wt} | `{tier}` | {dist_str} | {zone} | {match} |\n")

    lines += [
        "\n## Secondary observation — low-PPV positions (exploratory, not in bar)\n\n",
        "| Ca position | WT | Proximity tier | Distance |\n",
        "|---|---|---|---|\n",
    ]
    for pos, ppv_muts in sorted(fks1_caller.LOW_PPV_PANEL.items()):
        tier, dist = _lookup_tier(pos)
        wt = _wt_of(pos)
        dist_str = f"{dist:.2f} Å" if dist is not None else "N/A"
        lines.append(f"| {pos} | {wt} | `{tier}` | {dist_str} |\n")

    lines += [
        "\n## Summary\n\n",
        f"| Metric | Value |\n",
        f"|---|---|\n",
        f"| Core panel positions checked | {n_total} |\n",
        f"| Positions with contact-map entry | {n_mapped} |\n",
        f"| Positions in drug\\_contact or pocket\\_adjacent | {n_in_zone} |\n",
        f"| Fraction | {pct:.1f}% ({n_in_zone}/{n_mapped}) |\n",
        f"| Pre-registered bar | ≥ 80% |\n\n",
        f"## VERDICT: {verdict}\n\n",
    ]

    if verdict == "PASS":
        lines.append(
            "Core resistance panel positions cluster in the drug-contact / adjacent "
            "zone of the 9WZU cryo-EM structure. The contact map is structurally "
            "relevant to the known-resistance panel, anti-circularly validated.\n\n"
            "**Consequence (per preregistration):** implement `structural_verdict()` "
            "and wire into the pipeline.\n"
        )
    elif verdict == "FAIL":
        lines.append(
            "`structural` stays null for FKS1; proximity display remains Rung-2-only. "
            "Re-open only with a new structural source.\n"
        )
    else:
        lines.append(
            "Insufficient contact-map coverage of core panel positions. Treat as FAIL.\n"
        )

    RESULTS.write_text("".join(lines))
    print("".join(lines))
    print(f"VERDICT: {verdict}")
    return verdict == "PASS"


if __name__ == "__main__":
    ok = run()
    sys.exit(0 if ok else 1)
