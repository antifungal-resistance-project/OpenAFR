"""Source the run accessions behind the FKS1 benchmark paper (Misas et al., Microbiol Spectr
2025, PMC12323592, DOI 10.1128/spectrum.03147-24) WITHOUT the paywalled/CAPTCHA'd Table 1, and
verify each against ENA.

Why this exists. The concordance validation (work/PREREGISTRATION_fks1_concordance.md, engine
in openafr/concordance.py) is blocked on transcribing the paper's 100-row Table 1 into
data/earlywarning/recaller_sanity/fks1_benchmark_100.tsv. The paper's HTML is CAPTCHA-blocked and
an LLM scrape produced a self-contradictory row, so the rows cannot be auto-pulled. This script
recovers what CAN be sourced honestly, from machine APIs, with zero fabrication:

  1. Europe PMC annotations API text-mines accession numbers from the full text it cannot legally
     serve (the article is not open-access; it IS a US-Gov public-domain work). This yields the
     SRR run accessions + BioProjects the paper cites.
  2. ENA read_run confirms each run is ILLUMINA / WGS / PAIRED (the RUNBOOK's live-verify step)
     and maps it to its CDC isolate id (B-number) via sample title/alias.

What it CANNOT recover, and does not invent: the per-isolate FKS1 genotype (expected_panel),
echinocandin MIC, and R/S phenotype — the truth labels. Those live only in Table 1; BioSample
attributes carry isolate/geo/date/host but no AST/genotype (verified). So this emits a PARTIAL
skeleton with those columns marked PENDING, to be filled from the paper PDF (a public-domain
US-Gov work) and only then frozen at the graded fixture path. It deliberately writes to a
'*_partial.tsv' name, never the graded default, so it can't be mistaken for the frozen truth set.

Usage:
    python scripts/fks1_benchmark_accessions.py            # writes the partial skeleton + prints stats
"""
import json
import pathlib
import re
import sys
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "earlywarning" / "recaller_sanity" / "fks1_benchmark_srr_verified_partial.tsv"
PMCID = "PMC12323592"
CITATION = "PMC12323592;doi:10.1128/spectrum.03147-24"
ANNOT = ("https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds"
         f"?articleIds=PMC%3A{PMCID}&type=Accession%20Numbers&format=JSON")
ENA = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run"
       "&fields=run_accession,instrument_platform,library_strategy,library_layout,"
       "sample_accession,sample_title,sample_alias,study_accession&format=json")
CDC_ID = re.compile(r"^(USA-CDC-CAU-\S+|B\d{4,}\S*)$")


def _get(url, timeout=60):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.load(r)


def isolate_id(title, alias):
    for cand in (alias, title):
        if cand and CDC_ID.match(cand.strip()):
            return cand.strip()
    return (alias or title or "").strip() or "UNKNOWN"


def main():
    annots = _get(ANNOT)[0]["annotations"]
    srr = sorted({a["exact"] for a in annots if a["exact"].startswith("SRR")})
    prj = sorted({a["exact"] for a in annots if a["exact"].startswith("PRJ")})
    print(f"Europe PMC mined {len(srr)} SRR run accessions + {len(prj)} BioProjects: {prj}")

    rows = []
    for acc in srr:
        time.sleep(0.15)
        try:
            j = _get(ENA.format(acc=acc))
        except Exception as e:
            print(f"  {acc}: ENA error {e}", file=sys.stderr)
            continue
        if j:
            rows.append(j[0])

    bad = [r for r in rows if not (r.get("instrument_platform") == "ILLUMINA"
                                   and r.get("library_strategy") == "WGS"
                                   and r.get("library_layout") == "PAIRED")]
    print(f"ENA-verified {len(rows)}/{len(srr)} runs; non-ILLUMINA/WGS/PAIRED: {len(bad)}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write("# PARTIAL FKS1 benchmark skeleton — accessions sourced + ENA-verified, "
                 "TRUTH LABELS PENDING Table 1.\n")
        fh.write(f"# Source: {PMCID} (Misas et al., Microbiol Spectr 2025; US-Gov public domain). "
                 "Accessions text-mined via Europe PMC annotations API; run metadata from ENA.\n")
        fh.write(f"# Mined {len(srr)}/100 isolates (miner missed ~33 Table-1 cells). BioProjects: "
                 f"{','.join(prj)}\n")
        fh.write("# expected_panel/susceptibility/clade = PENDING: fill from the paper PDF Table 1, "
                 "then freeze at fks1_benchmark_100.tsv. DO NOT grade this file.\n")
        fh.write("# strain\tclade\trun_acc\texpected_panel\tsusceptibility\tnote\tcitation\n")
        for r in sorted(rows, key=lambda r: r["run_accession"]):
            strain = isolate_id(r.get("sample_title"), r.get("sample_alias"))
            note = (f"bioproject={r.get('study_accession')};sample={r.get('sample_accession')};"
                    f"ena_verified=ILLUMINA/WGS/PAIRED")
            fh.write(f"{strain}\tPENDING\t{r['run_accession']}\tPENDING\tPENDING\t{note}\t{CITATION}\n")
    print(f"wrote {OUT.relative_to(ROOT)} ({len(rows)} verified runs, labels PENDING)")


if __name__ == "__main__":
    main()
