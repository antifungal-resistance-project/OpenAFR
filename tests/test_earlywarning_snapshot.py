"""Tests for the early-warning snapshot store (openafr/earlywarning.py, issue #21).

These lock the four guarantees the emergence layer (#22) will lean on:
  * normalization keys isolates stably and never invents a resistance call,
  * a snapshot is byte-stable / idempotent (so git diffs mean something and a
    re-pull is a no-op),
  * baseline_as_of filters by NCBI visibility date, honestly excluding undated
    isolates instead of silently placing them,
  * diff surfaces exactly the newly-appeared isolates.
No network: everything runs against synthetic NCBI-shaped rows. The release-URL
helpers are exercised against a stubbed _get (still no network) so the release
parsing and its refuse-on-empty guarantee are covered without hitting NCBI.
"""
import pytest

from openafr import earlywarning as ew


def _raw(target_acc, run="SRR1", geo="Pakistan", created="2020-01-01",
         collected="2019", biosample="SAMN1", asm="NULL"):
    """A minimal NCBI-metadata-shaped row (only the columns we read)."""
    return {
        "target_acc": target_acc, "Run": run, "geo_loc_name": geo,
        "target_creation_date": created, "collection_date": collected,
        "biosample_acc": biosample, "asm_acc": asm, "lat_lon": "NULL",
        # NCBI columns the spike proved empty; present so normalize() is realistic.
        "AMR_genotypes": "NULL", "computed_types": "NULL",
    }


def test_normalize_keys_on_versionless_target_acc():
    recs, dropped = ew.normalize([_raw("PDT000585001.1")])
    assert dropped == 0
    assert len(recs) == 1
    assert recs[0]["isolate_key"] == "PDT000585001"
    assert recs[0]["target_acc"] == "PDT000585001.1"


def test_normalize_splits_country_and_blanks_nullish():
    recs, _ = ew.normalize([_raw("PDT1.1", geo="USA: California", asm="null")])
    r = recs[0]
    assert r["country"] == "USA"
    assert r["geo_loc_name"] == "USA: California"
    assert r["asm_acc"] == ""  # "null" -> blank, not the literal string


def test_normalize_never_fabricates_a_resistance_call():
    # The load-bearing honesty check: NCBI exposes no azole call, so erg11_call
    # must stay empty and resistance_source must say why (reads vs no reads).
    with_reads = ew.normalize([_raw("PDT1.1", run="SRR9")])[0][0]
    no_reads = ew.normalize([_raw("PDT2.1", run="NULL")])[0][0]
    assert with_reads["erg11_call"] == "" and no_reads["erg11_call"] == ""
    assert with_reads["resistance_source"] == "pending:sra-recaller"
    assert no_reads["resistance_source"] == "pending:no-reads"


def test_normalize_drops_unkeyable_rows_and_counts_them():
    recs, dropped = ew.normalize([_raw(""), _raw("NULL"), _raw("PDT3.1")])
    assert dropped == 2
    assert [r["isolate_key"] for r in recs] == ["PDT3"]


def test_normalize_dedupes_keeping_highest_version():
    recs, _ = ew.normalize([_raw("PDT9.1", geo="Old"), _raw("PDT9.2", geo="New")])
    assert len(recs) == 1
    assert recs[0]["target_acc"] == "PDT9.2"
    assert recs[0]["country"] == "New"


def test_snapshot_is_byte_stable_and_order_independent(tmp_path):
    a = ew.normalize([_raw("PDT2.1"), _raw("PDT1.1"), _raw("PDT3.1")])[0]
    b = ew.normalize([_raw("PDT3.1"), _raw("PDT1.1"), _raw("PDT2.1")])[0]
    p1, p2 = tmp_path / "a.tsv", tmp_path / "b.tsv"
    s1 = ew.write_snapshot(p1, a)
    s2 = ew.write_snapshot(p2, b)
    assert s1 == s2  # same isolates, any input order -> identical bytes
    assert p1.read_text() == p2.read_text()
    # round-trips back to the same records
    assert ew.read_snapshot(p1) == a


def test_baseline_as_of_filters_by_creation_date_and_flags_undated():
    recs, _ = ew.normalize([
        _raw("PDT1.1", created="2019-06-01"),
        _raw("PDT2.1", created="2021-06-01"),
        _raw("PDT3.1", created="NULL"),        # undated
        _raw("PDT4.1", created="2020"),        # partial -> not placeable
    ])
    in_base, undated = ew.baseline_as_of(recs, "2020-01-01")
    assert {r["isolate_key"] for r in in_base} == {"PDT1"}
    assert {r["isolate_key"] for r in undated} == {"PDT3", "PDT4"}


def test_diff_surfaces_new_removed_and_revised():
    old = ew.normalize([_raw("PDT1.1"), _raw("PDT2.1")])[0]
    new = ew.normalize([_raw("PDT2.2"), _raw("PDT3.1")])[0]  # PDT1 gone, PDT3 new, PDT2 bumped
    d = ew.diff_snapshots(old, new)
    assert [r["isolate_key"] for r in d["added"]] == ["PDT3"]
    assert d["removed"] == ["PDT1"]
    assert d["revised"] == ["PDT2"]


def test_index_upsert_is_idempotent_per_release(tmp_path):
    idx = tmp_path / "INDEX.tsv"
    base = {"release_tag": "PDG000000067.671", "pulled_at_utc": "t1",
            "n_isolates": "10", "n_dropped": "0", "sha256": "abc",
            "source_url": "u", "tool_version": "1"}
    ew.update_index(idx, base)
    ew.update_index(idx, {**base, "pulled_at_utc": "t2", "sha256": "def"})
    rows = ew.read_snapshot(idx)  # same TSV reader
    assert len(rows) == 1  # re-pull replaced, not appended
    assert rows[0]["sha256"] == "def"
    ew.update_index(idx, {**base, "release_tag": "PDG000000067.672"})
    rows = ew.read_snapshot(idx)
    assert [r["release_tag"] for r in rows] == [
        "PDG000000067.671", "PDG000000067.672"]  # sorted by release number


# --- release discovery URLs (stubbed _get, no network) ----------------------

def test_metadata_url_builds_the_tagged_release_path():
    tag, url = ew.metadata_url(671)
    assert tag == "PDG000000067.671"
    assert url == (
        "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Candidozyma_auris/"
        "PDG000000067.671/Metadata/PDG000000067.671.metadata.tsv"
    )


def test_latest_release_picks_the_highest_published_number(monkeypatch):
    # NCBI's directory listing links every release; latest_release must return the
    # numeric max, not the last-listed or first-listed one.
    listing = (
        '<a href="PDG000000067.669/">669</a>'
        '<a href="PDG000000067.671/">671</a>'
        '<a href="PDG000000067.670/">670</a>'
    ).encode("utf-8")
    monkeypatch.setattr(ew, "_get", lambda url: listing)
    assert ew.latest_release() == 671


def test_latest_release_refuses_when_no_releases_are_listed(monkeypatch):
    # An empty/HTML-shaped listing must raise, never silently return a bogus release.
    monkeypatch.setattr(ew, "_get", lambda url: b"<html>no releases here</html>")
    with pytest.raises(RuntimeError, match="could not list releases"):
        ew.latest_release()


def test_fetch_release_rows_parses_the_metadata_tsv(monkeypatch):
    # With release pinned, fetch_release_rows must return DictReader rows plus the
    # tag/url it fetched -- no network, no latest_release lookup.
    tsv = b"target_acc\tRun\nPDT1.1\tSRR1\nPDT2.1\tSRR2\n"
    monkeypatch.setattr(ew, "_get", lambda url: tsv)
    rows, tag, url = ew.fetch_release_rows(671)
    assert tag == "PDG000000067.671"
    assert url.endswith("PDG000000067.671.metadata.tsv")
    assert [r["target_acc"] for r in rows] == ["PDT1.1", "PDT2.1"]
