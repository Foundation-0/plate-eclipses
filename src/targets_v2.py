"""TARGETS-V2: expanded, ranked DASCH plate-search target list beyond the 31 ASAS-SN big dippers.

Every row comes from a catalogue query recorded below (source, ADQL, run date printed in the
output header). No row is hand-entered. Run from src/:

    SSL_CERT_FILE=/etc/ssl/cert.pem python3 targets_v2.py [--no-dasch] [--max-dasch N]

Sources (all queried live; failures are reported, never silently replaced):
  VSX     AAVSO VSX as mirrored by VizieR TAP table "B/vsx/vsx" (AAVSO direct site returned
          HTTP 403 from this machine on 2026-10-18, so the VizieR mirror is used and labelled).
          EA* types, max <= 13.0 mag, amplitude (min-max) >= 0.5 mag, and
          (Period IS NULL or Period >= 100 d).
  ASASSN  ASAS-SN Catalog of Variable Stars via VizieR TAP "II/366/catv2021":
          Type EA/EB, APASS B (BmagA) <= 13.5, Amp >= 0.5 mag, Per >= 30 d (the long tail of
          that catalogue; its bright deep EA/EB periods top out near 95 d).
  GAIA    Gaia DR3 gaiadr3.vari_eclipsing_binary joined to gaiadr3.gaia_source via ESA TAP:
          365 <= P=1/frequency <= 5000 d, G < 13.0, global_ranking > 0.5,
          derived_primary_ecl_depth > 0.25. Depth used for scoring is the geometric-model
          Gaussian depth, whose TAP description is explicitly "Magnitude depth".
          (--gaia-release dr4 switches the schema prefix for 2 Dec 2026; untested until release.)

DASCH coverage proxy: each deduped target is matched to DASCH APASS exactly as fetch.py does
(querycat 90", nearest; here accepted to 15" because VSX coordinates are lower precision, and
the separation is recorded). The target light curve alone (no comparison stars) is parsed with
core.parse/core.base; "normal" plates are ok detections within depth/2 of the raw median, and
the caught fraction reuses coverage.py's formula verbatim (fraction of 1890-1990 start epochs
whose window of length D contains a normal plate). This is a PROXY: no per-plate comparison
zero-point and no non-detections, unlike coverage.py proper. It ranks fields; it is not the
measurement coverage.py will produce once fetch.py has run on a target.

Score = caught_frac(duration_used) * min(depth_mag, 2.0). duration_used is measured only for
Gaia rows (derived_primary_ecl_duration * P); for VSX/ASAS-SN rows no duration is catalogued,
so a 30 d reference duration is used and flagged assumed. caught_30/caught_100 are always
reported for comparability. Rows with no DASCH proxy get score empty and sort last by a
brightness/depth prior (documented in results/targets_v2.md).

Exclusions: anything within 120" of an original 31 target (targets.T) is dropped and listed in
the run log (this is also a query sanity check: ASASSN-21co should appear there).

Tiers: A = VSX period unknown; B = VSX P>=365 d; B2 = VSX 100<=P<365 d;
C = ASAS-SN 30-95 d; D = Gaia DR3 long-P (sparse-sampling period, alias risk).

Outputs: ../results/targets_v2.csv and ../results/targets_v2.md
"""
import argparse, bisect, csv, io, json, math, os, sys, time, urllib.parse, urllib.request
from datetime import date

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from targets import T as ORIGINALS  # noqa: E402

VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
J0, J1 = 2411368.5, 2447892.5  # 1890.0 - 1990.0, as in coverage.py
UA = {"User-Agent": "foundation0-plate-eclipses targets_v2"}


def tap(base, query, timeout=180, tries=5):
    # POST (http GET to TAPVizieR reset repeatedly from this machine on 2026-10-18; POST/https worked)
    body = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}).encode()
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(base, data=body, headers={
                **UA, "Content-Type": "application/x-www-form-urlencoded"})
            text = urllib.request.urlopen(req, timeout=timeout).read().decode("utf-8", "replace")
            if "QUERY_STATUS\" value=\"ERROR" in text or text.lstrip().startswith("<?xml"):
                raise RuntimeError("TAP error: " + text[:400])
            return list(csv.DictReader(io.StringIO(text)))
        except Exception as e:  # noqa: BLE001 - reported to caller, source marked failed
            last = e
            time.sleep(2 ** (k + 1))
    raise RuntimeError(str(last))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def fetch_vsx(log):
    q = ("SELECT OID, Name, RAJ2000, DEJ2000, Type, max, n_max, min, n_min, Epoch, "
         "u_Period, l_Period, Period FROM \"B/vsx/vsx\" WHERE Type LIKE 'EA%' AND max <= 13.0 "
         "AND min - max >= 0.5 AND (Period IS NULL OR Period >= 100)")
    rows = tap(VIZIER_TAP, q)
    out = []
    for r in rows:
        period = fnum(r.get("Period"))
        depth = fnum(r.get("min")) - fnum(r.get("max"))
        tier = "A" if period is None else ("B" if period >= 365 else "B2")
        status = "period_unknown_in_vsx" if period is None else "catalogue_period"
        if (r.get("u_Period") or "").strip() == ":":
            status += ";period_flagged_uncertain"
        if ":" in (r.get("Type") or ""):
            status += ";type_uncertain"
        out.append(dict(name=(r.get("Name") or "").strip(), ra=fnum(r["RAJ2000"]),
                        dec=fnum(r["DEJ2000"]), source="VSX(B/vsx/vsx via VizieR TAP)",
                        source_id=f"OID {r.get('OID')}", type=(r.get("Type") or "").strip(),
                        mag_bright=fnum(r.get("max")), mag_band=(r.get("n_max") or "").strip(),
                        depth_mag=round(depth, 3), period_d=period, epoch_jd=fnum(r.get("Epoch")),
                        duration_d=None, period_status=status, tier=tier))
    log.append(f"VSX: {len(out)} rows")
    return out


def fetch_asassn(log):
    q = ("SELECT ID, \"ASASSN-V\" AS asassnv, RAJ2000, DEJ2000, Vmag, BmagA, Amp, Per, Type, "
         "Prob, HJD, GaiaDR3, TIC FROM \"II/366/catv2021\" WHERE (Type LIKE 'EA%' OR "
         "Type LIKE 'EB%') AND BmagA <= 13.5 AND Amp >= 0.5 AND Per >= 30")
    rows = tap(VIZIER_TAP, q)
    out = []
    for r in rows:
        out.append(dict(name=f"ASASSN-V {(r.get('asassnv') or '').strip()}",
                        ra=fnum(r["RAJ2000"]), dec=fnum(r["DEJ2000"]),
                        source="ASAS-SN VarCat(II/366/catv2021 via VizieR TAP)",
                        source_id=f"ID {r.get('ID')};GaiaDR3 {r.get('GaiaDR3')};TIC {r.get('TIC')}",
                        type=(r.get("Type") or "").strip(), mag_bright=fnum(r.get("BmagA")),
                        mag_band="B(APASS)", depth_mag=fnum(r.get("Amp")),
                        period_d=fnum(r.get("Per")), epoch_jd=fnum(r.get("HJD")),
                        duration_d=None,
                        period_status=f"asassn_catalogue_period;prob={r.get('Prob')}", tier="C"))
    log.append(f"ASAS-SN: {len(out)} rows")
    return out


def fetch_gaia(log, release="dr3"):
    s = f"gaia{release}"
    q = (f"SELECT eb.source_id, g.ra, g.dec, g.phot_g_mean_mag, g.phot_bp_mean_mag, "
         f"g.phot_rp_mean_mag, eb.frequency, eb.frequency_error, eb.global_ranking, "
         f"eb.geom_model_gaussian1_depth, eb.geom_model_gaussian2_depth, "
         f"eb.derived_primary_ecl_duration, eb.derived_primary_ecl_depth, "
         f"eb.reference_time FROM {s}.vari_eclipsing_binary AS eb JOIN {s}.gaia_source AS g "
         f"ON eb.source_id = g.source_id WHERE eb.frequency BETWEEN 0.0002 AND 0.00274 "
         f"AND g.phot_g_mean_mag < 13.0 AND eb.derived_primary_ecl_depth > 0.25 "
         f"AND eb.global_ranking > 0.5")
    rows = tap(GAIA_TAP, q)
    out = []
    for r in rows:
        freq = fnum(r["frequency"]); period = 1.0 / freq if freq else None
        depths = [d for d in (fnum(r.get("geom_model_gaussian1_depth")),
                              fnum(r.get("geom_model_gaussian2_depth"))) if d is not None]
        dur_phase = fnum(r.get("derived_primary_ecl_duration"))
        out.append(dict(name=f"Gaia {release.upper()} {r.get('source_id')}",
                        ra=fnum(r["ra"]), dec=fnum(r["dec"]),
                        source=f"Gaia {release.upper()} vari_eclipsing_binary (ESA TAP join gaia_source)",
                        source_id=f"source_id {r.get('source_id')}", type="EB(geometric model)",
                        mag_bright=fnum(r.get("phot_g_mean_mag")), mag_band="G",
                        depth_mag=round(max(depths), 3) if depths else None,
                        period_d=round(period, 3) if period else None,
                        epoch_jd=None,
                        duration_d=round(dur_phase * period, 3) if (dur_phase and period) else None,
                        period_status=("gaia_geometric_model_period_sparse_sampling_alias_risk;"
                                       f"global_ranking={r.get('global_ranking')};"
                                       f"freq_err={r.get('frequency_error')}"), tier="D"))
    log.append(f"Gaia {release}: {len(out)} rows")
    return out


def sep_arcsec(ra1, dec1, ra2, dec2):
    if None in (ra1, dec1, ra2, dec2):
        return 1e9
    dra = (ra1 - ra2) * math.cos(math.radians(dec1))
    return math.hypot(dra, dec1 - dec2) * 3600.0


def caught_frac(normal_jds, dur_d):
    """coverage.py lines 19-22, with bisect in place of numpy searchsorted."""
    if not normal_jds:
        return 0.0
    caught = 0; total = 0; t = J0
    while t < J1:
        i = bisect.bisect_left(normal_jds, t)
        if i < len(normal_jds) and normal_jds[i] <= t + dur_d:
            caught += 1
        total += 1; t += 1.0
    return caught / total


def dasch_proxy(row, log):
    """Fill dasch_* and caught_* fields from the target DASCH light curve alone. See docstring."""
    row.update(dasch_status="", dasch_ref="", dasch_sep_arcsec="", dasch_B="",
               dasch_num_matches="", n_normal_plates="", caught_30="", caught_100="",
               caught_at_duration="")
    try:
        from dasch import querycat, lightcurve
        from core import parse, base
    except Exception as e:  # noqa: BLE001
        row["dasch_status"] = f"import_failed:{e}"; return
    try:
        near = querycat(row["ra"], row["dec"], 90)
        for r in near:
            r["sep"] = math.hypot(float(r["dra_asec"]), float(r["ddec_asec"]))
        near.sort(key=lambda r: r["sep"])
        if not near or near[0]["sep"] > 15:
            row["dasch_status"] = "no_apass_match_within_15arcsec"; return
        t = near[0]
        row["dasch_ref"] = t.get("ref_text", "")
        row["dasch_sep_arcsec"] = round(float(t["sep"]), 2)
        row["dasch_B"] = t.get("stdmag", "")
        row["dasch_num_matches"] = t.get("num_matches", "")
        lightcurve(t)
        det, _lim = parse(t["ref_number"])
        med, _mad = base(det)
        if med is None:
            row["dasch_status"] = "fewer_than_20_ok_detections"; return
        depth = row.get("depth_mag") or 0.5
        normal = sorted(p["jd"] for p in det if p["ok"] and (p["m"] - med) < depth / 2)
        row["n_normal_plates"] = len(normal)
        durations = {30: "caught_30", 100: "caught_100"}
        for d, key in durations.items():
            row[key] = round(caught_frac(normal, d), 3)
        row["caught_at_duration"] = round(caught_frac(normal, max(row["duration_used_d"], 1)), 3)
        row["dasch_status"] = "proxy_ok" + (";sep_gt_5arcsec_match_uncertain"
                                            if float(t["sep"]) > 5 else "")
        if fnum(t.get("stdmag")) and fnum(t["stdmag"]) > 13.5 and t["stdmag"] != "99.9":
            row["dasch_status"] += ";apass_B_gt_13.5_faint_for_plates"
    except Exception as e:  # noqa: BLE001 - per-target failure must not kill the run
        row["dasch_status"] = f"failed:{type(e).__name__}:{str(e)[:120]}"
        log.append(f"DASCH proxy failed for {row['name']}: {e}")


FIELDS = ["rank", "tier", "name", "ra_deg", "dec_deg", "source", "source_id", "also_in",
          "type", "mag_bright", "mag_band", "depth_mag", "period_d", "period_status",
          "epoch_jd", "duration_used_d", "duration_basis", "dasch_status", "dasch_ref",
          "dasch_sep_arcsec", "dasch_B", "dasch_num_matches", "n_normal_plates",
          "caught_30", "caught_100", "caught_at_duration", "score", "score_basis"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-dasch", action="store_true", help="skip DASCH proxy (catalogue only)")
    ap.add_argument("--max-dasch", type=int, default=0, help="cap DASCH proxy calls (0 = all)")
    ap.add_argument("--gaia-release", default="dr3", choices=["dr3", "dr4"])
    args = ap.parse_args()
    run_date = date.today().isoformat()
    log = [f"run_date={run_date}"]
    rows, failed_sources = [], []
    for fn, label in ((fetch_vsx, "VSX"), (fetch_asassn, "ASAS-SN")):
        try:
            rows += fn(log)
        except Exception as e:  # noqa: BLE001
            failed_sources.append(f"{label}: {e}"); log.append(f"{label} FAILED: {e}")
    try:
        rows += fetch_gaia(log, args.gaia_release)
    except Exception as e:  # noqa: BLE001
        failed_sources.append(f"Gaia: {e}"); log.append(f"Gaia FAILED: {e}")

    # Exclude the original 31 (120" positional), log them as the query sanity check.
    kept, excluded = [], []
    for r in rows:
        hit = next((o[0] for o in ORIGINALS
                    if sep_arcsec(r["ra"], r["dec"], o[1], o[2]) < 120), None)
        (excluded.append(f"{r['name']} [{r['source']}] ~ original {hit}") if hit
         else kept.append(r))
    log.append(f"excluded as original-31 matches ({len(excluded)}): " + "; ".join(excluded))

    # Cross-source dedup at 60": keep first (VSX > ASAS-SN > Gaia order), record also_in.
    deduped = []
    for r in kept:
        dup = next((d for d in deduped
                    if sep_arcsec(r["ra"], r["dec"], d["ra"], d["dec"]) < 60), None)
        if dup:
            dup["also_in"] = (dup.get("also_in", "") + f";{r['source']}:{r['source_id']}").strip(";")
            if dup.get("duration_d") is None and r.get("duration_d") is not None:
                dup["duration_d"] = r["duration_d"]
        else:
            r["also_in"] = ""; deduped.append(r)
    log.append(f"after dedup: {len(deduped)} targets")

    for r in deduped:
        if r.get("duration_d"):
            r["duration_used_d"] = r["duration_d"]; r["duration_basis"] = "catalogue_measured"
        else:
            r["duration_used_d"] = 30; r["duration_basis"] = "assumed_reference_30d_not_catalogued"
        r.update(dasch_status="not_run", dasch_ref="", dasch_sep_arcsec="", dasch_B="",
                 dasch_num_matches="", n_normal_plates="", caught_30="", caught_100="",
                 caught_at_duration="", score="", score_basis="")

    if not args.no_dasch:
        n = 0
        for r in deduped:
            if args.max_dasch and n >= args.max_dasch:
                r["dasch_status"] = "skipped_max_dasch_cap"; continue
            dasch_proxy(r, log); n += 1
            print(f"[{n}/{len(deduped)}] {r['name']}: {r['dasch_status']} "
                  f"caught30={r.get('caught_30')}", flush=True)

    for r in deduped:
        depth = r.get("depth_mag")
        if r.get("caught_at_duration") != "" and depth:
            r["score"] = round(r["caught_at_duration"] * min(depth, 2.0), 4)
            r["score_basis"] = "caught_frac(duration_used)*min(depth,2) [DASCH proxy]"
        else:
            bright = r.get("mag_bright") or 14.5
            prior = round(min(depth or 0, 2.0) * max(0.0, min((14.5 - bright) / 3.0, 1.0)), 4)
            r["score"] = prior
            r["score_basis"] = "PRIOR ONLY (no DASCH proxy): min(depth,2)*clip((14.5-mag)/3,0,1)"

    tier_order = {"A": 0, "B": 1, "B2": 2, "C": 3, "D": 4}
    # proxy-scored rows first (test the data field, not the basis text: the PRIOR label
    # itself contains the words "DASCH proxy")
    deduped.sort(key=lambda r: (0 if r.get("caught_at_duration") != "" else 1,
                                -(r["score"] or 0), tier_order.get(r["tier"], 9)))
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
    csv_path = os.path.join(outdir, "targets_v2.csv")
    with open(csv_path, "w", newline="") as fh:
        fh.write(f"# TARGETS-V2 generated by src/targets_v2.py on {run_date}; "
                 f"see results/targets_v2.md for queries, scoring and claim labels\n")
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for i, r in enumerate(deduped, 1):
            w.writerow(dict(rank=i, ra_deg=r["ra"], dec_deg=r["dec"], **r))
    md_path = os.path.join(outdir, "targets_v2.md")
    write_md(md_path, deduped, log, failed_sources, run_date, args)
    print("\n".join(log))
    print(f"wrote {csv_path} ({len(deduped)} rows) and {md_path}")


def write_md(path, rows, log, failed_sources, run_date, args):
    tiers = {}
    for r in rows:
        tiers[r["tier"]] = tiers.get(r["tier"], 0) + 1
    n_proxy = sum(1 for r in rows if r.get("caught_at_duration") != "")
    lines = [f"# TARGETS-V2 — expanded ranked DASCH target list (generated {run_date})", "",
             "Generated by [`src/targets_v2.py`](../src/targets_v2.py); data in "
             "[`targets_v2.csv`](targets_v2.csv). Every row carries its catalogue source and "
             "source_id; no row is hand-entered. Claim labels: catalogue values are KNOWN "
             "(as catalogued, not re-derived here); DASCH caught fractions are VERIFIED only as "
             "the documented target-only proxy computation on fetched light curves; scores that "
             "use the 30 d reference duration or the brightness prior are CONJECTURE-rankings "
             "for ordering follow-up, not measured detectabilities.", "",
             f"- Rows: {len(rows)} (tiers: " +
             ", ".join(f"{k}={v}" for k, v in sorted(tiers.items())) + ")",
             f"- DASCH proxy computed: {n_proxy}/{len(rows)}",
             f"- Failed sources this run: {failed_sources or 'none'}", "",
             "## Queries (exact; rerun reproduces the catalogue side)", "",
             "- VSX (VizieR TAP `B/vsx/vsx`, mirror of AAVSO VSX): `Type LIKE 'EA%' AND "
             "max<=13.0 AND min-max>=0.5 AND (Period IS NULL OR Period>=100)`.",
             "- ASAS-SN (VizieR TAP `II/366/catv2021`): `Type EA/EB AND BmagA<=13.5 AND "
             "Amp>=0.5 AND Per>=30`.",
             f"- Gaia (ESA TAP `gaia{args.gaia_release}.vari_eclipsing_binary` JOIN "
             "`gaia_source`): `365<=1/frequency<=5000 d AND G<13.0 AND global_ranking>0.5 AND "
             "derived_primary_ecl_depth>0.25`.",
             "- Original 31 excluded positionally at 120\"; cross-source dedup at 60\".", "",
             "## Scoring", "",
             "`score = caught_frac(duration_used) * min(depth_mag, 2.0)`. caught_frac reuses "
             "[`src/coverage.py`](../src/coverage.py)'s definition on target-only normal plates "
             "(no comparison-star zero-point, no non-detections — a ranking proxy). "
             "duration_used is catalogue-measured only for Gaia rows; VSX/ASAS-SN rows use an "
             "assumed 30 d reference (flagged per row). Rows without a DASCH proxy carry a "
             "brightness/depth PRIOR score and sort after all proxy-scored rows.", "",
             "Proxy validation (2026-10-18, one star): running this proxy on ASASSN-21co "
             "(excluded from the list as an original) gives 555 normal plates and caught "
             "fractions 0.192/0.313/0.337 at 30/86/100 d, against the comparison-corrected "
             "[`coverage.txt`](coverage.txt) values 560 and 0.19/0.32/0.34. VERIFIED for that "
             "one field; fields with strong plate zero-point scatter may differ more.", "",
             "## Run log", "", "```"] + log + ["```", ""]
    lines += ["## What this list is not", "",
              "- Not a DASCH search: no eclipse is claimed for any row. Next step per target is "
              "the existing fetch.py/events.py/period.py pipeline with comparison stars.",
              "- VSX EA: entries (incl. NSV suspects) may be misclassified or spurious; the "
              "type_uncertain flag is per row. Period-known Tier B rows include well-studied "
              "systems — check the literature before claiming novelty for any plate event.",
              "- Gaia Tier D periods come from sparse sampling and carry alias risk; their "
              "model eclipse durations are short, so most rank low for plates by design.", ""]
    open(path, "w").write("\n".join(lines))


if __name__ == "__main__":
    main()
