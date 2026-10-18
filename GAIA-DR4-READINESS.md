# Gaia DR4 readiness: day-one TESS cross-check of DR4 transit and eclipsing-binary candidates

Plan for [astr0 AREAS item 2](https://github.com/Foundation-0/astr0/blob/main/AREAS.md)
("Gaia DR4 on day one, checked against TESS"), written 2026-10-18. It also feeds AREAS item 1:
long-period DR4 eclipsing binaries (EBs) that TESS cannot settle become DASCH targets via
[`src/targets_v2.py --gaia-release dr4`](src/targets_v2.py).

Claim labels follow this repo's convention: **KNOWN** = published/catalogued fact with source;
**VERIFIED** = finite check actually run (commands in this repo's history / targets_v2 run log);
**CONJECTURE** = expected but unconfirmed until the named check runs. Nothing here is proved in
the mathematical sense.

## Fixed dates and contents (KNOWN, to re-check before 25 Nov)

| Item | Value | Source |
|---|---|---|
| Gaia DR4 release | 2 Dec 2026 | [ESA Gaia release page](https://www.cosmos.esa.int/web/gaia/release), via astr0 AREAS item 2 |
| Gaia archive closure before release | 25 Nov 2026, 12:00 CET | same |
| DR4 size | about 400 TB, incl. epoch photometry for every source | [ESA Gaia DR4 page](https://www.cosmos.esa.int/web/gaia/dr4), via AREAS |
| DR4 planetary-transit candidates | 1,007 | same |
| DR4 eclipsing binaries | 2,392,500 | same |
| Existing local asset | TESS light curves already held/vetted for about 4 million stars (eumemic's pipeline) | astr0 AREAS item 2 |

These numbers are leads, not facts of record (AREAS says so itself): re-read both ESA pages in
the pre-25-Nov checklist below and update this table if they moved.

## The day-one question

For each DR4 transit candidate and each (filtered) DR4 EB: does the TESS light curve, folded on
Gaia's period and epoch, show the same event at the same phase, depth-consistent within
bandpass differences? Gaia samples sparsely (tens of epochs over years); TESS samples densely
over about 27 d sectors. Agreement at the Gaia-predicted phase is strong confirmation; a clean
TESS fold with coverage of the predicted event and no dip is refutation; no TESS coverage of
the phase is "undecided", and long-period undecided EBs route to DASCH (this repo).

## TIC <-> Gaia matching approach

Build the match table **before** 25 Nov against DR3, then re-point it at DR4 on release day.

1. **Start from TIC's own Gaia column.** TIC v8.x carries a Gaia source id per star, but it is a
   **Gaia DR2** id (KNOWN from the TIC definition, Stassun et al. 2019, AJ 158, 138; re-confirm
   the exact column name in the TIC build actually used locally before relying on it).
2. **Bridge DR2 -> DR3 with ESA's neighbourhood table.** `gaiadr3.dr2_neighbourhood`
   (`dr2_source_id`, `dr3_source_id`, plus angular distance and magnitude-difference columns)
   exists and answers on ESA TAP — VERIFIED 2026-10-18 (`SELECT TOP 1 ...` returned a row).
   Prefer the neighbourhood row flagged/selected as the best match; record rows where the DR2
   source splits into several DR3 sources (or merges) as match-ambiguous, never silently pick.
3. **Bridge DR3 -> DR4 the same way, name unconfirmed.** A `gaiadr4.dr3_neighbourhood` table
   with the same shape is expected (CONJECTURE — DR3 followed this pattern from DR2, but the DR4
   schema is not public pre-release). Day-one smoke test: `SELECT TOP 1 * FROM
   gaiadr4.dr3_neighbourhood`. If the name differs, list `TAP_SCHEMA.tables` for
   `gaiadr4.%neighbourhood%` before writing any join.
4. **Positional fallback for stars the ID bridge misses** (no TIC Gaia id, ambiguous split):
   2 arcsec cone on the Gaia position propagated to the TIC epoch with Gaia proper motion,
   nearest candidate with |G - Tmag_expected| sanity (TESS magnitude from TIC within about
   1.5 mag of a colour-based estimate from BP/RP), and exactly one candidate inside the cone;
   everything else goes to a manual-review list. Do **not** widen the cone to raise match
   rate: a wrong match folds the wrong light curve and manufactures a refutation.
5. **Reverse direction for DR4 candidates:** DR4 source_id -> (neighbourhood back-bridge or
   direct positional match to TIC at 2 arcsec with proper-motion propagation to the TIC epoch)
   -> TIC id -> local TESS light curve. Record match_method per candidate
   (`id_bridge`, `positional_unique`, `positional_ambiguous`, `no_match`) so day-one results
   can be filtered by match quality.
6. **Volume of the match job:** bridging is a table join, not a per-star service call. The
   neighbourhood tables are of order 1-2 billion rows in total (ESTIMATE from DR3 practice),
   so restrict with an uploaded/local list of the about 4 million TIC ids of interest (TAP
   upload join or a downloaded subset keyed by source_id) rather than downloading whole
   neighbourhood tables. Result to cache: one row per local TESS star
   (`tic_id, gaia_dr2_id, gaia_dr3_id, gaia_dr4_id, match_method, sep_arcsec`); expected local
   size well under 1 GB (ESTIMATE: 4e6 rows x about 100 bytes).

## Folding / vetting checker spec

To implement where the TESS light-curve cache lives (eumemic's pipeline), not in this repo —
this document is the spec the checker is built and reviewed against. Nothing below exists yet
(NOT DONE); treat the acceptance tests as part of the spec.

Input per candidate: Gaia source_id, period P and its uncertainty or frequency error, reference
epoch of primary eclipse/transit (convert Gaia's reference_time + phase columns to a BJD
epoch explicitly — DR3's `vari_eclipsing_binary` gives `reference_time` [BJD] plus model
phases, and `vari_planetary_transit` gives transit parameters; the epoch-conversion step is
where a phase-0.5 or UTC/BJD slip would hide, so it gets its own unit test), model depth and
duration, G/BP/RP magnitudes, match_method.

For each matched candidate with a TESS light curve:

1. **Phase-coverage gate first.** Compute which phases of the Gaia ephemeris the TESS sectors
   actually cover, widening the predicted event window by the accumulated period uncertainty
   between the Gaia epoch and the TESS sector epochs, plus the model duration. If the widened
   window is not covered: verdict `undecided_no_phase_coverage` (and if P > about 300 d and
   the target is bright, queue for DASCH via targets_v2). Never refute on uncovered phase.
2. **Fold and measure.** Detrend as the existing pipeline does, fold on P, and in parallel fold
   on P/2 and 2P (Gaia sparse sampling both misses factors of two and aliases). Measure in the
   predicted window vs out-of-event: median depth, SNR of a box/trapezoid matched filter with
   the Gaia duration, and the phase offset of the deepest bin.
3. **Verdicts (fixed vocabulary):**
   - `confirmed_tess_event`: dip in the predicted window, matched-filter SNR >= 7 (pipeline's
     existing detection threshold if different — one threshold, named in code, not tuned per
     candidate), depth within a factor of 3 of Gaia's model depth after no bandpass correction
     (wider tolerance is a follow-up question, not a day-one knob), no equally deep dip at the
     opposite phase for transits (that pattern marks an EB, verdict
     `tess_shows_eb_not_transit`), centroid/blend checks passed by the existing vetting.
   - `refuted_no_tess_event`: widened window fully covered, in-window median consistent with
     out-of-event within 3 sigma, and the 5-sigma depth upper limit below half the Gaia model
     depth. Report the upper limit, not just the label.
   - `period_alias_candidate`: event confirmed only at P/2 or 2P, or confirmed but phase
     offset larger than the widened window/2. Report the TESS-best period; do not overwrite
     Gaia's value in the output table.
   - `undecided_*`: `no_phase_coverage`, `no_tic_match`, `ambiguous_match`,
     `tess_lc_too_noisy` (5-sigma depth limit above the Gaia depth). Undecided is a result,
     not a failure bucket to be re-labelled.
4. **Odd/even and secondary checks for confirmations:** split events by cycle parity; a
   depth difference > 3 sigma flags an EB at 2P mis-folded, or a secondary eclipse. Mandatory
   for transit candidates before anyone says "planet" outside the project.
5. **Outputs:** one machine-readable table (candidate, match, verdict, measured depth/duration,
   SNR, phase coverage fraction, sectors used) plus a folded-light-curve figure per candidate
   not `undecided_no_phase_coverage`, generated without human clicking (bulk triage first,
   eyes second). Every figure stamps period, epoch source, and verdict so it can be checked
   standalone.
6. **Acceptance tests before 25 Nov (on DR3):**
   - Run the whole checker on DR3 `vari_planetary_transit` (table exists — VERIFIED 2026-10-18
     on ESA TAP) matched to TESS: recover the published TESS-confirmed subset behaviour
     qualitatively (known planets confirm; catalogue false positives do not all confirm).
   - Run on a sample of DR3 EBs with TESS-known periods (e.g. from the pipeline's own EB
     rejections): verdict distribution sane, alias rate measured, no crashes on missing epochs.
   - Negative control: fold TESS light curves on deliberately wrong periods (P x 1.037);
     `confirmed` rate must be consistent with the pipeline's false-alarm rate, else the
     checker is finding systematics, not events.
   - Epoch-conversion unit test on one hand-checked DR3 EB with a TESS eclipse at known BJD.

## Data volumes and day-one compute (ESTIMATES, sized to be wrong safely)

- **Transit candidates: 1,007.** Trivial once matched: at most 1,007 TESS light curves already
  on local disk, plus epoch photometry for 1,007 Gaia sources if wanted (DataLink per-source
  FITS/VOTable; order 10-100 KB each for tens-hundreds of epochs -> under about 100 MB).
  Whole job is minutes of compute; the work is vetting output by eye, budget person-hours not
  CPU-hours.
- **EBs: 2,392,500 catalogue rows — filter before touching light curves.** Day-one sensible
  subsets, in order: (a) G < 13 and global_ranking high (bright, reliable — the DR3 analogue
  query in targets_v2 with P >= 365 d returned only **3** rows, VERIFIED 2026-10-18, so the
  long-period bright corner is tiny); (b) all periods, G < about 14, matched to locally cached
  TESS curves (upper bound = overlap of 2.39M EBs with the 4M-star TESS cache; expect
  1e5-1e6, ESTIMATE); (c) the rest only if (a)-(b) finish clean. Folding is cheap per star
  (existing pipeline already phase-folds for vetting), so 1e6 folds is a batch job, not a
  research project — but do not promise (b) on day one until the DR3 dress rehearsal times it.
- **Epoch photometry:** fetch only for candidates whose verdict needs the Gaia light curve
  itself (refutations someone might dispute, aliases). Never mirror DR4 epoch photometry
  (that is the 400 TB). Cache fetched per-source files by source_id.
- **Archive load on 2 Dec will be extreme** (CONJECTURE from every previous Gaia release):
  day-one queries must be few, indexed (source_id lists, `TOP`-smoke-tested first), retrying
  with backoff, and every result cached to disk before analysis so a TAP outage mid-run loses
  nothing already fetched. Prefer the bulk-download copy of the small candidate tables
  (transit candidates, EB summary) over repeated TAP joins if ESA publishes them as files.

## Pre-25-Nov checklist (archive closes 25 Nov 12:00 CET)

Do these against DR3; anything unfinished on 25 Nov cannot be rehearsed until release day.

- [ ] Re-read ESA release + DR4 pages; update the dates/contents table above if changed.
- [ ] TIC -> Gaia DR3 match table built (ID bridge + positional fallback above), match_method
      recorded, ambiguous fraction reported; spot-check 20 matches by eye (position, G vs Tmag).
- [ ] DataLink epoch-photometry fetcher written and exercised on DR3 sources (tens of sources,
      incl. a transit candidate and an EB); cached-by-source_id layout in place.
- [ ] Folding checker implemented to the spec above; all four DR3 acceptance tests run and
      their outputs (incl. the negative-control confirmed-rate) written down in the pipeline's
      repo, linked from this file once they exist.
- [ ] Epoch conversion (Gaia reference_time + phase -> BJD of primary event) unit-tested.
- [ ] Dress rehearsal: full DR3 transit-candidate run end-to-end from TAP query to verdict
      table, timed; extrapolate the EB-subset timing from a 10,000-row DR3 EB sample.
- [ ] `src/targets_v2.py --gaia-release dr4` query string reviewed (schema prefix is the only
      change); accept that column names may drift — the script fails loudly per source rather
      than guessing.
- [ ] Decide the credit/announcement route in advance (e.g. RNAAS note / telegram / direct to
      Gaia DPAC contact) so a day-one confirmation is not delayed by process debate. Any
      external claim waits for a person to inspect figures (this repo's standing rule).
- [ ] Freeze: no pipeline changes between 25 Nov and first DR4 verdicts except the schema
      prefix and bug fixes with a DR3 regression run.

## Day-one run order (2 Dec)

1. Smoke-test `gaiadr4` TAP schema (`vari_eclipsing_binary`, `vari_planetary_transit`,
   neighbourhood table name) with `TOP 1` queries; fix names only.
2. Fetch and cache the transit-candidate table and the EB summary columns needed for
   filtering (ids, period/frequency, epoch/phase, depth, duration, ranking, magnitudes).
3. Extend the match table with `gaia_dr4_id`; match candidates to TIC/local TESS cache.
4. Transit candidates through the checker -> verdict table + figures -> person inspects.
5. EB subset (a) then (b) through the checker; long-period `undecided_no_phase_coverage` EBs
   with DASCH-bright magnitudes are appended as DASCH targets (targets_v2 regeneration or a
   dated supplement — never hand-edits to the ranked CSV).
6. Write down what was NOT checked (unmatched fraction, uncovered-phase fraction) in the same
   artefact as the confirmations.

## What this plan does NOT do

- It does not implement the checker or the fetcher; both are specified, not built.
- It does not touch DR4 astrometric/spectroscopic orbits, Solar-System, or classification
  products; the niche is the TESS photometric cross-check nobody else is positioned to run at
  full scale (AREAS item 2).
- It makes no claim about any individual DR4 candidate before the data exist.
