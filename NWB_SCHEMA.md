# NWB file schema

One `{bird}.nwb` file per bird, written by `make_NWB.py` (`initNWB`) to `nwbParams.saveDir`.
Reflects the actual code in `make_NWB.py` / `nwb_signals.py` / `nwb_units.py` — if you change what gets packed, update this file too.

Field names below are exact — copy/paste them when reading the file (e.g. `nwbfile.units['goodCell'][:]`, `nwbfile.processing['Sync Edges']`).

## Top-level file metadata

| what | where | notes |
|---|---|---|
| session description | `nwbfile.session_description` | from `nwbParams.animalInfo.sessionDescription` |
| session start time | `nwbfile.session_start_time` | first `nidq_fileCreateTime` in the offsets table (below), i.e. the first raw file's creation timestamp |
| identifier | `nwbfile.identifier` | random UUID, no meaning |
| subject | `nwbfile.subject` | `subject_id` (= birdName), `sex`, `date_of_birth`, `description`, `species` (hardcoded `"Taeniopygia guttata"`) |
| device | `nwbfile.devices['NPX 2.0 4 shank']` | serial number + a `description` dict string with `imDatFx_hw/pn/sn`, `imDatHs_pn/sn`, `imDatPrb_pn/sn` pulled from the first stream's `.ap.meta` |

## Acquisition (`nwbfile.acquisition[...]`) — raw signals

| name | contents | shape / rate | present when |
|---|---|---|---|
| `Mic` | raw mic signal, one NI analog channel (`nwbParams.micLine`) | 1D, `rate` = NI sample rate, unit `V` | always |
| `NIDigitalLines` | every raw NI digital line, one column per line | `(nSamples, nLines)`, `rate` = NI sample rate, unit `n/a`; `comments` lists the line numbers in column order | only if the bird has NI digital lines at all |

Note: the DAF signal (`nwbParams.dafLine`) itself is **not** stored raw — only its derived onset/offset intervals (`DAFs`, below).

## Electrodes (`nwbfile.electrodes.to_dataframe()`)

One row per recording channel, across every stream in `imecStreams` (probes and, if `sepShanks`, per-shank streams — all disambiguated by synthetic stream IDs).

| column | meaning |
|---|---|
| `x`, `y` | channel position (µm), from `imec{s}_probe.json` |
| `shank` | shank index (`kCoords` from the KS probe map) |
| `stream` | which imec stream this channel belongs to, e.g. `"imec0"` — **use this to recover which probe/stream a unit came from** (see Units, below) |
| `chanMapID` | Kilosort channel-map channel number |
| `group` | electrode group object, named `imec{s}` |
| `id` | unique across all streams (not reused per-stream) |

## Units (`nwbfile.units.to_dataframe()`)

One row per Kilosort cluster that survived to bombcell scoring, across all streams. `id` is a fresh sequential counter across streams — **it is not the original Kilosort cluster ID**, and there is no column storing the original per-stream cluster ID. If you need to trace a unit back to its stream, go through `electrodes` → the referenced electrode rows' `stream` column (a unit's top channels are always on one stream, so any of them tells you which).

| column | meaning |
|---|---|
| `spike_times` | spike times in **seconds, TPrime-aligned to the mic/NI timeline** (`st_fixed.npy`) — use this for anything time-aligned to audio/behavior |
| `spikeIDXes` | misleadingly named: this is actually the **raw, un-aligned** spike times in samples (`spike_times.npy`, not TPrime-corrected). Kept for debugging drift/TPrime issues, not for analysis |
| `electrodes` | indices into the electrodes table for this unit's top 10 channels by template power (`nTopChans=10`), ordered strongest-first |
| `waveforms` | raw extracted spike waveforms on those same top channels, pooled across all bombcell chunks, NaN-only spike slots dropped. Shape per unit: `(nSpikesKept, spikeWidth, nTopChans)` |
| `TemplateWF` | Kilosort's own template waveform on those top channels. Shape per unit: `(spikeWidth, nTopChans)` |
| `MedianPosition_X`, `MedianPosition_Y` | median spike position (from `spike_positions.npy`) |
| `goodCell` | bool — bombcell pass fraction across chunks ≥ `keep_BCthresh` (0.8) |
| `cellType` | always `"unknown"` (placeholder, never actually classified) |
| *(everything else)* | bombcell quality/ephys metrics, **averaged across all bombcell chunks** for that unit (`RPV_window_index`, `estimatedTauR`, `fractionRPVs_estimatedTauR`, and any all-NaN columns are dropped) |

Units present in `spike_clusters.npy` but missing from the bombcell output for that stream are silently skipped (counted internally as `nMissedUnits`, only visible in the run log, not in the NWB file).

## Processing modules (`nwbfile.processing[...]`)

### `"Sync Edges"`
One `DynamicTable` per stream, named `sync_edges_imec{stream}`, columns `imec_sync_sec` / `nidq_sync_sec` — the raw sync-pulse edge timestamps (seconds) from the same `.ap.xd*`/`.nidq.xd*` files TPrime uses for `-fromstream=`/`-tostream=`. If the two streams have different pulse counts, the shorter column is NaN-padded at the end (via pandas' default-index alignment, not necessarily a positional pairing).

### `"metadata"`
One `DynamicTable`, `catGT_offsets` — this **is** `pass1_offsets.csv` (`metaDF.csv` joined column-wise with each run's CatGT `ct_offsets.txt`: `smp_imap0`, `smp_nidq`, `sec_imap0`, `sec_nidq`), merged with pass-2/supercat offsets if the bird went through supercat. One row per raw sub-file (every run × g-index × t-index × probe). This is the master lookup table for converting between a raw file's local time and the pooled/concatenated recording timeline — includes the original `imec{s}_filename` columns and `nidq_fileCreateTime`.

### `"NCAF things"`
Only populated from files found under `rawDataDir`; any of the three below is skipped (not added) if no matching files exist for this bird.

| table | contents |
|---|---|
| `nCAF params` (TimeIntervals) | one row per `params.txt` found, matched to a start/stop time in the offsets table by run name: `Targ_Syll`, `Delay_1`, `nCAF_Mode` (`push_up`/`push_down`), `Targ_unit` |
| `OSS params` (DynamicTable) | raw, unparsed dump of every `params.txt` found (`param`, `value`, `run` columns) — the source data behind `nCAF params` |
| `syll Log` (DynamicTable) | raw concatenated contents of every `syll_log.txt` found |

`oss_input` directories found under `rawDataDir` are **not** copied into the NWB — only their paths are recorded in this module's text `description`.

### `"Syllable Labels"`
Only populated if `{outDir}/finalSylls.csv` and `motif.p` exist (offline syllable segmentation — not currently produced by any Snakemake rule, so this module is empty for every bird run through the pipeline today).

| table | contents |
|---|---|
| `Sylls` (TimeIntervals) | every segmented syllable: `syll_ID`, `warp_sc` (duration vs. that syllable's median) |
| `Motifs` (TimeIntervals) | syllables + inter-syllable gaps that form a detected motif: `motif_ID`, `syll_ID` (gaps labeled `gap_{a}_{b}`), `warp_sc` |

## Top-level time intervals (`nwbfile.intervals[...]`)

Not inside a processing module — added directly via `add_time_intervals`. Either is omitted entirely if empty (an empty `TimeIntervals` can't be written to NWB).

| name | contents |
|---|---|
| `DAFs` | onset/offset of the DAF/white-noise efference-copy signal, plus `vol` (mean dB in that window) |
| `Online Labels` | decoded online syllable-detect events from `nwbParams.onlineDetectMode`/`onlineDetectLines` — `syll_ID` = decoded label, interval is a fixed 0.01s marker at detection time (not a real duration) |

## Quick lookup — "where is X?"

- **Raw mic audio** → `nwbfile.acquisition['Mic']`
- **Spike times aligned to mic/behavior** → `nwbfile.units['spike_times']`
- **Which probe/stream a unit is from** → `nwbfile.units['electrodes'][i]` → row(s) in `nwbfile.electrodes` → `stream` column
- **Bombcell QC metrics for a unit** → `nwbfile.units.to_dataframe()`, any column not in the fixed list above (chunk-averaged)
- **Raw NI digital lines** → `nwbfile.acquisition['NIDigitalLines']` (+ `comments` for line order)
- **CatGT/TPrime sync pulse edges** → `nwbfile.processing['Sync Edges']['sync_edges_imec{s}']`
- **File-level timing offsets (which raw file, when, cat-time offset)** → `nwbfile.processing['metadata']['catGT_offsets']`
- **NCAF live-sorting params/targets** → `nwbfile.processing['NCAF things']`
- **Online syllable-detect labels** → `nwbfile.intervals['Online Labels']`
- **DAF (efference copy) timing** → `nwbfile.intervals['DAFs']`
- **Offline syllable/motif segmentation** → `nwbfile.processing['Syllable Labels']` (empty until that step is wired into the pipeline)

## Known gaps / gotchas

- `units['spikeIDXes']` is named misleadingly — it's raw, TPrime-uncorrected spike times, not indices.
- No column stores a unit's original per-stream Kilosort cluster ID; only the electrode-table trace-back (above) recovers the stream.
- `units['cellType']` is always `"unknown"` — never actually populated.
- Bombcell metrics on `units` are per-unit chunk-averages, not per-chunk — the per-chunk data lives on disk (`{ks_imec{s}}/bombcell/combinedBombOut.pkl`) but isn't packed into the NWB file.
