# snakeSort

![action potential snake](action_potential_snake_topdown_v2.png)

A Snakemake pipeline for Neuropixels preprocessing: pull a bird's recordings off the rig, run CatGT, concatenate runs with supercat, estimate motion with dredge, and spike-sort with Kilosort4.

## Pipeline

```
transferFiles (rclone)
  -> makeMetaDf (checkpoint)
  -> runCatGT          [per run]
  -> runSupercat
      -> backupAndDeleteData (raw data -> S3 Deep Archive)
      -> makePass1Offsets (metaDF.csv + CatGT ct_offsets.txt -> pass1_offsets.csv)
      -> doDredge        [per stream]
          -> runKilosort [per stream]
              -> doBombcell [per stream]
                  -> makeNWB (once every stream's doBombcell AND makePass1Offsets are done)
```

Bad channels are detected once per bird/stream (cached as `imec{stream}_bad_channels.json`) and reused by every downstream step instead of being recomputed.

`makeNWB` packs electrodes, units, bombcell QC metrics, and mic/DAF/syllable timeseries into a single `.nwb` file per bird in `nwbParams.saveDir`.

## Requirements

Besides the Python dependencies below, the pipeline shells out to a few external tools that aren't pip-installable:

- **CatGT** — vendored in `CatGT-linux/`, path set via `catGTparams.CATGTPATH` in `config.yaml`.
- **TPrime** — a Linux build, path set via `nwbParams.tprimePath`. Used by `makeNWB` to align spike times to the mic timeline; only invoked if `st_fixed.npy` doesn't already exist for a stream.
- **rclone** — must have a remote configured (matching `rigDir` in `config.yaml`) that can reach the recording rig.
- **AWS CLI** — must be configured with credentials that can write to the S3 bucket in `~/.aws/addr`.

`bombcell` is pulled in as a pip dependency, but from a fork/branch (`bombcellParams` below), not the upstream package — see `pyproject.toml`.

## Install

```bash
pip install -e .
```

This installs Snakemake and the Python dependencies (spikeinterface, kilosort, pandas, etc.) and makes the pipeline's helper scripts (`make_CatGT_args`, `make_DredgeFiles`, `runKilosort`, ...) importable.

## Configure

Everything the pipeline needs lives in `config.yaml`, one bird per run.

### Top level

| field | meaning |
|---|---|
| `birdName` | bird to process this run |
| `workingDir` | where `rawData/`, `pass_1/`, `pass_2/` get written |
| `rigDir` | rclone remote:path for the rig's recording share |
| `awsAddrPath` | path to a file containing the S3 bucket address `backupAndDeleteData` syncs raw data to |

### `catGTparams`

| field | meaning |
|---|---|
| `CATGTPATH` | path to CatGT's `runit.sh` |
| `filtmode` / `apfilter` / `zerofillmax` | passed straight through to CatGT |
| `sepShanks` | split each probe into 4 per-shank streams instead of one stream per probe |

### `bombcellParams`

| field | meaning |
|---|---|
| `chunkSizeSec` | window size `doBombcell` processes the recording in; `-1` means one chunk covering the whole recording instead of chunking |

### `nwbParams`

Consumed by `makeNWB` — everything here either can't be inferred from the pipeline's own outputs (mic/DAF line numbers, animal metadata) or names an external tool/output location.

| field | meaning |
|---|---|
| `micLine` | NI analog/digital line index carrying the mic signal |
| `dafLine` | NI line index carrying the DAF/efference copy signal |
| `onlineDetectMode` | `pulse`, `port`, or `none` — how online syllable-detect labels were encoded on the digital lines |
| `onlineDetectLines` | digital line(s) carrying those labels (single line for `pulse`, multiple for `port`) |
| `saveDir` | where the finished `{bird}.nwb` file gets written |
| `tprimePath` | path to a Linux TPrime binary (see Requirements) |
| `animalInfo.sessionDescription` | free-text session label, e.g. `"Lav69_pCAF"` |
| `animalInfo.sex` | `M` or `F` |
| `animalInfo.dateOfBirth` | ISO date, e.g. `"2026-03-30"` |
| `animalInfo.description` | free-text, e.g. parents/genotype |

`subject_id` isn't configured separately — it's always taken from `birdName`.

## Run

```bash
snakemake -c<N>          # run with N cores
snakemake -n              # dry run - see what would execute
```

Logs for each rule land in `{workingDir}/logs/`.
