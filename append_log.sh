#!/usr/bin/env bash
# Append a rule's finished log file to the consolidated master log under a
# timestamped header section. Used because Snakemake gives every rule its own
# log file (and CatGT writes its own CatGT.log on top of that) - this is the
# one place that stitches them all into a single readable timeline.
#
# Usage: append_log.sh <master_log> <section_title> <source_log>
set -euo pipefail

master_log="$1"
title="$2"
source_log="$3"

mkdir -p "$(dirname "$master_log")"

# doDredge/runKilosort fan out per-stream and may append concurrently, so
# serialize writes with a lock instead of letting them interleave.
(
    flock -x 200
    {
        echo ""
        echo "===== ${title} -- $(date '+%Y-%m-%d %H:%M:%S') ====="
        cat "$source_log"
    } >> "$master_log"
) 200>"${master_log}.lock"
