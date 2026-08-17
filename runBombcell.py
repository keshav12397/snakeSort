import argparse
from pathlib import Path

import bombcell as bc


def runBombcell(streamID, outDir, chunkSizeSec=600):
    outDir = Path(outDir)
    #produced by runKilosort.py's runKilosort4_KS
    ksDir = outDir / f'ks_imec{streamID}'
    if not ksDir.exists():
        raise FileNotFoundError(
            f"Missing kilosort output directory: {ksDir}\n"
            "Run runKilosort first or check that outDir is correct."
        )

    bombOut = ksDir / 'bombcell'
    bcParams, ephysParam = bc.makeBothParams(ksDir)
    bc.run_bombcellCHUNKEDALL(ksDir, bombOut, bcParams, ephysParam, chunkSizeSec=chunkSizeSec, calcEphys=True)


def main():
    parser = argparse.ArgumentParser(description='given a stream ID, run bombcell QC on the kilosort output')
    parser.add_argument('--streamID', type=str, required=True)
    parser.add_argument('--outDir', type=Path, required=True)
    parser.add_argument('--chunkSizeSec', type=int, default=600)
    args = parser.parse_args()

    runBombcell(args.streamID, args.outDir, args.chunkSizeSec)


if __name__ == "__main__":
    main()
