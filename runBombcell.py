import argparse
from pathlib import Path

import bombcell as bc


def runBombcell(streamID, outDir):
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
    bc.run_bombcellCHUNKEDALL(ksDir, bombOut, bcParams, ephysParam, calcEphys=True)


def main():
    parser = argparse.ArgumentParser(description='given a stream ID, run bombcell QC on the kilosort output')
    parser.add_argument('--streamID', type=str, required=True)
    parser.add_argument('--outDir', type=Path, required=True)
    args = parser.parse_args()

    runBombcell(args.streamID, args.outDir)


if __name__ == "__main__":
    main()
