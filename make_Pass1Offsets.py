import argparse
from pathlib import Path

import pandas as pd

from readSGLX import readOffsetsTxt


def buildPass1Offsets(metadf_path, pass1dir):
    metadf_path = Path(metadf_path)
    pass1dir = Path(pass1dir)

    #already sorted by day_ind0, day_ind1, g_ind, t_ind - make_MetaDF.py
    metaDF = pd.read_csv(metadf_path, index_col=0)

    holdOffsets = []
    for run_id, run_group in metaDF.groupby('run_id', sort=False):
        for g_ind, g_group in run_group.groupby('g_ind', sort=True):
            offsetsPath = pass1dir / run_id / f'catgt_{run_id}_g{g_ind}' / f'{run_id}_g{g_ind}_ct_offsets.txt'
            gOffsets = readOffsetsTxt(offsetsPath)
            #ct_offsets.txt rows are in the same t_ind order CatGT concatenated
            #them in, which matches g_group's own t_ind order since metaDF is
            #pre-sorted by t_ind within each (run_id, g_ind)
            assert len(gOffsets) == len(g_group), (
                f"{offsetsPath} has {len(gOffsets)} rows, expected {len(g_group)} to match metaDF"
            )
            holdOffsets.append(gOffsets)

    offsetsFrame = pd.concat(holdOffsets, axis=0).reset_index(drop=True)
    metaDF = metaDF.reset_index(drop=True)

    assert len(offsetsFrame) == len(metaDF), "offsets/metaDF row count mismatch after concatenation"

    return pd.concat([metaDF, offsetsFrame], axis=1)


def main():
    parser = argparse.ArgumentParser(description='build pass1_offsets.csv by joining metaDF.csv with each run/g_ind\'s CatGT ct_offsets.txt')
    parser.add_argument('--metadf', type=Path, required=True)
    parser.add_argument('--pass1dir', type=Path, required=True)
    parser.add_argument('--outPath', type=Path, required=True)
    args = parser.parse_args()

    fullDf = buildPass1Offsets(args.metadf, args.pass1dir)
    fullDf.to_csv(args.outPath, index=False)


if __name__ == "__main__":
    main()
