import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from kilosort.io import load_ops


def doTPrimeOneDir(ksDir: Path, imecDir: Path, niDir: Path, tprime_path: Path, fs: float = None):
    '''
    Given a kilosort output directory, convert spike times to seconds, get sync
    edges from the CatGT-generated imec/nidq xd files, and run TPrime to make
    the mic-aligned st_fixed.npy file.

    imecDir/niDir are passed explicitly (rather than derived from ksDir's
    parents) since snakeSort's kilosort output lives flat under pass_2/{bird}
    (ks_imec{stream}), not nested inside the CatGT probe folder.
    '''
    if fs is None:
        needOpsPath = ksDir / 'ops.npy'
        if needOpsPath.exists():
            ops = load_ops(ksDir / 'ops.npy')
            fs = ops['fs']
        else:
            raise FileNotFoundError("no sampling rate given and no ops.npy to read it from")

    st = np.load(ksDir / 'spike_times.npy')
    st_sec = st / fs
    st_sec = np.asarray(st_sec, dtype=np.float64)
    np.save(ksDir / 'st_sec.npy', st_sec)

    savefixedpath = ksDir / 'st_fixed.npy'
    if savefixedpath.exists():
        print(f'{ksDir.name} fixed st exists, deleting!')
        os.remove(savefixedpath)

    tprimestring = f"-events=5,{ksDir/'st_sec.npy'},{ksDir/'st_fixed.npy'}"

    frompathstream = [imecDir / x for x in os.listdir(imecDir) if '.ap.xd' in x]
    tostreampath = [niDir / x for x in os.listdir(niDir) if '.nidq.xd' in x]
    scoffsetpath = [niDir / x for x in os.listdir(niDir) if 'sc_offsets.txt' in x]

    if len(scoffsetpath) == 1:
        tprimeargs = {
            '-syncperiod=': '1.0',
            '-tostream=': str(tostreampath[0]),
            '-fromstream=': f'5,{frompathstream[0]}',
        }
    elif len(scoffsetpath) == 0:
        tprimeargs = {
            '-syncperiod=': '1.0',
            '-tostream=': str(tostreampath[0]),
            '-fromstream=': f'5,{frompathstream[0]}',
        }
    else:
        raise ValueError('something wrong with sc_offsets path')

    fullcommand = [str(tprime_path)] + [k + v for k, v in tprimeargs.items()] + [tprimestring]

    print("Running TPrime command:")
    print(" ".join(fullcommand))
    subprocess.run(fullcommand, check=True)

    if savefixedpath.exists():
        print('tprime worked!')
    else:
        print('fixed file not made!')


def combineAndScoreBombChunks(bombdir):
    '''
    Combines all ephys and quality metric dataframes and scores for how many chunks they pass and writes to a new
    pickle. Combines all waveforms too and writes to a new NPY. Use for input into NWB.
    '''
    chunkTimes = np.load(bombdir / 'chunkIDXs.npy')
    chunkOffsets = chunkTimes[:, 1]

    allchunks = np.sort([int(x.split('_')[-1]) for x in os.listdir(bombdir) if 'chunk_' in x])
    assert np.all(np.diff(allchunks) == 1), 'chunks not increasing by 1!'
    tempread = pd.read_parquet(bombdir / 'chunk_0/templates._bc_qMetrics.parquet')

    nUnits = len(tempread)
    nChunks = len(allchunks)

    paramsinfo = pd.read_parquet(bombdir / 'chunk_0/_bc_parameters._bc_qMetrics.parquet')
    nRawSpikes = int(paramsinfo['nRawSpikesToExtract'].iloc[0])
    spikeWidth = int(paramsinfo['spike_width'].iloc[0])
    nChans = int(paramsinfo['n_channels_rec'].iloc[0])

    holdWFarr = np.full((nChunks, nUnits, spikeWidth, nChans, nRawSpikes), np.nan, dtype=np.float64)

    holdallDFs = []
    for thing in allchunks:
        usedir = bombdir / f'chunk_{thing}'
        qm_par = usedir / 'templates._bc_qMetrics.parquet'
        qm_df = pd.read_parquet(qm_par)
        ephys_par = usedir / 'templates._bc_ephysProperties.parquet'
        ephys_df = pd.read_parquet(ephys_par)
        ephys_df.rename(columns={'unit_id': 'phy_clusterID'}, inplace=True)
        totaldf = pd.merge(
            qm_df,
            ephys_df,
            on="phy_clusterID",
            how="left",
        )
        holdallDFs.append(totaldf)
        allunits = totaldf['phy_clusterID'].tolist()

        for unitIDX, unit in enumerate(allunits):
            wfpath = usedir / f'RawWaveforms/Unit{unit}_RawSpikes.npy'
            if not wfpath.exists():
                continue
            else:
                wfArr = np.load(wfpath)
                nspikes = wfArr.shape[-1]
                holdWFarr[thing, unitIDX, :, :, :nspikes] = wfArr

    cattedDf = pd.concat(holdallDFs, axis=1, keys=[f'chunkOffset_{x}' for x in chunkOffsets])
    chunkCol = cattedDf.columns.get_level_values(0).unique()
    bcClass = np.column_stack([cattedDf[x]['bc_unitType'].values for x in chunkCol])
    pctViol = np.column_stack([cattedDf[x]['pct_violations'].values for x in chunkCol])
    boolPassArr = (bcClass == 'GOOD') & (pctViol < 1)
    passPct = boolPassArr.sum(axis=1) / nChunks

    diffThresh = [0.8, 0.9, 1.0]
    diffPasses = [np.where(passPct >= x)[0] for x in diffThresh]

    for t, p in zip(diffThresh, diffPasses):
        falseArr = np.zeros(nUnits)
        falseArr[p] = 1
        cattedDf[('PassPctThresh'), (t)] = falseArr

    cattedDf.to_pickle(bombdir / 'combinedBombOut.pkl')
    np.save(bombdir / 'combinedBombWFs.npy', holdWFarr)

    return cattedDf, holdWFarr


def fastGetUnitIDXes(arr):
    '''
    Fastest way to get indices of all unique values in array.
    '''
    order = np.argsort(arr)
    arr_sorted = arr[order]

    unique_vals, starts, counts = np.unique(
        arr_sorted,
        return_index=True,
        return_counts=True
    )

    idx_by_value = {
        val: order[start:start + count]
        for val, start, count in zip(unique_vals, starts, counts)
    }
    return idx_by_value
