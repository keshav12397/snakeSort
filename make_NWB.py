import argparse
import os
import pickle
from datetime import datetime
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from pynwb import NWBHDF5IO, NWBFile, TimeSeries, ProcessingModule
from pynwb.base import DynamicTable
from pynwb.epoch import TimeIntervals
from pynwb.file import Subject
from kilosort.io import load_probe

import readSGLX
import nwb_signals
import nwb_units
from make_DredgeFiles import findImecDir


def findPathFast(searchDir, target_name, kind='any'):
    if kind not in ['dir', 'file', 'any']:
        raise ValueError('kind must be one of "[dir,file,any]')

    stack = [Path(searchDir)]
    matches = []

    while stack:
        folder = stack.pop()
        try:
            with os.scandir(folder) as it:
                for entry in it:
                    if entry.is_dir():
                        if kind in ('dir', 'any') and entry.name == target_name:
                            matches.append(Path(entry.path))
                        stack.append(Path(entry.path))
                    elif entry.is_file():
                        if kind in ('file', 'any') and entry.name == target_name:
                            matches.append(Path(entry.path))
        except PermissionError:
            continue

    return matches


def get_param(outdf, key, default=0):
    vals = outdf.loc[outdf["param"] == key, "value"]
    if vals.empty:
        return default
    return vals.iloc[0]


def packDAFs(dafSignal, fs, micSeries):
    dafStart, dafEnd, dafVol = nwb_signals.get_feedbacks(dafSignal, 0)
    dafed = TimeIntervals(
        name="DAFs",
        description="LV generated white noise efference, bandpassed 1.5-8kHz",
    )
    dafed.add_column(name='vol', description="mean dB of DAF from the efference copy", data=np.array([], dtype=float))

    for n in tqdm(range(len(dafStart)), desc="Packing DAFs"):
        dafed.add_interval(start_time=dafStart[n] / fs, stop_time=dafEnd[n] / fs, vol=dafVol[n], timeseries=micSeries)
    return dafed


def packMotifs(syllDf, motifSeq, fs, micSeries):
    motifdict = nwb_signals.makeMotifs(syllDf, motifSeq, 0.3 * fs)

    motifed = TimeIntervals(
        name='Motifs', description='All syllables and gaps in a motif and the warping factor based on median duration')

    motifed.add_column('motif_ID', description='integer ID for motif')
    motifed.add_column('syll_ID', description='labeled syllable or gap')
    motifed.add_column('warp_sc', description='warping factor based on median syll/gap duration')

    allwarps = motifdict['warp_sc']
    allstarts = motifdict['starts']
    allstops = motifdict['stops']
    allmotifsIDs = motifdict['motif_ID']
    allsyllIDs = motifdict['syll_ID']

    print('total motifs found =', max(allmotifsIDs))
    print('1st and 99th percentile warping factors  = ', np.quantile(allwarps, [0.01, 0.99]))

    for n in tqdm(range(len(allwarps)), desc="Packing motifs"):
        motifed.add_interval(start_time=allstarts[n] / fs, stop_time=allstops[n] / fs, syll_ID=allsyllIDs[n], motif_ID=allmotifsIDs[n], warp_sc=allwarps[n], timeseries=micSeries)

    return motifed


def packAllSylls(syllDf, motifSeq, fs, micSeries):
    syllDict = nwb_signals.makeAllSylls(syllDf, motifSeq)

    allwarps = syllDict['warp_sc']
    allsyllIDs = syllDict['syll_ID']
    allstarts = syllDict['starts']
    allstops = syllDict['stops']
    n_sylls = len(allsyllIDs)
    n_segs = len(syllDf)
    pct_misssed = (n_segs - n_sylls) / n_segs

    print(f'*****\n{100 * pct_misssed} % of segments are unlabelled')
    vals, counts = np.unique(allsyllIDs, return_counts=True)
    print('Total Sylls = ')
    print(f"{'syll_ID':<20} {'Count':>10}")
    print("-" * 31)
    for v, c in zip(vals, counts):
        print(f"{str(v):<20} {c:>10}")

    sylled = TimeIntervals(name='Sylls', description='All syllables, even if not in a motif, and warping factor')
    sylled.add_column('syll_ID', description='labeled syllable ID')
    sylled.add_column('warp_sc', description='warping factor based on median syll/gap duration')

    for n in tqdm(range(len(allwarps)), desc="Packing syllables"):
        sylled.add_interval(start_time=allstarts[n] / fs, stop_time=allstops[n] / fs, syll_ID=allsyllIDs[n], warp_sc=allwarps[n], timeseries=micSeries)

    return sylled


def packOnlineLabels(syllLabels, syllStarts, micSeries, onlineDetectParams):
    vals, counts = np.unique(syllLabels, return_counts=True)
    print('Total Online detections = ')
    print(f"{'online_ID':<20} {'Count':>10}")
    print("-" * 31)
    for v, c in zip(vals, counts):
        print(f"{str(v):<20} {c:>10}")

    onlinesylled = TimeIntervals(name='Online Labels', description=f' {onlineDetectParams}')
    onlinesylled.add_column('syll_ID', description='labeled syllable ID', data=np.array([], dtype=int))

    for n in tqdm(range(len(syllLabels)), desc="Packing Online Labels"):
        onlinesylled.add_interval(start_time=syllStarts[n], stop_time=syllStarts[n] + 0.01, syll_ID=syllLabels[n], timeseries=micSeries)

    return onlinesylled


def packDeviceAndElectrodes(nwbfile: NWBFile, outDir: Path, imecStreams: list):
    #resolved per-stream (rather than one shared imecDir) - each synthetic
    #sepShanks stream ID already disambiguates across probes, so this handles
    #multi-probe birds for free
    firstImecDir = findImecDir(outDir, imecStreams[0])
    firstmetapath = next(firstImecDir.glob(f'*_tcat.imec{imecStreams[0]}.ap.meta'))
    usedevicemeta = readSGLX.readMeta(firstmetapath)
    keys_to_add = ['imDatFx_hw', 'imDatFx_pn', 'imDatFx_sn', 'imDatHs_pn', 'imDatHs_sn', 'imDatPrb_pn', 'imDatPrb_sn']
    addDict = {k: usedevicemeta[k] for k in keys_to_add}
    device_model = nwbfile.create_device_model(name='NPX 2.0 4 shank', manufacturer='Imec', model_number=addDict['imDatPrb_pn'])
    device = nwbfile.create_device(name='NPX 2.0 4 shank', description=str(addDict), serial_number=addDict['imDatPrb_sn'])
    nwbfile.add_electrode_column(name="shank", description="shank #")
    nwbfile.add_electrode_column(name='stream', description='imec stream ID')
    nwbfile.add_electrode_column(name='chanMapID', description='chan number from ks chan map')

    useid = 0
    for s in imecStreams:
        imecDir = findImecDir(outDir, s)
        metaPath = next(imecDir.glob(f'*_tcat.imec{s}.ap.meta'))
        meta = readSGLX.readMeta(metaPath)
        shank_id = str(s)[-1]
        electrode_group = nwbfile.create_electrode_group(name=f'imec{s}', description=f'imec stream split from shank {shank_id}', device=device, location='Area X')
        ref_type = meta['imChan0Ref'] + ' AND silver wire into brain'
        #probe json written flat by runKilosort.py, not nested under imecDir
        ksPrb = load_probe(outDir / f'imec{s}_probe.json')
        chanMap, xc, yc, kCoords, nChan = ksPrb.values()
        for i in range(nChan):
            nwbfile.add_electrode(
                id=useid,
                x=xc[i],
                y=yc[i],
                shank=int(kCoords[i]),
                group=electrode_group,
                reference=ref_type,
                stream=f'imec{s}',
                chanMapID=chanMap[i],
                location="Area X", enforce_unique_id=True
            )
            useid += 1


def packUnits(nwbfile: NWBFile, outDir: Path, imecStreams: list, tprime_path: Path, nTopChans=10, keep_BCthresh=0.8, DEBUG=False):
    nTotalUnits = 0
    nGoodUnits = 0
    nMissedUnits = 0
    nwb_id = 0
    for stream in imecStreams:
        ksDir = outDir / f'ks_imec{stream}'
        corrSTPath = ksDir / 'st_fixed.npy'
        if not corrSTPath.exists():  # run TPrime if corrected spike times don't exist
            imecDir = findImecDir(outDir, stream)
            niDir = imecDir.parent  # CatGT always writes the nidq tcat bin one level up from each probe's imec folder
            nwb_units.doTPrimeOneDir(ksDir, imecDir, niDir, tprime_path)

        corr_st = np.load(corrSTPath)  # TPrimed spike times aligned to mic stream
        raw_st = np.load(ksDir / 'spike_times.npy')
        spike_positions = np.load(ksDir / 'spike_positions.npy')  # (#nUnits,2)
        clu = np.load(ksDir / 'spike_clusters.npy')
        cluIdxDict = nwb_units.fastGetUnitIDXes(clu)

        templates = np.load(ksDir / 'templates.npy')  # (nUnits,spikewidth,nChans)
        chanMap = np.load(ksDir / 'channel_map.npy')
        chan_power = (templates ** 2).sum(axis=1)
        topChanIDX = np.argsort(chan_power, axis=-1)[:, -nTopChans:]
        topChanIDX = topChanIDX[:, ::-1]
        useTopChans = chanMap[topChanIDX]

        electTable = nwbfile.electrodes.to_dataframe()
        lookUpElectTable = electTable.reset_index().set_index(['stream', 'chanMapID'])['id']
        lookUpkeys = pd.MultiIndex.from_arrays([np.repeat(f'imec{stream}', useTopChans.size), useTopChans.ravel()], names=['stream', 'chanMapID'])
        unitElectTableIdx = lookUpElectTable.reindex(lookUpkeys).to_numpy().reshape(useTopChans.shape).astype(int)

        checkWFPath = ksDir / 'bombcell/combinedBombWFs.npy'
        checkQMsPath = ksDir / 'bombcell/combinedBombOut.pkl'

        if (not checkWFPath.exists()) or (not checkQMsPath.exists()):
            print('combining bombcell outputs for NWB')
            nwb_units.combineAndScoreBombChunks(ksDir / 'bombcell')

        allWFs = np.load(checkWFPath)  # (nChunks,nUnit,SpikeWidth, nChan,nSpikes)
        nWFchunks, nUnits, spikeWidth, nChan, nSpikes = allWFs.shape
        allQMs = pickle.load(open(checkQMsPath, 'rb'))
        checkUnit = allQMs[('PassPctThresh'), keep_BCthresh]
        nGoodUnits += checkUnit.sum()
        nTotalUnits += len(checkUnit)
        chunk_cols = [col for col in allQMs.columns if "chunk" in str(col[0])]
        chunk_df = allQMs.loc[:, chunk_cols].copy()
        chunk_numeric = chunk_df.apply(pd.to_numeric, errors="coerce")
        avg_df = chunk_numeric.T.groupby(level=1).mean().T
        avg_df = avg_df.set_index(avg_df['phy_clusterID'].astype(int))
        unUsedCols = ['RPV_window_index', 'estimatedTauR', 'fractionRPVs_estimatedTauR'] + avg_df.columns[avg_df.isna().all()].tolist()
        rdyDF = avg_df.drop(columns=unUsedCols)
        unitIDs = rdyDF.index.to_numpy()
        unit_to_row = {int(unit): i for i, unit in enumerate(unitIDs)}

        if DEBUG:
            print('clu max', clu.max())
            print('templates shape', templates.shape)
            print('bc WF shape', allWFs.shape)
            print('bc DF shape', allQMs.shape)

        assert templates.shape[0] == allWFs.shape[1], "templates and allWFs disagree on nUnits"
        assert len(checkUnit) == allWFs.shape[1], "checkUnit and allWFs disagree on nUnits"
        assert len(unitIDs) == allWFs.shape[1], "rdyDF and allWFs disagree on nUnits"

        columnsToAdd = ['goodCell', 'cellType'] + rdyDF.columns.tolist() + ['MedianPosition_X', 'MedianPosition_Y', 'TemplateWF', 'spikeIDXes']
        if stream == imecStreams[0]:
            for ok in columnsToAdd:
                if ok in ['TemplateWF', 'RawWFs', 'spikeIDXes']:
                    nwbfile.add_unit_column(name=ok, index=True, description='todo')
                else:
                    nwbfile.add_unit_column(name=ok, index=False, description='todo')

        for unit_id, spike_idx in tqdm(cluIdxDict.items(), desc='Adding Units'):
            if DEBUG:
                clu_units = set(map(int, cluIdxDict.keys()))
                bc_units = set(map(int, unitIDs))

                print("stream:", stream)
                print("n clu units:", len(clu_units))
                print("n bombcell/rdyDF units:", len(bc_units))
                print("overlap:", len(clu_units & bc_units))

            unit_id = int(unit_id)
            if unit_id not in unit_to_row:
                nMissedUnits += 1
                continue

            unit_row = unit_to_row[unit_id]

            keepSTCorr = np.sort(corr_st[spike_idx])
            keepSTRaw = np.sort(raw_st[spike_idx])
            keepSpikePos = spike_positions[spike_idx]
            medPos = np.nanmedian(keepSpikePos, axis=0)

            keepChanIDX = topChanIDX[unit_row, :]
            keepTemplateWFs = templates[unit_row, :, keepChanIDX]
            keepRawWfs = allWFs[:, unit_row, :, keepChanIDX, :].reshape(-1, spikeWidth, nTopChans)
            keepSpikeIDXes = ~np.isnan(keepRawWfs).all(axis=(1, 2))
            finalRawWFs = keepRawWfs[keepSpikeIDXes, :, :]

            metric_kwargs = {}
            for col in rdyDF.columns:
                val = rdyDF.loc[unit_id, col]
                metric_kwargs[col] = val

            nwbfile.add_unit(id=nwb_id,
                              electrodes=unitElectTableIdx[unit_row, :].astype(int),
                              spike_times=keepSTCorr,
                              cellType='unknown',
                              goodCell=bool(checkUnit[unit_row]),
                              MedianPosition_X=medPos[0], MedianPosition_Y=medPos[1],
                              TemplateWF=keepTemplateWFs, spikeIDXes=keepSTRaw, **metric_kwargs, waveforms=finalRawWFs)
            nwb_id += 1
    print(f'Done adding {nTotalUnits} units, you have {nGoodUnits} good Units and {nMissedUnits} were in clu array but not in the bombcell DFs')


def packNCAFthings(rawDataDir: Path, finalOffsetDF):
    '''
    Mess of a function to look for all relevant ncaf things and make a poorly formatted dynamic table to stick into
    processing lol. Also goes to a bit more effort to make a time interval with at minimum targ syll, delay 1, targ
    unit, and start and end times based on the offset table and the matching runName offset df data.
    '''
    thresholdmodemap = {0: 'push_up', 1: 'push_down'}

    paramName = 'params.txt'
    paramPaths = findPathFast(rawDataDir, paramName, kind='file')
    print(paramPaths)
    ncafparamTI: TimeIntervals = TimeIntervals(name='nCAF params', description='Output params used for Live Spike Sorting in this time interval')
    ncafparamTI.add_column('Targ_Syll', description='targeted Syllable index from the online detects')
    ncafparamTI.add_column('Delay_1', description='Delay in ms from syll onset to start counting spikes')
    ncafparamTI.add_column('nCAF_Mode', description='whether feedback is given below threshold(push up) or above threshold(push down)')
    ncafparamTI.add_column('Targ_unit', description='targeted template IDX from baseline spike sorting')

    finalParam = pd.DataFrame([])
    if len(paramPaths) > 0:
        print('Found nCAF param files at', paramPaths)
        holdparams = []
        for pp in paramPaths:
            outdf = pd.read_csv(pp, names=['param', 'value'], header=None)
            runOptions = [x for x in pp.parts if '_g' in x]
            if len(runOptions) == 0:
                runOptions = ['unknown']
            runName = runOptions[0]
            outdf['run'] = runName
            holdparams.append(outdf)
            useIDXes = finalOffsetDF.index[finalOffsetDF["file"].astype(str).str.contains(runName, regex=False)].tolist()
            if len(useIDXes) == 0:
                continue

            startIDX, endIDX = min(useIDXes), max(useIDXes)
            start_sec = float(finalOffsetDF.loc[startIDX, 'sec_nidq'])
            stop_sec = float(finalOffsetDF.loc[endIDX, 'sec_nidq'])

            usedelay1 = get_param(outdf, 'Delay 1')
            usetargsyll = get_param(outdf, 'Syllable Number List')
            rawthreshmode = get_param(outdf, 'Threshold mode')
            usethresholdmode = thresholdmodemap.get(int(rawthreshmode), 'unknown')
            useTargunit = get_param(outdf, 'Template Index List')

            ncafparamTI.add_interval(start_time=start_sec, stop_time=stop_sec, Targ_Syll=usetargsyll, Delay_1=usedelay1, nCAF_Mode=usethresholdmode, Targ_unit=useTargunit)

        if len(holdparams) > 0:
            finalParam = pd.concat(holdparams, axis=0)

    finalSyllLog = pd.DataFrame([])
    syllLogName = 'syll_log.txt'
    syllLogPaths = findPathFast(rawDataDir, syllLogName, kind='file')
    if len(syllLogPaths) > 0:
        print('found Syll Log paths at ', syllLogPaths)
        holdsylllogs = []
        for ss in syllLogPaths:
            onelog = pd.read_csv(ss, header=None)
            onelog.columns = [str(x) for x in onelog.columns]
            runOptions = [x for x in ss.parts if '_g' in x]
            if len(runOptions) == 0:
                runOptions = ['unknown']
            runName = runOptions[0]
            onelog['run'] = runName
            holdsylllogs.append(onelog)

        if len(holdsylllogs) > 0:
            finalSyllLog = pd.concat(holdsylllogs, axis=0)

    ossInputDir = 'oss_input'
    ossInputPaths = findPathFast(rawDataDir, ossInputDir, kind='dir')

    return ncafparamTI, finalParam, finalSyllLog, ossInputPaths


def initNWB(outDir: Path, animalInfo: dict, micLine: int, dafLine: int, onlineDetectParams: tuple, imecStreams: list, rawDataDir: Path, saveDir: Path, tprime_path: Path):
    #outDir doubles as the NWB-packing "masterDir" - it's pass_2/{bird}, and
    #its name is just the bird name in snakeSort's layout (no run/supercat
    #prefix to strip like the old Windows workflow's masterDir naming)
    outDir = Path(outDir)
    birdName = outDir.name

    niBinPath = next(outDir.rglob('*_tcat.nidq.bin'))
    offsetDfpath = outDir / 'pass1_offsets.csv'
    #only present for genuine multi-run supercat output (written directly at
    #outDir's top level by CatGT's -supercat mode) - exact filename suffix is
    #unverified since no multi-run bird has gone through this pipeline yet
    pass2offsetDfpath = next(outDir.glob('*offsets.txt'), None)

    offset1df = pd.read_csv(offsetDfpath)
    if pass2offsetDfpath is not None:
        offset2df = readSGLX.readOffsetsTxt(pass2offsetDfpath)
        finalOffsetDF = readSGLX.makeFinalOffsetDF(offset1df, offset2df)
    else:
        finalOffsetDF = offset1df

    expStart = pd.Timestamp(finalOffsetDF['nidq_fileCreateTime'].iloc[0]).to_pydatetime()

    nwbfile: NWBFile = NWBFile(
        session_description=animalInfo['session_description'],
        session_start_time=expStart,
        identifier=str(uuid4()),
    )

    if rawDataDir is not None:
        print('since you passed raw data dir im assuming you want to find nCAF params/inputs')

        nCAF_paramTI, fullParamDF, fullSyllLogDF, ossInputPaths = packNCAFthings(rawDataDir, finalOffsetDF)

        NCAF_mod = nwbfile.create_processing_module(
            name="NCAF things",
            description=f"params and inputs used for ncaf, ossinputs can be found = {ossInputPaths}"
        )

        NCAF_mod.add(nCAF_paramTI)
        if len(fullParamDF) > 0:
            paramTable = DynamicTable.from_dataframe(fullParamDF, name='OSS params')
            NCAF_mod.add(paramTable)
        if len(fullSyllLogDF) > 0:
            syllLogTable = DynamicTable.from_dataframe(fullSyllLogDF, name='syll Log')
            NCAF_mod.add(syllLogTable)

    _ = packDeviceAndElectrodes(nwbfile, outDir, imecStreams)
    _ = packUnits(nwbfile, outDir, imecStreams, tprime_path, keep_BCthresh=0.8, DEBUG=False)

    offsetTable = DynamicTable.from_dataframe(finalOffsetDF, name='catGT_offsets')
    metadata_mod: ProcessingModule = nwbfile.create_processing_module(
        name="metadata",
        description="catgt onset and offsets "
    )
    metadata_mod.add(offsetTable)

    subject: Subject = Subject(subject_id=animalInfo['subject_id'], sex=animalInfo['sex'], date_of_birth=animalInfo['date_of_birth'], description=animalInfo['description'], species='Taeniopygia  guttata')
    nwbfile.subject = subject

    maxSec = -1
    NI, _, ni_fs = nwb_signals.load_NI(niBinPath, [micLine, dafLine], 0, maxSec)
    mic = NI[micLine, :].ravel()
    mic_signal = TimeSeries(name='Mic', comments='raw mic signal ', data=mic, rate=ni_fs, unit='V', starting_time=0.0 / ni_fs)
    nwbfile.add_acquisition(mic_signal)

    print(f'*****\n Mic data added total duration = {int(NI.shape[-1]/ni_fs/60)} minutes')

    daf_signal = NI[dafLine, :].ravel()

    dafInter = packDAFs(daf_signal, ni_fs, mic_signal)
    _ = nwbfile.add_time_intervals(dafInter)
    print(f'*****\n {len(dafInter)} DAFs added')

    digArr, _, _ = nwb_signals.get_dig_traces(niBinPath, 0, maxSec, 'nidq')
    syllLabels, syllStarts = nwb_signals.getOnlineLabels(digArr, ni_fs, 0, onlineDetectParams)

    onlineLabelInter = packOnlineLabels(syllLabels, syllStarts, mic_signal, onlineDetectParams)
    _ = nwbfile.add_time_intervals(onlineLabelInter)

    syllabelMod = nwbfile.create_processing_module(
        name="Syllable Labels",
        description="Syllable and motif onset and offsets and labels"
    )

    expectedSyllPath = outDir / 'finalSylls.csv'
    if expectedSyllPath.exists():
        syllDf = pd.read_csv(expectedSyllPath)
        motif = pickle.load(open(outDir / 'motif.p', 'rb'))

        allSyllInter = packAllSylls(syllDf, motif, ni_fs, mic_signal)
        syllabelMod.add(allSyllInter)
        print('*****\n All syll data added!')

        motifInter = packMotifs(syllDf, motif, ni_fs, mic_signal)
        syllabelMod.add(motifInter)
        print('*****\n All motif data added!')

    with NWBHDF5IO(saveDir / f'{birdName}.nwb', 'w') as io:
        io.write(nwbfile)

    return nwbfile


def main():
    parser = argparse.ArgumentParser(description='pack a bird\'s pass_2 output (electrodes, units, bombcell QC, mic/DAF/syllable timeseries) into an NWB file')
    parser.add_argument('--outDir', type=Path, required=True)
    parser.add_argument('--rawDataDir', type=Path, required=True)
    parser.add_argument('--saveDir', type=Path, required=True)
    parser.add_argument('--tprimePath', type=Path, required=True)
    parser.add_argument('--streamIDs', type=str, required=True, help='comma-separated imec stream IDs')
    parser.add_argument('--micLine', type=int, required=True)
    parser.add_argument('--dafLine', type=int, required=True)
    parser.add_argument('--onlineDetectMode', choices=['pulse', 'port', 'none'], required=True)
    parser.add_argument('--onlineDetectLines', type=str, default='', help='comma-separated digital line indices')
    parser.add_argument('--subjectID', type=str, required=True)
    parser.add_argument('--sessionDescription', type=str, required=True)
    parser.add_argument('--sex', type=str, required=True)
    parser.add_argument('--dateOfBirth', type=str, required=True, help='ISO date, e.g. 2026-03-30')
    parser.add_argument('--description', type=str, required=True)
    args = parser.parse_args()

    imecStreams = [int(s) for s in args.streamIDs.split(",")]

    if args.onlineDetectMode == 'none':
        onlineDetectParams = (None, [])
    elif args.onlineDetectMode == 'pulse':
        onlineDetectParams = ('pulse', int(args.onlineDetectLines.split(",")[0]))
    else:
        onlineDetectParams = ('port', [int(x) for x in args.onlineDetectLines.split(",")])

    animalInfo = {
        'session_description': args.sessionDescription,
        'subject_id': args.subjectID,
        'sex': args.sex,
        'date_of_birth': datetime.fromisoformat(args.dateOfBirth),
        'description': args.description,
    }

    args.saveDir.mkdir(parents=True, exist_ok=True)

    initNWB(args.outDir, animalInfo, args.micLine, args.dafLine, onlineDetectParams, imecStreams, args.rawDataDir, args.saveDir, args.tprimePath)


if __name__ == "__main__":
    main()
