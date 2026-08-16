import argparse
import pickle
from pathlib import Path

import numpy as np
import spikeinterface.full as si
import spikeinterface.extractors as se
from spikeinterface.sortingcomponents.motion import interpolate_motion
from kilosort import run_kilosort
import kilosort.io as ksio

from make_DredgeFiles import loadBadChanIDs, findImecDir


def runKilosort4_KS(streamID, outDir):
    outDir = Path(outDir)
    imecDir = findImecDir(outDir, streamID)

    #produced by make_DredgeFiles.py's calcPeaksSI
    motion_path = outDir / f'imec{streamID}_full_motion.p'
    if not motion_path.exists():
        raise FileNotFoundError(
            f"Missing full_motion.p file: {motion_path}\n"
            "Run doDredge first or check that outDir is correct."
        )

    rec_start = se.read_spikeglx(imecDir, stream_id=f'imec{streamID}.ap')
    fs = rec_start.sampling_frequency

    with open(motion_path, 'rb') as f:
        motion = pickle.load(f)
    motiontimes = motion.temporal_bins_s[0]
    print(motiontimes.max())
    maxSamp = int((motiontimes - motiontimes.min()).max() * fs) ##looks rough but have to do bc the time vector doesnt start at 0
    print(maxSamp / fs)

    rec_trimmed = rec_start.frame_slice(None, maxSamp)

    #bad channels are precomputed once per bird/stream by getBadChansSI
    #(make_CatGT_args.py) and echoed to pass_2/{bird}/ - reuse that instead of
    #re-running channel detection here
    badChansPath = outDir / f'imec{streamID}_bad_channels.json'
    badChans = loadBadChanIDs(badChansPath)
    if len(badChans) > 0:
        rec_trimmed = rec_trimmed.remove_channels(badChans)

    print('Interpolating motion!')
    rec_corr = interpolate_motion(rec_trimmed, motion, spatial_interpolation_method='kriging', border_mode='remove_channels', dtype='float32')

    rec_final = si.scale(rec_corr, gain=200.0, dtype="int16")

    ksProbeDict = {
        'chanMap': np.arange(rec_final.get_num_channels()),
        'xc': rec_final.get_channel_locations()[:, 0],
        'yc': rec_final.get_channel_locations()[:, 1],
        'kcoords': rec_final.get_channel_groups(),
        'n_chan': rec_final.get_num_channels()
    }

    probe_path = outDir / f'imec{streamID}_probe.json'
    ksio.save_probe(ksProbeDict, probe_path)
    print(rec_final)
    print('Writing KS file to disk!')
    useBinName = f'imec{streamID}_postmotion.bin'
    print(useBinName)
    filename, N, c, s, fs, _ = ksio.spikeinterface_to_binary(
        rec_final, outDir, data_name=useBinName, dtype=np.int16,
        chunksize=120000, export_probe=False) ##exporting probe from SI bc the KS function doesnt give kcoords

    probe = ksio.load_probe(probe_path)
    settings = {'fs': fs, 'n_chan_bin': c, 'filename': filename, 'nblocks': 0, 'highpass_cutoff': 0.01, 'batch_size': 120000}

    ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = run_kilosort(
        settings=settings, probe=probe, clear_cache=False, do_CAR=False,
        results_dir=outDir / f'ks_imec{streamID}',
    )

    return outDir / f'ks_imec{streamID}'


def main():
    parser = argparse.ArgumentParser(description='given a stream ID, run kilosort4 on the motion-corrected recording')
    parser.add_argument('--streamID', type=str, required=True)
    parser.add_argument('--outDir', type=Path, required=True)
    args = parser.parse_args()

    runKilosort4_KS(args.streamID, args.outDir)


if __name__ == "__main__":
    main()
