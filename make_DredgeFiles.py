import argparse
import json
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless - this runs as a Snakemake job, not interactively
import matplotlib.pyplot as plt
import numpy as np
import spikeinterface.full as si
import spikeinterface.extractors as se
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from spikeinterface.sortingcomponents.peak_localization import localize_peaks
from spikeinterface.sortingcomponents.motion import estimate_motion


def loadBadChanIDs(badChansPath):
    #bad channels are precomputed once per bird/stream by getBadChansSI
    #(make_CatGT_args.py) and echoed to pass_2/{bird}/ - reuse that instead of
    #re-running channel detection here
    return json.loads(Path(badChansPath).read_text())["badIDs"]


def findImecDir(outDir, streamID):
    #mirrors getBadChansSI's own directory discovery in make_CatGT_args.py -
    #done here at actual run time (not as a Snakemake params function) since
    #supercat's output folder-naming isn't worth hardcoding, and this dir may not
    #exist yet when Snakemake evaluates params during DAG-building
    matches = [p for p in Path(outDir).rglob("*") if p.is_dir() and p.name.endswith(f"imec{streamID}")]
    return matches[0]


def getPeakFallOff(peaks, fs, peak_locs, bin_sec=3, n_depth_bins=5, ignore_um=150):
    peakTimes = peaks['sample_index'] / fs
    max_sec = peakTimes.max()
    t_bins = np.r_[np.arange(0, max_sec, bin_sec), max_sec]

    peakDepths = peak_locs['y']
    min_depth, max_depth = np.quantile(peakDepths, [0, 1])
    space_bins = np.linspace(min_depth + ignore_um, max_depth - ignore_um, n_depth_bins)
    binned, t_edges, d_edges = np.histogram2d(x=peakTimes, y=peakDepths, bins=[t_bins, space_bins])
    binned_z = (binned - binned.mean(0)) / binned.std(0)

    return binned_z, t_edges, d_edges


def calcPeaksSI(streamID, badChansPath, outDir):
    outDir = Path(outDir)
    imecDir = findImecDir(outDir, streamID)
    rec = se.read_spikeglx(imecDir, stream_id=f'imec{streamID}.ap')

    badChans = loadBadChanIDs(badChansPath)
    if len(badChans) > 0:
        rec = rec.remove_channels(badChans)
    noise_levels = si.get_noise_levels(rec, return_in_uV=False)
    print(rec)

    DETECT_KWARGS = {
        'peak_sign': 'neg',
        'detect_threshold': 8, #was 8
        'exclude_sweep_ms': 0.8, #was 0.8
        'radius_um': 80.0, #was 80
        'noise_levels': noise_levels,
        'method': 'locally_exclusive',
    }
    LOCALIZE_PEAKS_KWARGS = {
        'radius_um': 75.0,
        'max_distance_um': 150.0,
        'optimizer': 'minimize_with_log_penality',
        'enforce_decrease': True,
        'feature': 'ptp',
        'method': 'monopolar_triangulation',
    }
    ESTIMATE_MOTION_KWARGS = {
        'direction': 'y',
        'rigid': False,
        'win_shape': 'gaussian',
        'win_step_um': 100.0,
        'win_scale_um': 200.0,
        'win_margin_um': -100.0,
        'bin_s': 2,
        'time_horizon_s': 1000,
        'method': 'dredge_ap',
        'histogram_depth_smooth_um': 2,
        'histogram_time_smooth_s': 2,
        'max_disp_um': 50,
    }
    SI_JOB_KWARGS = dict(n_jobs=-1, chunk_duration="5s")

    peaks = detect_peaks(recording=rec, **DETECT_KWARGS, **SI_JOB_KWARGS)

    bad_idxes = [] ##clean up any weird peaks + im just gonna remove the top and bottom 5% by amplitude bc theyre prob noise
    for thing in ['sample_index', 'channel_index', 'amplitude']:
        out = np.where(~np.isfinite(peaks[thing]))[0]
        bad_idxes.extend(out)
        if thing == 'amplitude':
            finite_vals = peaks[thing][np.isfinite(peaks[thing])]
            minAmp, maxAmp = np.quantile(finite_vals, [0.05, 0.95])
            bad2 = np.where((peaks[thing] <= minAmp) | (peaks[thing] >= maxAmp))[0]
            bad_idxes.extend(bad2)

    peaks = np.delete(peaks, bad_idxes)
    np.save(outDir / f'peaks_imec{streamID}.npy', peaks)

    peak_locs = localize_peaks(recording=rec, peaks=peaks, **LOCALIZE_PEAKS_KWARGS, **SI_JOB_KWARGS)
    np.save(outDir / f'peak_locs_imec{streamID}.npy', peak_locs)

    fs = rec.sampling_frequency

    binnedPeaks, t_edges, d_edges = getPeakFallOff(peaks, fs, peak_locs)
    z_mean = binnedPeaks.mean(1)
    z_std = binnedPeaks.std(1)
    plt.plot(t_edges[1:], z_mean, linewidth=0.1)
    plt.fill_between(t_edges[1:], z_mean - z_std, z_mean + z_std, alpha=0.5)
    plt.xlim(t_edges[1], t_edges[-1])
    plt.xlabel('seconds')
    plt.title(f'imec{streamID}')
    plt.savefig(outDir / f'peakfalloff_imec{streamID}.png')
    plt.close()

    motion = estimate_motion(recording=rec, peaks=peaks, peak_locations=peak_locs, **ESTIMATE_MOTION_KWARGS)
    motion_path = outDir / f'imec{streamID}_full_motion.p'
    with open(motion_path, 'wb') as f:
        pickle.dump(motion, f)

    return motion_path


def main():
    parser = argparse.ArgumentParser(description='given a stream ID, detect+localize peaks w/ SI, and make dredge motion estimate')
    parser.add_argument('--streamID', type=str, required=True)
    parser.add_argument('--badChansPath', type=Path, required=True)
    parser.add_argument('--outDir', type=Path, required=True)
    args = parser.parse_args()

    calcPeaksSI(args.streamID, args.badChansPath, args.outDir)


if __name__ == "__main__":
    main()
