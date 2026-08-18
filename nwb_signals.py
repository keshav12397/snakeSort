import numpy as np
from scipy.ndimage import uniform_filter1d

from readSGLX import readMeta, SampRate, makeMemMapRaw, GainCorrectNI, ExtractDigital


def load_NI(binFullPath, chan, t0, t1):
    '''
    Just a copy to convert the NI data to an array

    Parameters
    ---------
    binFullPath: Path object to the NI stream .bin file
    chan: (list) of channels to extract, microphone is usually [0], efference is usually 1
    Returns
    -------
    convData: (array) microphone signal in Volts
    '''
    meta = readMeta(binFullPath)
    firstSamp = int(meta['firstSample'])
    sRate = SampRate(meta)
    tStart = t0
    tEnd = t1
    s0 = int(sRate * tStart)
    s1 = int(sRate * tEnd)
    rawData = makeMemMapRaw(binFullPath, meta)
    selectData = rawData[chan, s0:s1 + 1]
    convData = GainCorrectNI(selectData, chan, meta)
    return convData, firstSamp, sRate


def rising_and_falling_edge(x, thr, s0=0):
    x = np.asarray(x)

    above = x >= thr

    rising = np.flatnonzero(~above[:-1] & above[1:]) + 1
    falling = np.flatnonzero(above[:-1] & ~above[1:]) + 1

    # if signal starts already above threshold
    if above[0]:
        rising = np.r_[0, rising]

    # if signal ends still above threshold
    if above[-1]:
        falling = np.r_[falling, x.size]

    rising = rising + s0
    falling = falling + s0 + 1

    if rising.size == 0 or falling.size == 0:
        return np.empty((0, 2), dtype=int)

    return np.column_stack((rising, falling))


def rising_edge_times_digital(x, fs, s0, return_mode='seconds'):
    """
    x: 1D array of bool / {0,1} samples
    fs: sampling rate (Hz)
    s0: optional start sample offset
    returns: times (seconds or samples) where signal rises (0->1)
    """
    rising_idx = np.flatnonzero((x[1:] == 1) & (x[:-1] == 0)) + 1

    if return_mode == 'samples':
        return rising_idx
    if return_mode == 'seconds':
        times = (s0 + rising_idx) / fs
        return times.astype(np.float64)
    else:
        raise ValueError("return_mode must be one of ['samples','seconds']")


def pulses_to_sylls(pulses, gap_before_ms, on_window_ms):
    usegap = int(gap_before_ms * 40000 / 1000)
    usewindow = int(on_window_ms * 40000 / 1000)
    diffed = np.diff(pulses, prepend=0)

    syllOn = pulses[np.where(diffed >= usegap)[0]]

    hold_labels = []
    hold_starts = []
    for x in syllOn:
        start = x
        end = start + usewindow
        kept = pulses[(pulses >= x) & (pulses <= end)]
        n_sylls = kept.size
        first_time = kept[0]
        hold_labels.append(n_sylls)
        hold_starts.append(first_time)

    return np.array(hold_labels), np.array(hold_starts)


def decode_port_labels(digArr, lines, settle_samples=2, min_gap_samples=20):
    """
    Decode multi-line digital port labels.

    Parameters
    ----------
    digArr : array
        Shape: n_lines x n_samples
    lines : list-like
        Digital lines in bit order. First line is LSB, second line is next bit, etc.
    settle_samples : int
        Window (in samples) after event onset, inclusive, to take the max port value
        over - guards against reading a fixed offset that lands after a short pulse
        has already dropped back to 0.
    min_gap_samples : int
        Minimum gap between separate port events.

    Returns
    -------
    labels : ndarray
        Decoded integer labels.
    times_samples : ndarray
        Event onset times in absolute samples.
    """
    lines = np.asarray(lines)
    bits = digArr[lines, :].astype(bool)

    # weights: line 0 = 1, line 1 = 2, line 2 = 4, ...
    weights = 2 ** np.arange(len(lines))

    # integer port value at every sample
    port_value = bits.T @ weights

    # candidate event starts: port goes from 0 to nonzero
    event_idx = np.flatnonzero((port_value[1:] != 0) & (port_value[:-1] == 0)) + 1

    if event_idx.size == 0:
        return (
            np.array([], dtype=int),
            np.array([], dtype=int),
        )

    # enforce minimum gap between events
    keep = np.r_[True, np.diff(event_idx) >= min_gap_samples]
    event_idx = event_idx[keep]

    # take the max port value in the settle window after onset, rather than a
    # single fixed-offset sample, so a pulse narrower than settle_samples doesn't
    # get read after it's already dropped back to 0
    window_end = np.minimum(event_idx + settle_samples + 1, port_value.size)
    labels = np.array([port_value[s:e].max() for s, e in zip(event_idx, window_end)]).astype(int)
    times_samples = event_idx

    return labels, times_samples


def getOnlineLabels(digArr, fs, s0, onlineDetectParams):
    label_mode = onlineDetectParams[0]

    if label_mode == 'pulse':
        useline = onlineDetectParams[1]
        usearr = digArr[useline, :].ravel()
        alledges = rising_edge_times_digital(usearr, fs, s0, return_mode='samples')
        syllLabels, syllStarts = pulses_to_sylls(alledges, 10, 5)
        syllStarts = syllStarts / fs

    elif label_mode == 'port':
        useline = onlineDetectParams[1]
        if len(useline) <= 1:
            raise ValueError('port labels need > 1 digital line')
        useline = np.sort(useline)
        syllLabels, syllStarts = decode_port_labels(digArr, useline)
        syllStarts = syllStarts / fs

    elif label_mode is None:
        syllLabels = np.array([], dtype=float)
        syllStarts = np.array([], dtype=float)

    else:
        raise ValueError("onlineLabels must be one of [pulse, port, None]")

    return syllLabels, syllStarts


def get_dig_traces(usepath, t0, t1, streamType):
    meta = readMeta(usepath)
    if t1 == -1:
        t1 = float(meta['fileTimeSecs'])
    srate = SampRate(meta)

    firstSamp = int(t0 * srate)
    lastSamp = int(t1 * srate)

    if streamType == 'nidq':
        avail_lines = meta['niXDChans1'].split(":")
        dLineList = [x for x in range(int(avail_lines[0]), int(avail_lines[-1]) + 1)]
    elif streamType == 'imec':
        dLineList = [x for x in np.arange(16)]
    else:
        raise ValueError('stream type is wrong')

    dw = 0
    rawData = makeMemMapRaw(usepath, meta)
    digArr = ExtractDigital(rawData, firstSamp, lastSamp - 1, dw, dLineList, meta)

    return digArr, dLineList, firstSamp, srate


def get_feedbacks(arr, s0, th=0.02):
    smoothened = uniform_filter1d(np.abs(arr), size=100, mode="nearest")
    feedback_arr = rising_and_falling_edge(smoothened, thr=th, s0=s0)  # 2d arr of [start,stop]
    if feedback_arr.size == 0:
        return (
            np.array([], dtype=int),
            np.array([], dtype=int),
            np.array([], dtype=float)
        )

    hold_vols = []
    for x in feedback_arr:
        sub = arr[x[0]:x[1]]
        sub_abs = np.abs(sub)
        sub_abs = np.maximum(sub_abs, np.finfo(float).eps)
        sub_dB = 10 * np.log10(sub_abs)
        hold_vols.append(np.mean(sub_dB))

    return feedback_arr[:, 0], feedback_arr[:, 1], np.array(hold_vols)


def makeMotifs(syllDF, motifSseq, maxGapSamp):
    n_sylls = len(motifSseq)
    m_start = np.where(syllDF.final_cluster == motifSseq[0])[0]  # idxes starting w the first syll in motif

    holdlabels = []
    holdstarts = []
    holdends = []
    hold_motifID = []
    motifID = 0
    for x in m_start:
        end = x + n_sylls
        if end <= len(syllDF):  # stop before reaching the end of DF to avoid error
            sub = syllDF.iloc[x:end, :]
            seq = sub.final_cluster.values
            if np.all(seq == motifSseq):  # found valid seq
                syll_starts = sub.start.values
                syll_ends = sub.end.values
                gaps = syll_starts[1:] - syll_ends[:-1] - 1

                if (gaps <= maxGapSamp).all():  # no intersyllable gaps > maxGapSamp - final valid motif
                    for k in range(n_sylls):
                        hold_motifID.append(motifID)
                        holdlabels.append(motifSseq[k])
                        holdstarts.append(syll_starts[k])
                        holdends.append(syll_ends[k])
                        if k < n_sylls - 1:
                            hold_motifID.append(motifID)
                            holdlabels.append(f"gap_{motifSseq[k]}_{motifSseq[k+1]}")
                            holdstarts.append(syll_ends[k])  # half open intervals
                            holdends.append(syll_starts[k + 1])
                    motifID += 1

    allstarts = np.array(holdstarts)
    sort_idx = np.argsort(allstarts)  # sort all starts ascending

    allstarts = allstarts[sort_idx]
    alllabels = np.array(holdlabels)[sort_idx]
    allends = np.array(holdends)[sort_idx]
    allmotifs = np.array(hold_motifID)[sort_idx]

    alldurs = allends - allstarts + 1
    med_durs = {
        lab: np.median(alldurs[alllabels == lab])
        for lab in np.unique(alllabels)
    }

    warp_sc = np.array([
        dur / med_durs[lab]
        for dur, lab in zip(alldurs, alllabels)
    ])

    return {'motif_ID': allmotifs, 'syll_ID': alllabels, 'starts': allstarts, 'stops': allends, 'warp_sc': warp_sc}


def makeAllSylls(syllDF, motifSeq):
    required_cols = {'start', 'end', 'final_cluster'}

    missing = required_cols - set(syllDF.columns)
    if missing:
        raise ValueError(f'SyllDF missing columns {missing}')

    if not np.issubdtype(syllDF["start"].dtype, np.integer):
        raise TypeError("syllDF['start'] must be sample indices, not seconds")

    if not np.issubdtype(syllDF["end"].dtype, np.integer):
        raise TypeError("syllDF['end'] must be sample indices, not seconds")

    holdstarts = []
    holdstops = []
    holdwarps = []
    holdlabels = []

    for x in motifSeq:
        sub = syllDF[syllDF.final_cluster == x]
        holdlabels.extend(sub.final_cluster.values)
        holdstarts.extend(sub.start.values)
        holdstops.extend(sub.end.values)
        durations = sub.end.values - sub.start.values
        scale_sc = durations / np.median(durations)
        holdwarps.extend(scale_sc)

    allstarts = np.array(holdstarts)
    sort_idx = np.argsort(allstarts)
    allstarts = allstarts[sort_idx]
    alllabels = np.array(holdlabels)[sort_idx]
    allstops = np.array(holdstops)[sort_idx]
    allwarps = np.array(holdwarps)[sort_idx]

    return {'syll_ID': alllabels, 'starts': allstarts, 'stops': allstops, 'warp_sc': allwarps}
