
import argparse
import json
import pandas as pd
from pathlib import Path
import numpy as np
import spikeinterface.full as si
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre



def getBadChansSI(runDir,s,nSlices = 1,sliceSecs = 5):
    #bad channels shouldn't vary across runs for the same bird/stream, so compute
    #once per bird and reuse the cached result on every later call (incl. later
    #dredge/kilosort steps) instead of re-running detect_bad_channels each time.
    badchans_path = runDir / f"imec{s}_bad_channels.json"
    if badchans_path.exists():
        cached = json.loads(badchans_path.read_text())
        return cached["badIDs"], cached["badidx"], cached["catGTBad"]

    allIMdirs =  [p for p in runDir.rglob("*") if p.is_dir() and p.name.endswith(f'imec{s}')]
    
    streamID = f'imec{s}.ap'
    listofrecs = []
    for thing in allIMdirs:
        try:
            listofrecs.append(se.read_spikeglx(folder_path=thing,stream_id=streamID))
        except:
            print('failed',thing)
    fullrec =  si.concatenate_recordings(listofrecs) #there must be a better way to do this instead of concatenating everything 
    time_vec = fullrec.get_times()
    kms = np.arange(time_vec.min(),time_vec.max()-sliceSecs,sliceSecs)
    getStarts = np.random.choice(kms,size  = min(nSlices,len(kms)), replace=False)
    subrec = si.concatenate_recordings([fullrec.time_slice(x,x+sliceSecs) for x in getStarts])
    subrec = spre.highpass_filter(subrec,freq_min=300.0)
    subrec = spre.phase_shift(subrec)
    badIDs,badidx = spre.detect_bad_channels(subrec,method = 'coherence+psd') ##spike interface channel IDs are formatted ie imec0.ap#AP
    delStr =  f'{streamID}#AP'
    badnumbers = sorted([int(x.replace(delStr,'')) for x in badIDs])
    ##format string as catgt expects 
    if len(badnumbers)>0:
        catGTBad = f"{{{s};{','.join(str(x) for x in badnumbers)}}}"
    else:
        catGTBad =" " #empty string

    #persist so later calls (this bird/stream, or later dredge/kilosort steps)
    #can read these back in instead of recomputing
    badchans_path.write_text(json.dumps({
        "badIDs": np.asarray(badIDs).tolist(),
        "badidx": np.asarray(badidx).tolist(),
        "catGTBad": catGTBad,
    }))

    return badIDs,badidx,catGTBad


def main():
    parser = argparse.ArgumentParser(description = 'given the concatenated metadata file and the run , return ALL (gtlist,chnexcl,prb)')

    parser.add_argument('metadf', type = Path)
    parser.add_argument("--run", type = str, required = True)
    args = parser.parse_args()


    df = pd.read_csv(args.metadf)
    run_df = df[df["run_id"] == args.run]

    parts = []
    for g_ind, group in run_df.groupby("g_ind"):
        t0 = group["t_ind"].min()
        t1 = group["t_ind"].max()
        parts.append(f"{{{int(g_ind)},{int(t0)},{int(t1)}}}")


    gtlist = "".join(parts)

    #look at meta columns to get possible stream names 
    possible_streams = set([int(x.split('_')[0].replace("imec","")) for x in run_df.columns if "imec" in x])
    possible_streams = sorted(possible_streams)

    #find bad channels per stream 
    holdbadchans = []
    for s in possible_streams:
        runDir = args.metadf.parent
        _,_,badChans = getBadChansSI(runDir, s)
        holdbadchans.append(badChans)
    chnexcl = "".join(holdbadchans)

    if len(possible_streams)> 1:
        prb = ",".join(possible_streams) #this type errors but sorted() makes it a list 
    else: 
        prb = possible_streams[0]

    print(f"{gtlist}\t{chnexcl}\t{prb}")
    ##load the first recording into spikeinterface to find bad chans 


if __name__ == "__main__":
    main()