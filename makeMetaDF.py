import pandas as pd 
import numpy as np
from pathlib import Path
import os 
import re 
import argparse


def readMeta(binFullPath):
    metaName = binFullPath.stem + ".meta"
    metaPath = Path(binFullPath.parent / metaName)
    metaDict = {}
    if metaPath.exists():
        # print("meta file present")
        with metaPath.open() as f:
            mdatList = f.read().splitlines()
            # convert the list entries into key value pairs
            for m in mdatList:
                csList = m.split(sep='=')
                if csList[0][0] == '~':
                    currKey = csList[0][1:len(csList[0])]
                else:
                    currKey = csList[0]
                metaDict.update({currKey: csList[1]})
    else:
        print("no meta file")
    return(metaDict)

def get_GTS(filename):
    '''
    Given a typical imec ap bin get the g, t, and imec idx
    '''
    match = re.search(r'_g(\d+)_t(\d+)\.imec(\d+)', filename)
    if match is None:
        raise ValueError(f"Couldn't parse g/t/imec from filename: {filename}")

    g_val, t_val, imec_val = match.groups()
    return int(g_val), int(t_val), int(imec_val)
def get_GT_NI(filename):
    '''
    Given a typical nidqbin get the g and t index
    '''
    match = re.search(r'_g(\d+)_t(\d+)\.nidq.bin', filename)
    if match is None:
        raise ValueError(f"Couldn't parse g/t/imec from filename: {filename}")

    g_val, t_val  = match.groups()
    return int(g_val), int(t_val)
def get_day(name):
    '''
    Extract the day value after _d as an X.Y string, defaulting to X.0 if no decimal present
    '''
    match = re.search(r'_d(\d+)(?:\.(\d+))?', name)
    if match is None:
        raise ValueError(f"Couldn't parse day value from: {name}")

    whole = match.group(1)
    frac = match.group(2) if match.group(2) is not None else "0"
    return f"{whole}.{frac}"
def makeMetaDF(animalDir):
        
    requiredMetaKeys = ['fileSHA1', 'fileSizeBytes', 'fileTimeSecs','fileCreateTime','firstSample'] #bad files miss these fields idk
    minFileTimeSec=5 #skip files too short bc sync edges get messed up

    hold_immetas = []
    hold_nimetas = []
    for root,dir, files in os.walk(animalDir):
        root = Path(root)
        for file in files:
            if file.endswith('.ap.bin'):
                #in an imec path do the thing 
                meta = readMeta(root/file)

                #check for bad files 
                missing_keys = [k for k in requiredMetaKeys if k not in meta]
                if missing_keys :
                    print(f'Skipping {file}, missing meta key:{missing_keys}')
                    continue
                fileTimesec=float(meta['fileTimeSecs'])
                if fileTimesec<minFileTimeSec:
                    print(f'Skipping {file}, file too short')

                g_ind, t_ind, stream_ind  = get_GTS(file)
                stream_ind = f'imec{stream_ind}'
                toremove  = f'_g{g_ind}_t{t_ind}.{stream_ind}.ap.bin'

                run_id = file.replace(toremove,"")
                day_id = get_day(run_id)

                hold_immetas.append(
                    {
                    f'{stream_ind}_filename': file,
                    'run_id': run_id,
                    'day_id': day_id,
                    'g_ind': g_ind,
                    't_ind': t_ind,
                    f'{stream_ind}_fileCreateTime': meta['fileCreateTime'],
                    f'{stream_ind}_fileTimeSecs': meta['fileTimeSecs'],
                    f'{stream_ind}_firstSample':meta['firstSample']
                    })

            if file.endswith('nidq.bin'):
                meta = readMeta(root/file)
                #check for bad files 
                missing_keys = [k for k in requiredMetaKeys if k not in meta]
                if missing_keys :
                    print(f'Skipping {file}, missing meta key:{missing_keys}')
                    continue
                fileTimesec=float(meta['fileTimeSecs'])
                if fileTimesec<minFileTimeSec:
                    print(f'Skipping {file}, file too short')

                g_ind, t_ind  = get_GT_NI(file)
                stream_ind = "nidq"
                toremove =  f'_g{g_ind}_t{t_ind}.{stream_ind}.bin'
                run_id = file.replace(toremove,"")
                day_id = get_day(run_id)

                hold_nimetas.append(
                {
                f'{stream_ind}_filename': file,
                'run_id': run_id,
                'day_id': day_id,
                'g_ind': g_ind,
                't_ind': t_ind,
                f'{stream_ind}_fileCreateTime': meta['fileCreateTime'],
                f'{stream_ind}_fileTimeSecs': meta['fileTimeSecs'],
                f'{stream_ind}_firstSample':meta['firstSample']
                })

            
    ni_df = pd.DataFrame(hold_nimetas)
    if len(ni_df) == 0:
        raise ValueError(" no NIDQ files found!")
    im_df = pd.DataFrame(hold_immetas)
    if len(im_df) == 0:
        raise ValueError("no IMEC files found!")
    
    fullDf = ni_df.merge(im_df,on = ['run_id','day_id','g_ind','t_ind'])
    #now i need to split day ids to sort this shit 
    fullDf[["day_ind0","day_ind1"]] = fullDf.day_id.str.split('.',expand=True).astype(int)
    fullDf = fullDf.sort_values(by=["day_ind0", "day_ind1", "g_ind", "t_ind"]).reset_index(drop=True)
    timecols = [x for x in fullDf.columns if 'fileCreateTime' in x]
    fullDf[timecols] = fullDf[timecols].apply(pd.to_datetime)

    #check to make sure these are ascending in real time ! 
    for col in timecols:
        diffs = fullDf[col].diff().dropna()
        if not (diffs > pd.Timedelta(0)).all():
            bad_idx = diffs[diffs <= pd.Timedelta(0)].index
            raise ValueError(f"{col} is NOT strictly increasing at rows: {list(bad_idx)}")

    fullDf.to_csv(animalDir/'metaDF.csv')

def main():
    parser = argparse.ArgumentParser(description="Build metadata_df.csv for a bird's raw SGLX directory")
    parser.add_argument("animalDir", type=Path, help="Path to the bird's raw data directory")
    args = parser.parse_args()

    makeMetaDF(args.animalDir)

if __name__ == "__main__":
    main()