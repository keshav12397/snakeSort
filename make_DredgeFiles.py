import spikeinterface.extractors as se


DETECT_KWARGS= {
                                    'peak_sign': 'neg',
                                    'detect_threshold': 8, #was 8
                                    'exclude_sweep_ms':0.8, #was 0.8
                                    'radius_um': 80.0, #was 80
                                    'noise_levels': noise_levels,
                                    'method': 'locally_exclusive',
                                    }

def main():
    rec=se.read_spikeglx(imecDir,stream_id=f'imec{streamID}.ap')
    badChans = findDeadChansSI(rec)
    if len(badChans)>0:
        rec = rec.remove_channels(badChans)
        noise_levels = si.get_noise_levels(rec, return_in_uV=False)
        print(rec)