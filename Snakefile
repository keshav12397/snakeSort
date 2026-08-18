configfile: "config.yaml"
import pandas as pd


bird = config['birdName']
workingDir = config['workingDir']
rigDir = config['rigDir']
awsAddrPath = config['awsAddrPath']
masterLog = f"{bird}_master.log"  # relative - lands wherever `snakemake` is invoked from

#def checkInputJSON(configfile):


def all_targets(wildcards):
    #wildcards is unused - rule all has none, bird is pinned by config - but Snakemake
    #input-functions must accept one. Deferred like this (rather than called directly
    #in rule all's input) so it runs during DAG-building, after get_stream_ids_for_bird
    #and the makeMetaDf checkpoint it depends on are actually defined.
    return [
        f"{workingDir}/pass_2/{bird}/awsdone",
        f"{workingDir}/pass_2/{bird}/nwbdone",
        *expand(
            f"{workingDir}/pass_2/{bird}/imec{{stream}}_bcdone",
            stream=get_stream_ids_for_bird(bird).split(",")
        )
    ]

rule all:
    input:
        all_targets

rule transferFiles:
    # Move files from my rig computer to starling, this is configured as an SMB drive in rclone 
    params:
        src = f"{rigDir}/{{bird}}/" ##TODO have to remount when i switch rig computers
    output:
        flag = f"{workingDir}/rawData/{{bird}}/transfer_done"
    log:
        f"{workingDir}/logs/{{bird}}_rclone.log"
    shell:
        '''
        set -euo pipefail
        mkdir -p {workingDir}/rawData/{wildcards.bird}
        exec 3>&1 4>&2
        exec > >(tee {log}) 2>&1
        set -x
        rclone copy {params.src} {workingDir}/rawData/{wildcards.bird} --transfers 8 -v --stats 10s
        set +x
        exec 1>&3 2>&4
        bash append_log.sh {masterLog} "transferFiles bird={wildcards.bird}" {log}
        touch {output.flag}
        '''
checkpoint makeMetaDf:
    #Now done transferring, parse all data files and decide the order based on 
    #the _d value, gX and tX. also gets the probe numbers 
    input:
        flag = f"{workingDir}/rawData/{{bird}}/transfer_done",
        script = "make_MetaDF.py",
    output:
        metadf =  f"{workingDir}/rawData/{{bird}}/metaDF.csv"
    params:
        rawdir = f"{workingDir}/rawData/{{bird}}"
    log:
        f"{workingDir}/logs/{{bird}}_makeMetaDf.log"
    shell:
        '''
        set -euo pipefail
        exec 3>&1 4>&2
        exec > {log} 2>&1
        set -x
        python make_MetaDF.py {params.rawdir}
        set +x
        exec 1>&3 2>&4
        bash append_log.sh {masterLog} "makeMetaDf bird={wildcards.bird}" {log}
        '''

def getRunsForBird(bird):
    #parse out metadf to make arguements for catgt -
    checkpoint_output = checkpoints.makeMetaDf.get(bird=bird).output.metadf
    df = pd.read_csv(checkpoint_output)
    return df["run_id"].unique().tolist()


rule runCatGT:
    #run pass 1 cat gt
    input:
       metadf = f"{workingDir}/rawData/{{bird}}/metaDF.csv",
       script = "make_CatGT_args.py",
    params:
        catgtpath = config['catGTparams']['CATGTPATH'],
        startdir = f"{workingDir}/rawData/{{bird}}",
        dest = f"{workingDir}/pass_1/{{bird}}/{{run}}",
        apfilter = config['catGTparams']['apfilter'],
        filtmode = config['catGTparams']['filtmode'],
        zerofillmax =  config['catGTparams']['zerofillmax'],

        sepShanks = "true" if config['catGTparams'].get('sepShanks', True) else "false"
    output:
       flag = f"{workingDir}/pass_1/{{bird}}/{{bird}}_{{run}}_pass1done"
    log:
        f"{workingDir}/logs/{{bird}}_{{run}}_pass1catgt.log"
    shell:
        '''
        set -euo pipefail
        mkdir -p {params.dest}
        exec 3>&1 4>&2
        exec > {log} 2>&1
        set -x

        IFS=$'\t' read -r gtlist chnexcl prb < <(
        python make_CatGT_args.py {input.metadf} --run {wildcards.run})

        #build one -sepShanks=ip,ip0,ip1,ip2,ip3 flag per probe, and track the resulting
        #output stream numbers - these become the effective "probe" streams for everything
        #downstream (supercat's -prb, and eventually per-stream dredge/kilosort/etc).
        #scheme: probe p's 4 shanks land on streams (p+1)*100 .. (p+1)*100+3, so probes never collide.
        sepShankArgs=""
        streamIDs=""
        if [ "{params.sepShanks}" = "true" ]; then
            IFS=',' read -ra probe_arr <<< "$prb"
            for p in "${{probe_arr[@]}}"; do
                base=$(( (p+1)*100 ))
                s0=$base; s1=$((base+1)); s2=$((base+2)); s3=$((base+3))
                sepShankArgs="$sepShankArgs -sepShanks=$p,$s0,$s1,$s2,$s3"
                streamIDs="${{streamIDs:+$streamIDs,}}$s0,$s1,$s2,$s3"
            done
        else
            streamIDs="$prb"
        fi
        echo "$streamIDs" > {params.dest}/stream_ids.txt

        (
        cd {params.dest}
        bash {params.catgtpath} \\
        -dir={params.startdir} \\
        -run={wildcards.run} \\
        -prb=$prb \\
        -gtlist=$gtlist \\
        -chnexcl=$chnexcl \\
        -dest={params.dest} \\
        -apfilter={params.apfilter} \\
        -{params.filtmode} \\
        $sepShankArgs \\
        -zerofillmax={params.zerofillmax} \\
        -ap -ni -prb_fld -pass1_force_ni_ob_bin -out_prb_fld \\
        )

        cp {params.dest}/'CatGT.log' {params.dest}/{wildcards.bird}_{wildcards.run}_pass1CatGT.log

        set +x
        exec 1>&3 2>&4
        bash append_log.sh {masterLog} "runCatGT bird={wildcards.bird} run={wildcards.run}" {log}
        bash append_log.sh {masterLog} "runCatGT bird={wildcards.bird} run={wildcards.run} (CatGT.log)" {params.dest}/{wildcards.bird}_{wildcards.run}_pass1CatGT.log
        touch {output.flag}
        '''

def get_pass1_flags(wildcards):
    runs = getRunsForBird(wildcards.bird)
    return expand(
        f"{workingDir}/pass_1/{{bird}}/{{bird}}_{{run}}_pass1done",
        bird=wildcards.bird, run=runs
    )

def get_probes_for_bird(bird):
    #same "which imec streams are present" logic as make_CatGT_args.py's possible_streams,
    #but read straight from the checkpoint's metaDF.csv - doesn't require runCatGT to have
    #actually executed yet, since it's just column names, not anything CatGT computes.
    checkpoint_output = checkpoints.makeMetaDf.get(bird=bird).output.metadf
    df = pd.read_csv(checkpoint_output)
    possible_streams = set(
        int(col.split('_')[0].replace("imec", "")) for col in df.columns if "imec" in col
    )
    return sorted(possible_streams)

def get_stream_ids_for_bird(bird):
    #mirrors the -sepShanks stream-numbering runCatGT does in its shell block (same
    #(p+1)*100 + shank formula), computed here at DAG-build time so it can feed
    #params.streamIDs without waiting on runCatGT's stream_ids.txt to exist.
    #this same list is what future per-stream rules (dredge/kilosort/etc) will fan out over.
    probes = get_probes_for_bird(bird)
    if config['catGTparams'].get('sepShanks', True):
        streams = []
        for p in probes:
            base = (p + 1) * 100
            streams += [base, base + 1, base + 2, base + 3]
    else:
        streams = probes
    return ",".join(str(s) for s in streams)

def get_stream_ids(wildcards):
    return get_stream_ids_for_bird(wildcards.bird)

rule runSupercat:
    #EITHER 1) run supercat if multiple runs, OR 2) just copy data to a pass2 directory, this simplifies all the downstream steps by keeping it in a consistent directory 
    input:
        flags = get_pass1_flags,
        script = "make_SuperCat_args.py",
    output:
        flag = f"{workingDir}/pass_2/{{bird}}/supercat_done"
    params:
        catgtpath = config['catGTparams']['CATGTPATH'],
        pass1dir = f"{workingDir}/pass_1/{{bird}}",
        pass2dir = f"{workingDir}/pass_2/{{bird}}",
        streamIDs = get_stream_ids,
        runcount = lambda wc: len(getRunsForBird(wc.bird))
    log:
        f"{workingDir}/logs/{{bird}}_supercat.log" #this will only print something if youre just transferring files, supercat log is written in catgt.log 
    shell:
        '''
        set -euo pipefail
        mkdir -p {params.pass2dir}
        exec 3>&1 4>&2
        exec > {log} 2>&1
        set -x

        if [ {params.runcount} -gt 1 ]; then
            supercat_arg=$(python make_SuperCat_args.py {params.pass1dir})
            (
            cd {params.pass2dir}
            bash {params.catgtpath} \\
                -supercat=$supercat_arg \\
                -ap -ni -prb_fld -out_prb_fld \\
                -prb={params.streamIDs} \\
                -supercat_trim_edges \\
                -dest={params.pass2dir} \\
                )

        else
            echo "Single run — copying pass_1 output directly"
            cp -r {params.pass1dir}/. {params.pass2dir}/
        fi

        #carry the (already cross-run-validated) stream list forward so downstream
        #per-stream rules (dredge/kilosort/etc) can read it straight from pass_2
        echo "{params.streamIDs}" > {params.pass2dir}/stream_ids.txt

        set +x
        exec 1>&3 2>&4
        bash append_log.sh {masterLog} "runSupercat bird={wildcards.bird}" {log}
        touch {output.flag}
        '''

rule backupAndDeleteData:
    #nothing will the touch the raw data its now safe to back up to AWS and delete pass 1 to save space on the SSD. 
    input:
        flag = f"{workingDir}/pass_2/{bird}/supercat_done" 
    output:
        flag = f"{workingDir}/pass_2/{bird}/awsdone"
    log:
        f"/home/kbsuresh/.aws/logs/{bird}_awsTransferlog"
    shell:
        '''
        set -euo pipefail
        exec 3>&1 4>&2
        exec > >(tee -a {log}) 2>&1
        set -x
        aws s3 sync {workingDir}/rawData/{bird}  $(cat {awsAddrPath})/NpxRawRecordings/{bird}/ --storage-class DEEP_ARCHIVE --size-only --no-progress
        set +x
        exec 1>&3 2>&4
        bash append_log.sh {masterLog} "backupAndDeleteData bird={bird}" {log}
        touch {output.flag}
        '''
        #"rm -r {workingD}" #TODO add back in delete

def get_bad_chans_path(wildcards):
    #getBadChansSI (make_CatGT_args.py) echoes its cached bad-channels json here
    return f"{workingDir}/pass_2/{wildcards.bird}/imec{wildcards.stream}_bad_channels.json"

def get_dredge_flags(wildcards):
    streams = get_stream_ids(wildcards).split(",")
    return expand(
        f"{workingDir}/pass_2/{{bird}}/imec{{stream}}_full_motion.p",
        bird=wildcards.bird, stream=streams
    )

rule doDredge:
    #load data into spikeinterface, find+localize peaks and save dredge motion estimate to pickle.
    input:
        flag = f"{workingDir}/pass_2/{{bird}}/supercat_done",
        script = "make_DredgeFiles.py",
    params:
        badChansPath = get_bad_chans_path,
        outDir = f"{workingDir}/pass_2/{{bird}}",
    output:
        motion = f"{workingDir}/pass_2/{{bird}}/imec{{stream}}_full_motion.p"
    log:
        f"{workingDir}/logs/{{bird}}_imec{{stream}}_dredge.log"
    shell:
        '''
        set -euo pipefail
        exec 3>&1 4>&2
        exec > {log} 2>&1
        set -x
        python make_DredgeFiles.py --streamID {wildcards.stream} --badChansPath {params.badChansPath} --outDir {params.outDir}
        set +x
        exec 1>&3 2>&4
        bash append_log.sh {masterLog} "doDredge bird={wildcards.bird} stream={wildcards.stream}" {log}
        '''

def get_kilosort_flags(wildcards):
    streams = get_stream_ids(wildcards).split(",")
    return expand(
        f"{workingDir}/pass_2/{{bird}}/imec{{stream}}_ksdone",
        bird=wildcards.bird, stream=streams
    )

rule runKilosort:
    #now that I have my dredge motion estimate i can actually kilosort 
    input:
        motion = f"{workingDir}/pass_2/{{bird}}/imec{{stream}}_full_motion.p",
        #runKilosort.py imports loadBadChanIDs from make_DredgeFiles.py, so a
        #change to either script should count as staleness here too
        scripts = ["runKilosort.py", "make_DredgeFiles.py"],
    params:
        outDir = f"{workingDir}/pass_2/{{bird}}",
    output:
        flag = f"{workingDir}/pass_2/{{bird}}/imec{{stream}}_ksdone"
    log:
        f"{workingDir}/logs/{{bird}}_imec{{stream}}_kilosort.log"
    shell:
        '''
        set -euo pipefail
        exec 3>&1 4>&2
        exec > {log} 2>&1
        set -x
        python runKilosort.py --streamID {wildcards.stream} --outDir {params.outDir}
        set +x
        exec 1>&3 2>&4
        bash append_log.sh {masterLog} "runKilosort bird={wildcards.bird} stream={wildcards.stream}" {log}
        touch {output.flag}
        '''

def get_bombcell_flags(wildcards):
    streams = get_stream_ids(wildcards).split(",")
    return expand(
        f"{workingDir}/pass_2/{{bird}}/imec{{stream}}_bcdone",
        bird=wildcards.bird, stream=streams
    )

rule doBombcell:
    #run bombcell QC on the kilosort output for this stream
    input:
        flag = f"{workingDir}/pass_2/{{bird}}/imec{{stream}}_ksdone",
        script = "runBombcell.py",
    params:
        outDir = f"{workingDir}/pass_2/{{bird}}",
        chunkSizeSec = config['bombcellParams']['chunkSizeSec'],
    output:
        flag = f"{workingDir}/pass_2/{{bird}}/imec{{stream}}_bcdone"
    log:
        f"{workingDir}/logs/{{bird}}_imec{{stream}}_bombcell.log"
    shell:
        '''
        set -euo pipefail
        exec 3>&1 4>&2
        exec > {log} 2>&1
        set -x
        python runBombcell.py --streamID {wildcards.stream} --outDir {params.outDir} --chunkSizeSec {params.chunkSizeSec}
        set +x
        exec 1>&3 2>&4
        bash append_log.sh {masterLog} "doBombcell bird={wildcards.bird} stream={wildcards.stream}" {log}
        touch {output.flag}
        '''

rule makePass1Offsets:
    #join metaDF.csv with each run/g_ind's CatGT ct_offsets.txt into one pass1_offsets.csv - the piece initNWB needs to align pass-1 time bases to the mic timeline
    input:
        flags = get_pass1_flags,
        script = "make_Pass1Offsets.py",
    params:
        metadf = f"{workingDir}/rawData/{{bird}}/metaDF.csv",
        pass1dir = f"{workingDir}/pass_1/{{bird}}",
    output:
        f"{workingDir}/pass_2/{{bird}}/pass1_offsets.csv"
    log:
        f"{workingDir}/logs/{{bird}}_pass1Offsets.log"
    shell:
        '''
        set -euo pipefail
        exec 3>&1 4>&2
        exec > {log} 2>&1
        set -x
        python make_Pass1Offsets.py --metadf {params.metadf} --pass1dir {params.pass1dir} --outPath {output}
        set +x
        exec 1>&3 2>&4
        bash append_log.sh {masterLog} "makePass1Offsets bird={wildcards.bird}" {log}
        '''

rule makeNWB:
    #pack electrodes, units, bombcell QC, and mic/DAF/syllable timeseries into one NWB file per bird
    input:
        bcflags = get_bombcell_flags,
        offsets = f"{workingDir}/pass_2/{{bird}}/pass1_offsets.csv",
        scripts = ["make_NWB.py", "nwb_signals.py", "nwb_units.py", "make_DredgeFiles.py", "readSGLX.py"],
    params:
        outDir = f"{workingDir}/pass_2/{{bird}}",
        rawDataDir = f"{workingDir}/rawData/{{bird}}",
        streamIDs = get_stream_ids,
        micLine = config['nwbParams']['micLine'],
        dafLine = config['nwbParams']['dafLine'],
        onlineDetectMode = config['nwbParams']['onlineDetectMode'],
        onlineDetectLines = ",".join(str(x) for x in config['nwbParams'].get('onlineDetectLines', [])),
        saveDir = config['nwbParams']['saveDir'],
        tprimePath = config['nwbParams']['tprimePath'],
        sessionDescription = config['nwbParams']['animalInfo']['sessionDescription'],
        sex = config['nwbParams']['animalInfo']['sex'],
        dateOfBirth = config['nwbParams']['animalInfo']['dateOfBirth'],
        description = config['nwbParams']['animalInfo']['description'],
    output:
        flag = f"{workingDir}/pass_2/{{bird}}/nwbdone"
    log:
        f"{workingDir}/logs/{{bird}}_makeNWB.log"
    shell:
        '''
        set -euo pipefail
        exec 3>&1 4>&2
        exec > {log} 2>&1
        set -x
        python make_NWB.py --outDir {params.outDir} --rawDataDir {params.rawDataDir} --saveDir {params.saveDir} \\
            --tprimePath {params.tprimePath} --streamIDs {params.streamIDs} --micLine {params.micLine} --dafLine {params.dafLine} \\
            --onlineDetectMode {params.onlineDetectMode} --onlineDetectLines "{params.onlineDetectLines}" \\
            --subjectID {wildcards.bird} --sessionDescription "{params.sessionDescription}" --sex {params.sex} \\
            --dateOfBirth {params.dateOfBirth} --description "{params.description}"
        set +x
        exec 1>&3 2>&4
        bash append_log.sh {masterLog} "makeNWB bird={wildcards.bird}" {log}
        touch {output.flag}
        '''
