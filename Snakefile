configfile: "config.yaml"
import pandas as pd


bird = config['birdName'] 
workingDir = config['workingDir']
rigDir = config['rigDir']

def checkInputJSON(configfile):


def all_catgt_targets(wildcards):
    targets = []
    runs = getRunsForBird(bird)
    targets += expand(
        f"{workingDir}/pass_1/{{bird}}/{{bird}}_{{run}}_pass1done",
        bird=bird, run=runs
    )
    return targets

rule all:
    input:
        f"{workingDir}/pass_2/{bird}/awsdone"



rule transferFiles:
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
        rclone copy {params.src} {workingDir}/rawData/{wildcards.bird} --transfers 8 -v --stats 10s 2>&1 | tee {log}
        touch {output.flag}
        '''
checkpoint makeMetaDf:
    input:
        flag = f"{workingDir}/rawData/{{bird}}/transfer_done",
    output:
        metadf =  f"{workingDir}/rawData/{{bird}}/metaDF.csv"
    params:
        rawdir = f"{workingDir}/rawData/{{bird}}"
    shell:
        'python make_MetaDF.py {params.rawdir}'

def getRunsForBird(bird):
    #this gets the output as defined in the makeMetaDf function which is defined as a checkpoint- this is the path 
    checkpoint_output = checkpoints.makeMetaDf.get(bird=bird).output.metadf
    df = pd.read_csv(checkpoint_output)
    return df["run_id"].unique().tolist()


rule runCatGT:
    input:
       metadf = f"{workingDir}/rawData/{{bird}}/metaDF.csv"
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
        > {log} 2>&1
        cp {params.dest}/'CatGT.log' {params.dest}/{{bird}}_{{run}}_pass1CatGT.log
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

def get_stream_ids(wildcards):
    #mirrors the -sepShanks stream-numbering runCatGT does in its shell block (same
    #(p+1)*100 + shank formula), computed here at DAG-build time so it can feed
    #params.streamIDs without waiting on runCatGT's stream_ids.txt to exist.
    #this same list is what future per-stream rules (dredge/kilosort/etc) will fan out over.
    probes = get_probes_for_bird(wildcards.bird)
    if config['catGTparams'].get('sepShanks', True):
        streams = []
        for p in probes:
            base = (p + 1) * 100
            streams += [base, base + 1, base + 2, base + 3]
    else:
        streams = probes
    return ",".join(str(s) for s in streams)

rule runSupercat:
    input:
        flags = get_pass1_flags
    output:
        flag = f"{workingDir}/pass_2/{{bird}}/supercat_done"
    params:
        catgtpath = config['catGTparams']['CATGTPATH'],
        pass1dir = f"{workingDir}/pass_1/{{bird}}",
        pass2dir = f"{workingDir}/pass_2/{{bird}}",
        streamIDs = get_stream_ids,
        runcount = lambda wc: len(getRunsForBird(wc.bird))
    log:
        f"{workingDir}/logs/{{bird}}_supercat.log" #this log will not print anyhting thoi
    shell:
        '''
        set -euo pipefail
        mkdir -p {params.pass2dir}

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
                ) > {log} 2>&1

        else
            echo "Single run — copying pass_1 output directly" > {log}
            cp -r {params.pass1dir}/. {params.pass2dir}/
        fi

        #carry the (already cross-run-validated) stream list forward so downstream
        #per-stream rules (dredge/kilosort/etc) can read it straight from pass_2
        echo "{params.streamIDs}" > {params.pass2dir}/stream_ids.txt

        touch {output.flag}
        '''

rule backupAndDeleteData:
    ##nothing will the touch the raw data its now safe to transfer and delete pass 1 
    input:
        flag = f"{workingDir}/pass_2/{bird}/supercat_done" #"{workingDIR}/pass_supercat/{BIRD}/catgtdone"
    output:
        flag = f"{workingDir}/pass_2/{bird}/awsdone"
    log:
        f"/home/kbsuresh/.aws/logs/{bird}_awsTransferlog"
    shell:
        '''
        aws s3 sync {workingDir}/rawData/{bird}  $(cat /home/kbsuresh/.aws/addr)/NpxRawRecordings/{bird}/ --storage-class DEEP_ARCHIVE --size-only --no-progress |& tee -a {log}
        touch {output.flag}
        rm -r 
        '''
        #"rm -r {workingD}"

rule doDredge:
   input:

# rule DoKilosort:

# rule DoDartsort:

# rule doBombcell:

