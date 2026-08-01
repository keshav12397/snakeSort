configfile: "config.yaml"
import pandas as pd


bird = config['birdName']   # singular now, not a list


           # f"{workingDir}/rawData/{{bird}}/transfer_done", ##file rule will be {bird,stream,ksdir}

workingDir = config['workingDir']
rigDir = config['rigDir']


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
        all_catgt_targets


rule transferFiles:
    input:
        f"{rigDir}/{{bird}}/" ##TODO have to remount when i switch rig computers
    output:
        flag = f"{workingDir}/rawData/{{bird}}/transfer_done"
    log:
        f"{workingDir}/logs/{{bird}}_rsync.log"
    shell:
        '''
        set -euo pipefail
        mkdir -p {workingDir}/rawData/{wildcards.bird}
        rsync -avmh {input} {workingDir}/rawData/{wildcards.bird}/ 2>&1 | tee {log} #TODO change to rclone so its faster for many small files ? 
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
        'python makeMetaDF.py {params.rawdir}'

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
        
        sepShankFlag = "-sepShanks" if config['catGTparams'].get('sepShanks', True) else ""
    output:
       flag = f"{workingDir}/pass_1/{{bird}}/{{bird}}_{{run}}_pass1done"
    log:
        f"{workingDir}/logs/{{bird}}_{{run}}_pass1catgt.log"
    shell:
        '''
        set -euo pipefail
        mkdir {params.dest}
        IFS=$'\t' read -r gtlist chnexcl prb < <(
        python parse_catGTargs.py {input.metadf} --run {wildcards.run}) 

        bash {params.catgtpath} \\
        -dir={params.startdir} \\
        -run={wildcards.run} \\
        -prb=$prb \\
        -gtlist=$gtlist \\
        -chnexcl=$chnexcl \\
        -dest={params.dest} \\
        -apfilter={params.apfilter} \\
        -{params.filtmode} \\
        {params.sepShankFlag} \\
        -zerofillmax={params.zerofillmax} \\
        -ap -ni -prb_fld -pass1_force_ni_ob_bin -out_prb_fld \\
        > {log} 2>&1

        touch {output.flag}
        '''




rule makeSupercat:
    input:

# rule backupRawData:
#     ##nothing will the touch the raw data its now safe to back up to aws 

#     input:
#         "{workingDIR}/pass_supercat/{BIRD}/catgtdone"
#     output:
#         'awslog'
#     shell:
#         "aws sync {workingDIR}/rawData/{BIRD} XXXXX/NPXRECORDINGS/ "

# rule runMotionCorrection:
#     input:

# rule DoKilosort:

# rule DoDartsort:

# rule doBombcell:

