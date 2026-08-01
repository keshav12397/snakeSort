configfile: "config.yaml"

_birds_raw = config['birdList']
BIRDS = [_birds_raw] if isinstance(_birds_raw, str) else list(_birds_raw) ##if its one bird just stick it in a list idk so it doesnt break donwstream

           # f"{workingDir}/rawData/{{bird}}/transfer_done", ##file rule will be {bird,stream,ksdir}

workingDir = config['workingDir']
rigDir = config['rigDir']


rule all:
    input:
        expand(
            f"{workingDir}/rawData/{{bird}}/metaDF.csv",
            bird=BIRDS
        )


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


rule makeMetaDf:
    input:
        flag = f"{workingDir}/rawData/{{bird}}/transfer_done",
    output:
        metadf =  f"{workingDir}/rawData/{{bird}}/metaDF.csv"
    params:
        rawdir = f"{workingDir}/rawData/{{bird}}"
    shell:
        'python makeMetaDF.py {params.rawdir}'

rule runCatGT_andSupercat:
    input:
       metadf = f"{workingDIR}/rawData/{{bird}}/metadata_df.csv"
    output:
       flag = f"{workingDIR}/pass_supercat/{{bird}}/catgtdone"
    shell:
        '''
        set -euo pipefail
        'python -c catgtwrapper ,-catgt params, metadatadf'
        "touch  {workingDIR}/pass_supercat/{BIRD}/catgtdone"
        '''

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

