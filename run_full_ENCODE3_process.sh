#!/bin/bash -l 
# Run full ENCODE3 pipeline
export DBDIR=$HOME/data/cHMM/db
export BINDIR=$HOME/cHMM/bin
export DATADIR=$DBDIR/files
export CALLDIR=$DBDIR/calls
export TMP=$DATADIR/tmp

# For step1 - filtering unique reads 
export SEQDIR="/broad/compbio/anshul/projects/encode/rawdata/sequence"
export UMAPDIR="/broad/compbio/anshul/projects/umap" 
export CHBIN="/broad/compbio/anshul/projects/encode/preprocessing/segmentations/chromhmm/scripts"
# Get CHMM here: http://compbio.mit.edu/ChromHMM/ChromHMM.zip
export CHMM=$HOME/data/software/ChromHMM/ChromHMM.jar
# Older version in anshul's directory:
# export CHMM="/broad/compbio/anshul/projects/encode/preprocessing/segmentations/chromhmm/ChromHMM/ChromHMM.jar"

# -- Vars -- 
export MARKS=TRUE
export NUMSTATES=18
export sleep_time=10

cd $DBDIR # everything into data directory.
mkdir -p $DBDIR/Rout $DBDIR/out
mkdir -p $DATADIR $TMP $CALLDIR # make sure main dir exist

# Get experiment information:
R CMD BATCH --no-save --no-restore $BINDIR/get_experiments.R $DBDIR/Rout/output_get_experiments.Rout

# Get file links: 
types=(released_hg19 released_GRCh38 all)
for type in ${types[@]}
do
    echo $type
    qsub -cwd -j y -b y -V -N get_file_links_${type} \
        -o $DBDIR/out/get_file_links_${type}.out \
        -l mfree=2G "R CMD BATCH --no-save --no-restore \"--args type='${type}' marks.only='${MARKS}'\" $BINDIR/get_file_links.R Rout/output_get_file_links_${type}.Rout"
done

# Wait until get_file_links jobs are done!
sleep 60 # wait for jobs to get on SGE

# Bash trap while running these jobs:
jobstatus=$(qstat -u $USERNAME | grep "get_file")
while [ -n "$jobstatus" ] # while jobs running
do
    sleep $sleep_time
    jobstatus=$(qstat -u $USERNAME | grep "get_file")
done

# CHMM processing pipeline: 
for type in ${types[@]}
do 
    echo "Starting ChromHMM pipeline for data from ${type}"

    # Type specific directories: 
    export LDIR=$DBDIR/file_links/${type}
    export TYPE_DIR=$DATADIR/${type}
    export TCALL_DIR=$CALLDIR/${type}
    mkdir -p $TYPE_DIR
    mkdir -p $TCALL_DIR
    cd $TYPE_DIR

    NCELL=$( wc -l $LDIR/available_marks.tsv | awk '{print $1}' )
    JOB1=step1_ENCODE3_array_$type
    qsub -cwd -q long -l mfree=25G -t 1-$NCELL -N $JOB1 -j y -b y -V -r y -o $DBDIR/out/${JOB1}_${SGE_TASK_ID}.out $BINDIR/submit_ENCODE3_process_step1.sh

    JOB2=step2_ENCODE3_array_$type
    qsub -cwd -q long -l mfree=25G -t 1-$NCELL -N $JOB2 -hold_jid_ad $JOB1 -j y -b y -V -r y -o $DBDIR/out/${JOB2}_${SGE_TASK_ID}.out $BINDIR/submit_ENCODE3_process_step2.sh

    JOB3=step3_ENCODE3_array_$type
    qsub -cwd -l mfree=25G -t 1-$NCELL -N $JOB3 -hold_jid_ad $JOB2 -j y -b y -V -r y -o $DBDIR/out/${JOB3}_${SGE_TASK_ID}.out $BINDIR/submit_ENCODE3_process_step3.sh

done

