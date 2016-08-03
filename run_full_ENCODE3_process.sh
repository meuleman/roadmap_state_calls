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
# Older version in anshul's direcotry:
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
    LDIR=$DBDIR/file_links/${type}
    TYPE_DIR=$DATADIR/${type}
    TCALL_DIR=$CALLDIR/${type}
    mkdir -p $TYPE_DIR
    mkdir -p $TCALL_DIR
    cd $TYPE_DIR

    # Processing pipeline for each cell:
    while read -r q
    do
        eval $( echo $q | awk -F'\t' '{printf("cell=\"%s\";",$1)}' )
        echo ${cell}
        CELL_DIR=${TYPE_DIR}/${cell} # Keep all files at the $type/$cell level
        mkdir -p ${CELL_DIR}
        cd ${CELL_DIR}

        # For each epitope + DNase + WCE:
        step2_jobs="empty"
        while read -r p
        do 
            eval $( echo $p | awk -F'\t' '{printf("epitope=\"%s\";",$1)}' )
            echo "-- ${cell} + ${epitope}"

            # TODO ARRAY JOBS!
            step1_jobs="empty"
            while read -r repl
            do
                # STEP 0 & 1: Download data, filter, remove duplicates.
                eval $( echo $repl | awk -F',' '{a = $3; split(a,b,"/"); printf("link=\"%s\"; id=%s;",$3,b[5])}' )

                JOBNAME=step1_${cell}_${epitope}_${id}
                STEP1_FILE="${CELL_DIR}/${id}_${cell}_${epitope}.filt.nodup.srt.SE.map.tagAlign.gz"
                if [[ ! -s $STEP1_FILE ]]
                then
                    qsub -cwd -q long -l mfree=25G -N $JOBNAME -j y -b y -V -r y -o $DBDIR/out/$JOBNAME.out $BINDIR/code_ENCODE3_process_step1.sh $id $cell $epitope $link ${CELL_DIR}
                fi

                if [[ "${step1_jobs}" == "empty" ]] # Add jobs for holding list:
                then 
                    step1_jobs=${JOBNAME}
                else
                    step1_jobs=${step1_jobs},${JOBNAME}
                fi
            done < <(grep "${epitope},${cell}" $LDIR/${epitope}.csv)

            # STEP 2 -- Pool replicates!
            JOBNAME=step2_${cell}_${epitope} 
            STEP2_FILE="${CELL_DIR}/FINAL_${cell_type}_${epitope}.tagAlign.gz"
            if [[ ! -s $STEP2_FILE ]]
            then
                qsub -cwd -q long -l mfree=25G -hold_jid ${step1_jobs} -N $JOBNAME -o $DBDIR/out/$JOBNAME.out -j y -b y -V -r y $BINDIR/code_ENCODE3_process_step2.sh $cell $epitope ${CELL_DIR}
            fi

            if [[ "${step2_jobs}" == "empty" ]] 
            then 
                step2_jobs=${JOBNAME}
            else
                step2_jobs=${step2_jobs},${JOBNAME}
            fi

            # TODO Figure out if there are enough reads in the dataset!

        done < $DBDIR/epitopes # list of epitopes we are interested in

        # TODO hold, wait step 2 to finish
        CC_DIR=${TCALL_DIR}/${cell} # Calls directory
        mkdir -p ${CC_DIR}

        JOBNAME=step3_${cell}
        qsub -cwd -q long -l mfree=25G -hold_jid ${step2_jobs} -N $JOBNAME -o $DBDIR/out/$JOBNAME.out -j y -b y -V -r y $BINDIR/code_ENCODE3_process_step3.sh $cell ${CELL_DIR} ${CC_DIR}

    done < $LDIR/available_marks.tsv # list of available cell types for our marks 
done


