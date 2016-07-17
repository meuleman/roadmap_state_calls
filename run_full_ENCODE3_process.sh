#!/bin/bash -l 
# Run full ENCODE3 pipeline
export DBDIR=$HOME/data/cHMM/db
export BINDIR=$HOME/cHMM/bin
export DATADIR=$DBDIR/files
export TMP=$DATADIR/tmp

# -- Vars -- 
export MARKS=TRUE

cd $DBDIR # everything into data directory.
mkdir -p $DBDIR/Rout $DBDIR/out
mkdir -p $DATADIR $TMP # make sure they exist

# Get experiment information: 
R CMD BATCH --no-save --no-restore $BINDIR/get_experiments.R $DBDIR/Rout/output_get_experiments.Rout

# Get file links: 
types=(released_hg19 released_GRCh38 all)
for type in ${types[@]}
do
    echo $type
    qsub -cwd -j y -b y -V -N get_file_links_${type} \
        -o out/get_file_links_${type}.out \
        -l mem_free=2G "R CMD BATCH --no-save --no-restore \"--args type='${type}' marks.only='${MARKS}'\" $BINDIR/get_file_links.R Rout/output_get_file_links_${type}.Rout"
done

# TODO Wait until done! 


for type in ${types[@]}
do 
    echo "Starting ChromHMM pipeline for data from ${type}"

    # Type specific directories: 
    LDIR=$DBDIR/file_links/${type}
    TYPE_DIR=$DATADIR/${type}
    mkdir -p $TYPE_DIR
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
        while read -r p
        do 
            eval $( echo $p | awk -F'\t' '{printf("epitope=\"%s\";",$1)}' )
            echo "-- ${cell} + ${epitope}"

            # TODO make sure it handles replicates correctly.
            while read -r repl
            do
                # STEP 0 + 1 -- Download data and also
                # Filter and remove duplicates for each individual file:

                eval $( echo $repl | awk -F',' '{a = $3; split(a,b,"/"); printf("link=\"%s\"; id=%s;",$3,b[5])}' )

                JOBNAME=step1_${cell}_${epitope}_${id}
                qsub -cwd -q long -l m_mem_free=25G -N $JOBNAME -j y -b y -V -r y $BINDIR/code_ENCODE3_process_step1.sh $id $cell $epitope $link ${CELL_DIR}

            done < grep "${epitope},${cell}" $LDIR/${epitope}.csv 

            # TODO Hold until all previous step1 jobs for the cell type are done

            # STEP 2 -- Pool replicates!
            # TODO find phantopmeak quals software in pkg. 
            JOBNAME=step2_${cell}_${epitope} 
            qsub -cwd -q long -l m_mem_free=25G -N $JOBNAME -j y -b y -V -r y $BINDIR/code_ENCODE3_process_step2.sh $cell $epitope ${CELL_DIR}

            # TODO Figure out if there are enough reads in the dataset!

        done < $DBDIR/epitopes

        # STEP 0 - Get data we will need (bam format)



    done < $LDIR/available_marks.tsv

    # STEP 3
    # source $BINDIR/submit_ENCODE3_process_step3.sh
done


