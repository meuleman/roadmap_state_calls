#!/bin/bash -l 
# Run full ENCODE3 pipeline
export DBDIR=$HOME/data/cHMM/db
export BINDIR=$HOME/cHMM/bin

# -- Vars -- 
export MARKS=TRUE

cd $DBDIR # everything into data directory.
mkdir -p $DBDIR/Rout $DBDIR/out

# Get experiment information: 
R CMD BATCH --no-save --no-restore $BINDIR/get_experiments.R $DBDIR/Rout/output_get_experiments.Rout

# Get file links: 
types="released_hg19 released_hg38 preliminary proposed submitted all"
for type in ${types}; do
    qsub -cwd -j y -b y -V -N get_file_links_${type} \
        -o out/get_file_links_${type}.out \
        -l mem_free=2G "R CMD BATCH --no-save --no-restore \"--args type='${type}' marks.only='${MARKS}'\" \
        $BINDIR/get_file_links.R Rout/output_get_file_links_${type}.Rout"
done

# Wait until done!

# Find celltypes that have all data:  

# Run processing pipeline for each cell: 
# while read -r q
# do
    # eval $( echo $q | awk '{printf("CELL=%s;INFOFILE=%s;",$1,$2)}' )


    # STEP 1

    source $BINDIR/submit_ENCODE3_process_step1.sh


    # STEP 2
    source $BINDIR/submit_ENCODE3_process_step2.sh



    # STEP 3
    source $BINDIR/submit_ENCODE3_process_step3.sh

# done < cell_type_info


