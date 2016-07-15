#!/bin/bash -l 
# Run full ENCODE3 pipeline
export DBDIR=$HOME/data/cHMM/db
export BINDIR=$HOME/cHMM/bin
export DATADIR=$DBDIR/files

# -- Vars -- 
export MARKS=TRUE

cd $DBDIR # everything into data directory.
mkdir -p $DBDIR/Rout $DBDIR/out

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

# Wait until done! 


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

        # TODO Submit all of this as a qsub command {{{
        # For each epitope + DNase + WCE:
        while read -r p
        do 
            eval $( echo $p | awk -F'\t' '{printf("epitope=\"%s\";",$1)}' )
            echo "-- ${cell} + ${epitope}"

            # Download:
            eval $( grep "${epitope},${cell}" $LDIR/${epitope}.csv | awk -F',' '{a = $3; split(a,b,"/"); printf("link=\"%s\"; id=%s;",$3,b[5])}' )

            # TODO handle replicates!
            BAMFILE=${CELL_DIR}/${id}_${cell}_${epitope}.bam
            if [[ ! -s ${BAMFILE} ]] 
            then 
                echo "Downloading ${BAMFILE}"
                wget ${link} -o ${BAMFILE}.log -O ${BAMFILE} # BAM File
                samtools index $BAMFILE # Index 
            fi

            # STEP 1 =================================================================
            source $BINDIR/code_ENCODE3_process_step1.sh $id $cell $epitope ${CELL_DIR}



        done < $DBDIR/epitopes
        # }}}

        # STEP 0 - Get data we will need (bam format)


        # STEP 2
        # source $BINDIR/submit_ENCODE3_process_step2.sh


    done < $LDIR/available_marks.tsv

        # STEP 3
        # source $BINDIR/submit_ENCODE3_process_step3.sh
done


