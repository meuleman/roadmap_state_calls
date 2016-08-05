#!/bin/bash
# SGE Array for STEP 2 -- Pool replicates and do QC if possible!
cell=$( sed "${SGE_TASK_ID}q;d" $LDIR/available_marks.tsv )

# TODO implement serial and parallel for epitopes.
echo "STEP 2 for ${cell}"
CELL_DIR=${TYPE_DIR}/${cell} 
cd ${CELL_DIR}

# In serial:
# For each epitope + DNase + WCE:
IFS=$'\t'
while read epitope
do 
    echo "- STEP2 ${cell}_${epitope}"
    STEP2_FILE="${CELL_DIR}/FINAL_${cell}_${epitope}.tagAlign.gz"
    echo $STEP2_FILE
    if [[ ! -s ${STEP2_FILE} ]]
    then
        source $BINDIR/code_ENCODE3_process_step2.sh $cell $epitope ${CELL_DIR}
    fi

    # TODO Figure out if there are enough reads in the dataset!

done < $DBDIR/epitopes # list of epitopes we are interested in

