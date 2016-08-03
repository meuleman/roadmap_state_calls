#!/bin/bash
# SGE Array for STEP 0 & 1: Download data, filter, remove duplicates.
cell=$( sed "${SGE_TASK_ID}q;d" $LDIR/available_marks.tsv )
cell2=$(echo ${cell} | awk '{gsub(" ","_",$0); print $0}')

# TODO implement serial and parallel for epitopes.
echo "STEP 1 for ${cell}"
CELL_DIR=${TYPE_DIR}/${cell} # Keep all files at the $type/$cell level
mkdir -p ${CELL_DIR}
cd ${CELL_DIR}

# In serial:
# For each epitope + DNase + WCE:
IFS=$'\t'
while read epitope
do 
    echo "\n- STEP 1 for ${cell} + ${epitope}"

    # Read replicates:
    while read -r repl
    do
        eval $( echo $repl | awk -F',' '{a = $3; split(a,b,"/"); printf("link=\"%s\"; id=%s;",$3,b[5])}' )

        echo "-- Running ${cell}_${epitope}_${id}"
        STEP1_FILE="${CELL_DIR}/${id}_${cell}_${epitope}.filt.nodup.srt.SE.map.tagAlign.gz"
        if [[ ! -s ${STEP1_FILE} ]]
        then
            source $BINDIR/code_ENCODE3_process_step1.sh $id $cell $epitope $link ${CELL_DIR}
        fi

    done < <(grep "${epitope},${cell2}" $LDIR/${epitope}.csv)

done < $DBDIR/epitopes # list of epitopes we are interested in

