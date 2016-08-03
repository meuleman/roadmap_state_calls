#!/bin/bash
# SGE Array for STEP 3 -- Call states on processed data:
cell=$( sed "${SGE_TASK_ID}q;d" $LDIR/available_marks.tsv )

echo "\n- STEP3 for ${cell}"
CELL_DIR=${TYPE_DIR}/${cell} 
CC_DIR=${TCALL_DIR}/${cell} # State calls 
mkdir -p ${CC_DIR}
cd ${CELL_DIR}

# TODO put logic gate on this:
source $BINDIR/code_ENCODE3_process_step3.sh $cell ${CELL_DIR} ${CC_DIR}




