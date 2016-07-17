#!/bin/bash -l

###################################################
### HERE WE POOL ALL REPLICATES AND SUBSAMPLE TO 30M READS, ACCORDING TO ROADMAP
### WE ALSO RUN PHANTOMPEAKQUALTOOLS IN ORDER TO OBTAIN SOME QC STATS
###################################################

cell_type=$1;
epitope=$2;
CELL_DIR=$3

NREADS=30000000; # Final number of desired reads.

OFPREFIX="${CELL_DIR}/FINAL_${cell_type}_${epitope}";
FINAL_TA_FILE="${OFPREFIX}.tagAlign.gz"

if [[ ! -s ${FINAL_TA_FILE} ]]
then
    # =================================
    # Pool all tagAlign files for the same cell type and epitope, and subsample to 30M reads
    # ================================
    zcat ${CELL_DIR}/*_${cell_type}_${epitope}.filt.nodup.srt.SE.map.tagAlign.gz | \
        shuf -n ${NREADS} | sort -k1,1V -k2,2g | gzip -c > ${FINAL_TA_FILE}
fi 

if [[ ! -s ${OFPREFIX}.qc ]] 
then
    # =================================
    # Run 'run_spp_nodups.R' to obtain some QC stats on the final tagAlign file.
    # ================================
    # TODO FIND SOFTWARE!
    Rscript ../../../software/phantompeakqualtools/run_spp_nodups.R -rf -c=${FINAL_TA_FILE} -savp -out=${OFPREFIX}.qc

    # Sanity check to easily see if we were successful:
    cat ${OFPREFIX}.qc
fi

