#!/bin/bash -l

###################################################
### HERE WE POOL ALL REPLICATES AND SUBSAMPLE TO 30M READS, ACCORDING TO ROADMAP
### WE ALSO RUN PHANTOMPEAKQUALTOOLS IN ORDER TO OBTAIN SOME QC STATS
###################################################
SOFTWARE_DIR="/broad/compbio/meuleman/software"
SPP_EXISTS=0 # enable if you have it installed.

cell_type=$1;
epitope=$2;
CELL_DIR=$3

NREADS=30000000; # Final number of desired reads.

OFPREFIX="${CELL_DIR}/FINAL_${cell_type}_${epitope}";
FINAL_TA_FILE="${OFPREFIX}.tagAlign.gz"

# Check if gzipped file is empty:
if LC_ALL=C gzip -l ${FINAL_TA_FILE} | awk 'NR==2 {exit($2!=0)}'
then
    # =================================
    # Pool all tagAlign files for the same cell type and epitope, and subsample to 30M reads
    # ================================
    echo "Pool files and subsample to create ${FINAL_TA_FILE}"

    # Count number of reads left after filtering:
    REPL_READS=$( zcat ${CELL_DIR}/*_${cell_type}_${epitope}.filt.nodup.srt.SE.map.tagAlign.gz | wc -l )
    echo "There are ${REPL_READS} reads available from all replicates for ${cell_type} + ${epitope}" 
    if (( ${REPL_READS} >= $NREADS )) 
    then
        echo "Subsampling ${REPL_READS} to ${NREADS}"
        zcat ${CELL_DIR}/*_${cell_type}_${epitope}.filt.nodup.srt.SE.map.tagAlign.gz | shuf -n ${NREADS} | sort -k1,1V -k2,2g | gzip -c > ${FINAL_TA_FILE}
    else
        echo "Not enough reads! ${REPL_READS} > $NREADS"
    fi
fi 

# NOTE: need R library spp from http://compbio.med.harvard.edu/Supplements/ChIP-seq/ 
if [[ ! -s ${OFPREFIX}.qc && "${SPP_EXISTS}" == "1" ]] 
then
    # =================================
    # Run 'run_spp_nodups.R' to obtain some QC stats on the final tagAlign file.
    # ================================
    echo "QC stats on the final tagAlign file. Yields ${OFPREFIX}.qc"
    Rscript ${SOFTWARE_DIR}/phantompeakqualtools/run_spp_nodups.R -rf -c=${FINAL_TA_FILE} -savp -out=${OFPREFIX}.qc

    # Sanity check to easily see if we were successful:
    cat ${OFPREFIX}.qc
fi

