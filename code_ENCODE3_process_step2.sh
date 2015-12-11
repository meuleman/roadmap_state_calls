#!/bin/bash -l

##########################################################################################
### HERE WE POOL ALL REPLICATES AND SUBSAMPLE TO 30M READS, ACCORDING TO ROADMAP
### WE ALSO RUN PHANTOMPEAKQUALTOOLS IN ORDER TO OBTAIN SOME QC STATS
##########################################################################################

cell_type=$1;
epitope=$2;

NREADS=30000000; # Final number of desired reads.

OUTPUTDIR="ENCODE3_data/${cell_type}";
mkdir -p ${OUTPUTDIR};

OFPREFIX="${OUTPUTDIR}/FINAL_${cell_type}_${epitope}";

FINAL_TA_FILE="${OFPREFIX}.tagAlign.gz"

# =================================
# Pool all tagAlign files for the same cell type and epitope, and subsample to 30M reads
# ================================
zcat ${OUTPUTDIR}/*_${cell_type}_${epitope}.filt.nodup.srt.SE.map.tagAlign.gz | \
  shuf -n ${NREADS} | sort -k1,1V -k2,2g | gzip -c > ${FINAL_TA_FILE}

# =================================
# Run 'run_spp_nodups.R' to obtain some QC stats on the final tagAlign file.
# ================================
Rscript ../../../software/phantompeakqualtools/run_spp_nodups.R \
  -rf -c=${FINAL_TA_FILE} -savp -out=${OFPREFIX}.qc

# Sanity check to easily see if we were successful:
cat ${OFPREFIX}.qc

### Now we're ready for calling chromatin states using ChromHMM on ${FINAL_TA_FILE}.



