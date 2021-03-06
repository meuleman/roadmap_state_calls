#!/bin/bash -l

###################################################
### THIS CODE MAKES USE OF ANSHUL KUNDAJE'S ENCODE3 PROPOSAL FOR DATA PROCESSING.
### THE CURRENT VERSION OF THIS DOCUMENT CAN BE FOUND HERE:
### https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit
###################################################
# UMAPDIR="/broad/compbio/anshul/projects/umap"

id=$1;
cell=$2;
epitope=$3;
link=$4; # download link
CELL_DIR=$5; # cell type specific file directory

RLEN=36; # Desired read length, based on Roadmap data.
MAPQ_THRESH=30

# PREFIXES AND FILES: 
OFPREFIX=${CELL_DIR}/${id}_${cell}_${epitope}
RAW_BAM_FILE=${OFPREFIX}.bam

# Filtering files: 
FILT_BAM_PREFIX="${OFPREFIX}.filt.srt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc" # QC file for marking duplicates

# Final BAM files:
FINAL_BAM_PREFIX="${OFPREFIX}.filt.nodup.srt"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" # To be stored
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file
PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc" # quality stats

# Final TA files:
FINAL_TA_FILE="${FINAL_BAM_PREFIX}.SE.tagAlign.gz"
FINAL_TA_MAP_FILE="${FINAL_BAM_PREFIX}.SE.map.tagAlign.gz"

###################################################
### STEP 0 --- Download data:
###################################################
if LC_ALL=C gzip -l ${FINAL_TA_MAP_FILE} | awk 'NR==2 {exit($2!=0)}'  # Run if final TA file empty
then 
    echo "Empty/Missing final TagAlign file"

    if samtools view -F 0x904 -c ${RAW_BAM_FILE} | awk '{exit($1!=0)}'
    then 
        echo "Downloading ${RAW_BAM_FILE}"
        wget ${link} -o ${RAW_BAM_FILE}.log -O ${RAW_BAM_FILE} # BAM File
        echo "Indexing ${RAW_BAM_FILE}"
        samtools index ${RAW_BAM_FILE} # Index
        # TODO CHECK DOWNLOAD.
    fi

    ###################################################
    ### STEP 1B --- BASED ON ANSHUL'S ENCODE3 PROPOSAL
    ###################################################
    if samtools view -F 0x904 -c ${FILT_BAM_FILE} | awk '{exit($1!=0)}'
    then
        # =============================
        # Remove  unmapped, mate unmapped
        # not primary alignment, reads failing platform
        # Remove low MAPQ reads
        # ==================  
        echo "Filtering to get ${FILT_BAM_FILE}"
        samtools view -F 1805 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} | samtools sort - -T ${FILT_BAM_PREFIX}.tmp -o ${FILT_BAM_FILE}
        samtools view -H ${FILT_BAM_FILE} | grep SO

    # If still didnt work, likely Phred64
    # if samtools view -F 0x904 -c ${FILT_BAM_FILE} | awk '{exit($1!=0)}'
    #     echo "Didn't work -- converting to Phred64"
    #     # samtools view -h ${RAW_BAM_FILE} | | rescale_quals.py | samtools view -bS - > out.bam

    #     samtools view -F 1805 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} | samtools sort - -T ${FILT_BAM_PREFIX}.tmp -o ${FILT_BAM_FILE}
    #     samtools view -H ${FILT_BAM_FILE} | grep SO
    # fi

        # ========================
        echo "Marking duplicates for ${FILT_BAM_FILE}"
        # ======================
        TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
        MARKDUP="/seq/software/picard/1.802/bin/MarkDuplicates.jar";
        #MARKDUP="/seq/software/picard/current/bin/MarkDuplicates.jar"; # might break

        java -Xmx4G -jar ${MARKDUP} INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false # Crashes if low memory

        mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}
    fi

    # Make BAM file if it is empty
    if samtools view -F 0x904 -c ${FINAL_BAM_FILE} | awk '{exit($1!=0)}'
    then
        # ============================
        # Remove duplicates
        # Index final position sorted BAM
        # ============================
        echo "Removing duplicates to get ${FINAL_BAM_FILE}"
        samtools view -F 1804 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}

        echo "Indexing ${FINAL_BAM_FILE}"
        samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}

        samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}

        # =============================
        # Compute library complexity
        # =============================
        # Obtain unique count statistics

        # PBC File output:
        # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
        bedtools bamtobed -i ${FILT_BAM_FILE} | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | \
            awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\n",mt,m0,m1,m2}' > ${PBC_FILE_QC}
    fi

    # Remove temporary data if previous steps have non-zero output
    if samtools view -F 0x904 -c ${FINAL_BAM_FILE} | awk '{exit($1==0)}'
    then
        rm ${FILT_BAM_FILE}
        rm ${RAW_BAM_FILE} ${RAW_BAM_FILE}.bai
    fi

    ###################################################
    ### STEP 2A --- BASED ON ANSHUL'S ENCODE3 PROPOSAL
    ###################################################

    # ===================
    # Create tagAlign file
    # ===================
    if LC_ALL=C gzip -l ${FINAL_TA_FILE} | awk 'NR==2 {exit($2!=0)}'
    then
        bedtools bamtobed -i ${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | \
            awk 'BEGIN{FS="\t";OFS="\t"} $6=="+"{$3=$2+'"${RLEN}"';print $0} $6=="-"{$2=$3-'"${RLEN}"';print $0}' | \
            gzip -c > ${FINAL_TA_FILE}
    fi

    ###################################################
    ### mapFilterTagAlignFiles.sh
    ### THIS CODE IS BASED ON ANSHUL'S mapFilterTagAlignFiles.sh CODE.
    ### IT CALLS filterUniqueReads TO REMOVE ANY POSSIBLE READS IN 'UNMAPPABLE' LOCATIONS
    ###################################################

    # Set output and log file name (TODO not necessarily best place to define these)
    OFNAME=$(echo ${FINAL_TA_FILE} | sed -r -e 's/\.tagAlign\.gz$/.map.tagAlign/g')
    LOGFILE=$(echo ${FINAL_TA_FILE} | sed -r -e 's/\.tagAlign\.gz$/.map.tagAlign.logfile/g')

    if LC_ALL=C gzip -l ${OFNAME} | awk 'NR==2 {exit($2!=0)}'
    then
        export TMPLOC=${TMP}/tmp_${RANDOM}${RANDOM};
        mkdir -p ${TMPLOC};
        export MCR_CACHE_ROOT=${TMPLOC};
        $BINDIR/filterUniqueReads -s=${SDIR} -u=${UDIR} -v=${LOGFILE} ${FINAL_TA_FILE} | grep -v 'Warning' | gzip -c > ${OFNAME}.gz
        rm -rf ${TMPLOC}

        # NOTE THAT THIS SEEMS TO NEED AT LEAST _SOME_ KIND OF ACCESS TO A $DISPLAY.
        # SOMETIMES WHEN IT DOESN'T HAVE THIS, IT CRASHES WITHOUT WARNING.

        # Sanity check to see if we were successful:
        zcat ${OFNAME} | tail
    fi

    if LC_ALL=C gzip -l ${OFNAME} | awk 'NR==2 {exit($2==0)}'
    then
        rm ${FINAL_BAM_FILE}
        rm ${FINAL_TA_FILE}
    fi

fi

