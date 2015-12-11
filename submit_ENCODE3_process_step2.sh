#!/bin/bash

mkdir -p out;

cell_types=`ls ENCODE3_data`;
epitopes="WCE H3K4me1 H3K4me3 H3K27me3 H3K36me3 H3K9me3 H3K27ac";

for cell_type in ${cell_types}; do
  for epitope in ${epitopes}; do
    echo "$cell_type / $epitope";
      bsub -q priority -P compbiofolk -J ENCODE3_process_step2_${cell_type}_${epitope} \
           -oo out/output_ENCODE3_process_step2_${cell_type}_${epitope}.out \
           -R "rusage[mem=5]" ./code_ENCODE3_process_step2.sh ${cell_type} ${epitope};
  done;
done;


