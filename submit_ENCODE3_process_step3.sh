#!/bin/bash

mkdir -p out;

cell_types=`ls ENCODE3_data`;

for cell_type in ${cell_types}; do
  echo "$cell_type";
  bsub -q compbio-week -P compbiofolk -J ENCODE3_process_step3_${cell_type} \
       -oo out/output_ENCODE3_process_step3_${cell_type}.out \
       -R "rusage[mem=9]" ./code_ENCODE3_process_step3.sh ${cell_type}
done;


