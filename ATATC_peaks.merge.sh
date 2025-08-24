#this script used for the merge the peaks from different ATAC-seq biological replicates
#this analysis is running for Constant and Natural light conditions
#this work is running on HPC
#this work is done by Huifang Yuan

#!/bin/bash

REPS=(L1_S9.last_peaks.sorted.narrowPeak L2_S10.last_peaks.sorted.narrowPeak L3_S11.last_peaks.sorted.narrowPeak  L4_S12.last_peaks.sorted.narrowPeak L6_S13.last_peaks.sorted.narrowPeak L7_S14.last_peaks.sorted.narrowPeak L8_S15.last_peaks.sorted.narrowPeak L10_S16.last_peaks.sorted.narrowPeak) 

# 1. candidate pooled peaks
cat "${REPS[@]}" | sort -k1,1 -k2,2n | bedtools merge > Light_pooled_peaks.bed

# 2. ( bedtools pooled peak replicate with overlap 25% )
> Light_overlap_counts.bed

for rep in "${REPS[@]}"; do
  bedtools intersect -wo -a Light_pooled_peaks.bed -b $rep | \
  awk 'BEGIN{OFS="\t"} {
    lenA = $3 - $2;
    lenB = $12 - $11;
    overlap = $NF;
    if ((overlap / lenA >= 0.25) || (overlap / lenB >= 0.25)) {
      print $1, $2, $3
    }
  }' >> Light_temp25_overlap.bed
done

# 3. filter the peaks with at least 5 of 8 replicates support
sort Light_temp25_overlap.bed | uniq -c | awk '$1 >= 5 {print $2"\t"$3"\t"$4}' > 001.Light_filtered_25peaks_5of8.bed

## result file: 001.Light_filtered_25peaks_5of8.bed is at ~/01.ATAC_analysis_2025/R_all/02.peaks/001.light.peaks.all
## do the same analysis for the Natural samples

