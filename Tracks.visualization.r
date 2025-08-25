## this script is used to visuzalize the results of the ATAC peaks and coverage signals


################ for all 16 datasets 



library(BSgenome.Amphimedon.Auq)
library(Gviz)
library(rtracklayer)


setwd("01.ATAC_analysis_2025/R_all/")
## input all peaks
peakFiles <- c(
  "02.peaks/001.light.peaks.all/L1_S9.last_peaks.sorted.narrowPeak",
  "02.peaks/001.light.peaks.all/L2_S10.last_peaks.sorted.narrowPeak",
  "02.peaks/001.light.peaks.all/L3_S11.last_peaks.sorted.narrowPeak",
  "02.peaks/001.light.peaks.all/L4_S12.last_peaks.sorted.narrowPeak",
  "02.peaks/001.light.peaks.all/L6_S13.last_peaks.sorted.narrowPeak",
  "02.peaks/001.light.peaks.all/L7_S14.last_peaks.sorted.narrowPeak",
  "02.peaks/001.light.peaks.all/L8_S15.last_peaks.sorted.narrowPeak",
  "02.peaks/001.light.peaks.all/L10_S16.last_peaks.sorted.narrowPeak",
  "02.peaks/002.natural.peaks.all/N1_S1.last_peaks.sorted.narrowPeak",
  "02.peaks/002.natural.peaks.all/N2_S2.last_peaks.sorted.narrowPeak",
  "02.peaks/002.natural.peaks.all/N3_S3.last_peaks.sorted.narrowPeak",
  "02.peaks/002.natural.peaks.all/N5_S4.last_peaks.sorted.narrowPeak",
  "02.peaks/002.natural.peaks.all/N6_S5.last_peaks.sorted.narrowPeak",
  "02.peaks/002.natural.peaks.all/N8_S6.last_peaks.sorted.narrowPeak",
  "02.peaks/002.natural.peaks.all/N9_S7.last_peaks.sorted.narrowPeak",
  "02.peaks/002.natural.peaks.all/N10_S8.last_peaks.sorted.narrowPeak"
  #"02.peaks/pre_comp.peaks/precompetent_1.last_peaks.narrowPeak",
  #"02.peaks/pre_comp.peaks/precompetent_2.last_peaks.narrowPeak",
  #"02.peaks/pre_comp.peaks/precompetent_3.last_peaks.narrowPeak",
 # "02.peaks/pre_comp.peaks/competent_1.last_peaks.narrowPeak",
 # "02.peaks/pre_comp.peaks/competent_2.last_peaks.narrowPeak",
# "02.peaks/pre_comp.peaks/competent_3.last_peaks.narrowPeak"
)
options(ucscChromosomeNames=FALSE)
peakDataList <- lapply(peakFiles, function(file) {
  import(file, format = "bed", extraCols = c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer"))
})

## generate tracknames
trackNames <- paste0("Sample_", 1:16)  
dataTracks <- lapply(seq_along(peakDataList), function(i) {
  DataTrack(
    range = peakDataList[[i]],
    type = "histogram",        
    name = trackNames[i],      
    genome = "Aq3.1",
    data = "signalValue",     
    col.histogram = "blue", 
    fill.histogram = "blue"
    #ylim = c(0.5, 0.5) 
  )
})

## 
axisTrack <- GenomeAxisTrack(name = "AquScaffold_2492")

##
library(ChIPseeker)
library(GenomicFeatures)
library(ggplot2)
library(clusterProfiler)
library(GenomicRanges)
library(dplyr)
library(DOSE)
#library(topGO)

txdb <- makeTxDbFromGFF("01.result_16samples/03.Aqu3.1_Gene.checked.addHead.3.gff3")
grtrack <- GeneRegionTrack(txdb, 
                           chromosome = "AquScaffold_14", 
                           start = 704000, 
                           end = 714000, 
                           name = "Genes in Region")


plotTracks(grtrack)

### 
seqTrack <- SequenceTrack(Aqu3)
## region
chrom <- "AquScaffold_14"
start <- 704000 
end <- 714000 

## 
trackList <- c(list(axisTrack), grtrack,  dataTracks, list(seqTrack))
total_tracks <- 1 + 1 + length(dataTracks) +1
plotTracks(
  trackList,
  chromosome = chrom,
  from = start,
  to = end,
  background.title = "grey",
  col.title = "white",
  collapse = FALSE,
  sizes = rep(1, total_tracks)
)


###################### for consensus datasets

library(BSgenome.Amphimedon.Auq)
library(Gviz)
library(rtracklayer)
options(ucscChromosomeNames=FALSE)
library(ChIPseeker)
library(GenomicFeatures)
library(ggplot2)
library(clusterProfiler)
library(GenomicRanges)


### bed files
bed_files <- c(
  "02.peaks/002.natural.peaks.all/002.Natural_filtered_25peaks_5of8.bed",
  "02.peaks/001.light.peaks.all/001.Light_filtered_25peaks_5of8.bed"
)

track_names <- c("natural","constant")



track_list <- lapply(seq_along(bed_files), function(i) {
  bed <- import(bed_files[i], format = "bed")
  #strand(bed) <- "*"
  AnnotationTrack(
    range = bed,
    genome = "Aqu3.1",
    name = track_names[i],
    chromosome = "AquScaffold_1235" 
  )
})



###############

bw_files <- c(
  "01.result_16samples/001.natural_all.filtered.bw",
  "01.result_16samples/002.constant_all.filtered.bw"
)

bw_names <- c("Natural", "Light")


bw_tracks <- lapply(seq_along(bw_files), function(i) {
  DataTrack(
    range = bw_files[i],
    genome = "Aqu3.1",
    name = bw_names[i],
    type = "histogram",    
    chromosome = "AquScaffold_1235",
    ylim = c(1,20)
  )
})



##
axisTrack <- GenomeAxisTrack(name = "AquScaffold_1235")

txdb <- makeTxDbFromGFF("01.result_16samples/03.Aqu3.1_Gene.checked.addHead.3.gff3")

grtrack <- GeneRegionTrack(txdb, 
                           chromosome = "AquScaffold_1235", 
                           start = 528000, 
                           end = 546000, 
                           name = "Genes in Region")




#strand(grtrack@range) <- "*"
plotTracks(grtrack)

### 
seqTrack <- SequenceTrack(Aqu3)
## region
chrom <- "AquScaffold_1235"
start <- 528000
end <- 546000 

### 

 
trackList <- c(list(axisTrack), grtrack, track_list, bw_tracks, list(seqTrack))
total_tracks <- length(trackList)

plotTracks(
  trackList,
  chromosome = chrom,
  from = start,
  to = end,
  background.title = "grey",
  col.title = "white",
  collapse = FALSE,
  sizes = rep(1, total_tracks)
)



