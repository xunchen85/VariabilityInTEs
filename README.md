# Variability_In_TEs

## Introduction
Scripts included in this repository were used for the main analyses in the article of "Transposable elements are associated with the variable response to influenza infection" (In submission).

## inputs
Files under this "inputs" folder were used as the inputs for the R scripts. Zipped files should be unzipped before the running. Summary.5mC.enriched.table and Summary.expAndCentroid.enriched.table files can be shared upon request. They were not uploaded due to the limitation of file sizes but will be further deposited to zenodo.

## shell scripts
TEPeak_detection.sh shell script was used to detect TE instances that were overlapped with ATAC-seq or Chip-seq peaks. 
TEPeak_shuffle.sh was used to generate the expected distribution of peaks that are overlapped with each TE family.
These scripts are optimied based on the scripts wrote by Bordan et al. (https://github.com/lubogdan/ImmuneTE) Python scripts were written and used in these shell scripts. 

## R scripts
R scripts were used for most of the main analyses performed in the article. They can be ran one by one by following the orders. After the running, the scripts will also generate corresponding figures in pdf format.

