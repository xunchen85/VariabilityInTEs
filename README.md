# Variability_In_TEs

## 1) Introduction
Scripts included in this repository were used for the main analyses in the article of "Transposable elements are associated with the variable response to influenza infection" (bioRxiv, doi: https://doi.org/10.1101/2022.05.10.491101).

## 2) inputs
Files under this *"inputs"* folder were used as the inputs for the R scripts. Zipped files should be unzipped before the running. *Summary.5mC.enriched.table* and *Summary.expAndCentroid.enriched.table* files can be shared upon request. They were not uploaded due to the limitation of file sizes but will be further deposited to Zenodo.

## 3) shell and python scripts
&emsp;•	*1.TEPeak_detection.sh*: detect TE instances that were overlapped with ATAC-seq or Chip-seq peaks.

&emsp;•	*2.TEPeak_shuffle.sh*: generate the expected distribution of peaks that are overlapped with each TE family.<br />

The two shell scripts were optimized based on the scripts written by Bordan et al. (https://github.com/lubogdan/ImmuneTE).

&emsp;Python scripts were written and used in these shell scripts. *Convert_to_PeakSummit.py*: convert the peak regions to the peak summits in BED format. *Combined_TEenrichmentByTEfamily.py*: organize the number of actual/shuffled peaks-associated instances per TE family. 

## 4) R scripts
R scripts were used for most of the main analyses performed in the article. They can be ran one by one by following the orders. After the running, the scripts will also generate corresponding figures in pdf format.

&emsp;•	*3.EMC-PCA_analysis.R*: PCA analysis of TE and gene expression among individuals (Figure 1A).

&emsp;•	*4.EMC-diff_TE_exp_analysis.R*: differential gene/TE expression analysis between flu infected and non-infected samples. DESeq2 was used for the analysis.

&emsp;•	*5.EMC-diff_TE_exp_analysis_2.R*: prepare the vacano plot of TE differential expression results. The script will also plot the proportion of TE families per subclass that are up-/down-regulated upon influenza infection (Figure 1C).

&emsp;•	*6.EMC-TE_gene_correlation_analysis.R*: correlation analysis between each TE/gene and computed viral load among 38 samples. We computed the correlation at the basal expression levels, expression levels post infection, and expression fold change of each TE/gene versus viral load (Figure S1A, Figure S1B).

&emsp;•	*7.EMC-TE_gene_correlation_analysis_2.R*: prepare the scatter plots of log2FC and R-squared among all TE families. It was also used to prepare the distribution of R-squareds of each TE subclass at different conditions. Examples of specific TE families were also visualized. (Figure 1D, Figure 1E, Figure 1F, Figure 6A, and Figure S7A)

&emsp;•	*8.EMC-correlated_TE_permutationTest.R*: permutation analysis of the proportion of positively or negatively correlated TE families (with viral loads) per superfamily/subclass (Figure S1C).

&emsp;•	*9.EMC-TE_enrichment_analysis.R*: identify TE families with enhanced/reduced activity (ATAC-seq, various Chip-seq) (Figure 2D, Figure 2E, Figure 2F, Figure2G).

&emsp;•	*10.EMC-TEfamilies_clustering.R*: clustering analysis of TE enrichment levels among 35 flu-infected and non-infected samples separately (Figure 3C, Figure S3B).

&emsp;•	*11.EMC-TE_regulation_analysis.R*: study the TE regulation roles of nearby genes (Figure 4A, Figure 4B). 

&emsp;•	*12.EMC-TE_epigenetic_status.R*: analyze the epigenetic status of all accessible TEs from 36 families with enhanced accessibility (Figure 4C).

&emsp;•	*13.Determine_TE_peak_Regions.R*: detect "TE peak regions" along each consensus sequence.

&emsp;•	*14.EMC-TE_motif_analysis.R*: identification of candidate binding motifs per TE family (Figure 5D, Figure S5C, Figure S6B).

&emsp;•	*15.EMC-predictive_models.R*: develop predictive models for viral load post-infection (Figure 6F, Figure 6H, Figure 6I, Figure 6J, Figure S7B, Figure S7C, Figure S7D, Figure S7E, Figure S7F, Figure S7G, Figure S7I). 

## 5) Citation
Xun Chen, Alain Sarabia Pacis, Katherine A Aracena, Saideep Gona, Tony Kwan, Cristian Groza, Yen Lung Lin, Renata Helena Monteiro Sindeaux, Vania Yotova, Albena Pramatarova, Marie-Michelle Simon, Tomi M. Pastinen, Luis Barreiro, Guillaume Bourque: **Transposable elements are associated with the variable response to influenza infection**. bioRxiv 2022.05.10.491101; doi: https://doi.org/10.1101/2022.05.10.491101

