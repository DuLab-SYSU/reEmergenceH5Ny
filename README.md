# reEmergenceH5Ny
Code and data needed to reproduce the results in article "Spatio-temporal Genotype Replacement of H5N8 AIVs Contributed to H5N1 Emergence in 2021/2022 Panzootic"

According to the permissions of GISAID (https://gisaid.org/), the actucal sequences used in this work were deleted from this repository. To reproduce the results, you can download related sequences from GISAID with the meta-information provided in data dir. 

./data_processing.ipynb notebook including:
- data pre-processing for epidemiology data (HPAI outbreak data downloaded from Epress-i (empress-i.apps.fao.org), number of cases and deaths downloaded from WOAH(woah.org) & human cases downloaded from WHO (who.int)).  
- visualization for Figure 1.  
- quality control & multiple sequence alignemnt & drop no HPAIV & clade 2.3.4.4b filter for sequences data. 
- construct ML tree for HA gene of all 2.3.4.4b HPAI H5 viruses (built using IQTree and re-rooted using TreeTime). 
- pairwise sequence identity calculation for each gene of all 2.3.4.4b HPAI H5. 
- visualization for Figure 2.

./blastn_h5n1.ipynb notebook including:
- visualization for blastn analysis for H5N1 viruses isolated between 2021-2022 (ED fig 1, fig S1-S8). 

./tree_prcessing.ipynb notebook including:
- construct ML tree for each segment of all 2.3.4.4b HPAI H5N8 viruses isolated between 2019-12-31 to 2022.
- calculate MPD value for each internal node.
- identify well supported monophyletic distinct group using a depth-first visit method.
- annoated genotype information for each strain based on monophyletic group of each gene.
- describe the time, location and host distribution for each genotype.
- visulization for Figure 3.

./blast.ipynb notebook including:
- blastn analysis using H5N8 sequences as query.
- ML trees were constructed for H5N8 sequences and retrived top 500 hits using Fasttree. The most recent common ancestor for each monophyletic group was identified. Sequences within the larger clade subtending the second ancestor node with bootstrap value no less then 70% were selected.

./tree_beauti.ipynb notebook including:
- visualization for the time, location and host estimation of MRCA and MRCGA of each monophyletic group of G0 and G1 H5N8 genotypes (ED fig 2, ED fig 3).
- visualization for the early spread of G0 and G1 H5N8 genotypes (Fig4 A, C)
