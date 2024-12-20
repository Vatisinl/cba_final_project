# cba_final_project
Here we present code for our final project on Computational Biology of Aging course done by Mazalov Aleksei (PhD 1), Malygina Alexandra (MSc 2), Shipulina Eva (MSc 2). 
# ECM transcriptome dynamics during aging
## Introduction




## Results

The next stage of the project was to conduct WGCNA to find gene modules that change during ageing process. Due to technical reasons we used only samples from heart, liver and lung. 

Data were preprocessed: metadata was created based on GTEx Portal files: GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt and GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt; after sample filtration gene filtration was performed: min 70 counts per a sample; DESeq normalization was performed as well as variance stabilizing transformation. Then power of correlation was chosen for every tissue such minimum power that scale free topology model fit showed in R^2 more than 0.8. 

Then we constructed gene modules using blockwiseModules() with maxBlockSize = 7000 and mergeCutHeight = 0.4. If correlation between at least one pair of modules was higher than 0.7 mergeCutHeight was increased on 0.2 until there was no such modules anymore. As a result the most significant modules contained several thousands of genes. 

Then we took the most significant modules from each tissue and worked with them.

### Heart
```{figure} figs/dendro_plot_Heart.png
Separation of modules for heart samples.
```
We see several separated by colours modules. The pattern is good so we proceeded with this. 

```{figure} figs/significance_plot_Heart.png
Significant difference between age groups.
```
We saw that all three age groups significantly differs one from each other. Therefore there will be findings in modules changes.  

```{figure} figs/top_module_heatmap_plot_Heart.png
Separation of samples for heart tissues.
```
Then we constructed a plot for the most significant module which is ME1 and saw that, possibly, samples are not very homogeneous.

To identify biological meaning of ME1 we performed GO- and KEGG-enrichment. It showed a focus on neurogenesis (which can be expected since heart is innervated) but also on T-helper cells differentiation and Malaria. We do not now causes of death of tissue donors from GTEx portal due to ethical reasons. So it causes questions about the reliability of the approach based on this database in the context of ageing studies: if the person died being young, maybe there are important system-wide changes that differ this patient from the healthy ageing person. 

Then we intersected genes in ME1 with matrisome genes and got 143 genes (DCN, VCAN, COL9A2, ELN, COL11A1, NTN1, FBLN1, COL4A4, EPYC, COL16A1, LTBP4, COL9A3, PAPLN, MXRA5, SRPX, COMP, PCOLCE, AEBP1, OGN, ASPN, ECM2, LGI1, COL1A1, TECTA, VWA5A, MGP, COL12A1, COL9A1, SMOC2, IMPG1, SPARC, THBS4, EFEMP1, FN1, PRG4, MFAP2, LTBP2, TNN, TGFBI, CRISPLD1, FMOD, SRGN, TNFAIP6, MATN4, COL21A1, FGL2, COL5A1, LAMA5, EMILIN2, MATN2, POSTN, LAMC1, CHAD, IGFBPL1, EMILIN1, CILP, MMRN1, FRAS1, LUM, KERA, FBLN5, MFAP1, MFGE8, IGFBP4, COL6A1, NTN5, TINAGL1, FNDC7, DPT, HMCN1, ECM1, FBLN7, COL8A1, SLIT2, VWA5B2, HAPLN1, RSPO3, VWDE, IGFBP3, IGFBP1, RSPO2, CRIM1, SPOCK1, IGSF10, LGI2, LGI4, ABI3BP, ACAN, SPON2, CILP2, NTN3, MATN1, NTNG1,SNED1,COL6A3 ,IGFBP7 ,FBLN2 ,PCOLCE2 , ESM1 , BMPER, COL1A2 ,FNDC1 ,SBSPON ,CTHRC1 ,SVEP1, FBN1 , MFAP4, IGFBP6 ,LTBP3 ,TNXB ,COL3A1 ,COL4A3 ,RSPO1, THBS3 ,COL22A1 ,COL24A1 ,COL8A2 ,MMRN2 ,PODN ,BGN , TSKU ,COL18A1 ,SLIT3, NELL2, THBS2, EMID1, HAPLN4, COL14A1,COL4A5, AGRN, PRELP , RELN, NTNG2, COL27A1 ,COL13A1, MFAP5, SMOC1,COL5A2, COL15A1, COL6A6, COL28A1, BGLAP, SPON1). This intersection is not at random according to Chi-square test (p < 0.05). It would be interesting to cross them with literature sources devoted to ageing of heart tissues. 

### Liver
The most significant gene module for the Liver samples gave no GO and KEGG biological processes groups and no intersection with matrisome genes.

### Lung
Here we have the same idea of homogeneity for samples. 
```{figure} figs/top_module_heatmap_plot_Lung.png
Separation of samples for heart tissues.
```
All age groups significantly differs from each other. 
```{figure} figs/significance_plot_Lung.png
Significant changes between groups.
```

Also GO- and KEGG-enrichment gave more interesting results: during ageing somehow respiration, oxidative phosphorilation and mitochondrial gene expression change. It would be interesting to trace particulas genes of interest associated with such changes and to show their dynamics. Additionally, intersection with matrisome genes list gave 40 genes (IBSP, LAMC3, LAMA3, LAMC2, COL11A1, COL17A1, COL4A4, IGFALS, CHADL, COL20A1, LAMA1, SPOCK2, TECTA, FN1, SPP1, FGL2, LAMA5, MATN3, EMILIN2, KCP, CHAD, MFAP1, FBN3, HAPLN1, VWDE, DMP1, VWA2, COL4A3, RSPO1, LAMB2, NDNF, VWA1, CRELD2, GLDN, DMBT1, EYS, COL4A5, AGRN, NTNG2, LAMB3). It would be also interesting to find out are there protein coding genes related to cell adjusting cell connections in lungs or maintaining structure of alveoli. 

```{figure} figs/top_module_heatmap_plot_Lung.png
Separation of samples for heart tissues.
```

## Discussion

WGCN analysis did not show interesting findings which is expected. To find gene modules which can be interesting in terms of ageging we should use smaller step for gene network construction, more homogeneous samples. It would be interesting to follow up particular genes expression through age groups to understand dynamics.

## Credits
This text prepared by [Team member 1](https://linktoyourprofile/scholar/or/linkedin.com) ...

## References

```{bibliography}
:style: plain
:filter: docname in docnames
```
