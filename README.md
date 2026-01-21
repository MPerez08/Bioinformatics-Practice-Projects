# Analysis of MIS-C Gene Expression Profiles using PyDESeq2

## Project Overview
This repository contains a pipeline to identify differentially expressed genes (DEGs) in patients with Multisystem Inflammatory Syndrome (MIS) compared to a healthy control group.

## Goal
Identify significantly differentially expressed genes that may serve as biomarkers or provide insight into the inflammatory response for clincal diagnostic purposes

## Data Source
Dataset: E-CURD-149 (RNA-seq raw counts) 
Publication: Jackson HR, Miglietta L, Habgood-Coote D, D'Souza G, Shah P et al. (2023) 
" of Multisystem Inflammatory Syndrome in Children by a Whole-Blood Transcriptional Signature."

# Data Pre-processing & Quality Control
## Expression Threshold
My initial raw counts file contained too much baground noise coming from samples with 0 reads. 
To reduce my data skewing statistical models, I filtered the file to only keep genes with 1) numerical data and 2) At least 3 Sample IDs having a minimum of 10 raw counts or more for that gene.

```counts = counts.select_dtypes(include="number")``` & ```counts = counts[(counts >= 10).sum(axis=1) >= 3] ```

## Normalization
Data was log2-tranformed and normalized to account for differences in sequencing depth across samples. 

# Data Visualization Methods and Reasoning

## Principal Component Analysis (PCA)
Using PCA as visual method will allow me to see if there is an apparent batching effect or noticeable grouping between samples. My assumption is that the MIS samples will cluster together,
this would indicate their raw counts are correlated to the disease state and not simply biological noise

## Volcano PLot
Using a volcano plot visualization of all the genes will allow me to gauge the scale of 
change in the raw counts as well how statistically significant those changes were. My guess is that the genes that are most upregulated and down regulated will be directly a part of biological pathways or signaling groups involved in inflammation responses. 

## HeatMap
Using a heatmap visualization of all genes will allow me to gauge whether the up and down regulation of the top 20 most significant genes is consistent accross all sample ID's or if it's simply a scattered phenotype.
Secondly narrowing down the heatmap to only compare the 5 control vs 5 MIS samples sample will help clarify if the regulation of those top genes is consistent among the groupings of the disease state. 


# Visual Analysis
1.Principal Component Analysis (PCA)
![image alt](C:\Users\miles\Documents\Bioinformatics\Project Practice)
<img width="2711" height="1683" alt="PCA_MIS_vs_Control" src="https://github.com/user-attachments/assets/d99d2cf2-28ec-41f2-b8f6-3162b4159a8c" />

* Insight: We used PCA to visualize the high-dimensional data. PC1 accounts for [48.3%] of the variance. The clear separation     between MIS and Control groups suggests a distinct transcriptional signature for the syndrome.

## Exploratory Data Analysis: PCA
Before running differential expression, I performed Principal Component Analysis (PCA) to ensure the sample cluster by biological condition rather than technical noise.

> Interpretation: > PC1 captures 52% of the variance. We observe a clear separation between MIS and Control samples along the
> X-axis, suggesting.......

## Differential Expression Analysis 
I utilized PyDESeq2 to model gene expression. This approach accounts for the overdispersion 
typical of RNA-seq count data.

**Statistical Criteria:**
* Fold Change Threshold: |$|\log_2FC| > 1.0$ (2-fold change)
* Signicance Threshold Padj < 0.05 (FDR corrected)


```coding block ... ``` 















# Results Summary 
The volcano plot visalizes the differential expression results for all filtered genes across 
the 10 prioritized samples.

* Up-regulated (Red): Genes significantly higher in MIS patients.
* Down-regulated (Blue): Genes significantly lower in MIS patients.
* Key Findings: THe top 5 genes (annotated) show the highest statistical significance. Specifically, [GENE HERE] is a known inflammatory marker, validating the experimental results.


Down-Regulated 
[SYCP2L]: involved in negative regulation of programmed cell death (apoptosis). Crucial for female fertility, promotes                   survival of primordial oocytes

Up-Regulated
[SLX1A-SULT1A3]:  fusion gene (joined together and expressed as 1 unit)
           SLX1A:   normally involved in DNA repair, helps maintain genome stability, important when cells are under stress
           SULT1A3: involved in sulfanating neurotransmitters, drugs, & dietary compounds. Highly expressed in upper GI tract,                     brain and blood platlets
[CYP2A6]:   critical liver enzyme, known for metabolyzing nicotine, activated certain carcinogens (aflatoxin B1)
            lower expressed levels lead to slower nicotine metabolizing and lower chance of addiction.
            Also known for metabolizing other pharmaceuticals 
            
[ENSG00000266302]
[SLC12A1]: gene encoding instructions for making the NKCC2 protein, a crucial kidney transporter that reabsorbs sodium,                   potassium, and chloride back into the blood
