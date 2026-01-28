# Analysis of MIS-C Gene Expression Profiles using PyDESeq2

## Project Overview
This repository contains a pipeline to identify differentially expressed genes (DEGs) in patients with Multisystem Inflammatory Syndrome in Children (MIS-C) compared to a healthy control group.

## Goal
Identify significantly differentially expressed genes that may serve as biomarkers or provide insight into the inflammatory response for clinical diagnostic purposes

## Data Source
Dataset: E-CURD-149 (RNA-seq raw counts) 
Publication: Jackson HR, Miglietta L, Habgood-Coote D, D'Souza G, Shah P et al. (2023) 
" of Multisystem Inflammatory Syndrome in Children by a Whole-Blood Transcriptional Signature."

# Data Pre-processing & Quality Control
## Expression Threshold
My initial raw counts file contained too much background noise coming from samples with 0 reads. 
To reduce my data skewing statistical models, I filtered the file to only keep genes with 1) numerical data and 2) At least 3 Sample IDs having a minimum of 10 raw counts or more for that gene.

```counts = counts.select_dtypes(include="number")``` & ```counts = counts[(counts >= 10).sum(axis=1) >= 3] ```

## Normalization
Data was log2-tranformed and normalized to account for differences in sequencing depth across samples. 

# Technical Summary: The DESeq2 Model
In this project, we utilized DESeq2, a gold-standard Bioconductor-based tool (ported to Python via pydeseq2), to perform differential expression analysis.

**Why DESeq2?**
Standard t-tests often fail in RNA-seq because gene expression data is heteroscedastic (the variance changes with the mean). DESeq2 addresses this through a negative binomial distribution ideal for modeling gene expression counts with high variance. With internal normalization, it also accounts for differences in sequencing depth and RNA composition differences. This allows for a more accurate representation of all genes participating in the sequencing analysis.

# Data Visualization Methods and Reasoning

## Principal Component Analysis (PCA)
Using PCA as visual method will allow me to see if there is an apparent batching effect or noticeable grouping between samples. My assumption is that the MIS-C samples will cluster together,
this would indicate their raw counts are correlated to the disease state and not simply biological noise

## Volcano PLot
Using a volcano plot visualization of all the genes will allow me to gauge the scale of 
change in the raw counts as well how statistically significant those changes were. It's expected the genes most upregulated and down regulated will be directly a part of biological pathways or signaling groups involved in inflammation responses.

## HeatMap
Using a heatmap visualization of all genes will allow me to gauge whether the up and down regulation of the top 20 most significant genes is consistent across all sample ID's or if it's simply a scattered phenotype.
Secondly narrowing down the heatmap to only compare the 5 control vs 5 MIS-C samples sample will help clarify if the regulation of those top genes is consistent among the groupings of the disease state. 


# Visual Analysis
**1.Principal Component Analysis (PCA)**
<img width="2711" height="1683" alt="PCA_MIS_vs_Control" src="https://github.com/user-attachments/assets/c5341f65-8ee8-4718-b66c-375cca1d5e83" />

> **Insight**: The PCA plot reveals that **48.3% of the total variance (PC1)** is attributed to the difference between MIS-C and Control conditions. This high percentage indicates a robust experimental signal.

> **Top Driver Genes:** Through loadings analysis, we identified that genes such as **OLAH** and **CD177** among others are the primary contributors to PC1. These 10 genes represent the core "variable" separating the disease state from the healthy baseline.

> **Secondary Variation:** PC2 (11.1%) captures inter-individual variation. While the groups separate cleanly on PC1, the vertical spread in PC2 suggests minor heterogeneity within the patient cohorts, which is expected in clinical human samples.

**2. Differential Expression (Volcano Plot)** 
<img width="2552" height="1873" alt="Volcano_Plot_MIS_vs_Control" src="https://github.com/user-attachments/assets/1fd30760-786c-40b1-bedf-ac25989a4375" />

> **Significance:** We defined significant genes as those with a $|Log2FC| > 1.0$ and an adjusted $p$-value $< 0.05$.

> **Log2 Fold Change (X-axis):** Represents the magnitude of change. A value of 1.0 indicates a doubling of gene expression in MIS-C patients compared to controls.

> **Statistical Significance (Y-axis):** We used the $-\log_{10}$ of the adjusted p-value. The higher the point, the more confident we are that the difference isn't due to random chance.

> **Result:** We identified **4522** up-regulated and **4208** down-regulated genes when comparing how much the MIS-C patients differed relative to the baseline of the control group.

> **Interpretation:** Genes plotted in the top-right quadrant are significantly over-expressed in MIS-C patients. These represent the molecular "on-switches" of the syndrome. Conversely, genes in the top-left quadrant are suppressed in the disease state. A considerable amount of upregulation was observed relative to significantly downregulated genes. 

**3a. Expression Consistency (Heatmap)**
<img width="2924" height="3480" alt="Heat_Map_Top_20_Genes" src="https://github.com/user-attachments/assets/c2adb019-704c-429b-b238-e001b51c6745" />

> **Interpretation:** The heatmap utilizes Z-score normalization to highlight relative expression. We observe a highly consistent expression pattern across the top 20 DEGs, confirming that the identified biomarkers are representative of the MIS-C group. Interestingly

**3b. Top 20 significant genes compared across 5 normal patients and 5 MIS-C diagnosed**
<img width="2971" height="3539" alt="Heat_Map_MIS_vs_Control" src="https://github.com/user-attachments/assets/e94819ba-b1ee-4c6d-9aaf-b396cf020be9" />
> **Interpretation:** This heatmap showed a clear separation across all 20 genes in regards to how much upregulation was observed across all 5 patients sampled with diagnosed MIS-C vs the control normal group. This further emphasizes the correlation of these Top 20 genes as potential biomarkers for diagnostic test that looks for MIS specifically.

#Genes of Interest
**Key Biological Findings**
The top genes driving the variance in PC1 include:
** Innate Immune Drivers**
1. **OLAH:**: involved in fatty acid metabolism seen to be significantly elevated in life threatning viral infections with diseases such as H7N9, Covid-19, RSV, and seasonal influenza.[1]
   
2. **CD177**: surface protein involved in recruiting neutrophils thorough blood vessels to reach inflammed tissues. Golden marker for neutrophil activation. 

3. **ADAMTS2**: involved in inflammation mediation, fibrosis, and immune cell regulation. Plays a role in remodeling of of extracelluar matrix during injury. Has shown elevated level patterns in other inflammation related diseases like sytemic lupus erythematosus. [2] 

5 & 6. **DAAM2 / DAAM2-AS1** : These are involved in the Wnt signaling pathway, which regulates cell polarity and cytoskeleton organization. In this context, they likely represent the structural changes cells undergo to move or survive within inflamed tissue. Known to modulate oligodendrocyte differentiation, a pathway involved in producing immune sensors and antigen presenting cells. [3]

7. **CCNA1**: Often associated with extracellular signal-regulated kinase (ERK) signaling or abnormal cell proliferation/repair attempts in response to tissue injury. ERK1/2 signaling is rapidly activated in response to tissue injury helping stem cells transition from quiescence to proliferation, allowing for tissue self-repair.[4]
   
8. **LINC02943**:  In MIS-C patients, LINC02943 acts as an enhancer for the ETS2 transcription factor. Elevated expression of ETS2 in monocytes drives the "cytokine storm" characteristic of MIS-C, where the immune system overreacts following a SARS-CoV-2 infection. [5]

9. **NUPR1**: Often called "P8," this is a stress-induced protein. It acts as a "hub" in the cell nucleus to prevent cell death during extreme stress (like a cytokine storm). If this is high in patients, it indicates their cells are in a high-intensity "survival mode" due to metabolic or inflammatory stress.[8]
    
10. **KCNMA1-AS3**: the parent gene "KCNMA1" acts as a central molecule in the NF-kappa B-dependent inflammatory response of macrophages. [7] 

#Acknowledged Limitations
The raw data in counts obtained from all +30,000 genes had a limited sample group of normal patient samples to obtain expression counts from. This was attributed to the limited release of data from the original, comprehensive study found in the Journal of Pediatric Infectious Diseases Society from Oxford Journals that was then partially shared over to Expression Atlas for public downloading. 

#Future Direction
#**1. Functional Enrichment Analysis (GO & KEGG)**

Individual genes are a great start to understand what specific protein pathways are being targeted, but biological processes happen in groups and being able to connect which genes are involved in what pathways created a new depth of understanding in the disease pathology and immune responses of those diagnosed with the MIS-C.
   
   **How to do it:** Perform Gene Ontology (GO) enrichment and KEGG pathway analysis. 
> Gene Ontology is a structured and standardized representation of the biological knowledge that’s accrued from scientific research in genes specifically. GO helps define terms and connections of the genes to each other by their relational function. 
> **KEGG (Kyoto Encyclopedia of Genes and Genomes):** : another bioinformatic approach that identifies significantly overrepresented biological pathways by comparing them against curated database pathways. 

**GOAL:** Instead of listing off genes with up or down-regulation phenotypes, using both tools mentioned above would define what activation pathways are enriched or modulated by the significant genes identified from the disease group. 

#2. Machine Learning: Biomarker Discovery
Since the identified genes could be potentially used as biomarkers for diagnostic purposes, using algorithms to test model the predictive power of each marker could help narrow down and determine the best candidates for future testing developments. 
 **How to do it:** Set up a simple classifier to be used on the top 20 DEGs as features using Random Forest or SVM
> Random Forest: a learning algorithm that builds out multiple decision trees based on nodes (classifying questions/variables) about the dataset in order to create predictive numerical or classification results on future data entries
> SVM: another machine learning algorithm used for classification and regression tasks. IT tries to find the best boundary (known as the hyperplane) that separates different classes in the data. Very useful in binary classifications. 



**Sources**
1. High expression of oleoyl-ACP hydrolase underpins life-threatening respiratory viral diseases (https://pubmed.ncbi.nlm.nih.gov/39137778/)
2. Analysis of expression and its clinical significance of the ADAMTS-2 in systemic lupus erythematosus
https://pubmed.ncbi.nlm.nih.gov/39806074)
3. New insights into the immunologic role of oligodendrocyte lineage cells in demyelination diseases
(https://pmc.ncbi.nlm.nih.gov/articles/PMC9548433)
4. Role of the Extracellular Signal-Regulated Kinase 1/2 Signaling Pathway in Ischemia-Reperfusion Injury
(https://www.frontiersin.org/journals/physiology/articles/10.3389/fphys.2019.01038/full)
5. Rare genetic variants involved in multisystem inflammatory syndrome in children: a multicenter Brazilian cohort study
(https://pmc.ncbi.nlm.nih.gov/articles/PMC10426286)
6. The ferroptosis inhibitor NUPR1 coordinates the mitochondrial response to oxidative stress and cell metabolism during COPD pathogenesis in the lung
doi: ( https://doi.org/10.1101/2025.04.27.650854)
7.  KCNMA1 potassium calcium-activated channel subfamily M alpha 1
(https://www.ncbi.nlm.nih.gov/gene/3778)
8. Lactylation-Driven NUPR1 Promotes Immunosuppression of Tumor-Infiltrating Macrophages in Hepatocellular Carcinoma
    DOI: (10.1002/advs.202413095)
