# Analysis of MIS-C Gene Expression Profiles using PyDESeq2

## Project Overview

This project investigates the transcriptional differences between patients with Multisystem 
Inflammatory Syndrome (MIS) vs a Control group (normal healthy participants). The goal is to identify significantly
differentially expressed genes that may serve as biomarkers or provide insight into the inflammatory response. 

## Data Source
Dataset: E-CURD-149 (RNA-seq raw counts) 
Publication: Jackson HR, Miglietta L, Habgood-Coote D, D'Souza G, Shah P et al. (2023) 
" of Multisystem Inflammatory Syndrome in Children by a Whole-Blood Transcriptional Signature."

## Data Filtering Explained
My initial raw counts file contained much baground noise coming from samples with 0 reads. 
To reduce data being skewed, I filtered the file to only keep genes with at least 3 Sample IDs having a minimum of 10 raw counts or more 

```counts = counts.select_dtypes(include="number")
counts = counts[(counts >= 10).sum(axis=1) >= 3] ```


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

#Repositories Needed
## Web Application Used
Jupyer Notebook 

##Below are the libraries you'll need to download:

1) pandas: open sourced data analysis tool creating 2D DataFrame structures 
allowing for intelligent filtering and cleaning from seamless loading of a diverse set of 
sources inlcuding CSV, excel, SQL, JSON, and HDF5

2) numpy: fundamental package for computing multidimensional arrays packed with high level mathematical functions to analyze data

3) seaborn: Python data visualization providing high level interfaces and drawing informative statistical graphics with minimal code

4) sklearn(Scikit-Learn):  open sourced machine learning library that performs "non-deep learning" focusing on predictive data analysis
                           integrates well with Pandas Dataframes and NumPy arrays

5) pydeseq: Python version of the DESeq2 which is a specialized package designed to analyze high throughput RNA-Seq data 


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
