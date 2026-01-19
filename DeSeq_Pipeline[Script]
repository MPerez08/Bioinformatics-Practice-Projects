import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# 2. LOADING DATA
counts = pd.read_csv("E-CURD-149-raw-counts.tsv", sep="\t", index_col=0)
exp = pd.read_csv("E-CURD-149-experiment-design.tsv", sep="\t", index_col=0)

# 3. FILTERING
# Keep numeric columns and filter low-expression genes
counts = counts.select_dtypes(include="number")
counts = counts[(counts >= 10).sum(axis=1) >= 3] #will help reduce background noise later

# 4. CREATING A SUBSET OF SAMPLES ('Multisymptom Inflammatory Syndrome' vs 'Normal' Patients(Control)
normal = ('ERR11658908', 'ERR11658911', 'ERR11658920', 'ERR11658924', 'ERR11658937')
MIS = ('ERR11658903', 'ERR11658915', 'ERR11658916', 'ERR11658935', 'ERR11658918')
keep_ids = normal + MIS

fexp2 = exp.loc[exp.index.isin(keep_ids)] #filtering exp df for only the sample ID's in 'keep_ids'
fc = counts.loc[:, fexp2.index] #filtering counts df for only the sample ID's in 'fexp2'
fcx = fc.T # Transposing data to have sample ID's as indices 

# 5. METADATA SETUP
# Ensuring the condition list matches the 'keep_ids' list set
conditions = ['MIS','C','C','MIS','MIS','MIS','C','C','MIS','C']
metaData = pd.DataFrame({'Condition': conditions}, index=fcx.index) # Creates legend dataframe for future labeling

# 6. DESeq2 ANALYSIS
# Using design = "Condition" instead of design_factors to avoid DeprecationWarning
dds = DeseqDataSet(counts=fcx, metadata=metaData, design="~Condition") #Creating new object with filtered data for DESeq2 analysis
dds.deseq2()

# Initializing the Stats object
stat_res = DeseqStats(dds, contrast=("Condition", "MIS", "C"))

# Generating the results_df
stat_res.summary() 
res = stat_res.results_df 

# Saving stats report to a new file
res.to_csv("E-CURD-149_DESeq2_MIS_vs_Control.tsv", sep="\t")

# 7. NORMALIZATION & PCA PREP
# Incorporating normalized counts from dds
norm_counts = pd.DataFrame(dds.layers['normed_counts'], index=dds.obs_names, columns=dds.var_names)
log_counts = np.log2(norm_counts + 1) 

pca = PCA(n_components=2)
pca_result = pca.fit_transform(log_counts)

pca_df = pd.DataFrame(
    pca_result,
    index=log_counts.index,
    columns=["PC1", "PC2"]
)

#Joining Metadata
pca_df = pca_df.join(metaData)
