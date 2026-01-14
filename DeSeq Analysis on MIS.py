import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# 1. SETTING WORKING DIRECTORY (Safety check)
os.chdir(r"C:\Users\miles\Documents\Bioinformatics\Project Practice")

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

# 8. PLOTTING PCA
plt.figure(figsize=(8, 6))
for cond in pca_df["Condition"].unique():
    subset = pca_df[pca_df["Condition"] == cond]
    plt.scatter(subset["PC1"], subset["PC2"], label=cond, s=100, edgecolors='white')

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.legend()
plt.title("PCA of RNA-seq Samples")
plt.show()

# 9a. IDENTIFYING DRIVER GENES (Loadings) ---
loadings = pd.DataFrame(
    pca.components_.T,
    columns=['PC1', 'PC2'],
    index=pca_input.columns # This must match the index of pca_input
)

loadings['symbol'] = loadings.index.map(mapping_dict)

# Getting top 10 genes for PC1

top_pc1_genes = loadings['PC1'].abs().sort_values(ascending=False).head(10).index
print("\nTop genes driving PC1 variance:")

print(loadings.loc[top_pc1_genes, ['symbol', 'PC1']])

# 9b PLOT LOADINGS ---
def plot_top_loadings(loadings, pc='PC1', n=10):
    top_genes = loadings[pc].abs().sort_values(ascending=False).head(n).index
    plt.figure(figsize=(8,6))
    plt.barh(loadings.loc[top_genes, 'symbol'], loadings.loc[top_genes, pc], color='skyblue')
    plt.xlabel(f"Weight in {pc}")
    plt.title(f"Top {n} Genes Driving {pc}")
    plt.gca().invert_yaxis()
    plt.show()

plot_top_loadings(loadings, pc="PC1")

# 10. PLOTTING VOLCANO
# Setup thresholds
p_thresh = 0.05
lfc_thresh = 1.0 

# Filter for highlighting
res_clean = res.dropna(subset=['padj'])
up = res_clean[(res_clean['padj'] < p_thresh) & (res_clean['log2FoldChange'] > lfc_thresh)]
down = res_clean[(res_clean['padj'] < p_thresh) & (res_clean['log2FoldChange'] < -lfc_thresh)]

# Identifying Top 10 genes to label (by significance)
top_genes = res_clean.sort_values('padj').head(10)

# Visualization
plt.figure(figsize=(10, 7))

# Plotting all genes (Grey background)
plt.scatter(res_clean['log2FoldChange'], -np.log10(res_clean['padj']),
            c='grey', alpha=0.3, s=5, label='Not Significant')

# Plotting Up-regulated (Firebrick)
plt.scatter(up['log2FoldChange'], -np.log10(up['padj']),
            c='firebrick', s=15, label=f'Up-regulated (n={len(up)})')

# Plotting Down-regulated (Royalblue)
plt.scatter(down['log2FoldChange'], -np.log10(down['padj']),
            c='royalblue', s=15, label=f'Down-regulated (n={len(down)})')

# Adding Labels for Top Genes
for i, row in top_genes.iterrows():
    # Use symbol if mapping worked, otherwise fallback to index (ENSG ID)
    label = row['symbol'] if ('symbol' in row and pd.notnull(row['symbol'])) else i 
    plt.annotate(label, (row['log2FoldChange'], -np.log10(row['padj'])),
                 xytext=(5, 5), textcoords='offset points', 
                 fontsize=8, fontweight='bold', alpha=0.8)

# Adding Threshold Lines
plt.axhline(-np.log10(p_thresh), color='black', linestyle='--', lw=1)
plt.axvline(lfc_thresh, color='black', linestyle='--', lw=1)
plt.axvline(-lfc_thresh, color='black', linestyle='--', lw=1)

# Adding Labels and Title
plt.xlabel("$\log_2$ Fold Change (MIS vs Control)")
plt.ylabel("$-\log_{10}$ Adjusted P-value")

plt.title("Differential Gene Expression: MIS vs Control Samples")
plt.legend(loc='upper right')

# Adding a grid for better readability on GitHub
plt.grid(alpha=0.2)

plt.show()

# 11. CREATING HEATMAPS 

#DATA PREP
# 1. Get the indices (ENSG IDs) of the top 20 most significant genes (lowest padj)
#Filtering out NaNs first ensuring we valid p-values
top_20_indices = res.dropna(subset=['padj']).nsmallest(20, 'padj').index

# Extracting the genes from counts dataframe
heatmap_data = counts.loc[top_20_indices]

# Map the symbols for the Y-axis labels
heatmap_data.index = heatmap_data.index.map(mapping_dict)

# Creating a Heatmap for ALL SAMPLE ID's in Dataset

# Standardizing the data (Z-score) so we see relative changes across rows, 
# otherwise high-expression genes would wash out the whole plot.

g = sns.clustermap(heatmap_data,
                   z_score=0, # 0 indicates row-wise standardization
                   cmap='RdYlBu_r',
                   annot=False,
                   figsize=(10,12))

g.ax_heatmap.set_xticks(np.arange(len(heatmap_data.columns)) + 0.5) # 1. Force the x-axis to show a tick for every single column
g.ax_heatmap.set_xticklabels(heatmap_data.columns, rotation=90, fontsize=7)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)
g.ax_heatmap.set_title('Top 20 Most Significant Genes (Z-score Normalized)')  

# Creating a Heatmap for the 10 Sample IDs (MIS vs C)

# Geting the indices of the top 20 most significant genes
top_20_indices = res_clean.nsmallest(20, 'padj').index

# Genes set up as ROWS and Samples as COLUMNS for z_score=0 to work on genes
heatmap_data = log_counts[top_20_indices].T 

# Mapping the symbols for the Y-axis labels
heatmap_data.index = heatmap_data.index.map(mapping_dict)

# Defining Column Colors (MIS vs Control)
condition_colors = metaData['Condition'].map({'MIS': 'firebrick', 'C': 'royalblue'})

# CREATING HEATMAP
# z_score=0 now correctly standardizes each GENE (row) across all samples
g = sns.clustermap(heatmap_data, 
                   z_score=0, 
                   cmap='RdYlBu_r', 
                   col_colors=condition_colors, 
                   metric='correlation',
                   method='average',
                   annot=False, 
                   figsize=(10, 12))

# 6. Aesthetic adjustments
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=9)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)


# 7. Add Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='firebrick', label='MIS'),
                   Patch(facecolor='royalblue', label='Control')]
# We attach the legend to the color bar (ax_col_colors) or the figure
plt.legend(handles=legend_elements, bbox_to_anchor=(1.2, 1), 
           loc='upper left', title="Condition")

g.ax_heatmap.set_title('Top 20 DEGs: Z-score Normalized Expression', pad=20)
plt.show()
