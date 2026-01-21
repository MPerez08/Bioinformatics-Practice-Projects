# %%

## HEATMAP 1 ## 

#DATA PREP#
import seaborn as sns
# 1. Get the indices (ENSG IDs) of the top 20 most significant genes (lowest padj)
#Filtering out NaNs first ensuring we valid p-values
top_20_indices = res.dropna(subset=['padj']).nsmallest(20, 'padj').index

#2. Extracting the genes from counts dataframe
heatmap_data = counts.loc[top_20_indices]

#3. Map the symbols for the Y-axis labels
heatmap_data.index = heatmap_data.index.map(mapping_dict)

#CREATING HEATMAP#
# Standardizing the data (Z-score) so we see relative changes across rows, 
# otherwise high-expression genes would wash out the whole plot.

g = sns.clustermap(heatmap_data,
                   z_score=0, # 0 indicates row-wise standardization
                   cmap='RdYlBu_r',
                   annot=False,
                   figsize=(10,12))


g.ax_heatmap.set_xticks(np.arange(len(heatmap_data.columns)) + 0.5) # Forces the x-axis to show a tick for every single column
g.ax_heatmap.set_xticklabels(heatmap_data.columns, rotation=90, fontsize=7)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)
g.ax_heatmap.set_title('Top 20 Most Significant Genes (Z-score Normalized)')  

g.savefig("Heat_Map_Top_20_Genes.png",dpi=300, bbox_inches = 'tight')
