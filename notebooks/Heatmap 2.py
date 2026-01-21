# %%

## HEATMAP 2 ## 

# Get the indices of the top 20 most significant genes
top_20_indices = res_clean.nsmallest(20, 'padj').index

# IMPORTANT: We now switch to log_counts (normalized) NOT raw counts
# We want Genes as ROWS and Samples as COLUMNS for z_score=0 to work on genes
heatmap_data = log_counts[top_20_indices].T 

# Map the symbols for the Y-axis labels
heatmap_data.index = heatmap_data.index.map(mapping_dict)

# Define Column Colors (MIS vs Control)
condition_colors = metaData['Condition'].map({'MIS': 'firebrick', 'C': 'royalblue'})

# CREATING HEATMAP 2
# z_score=0 now correctly standardizes each GENE (row) across all samples

g = sns.clustermap(heatmap_data, 
                   z_score=0, 
                   cmap='RdYlBu_r', 
                   col_colors=condition_colors, 
                   metric='correlation',
                   method='average',
                   annot=False, 
                   figsize=(10, 12))

# Aesthetic adjustments
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=9)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)

# Adding a Legend MIS vs Control
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='firebrick', label='MIS'),
                   Patch(facecolor='royalblue', label='Control')]

# We attach the legend to the color bar (ax_col_colors) or the figure
plt.legend(handles=legend_elements, bbox_to_anchor=(1.2, 1), 
           loc='upper left', title="Condition")

g.ax_heatmap.set_title('Top 20 DEGs: Z-score Normalized Expression', pad=20)
g.savefig("Heat_Map_MIS_vs_Control.png",dpi=300, bbox_inches = 'tight')
plt.show()
