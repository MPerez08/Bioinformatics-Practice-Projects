# %%

# PLOTTING VOLCANO ---
# Setup thresholds
p_thresh = 0.05
lfc_thresh = 1.0 

# Filter for highlighting
res_clean = res.dropna(subset=['padj'])
up = res_clean[(res_clean['padj'] < p_thresh) & (res_clean['log2FoldChange'] > lfc_thresh)]
down = res_clean[(res_clean['padj'] < p_thresh) & (res_clean['log2FoldChange'] < -lfc_thresh)]

# Identify Top 10 genes to label (by significance)
top_genes = res_clean.sort_values('padj').head(10)

# Visualization
plt.figure(figsize=(10, 7))

# Plot all genes (Grey background)
plt.scatter(res_clean['log2FoldChange'], -np.log10(res_clean['padj']),
            c='grey', alpha=0.3, s=5, label='Not Significant')

# Plot Up-regulated (Firebrick)
plt.scatter(up['log2FoldChange'], -np.log10(up['padj']),
            c='firebrick', s=15, label=f'Up-regulated (n={len(up)})')

# Plot Down-regulated (Royalblue)
plt.scatter(down['log2FoldChange'], -np.log10(down['padj']),
            c='royalblue', s=15, label=f'Down-regulated (n={len(down)})')

# Add Labels for Top Genes
for i, row in top_genes.iterrows():
    # Use symbol if mapping worked, otherwise fallback to index (ENSG ID)
    label = row['symbol'] if ('symbol' in row and pd.notnull(row['symbol'])) else i 
    plt.annotate(label, (row['log2FoldChange'], -np.log10(row['padj'])),
                 xytext=(5, 5), textcoords='offset points', 
                 fontsize=8, fontweight='bold', alpha=0.8)

# Add Threshold Lines
plt.axhline(-np.log10(p_thresh), color='black', linestyle='--', lw=1)
plt.axvline(lfc_thresh, color='black', linestyle='--', lw=1)
plt.axvline(-lfc_thresh, color='black', linestyle='--', lw=1)

# Labels and Title
plt.xlabel("$\log_2$ Fold Change (MIS vs Control)")
plt.ylabel("$-\log_{10}$ Adjusted P-value")

plt.title("Differential Gene Expression: MIS vs Control Samples")
plt.legend(loc='upper left')

# Add a grid for better readability
plt.grid(alpha=0.2)

plt.savefig("Volcano_Plot_MIS_vs_Control.png",dpi=300, bbox_inches = 'tight')
print("Image Saved")

plt.show()
