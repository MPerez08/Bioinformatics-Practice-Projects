# %%
# 8. PLOTTING PCA
plt.figure(figsize=(8, 6))
for cond in pca_df["Condition"].unique():
    subset = pca_df[pca_df["Condition"] == cond]
    plt.scatter(subset["PC1"], subset["PC2"], label=cond, s=100, edgecolors='white')

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.legend()
plt.title("PCA of RNA-seq Samples")

# Saving Before Showing
plt.savefig("results/plots/pca_plot.png", dpi=300, bbox_inches='tight')
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
