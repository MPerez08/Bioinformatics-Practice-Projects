# %%
#This notebook requires the output from DeSeq_Pipeline.py

# 1. Cleaning up summary stats
res = stat_res.results_df
res_clean = res.dropna(subset=['padj'])

# 2. GENE MAPPING
mg = mygene.MyGeneInfo() #Initializing mygene client
ensembl_ids = res_clean.index.tolist() #Pulling list of Ensembl ID's
gene_info = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human') #Querying the Database to get Gene ID Symbol (Human Specific)

# Creating a Mapping Dictionary {ENSG_ID: Symbol}
mapping_dict = {item['query']: item.get('symbol', item['query']) for item in gene_info}

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
    index=log_counts.columns # This must match the index of pca_input
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
