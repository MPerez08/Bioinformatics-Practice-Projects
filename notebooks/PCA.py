# %%
#This notebook requires the output from DeSeq_Pipeline.py

# Cleaning up summary stats
res = stat_res.results_df
res_clean = res.dropna(subset=['padj'])

# GENE MAPPING
mg = mygene.MyGeneInfo() #Initializing mygene client
ensembl_ids = res_clean.index.tolist() #Pulling list of Ensembl ID's
gene_info = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human') #Querying the Database to get Gene ID Symbol (Human Specific)

# Creating a Mapping Dictionary {ENSG_ID: Symbol}
mapping_dict = {item['query']: item.get('symbol', item['query']) for item in gene_info}

# Defining clear labels for legend
label_map = {'C': 'Healthy Control', 'MIS': 'MIS Patient'}

# PLOTTING PCA
plt.figure(figsize=(8, 6))
for cond in pca_df["Condition"].unique():
    subset = pca_df[pca_df["Condition"] == cond]
    plt.scatter(subset["PC1"], subset["PC2"], 
                label=label_map.get(cond, cond), # Assigns the full name, 
                s=100, edgecolors='white')

# Labeling with variance explained
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")

# Title and Legend
plt.title("PCA of RNA-seq Samples: MIS vs Control", fontsize=14, pad=15)
plt.legend(title="Patient Group", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)
plt.grid(True, linestyle='--', alpha=0.3) # Grid to track coordinates

# Saving PCA Plot as a high-quality PNG
plt.savefig("PCA_MIS_vs_Control.png", dpi=300, bbox_inches='tight')


#Show Plot
plt.show()

# IDENTIFYING DRIVER GENES (Loadings)
loadings = pd.DataFrame(
    pca.components_.T,
    columns=['PC1', 'PC2'],
    index=log_counts.columns # This must match the index of pca_input
)

loadings['symbol'] = loadings.index.map(mapping_dict)

# Getting Top 10 genes for PC1
top_pc1_genes = loadings['PC1'].abs().sort_values(ascending=False).head(10).index
print("\nTop genes driving PC1 variance:")

print(loadings.loc[top_pc1_genes, ['symbol', 'PC1']])

# PLOT LOADINGS 
def plot_top_loadings(loadings, pc='PC1', n=10, save_name=None):
    top_genes = loadings[pc].abs().sort_values(ascending=False).head(n).index
    plt.figure(figsize=(8,6))
    
    # Plotting Code
    plt.barh(loadings.loc[top_genes, 'symbol'], loadings.loc[top_genes, pc], color='skyblue')
    plt.xlabel(f"Weight in {pc}")
    plt.title(f"Top {n} Genes Driving {pc}")
    plt.gca().invert_yaxis()
    
    #Saving Figure Prior To Showing
    if save_name is not None:
        plt.savefig(save_name, dpi=300, bbox_inches = 'tight')
        print(f"Plot saved as: {save_name}")
        
    plt.show()

plot_top_loadings(loadings, pc="PC1", save_name="PC1_Loadings.png")
