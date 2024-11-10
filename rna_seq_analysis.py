import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

# Reading gene count data from a CSV file
count_data = pd.read_csv('counts.csv', index_col=0)

# Descriptive data (e.g., condition: control and treated)
# Ensure that the order of the samples in the metadata matches the order of the columns in the gene count data
col_data = ['control', 'treated'] * (count_data.shape[1] // 2)

# Splitting the data into groups based on the condition
control_data = count_data.iloc[:, :len(col_data)//2]
treated_data = count_data.iloc[:, len(col_data)//2:]

# Performing a T-test to check for differences between the groups
p_values = []
log2_fold_changes = []

for gene in count_data.index:
    control_counts = control_data.loc[gene].values
    treated_counts = treated_data.loc[gene].values
    # Perform the T-test
    t_stat, p_val = ttest_ind(control_counts, treated_counts)
    p_values.append(p_val)
    
    # Calculate the relative change (log2 fold change)
    log2_fc = np.log2(np.mean(treated_counts) / np.mean(control_counts))
    log2_fold_changes.append(log2_fc)

# Create a DataFrame for the results
results = pd.DataFrame({
    'Gene': count_data.index,
    'log2FoldChange': log2_fold_changes,
    'pValue': p_values
})

# Adjust the results to obtain the adjusted p-values (FDR) using the Benjamini-Hochberg method
results['padj'] = results['pValue'] * len(results) / (np.argsort(results['pValue']) + 1)

# Extract significantly differentially expressed genes (padj < 0.05)
significant_genes = results[results['padj'] < 0.05]

# Print the significantly differentially expressed genes
print(significant_genes)

# Plot the Volcano plot
plt.figure(figsize=(10, 6))
sns.scatterplot(data=results, x='log2FoldChange', y=-np.log10(results['padj']), alpha=0.5)
plt.axhline(y=-np.log10(0.05), color='red', linestyle='--')
plt.axvline(x=0, color='gray', linestyle='--')
plt.title('Volcano Plot')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted P-value')
plt.show()
