import pyBigWig
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

# Define function to read bigWig data
def read_bigwig(file_path, chrom, start, end):
    bw = pyBigWig.open(file_path)
    data = bw.values(chrom, start, end)
    bw.close()
    return np.nan_to_num(data)  # Convert NaNs to zeros or other default value

# Example: Reading data from bigWig files (modify these paths and region accordingly)
control_file = 'GSM8608549_spleen_WT_LCMV_Qiazol_50069.bw'
treated_file = 'GSM8608549_spleen_WT_LCMV_Qiazol_50069.bw'  # Replace with treated file

# Specify the genomic region (example: chromosome 1, position 0-10000)
chrom = 'chr1'
start = 0
end = 10000

# Read data from bigWig files for control and treated samples
control_data = read_bigwig(control_file, chrom, start, end)
treated_data = read_bigwig(treated_file, chrom, start, end)

# Perform T-test between control and treated data
t_stat, p_value = ttest_ind(control_data, treated_data)

# Create a DataFrame for the results (for demonstration)
results = pd.DataFrame({
    'Gene': [f'{chrom}:{start}-{end}'],
    'log2FoldChange': [np.log2(np.mean(treated_data) / np.mean(control_data))],
    'pValue': [p_value]
})

# Adjust p-values (FDR)
results['padj'] = results['pValue'] * len(results) / (np.argsort(results['pValue']) + 1)

# Print the significantly differentially expressed genes
significant_genes = results[results['padj'] < 0.05]
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
