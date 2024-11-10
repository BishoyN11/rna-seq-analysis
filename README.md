# RNA-Seq Analysis

This project performs RNA-Seq data analysis using Python. The analysis includes reading RNA-Seq data files, processing the data, and visualizing gene expression patterns between different samples.

## How to use:

### Step 1: Prepare Your Data
The project requires BigWig files containing RNA-Seq data. You can find these files or use your own data.

### Step 2: Modify the Script to Specify Your Files
In the script `rna_seq_analysis.py`, you need to update the paths for the `control_file` and `treated_file`.

For example, change the following lines in the script to the correct paths to your files:

```python
# Define paths to BigWig files
control_file = 'path/to/control_file.bw'  # Replace with the path to your control file
treated_file = 'path/to/treated_file.bw'  # Replace with the path to your treated file

# Specify the genomic region (example: chromosome 1, position 0-10000)
chrom = 'chr1'
start = 0
end = 10000

control_file = 'GSM8608549_spleen_WT_LCMV_Qiazol_50069.bw'
treated_file = 'GSM8608550_spleen_WT_LCMV_Qiazol_50070.bw'

pip install -r requirements.txt

python rna_seq_analysis.py

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
control_file = 'GSM8608549_spleen_WT_LCMV_Qiazol_50069.bw'  # Replace with control file path
treated_file = 'GSM8608550_spleen_WT_LCMV_Qiazol_50070.bw'  # Replace with treated file path

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
