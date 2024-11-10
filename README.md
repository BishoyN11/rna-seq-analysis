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
If your files are in the same directory as the script, you can directly specify their names like this:

```python
control_file = 'GSM8608549_spleen_WT_LCMV_Qiazol_50069.bw'
treated_file = 'GSM8608550_spleen_WT_LCMV_Qiazol_50070.bw'
