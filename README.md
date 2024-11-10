# RNA-Seq Analysis

This project performs RNA-Seq data analysis using Python. The analysis includes reading RNA-Seq data files, processing the data, and visualizing gene expression patterns between different samples.

## How to use:

1. **Prepare Your Data**:
   - The project requires BigWig files containing RNA-Seq data (e.g., `GSM8608549_spleen_WT_LCMV_Qiazol_50069.bw`).
   - Place these files in the project directory or specify the correct path in the script.

2. **Modify the Script**:
   - In the script `rna_seq_analysis.py`, update the file paths:
     ```python
     control_file = 'path/to/control_file.bw'
     treated_file = 'path/to/treated_file.bw'
     chrom = 'chr1'  # Replace with the chromosome you want to analyze
     start = 0  # Start position of the region
     end = 10000  # End position of the region
     ```

3. **Run the Code**:
   - After installing dependencies, run the script:
     ```bash
     python rna_seq_analysis.py
     ```

4. **View Results**:
   - The results will show differentially expressed genes based on a T-test, and a Volcano plot will visualize the significance of gene expression changes.

## Prerequisites:

Before running the project, make sure you have Python installed along with the required libraries. You can install the dependencies using:

```bash
pip install -r requirements.txt
