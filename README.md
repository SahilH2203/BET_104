
# The Neighborhood Effect On The Angle Between Residues

This pipeline extracts the secondary structure context and calculates the signed 3D angle between the $C_\alpha$ -> Centroid vectors of adjacent residues surrounding an Arginine (X-R-X) within an Alpha Helix (HHH).

# Plot

<img width="3000" height="1800" alt="ss_profile_HHH_for_arg" src="https://github.com/user-attachments/assets/cfb6b646-1ad8-4646-904e-e9585d4d7dd7" />

## Repository Structure

* **`Snakefile`**: The core Snakemake workflow script that automatically orchestrates the entire pipeline, from unzipping data to generating the final plot.
* **`config.yaml`**: Contains customizable parameters (like the target amino acid and output file names) allowing you to change pipeline behaviors without editing the code.
* **`scripts/context_stride.py`**: Parses the STRIDE secondary structure output to find valid X-R-X motifs specifically inside alpha helices and extracts their sequence context.
* **`scripts/calc_angles.py`**: Reads the 3D coordinate data from the raw PDB files, calculates the center of mass for the relevant sidechains, and computes the spatial angles.
* **`scripts/plot_angles.py`**: Takes the computed angles spreadsheet and generates a categorized Kernel Density Estimate (KDE) plot based on the steric bulk of neighboring residues.
* **`result/`**: The output directory where the final compiled dataset (`angles.tsv`), the list of successfully processed PDBs (`valid_pdbs.txt`), and the visualization are saved.

## Instructions to Reproduce

1. Clone this repository.
2. Create a folder named `pdbs/` in the root directory and place all raw `.pdb` (or `.pdb.gz`) data files inside it.
3. Ensure you have the required dependencies installed (`snakemake`, `stride`, `biopython`, `numpy`, `pandas`, `seaborn`, `matplotlib`).
4. Execute the following single-line command:

```bash
snakemake --cores all --keep-going

## Output

The pipeline will automatically generate the `valid_runs.txt` file containing the PDB names of all successfully processed files, and the final Kernel Density Estimate plot.

