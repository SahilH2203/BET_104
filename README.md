
# Target Amino Acid Neighborhood Effect Pipeline

This pipeline extracts the secondary structure context and calculates the signed 3D angle between the $C_\alpha$ -> Centroid vectors of adjacent residues surrounding an Arginine (X-R-X) within an Alpha Helix (HHH).

# Plot

<img width="3000" height="1800" alt="ss_profile_HHH_for_arg" src="https://github.com/user-attachments/assets/cfb6b646-1ad8-4646-904e-e9585d4d7dd7" />

## Instructions to Reproduce

1. Clone this repository.
2. Create a folder named `pdbs/` in the root directory and place all raw `.pdb` (or `.pdb.gz`) data files inside it.
3. Ensure you have the required dependencies installed (`snakemake`, `stride`, `biopython`, `numpy`, `pandas`, `seaborn`, `matplotlib`).
4. Execute the following single-line command:

\`\`\`bash
snakemake --cores all
\`\`\`

## Output

The pipeline will automatically generate the `valid_runs.txt` file containing the PDB names of all successfully processed files, and the final Kernel Density Estimate plot shown below:
<img width="3000" height="1800" alt="ss_profile_HHH_for_arg" src="https://github.com/user-attachments/assets/c7dfc1bd-05dc-453a-bf7a-dafc737f4f90" />

