configfile: "config.yaml"

import glob
import os

# Pull variables from config.yaml, with safe fallbacks
AA = config.get("target_aa", "ARG")
STRIDE = config.get("stride_path", "stride")
PLOT_NAME = config.get("plot_name", "ss_profile_HHH_for_arg_with_valid_runs.png")
ANGLE_FILE = config.get("angle_file", "computed_angles.tsv")
VALID_RUNS = config.get("valid_runs", "valid_pdbs.txt")

# Find all PDBs dynamically
PDBS = [
    os.path.basename(f).replace(".pdb.gz", "")
    for f in glob.glob("pdbs/*.pdb.gz")
]

rule all:
    input:
        f"result/{PLOT_NAME}",
        f"result/{VALID_RUNS}"

rule unzip_pdb:
    input:
        pdb = "pdbs/{pdb}.pdb.gz"
    output:
        unzipped = temp("unzipped_pdbs/{pdb}.pdb")
    shell:
        "zcat {input.pdb} > {output.unzipped}"

rule run_stride:
    input:
        unzipped = "unzipped_pdbs/{pdb}.pdb"
    output:
        ss_out = "stride_out/{pdb}.ss.out"
    params:
        stride = STRIDE
    shell:
        "{params.stride} {input.unzipped} > {output.ss_out} 2>/dev/null || touch {output.ss_out}"

rule extract_context:
    input:
        ss = "stride_out/{pdb}.ss.out"
    output:
        tsv = "contexts/context_for_{aa}_in_{pdb}.tsv"
    shell:
        "python Script/context_stride.py {input.ss} {output.tsv} {wildcards.aa}"

rule calculate_angles:
    input:
        tsvs = expand("contexts/context_for_{aa}_in_{pdb}.tsv", aa=[AA], pdb=PDBS)
    output:
        angles = f"result/{ANGLE_FILE}",
        valid_pdbs = f"result/{VALID_RUNS}"
    params:
        aa = AA,
        contexts_dir = "contexts",
        pdbs_dir = "pdbs"
    shell:
        "python Script/calc_angles.py {params.contexts_dir} {output.angles} {params.aa} {params.pdbs_dir}"

rule plot_angles:
    input:
        angles = f"result/{ANGLE_FILE}"
    output:
        plot = f"result/{PLOT_NAME}"
    shell:
        "python Script/plot_angles.py {input.angles} {output.plot}"
