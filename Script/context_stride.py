import os
import re
import sys
import pandas as pd
import multiprocessing as mp
from functools import partial
from tqdm import tqdm

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

def process_single_stride(file_name, input_folder, output_folder, focus_residue):
    """Parses a single STRIDE file and saves the context TSV."""
    stride_input = os.path.join(input_folder, file_name)
    protein_id = file_name.split(".")[0]
    tsv_output = os.path.join(output_folder, f"{protein_id}.tsv")

    parsed_data = []

    try:
        with open(stride_input, 'r') as file_handle:
            for line in file_handle:
                if line.startswith("ASG"):
                    tokens = line.split()
                    if len(tokens) < 10:
                        continue

                    num_match = re.match(r"\d+", tokens[3])
                    if not num_match:
                        continue

                    parsed_data.append({
                        "aa_3": tokens[1],
                        "chain_id": tokens[2],
                        "seq_num": int(num_match.group()),
                        "index": int(tokens[4]),
                        "ss_char": tokens[5],
                        "ss_desc": tokens[6],
                        "ang_phi": float(tokens[7]),
                        "ang_psi": float(tokens[8]),
                        "surf_area": float(tokens[9]),
                    })
    except Exception:
        return

    struct_df = pd.DataFrame(parsed_data)
    final_rows = []

    if not struct_df.empty:
        for idx in range(1, len(struct_df) - 1):
            mid_res = struct_df.iloc[idx]
            
            if mid_res["aa_3"] != focus_residue:
                continue

            left_res = struct_df.iloc[idx - 1]
            right_res = struct_df.iloc[idx + 1]

            seq_triplet = (
                THREE_TO_ONE.get(left_res["aa_3"], "X") +
                THREE_TO_ONE.get(mid_res["aa_3"], "X") +
                THREE_TO_ONE.get(right_res["aa_3"], "X")
            )
            
            struct_triplet = left_res["ss_char"] + mid_res["ss_char"] + right_res["ss_char"]

            position_mapping = (
                f"{mid_res['chain_id']}:{left_res['seq_num']},"
                f"{mid_res['chain_id']}:{mid_res['seq_num']},"
                f"{mid_res['chain_id']}:{right_res['seq_num']}"
            )

            for res_block in [left_res, mid_res, right_res]:
                final_rows.append([
                    res_block["aa_3"],
                    res_block["chain_id"],
                    res_block["seq_num"],
                    res_block["index"],
                    res_block["ss_char"],
                    res_block["ss_desc"],
                    res_block["ang_phi"],
                    res_block["ang_psi"],
                    res_block["surf_area"],
                    THREE_TO_ONE.get(res_block["aa_3"], "X"),
                    seq_triplet,
                    struct_triplet,
                    protein_id,
                    position_mapping,
                ])

    if final_rows:
        export_df = pd.DataFrame(final_rows)
        export_df.to_csv(tsv_output, sep="\t", index=False, header=False)


def parse_stride_contexts_batch():
    stride_dir = sys.argv[1]
    out_dir = sys.argv[2]
    focus_residue = sys.argv[3]
    cores = int(sys.argv[4]) if len(sys.argv) > 4 else 12

    os.makedirs(out_dir, exist_ok=True)
    
    # Get all .ss.out files in the directory
    files = [f for f in os.listdir(stride_dir) if f.endswith(".ss.out")][:49000]
    #files = [f for f in os.listdir(stride_dir) if f.endswith(".ss.out")]

    print(f"Found {len(files)} STRIDE files. Processing on {cores} cores...")

    # Set up the worker function with fixed arguments
    worker = partial(
        process_single_stride, 
        input_folder=stride_dir, 
        output_folder=out_dir, 
        focus_residue=focus_residue
    )

    # Run multiprocessing pool with a progress bar
    with mp.Pool(cores) as pool:
        list(tqdm(
            pool.imap_unordered(worker, files, chunksize=100), 
            total=len(files), 
            desc="Extracting Contexts"
        ))

    print("Extraction complete!")


if __name__ == "__main__":
    parse_stride_contexts_batch()