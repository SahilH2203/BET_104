import gzip
import multiprocessing as mp
import os
import sys
from functools import partial

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from tqdm import tqdm


AA_VOLUME_CATEGORIES = {
    "K": "Large", "R": "Bulky", "D": "Intermediate", "E": "Large",
    "G": "Tiny", "A": "Tiny", "V": "Small", "L": "Intermediate", "I": "Intermediate",
    "M": "Large", "P": "Small", "S": "Small", "T": "Small", "N": "Intermediate",
    "Q": "Large", "C": "Small", "F": "Bulky", "Y": "Bulky", "W": "Bulky", "H": "Large",
}


def calc_spatial_angle(vec_a, vec_b, rot_axis, return_deg=True):
    vec_a = np.array(vec_a, dtype=float)
    vec_b = np.array(vec_b, dtype=float)
    rot_axis = np.array(rot_axis, dtype=float)

    vec_a /= np.linalg.norm(vec_a)
    vec_b /= np.linalg.norm(vec_b)
    rot_axis /= np.linalg.norm(rot_axis)

    cross_prod = np.cross(vec_a, vec_b)
    dot_prod = np.dot(vec_a, vec_b)

    angle_rad = np.arctan2(
        np.dot(cross_prod, rot_axis),
        dot_prod
    )

    return np.degrees(angle_rad) if return_deg else angle_rad


def parse_gzipped_pdb(prot_id, dir_path):
    file_loc = os.path.join(dir_path, f"{prot_id}.pdb.gz")
    if not os.path.exists(file_loc):
        return None
    struct_parser = PDBParser(QUIET=True)
    try:
        with gzip.open(file_loc, "rt") as archive:
            return struct_parser.get_structure(prot_id, archive)
    except Exception:
        return None


def locate_amino_acid(model_struct, target_chain, seq_id):
    try:
        for sub_model in model_struct:
            chain_data = sub_model[target_chain]
            for aa in chain_data:
                if aa.id[1] == seq_id:
                    return aa
    except Exception:
        return None
    return None


def fetch_alpha_carbon_coord(aa_obj):
    if aa_obj and "CA" in aa_obj:
        return aa_obj["CA"].get_coord()
    return None


def compute_sidechain_center(aa_obj):
    atom_pts = []
    for atm in aa_obj:
        if atm.get_name() not in ["N", "CA", "C", "O"]:
            atom_pts.append(atm.get_coord())
    if not atom_pts:
        return None
    return np.mean(atom_pts, axis=0)


def evaluate_local_geometry(tsv_file, input_dir, archive_dir, focus_res):
    """Parses context TSV, computes vectors, and returns structured data."""
    full_path = os.path.join(input_dir, tsv_file)

    if os.path.getsize(full_path) == 0:
        return ([], None)

    try:
        tsv_data = pd.read_csv(full_path, sep="\t", header=None)
    except Exception:
        return ([], None)

    if tsv_data.empty:
        return ([], None)

    prot_id = tsv_data.iloc[0, 12]
    parsed_3d = parse_gzipped_pdb(prot_id, archive_dir)
    if parsed_3d is None:
        return ([], None)

    extracted_metrics = []
    for block_idx in range(0, len(tsv_data), 3):
        try:
            row_left = tsv_data.iloc[block_idx]
            row_mid = tsv_data.iloc[block_idx + 1]
            row_right = tsv_data.iloc[block_idx + 2]
        except Exception:
            continue

        if row_mid[0] != focus_res:
            continue
        if row_mid[11] != "HHH":
            continue

        chain_letter = row_mid[1]
        idx_left = int(row_left[2])
        idx_mid = int(row_mid[2])
        idx_right = int(row_right[2])

        aa_left = locate_amino_acid(parsed_3d, chain_letter, idx_left)
        aa_mid = locate_amino_acid(parsed_3d, chain_letter, idx_mid)
        aa_right = locate_amino_acid(parsed_3d, chain_letter, idx_right)
        
        if not aa_left or not aa_mid or not aa_right:
            continue

        ca_left = fetch_alpha_carbon_coord(aa_left)
        ca_mid = fetch_alpha_carbon_coord(aa_mid)
        ca_right = fetch_alpha_carbon_coord(aa_right)
        
        sc_center_left = compute_sidechain_center(aa_left)
        sc_center_mid = compute_sidechain_center(aa_mid)

        if any(val is None for val in [ca_left, ca_mid, ca_right, sc_center_left, sc_center_mid]):
            continue

        vec_prior = sc_center_left - ca_left
        vec_center = sc_center_mid - ca_mid
        axis_rotation = ca_mid - ca_left
        
        if np.linalg.norm(axis_rotation) == 0:
            continue

        try:
            computed_angle = calc_spatial_angle(vec_prior, vec_center, axis_rotation)
        except Exception:
            continue

        neighbor_code = row_left[9]
        size_label = AA_VOLUME_CATEGORIES.get(neighbor_code, "Unknown")
        if size_label == "Unknown":
            continue

        extracted_metrics.append([prot_id, neighbor_code, size_label, computed_angle])

    return (extracted_metrics, prot_id if extracted_metrics else None)


def main():
    folder_contexts = sys.argv[1]
    dest_file = sys.argv[2]
    aa_query = sys.argv[3]
    folder_pdbs = sys.argv[4] if len(sys.argv) > 4 else "pdbs"
    num_cores = int(sys.argv[5]) if len(sys.argv) > 5 else 12

    tsv_list = sorted(f for f in os.listdir(folder_contexts) if f.endswith(".tsv"))

    mp_worker = partial(
        evaluate_local_geometry,
        input_dir=folder_contexts,
        archive_dir=folder_pdbs,
        focus_res=aa_query,
    )

    compiled_results = []
    successful_ids = set()

    with mp.Pool(num_cores) as compute_pool:
        for batch_rows, p_id in tqdm(
            compute_pool.imap_unordered(mp_worker, tsv_list, chunksize=16),
            total=len(tsv_list),
            desc="Calculating 3D Angles",
        ):
            if batch_rows:
                compiled_results.extend(batch_rows)
            if p_id is not None:
                successful_ids.add(str(p_id))

    final_table = pd.DataFrame(compiled_results, columns=["pdb", "left_aa", "size_class", "angle"])
    out_folder = os.path.dirname(dest_file) or "."
    os.makedirs(out_folder, exist_ok=True)
    final_table.to_csv(dest_file, sep="\t", index=False)

    tracker_file = os.path.join(out_folder, "valid_pdbs.txt")
    with open(tracker_file, "w") as txt_out:
        for p_id in sorted(successful_ids):
            txt_out.write(f"{p_id}\n")

    print(f"Angle data written to: {dest_file}")
    print(f"Processed PDB list written to: {tracker_file}")
    print(f"Total calculated angles: {len(final_table)}")
    print(f"Total valid structures: {len(successful_ids)}")


if __name__ == "__main__":
    main()