from collections import defaultdict
import torch
from torch.utils.data import Dataset
import pandas as pd
import os
from build3 import chothia_frameworks, chothia_cdrs, three_to_one, is_residue_in_range

fasta_root = "/Users/zacharycohen/Desktop/GitHub/SabDabCurate/sabdabfasta/all_chains"
structure_root = "/Users/zacharycohen/Desktop/GitHub/SabDabCurate/sabdab_structures/all_structures/chothia"

def parse_multiple_paired_remarks(text):
    result = {}
    for line in text.strip().splitlines():
        parts = line.strip().split()
        for p in parts:
            if "=" in p:
                field, val = p.split("=")
                field = field.strip()
                val = val.strip()
                # We only care about HCHAIN or LCHAIN assignments
                if field in ("HCHAIN", "LCHAIN"):
                    if val in result:
                        raise ValueError(f"Duplicate key: {val}")
                    if field not in result.values(): result[val] = field
    return result

def parse_single_ca_line(line):
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        raise ValueError("Line is not a valid ATOM/HETATM line: " + line)

    atom_name = line[12:16].strip()
    if atom_name != "CA":
        raise ValueError("Line is not a valid CA atom line: " + line)

    return {
        'seq': line[21].strip().upper(),          # chain ID
        'indice': line[22:26].strip(),
        'x': float(line[30:38].strip()),          # x coordinate
        'y': float(line[38:46].strip()),          # y coordinate
        'z': float(line[46:54].strip()),          # z coordinate
    }

class SabDabDatasetV2(Dataset):

    def load_fasta(self, filepath):
        sequences = {}
        id = filepath.split("/")[-1].split(".")[0]
        with open(filepath, 'r') as file:
            lines = file.readlines()
            for idx, line in enumerate(lines):
                if line.startswith('>'):
                    tag = line.strip().split('_')[1] # L or H
                    sequences[tag] = lines[idx+1].strip()
        return sequences, id

    def compute_distogram(self, coords, n_bins=64, start_a=2.0, end_a=22.0):
        coords_tensor = torch.tensor(coords)
        distances = torch.cdist(coords_tensor, coords_tensor)  # (N, N)
        bin_edges = torch.linspace(start_a, end_a, n_bins)
        bin_idx = torch.bucketize(distances, bin_edges) - 1
        bin_idx = torch.clamp(bin_idx, 0, n_bins - 1)

        distogram = torch.nn.functional.one_hot(bin_idx, num_classes=n_bins).float()

        return distogram
    
    def get_chothia_cdr(self, resnum):
        for key, val in chothia_cdrs.items():
            if is_residue_in_range(resnum, val[0], val[1]):
                return key
        return None
            

    def load_coords(self, filepath):
        resnums = defaultdict(list)
        coords = defaultdict(list)
        cdr_mask = defaultdict(list)
        keys = defaultdict(list)
        with open(filepath, 'r') as file:
            lines = file.readlines()
            key_lines = [line for line in lines if line.startswith("REMARK   5 PAIRED_HL") or line.startswith("REMARK   5 SINGLE")]
            keys = parse_multiple_paired_remarks("".join(key_lines))
            for line in lines: # ONLY PULLING 
                if line[13:15].strip() == "CA": # only timesaves a little
                    try:
                        res = parse_single_ca_line(line)
                    except ValueError as e:
                        continue
                    if res['seq'] in keys: 
                        resnum = res['indice']
                        coords[keys[res['seq']]].append((res['x'], res['y'], res['z']))
                        resnums[keys[res['seq']]].append(resnum)
                        if any(is_residue_in_range(resnum, start, end) for start, end in chothia_frameworks.values()):
                            cdr_mask[keys[res['seq']]].append(0)
                        else:
                            cdr = self.get_chothia_cdr(resnum)
                            if cdr:
                                cdr_mask[keys[res['seq']]].append(int(cdr.replace("H","").replace("L","")))
                            else:
                                cdr_mask[keys[res['seq']]].append(-1)
        return coords, cdr_mask

    def group_by_mask(self, coords, mask):
        grouped = defaultdict(list)
        for coord, m in zip(coords, mask):
            grouped[m].append(coord)
        return grouped

    def load_structure(self, filepath):
        coords, cdr_mask = self.load_coords(filepath)
        # assume cdr_mask is a dict of lists, where the keys are the same as coords and the values are lists of 0s and 1,2,3 for H1,H2,H3,L1,L2,L3

        distograms = {}
        for chain, coord_set in coords.items(): # chain is H or L
            mask = cdr_mask[chain]
            grouped_by_cdr = self.group_by_mask(coord_set, mask)
            for cdr, coord_set in grouped_by_cdr.items():
                if cdr in (1,2,3):
                    distogram = self.compute_distogram(coord_set).to(self.device)
                    distograms[f"{chain[0]}{cdr}"] = distogram
                    print(f"[Wrote] distogram for {chain[0]}{cdr}")
        # we want to return a set of distogram which is a dict of h1, h2, h3, l1, l2, l3:dgram mappings
        return distograms

    def __init__(self, fasta_dir, structure_dir, device='cpu'):
        self.fasta_dir = fasta_dir
        self.structure_dir = structure_dir
        self.data = []
        self.device = device
        fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith(".fasta")]

        for fasta_file in fasta_files:
            pdb_id = fasta_file.split(".")[0]
            fasta_path = os.path.join(fasta_dir, fasta_file)
            pdb_path = os.path.join(structure_dir, f"{pdb_id}.pdb")

            if not os.path.exists(pdb_path):
                print(pdb_path)
                print(f"[Warning] PDB file not found for {pdb_id}, skipping.")
                continue

            sequences, seq_id = self.load_fasta(fasta_path)


            try:
                distograms = self.load_structure(pdb_path)
            except Exception as e:
                print(f"[Error] Failed to load structure for {pdb_id}: {e}")
                # raise(e)

                continue
            distograms = {k.replace("HCHAIN", "H").replace("LCHAIN", "L"): v for k, v in distograms.items()}

            #[REMOVED] sanity check: chains match between FASTA and PDB
            # missing_chains = [ch for ch in sequences.keys() if ch not in distograms]
            # if missing_chains:
            #     print(f"[Warning] Missing distograms for chains {missing_chains} in {pdb_id}")
            #     print(distograms.keys())



            # for k, v in distograms.items(): # [REMOVED B/C ALREADY SLICING TO CHOTHIA] ASSUMPTION BASED: Slice because PDB includes extra residues
            #     if v.shape[0] != len(sequences[k]) or v.shape[1] != len(sequences[k]):
            #         print(f"[Trimmed] {pdb_id} {k}: {v.shape} -> ({len(sequences[k])}, {len(sequences[k])}, {v.shape[2]})")
            #         L = len(sequences[k])
            #         distograms[k] = v[:L, :L, :]


            self.data.append({
                "pdb_id": pdb_id,
                "seq_id": seq_id,
                "sequences": sequences,         
                "distograms": distograms        
            })

        print(f"[INFO] Loaded {len(self.data)} antibody structures.")
        

        

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx]

test_dataset = SabDabDatasetV2(fasta_root, structure_root) # DO WITH SMALL DATASET FIRST