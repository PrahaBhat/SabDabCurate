import pandas as pd
from pathlib import Path
from collections import OrderedDict

# --- 3-letter to 1-letter amino acid map ---
three_to_one = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
    'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
    'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
    'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
    'SEC':'U','PYL':'O'
}

# --- Chothia CDR boundaries ---
chothia_cdrs = {
    "H1": (26, 32),
    "H2": (52, 56),
    "H3": (95, 102),
    "L1": (24, 34),
    "L2": (50, 56),
    "L3": (89, 97),
}

chothia_frameworks = {
    "FR1": (1, 25),
    "FR2": (33, 51),
    "FR3": (57, 94),
    "FR4": (103, 113),
}


# --- Helper functions ---
def extract_sequence_with_positions(pdb_file):
    """
    Extract residues with numbers (including insertions) from ATOM records.
    Returns {chain_id: [(resnum, aa), ...]}
    """
    chain_residues = OrderedDict()
    seen = set()
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") and line[13:15].strip() == "CA":
                chain_id = line[21].strip()
                res_num = line[22:26].strip()  # includes insertions
                res_name = line[17:20].strip()
                res_id = (chain_id, res_num)
                if res_id not in seen:
                    seen.add(res_id)
                    if chain_id not in chain_residues:
                        chain_residues[chain_id] = []
                    chain_residues[chain_id].append((res_num, three_to_one.get(res_name, "X")))
    return chain_residues

def in_chothia_range(resnum_str, start, end):
    """Return True if residue number (with optional insertion) falls in range"""
    try:
        num = int(''.join([c for c in resnum_str if c.isdigit()]))
        return start <= num <= end
    except ValueError:
        return False

def mask_non_cdr_residues(residues, start, end):
    """Keep full length but dash out non-CDR residues"""
    return "".join(aa if in_chothia_range(resnum, start, end) else "-" for resnum, aa in residues)

def mask_cdr_residues(residues, cdr_boundaries):
    """Keep framework residues, dash out CDR residues"""
    return "".join(
        "-" if any(in_chothia_range(resnum, start, end) for start, end in cdr_boundaries) else aa
        for resnum, aa in residues
    )
def mask_non_fragment_residues(residues, fragment_ranges):
    """
    Keep full sequence length but dash out any residues
    not in the provided fragment ranges.
    """
    masked_seq = ""
    for resnum, aa in residues:
        num = int("".join([c for c in resnum if c.isdigit()]))
        if any(start <= num <= end for start, end in fragment_ranges):
            masked_seq += aa
        else:
            masked_seq += "-"
    return masked_seq

# === Settings ===
tsv_path = Path("sabdab_structures/sabdab_summary_all.tsv")
pdb_dir = Path("sabdab_structures/all_structures/chothia")
fasta_root = Path("fasta")

# Subdirectories
dirs = {
    "all_chains": fasta_root / "all_chains",
    "antigen": fasta_root / "antigen_only",
    "heavy": fasta_root / "heavy_only",
    "light": fasta_root / "light_only",
    "cdr": fasta_root / "cdrs",
    "framework_heavy": fasta_root / "frameworks/heavy",
    "framework_light": fasta_root / "frameworks/light",
}
for path in dirs.values():
    path.mkdir(parents=True, exist_ok=True)

for cdr in ["H1","H2","H3","L1","L2","L3"]:
    (dirs["cdr"] / cdr).mkdir(parents=True, exist_ok=True)

# === Read summary TSV ===
df = pd.read_csv(tsv_path, sep="\t")
unique_pdbs = df.drop_duplicates(subset=["pdb"])

# === Main Loop ===
for _, row in unique_pdbs.iterrows():
    pdb_id = row["pdb"].lower()
    pdb_file = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_file.exists():
        print(f"⚠️ {pdb_id}.pdb not found, skipping")
        continue

    chain_residues = extract_sequence_with_positions(pdb_file)

    h_chain = row.get("Hchain")
    l_chain = row.get("Lchain")
    antigen_chains = [c.strip() for c in str(row.get("antigen_chain","")).split("|") if c.strip().lower() != "na"]

    # --- Write full antibody FASTA ---
    with open(dirs["all_chains"] / f"{pdb_id}.fasta", "w") as f:
        if pd.notna(h_chain) and h_chain in chain_residues:
            seq = "".join(aa for _, aa in chain_residues[h_chain])
            f.write(f">{pdb_id}_H\n{seq}\n")
        if pd.notna(l_chain) and l_chain in chain_residues:
            seq = "".join(aa for _, aa in chain_residues[l_chain])
            f.write(f">{pdb_id}_L\n{seq}\n")

    # --- Heavy / Light FASTAs ---
    for chain, path, label in [(h_chain, dirs["heavy"], "H"), (l_chain, dirs["light"], "L")]:
        if pd.notna(chain) and chain in chain_residues:
            seq = "".join(aa for _, aa in chain_residues[chain])
            with open(path / f"{pdb_id}.fasta", "w") as f:
                f.write(f">{pdb_id}_{label}\n{seq}\n")

    # --- Antigen FASTA ---
    with open(dirs["antigen"] / f"{pdb_id}.fasta", "w") as f:
        for chain in antigen_chains:
            if chain in chain_residues:
                seq = "".join(aa for _, aa in chain_residues[chain])
                f.write(f">{pdb_id}_{chain}\n{seq}\n")

    # --- CDR FASTAs ---
    for chain, cdr_list, label in [(h_chain, ["H1","H2","H3"], "H"), (l_chain, ["L1","L2","L3"], "L")]:
        if pd.notna(chain) and chain in chain_residues:
            residues = chain_residues[chain]
            for cdr in cdr_list:
                start, end = chothia_cdrs[cdr]
                masked_seq = mask_non_cdr_residues(residues, start, end)
                with open(dirs["cdr"] / cdr / f"{pdb_id}_{cdr}.fasta", "w") as f:
                    f.write(f">{pdb_id}_{cdr}\n{masked_seq}\n")
                    
# --- Framework-only FASTAs (mask non-framework residues) ---
    for chain, label in [(h_chain, "H"), (l_chain, "L")]:
        if pd.notna(chain) and chain in chain_residues:
            residues = chain_residues[chain]

            # Chothia framework regions (FR1–FR4)
            framework_ranges = list(chothia_frameworks.values())

            masked_seq = mask_non_fragment_residues(residues, framework_ranges)

            out_dir = fasta_root / "framework_only" / ("heavy" if label == "H" else "light")
            out_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_dir / f"{pdb_id}_framework{label}.fasta"

            with open(out_path, "w") as f:
                f.write(f">{pdb_id}_framework{label}\n{masked_seq}\n")

    print(f"✅ Processed {pdb_id}")

print("\n🎉 All PDBs processed successfully!")
