from pathlib import Path

cdr_root = Path("/Users/prahladbhat/Downloads/sabdab_dataset/fasta/cdrs")
output_dir = Path("/Users/prahladbhat/Downloads/sabdab_dataset/fasta")

heavy_output = output_dir / "heavy_cdrs_merged.fasta"
light_output = output_dir / "light_cdrs_merged.fasta"

# Define CDRs per chain
cdr_map = {
    "H": ["H1","H2","H3"],
    "L": ["L1","L2","L3"]
}

def merge_cdrs_positions(chain_type, output_file):
    # collect all pdbs
    pdb_set = set()
    for cdr in cdr_map[chain_type]:
        for f in (cdr_root / cdr).glob("*.fasta"):
            pdb_id = f.stem.split("_")[0]
            pdb_set.add(pdb_id)

    with open(output_file, "w") as out_f:
        for pdb in sorted(pdb_set):
            # Initialize merged_seq as empty, will grow to full length automatically
            merged_seq = []
            # We'll assume all CDRs have the same length sequences with dashes in non-CDR positions
            # Use first existing CDR to get the length
            first_cdr_file = None
            for cdr in cdr_map[chain_type]:
                f = cdr_root / cdr / f"{pdb}_{cdr}.fasta"
                if f.exists():
                    first_cdr_file = f
                    break
            if first_cdr_file is None:
                continue
            with open(first_cdr_file) as f:
                seq_len = len("".join([l.strip() for l in f if not l.startswith(">")]))
                merged_seq = ["-"]*seq_len

            # Overlay each CDR onto merged_seq
            for cdr in cdr_map[chain_type]:
                cdr_file = cdr_root / cdr / f"{pdb}_{cdr}.fasta"
                if cdr_file.exists():
                    with open(cdr_file) as f:
                        cdr_seq = "".join([l.strip() for l in f if not l.startswith(">")])
                    for i, aa in enumerate(cdr_seq):
                        if aa != "-":
                            merged_seq[i] = aa

            # Write to FASTA
            header = f">{pdb}|CHAIN={chain_type}|TYPE=cdrs_merged|CDRs={','.join(cdr_map[chain_type])}|NAME=none"
            out_f.write(f"{header}\n{''.join(merged_seq)}\n\n")
    print(f"✅ Wrote {len(pdb_set)} sequences to {output_file}")

# Merge heavy and light
merge_cdrs_positions("H", heavy_output)
merge_cdrs_positions("L", light_output)
