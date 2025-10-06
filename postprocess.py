from pathlib import Path
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# --- Paths ---
tsv_path = Path("/Users/prahladbhat/Downloads/sabdab_dataset/sabdab_structures/sabdab_summary_all.tsv")
fasta_root = Path("/Users/prahladbhat/Downloads/sabdab_dataset/fasta")

# Standard FASTA directories
dirs = {
    "antigen": fasta_root / "antigen_only",
    "all_chains": fasta_root / "all_chains",
    "heavy": fasta_root / "heavy_only",
    "light": fasta_root / "light_only",
}

cdr_root = fasta_root / "cdrs"

# Output files
outputs = {
    "antigen": fasta_root / "all_antigens_combined.fasta",
    "all_chains": fasta_root / "all_chains_combined.fasta",
    "heavy": fasta_root / "heavy_only_combined.fasta",
    "light": fasta_root / "light_only_combined.fasta",
}


# Add framework directories
framework_dirs = {
    "framework_heavy": fasta_root / "framework_only" / "heavy",
    "framework_light": fasta_root / "framework_only" / "light",
}

# Add corresponding output paths
framework_outputs = {
    "framework_heavy": fasta_root / "framework_heavy_combined.fasta",
    "framework_light": fasta_root / "framework_light_combined.fasta",
}



cdr_types = ["H1","H2","H3","L1","L2","L3"]

# Read TSV
df = pd.read_csv(tsv_path, sep="\t")

# --- Metadata helper ---
def get_metadata(pdb_id, chain_type):
    try:
        tsv_match = df[df["pdb"].str.lower() == pdb_id]
        tsv_row = tsv_match.iloc[0] if not tsv_match.empty else None
        compound_lower = str(tsv_row.get("compound","")).lower() if tsv_row is not None else ""
        fragment = "FAB" if "fab" in compound_lower or "fragment" in compound_lower else (
                   "MAB" if "mab" in compound_lower or "monoclonal" in compound_lower else "None")
        nanobody = "nano" in compound_lower
        scfv = tsv_row.get("scfv", False) if tsv_row is not None else False
        source = tsv_row.get("organism","unknown") if tsv_row is not None else "unknown"
        antibody = True
        return fragment, nanobody, scfv, source, antibody, tsv_row
    except Exception as e:
        print(f"❌ Failed to get metadata for {pdb_id}: {e}")
        return "None", False, False, "unknown", True, None

# --- Process a single FASTA file (normal or CDR) ---
def process_fasta_file(fasta_file, seq_type, cdr=None):
    try:
        with open(fasta_file) as f:
            seq_lines = []
            header_line = None
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    header_line = line
                else:
                    seq_lines.append(line)
        seq = "".join(seq_lines)
        if not seq:
            print(f"⚠️ Empty sequence for {fasta_file}")
            return None

        # PDB ID from header or filename
        if header_line and ">" in header_line:
            pdb_id = header_line.split("_")[0].replace(">","").lower()
        else:
            pdb_id = fasta_file.stem.split("_")[0].lower()

        # Chain type: for normal, first letter of stem; for CDR, from CDR name
        chain_type = cdr[0] if cdr else fasta_file.stem[0].upper()

        fragment, nanobody, scfv, source, antibody, tsv_row = get_metadata(pdb_id, chain_type)

        # Antigen metadata
        if seq_type == "antigen" and tsv_row is not None:
            antigen_chains = [c.strip() for c in str(tsv_row.get("antigen_chain","")).split("|") if c.strip().lower() != "na"]
            antigen_names = [n.strip() for n in str(tsv_row.get("antigen_name","")).split("|") if n.strip()]
            if chain_type in antigen_chains:
                idx = antigen_chains.index(chain_type)
                name = antigen_names[idx]
                antigen_type = tsv_row.get("antigen_type","none")
            else:
                name = "none"
                antigen_type = "none"
        else:
            name = "none"
            antigen_type = "none"

        # CDR field
        cdr_field = cdr if cdr else "none"

        # Header
        header = (
            f">PDB={pdb_id}|CHAIN={chain_type}|TYPE={seq_type}|CDRs={cdr_field}|NAME={name}|SOURCE={source}"
            f"|SCFV={scfv}|ANTIBODY={antibody}|NANOBODY={nanobody}|FRAGMENT={fragment}|ANTIGEN_TYPE={antigen_type}"
        )

        return f"{header}\n{seq}\n\n"

    except Exception as e:
        print(f"❌ Failed processing {fasta_file}: {e}")
        return None

# --- Combine FASTAs (normal) ---
def combine_fasta_multithread(fasta_dir, output_file, seq_type, max_workers=8):
    print(f"\n--- Combining {seq_type} FASTAs into {output_file} ---")
    all_entries = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_fasta_file, f, seq_type): f for f in fasta_dir.glob("*.fasta")}
        for future in as_completed(futures):
            fasta_file = futures[future]
            try:
                entry = future.result()
                if entry:
                    all_entries.append(entry)
            except Exception as e:
                print(f"❌ Exception for {fasta_file}: {e}")
    with open(output_file, "w") as out_f:
        out_f.writelines(all_entries)
    print(f"✅ Wrote {len(all_entries)} sequences to {output_file}")

# --- Combine CDR FASTAs ---
def combine_cdrs_multithread(cdr_root, fasta_root, max_workers=8):
    for cdr in cdr_types:
        cdr_dir = cdr_root / cdr
        output_file = fasta_root / f"{cdr}_combined.fasta"
        print(f"\n--- Combining CDR {cdr} FASTAs into {output_file} ---")
        all_entries = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(process_fasta_file, f, "cdr", cdr): f for f in cdr_dir.glob("*.fasta")}
            for future in as_completed(futures):
                fasta_file = futures[future]
                try:
                    entry = future.result()
                    if entry:
                        all_entries.append(entry)
                except Exception as e:
                    print(f"❌ Exception for {fasta_file}: {e}")
        with open(output_file, "w") as out_f:
            out_f.writelines(all_entries)
        print(f"✅ Wrote {len(all_entries)} sequences to {output_file}")

def combine_framework_multithread(framework_dirs, framework_outputs, max_workers=8):
    for seq_type, dir_path in framework_dirs.items():
        output_file = framework_outputs[seq_type]
        print(f"\n--- Combining {seq_type} FASTAs into {output_file} ---")
        all_entries = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(process_fasta_file, f, seq_type): f for f in dir_path.glob("*.fasta")}
            for future in as_completed(futures):
                fasta_file = futures[future]
                try:
                    entry = future.result()
                    if entry:
                        all_entries.append(entry)
                except Exception as e:
                    print(f"❌ Exception for {fasta_file}: {e}")
        with open(output_file, "w") as out_f:
            out_f.writelines(all_entries)
        print(f"✅ Wrote {len(all_entries)} sequences to {output_file}")


# --- Run all ---
for seq_type, dir_path in dirs.items():
    combine_fasta_multithread(dir_path, outputs[seq_type], seq_type)

combine_cdrs_multithread(cdr_root, fasta_root)
combine_framework_multithread(framework_dirs, framework_outputs)

print("\n✅ All combined FASTAs (normal + CDR) completed.")
