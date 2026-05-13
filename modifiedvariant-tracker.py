import pandas as pd
from pathlib import Path
import pickle
from collections import defaultdict
import statistics

# === Paths ===
input_dir = Path("input1")
output_dir = Path("output1")
output_dir.mkdir(parents=True, exist_ok=True)

cache_file = output_dir / "variant_cache.pkl"
summary_file_pickle = output_dir / "summary_df.pkl"
summary_csv = output_dir / "variant_presence_summary.csv"
comparison_csv = output_dir / "variant_comparison_50percent_plus.csv"

# === Load cache ===
variant_dict = {}
processed_samples = set()

if cache_file.exists():
    with open(cache_file, "rb") as f:
        variant_dict = pickle.load(f)
    processed_samples = {Path(name).stem for v in variant_dict.values() for name in v.get("Sample Files", set())}
    print(f"📦 Loaded cached variant data ({len(processed_samples)} samples)")

# === Identify new files ===
all_files = list(input_dir.glob("*.csv"))
new_files = [f for f in all_files if Path(f).stem not in processed_samples]
total_files = len({Path(f).stem for f in all_files})
print(f"📁 Found {len(all_files)} files. Processing {len(new_files)} new...")

# === Process input files ===
for file in new_files:
    sample_name = file.name
    sample_base = Path(sample_name).stem

    if file.stat().st_size == 0:
        print(f"⚠️ Skipping empty file: {sample_name}")
        continue

    df = pd.read_csv(file)
    df.columns = df.columns.str.strip()

    required_cols = ['Chromosome', 'Position', 'End Position', 'Gene Symbol',
                     'Transcript Variant', 'Protein Variant']
    allele_fraction_col = next((col for col in df.columns if col.lower().strip().endswith("allele fraction")), None)
    classification_col = next((col for col in df.columns if col.lower().strip().endswith("classification")), None)
    gnomad_col = next((col for col in df.columns if "gnomad" in col.lower()), None)

    missing = [col for col in required_cols if col not in df.columns]
    if missing or not allele_fraction_col or not classification_col:
        print(f"⚠️ Skipping {sample_name}, missing: {missing}")
        continue

    df['Position'] = pd.to_numeric(df['Position'], errors='coerce')
    df['End Position'] = pd.to_numeric(df['End Position'], errors='coerce')
    df.dropna(subset=['Chromosome', 'Position', 'End Position', 'Gene Symbol'], inplace=True)

    for _, row in df.iterrows():
        chrom = str(row['Chromosome']).strip()
        pos = int(row['Position'])
        end = int(row['End Position'])
        gene = str(row['Gene Symbol']).strip()

        # Safe parsing for possibly non-string values
        transcript_str = str(row.get('Transcript Variant', '')).strip("; ")
        protein_str = str(row.get('Protein Variant', '')).strip("; ")
        transcripts = [t.strip() for t in transcript_str.split(";") if t.strip()]
        proteins = [p.strip() for p in protein_str.split(";") if p.strip()]
        gnomad = str(row.get(gnomad_col, "")).strip() if gnomad_col else ""

        max_len = max(len(transcripts), len(proteins))
        transcripts += [""] * (max_len - len(transcripts))
        proteins += [""] * (max_len - len(proteins))

        af = row.get(allele_fraction_col)
        classification = row.get(classification_col)

        for t, p in zip(transcripts, proteins):
            full_key = f"{chrom}_{pos}_{end}_{gene}_{t}_{p}"
            summary_key = f"{chrom}_{pos}_{end}_{gene}"

            if full_key not in variant_dict:
                variant_dict[full_key] = {
                    "Chromosome": chrom,
                    "Position": pos,
                    "End Position": end,
                    "Gene Symbol": gene,
                    "Transcript Variant": t,
                    "Protein Variant": p,
                    "Sample Files": set(),
                    "Allele Fractions": {},
                    "Classifications": {},
                    "gnomAD Frequency": gnomad,
                    "Summary Key": summary_key
                }

            variant_dict[full_key]["Sample Files"].add(sample_name)
            variant_dict[full_key]["Allele Fractions"][sample_name] = str(af).strip() if pd.notna(af) else ""
            variant_dict[full_key]["Classifications"][sample_name] = str(classification).strip() if pd.notna(classification) else ""

    processed_samples.add(sample_base)

# === Generate summary ===
summary_group = defaultdict(lambda: {
    "Transcript Variants": [],
    "Protein Variants": [],
    "Sample Files": set(),
    "Classifications": set(),
    "gnomAD Frequency": ""
})

for v in variant_dict.values():
    sk = v["Summary Key"]
    sg = summary_group[sk]
    sg.update({
        "Chromosome": v["Chromosome"],
        "Position": v["Position"],
        "End Position": v["End Position"],
        "Gene Symbol": v["Gene Symbol"],
        "gnomAD Frequency": v["gnomAD Frequency"] or sg["gnomAD Frequency"]
    })
    if v["Transcript Variant"]:
        sg["Transcript Variants"].append(v["Transcript Variant"])
    if v["Protein Variant"]:
        sg["Protein Variants"].append(v["Protein Variant"])
    sg["Sample Files"].update(v["Sample Files"])
    sg["Classifications"].update(c for c in v["Classifications"].values() if c)

summary_rows = []
for sk, info in summary_group.items():
    sample_ids = sorted(Path(s).stem for s in info["Sample Files"])
    presence_pct = (len(sample_ids) / total_files) * 100

    max_pairs = max(len(info["Transcript Variants"]), len(info["Protein Variants"]))
    paired_transcripts = info["Transcript Variants"] + [""] * (max_pairs - len(info["Transcript Variants"]))
    paired_proteins = info["Protein Variants"] + [""] * (max_pairs - len(info["Protein Variants"]))
    transcript_str = "; ".join(paired_transcripts)
    protein_str = "; ".join(paired_proteins)

    summary_rows.append({
        "Variant Key": sk,
        "Chromosome": info["Chromosome"],
        "Position": info["Position"],
        "End Position": info["End Position"],
        "Gene Symbol": info["Gene Symbol"],
        "Transcript Variants": transcript_str,
        "Protein Variants": protein_str,
        "Classifications": "; ".join(sorted(str(c) for c in info["Classifications"] if pd.notna(c))),
        "gnomAD Frequency": info["gnomAD Frequency"],
        "Files Containing Variant": len(sample_ids),
        "Total Samples": total_files,
        "Presence %": round(presence_pct, 2),
        "Sample Names": ", ".join(sample_ids)
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.sort_values(by="Presence %", ascending=False, inplace=True)
summary_df.reset_index(drop=True, inplace=True)
summary_df.to_csv(summary_csv, index=False)
summary_df.to_pickle(summary_file_pickle)
print(f"✅ Saved summary: {summary_csv}")

# === Comparison table for 50%+ ===
summary_to_keys = defaultdict(list)
for fk, v in variant_dict.items():
    summary_to_keys[v["Summary Key"]].append(fk)

def get_af_stats(keys):
    afs = []
    for k in keys:
        for af in variant_dict[k]["Allele Fractions"].values():
            try:
                afs.append(float(af))
            except:
                continue
    if afs:
        return round(statistics.mean(afs), 4), round(statistics.stdev(afs), 4) if len(afs) > 1 else 0.0
    return None, None

comparison_rows = []
for _, row in summary_df[summary_df["Presence %"] >= 50].iterrows():
    sk = row["Variant Key"]
    mean_af, sd_af = get_af_stats(summary_to_keys[sk])
    comparison_rows.append({
        **row,
        "Mean Allele Fraction": mean_af,
        "SD Allele Fraction": sd_af
    })

pd.DataFrame(comparison_rows).to_csv(comparison_csv, index=False)
print(f"📊 Saved comparison: {comparison_csv}")

# === Final cache ===
with open(cache_file, "wb") as f:
    pickle.dump(variant_dict, f)

print("💾 Variant cache updated.")