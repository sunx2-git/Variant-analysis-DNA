import sqlite3
import pickle
from pathlib import Path
import numpy as np

# === Paths ===
output_dir = Path("output1")
cache_file = output_dir / "variant_cache.pkl"
db_path = output_dir / "variant_data.db"

# === Load cache ===
with cache_file.open("rb") as f:
    variant_data = pickle.load(f)

print(f"📦 Loaded {len(variant_data)} variants from cache.")

# === Connect DB ===
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# === Reset tables ===
cursor.executescript("""
DROP TABLE IF EXISTS variants;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS allele_fractions;
""")

cursor.execute("""
CREATE TABLE variants (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol TEXT,
    chrom TEXT,
    position INTEGER,
    end_position INTEGER,
    transcript_variant TEXT,
    protein_variant TEXT,
    gnomad_frequency REAL,
    artifact_score REAL,
    artifact_flag INTEGER
)
""")

cursor.execute("""
CREATE TABLE samples (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_name TEXT UNIQUE
)
""")

cursor.execute("""
CREATE TABLE allele_fractions (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    variant_id INTEGER,
    sample_id INTEGER,
    allele_fraction REAL,
    classification TEXT,
    FOREIGN KEY (variant_id) REFERENCES variants(id),
    FOREIGN KEY (sample_id) REFERENCES samples(id)
)
""")

# === Canonical sample set (MATCHES STREAMLIT) ===
def normalize_sample(name):
    import re
    m = re.search(r"(\d{4}D)", name)
    return m.group(1) if m else name

all_samples = set()
for rec in variant_data.values():
    for s in rec.get("Sample Files", []):
        all_samples.add(normalize_sample(s))

total_samples = len(all_samples)
print(f"🧮 Total unique samples: {total_samples}")

# === Insert data ===
sample_id_map = {}
variant_inserted = 0
allele_inserted = 0

for record in variant_data.values():
    try:
        chrom = record["Chromosome"]
        pos = int(record["Position"])
        end_pos = int(record["End Position"])
        gene = record["Gene Symbol"]
        transcript = record["Transcript Variant"]
        protein = record["Protein Variant"]

        gnomad = record.get("gnomAD Frequency")
        gnomad = float(gnomad) if gnomad not in ("", None) else None

        # --- AF processing (NORMALIZED) ---
        af_values = []
        for v in record.get("Allele Fractions", {}).values():
            try:
                af = float(v)
                if af > 1:
                    af /= 100
                af_values.append(af)
            except:
                continue

        recurrence = len(set(normalize_sample(s) for s in record.get("Sample Files", [])))
        mean_af = np.mean(af_values) if af_values else 0.0
        sd_af = min(np.std(af_values), 1.0) if af_values else 0.0

        artifact_score = round((recurrence / total_samples) * (1 - sd_af), 4)
        artifact_flag = int(recurrence / total_samples >= 0.5)

        cursor.execute("""
            INSERT INTO variants (
                gene_symbol, chrom, position, end_position,
                transcript_variant, protein_variant, gnomad_frequency,
                artifact_score, artifact_flag
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            gene, chrom, pos, end_pos,
            transcript, protein, gnomad,
            artifact_score, artifact_flag
        ))

        variant_id = cursor.lastrowid
        variant_inserted += 1

        # --- Insert per-sample data ---
        for sample_name in record.get("Sample Files", []):
            sample_key = normalize_sample(sample_name)

            if sample_key not in sample_id_map:
                cursor.execute(
                    "INSERT OR IGNORE INTO samples (sample_name) VALUES (?)",
                    (sample_key,)
                )
                cursor.execute(
                    "SELECT id FROM samples WHERE sample_name = ?",
                    (sample_key,)
                )
                sample_id_map[sample_key] = cursor.fetchone()[0]

            af = record.get("Allele Fractions", {}).get(sample_name)
            try:
                af = float(af)
                if af > 1:
                    af /= 100
            except:
                af = None

            cl = str(record.get("Classifications", {}).get(sample_name, "")).strip()

            cursor.execute("""
                INSERT INTO allele_fractions
                (variant_id, sample_id, allele_fraction, classification)
                VALUES (?, ?, ?, ?)
            """, (variant_id, sample_id_map[sample_key], af, cl))

            allele_inserted += 1

    except Exception as e:
        print(f"⚠️ Skipped record → {e}")

conn.commit()
conn.close()

print(f"✅ Variants inserted: {variant_inserted}")
print(f"📥 Allele rows inserted: {allele_inserted}")
