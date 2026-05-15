import streamlit as st
import sqlite3
import pandas as pd

DB_PATH = "output1/variant_data.db"

# === Load Data from SQLite ===
@st.cache_data(ttl=3600)
def load_data():
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql("""
        SELECT
            v.id AS variant_id,
            v.gene_symbol,
            v.chrom,
            v.position,
            v.end_position,
            v.transcript_variant,
            v.protein_variant,
            v.gnomad_frequency,
            v.artifact_score,
            v.artifact_flag,
            af.sample_id,
            s.sample_name,
            af.allele_fraction,
            af.classification
        FROM variants v
        LEFT JOIN allele_fractions af ON v.id = af.variant_id
        LEFT JOIN samples s ON af.sample_id = s.id
    """, conn)
    conn.close()

    # Extract canonical sample ID
    df["normalized_sample_name"] = df["sample_name"].str.extract(r"(\d{4}D)")

    # Artifact labels
    df["artifact_label"] = df["artifact_flag"].map({1: "Likely Artifact", 0: "Not Artifact"})
    df["artifact_likelihood_percent"] = df["artifact_flag"] * 100

    return df

# === App Setup ===
st.set_page_config(page_title="Variant Viewer", layout="wide")
st.title("Variant Viewer")

df = load_data()
total_samples = df["normalized_sample_name"].nunique()

# === Sidebar Navigation ===
st.sidebar.header("Navigation")
view = st.sidebar.radio("Choose View", [
    "Variant → Sample Lookup",
    "Gene Summary",
    "Multi-Gene Sample Search",
    "Sample → Variant Lookup",
    "Artifact Scoring"
])

# ===================================================================
# VIEW 1: Variant → Sample Lookup
# ===================================================================
if view == "Variant → Sample Lookup":
    st.subheader("Variant → Sample Lookup")

    genes = sorted(df["gene_symbol"].dropna().unique())
    gene = st.selectbox("Select Gene", genes)
    df_gene = df[df["gene_symbol"] == gene]

    transcripts = sorted(df_gene["transcript_variant"].dropna().unique())
    transcript = st.selectbox("Select Transcript Variant", transcripts)

    proteins = sorted(df_gene[df_gene["transcript_variant"] == transcript]["protein_variant"].dropna().unique())
    protein = st.selectbox("Select Protein Variant", proteins)

    result = df_gene[
        (df_gene["transcript_variant"] == transcript) &
        (df_gene["protein_variant"] == protein)
    ][["sample_name", "normalized_sample_name", "allele_fraction", "classification", "gnomad_frequency"]]

    result = result.drop_duplicates(subset=["normalized_sample_name"])
    gnomad_freq = result["gnomad_frequency"].dropna().unique()
    gnomad_freq = gnomad_freq[0] if len(gnomad_freq) > 0 else "N/A"

    st.markdown(f"### Variant: {gene} / {transcript} / {protein}")
    st.markdown(f"- gnomAD Frequency: {gnomad_freq}")

    if not result.empty:
        result["allele_fraction"] = pd.to_numeric(result["allele_fraction"], errors="coerce").round(3)
        st.dataframe(result.drop(columns=["normalized_sample_name"]), use_container_width=True)

        csv = result.drop(columns=["normalized_sample_name"]).to_csv(index=False)
        st.download_button("Download CSV", data=csv, file_name="variant_sample_details.csv", mime="text/csv")
    else:
        st.info("No sample data for this variant.")

# ===================================================================
# VIEW 2: Gene Summary
# ===================================================================
elif view == "Gene Summary":
    st.subheader("Gene-Level Summary")
    total_samples = df["normalized_sample_name"].nunique()
    genes = sorted(df["gene_symbol"].dropna().unique())
    gene = st.selectbox("Select Gene", genes)
    df_gene = df[df["gene_symbol"] == gene]

    if df_gene.empty:
        st.warning("No variants found for this gene.")
    else:
        summary = (
            df_gene
            .drop_duplicates(subset=["variant_id", "normalized_sample_name"])
            .groupby(["variant_id", "transcript_variant", "protein_variant", "gnomad_frequency"], dropna=False)
            .agg(
                sample_count=("normalized_sample_name", "nunique"),
                mean_vaf=("allele_fraction", lambda x: pd.to_numeric(x, errors="coerce").mean()),
                sd_vaf=("allele_fraction", lambda x: pd.to_numeric(x, errors="coerce").std()),
                classifications=("classification", lambda x: "; ".join(sorted(set(str(c) for c in x if pd.notna(c)))))
            )
            .reset_index()
        )

        summary["presence_percent"] = round(summary["sample_count"] / total_samples * 100, 2)
        summary["mean_vaf"] = summary["mean_vaf"].round(3)
        summary["sd_vaf"] = summary["sd_vaf"].round(3)

        # --- SEARCH FUNCTIONALITY ---
        search_term = st.text_input(
            "🔍 Search variants (c., p., transcript)",
            placeholder="e.g. c.123A>G, p.G12D, NM_"
        ).strip().lower()

        if search_term:
            summary["search_hit"] = (
                summary["transcript_variant"].str.lower().str.contains(search_term, na=False) |
                summary["protein_variant"].str.lower().str.contains(search_term, na=False) |
                summary["classifications"].str.lower().str.contains(search_term, na=False)
            ).astype(int)
        else:
            summary["search_hit"] = 0

        # Sort hits to top
        summary = summary.sort_values(by=["search_hit", "sample_count"], ascending=[False, False])

        # Highlight matched rows
        def highlight_hits(row):
            return ["background-color: #fff3cd" if row["search_hit"] == 1 else "" for _ in row]

        st.markdown(f"### Summary for {gene} ({total_samples} samples total)")
        st.dataframe(
            summary.drop(columns=["variant_id"])
                   .style.apply(highlight_hits, axis=1),
            use_container_width=True
        )

        csv = summary.drop(columns=["variant_id"]).to_csv(index=False)
        st.download_button("Download Summary", data=csv, file_name="gene_variant_summary.csv", mime="text/csv")

        # Per-variant breakdown
        st.markdown("---")
        st.markdown("### Per-Variant Sample Breakdown")
        for _, row in summary.iterrows():
            label = f"{row['transcript_variant']} / {row['protein_variant']} (gnomAD: {row['gnomad_frequency'] or 'N/A'})"
            with st.expander(label):
                variant_df = df_gene[df_gene["variant_id"] == row["variant_id"]]
                variant_df = variant_df.drop_duplicates(subset=["normalized_sample_name"])
                variant_df["allele_fraction"] = pd.to_numeric(variant_df["allele_fraction"], errors="coerce").round(3)
                variant_df = variant_df.dropna(subset=["normalized_sample_name"])
                if not variant_df.empty:
                    st.dataframe(variant_df[["sample_name", "allele_fraction", "classification"]], use_container_width=True)
                    csv = variant_df[["sample_name", "allele_fraction", "classification"]].to_csv(index=False)
                    st.download_button("Download Variant Data", data=csv, file_name=f"{row['transcript_variant']}_{row['protein_variant']}_samples.csv", mime="text/csv")
                else:
                    st.info("No sample data available.")

# ===================================================================
# VIEW 3: Multi-Gene Sample Search
# ===================================================================
elif view == "Multi-Gene Sample Search":
    st.subheader("Multi-Gene Sample Search")
    gene_input = st.text_input("Enter one or more gene symbols (comma-separated)", "")
    if gene_input:
        genes = [g.strip().upper() for g in gene_input.split(",") if g.strip()]
        df_filtered = df[df["gene_symbol"].str.upper().isin(genes)].copy()
        df_filtered = df_filtered.drop_duplicates(subset=["normalized_sample_name", "gene_symbol"])
        if df_filtered.empty:
            st.warning("No variants found for the provided gene(s).")
        else:
            sample_gene_map = (
                df_filtered[["normalized_sample_name", "gene_symbol"]]
                .dropna()
                .drop_duplicates()
                .assign(present=1)
                .pivot_table(index="normalized_sample_name", columns="gene_symbol", values="present", fill_value=0)
                .reindex(columns=genes, fill_value=0)
            )
            st.markdown("### Sample Presence Matrix")
            display_matrix = sample_gene_map.replace({1: "✓", 0: ""})
            st.dataframe(display_matrix, use_container_width=True)

            gene_counts = sample_gene_map.sum().astype(int).reset_index()
            gene_counts.columns = ["Gene", "Sample Count"]
            st.markdown("### Sample Counts per Gene")
            st.dataframe(gene_counts, use_container_width=True)

            sample_gene_map["Gene_Count"] = sample_gene_map.sum(axis=1)
            overlap_summary = {
                "Samples with All Genes": (sample_gene_map["Gene_Count"] == len(genes)).sum(),
                "Samples with Any Gene": (sample_gene_map["Gene_Count"] >= 1).sum(),
            }
            if len(genes) == 2:
                only_1 = (sample_gene_map[genes[0]] == 1) & (sample_gene_map[genes[1]] == 0)
                only_2 = (sample_gene_map[genes[0]] == 0) & (sample_gene_map[genes[1]] == 1)
                overlap_summary[f"Samples with {genes[0]} Only"] = only_1.sum()
                overlap_summary[f"Samples with {genes[1]} Only"] = only_2.sum()

            st.markdown("### Sample Overlap Summary")
            st.dataframe(pd.DataFrame(list(overlap_summary.items()), columns=["Condition", "Sample Count"]), use_container_width=True)

            export_df = sample_gene_map.reset_index()
            st.download_button("Download Presence Matrix", data=export_df.to_csv(index=False), file_name="multi_gene_sample_matrix.csv", mime="text/csv")

# ===================================================================
# VIEW 4: Sample → Variant Lookup
# ===================================================================
elif view == "Sample → Variant Lookup":
    st.subheader("Sample → Variant Lookup")
    samples = sorted(df["normalized_sample_name"].dropna().unique())
    selected_sample = st.selectbox("Select Sample", samples)
    if selected_sample:
        sample_variants = df[df["normalized_sample_name"] == selected_sample].drop_duplicates(subset=["variant_id"])
        cols = [
            "gene_symbol",
            "chrom",
            "position",
            "end_position",
            "transcript_variant",
            "protein_variant",
            "classification",
            "allele_fraction",
            "gnomad_frequency"
        ]
        sample_variants = sample_variants[cols].sort_values(by=["gene_symbol", "position"])
        sample_variants["allele_fraction"] = pd.to_numeric(sample_variants["allele_fraction"], errors="coerce").round(3)
        sample_variants["gnomad_frequency"] = pd.to_numeric(sample_variants["gnomad_frequency"], errors="coerce").round(6)
        st.dataframe(sample_variants, use_container_width=True)
        csv = sample_variants.to_csv(index=False)
        st.download_button("Download Sample Variant Details", data=csv, file_name=f"{selected_sample}_variant_details.csv", mime="text/csv")

# ===================================================================
# VIEW 5: Artifact Scoring
# ===================================================================
elif view == "Artifact Scoring":
    st.subheader("Artifact Scoring Summary")
    artifact_df = (
        df.drop_duplicates(subset=["variant_id", "normalized_sample_name"])
        .groupby([
            "variant_id",
            "gene_symbol",
            "transcript_variant",
            "protein_variant",
            "gnomad_frequency",
            "artifact_score",
            "artifact_label",
            "artifact_likelihood_percent"
        ])
        .agg(sample_count=("normalized_sample_name", "nunique"))
        .reset_index()
    )
    artifact_df["variant_presence_percent"] = (artifact_df["sample_count"] / total_samples * 100).round(2)

    st.dataframe(artifact_df.drop(columns=["variant_id"]), use_container_width=True)
    csv = artifact_df.drop(columns=["variant_id"]).to_csv(index=False)
    st.download_button("Download Artifact Scoring CSV", data=csv, file_name="artifact_scoring.csv", mime="text/csv")
