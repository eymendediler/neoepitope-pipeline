import pandas as pd
import streamlit as st
from pathlib import Path

st.set_page_config(page_title="Neoepitope Candidate Browser", layout="wide")

st.title("Neoepitope Candidate Browser")
st.caption("MDT-supported and isoform-specific neoepitope candidates")

BASE = Path(__file__).resolve().parents[1]
MDT_FILE = BASE / "exports" / "top_mdt_supported_candidates.tsv"
NON_MDT_FILE = BASE / "exports" / "non_mdt_or_unannotated_candidates.tsv"

@st.cache_data
def load_data():
    mdt = pd.read_csv(MDT_FILE, sep="\t")
    non_mdt = pd.read_csv(NON_MDT_FILE, sep="\t")
    mdt["candidate_group"] = "MDT_supported"
    non_mdt["candidate_group"] = "non_MDT_or_isoform_specific"
    return mdt, non_mdt

mdt, non_mdt = load_data()

tab1, tab2 = st.tabs(["MDT-supported candidates", "Non-MDT / isoform-specific candidates"])

def filter_df(df, prefix):
    col1, col2, col3 = st.columns(3)

    with col1:
        gene = st.text_input("Gene search", "", key=f"{prefix}_gene")
    with col2:
        allele = st.text_input("HLA allele search", "", key=f"{prefix}_allele")
    with col3:
        peptide = st.text_input("Peptide search", "", key=f"{prefix}_peptide")

    out = df.copy()

    if gene:
        out = out[out["gene_name"].fillna("").str.contains(gene, case=False, na=False)]
    if allele:
        out = out[out["allele"].fillna("").str.contains(allele, case=False, na=False)]
    if peptide:
        out = out[out["peptide"].fillna("").str.contains(peptide, case=False, na=False)]

    return out

with tab1:
    st.subheader("MDT-supported neoepitope candidates")
    st.write(f"Total rows: {len(mdt)}")

    filtered = filter_df(mdt, "mdt")

    if "best_rank" in filtered.columns:
        filtered = filtered.sort_values("best_rank", ascending=True)

    st.dataframe(filtered, use_container_width=True)
    st.download_button(
        "Download filtered MDT-supported candidates",
        filtered.to_csv(sep="\t", index=False),
        file_name="filtered_mdt_supported_candidates.tsv"
    )

with tab2:
    st.subheader("Non-MDT / isoform-specific candidate pool")
    st.write(f"Total rows: {len(non_mdt)}")

    filtered = filter_df(non_mdt, "non_mdt")

    if "consensus_rank" in filtered.columns:
        filtered = filtered.sort_values("consensus_rank", ascending=True)

    st.dataframe(filtered, use_container_width=True)
    st.download_button(
        "Download filtered non-MDT candidates",
        filtered.to_csv(sep="\t", index=False),
        file_name="filtered_non_mdt_candidates.tsv"
    )
