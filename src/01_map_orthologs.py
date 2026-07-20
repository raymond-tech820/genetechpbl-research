"""
Step 1. Build a mouse -> human ortholog table.

Geneformer is human-pretrained and keyed on HUMAN Ensembl gene IDs (ENSG...).
Our data (TMS, GSE176206) is mouse gene symbols. This script pulls a
mouse-symbol -> human-Ensembl mapping from Ensembl BioMart and caches it to CSV,
so downstream steps can just read the CSV (BioMart can be slow/flaky, so we only
hit it once).

Run once:  python src/01_map_orthologs.py
Output:    data/mouse_to_human_orthologs.csv
"""

from pathlib import Path
import pandas as pd

OUT = Path("data")
OUT.mkdir(exist_ok=True)
OUT_CSV = OUT / "mouse_to_human_orthologs.csv"


def fetch_orthologs() -> pd.DataFrame:
    from pybiomart import Server

    server = Server(host="http://www.ensembl.org")
    mouse = (
        server.marts["ENSEMBL_MART_ENSEMBL"]
        .datasets["mmusculus_gene_ensembl"]
    )

    df = mouse.query(
        attributes=[
            "external_gene_name",              # mouse gene symbol
            "ensembl_gene_id",                 # mouse Ensembl id
            "hsapiens_homolog_ensembl_gene",   # human ortholog Ensembl id (ENSG...)
            "hsapiens_homolog_orthology_type", # one2one / one2many / many2many
        ]
    )
    df.columns = ["mouse_symbol", "mouse_ensembl", "human_ensembl", "orthology_type"]
    return df


def main():
    df = fetch_orthologs()

    # Keep only clean 1:1 orthologs, safest for scoring, avoids duplicating genes.
    # If you lose too many genes you can relax this later, but start strict.
    before = df["mouse_symbol"].nunique()
    df = df[df["orthology_type"] == "ortholog_one2one"].dropna(subset=["human_ensembl"])
    df = df.drop_duplicates(subset=["mouse_symbol"])
    after = df["mouse_symbol"].nunique()

    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV}  ({after} one2one orthologs kept of {before} mouse genes)")
    print(df.head())


if __name__ == "__main__":
    main()
