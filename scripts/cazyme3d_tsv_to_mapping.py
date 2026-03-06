#!/usr/bin/env python3
"""
Convert CAZyme3D full TSV (cazyme3d_full.tsv) to Foldseek mapping file (cazyme3d_to_cazy.tsv).

The mapping file must have:
  - One column for target ID (must match PDB filename without .pdb, i.e. the ID Foldseek uses).
  - One column for CAZy family (e.g. GH31, CBM50+GH23).

Your cazyme3d_full.tsv has: cazyid, family, uniprot_mapped, ...
  - If your PDB files are named like NP_652145.1.pdb (RefSeq/cazyid), use: --id-column cazyid
  - If your PDB files are named like A0A0G3V1Z8.pdb (UniProt), use: --id-column uniprot_mapped

Usage:
  python cazyme3d_tsv_to_mapping.py cazyme3d_full.tsv -o cazyme3d_to_cazy.tsv
  python cazyme3d_tsv_to_mapping.py cazyme3d_full.tsv -o cazyme3d_to_cazy.tsv --id-column uniprot_mapped
"""
import argparse
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(
        description="Convert CAZyme3D full TSV to target_id -> CAZy ID mapping for run_dbcan structure search."
    )
    ap.add_argument("input_tsv", help="Path to cazyme3d_full.tsv (tab-separated, with cazyid, family, uniprot_mapped, ...)")
    ap.add_argument("-o", "--output", default="cazyme3d_to_cazy.tsv", help="Output mapping file path (default: cazyme3d_to_cazy.tsv)")
    ap.add_argument(
        "--id-column",
        default="cazyid",
        choices=["cazyid", "uniprot_mapped"],
        help="Column to use as target ID (must match PDB filename without .pdb). Use uniprot_mapped if PDBs are named by UniProt (e.g. A0A0G3V1Z8.pdb).",
    )
    ap.add_argument("--family-column", default="family", help="Column name for CAZy family (default: family)")
    ap.add_argument(
        "--filter-missing-id",
        action="store_true",
        help="Drop rows where the chosen ID column is empty/NaN (e.g. uniprot_mapped NULL).",
    )
    args = ap.parse_args()

    try:
        import pandas as pd
    except ImportError:
        print("This script requires pandas: pip install pandas", file=sys.stderr)
        sys.exit(1)

    path = Path(args.input_tsv)
    if not path.exists():
        print(f"Error: input file not found: {path}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(path, sep="\t")
    if args.id_column not in df.columns:
        print(f"Error: column '{args.id_column}' not found. Available: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)
    if args.family_column not in df.columns:
        print(f"Error: column '{args.family_column}' not found. Available: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    out = df[[args.id_column, args.family_column]].copy()
    out.columns = ["target", "CAZy ID"]
    if args.filter_missing_id:
        out = out[out["target"].notna() & (out["target"].astype(str).str.strip() != "")]
    out = out.drop_duplicates(subset=["target"], keep="first")
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(out)} rows to {args.output}")


if __name__ == "__main__":
    main()
