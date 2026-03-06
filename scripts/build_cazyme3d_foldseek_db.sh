#!/usr/bin/env bash
# Build Foldseek database from a directory of CAZyme3D PDB files.
# Usage:
#   ./scripts/build_cazyme3d_foldseek_db.sh /path/to/pdb_dir [output_prefix] [index_tmp_dir]
# Example:
#   ./scripts/build_cazyme3d_foldseek_db.sh /data/cazyme3d_pdbs CAZyme3D /tmp/foldseek_idx

set -e
PDB_DIR="${1:?Usage: $0 <pdb_dir> [output_prefix] [index_tmp_dir]}"
DB_PREFIX="${2:-CAZyme3D}"
INDEX_TMP="${3:-/tmp/foldseek_index}"

if ! command -v foldseek &>/dev/null; then
  echo "Error: foldseek not found in PATH. Install from https://github.com/steineggerlab/foldseek" >&2
  exit 1
fi
if [[ ! -d "$PDB_DIR" ]]; then
  echo "Error: PDB directory not found: $PDB_DIR" >&2
  exit 1
fi

echo "Creating Foldseek DB from $PDB_DIR with prefix $DB_PREFIX ..."
foldseek createdb "$PDB_DIR" "$DB_PREFIX"
echo "Creating index (tmp: $INDEX_TMP) ..."
mkdir -p "$INDEX_TMP"
foldseek createindex "$DB_PREFIX" "$INDEX_TMP" --split-memory-limit 40
echo "Done. Copy all ${DB_PREFIX}* files and cazyme3d_to_cazy.tsv to your run_dbcan --db_dir."
