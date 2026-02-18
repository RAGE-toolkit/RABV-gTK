#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROFILE="h10n8_xml_test"
PUBLISH_DIR="${PROJECT_DIR}/test_out/h10n8_xml_test"
DB_NAME="h10n8_286285_test"
EMAIL="${VGTK_EMAIL:-}"
MAX_NORMALIZED_RF="${MAX_NORMALIZED_RF:-0.25}"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --email EMAIL             Email for NCBI calls used by downstream scripts
  --profile NAME            Nextflow profile (default: h10n8_xml_test)
  --publish-dir PATH        Output directory (default: test_out/h10n8_xml_test)
  --db-name NAME            Expected DB base name (default: h10n8_286285_test)
  --max-normalized-rf FLOAT Fail RF check if normalized RF exceeds this value (default: 0.25)
  -h, --help                Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --email)
            EMAIL="$2"
            shift 2
            ;;
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --publish-dir)
            PUBLISH_DIR="$2"
            shift 2
            ;;
        --db-name)
            DB_NAME="$2"
            shift 2
            ;;
        --max-normalized-rf)
            MAX_NORMALIZED_RF="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "${EMAIL}" ]]; then
    echo "--email (or VGTK_EMAIL) is required" >&2
    exit 1
fi

XML_DIR="${PROJECT_DIR}/test_data/h10n8_286285/GenBank-XML"
if [[ ! -d "${XML_DIR}" ]]; then
    echo "Missing fixture XML directory: ${XML_DIR}" >&2
    echo "Run scripts/fetch_h10n8_fixture.sh first." >&2
    exit 1
fi

nextflow run "${PROJECT_DIR}/vgtk-init.nf" \
    -profile "${PROFILE}" \
    --email "${EMAIL}" \
    --publish_dir "${PUBLISH_DIR}" \
    --db_name "${DB_NAME}" \
    -resume

DB_PATH="${PUBLISH_DIR}/${DB_NAME}.db"
if [[ ! -f "${DB_PATH}" ]]; then
    echo "Expected DB not found: ${DB_PATH}" >&2
    exit 1
fi

python "${PROJECT_DIR}/scripts/ValidateDbTree.py" \
    --db "${DB_PATH}" \
    --outdir "${PUBLISH_DIR}/tests/db_tree_validation"

python - "${DB_PATH}" <<'PY'
import sqlite3
import sys

db_path = sys.argv[1]
required_tables = [
    "meta_data",
    "sequences",
    "sequence_alignment",
    "features",
    "insertions",
    "host_taxa",
    "trees",
]

conn = sqlite3.connect(db_path)
try:
    cur = conn.cursor()
    for table in required_tables:
        cur.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name=?", (table,))
        if cur.fetchone()[0] != 1:
            raise SystemExit(f"Missing required table: {table}")
        cur.execute(f"SELECT COUNT(*) FROM {table}")
        n = cur.fetchone()[0]
        if n <= 0:
            raise SystemExit(f"Table {table} is empty")
    print("[info] DB schema/content checks passed")
finally:
    conn.close()
PY

EXPECTED_TREE="${PROJECT_DIR}/test_data/h10n8_286285/expected/iqtree.treefile"
ACTUAL_TREE="$(find "${PUBLISH_DIR}" -type f -name '*.treefile' | head -n 1 || true)"

if [[ -n "${ACTUAL_TREE}" && -f "${EXPECTED_TREE}" ]]; then
    python "${PROJECT_DIR}/scripts/ComputeRobinsonFoulds.py" \
        --tree-a "${EXPECTED_TREE}" \
        --tree-b "${ACTUAL_TREE}" \
        --max-normalized-rf "${MAX_NORMALIZED_RF}"
else
    echo "[warn] RF check skipped (expected baseline tree or actual tree missing)."
    echo "[warn] Expected baseline path: ${EXPECTED_TREE}"
fi

echo "H10N8 integration test completed successfully."