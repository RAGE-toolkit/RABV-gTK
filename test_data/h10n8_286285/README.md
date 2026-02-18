# H10N8 test fixture (taxid 286285)

This folder stores a frozen GenBank XML snapshot for integration testing of the influenza segmented workflow.

## Taxonomy

- Taxonomy ID: `286285`
- Name: `H10N8 subtype`

## Contents

- `GenBank-XML/` — downloaded XML batch files used as deterministic input for `--xml_dir`
- `manifest.tsv` — metadata for when/what was downloaded
- `checksums.sha256` — checksums for XML snapshot files
- `expected/iqtree.treefile` — optional baseline tree for RF regression checks

## Refreshing the XML snapshot

Run from repository root:

```bash
scripts/fetch_h10n8_fixture.sh --email your_email@example.com --clean
```

If the upstream taxonomy content changes, refresh this snapshot and review expected outputs.

## Running integration checks

```bash
scripts/run_h10n8_xml_test.sh --email your_email@example.com
```

If `expected/iqtree.treefile` exists, the runner computes Robinson-Foulds distance and enforces threshold checks.