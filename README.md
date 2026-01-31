# Universal Plasmid Designer

This tool constructs a functional plasmid for an unknown organism by:
- Identifying the origin of replication (ORI) using GC-skew analysis
- Assembling a new plasmid using modules specified in a user-provided design file

## Methodology

1. **ORI Identification**
   - The input genome is scanned using GC-skew analysis.
   - The minimum GC-skew position is assumed to correspond to the ORI.
   - A 600 bp region surrounding this position is extracted.

2. **Design-Based Assembly**
   - The plasmid is assembled using a white-list approach.
   - Only genetic elements listed in the design file are included.
   - Restriction sites or markers not listed (e.g., EcoRI in pUC19) are excluded by construction.

3. **Robustness**
   - Missing markers or enzymes are skipped with a warning.
   - The program does not crash on invalid entries.

## Test Case

The repository includes a test case using the plasmid **pUC19**:
- `pUC19.fa` – input genome
- `Design_pUC19.txt` – design specification
- `markers.tab` – marker database

EcoRI is present in the original pUC19 sequence but absent in the output plasmid, confirming correct design interpretation.

## Usage

```bash
python plasmid_designer_final.py \
  --input tests/pUC19.fa \
  --design tests/Design_pUC19.txt \
  --markers tests/markers.tab \
  --output Output.fa
