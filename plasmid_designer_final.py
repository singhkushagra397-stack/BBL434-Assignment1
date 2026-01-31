#!/usr/bin/env python3

"""
Universal Plasmid Designer
White-list assembly using GC-skew ORI detection
"""

import argparse

# -------------------------
# FASTA READER
# -------------------------
def load_fasta(path):
    sequence = ""
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip().upper()
    return sequence

# -------------------------
# GC-SKEW ORI FINDER
# -------------------------
def detect_ori_gc_skew(seq):
    skew = 0
    skew_profile = []

    for base in seq:
        if base == "G":
            skew += 1
        elif base == "C":
            skew -= 1
        skew_profile.append(skew)

    return skew_profile.index(min(skew_profile))

# -------------------------
# DESIGN FILE PARSER
# -------------------------
def parse_design(path):
    enzymes = []
    markers = []

    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            label, value = [x.strip() for x in line.split(",")]

            if "site" in label.lower():
                enzymes.append(value)
            else:
                markers.append(value)

    return enzymes, markers

# -------------------------
# MARKER DATABASE PARSER
# -------------------------
def load_marker_table(path):
    db = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or "," not in line:
                continue
            label, value = [x.strip() for x in line.split(",", 1)]
            db[label] = value
    return db

# -------------------------
# DATABASES
# -------------------------
RESTRICTION_DB = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG",
    "SphI": "GCATGC",
    "SalI": "GTCGAC",
    "XbaI": "TCTAGA",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SmaI": "CCCGGG",
    "NotI": "GCGGCCGC"
}

MARKER_SEQS = {
    "Ampicillin": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTG",
    "Kanamycin": "ATGAGCCATATTCAACGGGAAACGTCTTGCTCGAGGCC",
    "Chloramphenicol": "ATGGAGAAAAAAATCACTGGATATACCACCGTTGATATATCC",
    "Blue_White_Selection": "ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGAC"
}

BIOBRICK_SCAR = "TACTAGAG"

# -------------------------
# MAIN
# -------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--markers", required=True)
    parser.add_argument("--output", default="Output.fa")
    args = parser.parse_args()

    genome = load_fasta(args.input)

    ori_pos = detect_ori_gc_skew(genome)
    ori_region = genome[ori_pos - 300 : ori_pos + 300]

    enzymes, requested_markers = parse_design(args.design)
    marker_db = load_marker_table(args.markers)

    plasmid = ori_region

    # Add markers
    for marker in requested_markers:
        if marker in MARKER_SEQS:
            plasmid += BIOBRICK_SCAR + MARKER_SEQS[marker]
        else:
            print(f"WARNING: Marker '{marker}' not found, skipped")

    # Add MCS
    for enzyme in enzymes:
        if enzyme in RESTRICTION_DB:
            plasmid += BIOBRICK_SCAR + RESTRICTION_DB[enzyme]
        else:
            print(f"WARNING: Enzyme '{enzyme}' not found, skipped")

    with open(args.output, "w") as f:
        f.write(">Universal_Plasmid\n")
        for i in range(0, len(plasmid), 70):
            f.write(plasmid[i:i+70] + "\n")

    print("✅ Plasmid constructed successfully")
    print("ORI position:", ori_pos)
    print("Final length:", len(plasmid))

if __name__ == "__main__":
    main()
