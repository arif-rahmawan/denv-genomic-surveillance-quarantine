"""
step02b_refseq_subset.py
Subset gen dari RefSeq (prM, E, NS1) per serotipe → file terpisah _refseq.fasta
Simpan di: /Users/arifrahmawan/Downloads/REVISI 5/10_scripts/pipeline/
Jalankan : python3 step02b_refseq_subset.py
"""

from Bio import SeqIO
from Bio.Seq import Seq
import os

# ── KONFIGURASI ──────────────────────────────────────────────────────────────
BASE   = "/Users/arifrahmawan/Downloads/REVISI 5"
REFSEQ = f"{BASE}/01_raw_data/REFSEQ/REFSEQ_ALL.fasta"
OUTDIR = f"{BASE}/02_subsetting"

REFSEQ_ACC = {
    "DENV1": "NC_001477.1",
    "DENV2": "NC_001474.2",
    "DENV3": "NC_001475.2",
    "DENV4": "NC_002640.1",
}

# Koordinat gen (1-based, inclusive) — sama dengan step02
COORDS = {
    "DENV1": {"prM": (437, 934),  "E": (935, 2419),  "NS1": (2420, 3475)},
    "DENV2": {"prM": (439, 936),  "E": (937, 2421),  "NS1": (2422, 3477)},
    "DENV3": {"prM": (437, 934),  "E": (935, 2413),  "NS1": (2414, 3469)},
    "DENV4": {"prM": (441, 938),  "E": (939, 2423),  "NS1": (2424, 3479)},
}

# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  STEP02B — REFSEQ SUBSET (prM | E | NS1)")
    print("=" * 60)

    # Load semua RefSeq
    refseq_records = {}
    for rec in SeqIO.parse(REFSEQ, "fasta"):
        for sero, acc in REFSEQ_ACC.items():
            if acc in rec.id:
                refseq_records[sero] = rec
    
    print(f"\n  RefSeq dimuat: {len(refseq_records)} serotipe")
    for sero, rec in refseq_records.items():
        print(f"    {REFSEQ_ACC[sero]}: {len(rec.seq)} bp ✅")

    print()

    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        if sero not in refseq_records:
            print(f"  ❌ {sero}: RefSeq tidak ditemukan, skip")
            continue

        rec    = refseq_records[sero]
        seq    = str(rec.seq).upper()
        coords = COORDS[sero]
        outdir = os.path.join(OUTDIR, sero)
        os.makedirs(outdir, exist_ok=True)

        for gene, (start, end) in coords.items():
            subseq   = seq[start - 1:end]  # konversi ke 0-based
            out_path = os.path.join(outdir, f"{sero}_{gene}_refseq.fasta")

            new_rec = SeqIO.SeqRecord(
                Seq(subseq),
                id=REFSEQ_ACC[sero],
                description=f"{sero} {gene} RefSeq | {start}-{end} | {len(subseq)} bp"
            )
            SeqIO.write(new_rec, out_path, "fasta")
            print(f"  ✅ {sero}_{gene}_refseq.fasta → {len(subseq)} bp")

    print()
    print("=" * 60)
    print("  ✅ Selesai. RefSeq subset tersimpan di 02_subsetting/")
    print("=" * 60)

if __name__ == "__main__":
    main()
