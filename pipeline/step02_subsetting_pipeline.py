"""
step02_subsetting_pipeline.py
Pipeline subsetting lengkap: ekstraksi gen + QC panjang + QC translasi
Simpan di: /Users/arifrahmawan/Downloads/REVISI 5/10_scripts/pipeline/
Jalankan : python3 step02_subsetting_pipeline.py

Output   : /Users/arifrahmawan/Downloads/REVISI 5/02_subsetting/
"""

from Bio import SeqIO
from Bio.Seq import Seq
import os

# ── KONFIGURASI ──────────────────────────────────────────────────────────────
BASE       = "/Users/arifrahmawan/Downloads/REVISI 5/01_raw_data"
OUTDIR     = "/Users/arifrahmawan/Downloads/REVISI 5/02_subsetting"
REFSEQ     = f"{BASE}/REFSEQ/REFSEQ_ALL.fasta"

# Koordinat gen (1-based, inclusive) dari anotasi mat_peptide RefSeq NCBI
COORDS = {
    "DENV1": {"prM": (437, 934),  "E": (935, 2419),  "NS1": (2420, 3475)},
    "DENV2": {"prM": (439, 936),  "E": (937, 2421),  "NS1": (2422, 3477)},
    "DENV3": {"prM": (437, 934),  "E": (935, 2413),  "NS1": (2414, 3469)},
    "DENV4": {"prM": (441, 938),  "E": (939, 2423),  "NS1": (2424, 3479)},
}

# Accession RefSeq per serotipe
REFSEQ_ACC = {
    "DENV1": "NC_001477.1",
    "DENV2": "NC_001474.2",
    "DENV3": "NC_001475.2",
    "DENV4": "NC_002640.1",
}

# Posisi awal CDS di RefSeq (1-based) — untuk anchor offset detection
REFSEQ_CDS_START = {
    "DENV1": 97,
    "DENV2": 97,
    "DENV3": 97,
    "DENV4": 102,
}

# File isolat clean per serotipe
ISOLAT_FILES = {
    "DENV1": [
        f"{BASE}/INDONESIA/DENV1/IDN_DENV1_clean.fasta",
        f"{BASE}/REGIONAL/Thailand/DENV1/THAI_DENV1_clean.fasta",
        f"{BASE}/REGIONAL/Singapore/DENV1/SGP_DENV1_clean.fasta",
        f"{BASE}/REGIONAL/Philippine/DENV1/PHL_DENV1_clean.fasta",
        f"{BASE}/REGIONAL/Malaysia/DENV1/MYS_DENV1_clean.fasta",
    ],
    "DENV2": [
        f"{BASE}/INDONESIA/DENV2/IDN_DENV2_clean.fasta",
        f"{BASE}/REGIONAL/Thailand/DENV2/THAI_DENV2_clean.fasta",
        f"{BASE}/REGIONAL/Singapore/DENV2/SGP_DENV2_clean.fasta",
        f"{BASE}/REGIONAL/Philippine/DENV2/PHL_DENV2_clean.fasta",
        f"{BASE}/REGIONAL/Malaysia/DENV2/MYS_DENV2_clean.fasta",
    ],
    "DENV3": [
        f"{BASE}/INDONESIA/DENV3/IDN_DENV3_clean.fasta",
        f"{BASE}/REGIONAL/Thailand/DENV3/THAI_DENV3_clean.fasta",
        f"{BASE}/REGIONAL/Singapore/DENV3/SGP_DENV3_clean.fasta",
        f"{BASE}/REGIONAL/Philippine/DENV3/PHL_DENV3_clean.fasta",
    ],
    "DENV4": [
        f"{BASE}/INDONESIA/DENV4/IDN_DENV4_clean.fasta",
        f"{BASE}/REGIONAL/Thailand/DENV4/THAI_DENV4_clean.fasta",
        f"{BASE}/REGIONAL/Singapore/DENV4/SGP_DENV4_clean.fasta",
        f"{BASE}/REGIONAL/Philippine/DENV4/PHL_DENV4_clean.fasta",
        f"{BASE}/REGIONAL/Malaysia/DENV4/MYS_DENV4_clean.fasta",
    ],
}

# Panjang gen yang diharapkan (untuk QC)
EXPECTED_LENGTH = {
    "DENV1": {"prM": 498, "E": 1485, "NS1": 1056},
    "DENV2": {"prM": 498, "E": 1485, "NS1": 1056},
    "DENV3": {"prM": 498, "E": 1479, "NS1": 1056},
    "DENV4": {"prM": 498, "E": 1485, "NS1": 1056},
}

# ── FUNGSI ────────────────────────────────────────────────────────────────────

def load_refseq():
    """Load semua RefSeq dari file REFSEQ_ALL.fasta."""
    refseq = {}
    for r in SeqIO.parse(REFSEQ, "fasta"):
        for sero, acc in REFSEQ_ACC.items():
            if acc in r.id:
                refseq[sero] = str(r.seq).upper()
    return refseq


def detect_offset(isolat_seq, ref_seq, ref_cds_start):
    """
    Deteksi offset antara isolat dan RefSeq.
    Cari anchor 30-nt dari posisi CDS start RefSeq di dalam isolat.
    Return: offset integer, atau None jika tidak ditemukan.
    """
    anchor_start = ref_cds_start - 1
    anchor = ref_seq[anchor_start:anchor_start + 30]

    pos = isolat_seq.find(anchor)
    if pos == -1:
        anchor_short = anchor[:15]
        pos = isolat_seq.find(anchor_short)
        if pos == -1:
            return None

    return pos - anchor_start


def subset_gen(isolat_seq, sero, offset, coords):
    """Ekstrak gen dengan koordinat yang disesuaikan oleh offset."""
    results = {}
    for gene, (start, end) in coords.items():
        adj_start = start + offset - 1
        adj_end   = end + offset

        if adj_start < 0 or adj_end > len(isolat_seq):
            results[gene] = None
            continue

        subseq = isolat_seq[adj_start:adj_end]
        expected = end - start + 1
        results[gene] = subseq if len(subseq) == expected else None

    return results


def has_premature_stop(seq_str):
    """Cek apakah sekuens memiliki stop codon prematur."""
    try:
        protein = Seq(seq_str).translate()
        return str(protein[:-1]).count("*") > 0
    except Exception:
        return True


# ── TAHAP 1: SUBSETTING ───────────────────────────────────────────────────────

def run_subsetting(refseq):
    print("\n" + "─"*60)
    print("  TAHAP 1 — SUBSETTING GEN (Offset Dinamis)")
    print("─"*60)

    gene_records = {
        sero: {"prM": [], "E": [], "NS1": []}
        for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]
    }

    total_ok    = 0
    total_skip  = 0
    total_nooff = 0

    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        ref_seq       = refseq[sero]
        ref_cds_start = REFSEQ_CDS_START[sero]
        coords        = COORDS[sero]

        print(f"\n  {sero}")

        for fasta_path in ISOLAT_FILES[sero]:
            if not os.path.exists(fasta_path):
                continue

            records   = list(SeqIO.parse(fasta_path, "fasta"))
            file_ok   = 0
            file_skip = 0
            file_nooff = 0

            for rec in records:
                isolat_seq = str(rec.seq).upper().replace("-", "")

                offset = detect_offset(isolat_seq, ref_seq, ref_cds_start)
                if offset is None:
                    offset = 0
                    file_nooff += 1

                subseqs = subset_gen(isolat_seq, sero, offset, coords)
                all_ok  = all(v is not None for v in subseqs.values())

                if all_ok:
                    for gene, subseq in subseqs.items():
                        new_rec = SeqIO.SeqRecord(
                            Seq(subseq), id=rec.id, description=""
                        )
                        gene_records[sero][gene].append(new_rec)
                    file_ok += 1
                    total_ok += 1
                else:
                    file_skip += 1
                    total_skip += 1

            total_nooff += file_nooff
            fname = os.path.basename(fasta_path)
            print(f"    {fname}: ✅ OK={file_ok} | offset_gagal={file_nooff} | ❌ skip={file_skip}")

    return gene_records, total_ok, total_skip, total_nooff


# ── TAHAP 2: SIMPAN OUTPUT ────────────────────────────────────────────────────

def save_output(gene_records):
    print("\n" + "─"*60)
    print("  TAHAP 2 — MENYIMPAN OUTPUT")
    print("─"*60)

    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        sero_dir = os.path.join(OUTDIR, sero)
        os.makedirs(sero_dir, exist_ok=True)
        for gene in ["prM", "E", "NS1"]:
            records  = gene_records[sero][gene]
            if not records:
                continue
            out_path = os.path.join(sero_dir, f"{sero}_{gene}_subset.fasta")
            SeqIO.write(records, out_path, "fasta")
            print(f"  ✅ {sero}_{gene}_subset.fasta → {len(records)} sekuens")


# ── TAHAP 3: QC PANJANG ───────────────────────────────────────────────────────

def run_qc_panjang():
    print("\n" + "─"*60)
    print("  TAHAP 3 — QC PANJANG GEN")
    print("─"*60)

    all_ok = True
    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        for gene in ["prM", "E", "NS1"]:
            path     = os.path.join(OUTDIR, sero, f"{sero}_{gene}_subset.fasta")
            expected = EXPECTED_LENGTH[sero][gene]

            if not os.path.exists(path):
                print(f"  ❌ {sero}_{gene}: file tidak ditemukan")
                all_ok = False
                continue

            records = list(SeqIO.parse(path, "fasta"))
            lengths = set(len(r.seq) for r in records)

            if lengths == {expected}:
                print(f"  ✅ {sero}_{gene}: {len(records)} sekuens × {expected} bp — konsisten")
            else:
                print(f"  ❌ {sero}_{gene}: panjang tidak konsisten — {lengths}")
                all_ok = False

    return all_ok


# ── TAHAP 4: QC TRANSLASI ─────────────────────────────────────────────────────

def run_qc_translasi():
    print("\n" + "─"*60)
    print("  TAHAP 4 — QC TRANSLASI IN SILICO (filter stop codon prematur)")
    print("─"*60)

    total_buang = 0
    total_sisa  = 0
    laporan     = []

    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        # Cukup cek gen E sebagai representatif (1 gen = 1 isolat)
        # Filter dilakukan di semua gen sekaligus per isolat
        e_path = os.path.join(OUTDIR, sero, f"{sero}_E_subset.fasta")
        if not os.path.exists(e_path):
            continue

        # Kumpulkan ID yang harus dibuang (berdasarkan gen E)
        buang_ids = set()
        for r in SeqIO.parse(e_path, "fasta"):
            if has_premature_stop(str(r.seq)):
                buang_ids.add(r.id)

        # Terapkan filter ke semua 3 gen
        sero_buang = 0
        sero_sisa  = 0
        for gene in ["prM", "E", "NS1"]:
            path = os.path.join(OUTDIR, sero, f"{sero}_{gene}_subset.fasta")
            if not os.path.exists(path):
                continue
            records = list(SeqIO.parse(path, "fasta"))
            clean   = [r for r in records if r.id not in buang_ids]
            SeqIO.write(clean, path, "fasta")

        sero_buang = len(buang_ids)
        sero_sisa  = len([r for r in SeqIO.parse(e_path, "fasta")])

        total_buang += sero_buang
        total_sisa  += sero_sisa
        laporan.append((sero, sero_buang, sero_sisa))

        if sero_buang > 0:
            print(f"  ⚠️  {sero}: {sero_buang} sekuens dibuang → tersisa {sero_sisa}")
        else:
            print(f"  ✅ {sero}: tidak ada stop codon prematur → {sero_sisa} sekuens")

    return total_buang, total_sisa, laporan


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  STEP02 — SUBSETTING PIPELINE LENGKAP")
    print("  prM | E | NS1  ×  DENV1 | DENV2 | DENV3 | DENV4")
    print("=" * 60)

    # Buat/reset output directory
    os.makedirs(OUTDIR, exist_ok=True)

    # Load RefSeq
    print("\n  Loading RefSeq...")
    refseq = load_refseq()
    if not refseq:
        print("  ❌ RefSeq tidak ditemukan. Cek path REFSEQ.")
        return
    for sero, seq in refseq.items():
        print(f"    {REFSEQ_ACC[sero]}: {len(seq)} bp ✅")

    # Jalankan pipeline
    gene_records, t_ok, t_skip, t_nooff = run_subsetting(refseq)
    save_output(gene_records)
    qc_ok = run_qc_panjang()
    t_buang, t_sisa, laporan = run_qc_translasi()

    # Ringkasan akhir
    print("\n" + "=" * 60)
    print("  RINGKASAN AKHIR")
    print("=" * 60)
    print(f"  Subsetting berhasil  : {t_ok}")
    print(f"  Offset tidak terdeteksi (fixed): {t_nooff}")
    print(f"  Skip (koordinat gagal): {t_skip}")
    print(f"  QC panjang           : {'✅ Semua konsisten' if qc_ok else '❌ Ada masalah'}")
    print(f"  Stop codon dibuang   : {t_buang}")
    print(f"  Dataset final        : {t_sisa}")
    print()
    print("  Per serotipe:")
    for sero, buang, sisa in laporan:
        print(f"    {sero}: dibuang={buang} | final={sisa}")
    print()
    print(f"  Output tersimpan di : {OUTDIR}")
    print("=" * 60)
    print("  ✅ Pipeline selesai. File siap untuk MSA di UGENE.")


if __name__ == "__main__":
    main()
