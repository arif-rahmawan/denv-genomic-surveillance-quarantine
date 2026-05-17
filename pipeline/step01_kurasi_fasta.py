"""
step01_kurasi_fasta.py
Kurasi FASTA hasil download GISAID untuk REVISI 5
Output: laporan ringkas + FASTA bersih per file
"""
from Bio import SeqIO
import os

BASE = "/Users/arifrahmawan/Downloads/REVISI 5/01_raw_data"

FILES = {
    "IDN_DENV1": f"{BASE}/INDONESIA/DENV1/IDN_DENV1_raw.fasta",
    "IDN_DENV2": f"{BASE}/INDONESIA/DENV2/IDN_DENV2_raw.fasta",
    "IDN_DENV3": f"{BASE}/INDONESIA/DENV3/IDN_DENV3_raw.fasta",
    "IDN_DENV4": f"{BASE}/INDONESIA/DENV4/IDN_DENV4_raw.fasta",
    "THA_DENV1": f"{BASE}/REGIONAL/Thailand/DENV1/THAI_DENV1_raw.fasta",
    "THA_DENV2": f"{BASE}/REGIONAL/Thailand/DENV2/THAI_DENV2_raw.fasta",
    "THA_DENV3": f"{BASE}/REGIONAL/Thailand/DENV3/THAI_DENV3_raw.fasta",
    "THA_DENV4": f"{BASE}/REGIONAL/Thailand/DENV4/THAI_DENV4_raw.fasta",
    "SGP_DENV1": f"{BASE}/REGIONAL/Singapore/DENV1/SGP_DENV1_raw.fasta",
    "SGP_DENV2": f"{BASE}/REGIONAL/Singapore/DENV2/SGP_DENV2_raw.fasta",
    "SGP_DENV3": f"{BASE}/REGIONAL/Singapore/DENV3/SGP_DENV3_raw.fasta",
    "SGP_DENV4": f"{BASE}/REGIONAL/Singapore/DENV4/SGP_DENV4_raw.fasta",
    "PHL_DENV1": f"{BASE}/REGIONAL/Philippine/DENV1/PHL_DENV1_raw.fasta",
    "PHL_DENV2": f"{BASE}/REGIONAL/Philippine/DENV2/PHL_DENV2_raw.fasta",
    "PHL_DENV3": f"{BASE}/REGIONAL/Philippine/DENV3/PHL_DENV3_raw.fasta",
    "PHL_DENV4": f"{BASE}/REGIONAL/Philippine/DENV4/PHL_DENV4_raw.fasta",
    "MYS_DENV1": f"{BASE}/REGIONAL/Malaysia/DENV1/MYS_DENV1_raw.fasta",
    "MYS_DENV2": f"{BASE}/REGIONAL/Malaysia/DENV2/MYS_DENV2_raw.fasta",
    "MYS_DENV4": f"{BASE}/REGIONAL/Malaysia/DENV4/MYS_DENV4_raw.fasta",
}

MIN_LEN    = 500
MAX_N_PCT  = 0.05
total_raw  = 0
total_ok   = 0

print("=" * 60)
print("  STEP01 — KURASI FASTA REVISI 5")
print("=" * 60)

for label, path in FILES.items():
    if not os.path.exists(path):
        print(f"\n  ❌ {label}: FILE TIDAK DITEMUKAN")
        print(f"     {path}")
        continue

    records   = list(SeqIO.parse(path, "fasta"))
    seen_ids  = set()
    clean     = []

    buang_pendek = 0
    buang_n      = 0
    buang_duplik = 0

    for r in records:
        seq = str(r.seq).upper()

        # Duplikat
        if r.id in seen_ids:
            buang_duplik += 1
            continue
        seen_ids.add(r.id)

        # Terlalu pendek
        if len(seq) < MIN_LEN:
            buang_pendek += 1
            continue

        # N terlalu banyak
        n_pct = seq.count("N") / len(seq)
        if n_pct > MAX_N_PCT:
            buang_n += 1
            continue

        clean.append(r)

    # Simpan FASTA bersih
    out_path = path.replace("_raw.fasta", "_clean.fasta")
    SeqIO.write(clean, out_path, "fasta")

    total_raw += len(records)
    total_ok  += len(clean)
    buang_total = len(records) - len(clean)

    print(f"\n  ✅ {label}")
    print(f"     Raw: {len(records)} | Bersih: {len(clean)} | Dibuang: {buang_total}")
    if buang_duplik: print(f"     Duplikat   : {buang_duplik}")
    if buang_pendek: print(f"     Terlalu pendek: {buang_pendek}")
    if buang_n:      print(f"     N >5%      : {buang_n}")
    print(f"     Output     : {os.path.basename(out_path)}")

print()
print("=" * 60)
print(f"  TOTAL RAW   : {total_raw}")
print(f"  TOTAL BERSIH: {total_ok}")
print(f"  TOTAL BUANG : {total_raw - total_ok}")
print("=" * 60)
print("  ✅ Kurasi selesai. File _clean.fasta tersimpan.")
print("     Lanjut ke step02 (merge + RefSeq)")
