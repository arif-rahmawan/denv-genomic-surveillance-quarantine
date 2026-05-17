"""
validasi_t12_grid_profile.py
Validasi output Tahap 12 — GRID Profile Nukleotida
Arif Rahmawan | Poltekkes Kemenkes Bandung | 2025

Output log: 09_qc/validasi_hasil/validasi_t12_grid_profile_HASIL.txt
"""

import os
import csv
import datetime
from collections import Counter

try:
    from Bio import SeqIO
except ImportError:
    raise SystemExit("ERROR: pip3 install biopython --break-system-packages")

BASE_DIR  = "/Users/arifrahmawan/Downloads/REVISI 5"
GRID_DIR  = os.path.join(BASE_DIR, "04_grid_profile")
ALIGN_DIR = os.path.join(BASE_DIR, "03_alignment")
LOG_DIR   = os.path.join(BASE_DIR, "09_qc", "validasi_hasil")
LOG_FILE  = os.path.join(LOG_DIR, "validasi_t12_grid_profile_HASIL.txt")

SEROTYPES      = ["DENV1", "DENV2", "DENV3", "DENV4"]
GENES          = ["E", "NS1", "prM"]
MAX_FREQ_ERROR = 0.01

_lines = []

def log(msg=""):
    print(msg)
    _lines.append(str(msg))

def save_log():
    os.makedirs(LOG_DIR, exist_ok=True)
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_FILE, "w", encoding="utf-8") as f:
        f.write("Validasi dijalankan: " + ts + "\n\n")
        f.write("\n".join(_lines) + "\n")
    print("\n  Log tersimpan: " + LOG_FILE)

def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def main():
    total_pass = total_fail = total_warn = 0

    log("=" * 65)
    log("VALIDASI TAHAP 12 — GRID PROFILE NUKLEOTIDA")
    log("=" * 65)

    for sero in SEROTYPES:
        for gene in GENES:
            label     = sero + "_" + gene
            grid_file = os.path.join(GRID_DIR, sero, sero + "_" + gene + "_grid_profile.csv")
            aln_file  = os.path.join(ALIGN_DIR, sero, sero + "_" + gene + "_alignment.fasta")

            log("\n  [" + label + "]")

            # 1. File exist
            if not os.path.exists(grid_file):
                log("  GAGAL | File grid_profile.csv tidak ditemukan")
                total_fail += 1
                continue
            log("  LULUS | File grid_profile.csv ada")
            total_pass += 1

            rows  = load_csv(grid_file)
            n_pos = len(rows)

            # 2. Jumlah posisi = panjang alignment
            if os.path.exists(aln_file):
                aln_records = list(SeqIO.parse(aln_file, "fasta"))
                aln_len     = len(aln_records[0].seq) if aln_records else 0
                if n_pos == aln_len:
                    log("  LULUS | Jumlah posisi = panjang alignment (" + str(n_pos) + " bp)")
                    total_pass += 1
                else:
                    log("  GAGAL | Jumlah posisi " + str(n_pos) + " != panjang alignment " + str(aln_len))
                    total_fail += 1

            # 3. Frekuensi sum ≈ 1.0
            bad_sum = []
            for r in rows:
                try:
                    s = (float(r["freq_A"]) + float(r["freq_T"]) +
                         float(r["freq_G"]) + float(r["freq_C"]) + float(r["freq_gap"]))
                    if abs(s - 1.0) > MAX_FREQ_ERROR:
                        bad_sum.append(r["pos_alignment"])
                except (KeyError, ValueError):
                    pass
            if not bad_sum:
                log("  LULUS | Frekuensi A+T+G+C+gap ≈ 1.0 di semua posisi")
                total_pass += 1
            else:
                log("  GAGAL | " + str(len(bad_sum)) + " posisi frekuensi tidak sum ke 1")
                total_fail += 1

            # 4. Distribusi flag
            flags    = Counter(r["flag"] for r in rows)
            n_cons   = flags.get("conserved", 0)
            pct_cons = n_cons / n_pos * 100 if n_pos else 0
            if pct_cons >= 50:
                log("  LULUS | Conserved " + str(round(pct_cons, 1)) + "% (>=50%)")
                total_pass += 1
            else:
                log("  WARN  | Conserved hanya " + str(round(pct_cons, 1)) + "% — periksa alignment")
                total_warn += 1
            log("  INFO  | Flag: " + str(dict(flags)))

            # 5. Dominant nt valid
            bad_dom = [r["pos_alignment"] for r in rows
                       if r.get("flag") != "gap_dominant" and r.get("dominant_nt") == "-"]
            if not bad_dom:
                log("  LULUS | Dominant nt valid di semua posisi non-gap")
                total_pass += 1
            else:
                log("  GAGAL | " + str(len(bad_dom)) + " posisi non-gap dengan dominant='-'")
                total_fail += 1

    log("\n" + "=" * 65)
    log("RINGKASAN: LULUS=" + str(total_pass) + "  WARN=" + str(total_warn) + "  GAGAL=" + str(total_fail))
    if total_fail == 0:
        log("STATUS: TAHAP 12 VALID")
    else:
        log("STATUS: ADA KEGAGALAN — periksa detail di atas")
    log("=" * 65)

if __name__ == "__main__":
    main()
    save_log()
