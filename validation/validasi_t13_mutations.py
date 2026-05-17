"""
validasi_t13_mutations.py
Validasi output Tahap 13 — Identifikasi & Klasifikasi Mutasi Nukleotida
Arif Rahmawan | Poltekkes Kemenkes Bandung | 2025

Output log: 09_qc/validasi_hasil/validasi_t13_mutations_HASIL.txt
"""

import os
import csv
import datetime
from collections import Counter

BASE_DIR = "/Users/arifrahmawan/Downloads/REVISI 5"
MUT_DIR  = os.path.join(BASE_DIR, "05_mutations")
LOG_DIR  = os.path.join(BASE_DIR, "09_qc", "validasi_hasil")
LOG_FILE = os.path.join(LOG_DIR, "validasi_t13_mutations_HASIL.txt")

SEROTYPES = ["DENV1", "DENV2", "DENV3", "DENV4"]
GENES     = ["E", "NS1", "prM"]

VALID_CLASSES = {
    "conserved", "fixed_mutation", "variable_WT", "variable_mut",
    "hotspot_WT", "hotspot_mut", "insertion", "deletion", "gap_dominant"
}

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
    log("VALIDASI TAHAP 13 — IDENTIFIKASI & KLASIFIKASI MUTASI")
    log("=" * 65)

    # File global
    for fname in ["HOTSPOT_SUMMARY.csv", "MUTATION_SUMMARY.csv"]:
        fpath = os.path.join(MUT_DIR, fname)
        if os.path.exists(fpath):
            log("  LULUS | " + fname + " ada")
            total_pass += 1
        else:
            log("  GAGAL | " + fname + " tidak ditemukan")
            total_fail += 1

    for sero in SEROTYPES:
        for gene in GENES:
            label    = sero + "_" + gene
            mut_file = os.path.join(MUT_DIR, sero, gene + "_mutations.csv")
            nc_file  = os.path.join(MUT_DIR, sero, gene + "_mutations_nonconserved.csv")

            log("\n  [" + label + "]")

            # 1. File exist
            for f_path, f_name in [(mut_file, "mutations.csv"), (nc_file, "mutations_nonconserved.csv")]:
                if os.path.exists(f_path):
                    log("  LULUS | " + f_name + " ada")
                    total_pass += 1
                else:
                    log("  GAGAL | " + f_name + " tidak ditemukan")
                    total_fail += 1

            if not os.path.exists(mut_file):
                continue

            rows = load_csv(mut_file)

            # 2. Klasifikasi valid
            bad_cls = [r["pos_alignment"] for r in rows
                       if r.get("classification") not in VALID_CLASSES]
            if not bad_cls:
                log("  LULUS | Semua nilai classification valid")
                total_pass += 1
            else:
                log("  GAGAL | " + str(len(bad_cls)) + " posisi dengan classification tidak dikenal")
                total_fail += 1

            # 3. Conserved → dominant = WT
            wrong_cons = [r["pos_alignment"] for r in rows
                          if r.get("classification") == "conserved"
                          and r.get("nt_dominant") != r.get("nt_WT")
                          and r.get("nt_WT") not in ("-", "N")]
            if not wrong_cons:
                log("  LULUS | Posisi conserved: dominant_nt = nt_WT (konsisten)")
                total_pass += 1
            else:
                log("  GAGAL | " + str(len(wrong_cons)) + " posisi conserved tapi dominant != WT")
                total_fail += 1

            # 4. freq_dominant + freq_alternative <= 1
            bad_freq = [r["pos_alignment"] for r in rows
                        if float(r.get("freq_dominant", 0)) +
                           float(r.get("freq_alternative", 0)) > 1.01]
            if not bad_freq:
                log("  LULUS | freq_dominant + freq_alternative <= 1.0")
                total_pass += 1
            else:
                log("  GAGAL | " + str(len(bad_freq)) + " posisi frekuensi melebihi 1.0")
                total_fail += 1

            # 5. Fisher tidak pending
            pending = sum(1 for r in rows if r.get("p_fisher_raw") == "pending")
            if pending == 0:
                log("  LULUS | Kolom Fisher sudah terisi semua")
                total_pass += 1
            else:
                log("  WARN  | " + str(pending) + " posisi Fisher masih 'pending'")
                total_warn += 1

            # Statistik
            cls_counts = Counter(r["classification"] for r in rows)
            log("  INFO  | conserved=" + str(cls_counts.get("conserved", 0)) +
                " | fixed_mut=" + str(cls_counts.get("fixed_mutation", 0)) +
                " | hotspot=" + str(cls_counts.get("hotspot_WT", 0) + cls_counts.get("hotspot_mut", 0)) +
                " | variable=" + str(cls_counts.get("variable_WT", 0) + cls_counts.get("variable_mut", 0)))

    log("\n" + "=" * 65)
    log("RINGKASAN: LULUS=" + str(total_pass) + "  WARN=" + str(total_warn) + "  GAGAL=" + str(total_fail))
    if total_fail == 0:
        log("STATUS: TAHAP 13 VALID")
    else:
        log("STATUS: ADA KEGAGALAN — periksa detail di atas")
    log("=" * 65)

if __name__ == "__main__":
    main()
    save_log()
