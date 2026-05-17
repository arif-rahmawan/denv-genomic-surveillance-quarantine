"""
validasi_t17_dnds.py
Validasi output Tahap 17 — Analisis Tekanan Seleksi dN/dS
Arif Rahmawan | Poltekkes Kemenkes Bandung | 2025

Output log: 09_qc/validasi_hasil/validasi_t17_dnds_HASIL.txt
"""

import os
import csv
import datetime

BASE_DIR = "/Users/arifrahmawan/Downloads/REVISI 5"
DNDS_DIR = os.path.join(BASE_DIR, "15_dnds")
LOG_DIR  = os.path.join(BASE_DIR, "09_qc", "validasi_hasil")
LOG_FILE = os.path.join(LOG_DIR, "validasi_t17_dnds_HASIL.txt")

SEROTYPES = ["DENV1", "DENV2", "DENV3", "DENV4"]
GENES     = ["E", "NS1", "prM"]

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
    log("VALIDASI TAHAP 17 — ANALISIS TEKANAN SELEKSI dN/dS")
    log("Metode: Nei-Gojobori (1986) + koreksi Jukes-Cantor")
    log("=" * 65)

    # 1. Global summary
    global_file = os.path.join(DNDS_DIR, "DNDS_GLOBAL_SUMMARY.csv")
    if not os.path.exists(global_file):
        log("  GAGAL | DNDS_GLOBAL_SUMMARY.csv tidak ditemukan")
        log("\n" + "=" * 65)
        log("STATUS: GAGAL — jalankan step17_dnds.py terlebih dahulu")
        log("=" * 65)
        total_fail += 1
        return
    log("  LULUS | DNDS_GLOBAL_SUMMARY.csv ada")
    total_pass += 1

    global_rows = load_csv(global_file)

    # Cek semua mean omega < 1.0
    bad_omega_global = [r for r in global_rows
                        if r.get("mean_omega", "") not in ("", "Inf")
                        and float(r.get("mean_omega", 1)) >= 1.0]
    if not bad_omega_global:
        log("  LULUS | Semua mean omega < 1.0 (purifying selection terkonfirmasi)")
        total_pass += 1
    else:
        log("  GAGAL | " + str(len(bad_omega_global)) + " gen-serotipe dengan mean omega >= 1.0")
        total_fail += 1

    # Tabel global
    log("\n  Global Summary:")
    log("  " + "-" * 60)
    log("  " + "Label".ljust(16) + "mean_dN".rjust(10) + "mean_dS".rjust(10) +
        "mean_w".rjust(8) + "  Interpretasi")
    log("  " + "-" * 60)
    for r in global_rows:
        lbl = r.get("serotipe", "") + "_" + r.get("gen", "")
        log("  " + lbl.ljust(16) +
            str(r.get("mean_dN", "?")).rjust(10) +
            str(r.get("mean_dS", "?")).rjust(10) +
            str(r.get("mean_omega", "?")).rjust(8) +
            "  " + r.get("interpretation", "?"))

    # Per gen-serotipe
    for sero in SEROTYPES:
        for gene in GENES:
            label     = sero + "_" + gene
            dnds_file = os.path.join(DNDS_DIR, sero, gene + "_dnds.csv")

            log("\n  [" + label + "]")

            # 2. File ada
            if not os.path.exists(dnds_file):
                log("  GAGAL | dnds.csv tidak ditemukan")
                total_fail += 1
                continue
            log("  LULUS | dnds.csv ada")
            total_pass += 1

            rows = load_csv(dnds_file)
            if not rows:
                log("  GAGAL | File kosong")
                total_fail += 1
                continue

            # 3. Semua omega valid
            bad_w  = []
            w_vals = []
            for r in rows:
                v = r.get("omega", "")
                if v == "Inf":
                    w_vals.append(float("inf"))
                else:
                    try:
                        w_vals.append(float(v))
                    except ValueError:
                        bad_w.append(r.get("isolat_id", "?"))

            if not bad_w:
                log("  LULUS | Semua nilai omega valid (" + str(len(rows)) + " isolat)")
                total_pass += 1
            else:
                log("  GAGAL | " + str(len(bad_w)) + " isolat dengan omega tidak valid")
                total_fail += 1

            # 4. dN < dS untuk semua isolat
            not_purifying = []
            for r in rows:
                try:
                    dn = float(r.get("dN", 0))
                    ds = float(r.get("dS", 0))
                    if ds > 0 and dn >= ds:
                        not_purifying.append(r.get("isolat_id", "?"))
                except ValueError:
                    pass
            if not not_purifying:
                log("  LULUS | dN < dS untuk semua isolat")
                total_pass += 1
            else:
                log("  GAGAL | " + str(len(not_purifying)) + " isolat dN >= dS")
                total_fail += 1

            # 5. Range omega
            valid_w = [v for v in w_vals if v != float("inf")]
            if valid_w:
                mean_w = round(sum(valid_w) / len(valid_w), 4)
                min_w  = round(min(valid_w), 4)
                max_w  = round(max(valid_w), 4)
                if max_w < 1.0:
                    log("  LULUS | Semua omega < 1.0 (strong purifying)")
                    total_pass += 1
                else:
                    log("  WARN  | Ada omega >= 1.0 (max=" + str(max_w) + ") — periksa manual")
                    total_warn += 1
                log("  INFO  | mean_w=" + str(mean_w) + " | min=" + str(min_w) +
                    " | max=" + str(max_w) + " | n=" + str(len(rows)))

            # 6. Spot-check 3 isolat pertama
            log("  Spot-check 3 isolat pertama:")
            for r in rows[:3]:
                log("    " + str(r.get("isolat_id", "?"))[:35] +
                    " | dN=" + str(r.get("dN", "?")) +
                    " | dS=" + str(r.get("dS", "?")) +
                    " | omega=" + str(r.get("omega", "?")))

    log("\n" + "=" * 65)
    log("RINGKASAN: LULUS=" + str(total_pass) + "  WARN=" + str(total_warn) + "  GAGAL=" + str(total_fail))
    if total_fail == 0:
        log("STATUS: TAHAP 17 VALID")
    else:
        log("STATUS: ADA KEGAGALAN — periksa detail di atas")
    log("=" * 65)

if __name__ == "__main__":
    main()
    save_log()
