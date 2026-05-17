"""
validasi_t15_protein.py
Validasi output Tahap 15 — Translasi Protein + BLASTp Verifikasi
Arif Rahmawan | Poltekkes Kemenkes Bandung | 2025

Output log: 09_qc/validasi_hasil/validasi_t15_protein_HASIL.txt
"""

import os
import csv
import datetime

try:
    from Bio import SeqIO
except ImportError:
    raise SystemExit("ERROR: pip3 install biopython --break-system-packages")

BASE_DIR = "/Users/arifrahmawan/Downloads/REVISI 5"
PROT_DIR = os.path.join(BASE_DIR, "13_protein")
LOG_DIR  = os.path.join(BASE_DIR, "09_qc", "validasi_hasil")
LOG_FILE = os.path.join(LOG_DIR, "validasi_t15_protein_HASIL.txt")

SEROTYPES            = ["DENV1", "DENV2", "DENV3", "DENV4"]
GENES                = ["E", "NS1", "prM"]
MAX_BLASTP_FAIL_RATE = 0.15
IDN_TOTAL            = {"DENV1": 43, "DENV2": 193, "DENV3": 45, "DENV4": 105}

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
    log("VALIDASI TAHAP 15 — TRANSLASI PROTEIN + BLASTP VERIFIKASI")
    log("=" * 65)

    # PROTEIN_SUMMARY.csv
    summary_file = os.path.join(PROT_DIR, "PROTEIN_SUMMARY.csv")
    if os.path.exists(summary_file):
        log("  LULUS | PROTEIN_SUMMARY.csv ada")
        total_pass += 1
    else:
        log("  GAGAL | PROTEIN_SUMMARY.csv tidak ditemukan")
        total_fail += 1

    for sero in SEROTYPES:
        for gene in GENES:
            label      = sero + "_" + gene
            fasta_file = os.path.join(PROT_DIR, sero, gene + "_protein_IDN.fasta")
            aa_file    = os.path.join(PROT_DIR, sero, gene + "_mutations_aa.csv")

            log("\n  [" + label + "]")

            # 1. FASTA ada
            if not os.path.exists(fasta_file):
                log("  GAGAL | FASTA protein tidak ditemukan")
                total_fail += 1
                continue
            log("  LULUS | FASTA protein IDN ada")
            total_pass += 1

            prot_records = list(SeqIO.parse(fasta_file, "fasta"))
            n_pass  = len(prot_records)
            n_total = IDN_TOTAL.get(sero, 0)

            # 2. BLASTp fail rate
            fail_rate = (n_total - n_pass) / n_total if n_total > 0 else 0
            pct_fail  = round(fail_rate * 100, 1)
            if fail_rate <= MAX_BLASTP_FAIL_RATE:
                log("  LULUS | BLASTp fail rate " + str(pct_fail) + "% (<=" +
                    str(int(MAX_BLASTP_FAIL_RATE * 100)) + "%) — " +
                    str(n_pass) + "/" + str(n_total) + " lolos")
                total_pass += 1
            else:
                log("  WARN  | BLASTp fail rate tinggi: " + str(pct_fail) + "% (" +
                    str(n_total - n_pass) + " isolat gagal)")
                total_warn += 1

            # 3. Tidak ada stop codon
            stop_list = [r.id for r in prot_records if "*" in str(r.seq)]
            if not stop_list:
                log("  LULUS | Tidak ada stop codon di protein terverifikasi")
                total_pass += 1
            else:
                log("  GAGAL | " + str(len(stop_list)) + " isolat mengandung stop codon")
                total_fail += 1

            # 4. Panjang protein konsisten
            lengths = set(len(r.seq) for r in prot_records)
            if len(lengths) <= 2:
                log("  LULUS | Panjang protein konsisten: " + str(sorted(lengths)) + " AA")
                total_pass += 1
            else:
                log("  WARN  | Banyak varian panjang protein: " + str(sorted(lengths)))
                total_warn += 1

            # 5. Rasio NS/S < 1.0
            if os.path.exists(aa_file):
                aa_rows  = load_csv(aa_file)
                n_nonsyn = sum(1 for r in aa_rows if r.get("mutation_type") == "nonsynonymous")
                n_syn    = sum(1 for r in aa_rows if r.get("mutation_type") == "synonymous")
                ratio    = round(n_nonsyn / n_syn, 4) if n_syn > 0 else 0
                if ratio < 1.0:
                    log("  LULUS | Rasio NS/S = " + str(ratio) + " < 1.0 (purifying selection)")
                    total_pass += 1
                else:
                    log("  GAGAL | Rasio NS/S = " + str(ratio) + " >= 1.0 — tidak wajar")
                    total_fail += 1
                log("  INFO  | Nonsynonymous=" + str(n_nonsyn) + " | Synonymous=" + str(n_syn))

    log("\n" + "=" * 65)
    log("RINGKASAN: LULUS=" + str(total_pass) + "  WARN=" + str(total_warn) + "  GAGAL=" + str(total_fail))
    if total_fail == 0:
        log("STATUS: TAHAP 15 VALID")
    else:
        log("STATUS: ADA KEGAGALAN — periksa detail di atas")
    log("=" * 65)

if __name__ == "__main__":
    main()
    save_log()
