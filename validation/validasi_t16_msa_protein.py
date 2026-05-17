"""
validasi_t16_msa_protein.py
Validasi output Tahap 16 — MSA Protein per Gen-Serotipe
Arif Rahmawan | Poltekkes Kemenkes Bandung | 2025

Output log: 09_qc/validasi_hasil/validasi_t16_msa_protein_HASIL.txt
"""

import os
import csv
import datetime

try:
    from Bio import SeqIO
except ImportError:
    raise SystemExit("ERROR: pip3 install biopython --break-system-packages")

BASE_DIR    = "/Users/arifrahmawan/Downloads/REVISI 5"
PROTALN_DIR = os.path.join(BASE_DIR, "14_protein_alignment")
LOG_DIR     = os.path.join(BASE_DIR, "09_qc", "validasi_hasil")
LOG_FILE    = os.path.join(LOG_DIR, "validasi_t16_msa_protein_HASIL.txt")

SEROTYPES        = ["DENV1", "DENV2", "DENV3", "DENV4"]
GENES            = ["E", "NS1", "prM"]
MIN_PROTEIN_CONS = 0.85
REFSEQ_PREFIXES  = ("NP_", "YP_")

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
    log("VALIDASI TAHAP 16 — MSA PROTEIN PER GEN-SEROTIPE")
    log("=" * 65)

    # Summary
    summary_file = os.path.join(PROTALN_DIR, "PROTEIN_MSA_SUMMARY.csv")
    if os.path.exists(summary_file):
        log("  LULUS | PROTEIN_MSA_SUMMARY.csv ada")
        total_pass += 1
    else:
        log("  GAGAL | PROTEIN_MSA_SUMMARY.csv tidak ditemukan")
        total_fail += 1

    for sero in SEROTYPES:
        for gene in GENES:
            label     = sero + "_" + gene
            aln_file  = os.path.join(PROTALN_DIR, sero, gene + "_protein_alignment.fasta")
            cons_file = os.path.join(PROTALN_DIR, sero, gene + "_protein_conservation.csv")

            log("\n  [" + label + "]")

            # 1. FASTA alignment ada
            if not os.path.exists(aln_file):
                log("  GAGAL | protein_alignment.fasta tidak ditemukan")
                total_fail += 1
                continue
            log("  LULUS | protein_alignment.fasta ada")
            total_pass += 1

            # 2. CSV conservation ada
            if not os.path.exists(cons_file):
                log("  GAGAL | protein_conservation.csv tidak ditemukan")
                total_fail += 1
            else:
                log("  LULUS | protein_conservation.csv ada")
                total_pass += 1

            # 3. Panjang alignment konsisten
            records = list(SeqIO.parse(aln_file, "fasta"))
            lengths = set(len(r.seq) for r in records)
            if len(lengths) == 1:
                log("  LULUS | Panjang alignment konsisten (" +
                    str(list(lengths)[0]) + " AA, " + str(len(records)) + " records)")
                total_pass += 1
            else:
                log("  GAGAL | Panjang alignment tidak konsisten: " + str(sorted(lengths)))
                total_fail += 1

            # 4. RefSeq ada di alignment
            has_refseq = any(r.id.startswith(REFSEQ_PREFIXES) for r in records)
            if has_refseq:
                ref_id = next(r.id for r in records if r.id.startswith(REFSEQ_PREFIXES))
                log("  LULUS | RefSeq ditemukan: " + ref_id)
                total_pass += 1
            else:
                log("  GAGAL | RefSeq tidak ditemukan di alignment")
                total_fail += 1

            # 5. Konservasi >= 85%
            if os.path.exists(cons_file):
                cons_rows = load_csv(cons_file)
                n_total   = len(cons_rows)
                n_cons    = sum(1 for r in cons_rows if r.get("flag") == "conserved")
                pct_cons  = round(n_cons / n_total * 100, 1) if n_total else 0
                n_var     = sum(1 for r in cons_rows if r.get("flag") == "variable")
                n_hot     = sum(1 for r in cons_rows if r.get("flag") == "hotspot")
                if pct_cons >= MIN_PROTEIN_CONS * 100:
                    log("  LULUS | Konservasi protein " + str(pct_cons) +
                        "% (>=" + str(int(MIN_PROTEIN_CONS * 100)) + "%)")
                    total_pass += 1
                else:
                    log("  GAGAL | Konservasi protein " + str(pct_cons) +
                        "% < " + str(int(MIN_PROTEIN_CONS * 100)) + "%")
                    total_fail += 1
                log("  INFO  | conserved=" + str(n_cons) + " | variable=" + str(n_var) +
                    " | hotspot=" + str(n_hot) + " | total=" + str(n_total) + " AA")

    log("\n" + "=" * 65)
    log("RINGKASAN: LULUS=" + str(total_pass) + "  WARN=" + str(total_warn) + "  GAGAL=" + str(total_fail))
    if total_fail == 0:
        log("STATUS: TAHAP 16 VALID")
    else:
        log("STATUS: ADA KEGAGALAN — periksa detail di atas")
    log("=" * 65)

if __name__ == "__main__":
    main()
    save_log()
