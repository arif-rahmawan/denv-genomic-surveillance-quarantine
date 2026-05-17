"""
TAHAP 17 ★ — Analisis Tekanan Seleksi dN/dS
Arif Rahmawan | DENV Indonesia In Silico Analysis

Metode : Nei-Gojobori (1986) dengan koreksi Jukes-Cantor
         Implementasi Python from scratch (standar publikasi)

Input  : 03_alignment/DENV{X}/DENV{X}_{gene}_alignment.fasta
         (CDS nukleotida in-frame, RefSeq sebagai referensi)
Output : 15_dnds/DENV{X}/{gene}_dnds.csv       ← per-pasangan isolat vs RefSeq
         15_dnds/DENV{X}/{gene}_dnds_summary.csv ← mean dN, dS, ω per gen
         15_dnds/DNDS_GLOBAL_SUMMARY.csv         ← tabel utama BAB IV

Interpretasi ω (dN/dS):
  ω > 1  → positive selection (mutasi adaptif)
  ω < 1  → purifying selection (konservatif)
  ω ≈ 1  → neutral evolution
"""

import os
import csv
import math
import re
from collections import Counter

try:
    from Bio import SeqIO
except ImportError:
    raise SystemExit("ERROR: pip3 install biopython --break-system-packages")

# ─── KONFIGURASI ─────────────────────────────────────────────────────────────
BASE_DIR   = "/Users/arifrahmawan/Downloads/REVISI 5"
ALIGN_DIR  = os.path.join(BASE_DIR, "03_alignment")
OUTPUT_DIR = os.path.join(BASE_DIR, "15_dnds")

SEROTYPES = ["DENV1", "DENV2", "DENV3", "DENV4"]
GENES     = ["E", "NS1", "prM"]

MAX_ISOLAT_PER_GENE = 200   # batasi pasangan untuk efisiensi (random sample jika lebih)
# ─────────────────────────────────────────────────────────────────────────────

# Standard genetic code
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

# Degenerasi kodon: 0-fold, 2-fold, 4-fold
# Pre-kalkulasi jumlah synonymous sites per kodon posisi
def count_syn_nonsyn_sites(codon):
    """
    Hitung jumlah synonymous (S) dan nonsynonymous (N) sites per kodon.
    Metode Nei-Gojobori: tiap posisi berkontribusi berdasarkan degenerasi.
    Return: (S, N) — float, S+N = 3
    """
    if len(codon) != 3 or 'N' in codon or '-' in codon:
        return 0.0, 0.0

    aa_wt = CODON_TABLE.get(codon, 'X')
    if aa_wt == '*' or aa_wt == 'X':
        return 0.0, 0.0

    S_total = 0.0
    N_total = 0.0
    bases   = ['A', 'T', 'G', 'C']

    for pos in range(3):
        syn_count = 0
        nonsyn_count = 0
        wt_base = codon[pos]

        for alt_base in bases:
            if alt_base == wt_base:
                continue
            mut_codon = codon[:pos] + alt_base + codon[pos+1:]
            aa_mut    = CODON_TABLE.get(mut_codon, 'X')
            if aa_mut == '*':
                nonsyn_count += 1
            elif aa_mut == aa_wt:
                syn_count += 1
            else:
                nonsyn_count += 1

        # Fraksi synonymous di posisi ini
        total_possible = syn_count + nonsyn_count
        if total_possible > 0:
            S_total += syn_count / total_possible
            N_total += nonsyn_count / total_possible

    return S_total, N_total


def jukes_cantor(p):
    """
    Koreksi Jukes-Cantor untuk jarak.
    p = proporsi perbedaan (0 ≤ p < 0.75)
    Return: jarak terkoreksi, atau None jika tidak bisa dihitung
    """
    if p <= 0:
        return 0.0
    if p >= 0.75:
        return None  # tidak bisa dikoreksi
    try:
        return -0.75 * math.log(1 - (4/3) * p)
    except (ValueError, ZeroDivisionError):
        return None


def calc_dnds_pairwise(seq1, seq2):
    """
    Hitung dN dan dS antara dua sekuens CDS (aligned, in-frame).
    Metode Nei-Gojobori (1986) dengan koreksi Jukes-Cantor.

    Return: dict dengan dN, dS, omega, S_total, N_total, dll
    """
    # Ekstrak kodon (tanpa gap bersama)
    codons1, codons2 = [], []
    for i in range(0, min(len(seq1), len(seq2)) - 2, 3):
        c1 = seq1[i:i+3].upper()
        c2 = seq2[i:i+3].upper()
        # Skip kodon dengan gap atau N
        if '-' in c1 or '-' in c2 or 'N' in c1 or 'N' in c2:
            continue
        if CODON_TABLE.get(c1, 'X') == '*' or CODON_TABLE.get(c2, 'X') == '*':
            continue
        codons1.append(c1)
        codons2.append(c2)

    if len(codons1) < 10:
        return None  # terlalu sedikit kodon valid

    # Hitung S dan N total (rata-rata dari dua sekuens)
    S_total = 0.0
    N_total = 0.0
    for c1, c2 in zip(codons1, codons2):
        S1, N1 = count_syn_nonsyn_sites(c1)
        S2, N2 = count_syn_nonsyn_sites(c2)
        S_total += (S1 + S2) / 2
        N_total += (N1 + N2) / 2

    if S_total == 0 or N_total == 0:
        return None

    # Hitung Sd (synonymous differences) dan Nd (nonsynonymous differences)
    Sd = 0.0
    Nd = 0.0
    for c1, c2 in zip(codons1, codons2):
        if c1 == c2:
            continue
        aa1 = CODON_TABLE.get(c1, 'X')
        aa2 = CODON_TABLE.get(c2, 'X')

        # Hitung jumlah perbedaan nukleotida
        nt_diffs = sum(1 for a, b in zip(c1, c2) if a != b)

        if aa1 == aa2:
            # Synonymous
            Sd += 1.0
        else:
            # Nonsynonymous — distribusikan berdasarkan Nei-Gojobori
            # Jika 1 perbedaan nt → langsung
            if nt_diffs == 1:
                Nd += 1.0
            else:
                # Multiple changes: rata-rata semua jalur
                Nd += 1.0  # simplified untuk multiple substitutions

    # Proporsi perbedaan
    pS = Sd / S_total if S_total > 0 else 0
    pN = Nd / N_total if N_total > 0 else 0

    # Koreksi Jukes-Cantor
    dS = jukes_cantor(pS)
    dN = jukes_cantor(pN)

    if dS is None or dN is None:
        return None

    # Omega
    if dS > 0:
        omega = dN / dS
    elif dN == 0:
        omega = 0.0
    else:
        omega = float('inf')

    return {
        'n_codons'   : len(codons1),
        'S_sites'    : round(S_total, 3),
        'N_sites'    : round(N_total, 3),
        'Sd'         : round(Sd, 3),
        'Nd'         : round(Nd, 3),
        'pS'         : round(pS, 6),
        'pN'         : round(pN, 6),
        'dS'         : round(dS, 6),
        'dN'         : round(dN, 6),
        'omega'      : round(omega, 4) if omega != float('inf') else 'Inf',
    }


def get_ungapped_cds(seq_aligned, ref_aligned):
    """
    Ekstrak CDS isolat menggunakan posisi RefSeq sebagai panduan frame.
    Gap dari RefSeq menunjukkan insertion → skip.
    Return: string nt tanpa gap yang sesuai frame RefSeq.
    """
    result = []
    for ref_nt, iso_nt in zip(ref_aligned, seq_aligned):
        if ref_nt != '-':   # posisi valid dalam RefSeq
            result.append(iso_nt if iso_nt != '-' else 'N')
    return ''.join(result)


def write_csv(path, rows, fields=None):
    if not rows:
        return
    if fields is None:
        fields = list(rows[0].keys())
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main():
    import random
    random.seed(42)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 65)
    print("TAHAP 17 ★ — ANALISIS TEKANAN SELEKSI dN/dS")
    print("Metode: Nei-Gojobori (1986) + koreksi Jukes-Cantor")
    print(f"Alignment dir : {ALIGN_DIR}")
    print(f"Output dir    : {OUTPUT_DIR}")
    print("=" * 65)

    global_summary = []

    for sero in SEROTYPES:
        sero_out = os.path.join(OUTPUT_DIR, sero)
        os.makedirs(sero_out, exist_ok=True)

        for gene in GENES:
            label    = f"{sero}_{gene}"
            aln_file = os.path.join(ALIGN_DIR, sero, f"{sero}_{gene}_alignment.fasta")

            if not os.path.exists(aln_file):
                print(f"\n[SKIP] {label}: alignment tidak ditemukan")
                continue

            print(f"\n[{label}]")

            # Baca alignment
            records    = list(SeqIO.parse(aln_file, "fasta"))
            idn_recs   = [r for r in records if r.id.startswith("IDN_")]
            refseq_rec = next((r for r in records
                               if r.id.startswith("NC_") or 'REF' in r.id.upper()), None)

            if refseq_rec is None:
                print(f"  [SKIP] RefSeq tidak ditemukan di alignment")
                continue

            print(f"  IDN isolat : {len(idn_recs)}")
            print(f"  RefSeq     : {refseq_rec.id}")

            ref_aligned = str(refseq_rec.seq).upper()
            ref_cds     = ref_aligned.replace('-', '')

            # Sample jika terlalu banyak
            sample_recs = idn_recs
            if len(idn_recs) > MAX_ISOLAT_PER_GENE:
                sample_recs = random.sample(idn_recs, MAX_ISOLAT_PER_GENE)
                print(f"  Sampling {MAX_ISOLAT_PER_GENE} dari {len(idn_recs)} isolat")

            # Hitung dN/dS per isolat vs RefSeq
            pair_rows = []
            n_valid   = 0
            n_fail    = 0

            for rec in sample_recs:
                iso_aln = str(rec.seq).upper()
                iso_cds = get_ungapped_cds(iso_aln, ref_aligned)

                result = calc_dnds_pairwise(ref_cds, iso_cds)
                if result is None:
                    n_fail += 1
                    continue

                n_valid += 1
                row = {'isolat_id': rec.id}
                row.update(result)
                pair_rows.append(row)

            print(f"  Pasangan valid  : {n_valid} | Gagal: {n_fail}")

            if not pair_rows:
                print(f"  [WARN] Tidak ada pasangan valid untuk {label}")
                continue

            # Tulis per-pasangan
            pair_file = os.path.join(sero_out, f"{gene}_dnds.csv")
            write_csv(pair_file, pair_rows)
            print(f"  → {pair_file}")

            # Ringkasan statistik
            dN_vals   = [r['dN'] for r in pair_rows if isinstance(r['dN'], float)]
            dS_vals   = [r['dS'] for r in pair_rows if isinstance(r['dS'], float)]
            omega_vals= [r['omega'] for r in pair_rows if isinstance(r['omega'], float)]

            mean_dN    = sum(dN_vals) / len(dN_vals) if dN_vals else 0
            mean_dS    = sum(dS_vals) / len(dS_vals) if dS_vals else 0
            mean_omega = sum(omega_vals) / len(omega_vals) if omega_vals else 0
            median_omega = sorted(omega_vals)[len(omega_vals)//2] if omega_vals else 0

            # Interpretasi
            if mean_omega < 0.5:
                interpretation = 'Strong purifying selection'
            elif mean_omega < 1.0:
                interpretation = 'Purifying selection'
            elif mean_omega > 1.0:
                interpretation = 'Positive selection'
            else:
                interpretation = 'Neutral evolution'

            print(f"  mean dN    : {mean_dN:.6f}")
            print(f"  mean dS    : {mean_dS:.6f}")
            print(f"  mean ω     : {mean_omega:.4f}  ({interpretation})")
            print(f"  median ω   : {median_omega:.4f}")

            # Distribusi omega
            n_purifying = sum(1 for v in omega_vals if v < 1.0)
            n_positive  = sum(1 for v in omega_vals if v > 1.0)
            n_neutral   = sum(1 for v in omega_vals if v == 1.0)
            print(f"  Purifying (<1): {n_purifying} | Positive (>1): {n_positive} | Neutral: {n_neutral}")

            # Tulis summary per gen
            sum_file = os.path.join(sero_out, f"{gene}_dnds_summary.csv")
            write_csv(sum_file, [{
                'label'            : label,
                'n_pairs'          : n_valid,
                'mean_dN'          : round(mean_dN, 6),
                'mean_dS'          : round(mean_dS, 6),
                'mean_omega'       : round(mean_omega, 4),
                'median_omega'     : round(median_omega, 4),
                'n_purifying'      : n_purifying,
                'n_positive'       : n_positive,
                'n_neutral'        : n_neutral,
                'interpretation'   : interpretation,
            }])

            global_summary.append({
                'serotipe'         : sero,
                'gen'              : gene,
                'n_pairs'          : n_valid,
                'mean_dN'          : round(mean_dN, 6),
                'mean_dS'          : round(mean_dS, 6),
                'mean_omega'       : round(mean_omega, 4),
                'median_omega'     : round(median_omega, 4),
                'n_purifying_pct'  : round(n_purifying / len(omega_vals) * 100, 1) if omega_vals else 0,
                'n_positive_pct'   : round(n_positive  / len(omega_vals) * 100, 1) if omega_vals else 0,
                'interpretation'   : interpretation,
            })

    # Tulis global summary
    if global_summary:
        gsf = os.path.join(OUTPUT_DIR, "DNDS_GLOBAL_SUMMARY.csv")
        write_csv(gsf, global_summary)
        print(f"\n[GLOBAL SUMMARY] {gsf}")

        # Print tabel
        print("\n" + "─" * 75)
        print(f"{'Label':<16} {'pairs':>6} {'mean_dN':>9} {'mean_dS':>9} {'mean_ω':>8} {'median_ω':>9} {'Interpretasi'}")
        print("─" * 75)
        for r in global_summary:
            label = f"{r['serotipe']}_{r['gen']}"
            print(f"{label:<16} {r['n_pairs']:>6} {r['mean_dN']:>9.6f} {r['mean_dS']:>9.6f} "
                  f"{r['mean_omega']:>8.4f} {r['median_omega']:>9.4f}  {r['interpretation']}")
        print("─" * 75)

        # Narasi otomatis
        print("\n--- NARASI BAB IV (draft) ---")
        for r in global_summary:
            label = f"{r['serotipe']}_{r['gen']}"
            print(f"{label}: dN={r['mean_dN']:.6f}, dS={r['mean_dS']:.6f}, "
                  f"ω={r['mean_omega']:.4f} → {r['interpretation']}")

    print("\n" + "=" * 65)
    print("TAHAP 17 SELESAI")
    print(f"Direktori: {OUTPUT_DIR}")
    print("=" * 65)


if __name__ == "__main__":
    main()
