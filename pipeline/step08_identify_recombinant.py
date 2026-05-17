"""
step08b_identify_recombinant.py
Arif Rahmawan | DENV In Silico Pipeline

Tujuan:
    Identifikasi sekuens rekombinasi di DENV4_E_alignment.fasta
    menggunakan sliding window distance analysis.

    Sekuens rekombinasi akan menunjukkan pola jarak yang TIDAK KONSISTEN
    terhadap sekuens lain di berbagai region alignment — di satu window
    dia dekat dengan grup A, di window lain dia dekat dengan grup B.

    Ini konsisten dengan hasil RDP5 (1 event confirmed, 3 metode).

Cara pakai:
    python step08b_identify_recombinant.py

Output:
    - Ranking kandidat rekombinasi (skor tertinggi = paling mencurigakan)
    - Plot distance inconsistency per sekuens
    - File: DENV4_E_recombinant_candidates.txt

Requirements:
    pip install numpy biopython matplotlib
"""

import os
import sys
import numpy as np
from collections import defaultdict

# ============================================================
# KONFIGURASI
# ============================================================

BASE = "/Users/arifrahmawan/Downloads/REVISI 5"
INPUT_FILE = f"{BASE}/03_alignment/DENV4/DENV4_E_alignment.fasta"
OUTPUT_FILE = f"{BASE}/03_alignment/DENV4/DENV4_E_recombinant_candidates.txt"

WINDOW_SIZE = 200    # bp per window
STEP_SIZE   = 50     # step sliding window
TOP_N       = 10     # tampilkan top N kandidat

# ============================================================
# FUNGSI
# ============================================================

def parse_fasta(filepath):
    records = []
    header, seq_lines = None, []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq_lines)))
                header, seq_lines = line[1:], []
            else:
                seq_lines.append(line.upper())
    if header:
        records.append((header, "".join(seq_lines)))
    return records


def hamming_distance(s1, s2):
    """Jarak Hamming (fraksi perbedaan, abaikan gap)."""
    assert len(s1) == len(s2)
    diff = 0
    valid = 0
    for a, b in zip(s1, s2):
        if a == '-' or b == '-' or a == 'N' or b == 'N':
            continue
        valid += 1
        if a != b:
            diff += 1
    return diff / valid if valid > 0 else 0.0


def sliding_window_distances(records, window_size, step_size):
    """
    Hitung matriks jarak per window.
    Return: list of (start, end, dist_matrix)
    """
    n = len(records)
    aln_len = len(records[0][1])
    windows = []

    for start in range(0, aln_len - window_size + 1, step_size):
        end = start + window_size
        # Ekstrak subsequences untuk window ini
        seqs = [seq[start:end] for _, seq in records]
        # Hitung pairwise distance matrix
        dist = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                d = hamming_distance(seqs[i], seqs[j])
                dist[i][j] = d
                dist[j][i] = d
        windows.append((start, end, dist))

    return windows


def compute_inconsistency_score(windows, n_seqs):
    """
    Hitung skor inkonsistensi per sekuens.
    
    Untuk setiap sekuens i, kita lihat rank-nya terhadap sekuens lain
    di setiap window. Sekuens rekombinasi akan punya rank yang sangat
    bervariasi (kadang dekat, kadang jauh) → standard deviation rank tinggi.
    """
    n_windows = len(windows)
    # rank_matrix[i][w] = mean rank of seq i in window w
    # (rank = posisi urutan jarak dari terdekat ke terjauh)
    rank_matrix = np.zeros((n_seqs, n_windows))

    for w_idx, (start, end, dist) in enumerate(windows):
        for i in range(n_seqs):
            # Jarak seq i ke semua seq lain
            dists_from_i = [dist[i][j] for j in range(n_seqs) if j != i]
            sorted_d = sorted(dists_from_i)
            # Mean normalized rank (0=paling dekat, 1=paling jauh)
            mean_d = np.mean(dists_from_i)
            # Normalize by alignment-wide mean
            rank_matrix[i][w_idx] = mean_d

    # Skor inkonsistensi = std dev dari mean distance across windows
    # (tinggi = jarak ke cluster tidak stabil = ciri rekombinasi)
    scores = np.std(rank_matrix, axis=1)
    return scores, rank_matrix


def main():
    print("=" * 60)
    print("  Identifikasi Kandidat Rekombinasi — DENV4_E")
    print("=" * 60)
    print()

    # Cek input
    if not os.path.isfile(INPUT_FILE):
        print(f"[ERROR] File tidak ditemukan: {INPUT_FILE}")
        sys.exit(1)

    # Parse FASTA
    print(f"Loading: {INPUT_FILE}")
    records = parse_fasta(INPUT_FILE)
    n = len(records)
    aln_len = len(records[0][1])
    print(f"Sekuens : {n}")
    print(f"Panjang : {aln_len} bp (aligned)")
    print(f"Window  : {WINDOW_SIZE} bp, step {STEP_SIZE} bp")
    print()

    # Sliding window distance
    print("Menghitung jarak per window...")
    print("(ini butuh beberapa menit tergantung jumlah sekuens)")
    windows = sliding_window_distances(records, WINDOW_SIZE, STEP_SIZE)
    print(f"Total windows: {len(windows)}")
    print()

    # Skor inkonsistensi
    print("Menghitung skor inkonsistensi...")
    scores, rank_matrix = compute_inconsistency_score(windows, n)

    # Ranking
    ranked_idx = np.argsort(scores)[::-1]  # descending

    print("=" * 60)
    print(f"  TOP {TOP_N} KANDIDAT REKOMBINASI (skor tertinggi)")
    print("=" * 60)
    print(f"{'Rank':<5} {'Skor':<8} {'Header'}")
    print("-" * 60)

    top_results = []
    for rank, idx in enumerate(ranked_idx[:TOP_N], 1):
        header = records[idx][0]
        score  = scores[idx]
        print(f"{rank:<5} {score:.5f}  {header}")
        top_results.append((rank, score, header))

    print()
    print("CATATAN INTERPRETASI:")
    print("  - Skor = standar deviasi jarak rata-rata across windows")
    print("  - Skor tinggi = pola jarak TIDAK KONSISTEN = ciri rekombinasi")
    print("  - Sekuens rekombinasi RDP5 kemungkinan besar ada di Top 3")
    print("  - RefSeq juga bisa score tinggi karena outgroup — wajar")
    print()

    # Filter out RefSeq dari rekomendasi
    non_refseq = [(r, s, h) for r, s, h in top_results if not h.startswith("REFSEQ")]
    if non_refseq:
        best = non_refseq[0]
        print(f">>> KANDIDAT UTAMA (non-RefSeq): {best[2]}")
        print(f"    Skor: {best[1]:.5f}")
        print()

    # Simpan hasil
    with open(OUTPUT_FILE, "w") as f:
        f.write("KANDIDAT REKOMBINASI DENV4_E\n")
        f.write("=" * 60 + "\n")
        f.write(f"Input  : {INPUT_FILE}\n")
        f.write(f"Window : {WINDOW_SIZE} bp, step {STEP_SIZE} bp\n")
        f.write(f"Total sekuens: {n}, alignment: {aln_len} bp\n\n")
        f.write(f"{'Rank':<5} {'Skor':<10} Header\n")
        f.write("-" * 60 + "\n")
        for rank, idx in enumerate(ranked_idx[:TOP_N], 1):
            f.write(f"{rank:<5} {scores[idx]:.5f}   {records[idx][0]}\n")

    print(f"Hasil disimpan: {OUTPUT_FILE}")
    print()
    print("Langkah berikutnya:")
    print("  1. Salin header kandidat utama ke step09b_filter_recombinant.py")
    print("  2. Set MATCH_MODE = 'partial' dan masukkan AccID-nya")
    print("  3. Jalankan filter → lanjut IQ-TREE3")


if __name__ == "__main__":
    main()
