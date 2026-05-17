# RUNBOOK — TAHAP 12 & 13
# DENV Indonesia In Silico | Arif Rahmawan

## URUTAN EKSEKUSI

### Step A — Pastikan environment siap
pip install biopython pandas

### Step B — Cek bahwa alignment sudah ada (output Tahap 7)
ls "/Users/arifrahmawan/Downloads/REVISI 5/03_alignment/DENV1/"
# Harus ada: E_alignment.fasta, NS1_alignment.fasta, prM_alignment.fasta
# (dan sama untuk DENV2, DENV3, DENV4)

### Step C — Jalankan Tahap 12 (GRID Profile)
cd "/Users/arifrahmawan/Downloads/REVISI 5/10_scripts/pipeline"
python3 step12_grid_profile.py

# Output di: 04_grid_profile/
# ├── DENV1/E_grid_profile.csv
# ├── DENV1/NS1_grid_profile.csv
# ...dst (12 file CSV)
# └── summary_grid_ALL.csv

### Step D — Jalankan Tahap 13 (Identifikasi Mutasi)
python3 step13_identify_mutations.py

# Output di: 05_mutations/
# ├── DENV1/E_mutations.csv              ← semua posisi
# ├── DENV1/E_mutations_nonconserved.csv ← posisi variabel saja
# ...dst (24 file CSV)
# ├── HOTSPOT_SUMMARY.csv
# └── MUTATION_SUMMARY.csv

---

## KOLOM PENTING OUTPUT TAHAP 12

| Kolom           | Penjelasan                                |
|-----------------|-------------------------------------------|
| pos_alignment   | Posisi di alignment (1-based)             |
| A/T/G/C/gap     | Jumlah absolut per nukleotida             |
| freq_A/T/G/C    | Frekuensi relatif (0.0–1.0)              |
| shannon_entropy | Entropy Shannon (0=monomorphic, max≈2)    |
| dominant_nt     | Nukleotida paling sering                  |
| dominant_freq   | Frekuensi dominant nt                     |
| flag            | conserved / variable / hotspot / gap_dom  |

---

## KLASIFIKASI MUTASI (TAHAP 13)

| classification  | Kriteria                                              |
|-----------------|-------------------------------------------------------|
| conserved       | dominant ≥95%, dominant = WT                         |
| fixed_mutation  | dominant ≥95%, dominant ≠ WT (semua isolat bermutasi) |
| variable_WT     | dominant 80–95%, dominant = WT (5–20% alt)           |
| variable_mut    | dominant 80–95%, dominant ≠ WT                       |
| hotspot_WT      | dominant <80%, dominant = WT (>20% alt)              |
| hotspot_mut     | dominant <80%, dominant ≠ WT                         |
| insertion       | RefSeq gap → isolat ada nukleotida                   |
| deletion        | RefSeq ada nukleotida → isolat dominan gap           |
| gap_dominant    | >80% isolat gap (artefak/region tidak ter-cover)     |

---

## SETELAH TAHAP 11 (FISHER) SELESAI

Jalankan ulang Tahap 13:
python3 step13_identify_mutations.py

Kolom p_fisher_raw, p_bonferroni, sig_bonferroni akan terisi otomatis
dari file 06_fisher/DENV{X}/{gene}_fisher.csv

---

## TROUBLESHOOTING

Q: "RefSeq tidak ditemukan di alignment"
A: Cek header RefSeq di alignment file:
   python3 -c "
   from Bio import SeqIO
   for r in SeqIO.parse('03_alignment/DENV1/E_alignment.fasta', 'fasta'):
       print(r.id[:60])
       break  # hapus ini untuk lihat semua header
   "
   Lalu sesuaikan REFSEQ_PATTERN di step13 jika perlu.

Q: "Tidak ada isolat IDN_ ditemukan"
A: Cek apakah header sudah di-relabel dengan format IDN_{Kota}_{Tahun}_{AccID}
   Jika belum, jalankan ulang step03_merge_relabel.py (Tahap 5).
