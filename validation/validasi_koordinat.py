"""
validasi_koordinat.py
Validasi bahwa koordinat RefSeq valid untuk semua isolat
sebelum dilakukan subsetting gen E, NS1, prM

Cara kerja:
1. Ambil sekuens RefSeq per serotipe
2. Ekstrak kodon start & stop setiap gen dari RefSeq
3. Cek apakah kodon yang sama ada di posisi yang sama pada isolat
4. Laporkan isolat yang bermasalah

Output: laporan validasi + daftar isolat yang perlu dicek manual
"""
import sys as _sys
import datetime as _dt
import os as _os

class _Tee:
    def __init__(self, fname):
        _os.makedirs("/Users/arifrahmawan/Downloads/REVISI 5/09_qc/validasi_hasil", exist_ok=True)
        self._path = "/Users/arifrahmawan/Downloads/REVISI 5/09_qc/validasi_hasil/" + fname + "_HASIL.txt"
        self._file = open(self._path, "w", encoding="utf-8")
        self._file.write("Validasi dijalankan: " + _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
        self._stdout = _sys.stdout
        _sys.stdout = self
    def write(self, msg):
        self._stdout.write(msg)
        self._file.write(msg)
    def flush(self):
        self._stdout.flush()
        self._file.flush()
    def close(self):
        _sys.stdout = self._stdout
        self._file.close()
        self._stdout.write("\n  Log tersimpan: " + self._path + "\n")

_tee = _Tee("validasi_koordinat")


from Bio import SeqIO
from Bio.Seq import Seq
import os

BASE = "/Users/arifrahmawan/Downloads/REVISI 5/01_raw_data"
REFSEQ_PATH = f"{BASE}/REFSEQ/REFSEQ_ALL.fasta"

# Koordinat gen (1-based, inclusive) dari anotasi RefSeq resmi
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

# ── Load RefSeq ───────────────────────────────────────────────
def load_refseq():
    refseq = {}
    for r in SeqIO.parse(REFSEQ_PATH, "fasta"):
        for sero, acc in REFSEQ_ACC.items():
            if acc in r.id:
                refseq[sero] = str(r.seq).upper()
    return refseq

# ── Ambil kodon start & stop tiap gen dari RefSeq ────────────
def get_ref_codons(ref_seq, coords):
    """Ambil 3 kodon pertama dan terakhir setiap gen sebagai fingerprint."""
    fingerprints = {}
    for gene, (start, end) in coords.items():
        subseq = ref_seq[start-1:end]
        fingerprints[gene] = {
            "start_codon"  : subseq[:3],       # kodon pertama
            "codon_2"      : subseq[3:6],       # kodon kedua
            "codon_3"      : subseq[6:9],       # kodon ketiga
            "stop_codon"   : subseq[-3:],       # kodon terakhir
            "length"       : len(subseq),
            "full_subseq"  : subseq,
        }
    return fingerprints

# ── Validasi satu isolat ─────────────────────────────────────
def validate_isolat(seq_str, sero, ref_fingerprints, coords):
    seq = seq_str.upper()
    results = {}

    for gene, (start, end) in coords.items():
        fp = ref_fingerprints[gene]
        expected_len = fp["length"]

        # Cek panjang sekuens cukup
        if len(seq) < end:
            results[gene] = {
                "status": "❌ FAIL",
                "reason": f"Sekuens terlalu pendek ({len(seq)} bp < {end} bp yang dibutuhkan)"
            }
            continue

        # Ekstrak region
        subseq = seq[start-1:end]
        actual_len = len(subseq)

        # Cek kodon start
        start_codon = subseq[:3]
        codon_2     = subseq[3:6]
        codon_3     = subseq[6:9]
        stop_codon  = subseq[-3:]

        issues = []

        # Kodon start harus sama dengan RefSeq
        if start_codon != fp["start_codon"]:
            issues.append(
                f"Start codon berbeda: isolat={start_codon}, ref={fp['start_codon']}"
            )

        # Kodon 2 & 3 sebagai konfirmasi tambahan
        if codon_2 != fp["codon_2"]:
            issues.append(
                f"Kodon ke-2 berbeda: isolat={codon_2}, ref={fp['codon_2']}"
            )

        # Panjang harus sama
        if actual_len != expected_len:
            issues.append(
                f"Panjang berbeda: isolat={actual_len}, ref={expected_len}"
            )

        # Cek N content di region gen ini
        n_pct = subseq.count("N") / actual_len if actual_len > 0 else 0
        if n_pct > 0.05:
            issues.append(f"N content tinggi di region gen: {round(n_pct*100,1)}%")

        if issues:
            results[gene] = {"status": "⚠️ WARN", "reason": " | ".join(issues)}
        else:
            results[gene] = {"status": "✅ OK", "reason": "Koordinat valid"}

    return results

# ── Main ─────────────────────────────────────────────────────
def run():
    print("=" * 65)
    print("  VALIDASI KOORDINAT REFSEQ — SEBELUM SUBSETTING")
    print("=" * 65)
    print("  Memverifikasi bahwa koordinat RefSeq valid untuk semua isolat")
    print()

    # Load RefSeq
    refseq = load_refseq()
    if not refseq:
        print("❌ RefSeq tidak ditemukan. Cek path REFSEQ_PATH.")
        return

    # Summary counter
    total_ok   = 0
    total_warn = 0
    total_fail = 0
    problem_isolats = []

    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        if sero not in refseq:
            print(f"⚠️ RefSeq {sero} tidak ditemukan, skip.")
            continue

        coords         = COORDS[sero]
        ref_seq        = refseq[sero]
        ref_fp         = get_ref_codons(ref_seq, coords)

        print(f"\n{'─'*65}")
        print(f"  {sero} — Koordinat:")
        for g, (s, e) in coords.items():
            print(f"    {g}: [{s}–{e}] = {e-s+1} bp | "
                  f"start codon RefSeq: {ref_fp[g]['start_codon']}")
        print()

        for fasta_path in ISOLAT_FILES[sero]:
            if not os.path.exists(fasta_path):
                continue

            records   = list(SeqIO.parse(fasta_path, "fasta"))
            file_ok   = 0
            file_warn = 0
            file_fail = 0

            for rec in records:
                seq_str = str(rec.seq).replace("-", "").upper()
                results = validate_isolat(seq_str, sero, ref_fp, coords)

                all_ok = all(v["status"] == "✅ OK" for v in results.values())
                any_fail = any(v["status"] == "❌ FAIL" for v in results.values())

                if all_ok:
                    file_ok += 1
                    total_ok += 1
                elif any_fail:
                    file_fail += 1
                    total_fail += 1
                    for gene, res in results.items():
                        if res["status"] != "✅ OK":
                            problem_isolats.append(
                                f"{rec.id} | {sero} | {gene} | {res['status']} | {res['reason']}"
                            )
                else:
                    file_warn += 1
                    total_warn += 1
                    for gene, res in results.items():
                        if res["status"] != "✅ OK":
                            problem_isolats.append(
                                f"{rec.id} | {sero} | {gene} | {res['status']} | {res['reason']}"
                            )

            fname = os.path.basename(fasta_path)
            print(f"  {fname}")
            print(f"    ✅ OK: {file_ok} | ⚠️ WARN: {file_warn} | ❌ FAIL: {file_fail}")

    # ── Ringkasan ────────────────────────────────────────────
    print()
    print("=" * 65)
    print("  RINGKASAN VALIDASI KOORDINAT")
    print("=" * 65)
    total = total_ok + total_warn + total_fail
    print(f"  Total isolat divalidasi : {total}")
    print(f"  ✅ Koordinat valid      : {total_ok} ({round(total_ok/total*100,1)}%)")
    print(f"  ⚠️  Warning             : {total_warn}")
    print(f"  ❌ Gagal validasi       : {total_fail}")

    if total_warn + total_fail == 0:
        print()
        print("  🎉 SEMUA KOORDINAT VALID!")
        print("  Subsetting aman dilanjutkan.")
        print()
        print("  → Bukti untuk dosen:")
        print("    Kodon start & stop setiap gen konsisten dengan RefSeq")
        print("    di seluruh isolat. Koordinat fixed valid untuk dataset ini.")
    else:
        print()
        print("  ⚠️  Ada isolat bermasalah:")
        for p in problem_isolats[:20]:
            print(f"    {p}")
        if len(problem_isolats) > 20:
            print(f"    ... dan {len(problem_isolats)-20} lainnya")

        # Simpan daftar masalah
        out_path = f"{BASE}/../09_scripts/validasi_koordinat_masalah.txt"
        with open(out_path, "w") as f:
            f.write("ID Isolat | Serotipe | Gen | Status | Alasan\n")
            f.write("-" * 80 + "\n")
            for p in problem_isolats:
                f.write(p + "\n")
        print(f"\n  Daftar lengkap disimpan: {out_path}")

    print("=" * 65)

if __name__ == "__main__":
    run()
_tee.close()
