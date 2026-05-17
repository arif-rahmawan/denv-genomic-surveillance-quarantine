"""
validasi_refseq_subset.py
Validasi ketat hasil subset RefSeq: file exist, jumlah, panjang, header, translasi, komposisi
Simpan di: /Users/arifrahmawan/Downloads/REVISI 5/10_scripts/validasi/
Jalankan : python3 validasi_refseq_subset.py
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

_tee = _Tee("validasi_refseq_subset")



from Bio import SeqIO
from Bio.Seq import Seq
import os

# ── KONFIGURASI ──────────────────────────────────────────────────────────────
SUBSETDIR = "/Users/arifrahmawan/Downloads/REVISI 5/02_subsetting"

REFSEQ_ACC = {
    "DENV1": "NC_001477.1",
    "DENV2": "NC_001474.2",
    "DENV3": "NC_001475.2",
    "DENV4": "NC_002640.1",
}

EXPECTED_LENGTH = {
    "DENV1": {"prM": 498, "E": 1485, "NS1": 1056},
    "DENV2": {"prM": 498, "E": 1485, "NS1": 1056},
    "DENV3": {"prM": 498, "E": 1479, "NS1": 1056},
    "DENV4": {"prM": 498, "E": 1485, "NS1": 1056},
}

VALID_BASES = set("ATCGN")

# ── FUNGSI ────────────────────────────────────────────────────────────────────

def cek_file_exist(path):
    return os.path.exists(path)

def cek_jumlah(records):
    return len(records) == 1

def cek_panjang(seq, expected):
    return len(seq) == expected

def cek_header(rec_id, expected_acc):
    return expected_acc in rec_id

def cek_translasi(seq_str):
    try:
        protein = Seq(seq_str).translate()
        return str(protein[:-1]).count("*") == 0
    except Exception:
        return False

def cek_komposisi(seq_str):
    chars = set(seq_str.upper())
    invalid = chars - VALID_BASES
    return len(invalid) == 0, invalid

# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  VALIDASI REFSEQ SUBSET")
    print("  Cek: file | jumlah | panjang | header | translasi | komposisi")
    print("=" * 60)

    total_pass = 0
    total_fail = 0
    gagal_list = []

    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        print(f"\n  {sero} ({REFSEQ_ACC[sero]})")
        print("  " + "─" * 50)

        for gene in ["prM", "E", "NS1"]:
            label    = f"{sero}_{gene}_refseq.fasta"
            path     = os.path.join(SUBSETDIR, sero, label)
            expected = EXPECTED_LENGTH[sero][gene]
            status   = []
            lulus    = True

            # 1. File exist
            if not cek_file_exist(path):
                print(f"  ❌ {label}: FILE TIDAK DITEMUKAN")
                total_fail += 1
                gagal_list.append((label, "file tidak ditemukan"))
                continue

            records = list(SeqIO.parse(path, "fasta"))
            seq_str = str(records[0].seq).upper() if records else ""

            # 2. Jumlah sekuens
            if cek_jumlah(records):
                status.append("jumlah=✅")
            else:
                status.append(f"jumlah=❌({len(records)})")
                lulus = False

            # 3. Panjang
            if cek_panjang(seq_str, expected):
                status.append(f"panjang=✅({len(seq_str)}bp)")
            else:
                status.append(f"panjang=❌({len(seq_str)}bp, exp={expected}bp)")
                lulus = False

            # 4. Header
            if cek_header(records[0].id, REFSEQ_ACC[sero]):
                status.append("header=✅")
            else:
                status.append(f"header=❌({records[0].id})")
                lulus = False

            # 5. Translasi
            if cek_translasi(seq_str):
                status.append("translasi=✅")
            else:
                status.append("translasi=❌(stop codon prematur)")
                lulus = False

            # 6. Komposisi
            ok_komp, invalid = cek_komposisi(seq_str)
            if ok_komp:
                status.append("komposisi=✅")
            else:
                status.append(f"komposisi=❌(karakter invalid: {invalid})")
                lulus = False

            # Print baris hasil
            icon = "✅" if lulus else "❌"
            print(f"  {icon} {label}")
            print(f"     → {' | '.join(status)}")

            if lulus:
                total_pass += 1
            else:
                total_fail += 1
                gagal_list.append((label, " | ".join(status)))

    # Ringkasan
    print("\n" + "=" * 60)
    print("  RINGKASAN VALIDASI")
    print("=" * 60)
    print(f"  LULUS : {total_pass} / 12")
    print(f"  GAGAL : {total_fail} / 12")

    if gagal_list:
        print("\n  Detail yang GAGAL:")
        for nama, ket in gagal_list:
            print(f"    ❌ {nama}: {ket}")
    else:
        print("\n  ✅ Semua RefSeq subset valid. Siap digunakan untuk MSA.")

    print("=" * 60)

if __name__ == "__main__":
    main()
_tee.close()
