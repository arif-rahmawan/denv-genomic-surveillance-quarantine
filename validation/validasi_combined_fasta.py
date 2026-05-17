"""
validasi_combined_fasta.py
Validasi 12 file FASTA gabungan sebelum MSA di UGENE:
file exist | jumlah | RefSeq | panjang | duplikat | spasi di label
Simpan di: /Users/arifrahmawan/Downloads/REVISI 5/10_scripts/validasi/
Jalankan : python3 validasi_combined_fasta.py
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

_tee = _Tee("validasi_combined_fasta")



from Bio import SeqIO
import os

# ── KONFIGURASI ───────────────────────────────────────────────────────────────
ALIGNDIR = "/Users/arifrahmawan/Downloads/REVISI 5/03_alignment"

EXPECTED_COUNT = {
    "DENV1": 218,   # 217 isolat + 1 RefSeq
    "DENV2": 369,   # 368 isolat + 1 RefSeq
    "DENV3": 270,   # 269 isolat + 1 RefSeq
    "DENV4": 239,   # 238 isolat + 1 RefSeq
}

EXPECTED_LENGTH = {
    "DENV1": {"prM": 498,  "E": 1485, "NS1": 1056},
    "DENV2": {"prM": 498,  "E": 1485, "NS1": 1056},
    "DENV3": {"prM": 498,  "E": 1479, "NS1": 1056},
    "DENV4": {"prM": 498,  "E": 1485, "NS1": 1056},
}

REFSEQ_PREFIX = "REFSEQ_NC_"

# ── FUNGSI ────────────────────────────────────────────────────────────────────

def cek_jumlah(records, expected):
    return len(records) == expected, len(records)

def cek_refseq(records):
    """RefSeq harus ada di posisi pertama dengan prefix REFSEQ_NC_"""
    if not records:
        return False, "file kosong"
    first_id = records[0].id
    if first_id.startswith(REFSEQ_PREFIX):
        return True, first_id
    return False, f"posisi pertama: {first_id}"

def cek_panjang(records, expected):
    salah = [(r.id, len(r.seq)) for r in records if len(r.seq) != expected]
    return len(salah) == 0, salah

def cek_duplikat(records):
    ids      = [r.id for r in records]
    duplikat = [i for i in set(ids) if ids.count(i) > 1]
    return len(duplikat) == 0, duplikat

def cek_spasi(records):
    """Label tidak boleh ada spasi — UGENE/IQ-TREE akan truncate di spasi pertama"""
    bermasalah = [r.id for r in records if " " in r.id]
    return len(bermasalah) == 0, bermasalah

# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  VALIDASI COMBINED FASTA (siap MSA)")
    print("  Cek: file | jumlah | RefSeq | panjang | duplikat | spasi")
    print("=" * 60)

    total_pass = 0
    total_fail = 0
    gagal_list = []

    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        exp_count = EXPECTED_COUNT[sero]
        print(f"\n  {sero} (expected: {exp_count} sekuens per gen)")
        print("  " + "-" * 50)

        for gene in ["E", "NS1", "prM"]:
            label    = f"{sero}_{gene}_combined.fasta"
            path     = os.path.join(ALIGNDIR, sero, label)
            expected = EXPECTED_LENGTH[sero][gene]
            status   = []
            lulus    = True

            # 1. File exist
            if not os.path.exists(path):
                print(f"  [X] {label}: FILE TIDAK DITEMUKAN")
                total_fail += 1
                gagal_list.append((label, "file tidak ditemukan"))
                continue

            records = list(SeqIO.parse(path, "fasta"))

            # 2. Jumlah
            ok, count = cek_jumlah(records, exp_count)
            status.append(f"jumlah=OK({count})" if ok else f"jumlah=GAGAL({count})")
            if not ok: lulus = False

            # 3. RefSeq di posisi pertama
            ok, info = cek_refseq(records)
            status.append(f"refseq=OK" if ok else f"refseq=GAGAL({info})")
            if not ok: lulus = False

            # 4. Panjang
            ok, salah = cek_panjang(records, expected)
            status.append(f"panjang=OK({expected}bp)" if ok else f"panjang=GAGAL({len(salah)} salah)")
            if not ok: lulus = False

            # 5. Duplikat label
            ok, duplikat = cek_duplikat(records)
            status.append("duplikat=OK" if ok else f"duplikat=GAGAL({len(duplikat)})")
            if not ok: lulus = False

            # 6. Spasi di label
            ok, bermasalah = cek_spasi(records)
            status.append("spasi=OK" if ok else f"spasi=GAGAL({len(bermasalah)} label)")
            if not ok: lulus = False

            icon = "[v]" if lulus else "[X]"
            print(f"  {icon} {label}")
            print(f"     -> {' | '.join(status)}")

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
        print("\n  Detail GAGAL:")
        for nama, ket in gagal_list:
            print(f"    [X] {nama}: {ket}")
    else:
        print("\n  Semua file valid. Siap diinput ke UGENE untuk MSA.")

    print("=" * 60)

if __name__ == "__main__":
    main()
_tee.close()
