"""
validasi_subset_isolat.py
Validasi ketat hasil subset isolat: file exist, jumlah, panjang, duplikat, translasi, komposisi
IUPAC ambiguity code (Y, R, M, W, S, K, B, D, H, V) -> WARNING, bukan GAGAL
Simpan di: /Users/arifrahmawan/Downloads/REVISI 5/10_scripts/validasi/
Jalankan : python3 validasi_subset_isolat.py
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

_tee = _Tee("validasi_subset_isolat")



from Bio import SeqIO
from Bio.Seq import Seq
import os

# -- KONFIGURASI --------------------------------------------------------------
SUBSETDIR = "/Users/arifrahmawan/Downloads/REVISI 5/02_subsetting"

EXPECTED_COUNT = {
    "DENV1": 217,
    "DENV2": 368,
    "DENV3": 269,
    "DENV4": 238,
}

EXPECTED_LENGTH = {
    "DENV1": {"prM": 498, "E": 1485, "NS1": 1056},
    "DENV2": {"prM": 498, "E": 1485, "NS1": 1056},
    "DENV3": {"prM": 498, "E": 1479, "NS1": 1056},
    "DENV4": {"prM": 498, "E": 1485, "NS1": 1056},
}

VALID_BASES = set("ATCGN")
IUPAC_EXTRA = set("YRWSKBDHVM")   # kode ambigu valid IUPAC
VALID_ALL   = VALID_BASES | IUPAC_EXTRA

# -- FUNGSI -------------------------------------------------------------------

def cek_jumlah(records, expected):
    return len(records) == expected, len(records)

def cek_panjang(records, expected):
    salah = [r.id for r in records if len(r.seq) != expected]
    return len(salah) == 0, salah

def cek_duplikat(records):
    ids      = [r.id for r in records]
    duplikat = [i for i in ids if ids.count(i) > 1]
    return len(duplikat) == 0, list(set(duplikat))

def cek_translasi(records):
    gagal = []
    for r in records:
        try:
            protein = Seq(str(r.seq)).translate()
            if str(protein[:-1]).count("*") > 0:
                gagal.append(r.id)
        except Exception:
            gagal.append(r.id)
    return len(gagal) == 0, gagal

def cek_komposisi(records):
    """
    - Karakter di luar VALID_ALL -> ERROR (benar-benar rusak)
    - Karakter IUPAC_EXTRA       -> WARNING (ambigu tapi valid)
    """
    errors   = []
    warnings = []
    for r in records:
        chars = set(str(r.seq).upper())
        bad   = chars - VALID_ALL
        iupac = chars & IUPAC_EXTRA
        if bad:
            errors.append((r.id, bad))
        elif iupac:
            warnings.append((r.id, iupac))
    return errors, warnings

# -- MAIN ---------------------------------------------------------------------

def main():
    print("=" * 60)
    print("  VALIDASI SUBSET ISOLAT")
    print("  Cek: file | jumlah | panjang | duplikat | translasi | komposisi")
    print("  IUPAC ambiguity code -> WARNING (bukan GAGAL)")
    print("=" * 60)

    total_pass    = 0
    total_fail    = 0
    total_warning = 0
    gagal_list    = []
    warning_list  = []

    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        exp_count = EXPECTED_COUNT[sero]
        print(f"\n  {sero} (expected: {exp_count} sekuens per gen)")
        print("  " + "-" * 50)

        for gene in ["prM", "E", "NS1"]:
            label    = f"{sero}_{gene}_subset.fasta"
            path     = os.path.join(SUBSETDIR, sero, label)
            expected = EXPECTED_LENGTH[sero][gene]
            status   = []
            lulus    = True
            ada_warn = False

            # 1. File exist
            if not os.path.exists(path):
                print(f"  [X] {label}: FILE TIDAK DITEMUKAN")
                total_fail += 1
                gagal_list.append((label, "file tidak ditemukan"))
                continue

            records = list(SeqIO.parse(path, "fasta"))

            # 2. Jumlah
            ok, count = cek_jumlah(records, exp_count)
            status.append(f"jumlah=OK({count})" if ok else f"jumlah=GAGAL({count}, exp={exp_count})")
            if not ok: lulus = False

            # 3. Panjang
            ok, salah = cek_panjang(records, expected)
            status.append(f"panjang=OK({expected}bp)" if ok else f"panjang=GAGAL({len(salah)} salah)")
            if not ok: lulus = False

            # 4. Duplikat
            ok, duplikat = cek_duplikat(records)
            status.append("duplikat=OK" if ok else f"duplikat=GAGAL({len(duplikat)} ID ganda)")
            if not ok: lulus = False

            # 5. Translasi
            ok, gagal_trans = cek_translasi(records)
            status.append("translasi=OK" if ok else f"translasi=GAGAL({len(gagal_trans)} stop codon)")
            if not ok: lulus = False

            # 6. Komposisi
            errors, warnings = cek_komposisi(records)
            if errors:
                status.append(f"komposisi=GAGAL({len(errors)} karakter rusak)")
                lulus = False
            elif warnings:
                kode_unik = set(k for _, kode in warnings for k in kode)
                status.append(f"komposisi=WARN({len(warnings)} sekuens IUPAC: {kode_unik})")
                ada_warn = True
            else:
                status.append("komposisi=OK")

            # Print hasil
            if not lulus:
                icon = "[X]"
            elif ada_warn:
                icon = "[W]"
            else:
                icon = "[v]"

            print(f"  {icon} {label}")
            print(f"     -> {' | '.join(status)}")

            if lulus:
                total_pass += 1
                if ada_warn:
                    total_warning += 1
                    warning_list.append((label, warnings))
            else:
                total_fail += 1
                gagal_list.append((label, " | ".join(status)))

    # Ringkasan
    print("\n" + "=" * 60)
    print("  RINGKASAN VALIDASI")
    print("=" * 60)
    print(f"  LULUS   : {total_pass} / 12")
    print(f"  WARNING : {total_warning} / 12  (IUPAC ambiguity - tetap lulus)")
    print(f"  GAGAL   : {total_fail} / 12")

    if warning_list:
        print("\n  Detail WARNING (IUPAC - tidak perlu tindakan):")
        for nama, warns in warning_list:
            kode = set(k for _, kode in warns for k in kode)
            print(f"    [W] {nama}: {len(warns)} sekuens mengandung {kode}")

    if gagal_list:
        print("\n  Detail GAGAL:")
        for nama, ket in gagal_list:
            print(f"    [X] {nama}: {ket}")
    else:
        print("\n  Tidak ada kegagalan. Semua subset isolat siap untuk MSA.")

    print("=" * 60)

if __name__ == "__main__":
    main()
_tee.close()
