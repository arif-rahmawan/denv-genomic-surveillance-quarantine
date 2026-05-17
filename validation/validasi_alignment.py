"""
validasi_alignment.py
Tahap 7 — QC Alignment Nukleotida
Cek: jumlah sekuens, panjang konsisten, RefSeq posisi pertama,
     kolom 100% gap, delta kolom before/after
Simpan di: /Users/arifrahmawan/Downloads/REVISI 5/10_scripts/validasi/
Jalankan : python3 validasi_alignment.py
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

_tee = _Tee("validasi_alignment")



from Bio import SeqIO
import os

BASE      = "/Users/arifrahmawan/Downloads/REVISI 5"
ALIGNDIR  = f"{BASE}/03_alignment"
COMBINEDIR= f"{BASE}/03_alignment"

EXPECTED_COUNT = {
    "DENV1": 218, "DENV2": 369, "DENV3": 270, "DENV4": 239,
}

EXPECTED_RAW = {
    "DENV1": {"E": 1485, "NS1": 1056, "prM": 498},
    "DENV2": {"E": 1485, "NS1": 1056, "prM": 498},
    "DENV3": {"E": 1479, "NS1": 1056, "prM": 498},
    "DENV4": {"E": 1485, "NS1": 1056, "prM": 498},
}

REFSEQ_PREFIX = "REFSEQ_NC_"

def cek_jumlah(records, expected):
    return len(records) == expected, len(records)

def cek_panjang_konsisten(records):
    lengths = set(len(r.seq) for r in records)
    return len(lengths) == 1, lengths

def cek_refseq(records):
    if not records:
        return False, "kosong"
    return records[0].id.startswith(REFSEQ_PREFIX), records[0].id

def hitung_gap_100(records):
    if not records:
        return 0, []
    aln_len   = len(records[0].seq)
    gap_cols  = []
    for i in range(aln_len):
        kolom = [str(r.seq[i]).upper() for r in records]
        if all(c == '-' for c in kolom):
            gap_cols.append(i)
    return len(gap_cols), gap_cols

def main():
    print("=" * 70)
    print("  VALIDASI ALIGNMENT NUKLEOTIDA (Tahap 7 QC)")
    print("  Cek: jumlah | panjang | RefSeq | kolom 100% gap | delta")
    print("=" * 70)

    total_pass  = 0
    total_fail  = 0
    total_warn  = 0
    gagal_list  = []
    ringkasan   = []

    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        exp_count = EXPECTED_COUNT[sero]
        print(f"\n  {sero} (expected: {exp_count} sekuens)")
        print("  " + "-" * 62)

        for gene in ["E", "NS1", "prM"]:
            label    = f"{sero}_{gene}"
            aln_path = os.path.join(ALIGNDIR, sero, f"{sero}_{gene}_alignment.fasta")
            raw_len  = EXPECTED_RAW[sero][gene]
            status   = []
            lulus    = True
            ada_warn = False

            if not os.path.exists(aln_path):
                print(f"  [X] {label}: FILE TIDAK DITEMUKAN")
                total_fail += 1
                gagal_list.append((label, "file tidak ditemukan"))
                continue

            records = list(SeqIO.parse(aln_path, "fasta"))
            aln_len = len(records[0].seq) if records else 0
            delta   = aln_len - raw_len

            # 1. Jumlah
            ok, count = cek_jumlah(records, exp_count)
            status.append(f"jumlah=OK({count})" if ok else f"jumlah=GAGAL({count})")
            if not ok: lulus = False

            # 2. Panjang konsisten
            ok, lengths = cek_panjang_konsisten(records)
            status.append(f"panjang=OK({aln_len}bp)" if ok else f"panjang=GAGAL({lengths})")
            if not ok: lulus = False

            # 3. RefSeq posisi pertama
            ok, ref_id = cek_refseq(records)
            status.append("refseq=OK" if ok else f"refseq=GAGAL({ref_id})")
            if not ok: lulus = False

            # 4. Kolom 100% gap
            n_gap, gap_cols = hitung_gap_100(records)
            if n_gap == 0:
                status.append("gap100%=OK(0)")
            else:
                status.append(f"gap100%=WARN({n_gap} kolom)")
                ada_warn = True

            # 5. Delta kolom
            status.append(f"delta=+{delta}")

            # Print
            if not lulus:
                icon = "[X]"
            elif ada_warn:
                icon = "[W]"
            else:
                icon = "[v]"

            print(f"  {icon} {label}")
            print(f"     -> {' | '.join(status)}")
            if n_gap > 0:
                print(f"     -> posisi gap 100%: {gap_cols[:5]}{'...' if len(gap_cols)>5 else ''}")

            ringkasan.append({
                "label": label, "raw": raw_len, "aln": aln_len,
                "delta": delta, "gap100": n_gap, "lulus": lulus,
            })

            if lulus:
                total_pass += 1
                if ada_warn: total_warn += 1
            else:
                total_fail += 1
                gagal_list.append((label, " | ".join(status)))

    # Ringkasan tabel
    print("\n" + "=" * 70)
    print("  RINGKASAN QC ALIGNMENT")
    print("=" * 70)
    print(f"  {'Gen':<12} {'Raw (bp)':<12} {'Aligned (bp)':<15} {'Delta':<10} {'Gap 100%':<10} Status")
    print("  " + "-" * 65)
    for r in ringkasan:
        st = "LULUS" if r["lulus"] else "GAGAL"
        print(f"  {r['label']:<12} {r['raw']:<12} {r['aln']:<15} +{r['delta']:<9} {r['gap100']:<10} {st}")

    print(f"\n  LULUS   : {total_pass} / 12")
    print(f"  WARNING : {total_warn} / 12  (kolom gap 100%)")
    print(f"  GAGAL   : {total_fail} / 12")

    if total_fail == 0 and total_warn == 0:
        print("\n  Semua alignment valid. Tidak ada kolom 100% gap.")
        print("  -> Siap untuk Tahap 8 (RDP4) dan Tahap 9 (IQ-TREE).")
    elif total_fail == 0 and total_warn > 0:
        print("\n  Alignment valid. Ada kolom 100% gap — perlu dipertimbangkan")
        print("  apakah akan dihapus sebelum IQ-TREE (lihat posisi di atas).")
    else:
        print("\n  Ada kegagalan — cek detail di atas.")

    print("=" * 70)

if __name__ == "__main__":
    main()
_tee.close()
