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

_tee = _Tee("validasi_subsetting_refseq")

from Bio import SeqIO
from Bio.Seq import Seq
import os

BASE        = "/Users/arifrahmawan/Downloads/REVISI 5"
SUBSETDIR   = f"{BASE}/02_subsetting"
PROTEIN_REF = f"{BASE}/01_raw_data/REFSEQ/REFSEQ_PROTEIN_ALL.fasta"

PROTEIN_MAP = {
    "NP_722460.2":    ("DENV1", "E"),
    "NP_722461.1":    ("DENV1", "NS1"),
    "NP_733807.2":    ("DENV1", "prM"),
    "NP_739583.2":    ("DENV2", "E"),
    "NP_739584.2":    ("DENV2", "NS1"),
    "NP_739582.2":    ("DENV2", "prM"),
    "YP_001531168.2": ("DENV3", "E"),
    "YP_001531169.2": ("DENV3", "NS1"),
    "YP_001531166.1": ("DENV3", "prM"),
    "NP_740317.1":    ("DENV4", "E"),
    "NP_740318.1":    ("DENV4", "NS1"),
    "NP_740315.1":    ("DENV4", "prM"),
}

def load_protein_refseq():
    proteins = {}
    for rec in SeqIO.parse(PROTEIN_REF, "fasta"):
        acc = rec.id.split("|")[0].strip()
        if acc in PROTEIN_MAP:
            sero, gene = PROTEIN_MAP[acc]
            proteins[(sero, gene)] = str(rec.seq).upper()
    return proteins

def translate_nukleotida(nt_seq):
    try:
        protein = str(Seq(nt_seq).translate())
        if protein.endswith("*"):
            protein = protein[:-1]
        return protein
    except Exception:
        return None

def compare_detail(seq1, seq2):
    if len(seq1) != len(seq2):
        return 0, max(len(seq1), len(seq2)), 0.0, []
    total = len(seq1)
    mismatches = []
    match_count = 0
    for i, (a, b) in enumerate(zip(seq1, seq2)):
        if a == b:
            match_count += 1
        else:
            mismatches.append((i + 1, a, b))
    return match_count, total, (match_count / total) * 100, mismatches

def main():
    print("=" * 70)
    print("  VALIDASI SUBSETTING REFSEQ — PERBANDINGAN POSISI PER POSISI")
    print("=" * 70)
    proteins = load_protein_refseq()
    print(f"\n  Protein RefSeq NCBI dimuat: {len(proteins)} / 12\n")
    total_pass = 0
    total_fail = 0
    hasil = []
    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        print(f"  {'─'*62}")
        print(f"  {sero}")
        print(f"  {'─'*62}")
        for gene in ["E", "NS1", "prM"]:
            nt_path = os.path.join(SUBSETDIR, sero, f"{sero}_{gene}_refseq.fasta")
            label   = f"{sero}_{gene}"
            records = list(SeqIO.parse(nt_path, "fasta"))
            nt_seq       = str(records[0].seq).upper()
            protein_kita = translate_nukleotida(nt_seq)
            protein_ncbi = proteins[(sero, gene)]
            match, total_pos, identity, mismatches = compare_detail(protein_kita, protein_ncbi)
            lulus = identity == 100.0
            mid = len(protein_kita) // 2
            if lulus:
                print(f"  [v] {label}")
                print(f"       Posisi identik   : {match}/{total_pos} aa ({identity:.2f}%)")
                print(f"       Spot N-term (1-5): kita={protein_kita[:5]} | NCBI={protein_ncbi[:5]} OK")
                print(f"       Spot tengah      : kita={protein_kita[mid-2:mid+3]} | NCBI={protein_ncbi[mid-2:mid+3]} OK")
                print(f"       Spot C-term (-5) : kita={protein_kita[-5:]} | NCBI={protein_ncbi[-5:]} OK")
                total_pass += 1
            else:
                print(f"  [X] {label}: identitas {identity:.2f}% — {len(mismatches)} posisi berbeda")
                total_fail += 1
            hasil.append({"label": label, "nt_bp": len(nt_seq), "match": match, "total": total_pos, "identity": identity, "lulus": lulus})
            print()
    print("=" * 70)
    print("  RINGKASAN")
    print("=" * 70)
    print(f"  {'Gen':<12} {'nt (bp)':<10} {'Posisi Identik':<20} {'Identitas':<12} Status")
    print("  " + "-" * 62)
    for h in hasil:
        print(f"  {h['label']:<12} {h['nt_bp']:<10} {h['match']}/{h['total']:<16} {h['identity']:.2f}%      {'LULUS' if h['lulus'] else 'GAGAL'}")
    print(f"\n  LULUS : {total_pass} / 12  |  GAGAL : {total_fail} / 12")
    if total_fail == 0:
        print("\n  KESIMPULAN: Seluruh 12 gen-serotipe identik 100% pada")
        print("  setiap posisi AA. Subsetting nukleotida terbukti AKURAT.")
    print("=" * 70)

if __name__ == "__main__":
    main()
_tee.close()
