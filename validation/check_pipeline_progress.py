"""
step09_check_progress.py
Arif Rahmawan | DENV In Silico Pipeline

Tujuan:
    Monitor progress IQ-TREE3 — cek output mana yang sudah selesai,
    dan ekstrak model substitusi dari setiap .iqtree file.

Cara pakai:
    python step09_check_progress.py

Output:
    - Tabel status 12 file (+ DENV4_E_noRecomb)
    - Ringkasan model substitusi (untuk Tahap 10)
    - Copy-paste siap masuk ke dokumen BAB III/IV
"""

import os
import re
from pathlib import Path

BASE = "/Users/arifrahmawan/Downloads/REVISI 5"
OUT_DIR = f"{BASE}/07_phylogenetic"

# 12 file utama + 1 extra
RUNS = [
    ("DENV1", "DENV1_E"),
    ("DENV1", "DENV1_NS1"),
    ("DENV1", "DENV1_prM"),
    ("DENV2", "DENV2_E"),
    ("DENV2", "DENV2_NS1"),
    ("DENV2", "DENV2_prM"),
    ("DENV3", "DENV3_E"),
    ("DENV3", "DENV3_NS1"),
    ("DENV3", "DENV3_prM"),
    ("DENV4", "DENV4_E"),
    ("DENV4", "DENV4_NS1"),
    ("DENV4", "DENV4_prM"),
    ("DENV4", "DENV4_E_noRecomb"),   # run ekstra
]

EXPECTED_FILES = [".treefile", ".iqtree", ".log", ".contree", ".mldist"]


def extract_model_info(iqtree_file):
    """Ekstrak model terbaik dan log-likelihood dari .iqtree file."""
    info = {
        "best_model": "—",
        "log_likelihood": "—",
        "bootstrap_note": "—",
        "tree_length": "—",
    }
    if not os.path.isfile(iqtree_file):
        return info

    try:
        with open(iqtree_file, "r") as f:
            content = f.read()

        # Best-fit model
        m = re.search(r"Best-fit model according to BIC:\s+(\S+)", content)
        if m:
            info["best_model"] = m.group(1)

        # Log-likelihood
        m = re.search(r"Log-likelihood of the tree:\s+([-\d.]+)", content)
        if m:
            info["log_likelihood"] = m.group(1)

        # Bootstrap
        m = re.search(r"Number of bootstrap replicates:\s+(\d+)", content)
        if m:
            info["bootstrap_note"] = f"UFBoot {m.group(1)}x"

        # Tree length
        m = re.search(r"Total tree length \(sum of branch lengths\):\s+([\d.]+)", content)
        if m:
            info["tree_length"] = m.group(1)

    except Exception as e:
        info["best_model"] = f"[error: {e}]"

    return info


def main():
    print("=" * 70)
    print("  TAHAP 9 — STATUS IQ-TREE3")
    print("=" * 70)
    print("")

    results = []
    all_done = True

    for serotype, label in RUNS:
        prefix = os.path.join(OUT_DIR, serotype, label)
        file_status = {}
        for ext in EXPECTED_FILES:
            path = prefix + ext
            file_status[ext] = os.path.isfile(path)

        done = all(file_status.values())
        if not done:
            all_done = False

        model_info = extract_model_info(prefix + ".iqtree") if file_status[".iqtree"] else {}

        results.append({
            "label": label,
            "serotype": serotype,
            "done": done,
            "files": file_status,
            "model": model_info,
        })

    # ---- STATUS TABLE ----
    print(f"{'No':<4} {'Label':<22} {'Status':<12} {'Best Model':<25} {'Log-L':<15}")
    print("-" * 82)

    for i, r in enumerate(results, 1):
        status_str = "✅ SELESAI" if r["done"] else "⬜ BELUM"
        model_str = r["model"].get("best_model", "—") if r["done"] else "—"
        logl_str  = r["model"].get("log_likelihood", "—") if r["done"] else "—"
        label_str = r["label"]
        if "noRecomb" in label_str:
            label_str = f"  [EXTRA] {label_str}"
        print(f"{i:<4} {label_str:<22} {status_str:<12} {model_str:<25} {logl_str:<15}")

    print("")
    done_count = sum(1 for r in results if r["done"])
    total_count = len(results)
    print(f"Progress: {done_count}/{total_count} selesai")
    print("")

    # ---- RINGKASAN MODEL UNTUK TAHAP 10 ----
    completed = [r for r in results if r["done"] and "noRecomb" not in r["label"]]
    if completed:
        print("=" * 70)
        print("  RINGKASAN MODEL SUBSTITUSI (untuk Tahap 10 & Penulisan)")
        print("=" * 70)
        print("")
        for r in completed:
            m = r["model"]
            print(f"  {r['label']:<20} Model: {m.get('best_model','—'):<25} "
                  f"Tree length: {m.get('tree_length','—')}")
        print("")

        # Format tabel untuk dokumen
        print("--- Copy-paste ke BAB IV (format tabel) ---")
        print(f"{'Gen-Serotipe':<20} | {'Model Terbaik (BIC)':<25} | {'Log-Likelihood':<15} | {'Bootstrap':<12}")
        print("-" * 80)
        for r in completed:
            m = r["model"]
            print(f"{r['label']:<20} | {m.get('best_model','—'):<25} | "
                  f"{m.get('log_likelihood','—'):<15} | {m.get('bootstrap_note','—'):<12}")
        print("")

    # ---- CEK FILE YANG BELUM SELESAI ----
    not_done = [r for r in results if not r["done"]]
    if not_done:
        print("File yang belum selesai:")
        for r in not_done:
            missing = [ext for ext, exists in r["files"].items() if not exists]
            print(f"  - {r['label']}: missing {missing}")
        print("")
        print("Untuk melanjutkan, jalankan ulang step09_iqtree3_batch.sh")
        print("IQ-TREE3 tidak akan overwrite file yang sudah ada.")
        print("(Kalau mau jalankan ulang dari nol, hapus dulu output prefix)")
    else:
        print("Semua selesai! Lanjut ke Tahap 10 (dokumentasi model substitusi)")
        print("dan Tahap 11 (Fisher's Exact Test).")


if __name__ == "__main__":
    main()
