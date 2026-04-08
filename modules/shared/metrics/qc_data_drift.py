#!/usr/bin/env python3
"""
qc_data_drift.py — Insert size distribution drift detector for DeepVariant fine-tuning.

Compares current sample's insert size distribution (gatk_metrics.insert_size_metrics)
against a benchmark distribution (data used to train DV) using the KS test.
Outputs a plot and raises a WARNING signal if drift is detected.

Usage:
    python3 qc_data_drift.py \\
        --sample  gatk_metrics.insert_size_metrics \\
        --benchmark benchmark.gatk_metrics.insert_size_metrics \\
        --out     insert_size_drift \\
        [--ks-pvalue 0.05]

DeepVariant features monitored here: insert_size
"""

import argparse
import sys
import warnings
from typing import Tuple
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp


# ---------------------------------------------------------------------------
# Picard insert_size_metrics parser
# ---------------------------------------------------------------------------

def parse_insert_size_metrics(path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse a Picard CollectInsertSizeMetrics output file.

    Returns
    -------
    metrics : pd.DataFrame  — summary statistics (one row per orientation)
    histogram : pd.DataFrame — columns [insert_size, count]
    """
    with open(path) as fh:
        lines = fh.readlines()

    # Locate the two sections: METRICS and HISTOGRAM
    metrics_header_idx = None
    histogram_header_idx = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("MEDIAN_INSERT_SIZE"):
            metrics_header_idx = i
        if stripped.startswith("insert_size"):
            histogram_header_idx = i

    if metrics_header_idx is None:
        raise ValueError(f"Could not find METRICS section in {path}")
    if histogram_header_idx is None:
        raise ValueError(f"Could not find HISTOGRAM section in {path}")

    # Parse metrics table (header + data rows until blank line)
    metrics_rows = []
    for line in lines[metrics_header_idx:]:
        if not line.strip():
            break
        metrics_rows.append(line.rstrip('\n\r').split("\t"))
    metrics = pd.DataFrame(metrics_rows[1:], columns=metrics_rows[0])

    # Parse histogram table
    hist_rows = []
    for line in lines[histogram_header_idx:]:
        if not line.strip():
            break
        hist_rows.append(line.rstrip('\n\r').split("\t"))
    histogram = pd.DataFrame(hist_rows[1:], columns=hist_rows[0])
    histogram = histogram.rename(columns={histogram.columns[0]: "insert_size",
                                          histogram.columns[1]: "count"})
    histogram["insert_size"] = histogram["insert_size"].astype(int)
    histogram["count"] = histogram["count"].astype(float)

    return metrics, histogram


def histogram_to_array(histogram: pd.DataFrame, max_insert: int = 1500) -> np.ndarray:
    """Expand histogram counts into a flat array for KS test (capped at max_insert)."""
    hist = histogram[histogram["insert_size"] <= max_insert].copy()
    counts = hist["count"].astype(int).values
    sizes = hist["insert_size"].values
    return np.repeat(sizes, counts)


# ---------------------------------------------------------------------------
# KS test
# ---------------------------------------------------------------------------

def run_ks_test(arr_sample: np.ndarray, arr_benchmark: np.ndarray,
                pvalue_threshold: float) -> dict:
    stat, pvalue = ks_2samp(arr_sample, arr_benchmark)
    drift_detected = pvalue < pvalue_threshold
    return {"ks_statistic": stat, "pvalue": pvalue,
            "threshold": pvalue_threshold, "drift_detected": drift_detected}


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_distributions(hist_sample: pd.DataFrame, hist_benchmark: pd.DataFrame,
                       ks_result: dict, out_prefix: str, max_insert: int = 1500):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # --- Left: overlaid histograms (normalised) ---
    ax = axes[0]
    for hist, label, color in [
        (hist_benchmark, "Benchmark (DV train)", "#2196F3"),
        (hist_sample,    "Current sample",       "#F44336"),
    ]:
        h = hist[hist["insert_size"] <= max_insert]
        total = h["count"].sum()
        ax.plot(h["insert_size"], h["count"] / total,
                label=label, color=color, lw=1.5, alpha=0.8)
    ax.set_xlabel("Insert size (bp)")
    ax.set_ylabel("Frequency")
    ax.set_title("Insert size distribution")
    ax.legend(loc="upper left")

    # Annotate KS result (bottom-right, away from legend)
    drift_color = "#c62828" if ks_result["drift_detected"] else "#2e7d32"
    drift_text  = "DRIFT DETECTED — consider fine-tuning DV" \
                  if ks_result["drift_detected"] else "No significant drift"
    ax.text(0.97, 0.95,
            f"KS stat = {ks_result['ks_statistic']:.4f}\n"
            f"p-value = {ks_result['pvalue']:.2e}\n"
            f"{drift_text}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            color=drift_color,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=drift_color, lw=1.2))

    # --- Right: difference plot ---
    ax2 = axes[1]
    merged = pd.merge(
        hist_benchmark[hist_benchmark["insert_size"] <= max_insert]
                      .rename(columns={"count": "bench"}),
        hist_sample[hist_sample["insert_size"] <= max_insert]
                   .rename(columns={"count": "sample"}),
        on="insert_size", how="outer"
    ).fillna(0).sort_values("insert_size")
    bench_norm  = merged["bench"]  / merged["bench"].sum()
    sample_norm = merged["sample"] / merged["sample"].sum()
    diff = sample_norm - bench_norm
    colors = ["#F44336" if d > 0 else "#2196F3" for d in diff]
    ax2.bar(merged["insert_size"], diff, color=colors, width=1.0, alpha=0.7)
    ax2.axhline(0, color="black", lw=0.8)
    ax2.set_xlabel("Insert size (bp)")
    ax2.set_ylabel("Δ Frequency (sample − benchmark)")
    ax2.set_title("Distribution difference\n(red = sample higher, blue = benchmark higher)")

    plt.tight_layout()
    out_path = f"{out_prefix}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[qc_data_drift] Plot saved → {out_path}")


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_report(ks_result: dict, metrics_sample: pd.DataFrame,
                 metrics_benchmark: pd.DataFrame, out_prefix: str):
    out_path = f"{out_prefix}.txt"
    lines = [
        "=== Insert Size Drift Report ===",
        f"KS statistic : {ks_result['ks_statistic']:.6f}",
        f"p-value      : {ks_result['pvalue']:.4e}",
        f"Threshold    : {ks_result['threshold']}",
        f"Drift        : {'YES — fine-tune DeepVariant recommended' if ks_result['drift_detected'] else 'NO'}",
        "",
        "--- Benchmark metrics ---",
    ]
    lines += metrics_benchmark.to_string(index=False).splitlines()
    lines += ["", "--- Sample metrics ---"]
    lines += metrics_sample.to_string(index=False).splitlines()

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[qc_data_drift] Report saved → {out_path}")


# ---------------------------------------------------------------------------
# Signal
# ---------------------------------------------------------------------------

def raise_drift_signal(ks_result: dict):
    """Print a clear WARNING and exit with code 1 if drift is detected."""
    if ks_result["drift_detected"]:
        msg = (
            "\n"
            "╔══════════════════════════════════════════════════════════════╗\n"
            "║  DATA DRIFT DETECTED — DeepVariant fine-tuning recommended  ║\n"
            f"║  KS statistic = {ks_result['ks_statistic']:.4f}   "
            f"p-value = {ks_result['pvalue']:.2e}  (threshold={ks_result['threshold']})  ║\n"
            "╚══════════════════════════════════════════════════════════════╝\n"
        )
        print(msg, file=sys.stderr)
        sys.exit(1)
    else:
        print(
            f"[qc_data_drift] No significant drift detected "
            f"(KS p={ks_result['pvalue']:.4e} >= {ks_result['threshold']}). "
            "DeepVariant fine-tuning not required.",
            file=sys.stderr,
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect insert size distribution drift for DeepVariant QC."
    )
    parser.add_argument("--sample",    required=True,
                        help="Current sample: gatk_metrics.insert_size_metrics")
    parser.add_argument("--benchmark", required=True,
                        help="Benchmark: benchmark.gatk_metrics.insert_size_metrics")
    parser.add_argument("--out",       default="insert_size_drift",
                        help="Output prefix for .png and .txt (default: insert_size_drift)")
    parser.add_argument("--ks-pvalue", type=float, default=0.05,
                        dest="ks_pvalue",
                        help="KS test p-value threshold for drift (default: 0.05)")
    parser.add_argument("--max-insert", type=int, default=1500,
                        dest="max_insert",
                        help="Truncate histogram at this insert size (default: 1500)")
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"[qc_data_drift] Parsing sample    : {args.sample}")
    metrics_sample, hist_sample = parse_insert_size_metrics(args.sample)

    print(f"[qc_data_drift] Parsing benchmark  : {args.benchmark}")
    metrics_bench, hist_bench = parse_insert_size_metrics(args.benchmark)

    arr_sample = histogram_to_array(hist_sample, args.max_insert)
    arr_bench  = histogram_to_array(hist_bench,  args.max_insert)

    if arr_sample.size == 0 or arr_bench.size == 0:
        raise ValueError("Histogram arrays are empty — check input files.")

    print(f"[qc_data_drift] Running KS test (threshold p < {args.ks_pvalue})")
    ks_result = run_ks_test(arr_sample, arr_bench, args.ks_pvalue)

    plot_distributions(hist_sample, hist_bench, ks_result, args.out, args.max_insert)
    write_report(ks_result, metrics_sample, metrics_bench, args.out)
    raise_drift_signal(ks_result)


if __name__ == "__main__":
    main()