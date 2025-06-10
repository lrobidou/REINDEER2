import sys
import scipy.stats as st
import matplotlib.pyplot as plt
import numpy as np

def parse_subpart(part):
    if ":" in part:
        rang, ab = part.split(":")
        intab = 0 if ab == "*" else int(ab)
        low, high = rang.split("-")
        return [intab] * (int(high) - int(low) + 1)
    elif part == "*":
        return [0]
    else:
        print("aled, c'était pas un *")
        return []

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compare REINDEER and REINDEER2 abundance estimates.")
    parser.add_argument("reindeer_results", help="REINDEER output file")
    parser.add_argument("reindeer2_results", help="REINDEER2 output file")
    parser.add_argument("--genes_or_kmers", required=True, choices=["genes", "k-mers"],
                        help="Comparison unit (genes or k-mers)")
    parser.add_argument("--sparse_or_dense", required=True, choices=["sparse", "dense"],
                        help="REINDEER2 index type")
    parser.add_argument("--file_format", required=True, choices=["png", "pdf"],
                        help="Output file format")
    parser.add_argument("--log_transform", action="store_true",
                        help="Apply log2 transform to abundances")

    args = parser.parse_args()

    summary = {}

    # Parse REINDEER (original)
    with open(args.reindeer_results, "r") as res:
        next(res)
        for line in res:
            header, ab_parts = line.strip().split("\t", maxsplit=1)
            for color, ab_part in enumerate(ab_parts.split("\t")):
                all_vals = []
                if ',' in ab_part:
                    for sub_part in ab_part.split(','):
                        all_vals.extend(parse_subpart(sub_part))
                else:
                    all_vals.extend(parse_subpart(ab_part))
                vals = [v for v in all_vals if v > 0]
                if len(all_vals) - len(vals) > 0.5:
                    median = 0
                else:
                    vals.sort()
                    mid = len(vals) // 2
                    median = (vals[mid] + vals[mid-1]) // 2 if len(vals) % 2 == 0 else vals[mid]
                key = f"{header}___{color}"
                summary[key] = [median, 0]

    # Parse REINDEER2 (sparse/dense)
    with open(args.reindeer2_results, "r") as res:
        next(res)
        for line in res:
            header_part, color, ab_part = line.strip().split(",")
            ab = int(ab_part)
            header = header_part[1:].split(' ')[0]
            key = f"{header}___{color}"
            if key in summary:
                summary[key][1] = ab
            else:
                print(f"Warning: {key} not found in REINDEER results")

    abs_r, abs_r2, diffs = [], [], []

    for key in summary:
        ab_r, ab_r2 = summary[key]
        if ab_r2 == 0 and ab_r > 0:
            print(f"[FN] {key}: REINDEER={ab_r} vs REINDEER2={ab_r2}")

        abs_r.append(ab_r)
        abs_r2.append(ab_r2)

        diff = (ab_r - ab_r2) / ab_r if ab_r > 0 else (ab_r - ab_r2)
        diffs.append(diff)

    pearson, pvalue = st.pearsonr(abs_r, abs_r2)
    diff_mean = np.mean(diffs)
    diff_median = np.median(diffs)

    print("\n===== RESULTS =====")
    print(f"Correlations between {args.reindeer_results} and {args.reindeer2_results}\n")
    print(f"Min: {round(min(diffs), 2)}\nMax: {round(max(diffs), 2)}")
    print(f"Mean: {round(diff_mean, 2)}\nMedian: {round(diff_median, 2)}")
    print(f"Pearson: {round(pearson, 2)}\np-value: {round(pvalue, 4)}")

    # === Plotting ===
    x = np.array(abs_r2)
    y = np.array(abs_r)

    if args.log_transform:
        x = np.where(x == 0, 1, x)
        y = np.where(y == 0, 1, y)
        x_plot = np.log2(x)
        y_plot = np.log2(y)
    else:
        x_plot = x
        y_plot = y

    plt.figure(figsize=(8, 6))
    plt.scatter(x_plot, y_plot, s=5, alpha=0.5, color="black")

    # Linear regression line
    slope, intercept = np.polyfit(x_plot, y_plot, 1)
    plt.plot(x_plot, slope * x_plot + intercept, color="blue", linewidth=1.5, label="Linear fit")

    # Diagonal reference
    min_val = min(x_plot.min(), y_plot.min())
    max_val = max(x_plot.max(), y_plot.max())
    plt.plot([min_val, max_val], [min_val, max_val], color="red", linestyle="--", linewidth=1)

    # Annotate Pearson r
    plt.text(0.05, 0.95, f"r = {round(pearson, 4)}", transform=plt.gca().transAxes,
             fontsize=14, color="darkred", ha="left", va="top", fontweight="bold")

    # Labels and title
    log_suffix = " (log2)" if args.log_transform else ""
    plt.xlabel(f"REINDEER2 abundances (with {args.sparse_or_dense} index){log_suffix}", fontsize=14)
    plt.ylabel(f"REINDEER abundances{log_suffix}", fontsize=14)
    plt.title(f"Comparison for estimated abundances of {args.genes_or_kmers}", fontsize=16, fontweight="bold")

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, linestyle=":", linewidth=0.5)
    plt.tight_layout()

    # Save plot
    output_filename = f"rd1_rd2_{args.sparse_or_dense}_{args.genes_or_kmers}.{args.file_format}"
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved to {output_filename}")

