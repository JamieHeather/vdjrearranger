from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def generate_knee_plot(truth_tsv_path: Path, outdir: Path) -> Path:
    """
    Generates a log-log Barcode Rank vs UMI Count plot (knee plot) to visually assess the simulated separation
    between real cells and ambient background noise (analagous to the cellranger web summary plot).

    :param truth_tsv_path: str, path to the generated truth reads summary TSV (or TSV.gz).
    :param outdir: str, output directory where the resulting plot PNG will be saved.
    :return: str, full file path to the generated plot image.
    """
    df = pd.read_csv(truth_tsv_path, sep="\t")

    umi_counts = df.groupby("barcode")["umi"].nunique().reset_index()
    umi_counts = umi_counts.rename(columns={"umi": "umi_count"})

    # Rank barcodes descending by total unique UMIs
    umi_counts = umi_counts.sort_values("umi_count", ascending=False).reset_index(drop=True)
    umi_counts["rank"] = umi_counts.index + 1

    labels = df[["barcode", "clonotype"]].drop_duplicates(subset=["barcode"])
    plot_df = umi_counts.merge(labels, on="barcode", how="left")

    plot_df["cell_type"] = np.where(plot_df["clonotype"] == "NOISE", "Background", "Simulated Cell")

    plt.figure(figsize=(8, 6))

    bg_data = plot_df[plot_df["cell_type"] == "Background"]
    plt.plot(
        bg_data["rank"], bg_data["umi_count"],
        color="gray", alpha=0.5, label="Background (Noise)"
    )

    cell_data = plot_df[plot_df["cell_type"] == "Simulated Cell"]
    plt.scatter(
        cell_data["rank"], cell_data["umi_count"],
        color="blue", s=15, zorder=5, label="True Cells"
    )

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Barcodes (Ranked)")
    plt.ylabel("Unique UMI Count")
    plt.title("Simulated Barcode Rank Plot")
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.2)

    plot_path = outdir / "rank_umi_plot.png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    return plot_path
