import gzip
import shutil
import sys
import time
from pathlib import Path

import pandas as pd
import typer

from vdjrearranger.distributions import CountSampler
from vdjrearranger.fastq import FastqWriter
from vdjrearranger.logging_utils import write_run_log
from vdjrearranger.simulation import Simulator, upload_clonotypes
from vdjrearranger.plot import generate_knee_plot

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@app.command()
def main(
        input_tsv: str = typer.Option(..., help="Input TSV containing paired receptor sequences"),
        outdir: str = typer.Option(..., help="Output directory for FASTQs, logs, and truth files"),
        clonotype_col: str = typer.Option("TCR_name",
            help="Column containing clonotype identifiers"),
        chain1_col: str = typer.Option("TRA_nt",
            help="Column containing primary chain nucleotide sequences (e.g., TRA, IGL, TRG)"),
        chain2_col: str = typer.Option("TRB_nt",
            help="Column containing secondary chain nucleotide sequences (e.g., TRB, IGH, TRD)"),
        coverage_mode: str = typer.Option("tiling", help="Coverage distribution mode (random|normal|3prime|tiling)"),
        primers_tsv: str = typer.Option(None, help="Optional path to custom TSV overriding bundled inner primers"),
        lane: int = typer.Option(1, help="Illumina lane number used in FASTQ filenames"),
        sequencer: str = typer.Option("VDJREARRANGER", help="Sequencer ID used in generated FASTQ headers"),
        flowcell: str = typer.Option("SIMULATEDFLOWCELL", help="Flowcell ID used in generated FASTQ headers"),
        cells_mode: str = typer.Option("fixed",
            help="Distribution mode for cells per clonotype (fixed|normal|poisson|powerlaw)"),
        cells_value: int = typer.Option(1, help="Magnitude parameter for cells per clonotype"),
        cells_variance: float = typer.Option(2.0, help="Variance parameter for cells per clonotype"),
        umis_mode: str = typer.Option("fixed",
            help="Distribution mode for UMIs per cell (fixed|normal|poisson|powerlaw)"),
        umis_value: int = typer.Option(100, help="Magnitude parameter for UMIs per cell (across both chains)"),
        umis_variance: float = typer.Option(5.0, help="Variance parameter for UMIs per cell"),
        reads_mode: str = typer.Option("fixed",
            help="Distribution mode for reads per UMI (fixed|normal|poisson|powerlaw)"),
        reads_value: int = typer.Option(50, help="Magnitude parameter for reads per UMI"),
        reads_variance: float = typer.Option(2.0, help="Variance parameter for reads per UMI"),
        chain_ratio: float = typer.Option(0.5, help="Baseline expression ratio between chain 1 and chain 2"),
        chain_variance: float = typer.Option(0.0, help="Variance applied to the chain expression ratio per cell"),
        seed: int = typer.Option(1234, help="Random seed for deterministic runs"),
        chemistry: str = typer.Option("SC5P-R2", help="10x Genomics library chemistry profile"),
        sample_name: str = typer.Option("simulated_vdj", help="Sample prefix used for FASTQ generation"),
        noise_cells: int = typer.Option(100, help="Number of background noise cell barcodes"),
        noise_umis: int = typer.Option(1, help="Max UMIs per noise cell (ambient RNA)"),
        fragment_center: float = typer.Option(0.5,
            help="Proportional center position of fragment for normal mode (0.0 to 1.0)"),
        fragment_std: float = typer.Option(0.2,
            help="Proportional standard deviation of fragment center for normal mode"),
        fragment_dropoff: float = typer.Option(0.1,
            help="Proportional exponential scale parameter for 3prime mode dropoff"),
):
    """
    Main entry point for generating simulated single-cell V(D)J sequencing outputs.
    Initialises necessary directory, configures probability distributions,
    runs the simulation, and compiles end-of-run summaries.

    :param input_tsv: str, Path to input TSV.
    :param outdir: str, Path to output directory.
    :param clonotype_col: str, Name of clonotype column.
    :param chain1_col: str, Name of chain 1 sequence column.
    :param chain2_col: str, Name of chain 2 sequence column.
    :param primers_tsv: str | None, Path to custom primers TSV.
    :param lane: int, Illumina lane integer.
    :param sequencer: str, Sequencer ID for FASTQ.
    :param flowcell: str, Flowcell ID for FASTQ.
    :param cells_mode: str, Mode for cell distribution.
    :param cells_value: int, Magnitude for cell distribution.
    :param cells_variance: float, Variance for cell distribution.
    :param umis_mode: str, Mode for UMI distribution.
    :param umis_value: int, Magnitude for UMI distribution.
    :param umis_variance: float, Variance for UMI distribution.
    :param reads_mode: str, Mode for read distribution.
    :param reads_value: int, Magnitude for read distribution.
    :param reads_variance: float, Variance for read distribution.
    :param chain_ratio: float, Expression ratio between chains.
    :param chain_variance: float, Variance applied to expression ratio.
    :param seed: int, Random seed.
    :param chemistry: str, 10x chemistry profile.
    :param sample_name: str, Sample prefix for output FASTQs.
    :param noise_cells: int, Number of noise cells.
    :param noise_umis: int, UMIs per noise cell.
    :param coverage_mode: str, Coverage sampling algorithm to use.
    :param fragment_center: float, Proportional center for normal coverage.
    :param fragment_std: float, Proportional deviation for normal coverage.
    :param fragment_dropoff: float, Proportional dropoff for 3prime coverage.
    """
    start_time = time.time()

    outdir = Path(outdir)
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    with open(outdir / "run_command.txt", "w") as f:
        f.write(" ".join(sys.argv) + "\n")

    clonotypes = upload_clonotypes(
        input_tsv,
        clonotype_col,
        chain1_col,
        chain2_col,
    )

    sim = Simulator(
        clonotypes=clonotypes,
        chemistry=chemistry,
        noise_cells=noise_cells,
        noise_umis=noise_umis,
        coverage_mode=coverage_mode,
        fragment_center=fragment_center,
        fragment_std=fragment_std,
        fragment_dropoff=fragment_dropoff,
        chain1_name=chain1_col,
        chain2_name=chain2_col,
        sequencer=sequencer,
        flowcell=flowcell,
        primers_tsv=primers_tsv,
        cell_distribution=CountSampler(
            mode=cells_mode,
            value=cells_value,
            variance=cells_variance,
        ),
        umi_distribution=CountSampler(
            mode=umis_mode,
            value=umis_value,
            variance=umis_variance,
        ),
        reads_distribution=CountSampler(
            mode=reads_mode,
            value=reads_value,
            variance=reads_variance,
        ),
        chain_ratio=chain_ratio,
        chain_ratio_variance=chain_variance,
        seed=seed,
    )

    writer = FastqWriter(
        outdir=outdir,
        sample_name=sample_name,
        lane=lane,
    )

    truth_rows = []
    total_reads = 0

    for result in sim.generate_read_pairs():
        writer.write_pair(result["read_pair"])

        truth_rows.append({
            "read_id": result["read_pair"].read_id,
            "barcode": result["barcode"],
            "umi": result["umi"],
            "clonotype": result["clonotype"],
            "chain": result["chain"],
        })

        total_reads += 1

    writer.close()

    truth_path = outdir / "truth_reads.tsv.gz"
    df = pd.DataFrame(truth_rows)
    with gzip.open(truth_path, "wt") as f:
        df.to_csv(f, sep="\t", index=False)

    summary_path = outdir / "summary_stats.txt"
    with open(summary_path, "w") as f:
        f.write("=== Simulation Summary Statistics ===\n\n")

        if not df.empty:
            real_df = df[df["clonotype"] != "NOISE"]
            noise_df = df[df["clonotype"] == "NOISE"]

            real_reads = len(real_df)
            noise_reads = len(noise_df)
            real_cells = real_df["barcode"].nunique()

            f.write(f"Total Reads (All): {len(df)}\n")
            f.write(f"Total Reads (Real Cells): {real_reads}\n")
            f.write(f"Total Reads (Noise/Ambient): {noise_reads}\n")
            f.write(f"Total Estimated Cells: {real_cells}\n")

            mean_reads = real_reads / real_cells if real_cells > 0 else 0
            f.write(f"Mean Reads per Cell: {mean_reads:.2f}\n\n")

            if real_cells > 0:
                umi_counts = real_df.groupby(["barcode", "chain"])["umi"].nunique().reset_index()
                umi_pivot = umi_counts.pivot(index="barcode", columns="chain", values="umi").fillna(0)

                med_chain1 = umi_pivot[chain1_col].median() if chain1_col in umi_pivot.columns else 0
                med_chain2 = umi_pivot[chain2_col].median() if chain2_col in umi_pivot.columns else 0

                f.write(f"Median {chain1_col} UMIs per Cell: {med_chain1:.2f}\n")
                f.write(f"Median {chain2_col} UMIs per Cell: {med_chain2:.2f}\n\n")

                f.write("=== Cells per Clonotype ===\n")
                cells_per_clono = real_df.groupby("clonotype")["barcode"].nunique().sort_values(ascending=False)
                for clono, count in cells_per_clono.items():
                    f.write(f"{clono}: {count}\n")
        else:
            f.write("No reads generated.\n")

    typer.echo(f"Generated summary statistics at {summary_path}")

    try:
        plot_path = generate_knee_plot(truth_path, outdir)
        typer.echo(f"Generated rank-UMI plot at {plot_path}")
    except Exception as e:
        typer.echo(f"Warning: Could not generate rank-UMI plot: {e}")

    params_to_log = {
        "input_tsv": input_tsv,
        "outdir": str(outdir),
        "clonotype_col": clonotype_col,
        "chain1_col": chain1_col,
        "chain2_col": chain2_col,
        "primers_tsv": primers_tsv,
        "lane": lane,
        "sequencer": sequencer,
        "flowcell": flowcell,
        "cells_mode": cells_mode,
        "cells_value": cells_value,
        "cells_variance": cells_variance,
        "umis_mode": umis_mode,
        "umis_value": umis_value,
        "umis_variance": umis_variance,
        "reads_mode": reads_mode,
        "reads_value": reads_value,
        "reads_variance": reads_variance,
        "chain_ratio": chain_ratio,
        "chain_variance": chain_variance,
        "seed": seed,
        "chemistry": chemistry,
        "sample_name": sample_name,
        "noise_cells": noise_cells,
        "noise_umis": noise_umis,
        "coverage_mode": coverage_mode,
        "fragment_center": fragment_center,
        "fragment_std": fragment_std,
        "fragment_dropoff": fragment_dropoff,
    }

    write_run_log(outdir, params_to_log, start_time)


if __name__ == "__main__":
    app()
