import csv
import random
from itertools import count
from pathlib import Path
from importlib.resources import files

import numpy as np
import pandas as pd

from vdjrearranger.barcode import load_whitelist, random_umi, sample_barcodes
from vdjrearranger.chemistry import get_chemistry
from vdjrearranger.distributions import CountSampler
from vdjrearranger.models import Cell, Clonotype, ReadPair


class Simulator:
    """
    Core engine responsible for generating synthetic droplet-based V(D)J sequencing reads.
    Arranges data per the hierarchy of cells -> chains -> UMIs -> reads (and sequencing artifacts).
    """

    def __init__(
            self,
            clonotypes: list[Clonotype],
            chemistry: str = "SC5P-R2",
            cell_distribution: CountSampler = None,
            umi_distribution: CountSampler = None,
            reads_distribution: CountSampler = None,
            chain_ratio: float = 0.5,
            chain_ratio_variance: float = 0,
            seed: int = 1234,
            noise_cells: int = 100,
            noise_umis: int = 1,
            coverage_mode: str = "normal",
            fragment_center: float = 0.5,
            fragment_std: float = 0.2,
            fragment_dropoff: float = 0.1,
            chain1_name: str = "chain1",
            chain2_name: str = "chain2",
            sequencer: str = "VDJREARRANGER",
            flowcell: str = "SIMULATEDFLOWCELL",
            primers_tsv: str | None = None,
    ):
        """
        Initializes the simulation environment and mathematical distributions.

        :param clonotypes: list[Clonotype], List of Clonotype instances to be sampled.
        :param chemistry: str, Profile dictating barcode, UMI, and read lengths.
        :param cell_distribution: CountSampler, Sampler determining the number of cells per clonotype.
        :param umi_distribution: CountSampler, Sampler determining the number of unique transcripts per cell.
        :param reads_distribution: CountSampler, Sampler determining PCR duplicates per UMI.
        :param chain_ratio: float, Baseline expression division between Chain 1 and Chain 2.
        :param chain_ratio_variance: float, Cell-to-cell variance applied to the chain ratio.
        :param seed: int, Random seed for deterministic simulation.
        :param noise_cells: int, Number of background barcodes to assign ambient RNA.
        :param noise_umis: int, Maximum number of ambient UMIs per background droplet.
        :param coverage_mode: str, Coverage distribution mode ('random', 'normal', '3prime', 'tiling').
        :param fragment_center: float, Proportional center of the fragmentation window (0.0-1.0).
        :param fragment_std: float, Proportional standard deviation of the window.
        :param fragment_dropoff: float, Proportional exponential drop-off scale.
        :param chain1_name: str, Column identifier for the first chain.
        :param chain2_name: str, Column identifier for the second chain.
        :param sequencer: str, Sequencer ID string used in FASTQ headers.
        :param flowcell: str, Flowcell ID string used in FASTQ headers.
        :param primers_tsv: str | None, Optional external path to override bundled trimming primers.
        """
        self.clonotypes = clonotypes
        self.chemistry = get_chemistry(chemistry)
        self.cell_distribution = cell_distribution
        self.umi_distribution = umi_distribution
        self.reads_distribution = reads_distribution

        self.chain_ratio = chain_ratio
        self.chain_ratio_variance = chain_ratio_variance
        self.noise_cells = noise_cells
        self.noise_umis = max(1, noise_umis)

        self.coverage_mode = coverage_mode.lower()
        self.fragment_center = max(0.0, min(1.0, fragment_center))
        self.fragment_std = max(0.01, fragment_std)
        self.fragment_dropoff = max(0.01, fragment_dropoff)

        self.chain1_name = chain1_name
        self.chain2_name = chain2_name

        self.sequencer = sequencer
        self.flowcell = flowcell

        self.primers = self._load_primers(primers_tsv)

        random.seed(seed)
        np.random.seed(seed)

        self.whitelist = load_whitelist()
        self.quality_pool = self._precompute_quality_profiles(max_len=150)

    def _load_primers(self, tsv_path: str | None) -> dict[str, str]:
        """
        Loads inner sequence primers from a custom TSV or the bundled package data.

        :param tsv_path: str | None, Optional path to a custom TSV file.
        :return: A dictionary mapping chain IDs to their constant region primer sequences.
        """
        if tsv_path and Path(tsv_path).exists():
            df = pd.read_csv(tsv_path, sep="\t")
        else:
            bundled_path = files("vdjrearranger.data").joinpath("primers.tsv")
            try:
                df = pd.read_csv(bundled_path, sep="\t")
            except FileNotFoundError:
                return {}

        if "chain" not in df.columns or "primer_seq" not in df.columns:
            return {}
        return dict(zip(df["chain"].str.upper(), df["primer_seq"]))

    def _precompute_quality_profiles(self, max_len: int = 150) -> list[str]:
        """
        Pre-generates a pool of simulated quality score strings to reduce runtime overhead.

        :param max_len: int, Maximum read length to generate qualities for.
        :return: A list of pre-formatted ASCII Phred+33 strings.
        """
        profiles = []
        for _ in range(1000):
            is_noise = random.random() < 0.1
            start_q = 40
            end_q = 20 if is_noise else 30
            quals = []
            for i in range(max_len):
                progress = i / max(1, max_len - 1)
                mean_q = start_q - (start_q - end_q) * progress
                final_q = max(0, min(41, int(mean_q)))
                quals.append(chr(final_q + 33))
            profiles.append("".join(quals))
        return profiles

    def _generate_quality_scores(self, length: int) -> str:
        """
        Retrieves a pre-computed quality score string cut to the target length.

        :param length: int, Desired length of the output string.
        :return: Phred+33 quality string.
        """
        return random.choice(self.quality_pool)[:length]

    def build_cells(self) -> list[Cell]:
        """
        Constructs the physical representation of the droplets by assigning barcodes.

        :return: A list of configured Cell objects.
        """
        all_cells_list = []
        for clonotype in self.clonotypes:
            n_cells = self.cell_distribution.sample()
            for _ in range(n_cells):
                all_cells_list.append(clonotype.name)

        all_cells_list.extend(["NOISE"] * self.noise_cells)
        barcodes = sample_barcodes(len(all_cells_list), self.whitelist)

        return [Cell(barcode=barcodes[i], clonotype_name=name) for i, name in enumerate(all_cells_list)]

    def trim_primer(self, chain_name: str, seq: str) -> str:
        """
        Identifies and removes the constant region internal primer sequence from a transcript.

        :param chain_name: str, The identifier of the chain to determine which primer to look for.
        :param seq: str, The full nucleotide sequence of the chain.
        :return: The trimmed nucleotide sequence.
        """
        primer = None
        upper_chain = chain_name.upper()

        for key in self.primers:
            if key in upper_chain:
                primer = self.primers[key]
                break

        if primer:
            primer_idx = seq.find(primer)
            if primer_idx != -1:
                return seq[:primer_idx + len(primer)]
        return seq

    def get_tiling_starts(self, seq_len: int, read_len: int, jitter: int = 0) -> list[int]:
        """
        Calculates starting indices for overlapping tiles in both 5'->3' and 3'->5' directions.

        :param seq_len: int, Length of the target nucleotide sequence.
        :param read_len: int, Length of the sequence read (tile size).
        :param jitter: int, Random offset to apply to the starting positions to vary coverage.
        :return: A list of unique integer starting indices.
        """
        max_start = seq_len - read_len
        if max_start <= 0:
            return [0]

        step = max(1, read_len // 2)

        # 5' to 3' forward pass
        forward_starts = list(range(jitter, max_start + 1, step))
        if not forward_starts or forward_starts[-1] + read_len < seq_len:
            forward_starts.append(max_start)

        # 3' to 5' backward pass
        backward_starts = list(range(max_start - jitter, -1, -step))
        if not backward_starts or backward_starts[-1] > 0:
            backward_starts.append(0)

        combined = forward_starts + backward_starts

        # Deduplicate coordinates while preserving the bidirectional order
        seen = set()
        deduped = []
        for s in combined:
            s_clamped = max(0, min(max_start, s))
            if s_clamped not in seen:
                seen.add(s_clamped)
                deduped.append(s_clamped)

        return deduped

    def fragment_sequence(self, seq: str, start_idx: int | None = None) -> str:
        """
        Simulates library fragmentation by extracting a targeted window from the transcript.

        :param seq: str, Parent nucleotide sequence.
        :param start_idx: int | None, Optional forced starting index (used exclusively for tiling mode).
        :return: A fragmented, oriented read sequence tailored to the chemistry.
        """
        read_len = min(self.chemistry.r2_length, len(seq))
        max_start = len(seq) - read_len

        if max_start <= 0:
            start = 0
        elif start_idx is not None:
            # Force the specific start index but strictly clamp it to valid biological bounds
            start = max(0, min(max_start, start_idx))
        else:
            if self.coverage_mode == "random":
                start = random.randint(0, max_start)

            elif self.coverage_mode == "3prime":
                start = -1
                attempts = 0
                scale_val = self.fragment_dropoff * max_start
                while (start < 0 or start > max_start) and attempts < 10:
                    shift = int(np.random.exponential(scale=scale_val))
                    start = max_start - shift
                    attempts += 1

                if start < 0: start = 0
                if start > max_start: start = max_start

            else:
                start = -1
                attempts = 0
                center_val = self.fragment_center * max_start
                std_val = self.fragment_std * max_start
                while (start < 0 or start > max_start) and attempts < 10:
                    start = int(np.random.normal(center_val, std_val))
                    attempts += 1

                if start < 0: start = 0
                if start > max_start: start = max_start

        insert = seq[start: start + read_len]

        if random.random() < 0.01:
            output_insert = insert
        else:
            rc_map = str.maketrans("ACGTN", "TGCAN")
            output_insert = insert.translate(rc_map)[::-1]

        return output_insert

    def _create_read_pair(
            self,
            read_counter: count,
            barcode: str,
            umi: str,
            seq: str,
            chain_name: str,
            cell_name: str,
            start_idx: int | None = None
    ) -> dict:
        """
        Utility method to finalize read simulation and package it into a dictionary format.

        :param read_counter: count, Global iterator generating sequential integer read IDs.
        :param barcode: str, The cell barcode string.
        :param umi: str, The UMI string.
        :param seq: str, The transcript sequence assigned to this iteration.
        :param chain_name: str, The originating chain column name.
        :param cell_name: str, The unique identifier for the host cell/clonotype.
        :param start_idx: int | None, Optional forced starting index for fragmentation.
        :return: Dictionary containing the populated ReadPair object and metadata.
        """
        fragment = self.fragment_sequence(seq, start_idx=start_idx)
        read_number = next(read_counter)
        read_id = f"{self.sequencer}:1:{self.flowcell}:1:1101:{20000 + read_number}:{1000 + read_number}"

        return {
            "read_pair": ReadPair(
                read_id=read_id,
                r1_seq=barcode + umi,
                r2_seq=fragment,
                r1_qual=self._generate_quality_scores(len(barcode + umi)),
                r2_qual=self._generate_quality_scores(len(fragment)),
            ),
            "barcode": barcode,
            "umi": umi,
            "clonotype": cell_name,
            "chain": chain_name,
        }

    def generate_read_pairs(self):
        """
        Iterates through the physical hierarchy yielding read dictionaries.

        :return: Generator yielding paired read dictionaries.
        """
        read_counter = count(1)
        cells = self.build_cells()
        clonotype_lookup = {c.name: c for c in self.clonotypes}

        for cell in cells:
            barcode = cell.barcode[: self.chemistry.barcode_length]
            is_noise_cell = (cell.clonotype_name == "NOISE")

            # 1. Determine which sequences and how many UMIs to generate for this cell
            if is_noise_cell:
                random_clonotype = random.choice(self.clonotypes)
                chain_name = self.chain1_name if random.random() < 0.5 else self.chain2_name
                raw_seq = random_clonotype.chain1 if chain_name == self.chain1_name else random_clonotype.chain2

                n_umis = random.randint(1, self.noise_umis)
                chains_to_process = [(chain_name, raw_seq, n_umis)]
            else:
                clonotype = clonotype_lookup[cell.clonotype_name]
                ratio = np.random.normal(self.chain_ratio, self.chain_ratio_variance)
                ratio = min(max(ratio, 0.05), 0.95)

                total_umis = self.umi_distribution.sample()
                c1_umis = max(2, int(total_umis * ratio))
                c2_umis = max(2, total_umis - c1_umis)

                chains_to_process = [
                    (self.chain1_name, clonotype.chain1, c1_umis),
                    (self.chain2_name, clonotype.chain2, c2_umis)
                ]

            # 2. Process each assigned chain and determine read start positions
            for chain_name, raw_seq, target_umis in chains_to_process:
                seq = self.trim_primer(chain_name, raw_seq)
                starts = []

                if self.coverage_mode == "tiling":
                    # Generate one complete bidirectional pass
                    base_starts = self.get_tiling_starts(len(seq), self.chemistry.r2_length, jitter=0)
                    starts.extend(base_starts)

                    # If target UMIs exceeds base tiles, loop with random jitter until threshold is met
                    while len(starts) < target_umis:
                        jitter = random.randint(1, max(1, self.chemistry.r2_length // 2))
                        starts.extend(self.get_tiling_starts(len(seq), self.chemistry.r2_length, jitter=jitter))

                    # Slice to ensure exact target count, but never output less than one complete tile set
                    final_count = max(len(base_starts), target_umis)
                    starts = starts[:final_count]
                else:
                    # Non-tiling modes allow dynamic random sampling later
                    starts = [None] * target_umis

                # 3. Generate final reads for every determined position
                for start_idx in starts:
                    umi = random_umi()[: self.chemistry.umi_length]
                    n_reads = 1 if is_noise_cell else max(2, self.reads_distribution.sample())

                    for _ in range(n_reads):
                        yield self._create_read_pair(
                            read_counter, barcode, umi, seq, chain_name, cell.clonotype_name, start_idx=start_idx
                        )


def upload_clonotypes(
        tsv_path: str,
        clonotype_col: str,
        chain1_col: str,
        chain2_col: str,
) -> list[Clonotype]:
    """
    Parses the input TSV file to extract sequence pairs into an internal representation.

    :param tsv_path: str, Path to the clonotype TSV.
    :param clonotype_col: str, Column header containing the core sequence name.
    :param chain1_col: str, Column header for the primary sequence chain.
    :param chain2_col: str, Column header for the secondary sequence chain.
    :return: A list of populated Clonotype dataclasses.
    """
    df = pd.read_csv(tsv_path, sep="\t", quoting=csv.QUOTE_NONE, dtype=str)
    clonotypes = []
    for i, row in df.iterrows():
        base_name = str(row[clonotype_col])
        unique_name = f"{base_name}_row{i}"
        clonotypes.append(Clonotype(
            name=unique_name,
            chain1=str(row[chain1_col]),
            chain2=str(row[chain2_col]),
        ))
    return clonotypes
