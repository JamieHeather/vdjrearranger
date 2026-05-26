from dataclasses import dataclass


@dataclass
class Clonotype:
    """
    Represents a single unique paired-chain sequence variant.

    :param name: str, unique string identifier for the clonotype.
    :param chain1: str, nt sequence for the primary/light chain (e.g. TRA, IGL, TRG).
    :param chain2: str, nt sequence for the secondary/heavy chain (e.g. TRB, IGH, TRD).
    """
    name: str
    chain1: str
    chain2: str


@dataclass
class Cell:
    """
    Represents an individual simulated cell or droplet.

    :param barcode: str, 10x-style cell barcode associated with this droplet.
    :param clonotype_name: str, identifier of the clonotype contained within.
    """
    barcode: str
    clonotype_name: str


@dataclass
class ReadPair:
    """
    Represents a single simulated paired-end sequence read.

    :param read_id: str, formatted FASTQ sequence identifier.
    :param r1_seq: str, nt sequence for Read 1 (e.g. barcode + UMI for SC5P-R2).
    :param r2_seq: str, nt sequence for Read 2 (e.g. transcript fragment for SC5P-R2).
    :param r1_qual: str, phred+33 quality string for Read 1.
    :param r2_qual: str, phred+33 quality string for Read 2.
    """
    read_id: str
    r1_seq: str
    r2_seq: str
    r1_qual: str
    r2_qual: str
