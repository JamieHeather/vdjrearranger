import gzip
import random
from importlib.resources import files


def load_whitelist() -> list[str]:
    """
    Loads the standard 10x 737K cell barcode whitelist bundled within the data module.

    :return: list of valid barcode strings.
    """
    bundled_path = files("vdjrearranger.data").joinpath(
        "737K-august-2016.txt.gz"
    )

    with gzip.open(bundled_path, "rt") as handle:
        return [line.strip() for line in handle]


def sample_barcodes(n: int, whitelist: list[str]) -> list[str]:
    """
    Randomly selects a unique set of barcodes from the whitelist population.

    :param n: int, number of unique barcodes to sample.
    :param whitelist: list of str, population of valid barcodes to draw from.
    :return: list of uniquely sampled barcodes.
    """
    return random.sample(whitelist, n)


def random_umi(length: int = 10) -> str:
    """
    Generates a random sequence of nucleotides to act as a Unique Molecular Identifier.

    :param length: int, length of the sequence to generate (default=10)
    :return: A random string composed of A, C, G, and T.
    """
    return "".join(random.choices("ACGT", k=length))
