from dataclasses import dataclass


@dataclass(frozen=True)
class Chemistry:
    """
    Defines the structural properties of a single-cell sequencing library format.

    :param name: Identifier for the chemistry type (e.g., 'SC5P-R2').
    :param barcode_length: Expected length of the cell barcode in nucleotides.
    :param umi_length: Expected length of the Unique Molecular Identifier.
    :param r1_tail: Optional static sequence appended to Read 1.
    :param r1_length: Total simulated length of Read 1.
    :param r2_length: Total simulated length of Read 2 (the transcript).
    """
    name: str
    barcode_length: int
    umi_length: int
    r1_tail: str
    r1_length: int
    r2_length: int


CHEMISTRIES = {
    "SC5P-R2": Chemistry(
        name="SC5P-R2",
        barcode_length=16,
        umi_length=10,
        r1_tail="",
        r1_length=26,
        r2_length=91,
    )
}


def get_chemistry(name: str) -> Chemistry:
    """
    Retrieves the configuration for a given library chemistry.

    :param name: str, ame of the chemistry to look up.
    :return: A populated Chemistry dataclass instance.
    :raises ValueError: If the requested chemistry is not found in the CHEMISTRIES mapping.
    """
    if name not in CHEMISTRIES:
        raise ValueError(
            f"Unknown chemistry '{name}'. "
            f"Available: {list(CHEMISTRIES.keys())}"
        )
    return CHEMISTRIES[name]
