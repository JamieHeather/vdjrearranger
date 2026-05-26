import gzip
from pathlib import Path
from vdjrearranger.models import ReadPair


class FastqWriter:
    """
    Handles sequential writing of paired-end sequences to standard gzipped FASTQ files.
    """

    def __init__(self, outdir: Path, sample_name: str, lane: int):
        """
        Initializes output file streams.

        :param outdir: str, directory path where the FASTQ files will be created.
        :param sample_name: str, prefix used for the FASTQ filenames.
        :param lane: int, Illumina lane integer to insert into the filename formatting.
        """
        r1_name = f"{sample_name}_S1_L{lane:03d}_R1_001.fastq.gz"
        r2_name = f"{sample_name}_S1_L{lane:03d}_R2_001.fastq.gz"

        self.r1_path = outdir / r1_name
        self.r2_path = outdir / r2_name

        self.r1 = gzip.open(self.r1_path, "wt")
        self.r2 = gzip.open(self.r2_path, "wt")

    def write_pair(self, pair: ReadPair):
        """
        Formats and writes a single ReadPair object to the open FASTQ handles.

        :param pair: A populated ReadPair dataclass instance.
        """
        self.r1.write(
            f"@{pair.read_id} 1:N:0:GTCCTAAC\n"
            f"{pair.r1_seq}\n"
            f"+\n"
            f"{pair.r1_qual}\n"
        )

        self.r2.write(
            f"@{pair.read_id} 2:N:0:GTCCTAAC\n"
            f"{pair.r2_seq}\n"
            f"+\n"
            f"{pair.r2_qual}\n"
        )

    def close(self):
        """
        Safely flushes and closes the underlying gzip file handles.
        """
        self.r1.close()
        self.r2.close()
