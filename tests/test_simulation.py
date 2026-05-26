import pandas as pd
from vdjrearranger.simulation import Simulator, upload_clonotypes


def test_upload_clonotypes(tmp_path):
    mock_tsv = tmp_path / "mock_clonotypes.tsv"
    df = pd.DataFrame({
        "clono_name": ["cell_A", "cell_B"],
        "chain1_nt": ["ATGC", "CGTA"],
        "chain2_nt": ["TTTT", "GGGG"]
    })
    df.to_csv(mock_tsv, sep="\t", index=False)

    clonotypes = upload_clonotypes(
        tsv_path=str(mock_tsv),
        clonotype_col="clono_name",
        chain1_col="chain1_nt",
        chain2_col="chain2_nt"
    )

    assert len(clonotypes) == 2
    assert clonotypes[0].name == "cell_A_row0"
    assert clonotypes[0].chain1 == "ATGC"
    assert clonotypes[1].chain2 == "GGGG"


def test_trim_primer():
    sim = Simulator(clonotypes=[], chemistry="SC5P-R2")
    sim.primers = {"TRA": "CGTGTACC"}

    raw_seq = "AAAAACGTGTACCGCGCGC"

    trimmed = sim.trim_primer("TRA_nt", raw_seq)
    assert trimmed == "AAAAACGTGTACC"

    untrimmed = sim.trim_primer("TRB_nt", raw_seq)
    assert untrimmed == "AAAAACGTGTACCGCGCGC"


def test_fragment_sequence_bounds():
    # Verify sequence lengths enforce chemistry limits
    sim = Simulator(clonotypes=[], chemistry="SC5P-R2", coverage_mode="random")
    long_seq = "A" * 500
    fragment = sim.fragment_sequence(long_seq)

    # The SC5P-R2 chemistry limits R2 to 91 bases
    assert len(fragment) == 91

    # Extremely short sequences shouldn't break the clamping
    short_seq = "A" * 50
    short_fragment = sim.fragment_sequence(short_seq)
    assert len(short_fragment) == 50


def test_tiling_starts_bidirectional():
    sim = Simulator(clonotypes=[], chemistry="SC5P-R2")

    # 200bp sequence, 91bp read length.
    # Step = 91 // 2 = 45. Max start = 200 - 91 = 109.
    starts = sim.get_tiling_starts(seq_len=200, read_len=91, jitter=0)

    # Must include perfect 5' edge
    assert 0 in starts

    # Must include perfect 3' edge
    assert 109 in starts

    # Should contain bidirectional intermediate steps
    assert 45 in starts  # First forward step
    assert 64 in starts  # First backward step (109 - 45)


def test_tiling_starts_with_jitter():
    sim = Simulator(clonotypes=[], chemistry="SC5P-R2")

    jitter_starts = sim.get_tiling_starts(seq_len=200, read_len=91, jitter=10)

    # Forward pass should start offset by the jitter
    assert 10 in jitter_starts

    # Backward pass should start offset by the jitter from the end
    assert (109 - 10) in jitter_starts
