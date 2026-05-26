from vdjrearranger.barcode import random_umi, sample_barcodes


def test_random_umi_default_length():
    umi = random_umi()
    assert len(umi) == 10
    assert all(base in "ACGT" for base in umi)


def test_random_umi_custom_length():
    umi = random_umi(15)
    assert len(umi) == 15


def test_sample_barcodes():
    mock_whitelist = ["AAAA", "CCCC", "GGGG", "TTTT"]
    sampled = sample_barcodes(2, mock_whitelist)

    assert len(sampled) == 2
    # Ensure they are unique
    assert len(set(sampled)) == 2
    # Ensure they actually came from the whitelist
    assert all(b in mock_whitelist for b in sampled)
