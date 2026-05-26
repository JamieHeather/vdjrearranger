import pytest
from vdjrearranger.chemistry import get_chemistry, Chemistry


def test_get_chemistry_valid():
    chem = get_chemistry("SC5P-R2")
    assert isinstance(chem, Chemistry)
    assert chem.barcode_length == 16
    assert chem.umi_length == 10
    assert chem.r2_length == 91


def test_get_chemistry_invalid():
    with pytest.raises(ValueError, match="Unknown chemistry"):
        get_chemistry("NONEXISTENT-CHEMISTRY")
