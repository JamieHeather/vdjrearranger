import pytest
from vdjrearranger.distributions import CountSampler


def test_fixed_sampler():
    sampler = CountSampler("fixed", 5)
    assert sampler.sample() == 5


def test_normal_sampler_bounds():
    # Mean of -10 should still be clamped to 1
    sampler = CountSampler("normal", -10, 1)
    assert sampler.sample() == 1


def test_poisson_sampler():
    sampler = CountSampler("poisson", 100)
    val = sampler.sample()
    assert isinstance(val, int)
    assert val >= 1


def test_invalid_mode():
    sampler = CountSampler("unknown_mode", 10)
    with pytest.raises(ValueError, match="Unsupported mode"):
        sampler.sample()
