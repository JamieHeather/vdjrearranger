import numpy as np


class CountSampler:
    """
    Provides a standardized interface for sampling counts from various statistical distributions.
    Used to model cell population sizes, UMI capture rates, and PCR duplicate read counts.
    """

    def __init__(self, mode: str, value: int, variance: float = 0):
        """
        Initialises the sampler with specific distributional parameters.

        :param mode: str, type of distribution ('fixed', 'poisson', 'normal', 'powerlaw').
        :param value: int, the primary abundance parameter (mean, fixed value, or scale factor).
        :param variance: int, secondary spread parameter (standard deviation or alpha for powerlaw).
        """
        self.mode = mode
        self.value = value
        self.variance = variance

    def sample(self) -> int:
        """
        Draws a single integer sample from the configured distribution.
        Values are strictly clamped to a minimum of 1.

        :return: A sampled integer count.
        :raises ValueError: If an unsupported mode string is provided during instantiation.
        """
        if self.mode == "fixed":
            return max(1, int(self.value))

        if self.mode == "poisson":
            return max(1, int(np.random.poisson(self.value)))

        if self.mode == "normal":
            return max(
                1,
                int(np.random.normal(self.value, self.variance))
            )

        if self.mode == "powerlaw":
            alpha = self.variance if self.variance > 0 else 1.5
            val = np.random.pareto(alpha) * self.value
            return max(1, int(round(val)))

        raise ValueError(f"Unsupported mode: {self.mode}")
