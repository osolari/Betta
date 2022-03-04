import numpy
import numpy as np
from scipy.stats import pareto
from scipy.stats._distn_infrastructure import rv_continuous

from typing import List


class ConditionalMonteCarlo:
    def __init__(self, random_variables: List[rv_continuous], random_state: int = None):
        """
        To perform Conditional Monte-Carlo estimation of the large-deviance probability
        :param random_variables: a list of instances of rv_continuous. each instance
                must have .pdf(), .cdf(), .rv() methods implemented.
        :param random_state: for reproducible results
        """

        self.random_variables = random_variables
        self.random_state = random_state

    @classmethod
    def construct_pareto_distributions(
        cls,
        shapes: List[float] = None,
        scales: List[float] = None,
        num_rvs: int = 10,
    ):
        """
        A class method for constructing an instance with num_rvs pareto
        distributed random variables with shapes and scales.
        :param shapes: list of shapes of pareto random variables
        :param scales: list of scales of pareto random variables
        :param num_rvs: number of random variables
        :return:
        """

        if (shapes is None and scales is None) and (num_rvs is None):
            raise ValueError(f"Either shapes and scales or num_rvs must be specified")

        if shapes is None:
            shapes = np.linspace(1, 2, num_rvs, endpoint=False)
        if scales is None:
            scales = np.repeat(1, num_rvs)

        random_variables = [
            pareto(b=shape, scale=scale) for shape, scale in zip(shapes, scales)
        ]

        return cls(random_variables=random_variables)

    def sample(self, num_samples: int = 100) -> numpy.ndarray:
        """
        Draws num_samples samples from random_variables
        :param num_samples: number of samples per random variable
        :return: a len(self.random_variables) x num_samples sample matrix
        """

        return numpy.array([rv.rvs(num_samples) for rv in self.random_variables])

    def conditional_monte_carlo(self, sample_matrix: numpy.ndarray, c: int):
        """
        Performs Conditional Monte-Carlo estimation of large deviation probability
        greater than c
        :param sample_matrix: len(self.random_variables) x num_samples sample matrix
        :param c: int
        :return:
        """

        num_rvs = sample_matrix.shape[0]
        num_samples = sample_matrix.shape[1]

        if num_rvs != len(self.random_variables):
            raise ValueError(
                f"Expecting samples from {len(self.random_variables)},"
                f" received {num_rvs}"
            )

        Z = []
        for i, rv in enumerate(self.random_variables):
            ix_i = list(range(i)) + list(range(i + 1, num_samples))

            z_i = rv.sf(
                np.max(
                    np.concatenate(
                        (
                            c - np.sum(sample_matrix[ix_i, :], 0),
                            np.max(sample_matrix[ix_i, :], 0),
                        )
                    ),
                    0,
                )
            )

            Z.append(z_i)

        z = numpy.concatenate(Z).sum(axis=0)
        z_bar = z.mean()

        return z, z_bar


class CMC:
    def __init__(self, shape=None, scale=None, c=100, N=10, n=100, K=50):

        if shape is None:
            self.shape = np.linspace(1, 2, N, endpoint=False)
        else:
            self.shape = shape

        if scale is None:
            self.scale = np.repeat(1, len(shape))
        else:
            self.scale = scale

        assert isinstance(n, int), "n must be an integer."

        self.c = c
        self.N = len(self.shape)
        self.n = n
        self.K = K

    def take_sample(self):

        self.sample = np.asmatrix(
            [
                np.random.pareto(self.shape[i], self.n) * self.scale[i]
                for i in range(len(self.shape))
            ]
        )

    def cmc(self):

        z = 0

        def compZi(i, sample, c, shapei, scalei):

            ix = range(i) + range(i + 1, sample.shape[0])
            return pareto.sf(
                np.max(
                    np.concatenate(
                        (c - np.sum(sample[ix, :], 0), np.max(sample[ix, :], 0))
                    ),
                    0,
                ),
                shapei,
                scale=scalei,
            )

        self.Z = np.sum(
            np.concatenate(
                [
                    compZi(i, self.sample, self.c, self.shape[i], self.scale[i])
                    for i in range(self.sample.shape[0])
                ]
            ),
            0,
        )
        self.Zbar = np.mean(self.Z)

    def mc(self):

        self.muHat = (np.sum(self.sample, 0) > self.c).astype(int)
        self.muBar = np.mean(self.muHat)

    def compRelEff(self, mu, kappa=0.01):
        def mccmc(obj):
            obj.take_sample()
            obj.cmc()
            obj.mc()

            return (obj.muBar, obj.Zbar)

        self.muHatZbars = np.asmatrix([mccmc(self) for i in range(self.K)])

        self.relEff = np.log(
            np.mean(np.abs(self.muHatZbars[:, 1] - mu) > kappa * mu)
        ) - np.log(np.mean(np.abs(self.muHatZbars[:, 0] - mu) > kappa * mu))
