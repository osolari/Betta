import numpy
from cmc import ConditionalMonteCarlo

import logging

logger = logging.getLogger("cmc")


def monte_carlo_estimate(sample_matrix: numpy.ndarray, c: int):
    """

    :param sample_matrix:
    :param c:
    :return:
    """

    mu_hat = (sample_matrix.sum(axis=0) > c).astype(int)
    mu_bar = mu_hat.mean()

    return mu_hat, mu_bar


def compute_relative_efficiency(
    mu, num_samples: int = 100, c: int = 100, k: int = 50, kappa: float = 0.01
):
    """

    :param mu:
    :param num_samples:
    :param c:
    :param k:
    :param kappa:
    :return:
    """

    z_bar = []
    mu_bar = []
    for _ in range(k):
        cmc = ConditionalMonteCarlo.construct_pareto_distributions()
        sample_matrix = cmc.sample(num_samples=num_samples)
        z_bar.append(cmc.conditional_monte_carlo(sample_matrix=sample_matrix, c=c)[1])
        mu_bar.append(monte_carlo_estimate(sample_matrix=sample_matrix, c=c)[1])

    z_bar = numpy.array(z_bar)
    mu_bar = numpy.array(mu_bar)

    relative_efficiency = numpy.log(
        (numpy.abs(z_bar - mu) > kappa * mu).mean()
    ) - numpy.log((numpy.abs(mu_bar - mu) > kappa * mu).mean())

    return relative_efficiency


def compute_mu(shapes, c):
    """

    :param shapes:
    :param c:
    :return:
    """

    cmc = ConditionalMonteCarlo.construct_pareto_distributions(shapes=shapes)
    sample_matrix = cmc.sample(num_samples=1e7)
    _, z_bar = cmc.conditional_monte_carlo(sample_matrix=sample_matrix, c=c)

    return z_bar


def compute_relative_efficiency_curve(
    shapes=numpy.linspace(1, 2, 10, endpoint=False),
    num_samples=numpy.linspace(2, 2000, 50).astype(int),
    kappa=0.01,
    k=20,
    c=100,
):

    """

    :param shapes:
    :param num_samples:
    :param kappa:
    :param k:
    :param c:
    :return:
    """

    logger.info("Computing Mu")
    mu = compute_mu(shapes=shapes, c=c)
    logger.info(f"Estimated mu {mu}")

    relative_efficiencies = numpy.zeros((k, len(num_samples)))
    for i in range(len(num_samples)):

        logger.info(f"Computing relative efficiency for {num_samples[i]} samples.")
        tmp_cmc = ConditionalMonteCarlo.construct_pareto_distributions(shapes=shapes)

        for j in range(k):
            relative_efficiencies[j, i] = compute_relative_efficiency(
                mu, num_samples=num_samples[i], c=c, k=k, kappa=kappa
            )

    return num_samples, relative_efficiencies, shapes, mu, c


def compute_relative_efficiency_variable_shape(kappa: int = 0.005, num_jobs: int = 1):
    """

    :param kappa:
    :param num_jobs:
    :return:
    """
    shape_factors = [numpy.linspace(1, e, 10, endpoint=False) for e in range(2, 6)]

    return [
        compute_relative_efficiency_curve(shape, kappa=kappa) for shape in shape_factors
    ]


def compute_relative_efficiency_variiable_c(kappa=5e-3, num_jobs: int = 1):
    """

    :param kappa:
    :param num_jobs:
    :return:
    """

    c = numpy.linspace(100, 1000, 10, endpoint=True)
    logger.info(f"Computing relative efficiency with bounds {c}")
    return [
        compute_relative_efficiency_curve(
            shapes=numpy.linspace(1, 3, 10, endpoint=False), c=ci, kappa=kappa
        )
        for ci in c
    ]


# def compute_relative_efficiency_variiable_minimum_shape
# K
# def compRelEffVarShapeVarMin(
#     kappa=0.005,
#     c=100,
#     ALPHA=[
#         (np.float(e) / 2 + np.linspace(0, 1, 10, endpoint=False)) for e in range(2, 6)
#     ],
# ):
#
#     return [compRelEffCurve(shape, kappa=kappa) for shape in ALPHA]
