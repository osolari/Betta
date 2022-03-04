import numpy

from cmc import ConditionalMonteCarlo


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
