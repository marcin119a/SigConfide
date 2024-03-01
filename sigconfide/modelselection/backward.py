import numpy as np
from sigconfide.estimates.bootstrap import bootstrapSigExposures
from sigconfide.estimates.standard import findSigExposures

from sigconfide.decompose.qp import decomposeQP
from sigconfide.utils.utils import is_wholenumber


def compute_p_value(exposures, threshold=0.01):
    """
    Calculate the proportion of exposures greater than a specified threshold for each signature.

    This function computes the p-value for each signature by determining the fraction of bootstrap exposures that exceed a given threshold. It's used to assess the significance of each signature's contribution.

    :param exposures: A numpy array of signature exposures obtained from bootstrap resampling.
    :type exposures: numpy.ndarray
    :param threshold: The threshold above which exposures are considered significant. Defaults to 0.01.
    :type threshold: float, optional

    :returns: A numpy array of p-values for each signature, indicating the proportion of non-significant exposures.
    :rtype: numpy.ndarray
    """

    grater_than_threshold = exposures > threshold

    return 1 - grater_than_threshold.sum(axis=1) / grater_than_threshold.shape[1]


def bootstraped_patient(m, mutation_count, R):
    """
    Generate a bootstrap distribution of mutation profiles for a patient/sample.

    This function creates a set of bootstrap replicates of the patient's mutation profile by resampling with replacement. It's used to simulate the variability in the mutation profile and assess the stability of the signature exposures.

    :param m: The observed mutation profile vector for a patient/sample.
    :type m: numpy.ndarray
    :param mutation_count: The total number of mutations to be resampled in each bootstrap replicate. If not provided, it will be calculated from 'm'.
    :type mutation_count: int
    :param R: The number of bootstrap replicates to generate.
    :type R: int

    :raises ValueError: If 'mutation_count' is not specified and 'm' does not contain integer counts.

    :returns: A matrix of bootstrap replicates of the patient's mutation profile.
    :rtype: numpy.ndarray
    """
    K = len(m)

    if mutation_count is None:
        if all(is_wholenumber(val) for val in m):
            mutation_count = int(m.sum())
        else:
            raise ValueError("Please specify the parameter 'mutation_count' in the function call or provide mutation counts in parameter 'm'.")
    m = m / np.sum(m)


    def bootstrap_sample(m, mutation_count, K):
        mutations_sampled = np.random.choice(K, size=mutation_count, p=m)
        return np.bincount(mutations_sampled, minlength=K) / mutation_count

    M = np.column_stack([bootstrap_sample(m, mutation_count, K) for _ in range(R)])

    return M


def backward_elimination(
    m, P, R, threshold, mutation_count, significance_level, decomposition_method=decomposeQP
):
    """
    Perform backward elimination to identify the most significant signatures contributing to a patient's mutation profile.

    Starting with all signatures, this function iteratively removes the least significant signature until all remaining signatures have a p-value below a specified significance level. It uses bootstrap resampling to assess the significance of each signature's contribution.

    :param m: The observed mutation profile vector for a patient/sample.
    :type m: numpy.ndarray
    :param P: The signature profile matrix.
    :type P: numpy.ndarray
    :param R: The number of bootstrap replicates used for significance testing.
    :type R: int
    :param threshold: The threshold used to determine if a signature's exposure is significant in the bootstrap samples.
    :type threshold: float
    :param mutation_count: The total number of mutations in the patient's profile. Used for bootstrap resampling.
    :type mutation_count: int
    :param significance_level: The p-value threshold below which a signature is considered significant.
    :type significance_level: float
    :param decomposition_method: The method used to decompose the mutation profile into signature exposures. Defaults to 'decomposeQP'.
    :type decomposition_method: function, optional

    :returns: A tuple containing the indices of the significant signatures, the exposures and errors from the final bootstrap resampling, and the exposures and errors from decomposing the original mutation profile.
    :rtype: tuple(numpy.ndarray, tuple(numpy.ndarray, numpy.ndarray), tuple(numpy.ndarray, numpy.ndarray))
    """

    best_columns = np.arange(P.shape[1])
    P_temp = P

    M = bootstraped_patient(m, mutation_count, R)

    while True:
        changed = False

        exposures, errors = findSigExposures(
            M, P_temp, decomposition_method=decomposition_method
        )
        p_values = compute_p_value(exposures, threshold=threshold)

        max_p_value = p_values.max()
        if max_p_value > significance_level:
            indices_with_max = np.where(p_values == max_p_value)[0]
            max_p_var = np.random.choice(indices_with_max)
            best_columns = np.delete(best_columns, max_p_var)
            P_temp = P[:, best_columns]

            changed = True

        if not changed:
            break

    return (
        best_columns,
        bootstrapSigExposures(
            m,
            P_temp,
            mutation_count=mutation_count,
            R=R,
            decomposition_method=decomposition_method,
        ),
        findSigExposures(
            m.reshape(-1, 1), P_temp, decomposition_method=decomposition_method
        ),
    )


