import numpy as np
from sigconfide.utils.utils import is_wholenumber, FrobeniusNorm
from sigconfide.decompose.qp import decomposeQP

def bootstrapSigExposures(m, P, R, mutation_count=None, decomposition_method=decomposeQP):
    """
    Obtain the bootstrap distribution of signature exposures for a tumor sample.

    This function allows obtaining the bootstrap distribution of the signature exposures for
    a specific tumor sample using a specified decomposition method.

    Parameters:
        m (numpy.ndarray): Observed tumor profile vector for a patient/sample.
            It should have a shape of (96, 1) and can represent mutation counts or mutation probabilities.
        P (numpy.ndarray): Signature profile matrix with a shape of (96, N),
            where N is the number of signatures (e.g., COSMIC: N=30).
        R (int): The number of bootstrap replicates.
        mutation_count (int, optional): If 'm' is a vector of counts, then 'mutation_count' equals
            the summation of all the counts. If 'm' is probabilities, 'mutation_count' must be specified.
        decomposition_method (function, optional): The method selected to get the optimal solution.
            It should be a function. Default is 'decomposeQP'.

    Returns:
        tuple: A tuple containing two numpy arrays.
            - exposures (numpy.ndarray): Matrix of signature exposures for each bootstrap replicate (column).
            - errors (numpy.ndarray): Estimation error for each bootstrap replicate (Frobenius norm).

    Raises:
        ValueError: If the length of vector 'm' and the number of rows of matrix 'P' do not match,
            if 'P' has less than 2 columns, if 'mutation_count' is not specified and 'm' does not contain counts.

    Examples:
        bootstrapSigExposures(tumorBRCA[:, 1], signaturesCOSMIC, 100, 2000, decomposeQP)
        sigsBRCA = [1, 2, 3, 5, 6, 8, 13, 17, 18, 20, 26, 30]
        bootstrapSigExposures(tumorBRCA[:, 1], signaturesCOSMIC[:, sigsBRCA], 10, 1000, decomposeQP)
    """

    if len(m) != P.shape[0]:
        raise ValueError("Length of vector 'm' and number of rows of matrix 'P' must be the same.")
    if m.shape[0] != P.shape[0]:
        raise ValueError("Elements of vector 'm' and rows of matrix 'P' must have the same names (mutations types).")
    #if P.shape[1] == 1:
    #    raise ValueError("Matrices 'P' must have at least 2 columns (signatures).")

    # If 'mutation_count' is not specified, 'm' has to contain counts
    if mutation_count is None:
        if all(is_wholenumber(val) for val in m):
            mutation_count = int(m.sum())
        else:
            raise ValueError("Please specify the parameter 'mutation_count' in the function call or provide mutation counts in parameter 'm'.")

    # Normalize m to be a vector of probabilities.
    m = m / np.sum(m)

    # Find optimal solutions using provided decomposition method for each bootstrap replicate
    # Matrix of signature exposures per replicate (column)
    K = len(m)  # number of mutation types

    def bootstrap_sample(m, mutation_count, K):
        mutations_sampled = np.random.choice(K, size=mutation_count, p=m)
        return np.bincount(mutations_sampled, minlength=K) / mutation_count

    exposures = np.column_stack([
        decomposition_method(bootstrap_sample(m, mutation_count, K), P) for _ in range(R)
    ])
    exposures = exposures / np.sum(exposures, axis=0)  # Normalize exposures

    # Compute estimation error for each replicate/trial (Frobenius norm)
    # G x R
    errors = np.vectorize(lambda i: FrobeniusNorm(m, P, exposures[:, i]))(range(exposures.shape[1]))

    return exposures, errors
