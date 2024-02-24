import numpy as np
from sigconfide.decompose.qp import decomposeQP
from sigconfide.utils.utils import FrobeniusNorm
def findSigExposures(M, P, decomposition_method=decomposeQP):
    """
     Find signature exposures for tumor profiles using specified decomposition method.

     This function allows obtaining the optimal solution by specifying quadratic programming to solve the optimization problem.

     Parameters:
         M (numpy.ndarray): Observed tumor profile matrix for all patients/samples.
             It should have a shape of (96, G), where G is the number of patients.
             Each column can represent mutation counts or mutation probabilities, and
             each column will be normalized to sum up to 1.
         P (numpy.ndarray): Signature profile matrix with a shape of (96, N),
             where N is the number of signatures (e.g., COSMIC: N=30).
         decomposition_method (function, optional): The method selected to get the
             optimal solution. It should be a function. Default is 'decomposeQP'.

     Returns:
         tuple: A tuple containing two numpy arrays.
             - exposures (numpy.ndarray): Matrix of signature exposures per sample/patient (column).
             - errors (numpy.ndarray): Estimation error for each sample/patient (Frobenius norm).

     Raises:
         ValueError: If 'M' and 'P' do not have the same number of rows (mutations types),
             or if 'P' has less than 2 columns, or if 'decomposition_method' is not a function.

     Examples:
         E1 = findSigExposures(tumorBRCA, signaturesCOSMIC, decomposeQP)
         sigsBRCA = [1, 2, 3, 5, 6, 8, 13, 17, 18, 20, 26, 30]
         E2 = findSigExposures(tumorBRCA, signaturesCOSMIC[:, sigsBRCA], decomposeQP)
         E3 = findSigExposures(np.round(tumorBRCA * 10000), signaturesCOSMIC, decomposeQP)
     """
    # Process and check function parameters
    # M, P
    if M.shape[0] != P.shape[0]:
        raise ValueError("Matrices 'M' and 'P' must have the same number of rows (mutations types).")

    if P.shape[1] == 1:
        raise ValueError("Matrices 'P' must have at least 2 columns (signatures).")

    # decomposition.method
    if not callable(decomposition_method):
        raise ValueError("Parameter 'decomposition_method' must be a function.")

    # Normalize M by column (just in case it is not normalized)
    M = M / M.sum(axis=0)

    # Find solutions
    # Matrix of signature exposures per sample/patient (column)
    exposures = np.apply_along_axis(decomposition_method, 0, M, P)

    # Compute estimation error for each sample/patient (Frobenius norm)
    errors = np.vectorize(lambda i: FrobeniusNorm(M[:, i], P, exposures[:, i]))(range(M.shape[1]))

    return exposures, errors
