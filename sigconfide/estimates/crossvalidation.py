import numpy as np
from sigconfide.utils.utils import FrobeniusNorm
from sigconfide.decompose.qp import decomposeQP


def crossValidationSigExposures(m, P, fold_size, shuffle=True, decomposition_method=decomposeQP):
    """
    Perform cross-validation to estimate signature exposures for a tumor sample.

    This function performs cross-validation to estimate signature exposures for a specific tumor sample
    using a specified decomposition method.

    Parameters:
        m (numpy.ndarray): Observed tumor profile vector for a patient/sample.
            It should have a shape of (n, 1), where n is the number of mutations.
        P (numpy.ndarray): Signature profile matrix with a shape of (n, N),
            where N is the number of signatures.
        fold_size (int): The number of cross-validation size.
        shuffle (bool): Change the order of mutations
        decomposition_method (function, optional): The method selected to get the optimal solution.
            It should be a function. Default is 'decomposeQP'.

    Returns:
        tuple: A tuple containing two numpy arrays.
            - fold_exposures (numpy.ndarray): Matrix of signature exposures for each cross-validation fold (column).
            - errors (numpy.ndarray): Estimation error for each fold (Frobenius norm).

    Raises:
        ValueError: If the length of vector 'm' and the number of rows of matrix 'P' do not match,
            if 'P' has less than 2 columns.

    Examples:
        num_folds = 5
        sigsBRCA = [1, 2, 3, 5, 6, 8, 13, 17, 18, 20, 26, 30]
        fold_exposures, errors = crossValidationSigExposures(tumorBRCA[:, 1], signaturesCOSMIC[:, sigsBRCA], num_folds=5, decomposeQP)
    """
    # Process and check function parameters
    P = np.array(P)

    if len(m) != P.shape[0]:
        raise ValueError("Length of vector 'm' and number of rows of matrix 'P' must be the same.")

    if P.shape[1] == 1:
        raise ValueError("Matrices 'P' must have at least 2 columns (signatures).")

    m = m / np.sum(m)

    if shuffle:
        permutation_indices = np.random.permutation(len(m))
        m = m[permutation_indices]
        P = P[permutation_indices,:]

    folds = [m[i:i + fold_size] for i in range(0, (len(m) - (len(m) % fold_size)), fold_size)]

    # Handle the remaining elements that do not fit in full folds
    if len(m) % fold_size != 0:
        last_fold = m[-(len(m) % fold_size):]
        folds.append(last_fold)

    def calculate_fold_exposures(i, num_folds):
        fold = np.concatenate([folds[j] if j != i else [0] * len(folds[i] + 1) for j in range(num_folds)])
        normalized_fold = fold / fold.sum()
        return normalized_fold

    # Perform cross-validation for each replicate
    fold_exposures = np.column_stack([
        decomposition_method(calculate_fold_exposures(i, len(folds)), P)
        for i in range(len(folds))
    ])
    fold_exposures = fold_exposures / np.sum(fold_exposures, axis=0)

    # Compute estimation error for each replicate/trial (Frobenius norm)
    errors = np.vectorize(lambda i: FrobeniusNorm(m, P, fold_exposures[:, i]))(range(fold_exposures.shape[1]))

    return fold_exposures, errors

