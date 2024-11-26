import numpy as np
from sigconfide.estimates.standard import findSigExposures
from sigconfide.decompose.qp import decomposeQP
from sigconfide.utils.utils import is_wholenumber
from sigconfide.modelselection.backward import compute_p_value, bootstraped_patient

# Funkcje compute_p_value i bootstraped_patient pozostają bez zmian

import numpy as np
from sigconfide.estimates.standard import findSigExposures
from sigconfide.decompose.qp import decomposeQP
from sigconfide.utils.utils import is_wholenumber
from sigconfide.modelselection.backward import compute_p_value, bootstraped_patient

# Funkcje compute_p_value i bootstraped_patient pozostają bez zmian

def hybrid_selection(
    m, P, R, threshold, mutation_count, significance_level, decomposition_method=decomposeQP
):
    """
    Perform a hybrid selection process: start with backward elimination and follow up with forward selection.
    
    This function first performs backward elimination to remove non-significant signatures, then re-evaluates 
    removed signatures to add back any that significantly contribute to the mutation profile.

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

    :returns: A tuple containing the indices of the significant signatures and the exposures and errors from decomposing the original mutation profile.
    :rtype: tuple(numpy.ndarray, tuple(numpy.ndarray, numpy.ndarray))
    """
    # Step 1: Backward Elimination
    best_columns = np.arange(P.shape[1])
    P_temp = P
    M = bootstraped_patient(m, mutation_count, R)

    removed_columns = []
    _, errors = findSigExposures(m.reshape(-1, 1), P[:, best_columns], decomposition_method=decomposition_method)
    current_error = errors[0]
    bs = []
    while True:
        changed = False
        exposures, errors = findSigExposures(M, P[:, best_columns], decomposition_method=decomposition_method)
        p_values = compute_p_value(exposures, threshold=threshold)

        max_p_value = p_values.max()
        if max_p_value > significance_level:
            indices_with_max = np.where(p_values == max_p_value)[0]
            max_p_var = np.random.choice(indices_with_max)
            candidate_columns = np.delete(best_columns, max_p_var)
            _, new_error = findSigExposures(m.reshape(-1, 1), P[:, candidate_columns], decomposition_method=decomposition_method)

            if new_error[0] - current_error <= 0.01:
                bs.append(current_error)
                current_error = new_error[0]
                removed_columns.append(best_columns[max_p_var])
                best_columns = candidate_columns
                changed = True

        if not changed:
            break


    return (
        best_columns,
        findSigExposures(m.reshape(-1, 1), P[:, best_columns], decomposition_method=decomposition_method),
    )