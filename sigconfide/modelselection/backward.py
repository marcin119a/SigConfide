import numpy as np
from estimages.bootstrap import bootstrapSigExposures
from estimates.standard import findSigExposures

from decompose.qp import decomposeQP
from utils.utils import is_wholenumber
def compute_p_value(exposures, threshold=0.01):
    grater_than_threshold = exposures > threshold

    return 1 - grater_than_threshold.sum(axis=1) / grater_than_threshold.shape[1]


def bootstraped_patient(m, mutation_count, R):
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
    best_columns = np.arange(P.shape[1])
    P_temp = P
    model_error, _ = findSigExposures(
        m.reshape(-1, 1), P_temp, decomposition_method
    )
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


