import quadprog
import numpy as np

def decomposeQP(m, P):
    # N: how many signatures are selected
    N = P.shape[1]
    # G: matrix appearing in the quadratic programming objective function
    G = np.dot(P.T, P).astype(float)
    # C: matrix constraints under which we want to minimize the quadratic programming objective function.
    C = np.column_stack([np.ones(N), np.eye(N)]).astype(float)
    # b: vector containing the values of b_0.
    b = np.array([1] + [0]*N).astype(float)
    # d: vector appearing in the quadratic programming objective function
    d = np.dot(m.T, P).astype(float)

    # Solve quadratic programming problem
    out = quadprog.solve_qp(G, d, C, b, meq=1)

    # Some exposure values are negative, but very close to 0
    # Change these negative values to zero and renormalize
    exposures = out[0]
    exposures[exposures < 0] = 0
    exposures /= sum(exposures)

    # return the exposures
    return exposures