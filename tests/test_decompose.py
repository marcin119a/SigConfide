import numpy as np
import unittest
from scipy.optimize import minimize

from sigconfide.decompose.qp import decomposeQP
from sigconfide.utils.utils import load_and_process_data, FrobeniusNorm


def decomposeQPScipy(m, P):
    N = P.shape[1]

    def objective(E):
        return FrobeniusNorm(m, P, E)

    constraints = ({'type': 'eq', 'fun': lambda E: np.array([np.sum(E) - 1])})
    bounds = [(0, None)] * N

    out = minimize(
        objective, x0=np.zeros(N), method='SLSQP', bounds=bounds, constraints=constraints)


    exposures = out.x
    exposures[exposures < 0] = 0
    exposures = exposures / np.sum(exposures)

    return exposures

class TestDecomposeQP(unittest.TestCase):
    def test_decomposeQP(self):
        # Example input data
        m = np.array([0.5, 0.3, 0.2])
        P = np.array([[0.2, 0.3, 0.5], [0.1, 0.4, 0.5], [0.3, 0.1, 0.6]])

        exposuresFast = decomposeQP(m, P)

        exposuresSlow = decomposeQPScipy(m, P)

        np.testing.assert_array_almost_equal(exposuresFast, exposuresSlow, decimal=3)

    def test_real_data_decomposeQP(self):
        first_patient, signaturesCOSMIC = (
            load_and_process_data(patient_index=0,
                                  mutational_profiles='data/tumorBRCA.csv',
                                  predf_mutational_signatures='data/signaturesCOSMIC.csv'))

        exposuresFast = decomposeQP(first_patient, signaturesCOSMIC)
        exposuresSlow = decomposeQPScipy(first_patient, signaturesCOSMIC)

        np.testing.assert_array_almost_equal(exposuresFast, exposuresSlow, decimal=1)


if __name__ == '__main__':
    unittest.main()
