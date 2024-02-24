import unittest
from sigconfide.estimates.bootstrap import bootstrapSigExposures
from sigconfide.estimates.crossvalidation import crossValidationSigExposures
from sigconfide.estimates.standard import findSigExposures
from sigconfide.utils.utils import load_and_process_data
import numpy as np


class TestEstimateExposures(unittest.TestCase):
    def test_findSigExposures(self):
        M = np.array([[0.5, 0.3, 0.2], [0.9, 0.05, 0.05], [0.7, 0.1, 0.2]])
        P = np.array([[0.2, 0.3, 0.5], [0.1, 0.4, 0.5], [0.3, 0.1, 0.6]])

        exposures, errors = findSigExposures(M, P)
        #data obtaining from R code
        expected_exposures = np.array([
            [0.2007233, 0.4317862, 0.6437908],
            [0.4755877, 0.2883263, 0.0000000],
            [0.3236890, 0.2798875, 0.3562092]
        ])
        expected_errors = np.array([
            0.1245907, 0.4137049, 0.1939068]
        )

        np.testing.assert_array_almost_equal(exposures, expected_exposures, decimal=7)
        np.testing.assert_array_almost_equal(errors, expected_errors, decimal=7)

    def test_findSigExposuresReal(self):
        profile, signatures = load_and_process_data(None,
                                                         'data/tumorBRCA.csv',
                                                         'data/signaturesCOSMIC.csv')

        exposures, errors = findSigExposures(profile, signatures)
        expected_exposures = np.genfromtxt('data/R_exposures.csv', delimiter=',', skip_header=1)
        expected_exposures = np.delete(expected_exposures, 0, axis=1).squeeze()
        expected_errors = np.genfromtxt('data/R_errors.csv', delimiter=',', skip_header=1)
        expected_errors = np.delete(expected_errors, 0, axis=1).squeeze()

        np.testing.assert_array_almost_equal(exposures, expected_exposures, decimal=7)
        np.testing.assert_array_almost_equal(errors, expected_errors, decimal=7)
class TestBootstrapSigExposures(unittest.TestCase):
    def test_bootstrap_sample(self):
        m = np.array([0.5, 0.3, 0.2])
        P = np.array([[0.2, 0.3, 0.5], [0.1, 0.4, 0.5], [0.3, 0.1, 0.6]])

        expected_exposures = np.array(
            [[0.1025316, 0.2443038, 0.0481013],
            [0.664557, 0.4797468, 0.7065823],
            [0.2329114, 0.2759494, 0.2453165]]
        )

        np.random.seed(42)

        exposures, errors = bootstrapSigExposures(m, mutation_count=100, R=3, P=P)

        np.testing.assert_array_almost_equal(exposures, expected_exposures, decimal=7)


class TestCrossValidationSigExposures(unittest.TestCase):

    def test_crossvalidation(self):
        m = np.array([0.5, 0.3, 0.2])
        P = np.array([[0.2, 0.3, 0.5], [0.1, 0.4, 0.5], [0.3, 0.1, 0.6]])
        np.random.seed(42)
        expected_exposures = np.array([
            [0.0886076, 0.6764706, 0.],
            [0.5594937, 0., 0.9583333],
            [0.3518987, 0.3235294, 0.0416667]]
        )

        exposures, errors = crossValidationSigExposures(m, P, fold_size=1)
        np.testing.assert_array_almost_equal(exposures, expected_exposures, decimal=7)

if __name__ == '__main__':
    unittest.main()
