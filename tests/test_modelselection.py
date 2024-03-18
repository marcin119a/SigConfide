import numpy as np
import unittest

from sigconfide.modelselection.backward import compute_p_value
from sigconfide.modelselection.backward import bootstraped_patient
class TestComputePValue(unittest.TestCase):

    def test_compute_p_value(self):
        # Create a mock exposures array
        exposures = np.array([
            [0.02, 0.03, 0.00, 0.01],  # 2 out of 4 exposures are greater than threshold
            [0.05, 0.02, 0.03, 0.04],  # All exposures are greater than threshold
            [0.00, 0.00, 0.00, 0.00]   # No exposures are greater than threshold
        ])

        # Expected p-values
        expected = np.array([
            0.5,  # 50% of the exposures are not significant
            0.0,  # 0% of the exposures are not significant
            1.0   # 100% of the exposures are not significant
        ])

        # Compute the p-values
        actual = compute_p_value(exposures)

        # Assert that the actual p-values match the expected p-values
        np.testing.assert_array_almost_equal(actual, expected, decimal=2)


class TestBootstrapedPatient(unittest.TestCase):

    def test_bootstraped_patient_valid_input(self):
        m = np.array([10, 20, 30])
        mutation_count = 60  # Explicit mutation count
        R = 10  # Number of replicates
        bootstrap_replicates = bootstraped_patient(m, mutation_count, R)

        self.assertEqual(bootstrap_replicates.shape, (len(m), R))  # Check shape of the output

    def test_bootstraped_patient_no_mutation_count(self):
        m = np.array([10, 20, 30])
        R = 10
        bootstrap_replicates = bootstraped_patient(m, None, R)

        self.assertEqual(bootstrap_replicates.shape, (len(m), R))  # Should compute mutation_count automatically

    def test_bootstraped_patient_non_integer_values(self):
        m = np.array([10.5, 20.2, 30.3])  # Non-integer values
        R = 10

        with self.assertRaises(ValueError):
            bootstraped_patient(m, None, R)  # Should raise ValueError

