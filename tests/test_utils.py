import unittest
from sigconfide.utils.utils import *
import os
class TestUtils(unittest.TestCase):
    def test_detect_format(self):
        """Test the detection of data formats."""
        csv_line = 'C>A,ACA,66,83,88,78,90,46,73'
        tsv_line = 'AC>AA\t66\t83\t88\t78\t90\t46\t73'
        mutated_tsv_line = 'A[C>A]A\t66\t83\t88\t78\t90\t46\t73'  # Example same as TSV for illustration
        unknown_line = 'C>A ACA 66 83 88 78 90 46 73'

        # Test CSV format detection
        self.assertEqual(detect_format(csv_line), 'CSV Format', "Should detect CSV Format")

        # Test TSV format detection
        self.assertEqual(detect_format(tsv_line), 'TSV Format', "Should detect TSV Format")

        # Test Mutated TSV format detection
        self.assertEqual(detect_format(mutated_tsv_line), 'Mutated TSV Format', "Should detect Mutated TSV Format")

        # Test unknown format detection
        self.assertEqual(detect_format(unknown_line), 'Unknown Format', "Should detect Unknown Format")

    def test_detect_format(self):
        """Test the detection of data formats."""
        with open('data/format_2.dat', 'r') as file:
            csv_line = file.readline().strip()

        with open('data/format_1.dat', 'r') as file:
            tsv_line = file.readline().strip()

        # Assume detect_format function exists and performs format detection
        # Test CSV format detection
        self.assertEqual(detect_format(csv_line), 'CSV Format', "Should detect CSV Format")

        # Test TSV format detection
        self.assertEqual(detect_format(tsv_line), 'TSV Format', "Should detect TSV Format")

    def test_frobenius_norm(self):
        M = np.array([[1, 2], [4, 5]])
        P = np.array([[1, 1], [2, 8]])
        E = np.array([[2, 0], [1, 1]])

        result = FrobeniusNorm(M, P, E)

        # Define the expected Frobenius norm value
        expected_result = 8.83176

        # Check if the result matches the expected result
        self.assertAlmostEqual(result, expected_result, places=5)
    def test_valid_data_files(self):
        """Test loading and processing with valid data files."""
        profile, signatures = load_and_process_data(0, "data/tumorBRCA.csv", "data/signaturesCOSMIC.csv")
        self.assertIsNotNone(profile, "Profile should not be None")
        self.assertIsNotNone(signatures, "Signatures should not be None")

    def test_nonexistent_files(self):
        """Test handling of non-existent files."""
        with self.assertRaises(FileNotFoundError):
            load_and_process_data(0, "nonexistent.csv", "nonexistent.csv")

    def test_empty_files(self):
        """Test handling of empty files."""
        # Create empty file
        open('data/empty.csv', 'w').close()
        with self.assertRaises(ValueError):
            load_and_process_data(0, "data/empty.csv", "data/empty.csv")
        # Cleanup
        os.remove('data/empty.csv')


    def test_is_wholenumber(self):
        # Test with whole numbers
        self.assertTrue(is_wholenumber(2), "2 should be identified as a whole number")
        self.assertTrue(is_wholenumber(0), "0 should be identified as a whole number")
        self.assertTrue(is_wholenumber(-3), "-3 should be identified as a whole number")

        # Test with a number very close to whole number (within tolerance)
        self.assertTrue(is_wholenumber(2.000000000000001),
                        "2.000000000000001 should be identified as a whole number")

        # Test with non-whole numbers
        self.assertFalse(is_wholenumber(3.5), "3.5 should not be identified as a whole number")
        self.assertFalse(is_wholenumber(-1.2), "-1.2 should not be identified as a whole number")

        # Test with numbers close to whole numbers within the default tolerance
        self.assertTrue(is_wholenumber(2.0000000000000001),
                        "2.0000000000000001 should be identified as a whole number within default tolerance")
        self.assertTrue(is_wholenumber(-3.0000000000000001),
                        "-3.0000000000000001 should be identified as a whole number within default tolerance")
