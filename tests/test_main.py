from sigconfide.modelselection.analyzer import cosmic_fit

import unittest
import os
import shutil
current_dir = os.path.dirname(os.path.abspath(__file__))

def remove_folder(folder_path):
    try:
        shutil.rmtree(folder_path)
    except OSError as error:
        print({error})
class TestCosmicFit(unittest.TestCase):

    def test_cosmic_fit_generates_csv(self):
        output_dir = 'output'
        expected_output_filename = 'Assignment_Solution_Activities.csv'
        expected_output_path = os.path.join(output_dir, expected_output_filename)

        if os.path.exists(expected_output_path):
            os.remove(expected_output_path)

        cosmic_fit(os.path.join(current_dir, 'data', 'reduced_data.dat'),  output_dir, cosmic_version=2.0)

        self.assertTrue(os.path.exists(expected_output_path), "The CSV file was not generated.")
        remove_folder(output_dir)

