from sigconfide.modelselection.backward import backward_elimination
from sigconfide.utils.utils import load_samples_file, load_signatures_file
from sigconfide.utils import utils
import numpy as np
import os
import sys

module_path = os.path.dirname(utils.__file__)

versions = {
    1.0: f'{module_path}/data/COSMIC_v1_SBS_GRCh37.txt',
    2.0: f'{module_path}/data/COSMIC_v2_SBS_GRCh37.txt',
    3.0: f'{module_path}/data/COSMIC_v3_SBS_GRCh37.txt',
    3.1: f'{module_path}/data/COSMIC_v3.1_SBS_GRCh37.txt',
    3.4: f'{module_path}/data/COSMIC_v3.4_SBS_GRCh37.txt'
}

def process_sample(args):
    i, col, COSMIC, threshold, mutation_count, R, significance_level = args
    try:
        best_columns, estimation_exposures = backward_elimination(
            col, COSMIC, threshold=threshold, mutation_count=mutation_count, R=R, significance_level=significance_level
        )
        return (i, best_columns, estimation_exposures)
    except Exception as e:
        print(f"Error processing sample {i}: {e}")
        return (i, None, None)

def cosmic_fit(samples_file, output_folder, threshold=0.01,
               mutation_count=None, R=100, significance_level=0.01, cosmic_version=3.4,
               drop_zeros_columns=False):
    """
     Analyzes genetic samples to determine their fit against known COSMIC (Catalogue Of Somatic Mutations In Cancer) signatures.
     This function processes each sample from a given input file to identify the best fitting mutational signatures based on the COSMIC database, and then exports the exposure estimates to a specified output folder.

     Parameters:
     - samples_file (str): Path to the file containing the genetic sample data to be analyzed.
     - output_folder (str): Path to the folder where the output files containing the assignment of samples to mutational signatures will be saved.
     - threshold (float, optional): The threshold value used to determine the fit of a sample to a signature. Default is 0.01.
     - mutation_count (int, optional): The total number of mutations to consider in the analysis. If None, the function will use all mutations available in the samples. Default is None.
     - R (int, optional): The number of iterations for the fitting algorithm. Higher values increase accuracy but also computational time. Default is 100.
     - significance_level (float, optional): The statistical significance level used in the fitting process. Default is 0.01.
     - cosmic_version (float, optional): The version of the COSMIC mutational signatures to use. Default is 3.4.
     - drop_zeros_columns (bool, optional): If True, columns with all zero values in the output matrix will be removed. Default is False.

     Returns:
     - None. The function saves the results in the specified output folder as "Assignment_Solution_Activities.csv". If drop_zeros_columns is True, columns containing only zeros will be excluded from the output.

     Note:
     - The function requires numpy for matrix operations and assumes the availability of `load_signatures_file`, `load_samples_file`, `process_sample`, and `utils.create_folder_if_not_exists` utility functions.
     - The output CSV file will contain the names of the signatures as the first row and the names of the samples as the first column. The rest of the matrix represents the estimated exposures of each sample to each signature.
     """
    if isinstance(cosmic_version, float):
        COSMIC, names_signatures = load_signatures_file(versions[cosmic_version])
    if isinstance(cosmic_version, str):
        COSMIC, names_signatures = load_signatures_file(cosmic_version)
    samples, names_patients = load_samples_file(samples_file)

    output = np.zeros((samples.shape[1], len(names_signatures)))
    output = np.vstack([names_signatures, output])

    for i in range(samples.shape[1]):
        _, best_columns, estimation_exposures = process_sample(
            (i, samples[:, i], COSMIC, threshold, mutation_count, R, significance_level))

        if best_columns is not None:
                for ind, col in enumerate(best_columns):
                    output[i + 1, col + 1] = estimation_exposures[0].squeeze()[ind]
                output[i + 1, 0] = names_patients[i]

                percent = (i + 1) / samples.shape[1]
                sys.stdout.write('\r')
                sys.stdout.write("[%-20s] %d%%" % ('=' * int(20 * percent), 100 * percent))
                sys.stdout.flush()
    utils.create_folder_if_not_exists(output_folder)
    if not drop_zeros_columns:
        np.savetxt(output_folder + "/Assignment_Solution_Activities.csv", output, delimiter=",", fmt='%s')
    else:
        is_non_zero_column = np.array([i != 0 if i > 0 else True for i in range(output.shape[1])])
        for col in range(1, output.shape[1]):
            if np.all(output[1:, col].astype(np.float32) == 0.0):
                is_non_zero_column[col] = False
        filtered_output = output[:, is_non_zero_column]
        np.savetxt(output_folder + "/Assignment_Solution_Activities.csv", filtered_output, delimiter=',', fmt='%s')
