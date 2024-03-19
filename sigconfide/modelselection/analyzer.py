from sigconfide.modelselection.backward import backward_elimination
from sigconfide.utils.utils import load_samples_file, load_signatures_file
from sigconfide.utils import utils
import numpy as np
import os
import sys
from concurrent.futures import ThreadPoolExecutor

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

def cosmic_fit(samples_file, output_file, threshold=0.01,
               mutation_count=None, R=100, significance_level=0.01, cosmic_version=3.4,
               drop_zeros_columns=False):
    COSMIC, names_signatures = load_signatures_file(versions[cosmic_version])
    samples, names_patients = load_samples_file(samples_file)

    output = np.zeros((samples.shape[1], len(names_signatures)))
    output = np.vstack([names_signatures, output])

    args_list = [(i, samples[:, i], COSMIC, threshold, mutation_count, R, significance_level) for i in range(samples.shape[1])]

    with ThreadPoolExecutor() as executor:
        for i, best_columns, estimation_exposures in executor.map(process_sample, args_list):
            if best_columns is not None:
                for ind, col in enumerate(best_columns):
                    output[i + 1, col + 1] = estimation_exposures[0].squeeze()[ind]
                output[i + 1, 0] = names_patients[i]

                percent = (i + 1) / samples.shape[1]
                sys.stdout.write('\r')
                sys.stdout.write("[%-20s] %d%%" % ('=' * int(20 * percent), 100 * percent))
                sys.stdout.flush()
    if not drop_zeros_columns:
        np.savetxt(output_file + "/Assignment_Solution_Activities.csv", output, delimiter=",", fmt='%s')
    else:
        is_non_zero_column = np.array([i != 0 if i > 0 else True for i in range(output.shape[1])])
        for col in range(1, output.shape[1]):
            if np.all(output[1:, col].astype(np.float32) == 0.0):
                is_non_zero_column[col] = False
        filtered_output = output[:, is_non_zero_column]
        np.savetxt(output_file + "/Assignment_Solution_Activities.csv", filtered_output, delimiter=',', fmt='%s')
