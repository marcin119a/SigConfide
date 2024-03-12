from sigconfide.modelselection.backward import backward_elimination
from sigconfide.utils.utils import load_samples_file, load_signatures_file
from sigconfide.utils import utils
import numpy as np
import os

module_path = os.path.dirname(utils.__file__)

versions = {3.4: f'{module_path}/data/COSMIC_v3_SBS_GRCh37.txt'}
def cosmic_fit(samples_file, output_file, threshold=0.01, mutation_count=None, R=100, significance_level=0.01, cosmic_version=3.4):
    COSMIC, names_signatures = load_signatures_file(versions[cosmic_version])

    samples, names_patients = load_samples_file(samples_file)
    output = np.zeros((samples.shape[1], len(names_signatures)))
    output = np.vstack([names_signatures, output])
    import sys

    for i in range(samples.shape[1]):
        col = samples[:, i]
        try:
            best_columns, estimation_exposures = backward_elimination(col, COSMIC, threshold=threshold, mutation_count=mutation_count, R=R, significance_level=significance_level)

            for ind, col in enumerate(best_columns):
                output[i+1, col+1] = estimation_exposures[0].squeeze()[ind]
            output[i+1, 0] = names_patients[i]
            percent = (i + 1) / samples.shape[1]
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%%" % ('=' * int(20 * percent), 100 * percent))
            sys.stdout.flush()

        except Exception as e:
            print(e)
            continue
    np.savetxt(output_file + "Assignment_Solution_Activities.csv", output, delimiter=",", fmt='%s')


