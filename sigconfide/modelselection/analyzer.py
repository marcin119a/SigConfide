from sigconfide.modelselection.backward import backward_elimination
from sigconfide.utils.utils import load_samples_file
from sigconfide.utils import utils
import numpy as np
import os

module_path = os.path.dirname(utils.__file__)

versions = {3.4: f'{module_path}/data/COSMIC_v3_SBS_GRCh37.txt'}
def cosmic_fit(samples_file, output_file, threshold=0.01, mutation_count=None, R=100, significance_level=0.01, cosmic_version=3.4):
    COSMIC = np.genfromtxt(versions[cosmic_version], delimiter='\t', skip_header=1)
    names_signatures = np.genfromtxt(versions[cosmic_version], delimiter='\t', skip_header=0, max_rows=1, dtype=str)[1:]
    names_signatures = np.insert(names_signatures, 0, 'Samples', axis=0)
    COSMIC = np.delete(COSMIC, 0, axis=1)

    samples = load_samples_file(samples_file)
    print(samples)
    output = np.zeros((samples.shape[1], len(names_signatures)), dtype=COSMIC.dtype)
    output = np.vstack([names_signatures, output])
    import sys

    for i in range(samples.shape[1]):
        col = samples[:, i]
        try:
            best_columns, estimation_exposures = backward_elimination(col, COSMIC, threshold=threshold, mutation_count=mutation_count, R=R, significance_level=significance_level)

            for ind, col in enumerate(best_columns):
                output[i+1, col] = estimation_exposures[0].squeeze()[ind]
            percent = (i + 1) / samples.shape[1]
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%%" % ('=' * int(20 * percent), 100 * percent))
            sys.stdout.flush()

        except Exception as e:
            print(e)
            continue
    np.savetxt(output_file + "Assignment_Solution_Activities.csv", output, delimiter=",", fmt='%s')


