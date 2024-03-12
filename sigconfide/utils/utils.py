import numpy as np

def FrobeniusNorm(M, P, E):
    return np.sqrt(np.sum((M - np.dot(P, E))**2))


def is_wholenumber(x, tol=1e-15):
    return np.abs(x - np.round(x)) < tol

def load_and_process_data(patient_index, mutational_profiles, predf_mutational_signatures):
    # Load and process the first file
    profiles = np.genfromtxt(mutational_profiles, delimiter=',', skip_header=1)
    if profiles.size == 0:
        raise ValueError(f"Empty data in {mutational_profiles}")
    profiles = np.delete(profiles, 0, axis=1)

    profile = profiles
    if patient_index is not None:
        profile = profiles[:, patient_index]

    # Load and process the second file
    signatures = np.genfromtxt(predf_mutational_signatures, delimiter=',', skip_header=1)
    if signatures.size == 0:
        raise ValueError(f"Empty data in {predf_mutational_signatures}")
    signatures = np.delete(signatures, 0, axis=1)

    # Return processed data
    return profile, signatures


def detect_format(line):
    # Check if the line contains tabs - this suggests TSV format
    if '\t' in line:
        # Additionally, check if square brackets are present
        if '[' in line and (']' in line):
            return 'Mutated TSV Format', '\t'
        return 'TSV Format', '\t'

    # Check if the line contains commas - this suggests CSV format
    elif ',' in line:
        return 'CSV Format', ','

    # If none of the above were detected, return an unknown value
    return 'Unknown Format', None

def load_samples_file(file_name):
    with open(file_name, 'r') as file:
        csv_line = file.readline().strip()
        format, sep = detect_format(csv_line)
        file.seek(0)
        samples = np.genfromtxt(file, delimiter=sep, skip_header=1)

    if format == 'TSV Format':
        samples = np.delete(samples, 0, axis=1)

    if format == 'CSV Format':
        samples = np.delete(samples, [0, 1], axis=0)

    if format == 'Unknown Format':
        raise ValueError('Unknown Format')

    return samples

