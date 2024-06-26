import numpy as np
import os

def FrobeniusNorm(M, P, E):
    return np.sqrt(np.sum((M - np.dot(P, E))**2))


def is_wholenumber(x, tol=1e-15):
    return np.abs(x - np.round(x)) < tol

def detect_format(line):
    # Check if the line contains tabs - this suggests TSV format
    if '\t' in line:
        # Additionally, check if square brackets are present
        if '[' in line and (']' in line):
            return 'Mutated TSV Format', '\t'
        return 'TSV Format', '\t'

    # Check if the line contains commas - this suggests CSV format
    elif ',' in line:
        if ('[' in line) and (']' in line):
            return 'Mutated CSV Format', ','
        return 'CSV Format', ','

    # If none of the above were detected, return an unknown value
    return 'Unknown Format', None

def load_samples_file(file_name):
    with open(file_name, 'r') as file:
        csv_line = ''.join(file.readlines()).strip()
        format, sep = detect_format(csv_line)
        file.seek(0)
    if format == 'Unknown Format':
        raise ValueError('Unknown Format')

    samples = np.genfromtxt(file_name, delimiter=sep, skip_header=1)
    patient_names = np.genfromtxt(file_name, delimiter=sep, skip_header=0, max_rows=1, dtype=str).squeeze()

    if format == 'TSV Format' or format == 'Mutated TSV Format':
        samples = np.delete(samples, 0, axis=1)
        patients_names = patient_names[1:]

    if format == 'CSV Format' or format == 'Mutated CSV Format':
        samples = np.delete(samples, [0, 1], axis=1)
        patients_names = patient_names[2:]

    return samples, patients_names

def load_signatures_file(file_name):
    with open(file_name, 'r') as file:
        csv_line = ''.join(file.readlines()).strip()
        format, sep = detect_format(csv_line)
        file.seek(0)
    signatures = np.genfromtxt(file_name, delimiter=sep, skip_header=1)
    names_signatures = np.genfromtxt(file_name, delimiter=sep, skip_header=0, max_rows=1, dtype=str)[1:]
    names_signatures = np.insert(names_signatures, 0, 'Samples', axis=0)
    signatures = np.delete(signatures, 0, axis=1)

    return signatures, names_signatures

def create_folder_if_not_exists(folder_path):
    try:
        os.makedirs(folder_path, exist_ok=True)
    except OSError as error:
        print(f"Błąd podczas tworzenia katalogu {folder_path}: {error}")

