# SigConfide

SigConfide facilitates the allocation of established mutational signatures to specific samples and individual somatic mutations. 
This tool reanalyzes various reference mutational signatures, such as those from [COSMIC](https://cancer.sanger.ac.uk/signatures/), in addition to personalized signature databases. 
The process of reanalyzing known mutational signatures involves a numerical optimization method, 
which not only discerns the array of active mutational signatures within a given sample 
but also measures the quantity of mutations attributed to each detected signature within that sample.
In addition, it performs a bootstrapping approach for a single sample and selects the signatures most relevant to it. 
In addition, the tool can be run on an entire sample of patients. 

# Installation

In this version install the current package from github:
```python
pip install -e git+https://github.com/marcin119a/sigconfide.git#egg=sigconfide
```

# Running 


#### COSMIC Fit Function Description

The `cosmic_fit` function is designed to analyze genetic samples and evaluate their alignment with known COSMIC (Catalogue Of Somatic Mutations In Cancer) signatures. 
It meticulously processes each sample provided in the input file to pinpoint the most compatible mutational signatures from the COSMIC database. Subsequently,
it calculates the exposure estimates for these signatures and outputs the results into a designated folder.

### Parameters

The function's behavior can be customized through the following parameters:

| Parameter           | Type         | Description                                                                    | Default |
|---------------------|--------------|--------------------------------------------------------------------------------|---------|
| `samples_file`      | str          | Path to the file with genetic sample data.                                     | -       |
| `output_folder`     | str          | Path to the directory for saving output files.                                 | -       |
| `threshold`         | float        | Threshold for matching a sample to a signature.                                | 0.01    |
| `mutation_count`    | int          | Number of mutations to analyze. If `None`, all available mutations in the samples are used. | None    |
| `R`                 | int          | Iteration count for the fitting algorithm, affecting both accuracy and computation time. | 100     |
| `significance_level`| float        | Statistical significance level for the fitting process.                        | 0.01    |
| `cosmic_version`    | float or str | Version of the COSMIC mutational signatures database (2.0, 3.4, 3.1, 3.0) to use or signatures file. | 3.4     |
| `drop_zeros_columns`| bool         | If `True`, excludes columns with zero values from the output matrix.           | False   |

### Output

The function does not return any values but instead writes the analysis results to a CSV file named `Assignment_Solution_Activities.csv` in the `output_folder`.
The CSV file's first row lists the signatures, the first column lists the sample names, and the subsequent cells contain the estimated exposure levels.


### Examples