import sys
import numpy as np
import pandas as pd
from SparCC import SparCC

def run_sparcc(input_file, output_file):
    """
    Runs SparCC on the input dataset and outputs the covariance matrix.

    Args:
        input_file (str): The path to the CSV file containing the input sequence count data.
        output_file (str): The path to the CSV file where the covariance matrix will be saved.

    Returns:
        None. The covariance matrix is saved to `output_file`.

    Example:
        To run SparCC on a dataset:

        ```bash
        python sparcc_script.py input.csv output.csv
        ```

    Notes:
        This script requires SparCC and its dependencies to be installed in a Python 2.7 environment.
    """
    data = pd.read_csv(input_file, index_col=0)
    sparcc = SparCC(data)
    cov = sparcc.covariance_
    pd.DataFrame(cov, index=data.columns, columns=data.columns).to_csv(output_file)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    run_sparcc(input_file, output_file)
