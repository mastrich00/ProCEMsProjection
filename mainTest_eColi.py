import numpy as np
import copy
import subprocess
from zipfile import ZipFile 
import re
import os
from datetime import datetime
import time
#import efmtool
#from projection.marianneBianca import get_blocked_reactions, rm_reactions, split_all_reversible_reactions, indicate_exchange
from projection.projection import runMarashiWithPolco, runMarashiWithMPLRS, runFELWithPolco
#from gmpy2 import mpq, mpfr
from argparse import ArgumentParser, ArgumentTypeError
from projection.logging_config import logger

polcoPath = "polco.jar"
mplrsPath = "mplrs"
import pandas as pd

matrix = np.loadtxt("ecmtool/test_network.txt", delimiter=',')
print(matrix)
input_indices = [15, 21, 25, 27, 31, 33, 35]
output_indices = [12, 13, 14, 16, 17, 18, 19, 20, 22, 23, 24, 26, 28, 29, 30, 32, 34, 36, 37, 38, 39]
rowIndex = [15, 21, 25, 27, 31, 33, 35, 12, 13, 14, 16, 17, 18, 19, 20, 22, 23, 24, 26, 28, 29, 30, 32, 34, 36, 37, 38, 39]

# Get the number of rows in the matrix
num_rows = matrix.shape[0]

# Create new columns for input indices
for index in input_indices:
    new_column = np.zeros((num_rows, 1))  # Column of zeros
    new_column[index, 0] = 1             # Set the entry for the input index to 1
    matrix = np.hstack((matrix, new_column))  # Add the new column to the matrix

# Create new columns for output indices
for index in output_indices:
    new_column = np.zeros((num_rows, 1))  # Column of zeros
    new_column[index, 0] = -1            # Set the entry for the output index to -1
    matrix = np.hstack((matrix, new_column))  # Add the new column to the matrix

# Get the indices of the newly created columns
new_column_indices = list(range(matrix.shape[1] - len(input_indices) - len(output_indices), matrix.shape[1]))

# Return the new column indices
print(new_column_indices, matrix.shape)
print(matrix)


## Find column indices where any value in the specified rows is non-zero
#column_indices = np.where(np.any(matrix[rowIndex, :] != 0, axis=0))[0]
#column_indices = list(column_indices)
## Output the result
#print("Column indices with non-zero values in specified rows:", column_indices)

proCEMs, efps = runMarashiWithPolco(matrix, new_column_indices, "testResults/polco_ecoli/", 8, False, True)
#proCEMs, efps = runMarashiWithMPLRS(matrix, column_indices, "testResults/mplrs_ecoli/", mplrsPath, 4, False, True)
# proCEMs, efps = runMarashiWithMPLRS(inpMatrix, inpDims, "testResults/mplrs/", mplrsPath)
# proCEMs, efps = runFELWithPolco(inpMatrix, inpDims, "~/secondDraft/testResults/fel/", mplrsPath)
#logger.info(f"proCEMS: {len(proCEMs)}")
#logger.info(f"efps: {len(efps)}")
