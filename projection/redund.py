import subprocess
import numpy as np
import os
from math import gcd
from functools import reduce
from fractions import Fraction
from decimal import Decimal, getcontext

# Function to simplify and round each float to significant digits
def simplify_and_round(value):
    return round(value, 6)

# Function to apply to every element of the numpy array
def process_matrix(matrix):
    return np.round(matrix,6)
    # Vectorize the function to apply it element-wise
    #vectorized_function = np.vectorize(simplify_and_round)
    
    # Apply the function to the entire matrix
    #return vectorized_function(matrix)

# Function to compute GCD of a list
def row_gcd(row):
    # Convert floats to fractions
    fractions = [Fraction(x).limit_denominator() for x in row]
    # Find the least common multiple of the denominators
    denominators = [frac.denominator for frac in fractions]
    common_denominator = reduce(gcd, denominators)
    # Convert fractions to integers with the common denominator
    numerators = [frac.numerator * (common_denominator // frac.denominator) for frac in fractions]
    # Find the GCD of the numerators
    common_numerator_gcd = reduce(gcd, numerators)
    # Compute the greatest common divisor of the row as a float
    return common_numerator_gcd / common_denominator

# Function to check if any number in the row exceeds 7 digits
def exceeds_seven_digits(row):
    return np.any(np.abs(num) >= 10**6 for num in row)

# Function to scale down a row so all numbers are less than 7 digits
def scale_row(row):
    max_val = np.max(np.abs(row))
    scaling_factor = 10**(len(str(int(max_val))) - 6)
    #return (row / scaling_factor).astype(int)
    return row / scaling_factor
# Normalize matrix by GCD and ensure no number exceeds 7 digits
# def normalize_matrix_by_gcd_and_scale(matrix):
#     # matrix = np.round(matrix, decimals=6)
#     result = np.zeros_like(matrix, dtype=object)
#     for i, row in enumerate(matrix):
#         # Normalize the row by its GCD
#         #row_gcd_value = row_gcd(row)
#         #if row_gcd_value > 0:  # Avoid division by zero
#         #    normalized_row = row / row_gcd_value
#         #else:
#         #    normalized_row = row
#         normalized_row = row
#         # Check if any number in the row exceeds 7 digits
#         # if exceeds_seven_digits(normalized_row):
#         #     normalized_row = scale_row(normalized_row)

#         result[i] = normalized_row
#     print(result)
#     return result

# def normalize_matrix_by_gcd_and_scale(row, threshold=1e7):
#     """
#     Normalizes the row if any element in the row has an absolute value
#     greater than or equal to the threshold. The normalization divides every
#     element in the row by the maximum absolute value in the row, preserving
#     the ratios. The resulting row is rounded to 6 decimal places.
    
#     Parameters:
#       row: 1D NumPy array of numbers.
#       threshold: the threshold value (default 1e7) to decide if normalization is needed.
      
#     Returns:
#       A normalized (or unmodified) row as a 1D NumPy array.
#     """
#     # Check if any element meets the threshold.
#     if np.any(np.abs(row) >= threshold):
#         # Avoid division by zero by ensuring the maximum is nonzero.
#         max_val = np.max(np.abs(row))
#         if max_val != 0:
#             normalized_row = row / max_val
#         else:
#             normalized_row = row
#         # Round the normalized values to 6 decimal places.
#         return np.round(normalized_row, 6)
#     else:
#         # Optionally, round the row even if not normalized.
#         print(row)
#         return np.round(row, 6)

def round_fraction_to_fraction(frac, ndigits=4):
    """
    Rounds a Fraction to a given number of decimal places,
    returning a Fraction that approximates the original value.
    """
    multiplier = 10 ** ndigits
    # Multiply the fraction, round the result, then create a new Fraction.
    return Fraction(round(frac * multiplier), multiplier)

def normalize_matrix_by_gcd_and_scale(row, threshold=1e15):
    """
    Normalizes the row if any element in the row has an absolute value
    greater than or equal to the threshold. The normalization divides every
    element in the row by the maximum absolute value in the row, preserving
    the ratios. The resulting row is rounded to 6 decimal places.
    
    Parameters:
      row: 1D NumPy array of numbers.
      threshold: the threshold value (default 1e7) to decide if normalization is needed.
      
    Returns:
      A normalized (or unmodified) row as a 1D NumPy array.
    """
    # Check if any element meets the threshold.
    if np.any(np.abs(row) >= threshold):
        # Avoid division by zero by ensuring the maximum is nonzero.
        max_val = np.max(np.abs(row))
        if max_val != 0:
            normalized_row = row / max_val
        else:
            normalized_row = row
        # Round the normalized values to 6 decimal places.
        #normalized_row = np.array([round_fraction_to_fraction(frac) for frac in normalized_row], dtype=object)
        return normalized_row
    else:
        # Optionally, round the row even if not normalized.
        #row = np.array([round_fraction_to_fraction(frac) for frac in row], dtype=object)
        return row
    
def removeRowsWithZeros(matrix: np.matrix, verbose = True):
    '''
    deletes rows from a matrix that contain only zeros
    '''
    rowsToDelete = []
    numRows = matrix.shape[0]
    # numCols = matrix.shape[1]
    for row in range(numRows):
        toDelete = True
        for value in np.nditer(matrix[row,:]):
            if abs(value) >= 10**-6:
                toDelete = False
                break
        if toDelete:
            rowsToDelete.append(row)

    if len(rowsToDelete) > 0:
        rowsToKeep = list(range(numRows))
        for i in rowsToDelete:
            rowsToKeep.remove(i)
        # rowsToKeep = list(range(numRows)) #.remove(rowsToDelete)
        # numRows -= len(rowsToDelete)
        # resMatrix = np.delete(matrix,rowsToDelete) # .reshape(numRows, numCols)
        resMatrix = matrix[rowsToKeep, :]
    if verbose:
        print(resMatrix)
    return resMatrix

def redund(n_cores, path_mplrs, inputFilePath, outFilePath, verbose=True):
    """
    Performs mplrs redund on given .ine file.
    """
    with open(inputFilePath, 'a') as file:
        file.write('redund 0 0\n')
    #original_cmd = "mpirun -np 3 /opt/lrslib/v072/mplrs -redund ./h_representation.ine > input_postredund.ine"
    # cmd = ["mpirun", "-np", str(n_cores), path_mplrs, "-minrep", inputFilePath, outFilePath, "-lastp", str(1000), "-lastrows", str(10000), "-rows", str(500), "-time", str(172800)]
    cmd = ["mpirun", "-np", str(n_cores), path_mplrs, inputFilePath, outFilePath]
    if verbose:
        subprocess.run(cmd)
    else:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)


def redundIterative(n_cores, path_mplrs, inputFilePath, outFilePath, tempFolder, chunkSize=10000, verbose=True):
    """
    Performs mplrs redund iteratively on the given .ine file and combines the results.
    """
    
    # Read the original file
    with open(inputFilePath, 'r') as f:
        lines = f.readlines()

    # Find the 'begin' section and the number of rows and columns
    start_index = lines.index('begin\n') + 1
    header_line = lines[start_index].strip()  # Get the first line after 'begin'
    
    # Extract the number of rows and columns from the header line
    num_rows, num_columns, _ = header_line.split()  # Splits into 3 parts
    num_rows = int(num_rows)
    num_columns = int(num_columns)
    
    # Extract rows data (after the 'begin' and header line, before 'end')
    end_index = lines.index('end\n')
    rows = lines[start_index + 1:end_index]

    # Split the rows into chunks of 'chunk_size' rows
    chunks = [rows[i:i + chunkSize] for i in range(0, len(rows), chunkSize)]

    # List to store paths of the output files for combining later
    temp_out_files = []

    # Write each chunk into a separate file
    for idx, chunk in enumerate(chunks, 1):
        # Write the chunk to the new file
        tempRedundFile = os.path.join(tempFolder, f"tempRedund_{idx}.ine")
        outTempRedundFile = os.path.join(tempFolder, f"outTempRedund_{idx}.ine")
        with open(tempRedundFile, 'w') as out_f:
            out_f.write("poly\n")
            out_f.write("H-representation\n")
            out_f.write("begin\n")
            out_f.write(f"{len(chunk)} {num_columns} rational\n")  # Write the number of columns dynamically
            out_f.writelines(chunk)
            out_f.write("end\n")
        
        # Run the redund process
        redund(n_cores, path_mplrs, tempRedundFile, outTempRedundFile, verbose=True)
        
        # Add the output file path to the list
        temp_out_files.append(outTempRedundFile)
        
        print(f"Written: {outTempRedundFile}")

    num_new_rows = 0
    for temp_file in temp_out_files:
        with open(temp_file, 'r') as f:
            lines = f.readlines()
    
        # Find the 'begin' section and the number of rows and columns
        start_index = lines.index('begin\n') + 1
        header_line = lines[start_index].strip()  # Get the first line after 'begin'
        
        # Extract the number of rows and columns from the header line
        num_rows, num_columns, _ = header_line.split()  # Splits into 3 parts
        num_rows = int(num_rows)
        num_new_rows += num_rows
    
    # Combine the results into the final output file
    with open(outFilePath, 'w') as final_out_file:
        final_out_file.write("poly\n")
        final_out_file.write("H-representation\n")
        final_out_file.write("begin\n")
        
        final_out_file.write(f"{num_new_rows} {num_columns} rational\n")  # Write header to the final file
        
        # Write rows from each temp file
        for temp_file in temp_out_files:
            # Read the original file
            with open(temp_file, 'r') as f:
                lines = f.readlines()
        
            # Find the 'begin' section and the number of rows and columns
            start_index = lines.index('begin\n') + 1
            header_line = lines[start_index].strip()  # Get the first line after 'begin'
            
            # Extract the number of rows and columns from the header line
            num_rows, num_columns, _ = header_line.split()  # Splits into 3 parts
            num_rows = int(num_rows)
            num_columns = int(num_columns)
            
            # Extract rows data (after the 'begin' and header line, before 'end')
            end_index = lines.index('end\n')
            rows = lines[start_index + 1:end_index]
            final_out_file.writelines(rows)
        
        final_out_file.write("end\n")
    
    print(f"Final combined output written to: {outFilePath}")





    