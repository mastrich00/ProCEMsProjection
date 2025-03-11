import subprocess
import numpy as np
import os
from math import gcd
from functools import reduce
from fractions import Fraction
from decimal import Decimal, getcontext

#def round_fraction_to_fraction(frac, ndigits=4):
#    """
#    Rounds a Fraction to a given number of decimal places,
#    returning a Fraction that approximates the original value.
#    """
#    multiplier = 10 ** ndigits
#    # Multiply the fraction, round the result, then create a new Fraction.
#    return Fraction(round(frac * multiplier), multiplier)

#def normalize_matrix_by_gcd_and_scale(row, threshold=1e15):
#    """
#    Normalizes the row if any element in the row has an absolute value
#    greater than or equal to the threshold. The normalization divides every
#    element in the row by the maximum absolute value in the row, preserving
#    the ratios. The resulting row is rounded to 6 decimal places.
#    
#    Parameters:
#      row: 1D NumPy array of numbers.
#      threshold: the threshold value (default 1e7) to decide if normalization is needed.
#      
#    Returns:
#      A normalized (or unmodified) row as a 1D NumPy array.
#    """
#    # Check if any element meets the threshold.
#    if np.any(np.abs(row) >= threshold):
#        # Avoid division by zero by ensuring the maximum is nonzero.
#        max_val = np.max(np.abs(row))
#        if max_val != 0:
#            normalized_row = row / max_val
#        else:
#            normalized_row = row
#        # Round the normalized values to 6 decimal places.
#        #normalized_row = np.array([round_fraction_to_fraction(frac) for frac in normalized_row], dtype=object)
#        return normalized_row
#    else:
#        # Optionally, round the row even if not normalized.
#        #row = np.array([round_fraction_to_fraction(frac) for frac in row], dtype=object)
#        return row

def redund(n_cores, path_mplrs, inputFilePath, outFilePath, verbose=True):
    """
    Performs mplrs redund on given .ine file.
    """
    with open(inputFilePath, 'a') as file:
        file.write('redund 0 0\n')
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





    