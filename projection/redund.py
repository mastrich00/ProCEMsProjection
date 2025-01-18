import subprocess
import numpy as np
import os

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





    