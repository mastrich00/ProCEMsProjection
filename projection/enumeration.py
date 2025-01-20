import numpy as np
import subprocess
from zipfile import ZipFile 
import re
import os
from dotenv import load_dotenv
from projection.logging_config import logger
from fractions import Fraction

if not os.getenv("POLCO_PATH"):
    load_dotenv("env") # load env
    
POLCO_PATH = os.getenv("POLCO_PATH")

def as_fraction(number, approximation=None):
    """
    Takes integers or floats and converts them to rational strings.
    
    Parameters
    ----------
    number : int, float
        Give number for conversion to a rational.
    approximation : int default None
        Sets the border for number a denominator can have.
        If border of 1e12 is set, you get a rational unequal to zero for numbers down to 1e-12 and zero for numbers lower than that.
        Approximation default is 1e6. Keep attention. Less strict approximation borders can lead to huge numbers which slow down the calculation intensely.
    
    Returns
    -------
    str
        A rational or integer number as a string.
    
    Examples
    --------
    >>> into_fractions(5)
    '5'
    >>> into_fractions(0)
    '0'
    >>> into_fractions(0.5)
    '1/2'
    >>> into_fractions(1/6)
    '1/6'
    >>> into_fractions(1/6, approximation=1e20)
    '6004799503160661/36028797018963968'
    >>> into_fractions(1/1e6)
    '1/1000000'
    >>> into_fractions(1/1e7)
    '0'
    """
    if approximation:    
        return str(Fraction(number).limit_denominator(int(float(approximation))))
    else:
        return str(Fraction(number).limit_denominator())

def doubleDescriptionINE(ineFile: str, outputDir: str, outputFile: str, stage = "", verbose=True, numThreads=28):
    '''
    - performs polco's double description method on the system specified by the two input matrices: 
        - inequality matrix: ${iqMatrix} [NOT optional]
        - equality matrix: ${eqMatrix} [optional]
    - input matrices are automatically converted to a text-file, specified under ${outputDir}/${iqFile} and ${outputDir}/${eqFile}
    - polcos rayMatrix will be unzipped, fractions are replaced and final matrix is written to ${outputDir}/${outputFile}
    - ${stage} is used for logging purposes
    '''
       
    if(".zip" not in outputFile):
        outputFile = outputFile + ".zip"
    outputFilePath = os.path.join(outputDir, outputFile)
    logPath = outputDir + stage + "_polco.log"

    cmd = ["java", "-Xms10g", "-Xmx15g", "-jar", POLCO_PATH, "-kind", "cdd", "-in", ineFile, "-out", "zip", outputFilePath, "-maxthreads", str(numThreads)]

    if verbose:
        print(cmd)
    if verbose:
        with open(logPath, 'w') as logFile:
            subprocess.run(cmd, stdout=logFile)
    else:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)
    
    with ZipFile(outputFilePath, 'r') as zObject: 
        zObject.extract("R", path=outputDir) 
        zObject.close()

    rOutPath = os.path.join(outputDir, stage + "_R")
    os.rename(os.path.join(outputDir,"R"), rOutPath) 
    
    raysWithOutFractionsFilePath = os.path.join(outputDir, stage + "_RwoFractions")
    if os.path.exists(raysWithOutFractionsFilePath):
        os.remove(raysWithOutFractionsFilePath)
    raysWithOutFractionsFile = open(raysWithOutFractionsFilePath, "a")
        
    with open(rOutPath, "r") as file: # replace fractions with decimal numbers
        for line in file.readlines():
            results = re.findall(r"(\d+)\/(\d+)", line)
            for firstNumber,secondNumber in results:
                decimalNumber = float(firstNumber) / float(secondNumber) 
                line = line.replace(str(firstNumber) + "/" + str(secondNumber), str(decimalNumber)) # TODO: replace string with regex match
                
            raysWithOutFractionsFile.write(line)
    raysWithOutFractionsFile.close()        
    rayMatrix = np.loadtxt(raysWithOutFractionsFilePath)
    
    # rayMatrix = rayMatrix.transpose()
    if(verbose):
        logger.info("resulting Rays:")
        logger.info(rayMatrix)
    return rayMatrix

def doubleDescription(iqMatrix: np.matrix, iqFile: str, outputDir: str, outputFile: str, stage = "", eqMatrix = np.zeros(0), eqFile = "", verbose=True, numThreads=28):
    '''
    - performs polco's double description method on the system specified by the two input matrices: 
        - inequality matrix: ${iqMatrix} [NOT optional]
        - equality matrix: ${eqMatrix} [optional]
    - input matrices are automatically converted to a text-file, specified under ${outputDir}/${iqFile} and ${outputDir}/${eqFile}
    - polcos rayMatrix will be unzipped, fractions are replaced and final matrix is written to ${outputDir}/${outputFile}
    - ${stage} is used for logging purposes
    '''
    if(eqFile != ""):
        eqFilePath = os.path.join(outputDir, eqFile)
        
        logger.info("Equality Matrix:")
        logger.info(eqMatrix)
        np.savetxt(eqFilePath, eqMatrix, fmt='%.6f', delimiter="\t") # TODO: fmt aendern
    
    iqFilePath = os.path.join(outputDir, iqFile)
    logger.info("Inequality Matrix:")
    logger.info(iqMatrix)
    np.savetxt(iqFilePath, iqMatrix, fmt='%.6f', delimiter="\t")
    
    if(".zip" not in outputFile):
        outputFile = outputFile + ".zip"
    outputFilePath = os.path.join(outputDir, outputFile)
    
    # date = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    logPath = outputDir + stage + "_polco.log"
    # TODO: make polco path relative
    #cmd = ["java", "-jar", POLCO_PATH, "-kind", "text", "-iq", iqFilePath]
    # if(eqFile != ""):
    #     cmd.append("-eq " + eqFilePath)

    # cmd.extend(["-out", "zip", outputFilePath])
    # cmd = ["java", "-jar", POLCO_PATH, "-kind", "text", "-iq", iqFilePath]
    if(eqFile != ""):
        cmd = ["java", "-Xms10g", "-Xmx13g", "-jar", POLCO_PATH, "-kind", "text", "-eq", eqFilePath, "-iq", iqFilePath, "-out", "zip", outputFilePath, "-maxthreads", str(numThreads)]
    else:
        cmd = ["java", "-Xms10g", "-Xmx13g", "-jar", POLCO_PATH, "-kind", "text", "-iq", iqFilePath, "-out", "zip", outputFilePath, "-maxthreads", str(numThreads)]

    
    if verbose:
        print(cmd)
    if verbose:
        with open(logPath, 'w') as logFile:
            subprocess.run(cmd, stdout=logFile)
    else:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)
    
    with ZipFile(outputFilePath, 'r') as zObject: 
        zObject.extract("R", path=outputDir) 
        zObject.close()

    rOutPath = os.path.join(outputDir, stage + "_R")
    os.rename(os.path.join(outputDir,"R"), rOutPath) 
    
    raysWithOutFractionsFilePath = os.path.join(outputDir, stage + "_RwoFractions")
    if os.path.exists(raysWithOutFractionsFilePath):
        os.remove(raysWithOutFractionsFilePath)
    raysWithOutFractionsFile = open(raysWithOutFractionsFilePath, "a")
        
    with open(rOutPath, "r") as file: # replace fractions with decimal numbers
        for line in file.readlines():
            results = re.findall(r"(\d+)\/(\d+)", line)
            for firstNumber,secondNumber in results:
                decimalNumber = float(firstNumber) / float(secondNumber) 
                line = line.replace(str(firstNumber) + "/" + str(secondNumber), str(decimalNumber)) # TODO: replace string with regex match
                
            raysWithOutFractionsFile.write(line)
    raysWithOutFractionsFile.close()
    rayMatrix = np.loadtxt(raysWithOutFractionsFilePath)
    
    # rayMatrix = rayMatrix.transpose()
    if(verbose):
        logger.info("resulting Rays:")
        logger.info(rayMatrix)
    return rayMatrix

def convertMatrixToHRepresentation(matrix: np.matrix, filePath: str):
    '''
    takes a numpy matrix and converts it to its H-Representation
    '''
    if os.path.exists(filePath):
        os.remove(filePath)

    with open(filePath, "w") as file:
        file.write("poly\n")
        file.write("H-representation\n")
        file.write("begin\n")
        line = str(matrix.shape[0]) + " " + str((matrix.shape[1]+1)) + " rational\n" 
        file.write(line)
        print(matrix.shape)
        for i in range(matrix.shape[0]):
            line = "0" # no coefficient, TODO: check if zero
            #for value in list(matrix[i,:]):
            #print(matrix[i,0])
            for j in range(matrix.shape[1]):
                #print(matrix[i][j])
                line += " "
                value = matrix[i,j]
                if value == int(value): 
                    insertValue = int(value) 
                else:
                    insertValue = as_fraction(value)#Fraction(float(value)).limit_denominator()
                line += str(insertValue)
            line += "\n"
            file.write(line)
        file.write("end\n")
        file.close()

def convertMatrixToVRepresentation(matrix: np.matrix, filePath: str):
    '''
    takes a numpy matrix and converts it to its V-Representation
    '''
    if os.path.exists(filePath):
        os.remove(filePath)

    with open(filePath, "w") as file:
        file.write("poly\n")
        file.write("V-representation\n")
        file.write("begin\n")
        line = str(matrix.shape[0] + 1) + " " + str((matrix.shape[1]+1)) + " rational\n" 
        file.write(line)
        line = "1" + " 0" * matrix.shape[1]
        file.write(line)
        file.write("\n")
        print(matrix.shape)
        for i in range(matrix.shape[0]):
            line = "0" # no coefficient, TODO: check if zero
            #for value in list(matrix[i,:]):
            #print(matrix[i,0])
            for j in range(matrix.shape[1]):
                #print(matrix[i][j])
                line += " "
                value = matrix[i,j]
                if value == int(value): 
                    insertValue = int(value) 
                else:
                    insertValue = as_fraction(value)#Fraction(float(value)).limit_denominator()
                line += str(insertValue)
            line += "\n"
            file.write(line)
        file.write("end\n")
        file.close()

def mplrs_conversion(n_cores, path_mplrs: str, hFilePath: str, vFilePath: str, verbose=True):
    """Convert H-representation to V-representation."""
    cmd = ["mpirun", "-np", str(n_cores), path_mplrs, hFilePath, vFilePath]
    logger.info('Run conversion')
    if verbose:
        subprocess.run(cmd)
    else:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)

def convertHtoVrepresentation(hFilePath: str, vFilePath: str, verbose = True):
    '''alternative using only lrs'''
    cmd = ["lrs", hFilePath, vFilePath]
    if verbose:
        subprocess.run(cmd)
    else:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)

def mplrs_project(smatrix, ex_reactions, n_cores, path_mplrs, postRedundFile, outProjectedFile, rows=60, lastp=10, lastrows=10, verbose=True):
    """
    Performs mplrs project.
    """
    # calculate n times
    # if we want the last H-representation
    n = smatrix.shape[1] - len(ex_reactions)
    print(f'Start projection of postredund file.')
    
    # mplrs project
    # cmd = ["mpirun", "-np", str(n_cores), path_mplrs, postRedundFile, outProjectedFile, "-rows", str(rows),  "-lastp", str(lastp), "-lastrows", str(lastrows)]
    cmd = ["mpirun", "-np", str(n_cores), path_mplrs, "-fel", postRedundFile, outProjectedFile]
    for i in range(0, n):
        logger.info(f'Projection-Run {i + 1}/{n}')
        if verbose:
            subprocess.run(cmd)
        else:
            subprocess.run(cmd, stdout=subprocess.DEVNULL)
        os.replace(outProjectedFile, postRedundFile)
    os.replace(postRedundFile, outProjectedFile)
    
# def getMatrixFromHrepresentation(hFilePath: str):
#     '''
#     builds a matrix from a H-Representation
#     '''
#     beginOfRays = False
#     rayMatrix = []
#     with open(hFilePath, "r") as file:
#         for line in file.readlines():
#             if "end" in line:
#                 break
#             if "begin" in line:
#                 beginOfRays = True
#                 continue
#             if beginOfRays:
#                 if "*" in line:
#                     continue
#                 matches = re.findall("-?\d+", line) # get all numbers from line
#                 if matches[0] != "0":
#                     continue
#                 row = []
#                 matches.pop(0) # remove first number that indicates rhs
#                 for match in matches:
#                     row.append(match)
#                 # print(row)
#                 rayMatrix.append(row)
#         file.close()

#     # print(np.matrix(rayMatrix, dtype="float64"))    
#     return np.matrix(rayMatrix, dtype="float64") #.transpose() # transpose since rows are rays and we want columns to be rays

# def getRaysFromVrepresentation(vFilePath: str):
#     '''
#     builds a matrix where the columns are the rays that are contained in a V-Representation
#     '''
#     beginOfRays = False
#     listRayMatrix = []
#     with open(vFilePath, "r") as file:
#         for line in file.readlines():
#             if "end" in line:
#                 break
#             if "begin" in line:
#                 beginOfRays = True
#                 continue
#             if beginOfRays:
#                 if "*" in line:
#                     continue
#                 matches = re.findall("-?\d+", line) # get all numbers from line
#                 if matches[0] != "0":
#                     continue
#                 row = []
#                 matches.pop(0) # remove first number that indicates if line is vertex or ray
#                 for match in matches:
#                     row.append(match)
#                 # print(row)
#                 listRayMatrix.append(row)
#         file.close()

#     # print(np.matrix(listRayMatrix, dtype="float64")) 
#     print(listRayMatrix)
#     print(len(listRayMatrix))
#     rayMatrix = np.array(listRayMatrix, dtype="float64") #.transpose() # transpose since rows are rays and we want columns to be rays
#     #print(rayMatrix)
#     return rayMatrix

def getMatrixFromHrepresentation(hFilePath: str):
    '''
    Builds a matrix from an H-Representation, properly handling rational numbers (numerator/denominator).
    '''
    beginOfRays = False
    rayMatrix = []

    with open(hFilePath, "r") as file:
        for line in file.readlines():
            if "end" in line:
                break
            if "begin" in line:
                beginOfRays = True
                continue
            if beginOfRays:
                if "*" in line:
                    continue

                # Match rational numbers in the form of numerator/denominator or integers
                matches = re.findall(r"-?\d+(?:/\d+)?", line)

                if not matches or matches[0] != "0":
                    continue

                row = []
                matches.pop(0)  # Remove the first number (rhs indicator)

                for match in matches:
                    # Convert each match to a Fraction for exact arithmetic
                    row.append(float(Fraction(match)))

                rayMatrix.append(row)

    # Convert the ray matrix to a NumPy matrix (float64 for compatibility with numerical computations)
    return np.array(rayMatrix, dtype="float64")

def getRaysFromVrepresentation(vFilePath: str):
    '''
    builds a matrix where the columns are the rays that are contained in a V-Representation
    '''
    beginOfRays = False
    listRayMatrix = []
    with open(vFilePath, "r") as file:
        for line in file.readlines():
            if "end" in line:
                break
            if "begin" in line:
                beginOfRays = True
                continue
            if beginOfRays:
                if "*" in line:
                    continue
                # Updated regex to capture rational numbers and integers
                matches = re.findall(r"-?\d+(?:/\d+)?", line)
                if not matches or matches[0] != "0":
                    continue
                row = []
                matches.pop(0)  # Remove the first number (indicates if line is vertex or ray)
                for match in matches:
                    # Convert to Fraction and append
                    row.append(float(Fraction(match)))
                listRayMatrix.append(row)
        file.close()

    # Convert to NumPy array
    rayMatrix = np.array(listRayMatrix, dtype="float64")  # Transpose if needed
    return rayMatrix


# def convertEqualitiesAndInequalities2hRep(eqMatrix: np.matrix, iqMatrix: np.matrix, hFilePath: str):
#     '''
#     converts an equality matrix and an inequality matrix to the corresponding H-representation. Output is written to ${hFilePath}
#     '''
#     if os.path.exists(hFilePath):
#         os.remove(hFilePath)

#     rows = eqMatrix.shape[0] + iqMatrix.shape[0]
#     #rows = eqMatrix.shape[0] * 2 + iqMatrix.shape[0]
#     if eqMatrix.shape[1] != iqMatrix.shape[1]:
#         raise Exception("Matrices must have the same column dimensions")
    
#     with open(hFilePath, "w") as file:
#         file.write("poly\n")
#         file.write("H-representation\n")
#         linearities = np.array(range(eqMatrix.shape[0]))+1
#         file.write(f"linearity {eqMatrix.shape[0]} {' '.join(map(str, linearities))}\n")
#         file.write("begin\n")
#         line = str(rows) + " " + str((eqMatrix.shape[1]+1)) + " rational\n" 
#         file.write(line)
#         for i in range(eqMatrix.shape[0]): # convert equality matrix
#             line = "0"
#             for j in range(eqMatrix.shape[1]):
#                 value = eqMatrix[i,j]
#                 line += " " 
#                 if value == int(value): 
#                     insertValue = int(value) 
#                 else:
#                     insertValue = Fraction(float(value)).limit_denominator()
#                 line += str(insertValue)
#             line += "\n"
#             file.write(line)
            
#         for i in range(iqMatrix.shape[0]): # convert inequality matrix
#             line = "0"
#             for j in range(iqMatrix.shape[1]):
#                 value = iqMatrix[i,j]
#                 line += " " 
#                 if value == int(value): 
#                     insertValue = int(value) 
#                 else:
#                     insertValue = Fraction(float(value)).limit_denominator()
#                 line += str(insertValue)
#             line += "\n"
#             file.write(line)

#         file.write("end\n")
#         file.close()

# def convertEqualities2hRep(eqMatrix: np.matrix, hFilePath: str):
#     '''
#     converts an equality matrix to the corresponding H-representation. Output is written to ${hFilePath}
#     '''
#     if os.path.exists(hFilePath):
#         os.remove(hFilePath)

#     rows = eqMatrix.shape[0]
    
#     with open(hFilePath, "w") as file:
#         file.write("poly\n")
#         file.write("H-representation\n")
#         linearities = np.array(range(eqMatrix.shape[0]))+1
#         file.write(f"linearity {eqMatrix.shape[0]} {' '.join(map(str, linearities))}\n")
#         file.write("begin\n")
#         line = str(rows) + " " + str((eqMatrix.shape[1]+1)) + " rational\n" 
#         file.write(line)
#         for i in range(eqMatrix.shape[0]): # convert equality matrix
#             line = "0"
#             for j in range(eqMatrix.shape[1]):
#                 value = eqMatrix[i,j]
#                 line += " " 
#                 if value == int(value): 
#                     insertValue = int(value) 
#                 else:
#                     insertValue = Fraction(float(value)).limit_denominator()
#                 line += str(insertValue)
#             line += "\n"
#             file.write(line)
    

#         file.write("end\n")
#         file.close()
    
def convertEqualitiesAndInequalities2hRep(eqMatrix: np.matrix, iqMatrix: np.matrix, hFilePath: str):
    '''
    converts an equality matrix and an inequality matrix to the corresponding H-representation. Output is written to ${hFilePath}
    '''
    if os.path.exists(hFilePath):
        os.remove(hFilePath)

    rows = eqMatrix.shape[0] * 2 + iqMatrix.shape[0]
    if eqMatrix.shape[1] != iqMatrix.shape[1]:
        raise Exception("Matrices must have the same column dimensions")
    
    with open(hFilePath, "w") as file:
        file.write("poly\n")
        file.write("H-representation\n")
        file.write("begin\n")
        line = str(rows) + " " + str((eqMatrix.shape[1]+1)) + " rational\n" 
        file.write(line)
        for i in range(eqMatrix.shape[0]): # convert equality matrix
            line = "0"
            for j in range(eqMatrix.shape[1]):
                value = eqMatrix[i,j]
                line += " " 
                if value == int(value): 
                    insertValue = int(value) 
                else:
                    insertValue = as_fraction(value)#Fraction(float(value)).limit_denominator()
                line += str(insertValue)
            line += "\n"
            file.write(line)

            line = "0"
            for j in range(eqMatrix.shape[1]):
                value = -eqMatrix[i,j]
                line += " " 
                if value == int(value): 
                    insertValue = int(value) 
                else:
                    insertValue = as_fraction(value)#Fraction(float(value)).limit_denominator()
                line += str(insertValue)
            line += "\n"
            file.write(line)
        
        for i in range(iqMatrix.shape[0]): # convert inequality matrix
            line = "0"
            for j in range(iqMatrix.shape[1]):
                value = iqMatrix[i,j]
                line += " " 
                if value == int(value): 
                    insertValue = int(value) 
                else:
                    insertValue = as_fraction(value)#Fraction(float(value)).limit_denominator()
                line += str(insertValue)
            line += "\n"
            file.write(line)

        file.write("end\n")
        file.close()

def convertEqualities2hRep(eqMatrix: np.matrix, hFilePath: str):
    '''
    converts an equality matrix to the corresponding H-representation. Output is written to ${hFilePath}
    '''
    if os.path.exists(hFilePath):
        os.remove(hFilePath)

    rows = eqMatrix.shape[0] * 2
    
    with open(hFilePath, "w") as file:
        file.write("poly\n")
        file.write("H-representation\n")
        file.write("begin\n")
        line = str(rows) + " " + str((eqMatrix.shape[1]+1)) + " rational\n" 
        file.write(line)
        for i in range(eqMatrix.shape[0]): # convert equality matrix
            line = "0"
            for j in range(eqMatrix.shape[1]):
                value = eqMatrix[i,j]
                line += " " 
                if value == int(value): 
                    insertValue = int(value) 
                else:
                    insertValue = as_fraction(value)#Fraction(float(value)).limit_denominator()
                line += str(insertValue)
            line += "\n"
            file.write(line)

            line = "0"
            for j in range(eqMatrix.shape[1]):
                value = -eqMatrix[i,j]
                line += " " 
                if value == int(value): 
                    insertValue = int(value) 
                else:
                    insertValue = as_fraction(value)#Fraction(float(value)).limit_denominator()
                line += str(insertValue)
            line += "\n"
            file.write(line)
    

        file.write("end\n")
        file.close()