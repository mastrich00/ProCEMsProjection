import numpy as np
import subprocess
import time
from datetime import datetime
from projection.helpers import printDatetime
from projection.marashiProjectionHelpers import sortStoicMatrix, buildInterestedMatrix, buildEliminationMatrix, buildProjectedConeMatrix, \
    calculateEFPsFromProCEMs, logTimeToCSV
from projection.enumeration import doubleDescription, doubleDescriptionINE, convertMatrixToHRepresentation, convertMatrixToVRepresentation, getMatrixFromHrepresentation,convertEqualitiesAndInequalities2hRep, \
    getRaysFromVrepresentation, convertHtoVrepresentation, mplrs_conversion, convertEqualities2hRep
from projection.redund import redund, redundIterative, normalize_matrix_by_gcd_and_scale, process_matrix
from projection.logging_config import logger
import os
from fractions import Fraction
if not os.getenv("POLCO_PATH"):
    load_dotenv("env") # load env
    
POLCO_PATH = os.getenv("POLCO_PATH")
mplrsPath = MPLRS_PATH = os.getenv("MPLRS_PATH")
mplrsPath = MPLRS_PATH = "mplrsV7_2"
import os
import time

def round_sig(x, sig=15):
    """
    Round a number to 'sig' significant figures.
    Returns 0.0 if x is zero.
    """
    if x == 0:
        return 0.0
    # Compute the number of digits after the decimal point required
    n_digits = int(sig - np.floor(np.log10(abs(x))) - 1)
    return round(x, n_digits)

def scale_and_round_row_np(row, sig=15):
    """
    Given a row (list or NumPy array) representing an inequality,
    scale it so that the maximum absolute value is about 10^(sig)
    (i.e. roughly sig digits) and then round each entry to sig
    significant figures.
    """
    # Convert the row to a NumPy array (float64 for precision)
    row = np.array(row, dtype=np.float64)
    
    # Identify nonzero entries (to avoid log10(0))
    nonzero = row != 0
    if not np.any(nonzero):
        return row  # if all entries are zero, no scaling is needed.
    
    # Find the maximum absolute nonzero value
    max_val = np.max(np.abs(row[nonzero]))
    
    # Compute scaling factor:
    # We want to bring the order of magnitude of max_val to about 10^sig.
    factor = 10 ** (np.ceil(np.log10(max_val)) - sig)
    
    # Scale the row by the factor
    scaled_row = row / factor
    
    # Vectorize the round_sig function to apply it elementwise
    vec_round_sig = np.vectorize(round_sig)
    rounded_row = vec_round_sig(scaled_row, sig)
    
    return rounded_row

def truncate_float(x, digits):
    factor = 10 ** digits
    return np.round(x * factor) / factor
                
def normalize_row_if_needed(row, digit_threshold=12, trunc_digits=3):
    """
    Normalizes the row so that its elements sum to 1 only if at least one element 
    in the row has more than norm_threshold digits in its numerator or denominator.
    
    After normalization, if any normalized value has a numerator or denominator with 
    more than digit_threshold digits, then each value is reduced as follows: its 
    floating-point value is truncated to trunc_digits digits after the decimal point 
    (i.e. any extra decimals are cut off) and then converted back to a Fraction.
    """
    # Decide if normalization is needed based on the original fractions.
    needs_normalization = any(
        len(str(abs(x.numerator))) > digit_threshold or len(str(abs(x.denominator))) > digit_threshold
        for x in row
    )
    
    if needs_normalization:            
            # Truncate each normalized element to trunc_digits decimal places,
            # then convert back to a Fraction.
            normalized_truncated = [
                Fraction(truncate_float(float(val), trunc_digits))
                for val in row
            ]
            return np.array(normalized_truncated, dtype=object)
    else:
        # If no element requires normalization, return the row unchanged.
        return row
            
def get_sorted_column_indices_from_array(arr, lenOriginalReactions):
    """
    Given a 2D numpy array, ignores the first 24 columns and returns a list of 
    original column indices (of the remaining columns) ordered by descending sum 
    of the absolute values of their entries.

    Parameters:
        arr (np.array): A 2D numpy array.

    Returns:
        List[int]: Sorted original column indices (ignoring the first 24 columns) 
                   ordered by descending sum of absolute values.
    """
    # Verify that the array has more than 24 columns.
    if arr.shape[1] <= lenOriginalReactions:
        raise ValueError(f"The array must have more than {lenOriginalReactions} columns.")
    
    # Select columns after the first 24.
    remaining_columns = arr[:, lenOriginalReactions:]
    
    # Compute the sum of absolute values for each of the remaining columns.
    abs_sums = np.sum(np.abs(remaining_columns), axis=0)
    
    # Get the indices that would sort these sums in descending order.
    # np.argsort sorts in ascending order by default, so we sort the negative values to reverse the order.
    sorted_order = np.argsort(-abs_sums)
    
    # Adjust the indices to match the original array by adding the offset (24).
    original_indices = sorted_order + lenOriginalReactions
    
    return original_indices.tolist()

def runMarashiWithMPLRSSubsets(stoicMatrix, projectOntoDimensions, outputDir, numThreads=8, stopAfterProjection=True, verbose = True, iteration=0, originalProjectionReactions=[], logTimesFile="times.csv", reversibleList=[]):
    time_setup_stoic_start = time.time()
    printDatetime("Start:", datetime.now())
    sortedStoicMatrix = stoicMatrix
    convertMatrixToHRepresentation(sortedStoicMatrix, os.path.join(outputDir, "stoicMatrix.ine"))
    if iteration == 0:
        p = len(originalProjectionReactions)
        q = sortedStoicMatrix.shape[1] - p
        convertEqualitiesAndInequalities2hRep(sortedStoicMatrix, -np.identity(p+q)[reversibleList],
            os.path.join(outputDir, f"{iteration}_stoicMatrix.ine")
        )
        # convertEqualitiesAndInequalities2hRep(sortedStoicMatrix, 
        #     np.hstack([
        #         np.vstack([
        #             -np.identity(p),
        #             np.zeros((q,p))]),
        #         np.vstack([
        #             np.zeros((p,q)),
        #             -np.identity(q)]),
        #     ]),
        #     os.path.join(outputDir, f"{iteration}_stoicMatrix.ine")
        # )
        # convertEqualitiesAndInequalities2hRep(stoicMatrix[:,:len(originalProjectionReactions)], np.identity(stoicMatrix[:,:len(originalProjectionReactions)].shape[1]), os.path.join(outputDir, f"{iteration}_stoicMatrixInterested.ine"))
        # redund(numThreads,mplrsPath, os.path.join(outputDir, f"{iteration}_stoicMatrixElimination.ine"),os.path.join(outputDir, f"{iteration}_redund_stoicMatrixElimination.ine"))
        redund(numThreads,mplrsPath, os.path.join(outputDir, f"{iteration}_stoicMatrix.ine"),os.path.join(outputDir, f"{iteration}_stoicMatrix_redund.ine"))
        sortedStoicMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, f"{iteration}_stoicMatrix_redund.ine"))
        # if verbose:
        #     logger.info("InterestedMatrix:")
        #     logger.info(interestedMatrix)
        # colsElimination = eliminationMatrix.shape[1] - (sortedStoicMatrix.shape[1] - len(projectOntoDimensions))
        # interestedMatrix = np.hstack([interestedMatrix, eliminationMatrix[:,:colsElimination]])
        # eliminationMatrix = eliminationMatrix[:,colsElimination:]
        interestedMatrix = sortedStoicMatrix[:,:len(projectOntoDimensions)]
        eliminationMatrix = sortedStoicMatrix[:,len(projectOntoDimensions):]
    else:
        interestedMatrix = sortedStoicMatrix[:, :len(projectOntoDimensions)]
        eliminationMatrix = sortedStoicMatrix[:, len(projectOntoDimensions):]
    eliminationMatrix = eliminationMatrix.transpose()
    print(f"shape el: {eliminationMatrix.shape}")
    print(f"shape int: {interestedMatrix.shape}")
    time_setup_stoic_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Setup Elimination-Matrix", time_setup_stoic_end - time_setup_stoic_start)
    
    time_projectioncone_redund_start = time.time()
    hEliminationInePath = os.path.join(outputDir, "elimination_H.ine")
    hEliminationRedundInePath = os.path.join(outputDir, "elimination_redund_H.ine")
    vEliminationInePath = os.path.join(outputDir, "elimination_V.ine")
    vEliminationRedundInePath = os.path.join(outputDir, "elimination_redund_V.ine")
    convertEqualitiesAndInequalities2hRep(eliminationMatrix, np.identity(eliminationMatrix.shape[1]), hEliminationInePath) # build W as ine file
    printDatetime("Finished creating eliminiation_H.ine:", datetime.now())
    redund(numThreads,mplrsPath,hEliminationInePath,hEliminationRedundInePath) # remove redundancy
    time_projectioncone_redund_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Redund Projection-Cone", time_projectioncone_redund_end - time_projectioncone_redund_start)

    time_projectioncone_conversion_start = time.time()
    mplrs_conversion(numThreads, mplrsPath, hEliminationRedundInePath, vEliminationInePath) # convert to v-representation
    # mplrs_conversion(numThreads, mplrsPath, hEliminationInePath, vEliminationInePath) # convert to v-representation
    time_projectioncone_conversion_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Enumeration Projection-Cone", time_projectioncone_conversion_end - time_projectioncone_conversion_start)

    time_rays_redund_start = time.time()
    printDatetime("Finished mplrs conversion:", datetime.now())    
    rayMatrix = getRaysFromVrepresentation(vEliminationInePath)
    time_rays_redund_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Redund Rays", time_rays_redund_end - time_rays_redund_start)
    logger.info(f"Raymatrix - Shape: {rayMatrix.shape}")
    if verbose:
        logger.info("RayMatrix:")
        logger.info(rayMatrix)

    if verbose:
        logger.info("InterestedMatrix:")
        logger.info(interestedMatrix)
        
    time_build_projectedcone_start = time.time()
    projectedConeMatrix = buildProjectedConeMatrix(rayMatrix, interestedMatrix)
    logger.info(f"projectedConeMatrix - Shape: {projectedConeMatrix.shape}")
    time_build_projectedcone_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Build projected Cone", time_build_projectedcone_end - time_build_projectedcone_start)

    time_redund_projectedcone_start = time.time()
    hProjectedConeInePath = os.path.join(outputDir, f"{iteration}_projectedCone_H.ine")
    #projectedConeMatrix = normalize_matrix_by_gcd_and_scale(projectedConeMatrix)
    # projectedConeMatrix = np.array([normalize_row_if_needed(row) for row in projectedConeMatrix], dtype=object)
    convertMatrixToHRepresentation(projectedConeMatrix, hProjectedConeInePath)
    printDatetime("Finished converting projected Cone Matrix to H-Rep:", datetime.now())
    if projectedConeMatrix.shape[0] > 300000:
        print(f"Error: shape={projectedConeMatrix.shape}")
        exit(0)

    hRedundProjectedConeInePath = os.path.join(outputDir, "projectedConeMatrix_H_redund.ine")
    if not os.path.exists(os.path.join(outputDir,"temp")):
        os.mkdir(os.path.join(outputDir,"temp"))
    # chunkSize = 10000
    #logger.info(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} with chunkSize {chunkSize}")
    #redundIterative(numThreads, mplrsPath, hProjectedConeInePath, os.path.join(outputDir, "tempOutprojectedCone_H.ine"), os.path.join(outputDir,"temp"), chunkSize=chunkSize, verbose=True)
    #logger.info(f"End: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    redund(numThreads, mplrsPath, hProjectedConeInePath, hRedundProjectedConeInePath)
    printDatetime("Finished redund step:", datetime.now())
    matrix = getMatrixFromHrepresentation(hRedundProjectedConeInePath)
    # matrix = getMatrixFromHrepresentation(hProjectedConeInePath)
    time_redund_projectedcone_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Redund projected Cone", time_redund_projectedcone_end - time_redund_projectedcone_start)
    if verbose:
        logger.info("projectedConeMatrix:")
        logger.info(projectedConeMatrix)
        logger.info(f"Shape: {projectedConeMatrix.shape}")
    if stopAfterProjection:
        return matrix, None

    time_enumerate_projectedcone_start = time.time()
    vProjectedConeInePath = os.path.join(outputDir, "projectedConeMatrix_V.ine")
    #!!! mplrs_conversion(numThreads,mplrsPath, hProjectedConeInePath, vProjectedConeInePath)
    mplrs_conversion(numThreads,mplrsPath, hRedundProjectedConeInePath, vProjectedConeInePath)
    # convertHtoVrepresentation(hRedundProjectedConeInePath, vProjectedConeInePath)
    printDatetime("Finished H-Rep to V-Rep conversion:", datetime.now())
    proCEMs = getRaysFromVrepresentation(vProjectedConeInePath)
    time_enumerate_projectedcone_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Enumerate projected Cone", time_enumerate_projectedcone_end - time_enumerate_projectedcone_start)
    if verbose:
        print("\n#########\nBefore Redund - proCEMs:\n", proCEMs)
        print("\n#########\nResulting proCEMs-Shape:\n", proCEMs.shape)

    time_redund_procems_start = time.time()
    vProCEMsPath = os.path.join(outputDir,"proCEMs_V.ine")
    redundvProCEMsPath = os.path.join(outputDir,"redund_proCEMs_V.ine")
    convertMatrixToVRepresentation(proCEMs, vProCEMsPath)
    if proCEMs.shape[0] > 300000:
        print(f"Error: shape={projectedConeMatrix.shape}")
        exit(0)
    redund(numThreads,mplrsPath,vProCEMsPath,redundvProCEMsPath)
    printDatetime("Finished redund step for proCEMs:", datetime.now())
    proCEMs = getRaysFromVrepresentation(redundvProCEMsPath)
    time_redund_procems_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Redund proCEMs", time_redund_procems_end - time_redund_procems_start)
    if verbose:
        print("\n#########\nAfter Redund - proCEMs:\n", proCEMs)
        print("\n#########\nResulting proCEMs-Shape:\n", proCEMs.shape)

    # time_calc_efps_start = time.time()
    # efps = calculateEFPsFromProCEMs(proCEMs)
    # time_calc_efps_end = time.time()
    # logTimeToCSV(logTimesFile, f"Iter {iteration}", "Calc EFPs", time_calc_efps_end - time_calc_efps_start)
    # printDatetime("Finished calculating efps:", datetime.now())
    #if verbose:
    #    print("EFPs:", efps)
    efps = []
    return proCEMs, efps

def runMarashiWithPolcoSubsets(stoicMatrix, projectOntoDimensions, outputDir, numThreads=8, stopAfterProjection=True, verbose = True, iteration=0, originalProjectionReactions=[], logTimesFile="times.csv", reversibleList=[]):
    time_setup_stoic_start = time.time()
    printDatetime("Start: ", datetime.now())
    #sortedStoicMatrix = sortStoicMatrix(stoicMatrix, projectOntoDimensions)
    sortedStoicMatrix = stoicMatrix
    # sortReactions = get_sorted_column_indices_from_array(sortedStoicMatrix, len(originalProjectionReactions))
    # sortReactions = list(range(len(originalProjectionReactions))) + sortReactions
    # #print(sortReactions)
    # sortedStoicMatrix = sortedStoicMatrix[:,sortReactions]
    np.savetxt(os.path.join(outputDir, 'array.csv'), sortedStoicMatrix, delimiter=',', fmt='%d')

    convertMatrixToHRepresentation(sortedStoicMatrix, os.path.join(outputDir, "stoicMatrix.ine"))
    if iteration == 0:
        p = len(originalProjectionReactions)
        q = sortedStoicMatrix.shape[1] - p
        convertEqualitiesAndInequalities2hRep(sortedStoicMatrix, -np.identity(p+q)[reversibleList],
            os.path.join(outputDir, f"{iteration}_stoicMatrix.ine")
        )
        # convertEqualitiesAndInequalities2hRep(sortedStoicMatrix, 
        #     np.hstack([
        #         np.vstack([
        #             -np.identity(p),
        #             np.zeros((q,p))]),
        #         np.vstack([
        #             np.zeros((p,q)),
        #             -np.identity(q)]),
        #     ]),
        #     os.path.join(outputDir, f"{iteration}_stoicMatrix.ine")
        # )
        # convertEqualitiesAndInequalities2hRep(stoicMatrix[:,:len(originalProjectionReactions)], np.identity(stoicMatrix[:,:len(originalProjectionReactions)].shape[1]), os.path.join(outputDir, f"{iteration}_stoicMatrixInterested.ine"))
        # redund(numThreads,mplrsPath, os.path.join(outputDir, f"{iteration}_stoicMatrixElimination.ine"),os.path.join(outputDir, f"{iteration}_redund_stoicMatrixElimination.ine"))
        redund(numThreads,mplrsPath, os.path.join(outputDir, f"{iteration}_stoicMatrix.ine"),os.path.join(outputDir, f"{iteration}_stoicMatrix_redund.ine"))
        sortedStoicMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, f"{iteration}_stoicMatrix_redund.ine"))
        # if verbose:
        #     logger.info("InterestedMatrix:")
        #     logger.info(interestedMatrix)
        # colsElimination = eliminationMatrix.shape[1] - (sortedStoicMatrix.shape[1] - len(projectOntoDimensions))
        # interestedMatrix = np.hstack([interestedMatrix, eliminationMatrix[:,:colsElimination]])
        # eliminationMatrix = eliminationMatrix[:,colsElimination:]
        interestedMatrix = sortedStoicMatrix[:,:len(projectOntoDimensions)]
        eliminationMatrix = sortedStoicMatrix[:,len(projectOntoDimensions):]
    else:
        interestedMatrix = sortedStoicMatrix[:, :len(projectOntoDimensions)]
        eliminationMatrix = sortedStoicMatrix[:, len(projectOntoDimensions):]
    print(f"shape el: {eliminationMatrix.shape}")
    print(f"shape int: {interestedMatrix.shape}")
    #convertMatrixToHRepresentation(-eliminationMatrix.transpose(), os.path.join(outputDir, "tempRedundSort.ine"))
    #redund(numThreads,mplrsPath,os.path.join(outputDir, "tempRedundSort.ine"),os.path.join(outputDir, "outTempRedundSort.ine")) # remove redundancy
    #eliminationMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, "outTempRedundSort.ine"))
    eliminationMatrix = eliminationMatrix.transpose()
    time_setup_stoic_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Setup Elimination-Matrix", time_setup_stoic_end - time_setup_stoic_start)

    time_projectioncone_conversion_start = time.time()
    eqFile = 'file.eq'
    iqFileIdentity = 'identity.iq'
    stage = 'projectionCone'
    outputFile = stage + '_out.zip'
    if verbose:
        logger.info("Start Double Description:")
    # enumerate rays of W:
    #rayMatrix = doubleDescription(-np.vstack([eliminationMatrix, np.identity(eliminationMatrix.shape[1])]), iqFileIdentity, outputDir, outputFile, stage, numThreads=numThreads)
    #!!! rayMatrix = doubleDescription(np.identity(eliminationMatrix.shape[1]), iqFileIdentity, outputDir, outputFile, stage, eliminationMatrix, eqFile, numThreads=numThreads)
    convertEqualitiesAndInequalities2hRep(eliminationMatrix, np.identity(eliminationMatrix.shape[1]), os.path.join(outputDir, "doubleDescription.ine"))
    # convertEqualities2hRep(eliminationMatrix, os.path.join(outputDir, "doubleDescription.ine"))
    #redund(numThreads,mplrsPath,os.path.join(outputDir, "doubleDescription.ine"),os.path.join(outputDir, "outDoubleDescription.ine")) # remove redundancy
    #exit(0)
    rayMatrix = doubleDescriptionINE(os.path.join(outputDir, "doubleDescription.ine"), outputDir, outputFile, stage, numThreads=numThreads)
    time_projectioncone_conversion_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Enumeration Projection-Cone", time_projectioncone_conversion_end - time_projectioncone_conversion_start)
    
    # print(f"Before redund - elimination rays: {rayMatrix.shape}")
    # time_rays_redund_start = time.time()
    # elimination_rays_V_path = os.path.join(outputDir, "elimination_rays_V.ine")
    # elimination_rays_V_Redund_path = os.path.join(outputDir, "elimination_rays_V_redund.ine")
    # convertMatrixToVRepresentation(rayMatrix, elimination_rays_V_path)
    # redund(numThreads,mplrsPath,elimination_rays_V_path,elimination_rays_V_Redund_path)
    # rayMatrix = getRaysFromVrepresentation(elimination_rays_V_Redund_path)
    # print(f"After redund - elimination rays: {rayMatrix.shape}")
    logger.info(f"Raymatrix - Shape: {rayMatrix.shape}")
    printDatetime("Finished first double description step:", datetime.now())
    # rayMatrix = rayMatrix.T
    if verbose:
        logger.info("RayMatrix:")
        logger.info(rayMatrix)
    # rayMatrix[:,[1,3]] = rayMatrix[:,[3,1]]
    # print("RayMatrix:\n", rayMatrix)
    # time_rays_redund_end = time.time()
    # logTimeToCSV(logTimesFile, f"Iter {iteration}", "Redund Rays", time_rays_redund_end - time_rays_redund_start)
    
    # projectedConeMatrix = buildProjectedConeMatrix(rayMatrix.transpose(), interestedMatrix)
    time_build_projectedcone_start = time.time()
    printDatetime("Calculating projected cone matrix:", datetime.now())
    projectedConeMatrix = buildProjectedConeMatrix(rayMatrix, interestedMatrix) # build condition of projected cone
    printDatetime("Finished calculating projected cone matrix:", datetime.now())
    time_build_projectedcone_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Build projected Cone", time_build_projectedcone_end - time_build_projectedcone_start)

    #projectedConeMatrix = process_matrix(projectedConeMatrix)
    #if iteration == 0:
    #    projectedConeMatrix *= -1
    logger.info(f"projectedConeMatrix - Shape: {projectedConeMatrix.shape}")
    # projectedConeMatrix = removeRowsWithZeros(projectedConeMatrix)

    time_redund_projectedcone_start = time.time()
    ineFile = os.path.join(outputDir, "projectedCone_H.ine")
    
    # projectedConeMatrix = normalize_matrix_by_gcd_and_scale(projectedConeMatrix)
    #projectedConeMatrix = np.array([normalize_matrix_by_gcd_and_scale(row) for row in projectedConeMatrix])    
    # Apply conditional normalization to each row of the matrix.
    # projectedConeMatrix = np.array([normalize_row_if_needed(row) for row in projectedConeMatrix], dtype=object
    # scale_and_round_row_np
    # projectedConeMatrix = np.array([scale_and_round_row_np(row) for row in projectedConeMatrix], dtype=object)
    # print(projectedConeMatrix)
    
    convertMatrixToHRepresentation(projectedConeMatrix, ineFile)
    printDatetime("Finished converting projected Cone Matrix to H-Rep:", datetime.now())
    if projectedConeMatrix.shape[0] > 300000:
        print(f"Error: shape={projectedConeMatrix.shape}")
        exit(0)
    #projectedConeMatrix[[2,3],:]= projectedConeMatrix[[3,2],:]
    #projectedConeMatrix = removeRedundancy(projectedConeMatrix)
    if not os.path.exists(os.path.join(outputDir, "temp")):
        os.mkdir(os.path.join(outputDir, "temp"))
    #chunkSize = 200
    #logger.info(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} with chunkSize {chunkSize}")
    #redundIterative(numThreads, mplrsPath, ineFile, os.path.join(outputDir, "tempFirstOutprojectedCone_H.ine"), os.path.join(outputDir, "temp"), chunkSize=chunkSize, verbose=True)
    #chunkSize = 600
    #redundIterative(numThreads, mplrsPath, os.path.join(outputDir, "tempFirstOutprojectedCone_H.ine"), os.path.join(outputDir, "tempOutprojectedCone_H.ine"), os.path.join(outputDir, "temp"), chunkSize=chunkSize, verbose=True)
    #logger.info(f"End: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    #redund(numThreads,mplrsPath,os.path.join(outputDir, "tempOutprojectedCone_H.ine"),os.path.join(outputDir, f"{iteration}_projectedConeMatrix_final_H.ine")) # remove redundancy
    redund(numThreads,mplrsPath,ineFile,os.path.join(outputDir, f"{iteration}_projectedConeMatrix_final_H.ine")) # remove redundancy
    printDatetime("Finished redund step:", datetime.now())

    projectedConeMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, f"{iteration}_projectedConeMatrix_final_H.ine"))
    time_redund_projectedcone_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Redund projected Cone", time_redund_projectedcone_end - time_redund_projectedcone_start)
    
    #projectedConeMatrix = getMatrixFromHrepresentation(ineFile)
    if verbose:
        logger.info("projectedConeMatrix:")
        logger.info(projectedConeMatrix)
        logger.info(f"Shape: {projectedConeMatrix.shape}")
        
    if stopAfterProjection:
        # projectedConeMatrix = np.vstack([projectedConeMatrix, -projectedConeMatrix])
        return projectedConeMatrix, None
    time_enumerate_projectedcone_start = time.time()
    eqFileIdentity = 'identity.eq'
    iqFile = 'file.iq'
    stage = 'proCEMs'
    outputFile = stage + '_out.zip'
    if verbose:
        printDatetime("Start Double Description:", datetime.now())
    proCEMs = doubleDescription(projectedConeMatrix, iqFile, outputDir, outputFile, stage, numThreads=numThreads)
    time_enumerate_projectedcone_end= time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Enumerate projected Cone", time_enumerate_projectedcone_end - time_enumerate_projectedcone_start)
    printDatetime("Finished second double description step:", datetime.now())
    # if verbose:
    #     print("\n#########\nBefore redund - proCEMs:\n", proCEMs)
    #     print("Shape: ", proCEMs.shape)
    time_redund_procems_start = time.time()
    convertMatrixToVRepresentation(proCEMs, os.path.join(outputDir, "redund_proCEMs_V.ine"))
    if proCEMs.shape[0] > 300000:
        print(f"Error: shape={proCEMs.shape}")
        exit(0)
    redund(numThreads,mplrsPath,os.path.join(outputDir, "redund_proCEMs_V.ine"),os.path.join(outputDir, "outRedundProCEMs_V.ine"))
    printDatetime("Finished redund step for proCEMs:", datetime.now())
    proCEMs = getRaysFromVrepresentation(os.path.join(outputDir, "outRedundProCEMs_V.ine"))
    time_redund_procems_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Redund proCEMs", time_redund_procems_end - time_redund_procems_start)

    
    ## proCEMs= proCEMs.T
    # if verbose:
    #     print("\n#########\nAfter redund - proCEMs:\n", proCEMs)
    #     print("Shape: ", proCEMs.shape)
    efps = []
    # efps = calculateEFPsFromProCEMs(proCEMs)
    # printDatetime("Finished calculating efps:", datetime.now())
    # if verbose:
    #     # print("EFPs:", efps)
    #     print("EFPs-length:", len(efps))
    return proCEMs, efps

def eliminate(ineFilePath: str, dimensionsToEliminate: list, outFilePath: str, mplrsPath: str, numProc= 4, verbose = True, logEliminations = False):
    """
    Eliminates given dimensions in .ine file. Writes output to 'outFilePath'
    Dimension index starts at 1, not zero.
    """
    if len(dimensionsToEliminate) == 0:
        if verbose:
            logger.info("No dimensions to eliminate!")
        return 0
    
    file = open(ineFilePath, "a")
    cmd = "eliminate " + str(len(dimensionsToEliminate))
    for index in dimensionsToEliminate:
        cmd = cmd + " " + str(index)
    if verbose:
        logger.info(f"Command: {cmd}")
    cmd += "\n"
    file.write(cmd)
    if verbose:
        file.write("verbose")
    file.close()

    tmpFilePath = ineFilePath + "_tmp.ine"
    if logEliminations:
        subprocess.run(["cp",ineFilePath, tmpFilePath + "1"])
    else:
        subprocess.run(["cp",ineFilePath, tmpFilePath])

    for i in range(1, len(dimensionsToEliminate)):
        if verbose:
            logger.info("Running elimination step " + str(i) + "/" + str(len(dimensionsToEliminate)))
        if logEliminations:
            cmd = ["mpirun", "-np", str(numProc), mplrsPath, "-fel", tmpFilePath+str(i), outFilePath+str(i)] # TODO: check if -fel necessary
        else:
            cmd = ["mpirun", "-np", str(numProc), mplrsPath, "-fel", tmpFilePath, outFilePath]
        subprocess.run(cmd)

        if logEliminations:
            subprocess.run(["cp", outFilePath+str(i), tmpFilePath+str(i+1)])
        else:
            subprocess.run(["cp", outFilePath, tmpFilePath])

    if verbose:
        logger.info("Running elimination step " + str(len(dimensionsToEliminate)) + "/" + str(len(dimensionsToEliminate)))
    if logEliminations:
        cmd = ["mpirun", "-np", str(numProc), mplrsPath, "-fel", tmpFilePath+str(len(dimensionsToEliminate)), outFilePath+str(len(dimensionsToEliminate))]
    else:
        cmd = ["mpirun", "-np", str(numProc), mplrsPath, "-fel", tmpFilePath, outFilePath]
    subprocess.run(cmd)
    
    subprocess.run(["cp", outFilePath+str(len(dimensionsToEliminate)), outFilePath])

def runMarashiWithPolcoIterativeV2(stoicMatrix, projectOntoDimensions, outputDir, numThreads=8, stopAfterProjection=True, verbose = True):
    print(f"Number projection dims: {len(projectOntoDimensions)}")
    print(f"shape stoic matrix: {stoicMatrix.shape}")
    printDatetime("Start: ", datetime.now())
    #sortedStoicMatrix = sortStoicMatrix(stoicMatrix, projectOntoDimensions)
    sortedStoicMatrix = stoicMatrix
    print(f"stoicmatrix: {sortedStoicMatrix.shape}")
    interestedMatrix = buildInterestedMatrix(sortedStoicMatrix, projectOntoDimensions) # build G
    if verbose:
        logger.info("InterestedMatrix:")
        logger.info(interestedMatrix)
        print(f"InterestedMatrix: {interestedMatrix.shape}")
        
    eliminationMatrix = buildEliminationMatrix(sortedStoicMatrix, projectOntoDimensions) # build H
    #convertMatrixToHRepresentation(-eliminationMatrix.transpose(), os.path.join(outputDir, "tempRedundSort.ine"))
    #redund(numThreads,mplrsPath,os.path.join(outputDir, "tempRedundSort.ine"),os.path.join(outputDir, "outTempRedundSort.ine")) # remove redundancy
    #eliminationMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, "outTempRedundSort.ine"))
    eliminationMatrix = eliminationMatrix.transpose()
    eqFile = 'file.eq'
    iqFileIdentity = 'identity.iq'
    stage = 'projectionCone'
    outputFile = stage + '_out.zip' 
    if verbose:
        logger.info("Start Double Description:")
    # enumerate rays of W:
    #eliminationMatrixT = eliminationMatrix.transpose()
    eliminationMatrixT = eliminationMatrix
    print(f"elMatrixT: {eliminationMatrixT.shape}")
    identityEliminationMatrixT = np.identity(eliminationMatrix.shape[1])
    i = stoicMatrix.shape[1] - len(projectOntoDimensions)
    stepSize = min(5, i)
    print(f"i={i}, stepSize={stepSize}")
    boolStop = False
    #matrix = np.empty((0, 409), float)
    # matrix = np.empty((0, 336), float)
    #projectedConeMatrix = interestedMatrix
    while i > 0:
        print(f"Iter: {i}")
        if (i - stepSize < stepSize):
            #convertEqualitiesAndInequalities2hRep(eliminationMatrixT[0:i,:], identityEliminationMatrixT, os.path.join(outputDir, "doubleDescription.ine"))
            print(f"elRows: 0:{i}")
            tempEliminationMatrix = eliminationMatrixT[0:i,:]
            boolStop = True
        else:
            print(f"elRows: {i-stepSize}:{i}")
            tempEliminationMatrix = eliminationMatrixT[i-stepSize:i,:]
            #convertEqualitiesAndInequalities2hRep(eliminationMatrixT[i-5:i,:], identityEliminationMatrixT, os.path.join(outputDir, "doubleDescription.ine"))
        #redund(numThreads,mplrsPath,os.path.join(outputDir, "doubleDescription.ine"),os.path.join(outputDir, "outDoubleDescription.ine")) # remove redundancy
        
        #rayMatrix = doubleDescriptionINE(os.path.join(outputDir, "outDoubleDescription.ine"), outputDir, outputFile, stage, numThreads=numThreads)
        rayMatrix = doubleDescription(identityEliminationMatrixT, iqFileIdentity, outputDir, outputFile, stage, tempEliminationMatrix, eqFile, numThreads=numThreads)
        logger.info(f"Raymatrix - Shape: {rayMatrix.shape}")
        print(f"Raymatrix - Shape: {rayMatrix.shape}")
        convertMatrixToVRepresentation(rayMatrix, os.path.join(outputDir, "rayMatrixTest.ine"))
        redund(numThreads,mplrsPath,os.path.join(outputDir, "rayMatrixTest.ine"),os.path.join(outputDir, f"{i}_redundRayMatrixTest.ine")) # remove redundancy
        printDatetime(f"Finished Raymatrix step {i-stepSize}-{i}:", datetime.now())
        matrix = np.vstack([matrix, getRaysFromVrepresentation(os.path.join(outputDir, f"{i}_redundRayMatrixTest.ine"))])
        #matrix = np.vstack([matrix, rayMatrix])
        #projectedConeMatrix = buildProjectedConeMatrix(rayMatrix, projectedConeMatrix)
        i = i - stepSize
        if boolStop:
            break
    printDatetime("Finished first double description step:", datetime.now())
    # rayMatrix = rayMatrix.T
    #logger.info(matrix)
    #print(matrix)
    convertMatrixToVRepresentation(matrix, os.path.join(outputDir, f"matrix.ine"))
    rayMatrix = matrix
    if verbose:
        logger.info("RayMatrix:")
        logger.info(rayMatrix)
        print(f"rayMatrix: {rayMatrix.shape}")
    ## rayMatrix[:,[1,3]] = rayMatrix[:,[3,1]]
    ## print("RayMatrix:\n", rayMatrix)

    # projectedConeMatrix = buildProjectedConeMatrix(rayMatrix.transpose(), interestedMatrix)
    projectedConeMatrix = -buildProjectedConeMatrix(rayMatrix, interestedMatrix) # build condition of projected cone
    #logger.info(f"projectedConeMatrix - Shape: {projectedConeMatrix.shape}")
    # projectedConeMatrix = removeRowsWithZeros(projectedConeMatrix)
    ineFile = os.path.join(outputDir, "projectedCone_H.ine")
    convertMatrixToHRepresentation(projectedConeMatrix, ineFile)
    printDatetime("Finished converting projected Cone Matrix to H-Rep:", datetime.now())

    #projectedConeMatrix[[2,3],:]= projectedConeMatrix[[3,2],:]
    #projectedConeMatrix = removeRedundancy(projectedConeMatrix)
    if not os.path.exists(os.path.join(outputDir, "temp")):
        os.mkdir(os.path.join(outputDir, "temp"))
    chunkSize = 10000
    #logger.info(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} with chunkSize {chunkSize}")
    #redundIterative(numThreads, mplrsPath, ineFile, os.path.join(outputDir, "tempOutprojectedCone_H.ine"), os.path.join(outputDir, "temp"), chunkSize=chunkSize, verbose=True)
    #logger.info(f"End: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    redund(numThreads,mplrsPath,ineFile,os.path.join(outputDir, "projectedConeMatrix_final_H.ine")) # remove redundancy
    printDatetime("Finished redund step:", datetime.now())

    projectedConeMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, "projectedConeMatrix_final_H.ine"))
    if verbose:
        logger.info("projectedConeMatrix:")
        logger.info(projectedConeMatrix)
        logger.info(f"Shape: {projectedConeMatrix.shape}")
        
    if stopAfterProjection:
        return None, None

    eqFileIdentity = 'identity.eq'
    iqFile = 'file.iq'
    stage = 'proCEMs'
    outputFile = stage + '_out.zip'
    if verbose:
        printDatetime("Start Double Description:", datetime.now())

    
    #proCEMs = doubleDescriptionINE(os.path.join(outputDir, "projectedConeMatrix_final_H.ine"), outputDir, outputFile, stage, numThreads=numThreads)
    #logger.info(f"proCEMs - Shape: {proCEMs.shape}")
    #print(f"proCEMs - Shape: {proCEMs.shape}")
    print(projectedConeMatrix)
    proCEMs = doubleDescription(projectedConeMatrix, iqFile, outputDir, outputFile, stage, numThreads=numThreads)
    logger.info(f"proCEMs - Shape: {proCEMs.shape}")
    print(f"proCEMs - Shape: {proCEMs.shape}")
    printDatetime("Finished second double description step:", datetime.now())
    # proCEMs= proCEMs.T
    if verbose:
        print("\n#########\nResulting proCEMs:\n", proCEMs)
    efps = calculateEFPsFromProCEMs(proCEMs)
    printDatetime("Finished calculating efps:", datetime.now())
    if verbose:
        print("EFPs:", efps)
        print("EFPs-length:", len(efps))
    return proCEMs, efps
    
def runMarashiWithPolcoIterative(stoicMatrix, projectOntoDimensions, outputDir, numThreads=8, stopAfterProjection=True, verbose = True):
    print(f"Number projection dims: {len(projectOntoDimensions)}")
    print(f"shape stoic matrix: {stoicMatrix.shape}")
    printDatetime("Start: ", datetime.now())
    #sortedStoicMatrix = sortStoicMatrix(stoicMatrix, projectOntoDimensions)
    sortedStoicMatrix = stoicMatrix
    print(f"stoicmatrix: {sortedStoicMatrix.shape}")
    interestedMatrix = buildInterestedMatrix(sortedStoicMatrix, projectOntoDimensions) # build G
    if verbose:
        logger.info("InterestedMatrix:")
        logger.info(interestedMatrix)
        print(f"InterestedMatrix: {interestedMatrix.shape}")
        
    eliminationMatrix = buildEliminationMatrix(sortedStoicMatrix, projectOntoDimensions) # build H
    #convertMatrixToHRepresentation(-eliminationMatrix.transpose(), os.path.join(outputDir, "tempRedundSort.ine"))
    #redund(numThreads,mplrsPath,os.path.join(outputDir, "tempRedundSort.ine"),os.path.join(outputDir, "outTempRedundSort.ine")) # remove redundancy
    #eliminationMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, "outTempRedundSort.ine"))
    eliminationMatrix = eliminationMatrix.transpose()
    eqFile = 'file.eq'
    iqFileIdentity = 'identity.iq'
    stage = 'projectionCone'
    outputFile = stage + '_out.zip' 
    if verbose:
        logger.info("Start Double Description:")
    # enumerate rays of W:
    #eliminationMatrixT = eliminationMatrix.transpose()
    eliminationMatrixT = eliminationMatrix
    print(f"elMatrixT: {eliminationMatrixT.shape}")
    identityEliminationMatrixT = np.identity(eliminationMatrix.shape[1])
    i = stoicMatrix.shape[1] - len(projectOntoDimensions)
    stepSize = min(5, i)
    print(f"i={i}, stepSize={stepSize}")
    boolStop = False
    matrix = np.empty((0, 30), float) 
    #matrix = np.empty((0, 409), float)
    # matrix = np.empty((0, 336), float)
    #projectedConeMatrix = interestedMatrix
    while i > 0:
        print(f"Iter: {i}")
        if (i - stepSize < stepSize):
            #convertEqualitiesAndInequalities2hRep(eliminationMatrixT[0:i,:], identityEliminationMatrixT, os.path.join(outputDir, "doubleDescription.ine"))
            print(f"elRows: 0:{i}")
            tempEliminationMatrix = eliminationMatrixT[0:i,:]
            boolStop = True
        else:
            print(f"elRows: {i-stepSize}:{i}")
            tempEliminationMatrix = eliminationMatrixT[i-stepSize:i,:]
            #convertEqualitiesAndInequalities2hRep(eliminationMatrixT[i-5:i,:], identityEliminationMatrixT, os.path.join(outputDir, "doubleDescription.ine"))
        #redund(numThreads,mplrsPath,os.path.join(outputDir, "doubleDescription.ine"),os.path.join(outputDir, "outDoubleDescription.ine")) # remove redundancy
        
        #rayMatrix = doubleDescriptionINE(os.path.join(outputDir, "outDoubleDescription.ine"), outputDir, outputFile, stage, numThreads=numThreads)
        rayMatrix = doubleDescription(identityEliminationMatrixT, iqFileIdentity, outputDir, outputFile, stage, tempEliminationMatrix, eqFile, numThreads=numThreads)
        logger.info(f"Raymatrix - Shape: {rayMatrix.shape}")
        print(f"Raymatrix - Shape: {rayMatrix.shape}")
        convertMatrixToVRepresentation(rayMatrix, os.path.join(outputDir, "rayMatrixTest.ine"))
        redund(numThreads,mplrsPath,os.path.join(outputDir, "rayMatrixTest.ine"),os.path.join(outputDir, f"{i}_redundRayMatrixTest.ine")) # remove redundancy
        printDatetime(f"Finished Raymatrix step {i-stepSize}-{i}:", datetime.now())
        matrix = np.vstack([matrix, getRaysFromVrepresentation(os.path.join(outputDir, f"{i}_redundRayMatrixTest.ine"))])
        #matrix = np.vstack([matrix, rayMatrix])
        #projectedConeMatrix = buildProjectedConeMatrix(rayMatrix, projectedConeMatrix)
        i = i - stepSize
        if boolStop:
            break
    printDatetime("Finished first double description step:", datetime.now())
    # rayMatrix = rayMatrix.T
    #logger.info(matrix)
    #print(matrix)
    convertMatrixToVRepresentation(matrix, os.path.join(outputDir, f"matrix.ine"))
    rayMatrix = matrix
    if verbose:
        logger.info("RayMatrix:")
        logger.info(rayMatrix)
        print(f"rayMatrix: {rayMatrix.shape}")
    ## rayMatrix[:,[1,3]] = rayMatrix[:,[3,1]]
    ## print("RayMatrix:\n", rayMatrix)

    # projectedConeMatrix = buildProjectedConeMatrix(rayMatrix.transpose(), interestedMatrix)
    projectedConeMatrix = -buildProjectedConeMatrix(rayMatrix, interestedMatrix) # build condition of projected cone
    #logger.info(f"projectedConeMatrix - Shape: {projectedConeMatrix.shape}")
    # projectedConeMatrix = removeRowsWithZeros(projectedConeMatrix)
    ineFile = os.path.join(outputDir, "projectedCone_H.ine")
    convertMatrixToHRepresentation(projectedConeMatrix, ineFile)
    printDatetime("Finished converting projected Cone Matrix to H-Rep:", datetime.now())

    #projectedConeMatrix[[2,3],:]= projectedConeMatrix[[3,2],:]
    #projectedConeMatrix = removeRedundancy(projectedConeMatrix)
    if not os.path.exists(os.path.join(outputDir, "temp")):
        os.mkdir(os.path.join(outputDir, "temp"))
    chunkSize = 10000
    #logger.info(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} with chunkSize {chunkSize}")
    #redundIterative(numThreads, mplrsPath, ineFile, os.path.join(outputDir, "tempOutprojectedCone_H.ine"), os.path.join(outputDir, "temp"), chunkSize=chunkSize, verbose=True)
    #logger.info(f"End: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    redund(numThreads,mplrsPath,ineFile,os.path.join(outputDir, "projectedConeMatrix_final_H.ine")) # remove redundancy
    printDatetime("Finished redund step:", datetime.now())

    projectedConeMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, "projectedConeMatrix_final_H.ine"))
    if verbose:
        logger.info("projectedConeMatrix:")
        logger.info(projectedConeMatrix)
        logger.info(f"Shape: {projectedConeMatrix.shape}")
        
    if stopAfterProjection:
        return None, None

    eqFileIdentity = 'identity.eq'
    iqFile = 'file.iq'
    stage = 'proCEMs'
    outputFile = stage + '_out.zip'
    if verbose:
        printDatetime("Start Double Description:", datetime.now())

    
    #proCEMs = doubleDescriptionINE(os.path.join(outputDir, "projectedConeMatrix_final_H.ine"), outputDir, outputFile, stage, numThreads=numThreads)
    #logger.info(f"proCEMs - Shape: {proCEMs.shape}")
    #print(f"proCEMs - Shape: {proCEMs.shape}")
    print(projectedConeMatrix)
    proCEMs = doubleDescription(projectedConeMatrix, iqFile, outputDir, outputFile, stage, numThreads=numThreads)
    logger.info(f"proCEMs - Shape: {proCEMs.shape}")
    print(f"proCEMs - Shape: {proCEMs.shape}")
    printDatetime("Finished second double description step:", datetime.now())
    # proCEMs= proCEMs.T
    if verbose:
        print("\n#########\nResulting proCEMs:\n", proCEMs)
    efps = calculateEFPsFromProCEMs(proCEMs)
    printDatetime("Finished calculating efps:", datetime.now())
    if verbose:
        print("EFPs:", efps)
        print("EFPs-length:", len(efps))
    return proCEMs, efps

def runMarashiWithPolco(stoicMatrix, projectOntoDimensions, outputDir, numThreads=8, stopAfterProjection=True, verbose = True):
    printDatetime("Start: ", datetime.now())
    sortedStoicMatrix = stoicMatrix#sortStoicMatrix(stoicMatrix, projectOntoDimensions)
    #convertEqualities2hRep(sortedStoicMatrix, os.path.join(outputDir, "tempRedundSort.ine"))
    #redund(numThreads,mplrsPath,os.path.join(outputDir, "tempRedundSort.ine"),os.path.join(outputDir, "outTempRedundSort.ine")) # remove redundancy
    #sortedStoicMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, "outTempRedundSort.ine"))
    convertMatrixToHRepresentation(sortedStoicMatrix, os.path.join(outputDir, "stoicMatrix.ine"))
    eliminationMatrix = buildEliminationMatrix(sortedStoicMatrix, projectOntoDimensions) # build H
    convertMatrixToHRepresentation(-eliminationMatrix.transpose(), os.path.join(outputDir, "tempRedundSort.ine"))
    redund(numThreads,mplrsPath,os.path.join(outputDir, "tempRedundSort.ine"),os.path.join(outputDir, "outTempRedundSort.ine")) # remove redundancy
    eliminationMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, "outTempRedundSort.ine"))
    
    eqFile = 'file.eq'
    iqFileIdentity = 'identity.iq'
    stage = 'projectionCone'
    outputFile = stage + '_out.zip'
    if verbose:
        logger.info("Start Double Description:")
    # enumerate rays of W:
    rayMatrix = doubleDescription(-np.vstack([eliminationMatrix, np.identity(eliminationMatrix.shape[1])]), iqFileIdentity, outputDir, outputFile, stage, numThreads=numThreads)
    #rayMatrix = doubleDescription(np.identity(eliminationMatrix.shape[1]), iqFileIdentity, outputDir, outputFile, stage, eliminationMatrix, eqFile, numThreads=numThreads)
    #convertEqualitiesAndInequalities2hRep(eliminationMatrix, np.identity(eliminationMatrix.shape[1]), os.path.join(outputDir, "doubleDescription.ine"))
    #redund(numThreads,mplrsPath,os.path.join(outputDir, "doubleDescription.ine"),os.path.join(outputDir, "outDoubleDescription.ine")) # remove redundancy
    
    #rayMatrix = doubleDescriptionINE(os.path.join(outputDir, "outDoubleDescription.ine"), outputDir, outputFile, stage, numThreads=numThreads)
    logger.info(f"Raymatrix - Shape: {rayMatrix.shape}")
    printDatetime("Finished first double description step:", datetime.now())
    # rayMatrix = rayMatrix.T
    if verbose:
        logger.info("RayMatrix:")
        logger.info(rayMatrix)
    # rayMatrix[:,[1,3]] = rayMatrix[:,[3,1]]
    # print("RayMatrix:\n", rayMatrix)
    interestedMatrix = buildInterestedMatrix(sortedStoicMatrix, projectOntoDimensions) # build G
    if verbose:
        logger.info("InterestedMatrix:")
        logger.info(interestedMatrix)
    # projectedConeMatrix = buildProjectedConeMatrix(rayMatrix.transpose(), interestedMatrix)
    projectedConeMatrix = -buildProjectedConeMatrix(rayMatrix, interestedMatrix) # build condition of projected cone
    logger.info(f"projectedConeMatrix - Shape: {projectedConeMatrix.shape}")
    # projectedConeMatrix = removeRowsWithZeros(projectedConeMatrix)
    ineFile = os.path.join(outputDir, "projectedCone_H.ine")
    convertMatrixToHRepresentation(projectedConeMatrix, ineFile)
    printDatetime("Finished converting projected Cone Matrix to H-Rep:", datetime.now())

    #projectedConeMatrix[[2,3],:]= projectedConeMatrix[[3,2],:]
    #projectedConeMatrix = removeRedundancy(projectedConeMatrix)
    if not os.path.exists(os.path.join(outputDir, "temp")):
        os.mkdir(os.path.join(outputDir, "temp"))
    chunkSize = 10000
    logger.info(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} with chunkSize {chunkSize}")
    redundIterative(numThreads, mplrsPath, ineFile, os.path.join(outputDir, "tempOutprojectedCone_H.ine"), os.path.join(outputDir, "temp"), chunkSize=chunkSize, verbose=True)
    logger.info(f"End: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    redund(numThreads,mplrsPath,os.path.join(outputDir, "tempOutprojectedCone_H.ine"),os.path.join(outputDir, "projectedConeMatrix_final_H.ine")) # remove redundancy
    printDatetime("Finished redund step:", datetime.now())

    projectedConeMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, "projectedConeMatrix_final_H.ine"))
    if verbose:
        logger.info("projectedConeMatrix:")
        #logger.info(projectedConeMatrix)
        logger.info(f"Shape: {projectedConeMatrix.shape}")
        print("Shape: ", projectedConeMatrix.shape)
    if stopAfterProjection:
        return None, None

    eqFileIdentity = 'identity.eq'
    iqFile = 'file.iq'
    stage = 'proCEMs'
    outputFile = stage + '_out.zip'
    if verbose:
        printDatetime("Start Double Description:", datetime.now())
    proCEMs = doubleDescription(projectedConeMatrix, iqFile, outputDir, outputFile, stage, numThreads=numThreads)
    printDatetime("Finished second double description step:", datetime.now())
    if verbose:
        # print("\n#########\nBefore redund - proCEMs:\n", proCEMs)
        print("proCEMs:\n", proCEMs)
        print("projected Matrix-Shape: ", proCEMs.shape)
    # convertMatrixToVRepresentation(proCEMs, os.path.join(outputDir, "redund_proCEMs_V.ine"))
    # redund(numThreads,mplrsPath,os.path.join(outputDir, "redund_proCEMs_V.ine"),os.path.join(outputDir, "outRedundProCEMs_V.ine"))
    # printDatetime("Finished redund step for proCEMs:", datetime.now())
    # proCEMs = getRaysFromVrepresentation(os.path.join(outputDir, "outRedundProCEMs_V.ine"))
    ## proCEMs= proCEMs.T
    # if verbose:
        # print("\n#########\nAfter redund - proCEMs:\n", proCEMs)
        # print("Shape: ", proCEMs.shape)
    efps = calculateEFPsFromProCEMs(proCEMs)
    printDatetime("Finished calculating efps:", datetime.now())
    if verbose:
        print("EFPs:", efps)
        print("EFPs-length:", len(efps))
    return proCEMs, efps

def runMarashiWithMPLRS(stoicMatrix, projectOntoDimensions, outputDir, mplrsPath, numThreads = 8, stopAfterProjection=True, verbose = True):
    printDatetime("Start:", datetime.now())
    sortedStoicMatrix = sortStoicMatrix(stoicMatrix, projectOntoDimensions)
    eliminationMatrix = buildEliminationMatrix(sortedStoicMatrix, projectOntoDimensions) # build H
    
    hEliminationInePath = os.path.join(outputDir, "elimination_H.ine")
    hEliminationRedundInePath = os.path.join(outputDir, "elimination_redund_H.ine")
    vEliminationInePath = os.path.join(outputDir, "elimination_V.ine")
    convertEqualitiesAndInequalities2hRep(eliminationMatrix.transpose(), np.identity(eliminationMatrix.shape[0]), hEliminationInePath) # build W as ine file
    printDatetime("Finished creating eliminiation_H.ine:", datetime.now())
    redund(numThreads,mplrsPath,hEliminationInePath,hEliminationRedundInePath) # remove redundancy
    mplrs_conversion(numThreads, mplrsPath, hEliminationRedundInePath, vEliminationInePath) # convert to v-representation
    printDatetime("Finished mplrs conversion:", datetime.now())
    rayMatrix = getRaysFromVrepresentation(vEliminationInePath)
    logger.info(f"Raymatrix - Shape: {rayMatrix.shape}")
    if verbose:
        logger.info("RayMatrix:")
        logger.info(rayMatrix)

    interestedMatrix = buildInterestedMatrix(sortedStoicMatrix, projectOntoDimensions) # build G
    if verbose:
        logger.info("InterestedMatrix:")
        logger.info(interestedMatrix)
    projectedConeMatrix = -buildProjectedConeMatrix(rayMatrix, interestedMatrix)
    logger.info(f"projectedConeMatrix - Shape: {projectedConeMatrix.shape}")

    hProjectedConeInePath = os.path.join(outputDir, "projectedCone_H.ine")
    convertMatrixToHRepresentation(projectedConeMatrix, hProjectedConeInePath)
    printDatetime("Finished converting projected Cone Matrix to H-Rep:", datetime.now())

    hRedundProjectedConeInePath = os.path.join(outputDir, "projectedConeMatrix_H_redund.ine")
    if not os.path.exists(os.path.join(outputDir,"temp")):
        os.mkdir(os.path.join(outputDir,"temp"))
    chunkSize = 10000
    logger.info(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} with chunkSize {chunkSize}")
    redundIterative(numThreads, mplrsPath, hProjectedConeInePath, os.path.join(outputDir, "tempOutprojectedCone_H.ine"), os.path.join(outputDir,"temp"), chunkSize=chunkSize, verbose=True)
    logger.info(f"End: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    redund(numThreads,mplrsPath,os.path.join(outputDir, "tempOutprojectedCone_H.ine"),hRedundProjectedConeInePath) # remove redundancy

    #redund(numThreads, mplrsPath, hProjectedConeInePath, hRedundProjectedConeInePath)
    printDatetime("Finished redund step:", datetime.now())

    if verbose:
        logger.info("projectedConeMatrix:")
        logger.info(projectedConeMatrix)
        logger.info(f"Shape: {projectedConeMatrix.shape}")
    if stopAfterProjection:
        return None, None

    vProjectedConeInePath = os.path.join(outputDir, "projectedConeMatrix_V.ine")
    convertHtoVrepresentation(hRedundProjectedConeInePath, vProjectedConeInePath)
    printDatetime("Finished H-Rep to V-Rep conversion:", datetime.now())
    proCEMs = getRaysFromVrepresentation(vProjectedConeInePath)
    if verbose:
        print("\n#########\nResulting proCEMs:\n", proCEMs)
    efps = calculateEFPsFromProCEMs(proCEMs)
    printDatetime("Finished calculating efps:", datetime.now())
    if verbose:
        print("EFPs:", efps)
    return proCEMs, efps

def runFELWithPolco(stoicMatrix, projectOntoDimensions, outputDir, mplrsPath, mplrsProcesses=4, verbose=True):
    printDatetime("Start:", datetime.now())
    sortedStoicMatrix = sortStoicMatrix(stoicMatrix, projectOntoDimensions)
    hBeforeEliminationInePath = outputDir + "/beforeElimination_H.ine"
    convertEqualitiesAndInequalities2hRep(sortedStoicMatrix, np.identity(stoicMatrix.shape[1]), hBeforeEliminationInePath)
    printDatetime("Finished creating beforeElimination_H.ine:", datetime.now())
    eliminateDimensions = []
    for i in range(len(projectOntoDimensions), stoicMatrix.shape[1]):
        eliminateDimensions.append(i+1)
    hProjectedInePath = outputDir + "/projected_H.ine"
    eliminate(hBeforeEliminationInePath, eliminateDimensions, hProjectedInePath, mplrsPath, mplrsProcesses, verbose)
    printDatetime("Finished FEL elimination:", datetime.now())

    projectedMatrix = getMatrixFromHrepresentation(hProjectedInePath)
    if verbose:
        print("Projected Matrix:\n", projectedMatrix)

    iqFile = 'file.iq'
    stage = 'Rays'
    outputFile = stage + '_out.zip'
    if verbose:
        print("\nStart Double Description:")
    proCEMs = doubleDescription(projectedMatrix, iqFile, outputDir, outputFile, stage)
    printDatetime("Finished double description step:", datetime.now())
    if verbose:
        print("\n#########\nResulting Rays:\n", proCEMs)
    efps = calculateEFPsFromProCEMs(proCEMs)
    printDatetime("Finished calculating efps:", datetime.now())
    if verbose:
        print("EFPs:", efps)
    
    return proCEMs, efps