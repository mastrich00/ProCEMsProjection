import numpy as np
import subprocess
import time
from datetime import datetime
from projection.helpers import printDatetime
from projection.marashiProjectionHelpers import sortStoicMatrix, buildInterestedMatrix, buildEliminationMatrix, buildProjectedConeMatrix, \
    calculateEFPsFromProCEMs, logTimeToCSV
from projection.enumeration import doubleDescription, doubleDescriptionINE, convertMatrixToHRepresentation, convertMatrixToVRepresentation, getMatrixFromHrepresentation,convertEqualitiesAndInequalities2hRep, \
    getRaysFromVrepresentation, convertHtoVrepresentation, mplrs_conversion, convertEqualities2hRep
from projection.redund import redund, redundIterative
from projection.logging_config import logger
import os
from fractions import Fraction
import os
import time
import shutil

#import mpmath
#mpmath.mp.dps = 100 

#def has_only_2_and_5_factors(n):
#    """Check if the denominator has only 2 and/or 5 as prime factors."""
#    while n % 2 == 0:
#        n //= 2
#    while n % 5 == 0:
#        n //= 5
#    return n == 1  # If 1 remains, it had only 2 and 5 as factors

#def fraction_to_safe_float(frac):
#    """Convert fraction to float only if it has no repeating decimals."""
#    if has_only_2_and_5_factors(frac.denominator):
#        #print(frac)
#        return mpmath.mpf(frac.numerator) / mpmath.mpf(frac.denominator)
#        # return float(frac)  # Convert to double precision
#    return frac  # Keep as fraction

#def float_to_fraction(value):
#    """Convert float back to fraction if it was originally converted."""
#    if isinstance(value, float):
#        return Fraction(value).limit_denominator()  # Convert back to fraction
#    return value  # Keep original fractions unchanged

def runMarashiWithMPLRSSubsets(stoicMatrix, projectOntoDimensions, outputDir, numThreads=8, stopAfterProjection=True, verbose = True, iteration=0, originalProjectionReactions=[], logTimesFile="times.csv", reversibleList=[], MPLRS_PATH="mplrs", boolCalcEFPs=False):
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
        redund(numThreads,MPLRS_PATH, os.path.join(outputDir, f"{iteration}_stoicMatrix.ine"),os.path.join(outputDir, f"{iteration}_stoicMatrix_redund.ine"))
        sortedStoicMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, f"{iteration}_stoicMatrix_redund.ine"))
        
    interestedMatrix = sortedStoicMatrix[:,:len(projectOntoDimensions)]
    eliminationMatrix = sortedStoicMatrix[:,len(projectOntoDimensions):]
    eliminationMatrix = eliminationMatrix.transpose()
    print(f"shape elimination: {eliminationMatrix.shape}")
    print(f"shape interest: {interestedMatrix.shape}")
    time_setup_stoic_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Setup Elimination-Matrix", time_setup_stoic_end - time_setup_stoic_start)
    
    time_projectioncone_redund_start = time.time()
    hEliminationInePath = os.path.join(outputDir, "elimination_H.ine")
    hEliminationRedundInePath = os.path.join(outputDir, "elimination_redund_H.ine")
    vEliminationInePath = os.path.join(outputDir, "elimination_V.ine")
    convertEqualitiesAndInequalities2hRep(eliminationMatrix, np.identity(eliminationMatrix.shape[1]), hEliminationInePath) # build W as ine file
    printDatetime("Finished creating eliminiation_H.ine:", datetime.now())
    redund(numThreads,MPLRS_PATH,hEliminationInePath,hEliminationRedundInePath) # remove redundancy
    time_projectioncone_redund_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Redund Projection-Cone", time_projectioncone_redund_end - time_projectioncone_redund_start)

    time_projectioncone_conversion_start = time.time()
    mplrs_conversion(numThreads, MPLRS_PATH, hEliminationRedundInePath, vEliminationInePath) # convert to v-representation
    time_projectioncone_conversion_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Enumeration Projection-Cone", time_projectioncone_conversion_end - time_projectioncone_conversion_start)

    vEliminationBeforeRedundInePath = os.path.join(outputDir, "elimination_before_redund_V.ine")
    vEliminationRedundInePath = os.path.join(outputDir, "elimination_redund_V.ine")
    eliminationRayMatrix = getRaysFromVrepresentation(vEliminationInePath)
    convertMatrixToVRepresentation(eliminationRayMatrix, vEliminationBeforeRedundInePath)
    redund(numThreads,MPLRS_PATH,vEliminationBeforeRedundInePath,vEliminationRedundInePath)

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
    convertMatrixToHRepresentation(projectedConeMatrix, hProjectedConeInePath)
    printDatetime("Finished converting projected Cone Matrix to H-Rep:", datetime.now())

#    if projectedConeMatrix.shape[0] > 300000: # protect server from crashing
#        print(f"Error: shape={projectedConeMatrix.shape}")
#        exit(0)

    hRedundProjectedConeInePath = os.path.join(outputDir, "projectedConeMatrix_H_redund.ine")
    if not os.path.exists(os.path.join(outputDir,"temp")):
        os.mkdir(os.path.join(outputDir,"temp"))

    redund(numThreads, MPLRS_PATH, hProjectedConeInePath, hRedundProjectedConeInePath)
    printDatetime("Finished redund step:", datetime.now())
    matrix = getMatrixFromHrepresentation(hRedundProjectedConeInePath)
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
    mplrs_conversion(numThreads,MPLRS_PATH, hRedundProjectedConeInePath, vProjectedConeInePath)
    printDatetime("Finished H-Rep to V-Rep conversion:", datetime.now())
    proCEMs = getRaysFromVrepresentation(vProjectedConeInePath)
    time_enumerate_projectedcone_end = time.time()
    logTimeToCSV(logTimesFile, f"Enum. ProCEMs", "Enumerating ProCEMs", time_enumerate_projectedcone_end - time_enumerate_projectedcone_start)
    if verbose:
        print("\n#########\nBefore Redund - proCEMs:\n", proCEMs)
        print("\n#########\nResulting proCEMs-Shape:\n", proCEMs.shape)

    time_redund_procems_start = time.time()
    vProCEMsPath = os.path.join(outputDir,"proCEMs_V.ine")
    redundvProCEMsPath = os.path.join(outputDir,"redund_proCEMs_V.ine")
    convertMatrixToVRepresentation(proCEMs, vProCEMsPath)
    #if proCEMs.shape[0] > 300000:
    #    print(f"Error: shape={projectedConeMatrix.shape}")
    #    exit(0)
    redund(numThreads,MPLRS_PATH,vProCEMsPath,redundvProCEMsPath)
    printDatetime("Finished redund step for proCEMs:", datetime.now())
    proCEMs = getRaysFromVrepresentation(redundvProCEMsPath)
    time_redund_procems_end = time.time()
    logTimeToCSV(logTimesFile, f"Enum. ProCEMs", "Redund proCEMs", time_redund_procems_end - time_redund_procems_start)
    if verbose:
        print("\n#########\nAfter Redund - proCEMs:\n", proCEMs)
        print("\n#########\nResulting proCEMs-Shape:\n", proCEMs.shape)

    efps = []
    if boolCalcEFPs:
        time_calc_efps_start = time.time()
        efps = calculateEFPsFromProCEMs(proCEMs)
        time_calc_efps_end = time.time()
        logTimeToCSV(logTimesFile, f"Iter {iteration}", "Calc EFPs", time_calc_efps_end - time_calc_efps_start)
        printDatetime("Finished calculating efps:", datetime.now())
        if verbose:
            print("Number of EFPs:", len(efps))

    parent_folder = os.path.dirname(os.path.abspath(outputDir))
    shutil.copy(redundvProCEMsPath, os.path.join(parent_folder, "proCEMs.ine"))
    
    return proCEMs, efps

def runMarashiWithPolcoSubsets(stoicMatrix, projectOntoDimensions, outputDir, numThreads=8, stopAfterProjection=True, verbose = True, iteration=0, originalProjectionReactions=[], logTimesFile="times.csv", reversibleList=[], MPLRS_PATH="mplrs", POLCO_PATH="polco.jar", safetyThresholdRAM=0, boolCalcEFPs=False):
    time_setup_stoic_start = time.time()
    printDatetime("Start: ", datetime.now())
    sortedStoicMatrix = stoicMatrix
    
    convertMatrixToHRepresentation(sortedStoicMatrix, os.path.join(outputDir, "stoicMatrix.ine"))
    if iteration == 0:
        p = len(originalProjectionReactions)
        q = sortedStoicMatrix.shape[1] - p
        convertEqualitiesAndInequalities2hRep(sortedStoicMatrix, -np.identity(p+q)[reversibleList],
            os.path.join(outputDir, f"{iteration}_stoicMatrix.ine")
        )
        redund(numThreads,MPLRS_PATH, os.path.join(outputDir, f"{iteration}_stoicMatrix.ine"),os.path.join(outputDir, f"{iteration}_stoicMatrix_redund.ine"))
        sortedStoicMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, f"{iteration}_stoicMatrix_redund.ine"))   
    else:
        interestedMatrix = sortedStoicMatrix[:, :len(projectOntoDimensions)]
        eliminationMatrix = sortedStoicMatrix[:, len(projectOntoDimensions):]
    interestedMatrix = sortedStoicMatrix[:,:len(projectOntoDimensions)]
    eliminationMatrix = sortedStoicMatrix[:,len(projectOntoDimensions):]
    print(f"shape el: {eliminationMatrix.shape}")
    print(f"shape int: {interestedMatrix.shape}")
    #convertMatrixToHRepresentation(-eliminationMatrix.transpose(), os.path.join(outputDir, "tempRedundSort.ine"))
    #redund(numThreads,MPLRS_PATH,os.path.join(outputDir, "tempRedundSort.ine"),os.path.join(outputDir, "outTempRedundSort.ine")) # remove redundancy
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
    convertEqualitiesAndInequalities2hRep(eliminationMatrix, np.identity(eliminationMatrix.shape[1]), os.path.join(outputDir, "doubleDescription.ine"))
    #redund(numThreads,MPLRS_PATH,os.path.join(outputDir, "doubleDescription.ine"),os.path.join(outputDir, "outDoubleDescription.ine")) # remove redundancy
    rayMatrix = doubleDescriptionINE(os.path.join(outputDir, "doubleDescription.ine"), outputDir, outputFile, stage, numThreads=numThreads, POLCO_PATH=POLCO_PATH)
    time_projectioncone_conversion_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Enumeration Projection-Cone", time_projectioncone_conversion_end - time_projectioncone_conversion_start)
    
    # print(f"Before redund - elimination rays: {rayMatrix.shape}")
    # time_rays_redund_start = time.time()
    # elimination_rays_V_path = os.path.join(outputDir, "elimination_rays_V.ine")
    # elimination_rays_V_Redund_path = os.path.join(outputDir, "elimination_rays_V_redund.ine")
    # convertMatrixToVRepresentation(rayMatrix, elimination_rays_V_path)
    # redund(numThreads,MPLRS_PATH,elimination_rays_V_path,elimination_rays_V_Redund_path)
    # rayMatrix = getRaysFromVrepresentation(elimination_rays_V_Redund_path)
    # print(f"After redund - elimination rays: {rayMatrix.shape}")
    logger.info(f"Raymatrix - Shape: {rayMatrix.shape}")
    printDatetime("Finished first double description step:", datetime.now())
    if verbose:
        logger.info("RayMatrix:")
        logger.info(rayMatrix)
    # time_rays_redund_end = time.time()
    # logTimeToCSV(logTimesFile, f"Iter {iteration}", "Redund Rays", time_rays_redund_end - time_rays_redund_start)
    
    time_build_projectedcone_start = time.time()
    printDatetime("Calculating projected cone matrix:", datetime.now())
    projectedConeMatrix = buildProjectedConeMatrix(rayMatrix, interestedMatrix) # build condition of projected cone
    printDatetime("Finished calculating projected cone matrix:", datetime.now())
    time_build_projectedcone_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Build projected Cone", time_build_projectedcone_end - time_build_projectedcone_start)
    logger.info(f"projectedConeMatrix - Shape: {projectedConeMatrix.shape}")

    time_redund_projectedcone_start = time.time()
    ineFile = os.path.join(outputDir, "projectedCone_H.ine")
    convertMatrixToHRepresentation(projectedConeMatrix, ineFile)
    printDatetime("Finished converting projected Cone Matrix to H-Rep:", datetime.now())
    if safetyThresholdRAM > 0 and projectedConeMatrix.shape[0] > safetyThresholdRAM:
        print(f"Error: too many rays (RAM Safety Threshold). Shape={projectedConeMatrix.shape}")
        exit(0)
    redund(numThreads,MPLRS_PATH,ineFile,os.path.join(outputDir, f"{iteration}_projectedConeMatrix_final_H.ine")) # remove redundancy
    printDatetime("Finished redund step:", datetime.now())
    projectedConeMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, f"{iteration}_projectedConeMatrix_final_H.ine"))
    time_redund_projectedcone_end = time.time()
    logTimeToCSV(logTimesFile, f"Iter {iteration}", "Redund projected Cone", time_redund_projectedcone_end - time_redund_projectedcone_start)
            
    if stopAfterProjection:
        return projectedConeMatrix, None
        
    time_enumerate_projectedcone_start = time.time()
    eqFileIdentity = 'identity.eq'
    iqFile = 'file.iq'
    stage = 'proCEMs'
    outputFile = stage + '_out.zip'
    if verbose:
        printDatetime("Start Double Description:", datetime.now())
    convertMatrixToHRepresentation(projectedConeMatrix, os.path.join(outputDir, "doubleDescriptionProjectedCone.ine"))
    proCEMs = doubleDescriptionINE(os.path.join(outputDir, "doubleDescriptionProjectedCone.ine"), outputDir, outputFile, stage, numThreads=numThreads,POLCO_PATH=POLCO_PATH)
    time_enumerate_projectedcone_end= time.time()
    logTimeToCSV(logTimesFile, f"Enum. ProCEMs", "Enumerating ProCEMs", time_enumerate_projectedcone_end - time_enumerate_projectedcone_start)
    printDatetime("Finished second double description step:", datetime.now())

    time_redund_procems_start = time.time()
    convertMatrixToVRepresentation(proCEMs, os.path.join(outputDir, "redund_proCEMs_V.ine"))
    if safetyThresholdRAM > 0 and proCEMs.shape[0] > safetyThresholdRAM:
        print(f"Error: shape={proCEMs.shape}")
        exit(0)

    redundvProCEMsPath = os.path.join(outputDir, "outRedundProCEMs_V.ine")
    redund(numThreads,MPLRS_PATH, os.path.join(outputDir, "redund_proCEMs_V.ine"), redundvProCEMsPath)
    printDatetime("Finished redund step for proCEMs:", datetime.now())
    proCEMs = getRaysFromVrepresentation(redundvProCEMsPath)
    time_redund_procems_end = time.time()
    logTimeToCSV(logTimesFile, f"Enum. ProCEMs", "Redund proCEMs", time_redund_procems_end - time_redund_procems_start)


    efps = []
    if boolCalcEFPs:
        time_calc_efps_start = time.time()
        efps = calculateEFPsFromProCEMs(proCEMs)
        time_calc_efps_end = time.time()
        logTimeToCSV(logTimesFile, f"Iter {iteration}", "Calc EFPs", time_calc_efps_end - time_calc_efps_start)
        printDatetime("Finished calculating efps:", datetime.now())
        if verbose:
            print("Number of EFPs:", len(efps))

    parent_folder = os.path.dirname(os.path.abspath(outputDir))
    shutil.copy(redundvProCEMsPath, os.path.join(parent_folder, "proCEMs.ine"))
                                         
    return proCEMs, efps

def eliminate(ineFilePath: str, dimensionsToEliminate: list, outFilePath: str, MPLRS_PATH: str, numProc= 4, verbose = True, logEliminations = False):
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
            cmd = ["mpirun", "-np", str(numProc), MPLRS_PATH, "-fel", tmpFilePath+str(i), outFilePath+str(i)] # TODO: check if -fel necessary
        else:
            cmd = ["mpirun", "-np", str(numProc), MPLRS_PATH, "-fel", tmpFilePath, outFilePath]
        subprocess.run(cmd)

        if logEliminations:
            subprocess.run(["cp", outFilePath+str(i), tmpFilePath+str(i+1)])
        else:
            subprocess.run(["cp", outFilePath, tmpFilePath])

    if verbose:
        logger.info("Running elimination step " + str(len(dimensionsToEliminate)) + "/" + str(len(dimensionsToEliminate)))
    if logEliminations:
        cmd = ["mpirun", "-np", str(numProc), MPLRS_PATH, "-fel", tmpFilePath+str(len(dimensionsToEliminate)), outFilePath+str(len(dimensionsToEliminate))]
    else:
        cmd = ["mpirun", "-np", str(numProc), MPLRS_PATH, "-fel", tmpFilePath, outFilePath]
    subprocess.run(cmd)
    
    subprocess.run(["cp", outFilePath+str(len(dimensionsToEliminate)), outFilePath])


def runFELWithPolco(stoicMatrix, projectOntoDimensions, outputDir, MPLRS_PATH, mplrsProcesses=4, verbose=True):
    printDatetime("Start:", datetime.now())
    sortedStoicMatrix = sortStoicMatrix(stoicMatrix, projectOntoDimensions)
    hBeforeEliminationInePath = outputDir + "/beforeElimination_H.ine"
    convertEqualitiesAndInequalities2hRep(sortedStoicMatrix, np.identity(stoicMatrix.shape[1]), hBeforeEliminationInePath)
    printDatetime("Finished creating beforeElimination_H.ine:", datetime.now())
    eliminateDimensions = []
    for i in range(len(projectOntoDimensions), stoicMatrix.shape[1]):
        eliminateDimensions.append(i+1)
    hProjectedInePath = outputDir + "/projected_H.ine"
    eliminate(hBeforeEliminationInePath, eliminateDimensions, hProjectedInePath, MPLRS_PATH, mplrsProcesses, verbose)
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