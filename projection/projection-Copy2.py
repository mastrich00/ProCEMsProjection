import numpy as np
import subprocess
from datetime import datetime
from projection.helpers import printDatetime
from projection.marashiProjectionHelpers import sortStoicMatrix, buildInterestedMatrix, buildEliminationMatrix, buildProjectedConeMatrix, \
    calculateEFPsFromProCEMs
from projection.enumeration import doubleDescription, doubleDescriptionINE, convertMatrixToHRepresentation, convertMatrixToVRepresentation, getMatrixFromHrepresentation,convertEqualitiesAndInequalities2hRep, \
    getRaysFromVrepresentation, convertHtoVrepresentation, mplrs_conversion, convertEqualities2hRep
from projection.redund import redund, redundIterative
from projection.logging_config import logger
import os
if not os.getenv("POLCO_PATH"):
    load_dotenv("env") # load env

POLCO_PATH = os.getenv("POLCO_PATH")
mplrsPath = MPLRS_PATH = os.getenv("MPLRS_PATH")
mplrsPath = MPLRS_PATH = "mplrsV7_2"
import os

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

def runMarashiWithPolcoIterative(stoicMatrix, projectOntoDimensions, outputDir, numThreads=8, stopAfterProjection=True, verbose = True):
    print(f"Number projection dims: {len(projectOntoDimensions)}")
    print(f"shape stoic matrix: {stoicMatrix.shape}")
    printDatetime("Start: ", datetime.now())
    sortedStoicMatrix = sortStoicMatrix(stoicMatrix, projectOntoDimensions)
    print(f"stoicmatrix: {sortedStoicMatrix.shape}")
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
    stepSize = min(4, i)
    print(f"i={i}, stepSize={stepSize}")
    boolStop = False
    matrix = np.empty((0, 30), float) 
    #matrix = np.empty((0, 409), float)
    # matrix = np.empty((0, 336), float)
    while i > 0:
        print(f"Iter: {i}")
        if (i - stepSize <= stepSize):
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
        #convertMatrixToHRepresentation(rayMatrix[:,:], os.path.join(outputDir, "rayMatrixTest.ine"))
        #redund(numThreads,mplrsPath,os.path.join(outputDir, "rayMatrixTest.ine"),os.path.join(outputDir, f"{i}_redundRayMatrixTest.ine")) # remove redundancy
        printDatetime(f"Finished Raymatrix step {i-stepSize}-{i}:", datetime.now())
        #matrix = np.vstack([matrix, getMatrixFromHrepresentation(os.path.join(outputDir, f"{i}_redundRayMatrixTest.ine"))])
        matrix = np.vstack([matrix, rayMatrix])
        i = i - stepSize
        if boolStop:
            break
    printDatetime("Finished first double description step:", datetime.now())
    # rayMatrix = rayMatrix.T
    logger.info(matrix)
    print(matrix)
    convertMatrixToHRepresentation(matrix, os.path.join(outputDir, f"{i}_redundRayMatrixTest.ine"))
    rayMatrix = matrix
    if verbose:
        logger.info("RayMatrix:")
        logger.info(rayMatrix)
        print(f"rayMatrix: {rayMatrix.shape}")
    # rayMatrix[:,[1,3]] = rayMatrix[:,[3,1]]
    # print("RayMatrix:\n", rayMatrix)
    interestedMatrix = buildInterestedMatrix(sortedStoicMatrix, projectOntoDimensions) # build G
    if verbose:
        logger.info("InterestedMatrix:")
        logger.info(interestedMatrix)
        print(f"InterestedMatrix: {interestedMatrix.shape}")
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
    sortedStoicMatrix = sortStoicMatrix(stoicMatrix, projectOntoDimensions)
    
    convertEqualities2hRep(sortedStoicMatrix, os.path.join(outputDir, "tempRedundSort.ine"))
    redund(numThreads,mplrsPath,os.path.join(outputDir, "tempRedundSort.ine"),os.path.join(outputDir, "outTempRedundSort.ine")) # remove redundancy
    sortedStoicMatrix = getMatrixFromHrepresentation(os.path.join(outputDir, "outTempRedundSort.ine"))

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
    # rayMatrix = doubleDescription(np.identity(eliminationMatrix.shape[1]), iqFileIdentity, outputDir, outputFile, stage, eliminationMatrix, eqFile, numThreads=numThreads)
    convertEqualitiesAndInequalities2hRep(eliminationMatrix, np.identity(eliminationMatrix.shape[1]), os.path.join(outputDir, "doubleDescription.ine"))
    redund(numThreads,mplrsPath,os.path.join(outputDir, "doubleDescription.ine"),os.path.join(outputDir, "outDoubleDescription.ine")) # remove redundancy
    
    rayMatrix = doubleDescriptionINE(os.path.join(outputDir, "outDoubleDescription.ine"), outputDir, outputFile, stage, numThreads=numThreads)
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
    proCEMs = doubleDescription(projectedConeMatrix, iqFile, outputDir, outputFile, stage, numThreads=numThreads)
    printDatetime("Finished second double description step:", datetime.now())
    # if verbose:
    #     print("\n#########\nBefore redund - proCEMs:\n", proCEMs)
    #     print("Shape: ", proCEMs.shape)
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
