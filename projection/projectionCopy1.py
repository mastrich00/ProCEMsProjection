import numpy as np
import subprocess
from datetime import datetime
from projection.helpers import printDatetime
from projection.marashiProjectionHelpers import sortStoicMatrix, buildInterestedMatrix, buildEliminationMatrix, buildProjectedConeMatrix, \
    calculateEFPsFromProCEMs
from projection.enumeration import doubleDescription, convertMatrixToHRepresentation, getMatrixFromHrepresentation,convertEqualitiesAndInequalities2hRep, \
    getRaysFromVrepresentation, convertHtoVrepresentation, mplrs_conversion
from projection.redund import redund
polcoPath = "polco.jar"
mplrsPath = "mplrs"
import os

def eliminate(ineFilePath: str, dimensionsToEliminate: list, outFilePath: str, mplrsPath: str, numProc= 4, verbose = True, logEliminations = False):
    """
    Eliminates given dimensions in .ine file. Writes output to 'outFilePath'
    Dimension index start at 1, not zero.
    """
    if len(dimensionsToEliminate) == 0:
        if verbose:
            print("No dimensions to eliminate!")
        return 0
    
    file = open(ineFilePath, "a")
    cmd = "eliminate " + str(len(dimensionsToEliminate))
    for index in dimensionsToEliminate:
        cmd = cmd + " " + str(index)
    if verbose:
        print("Command: " + cmd)
    file.write(cmd)
    if verbose:
        file.write("\nverbose")
    file.close()

    tmpFilePath = ineFilePath + "_tmp.ine"
    if logEliminations:
        subprocess.run(["cp",ineFilePath, tmpFilePath + "1"])
    else:
        subprocess.run(["cp",ineFilePath, tmpFilePath])

    for i in range(1, len(dimensionsToEliminate)):
        if verbose:
            print("Running elimination step " + str(i) + "/" + str(len(dimensionsToEliminate)))
        if logEliminations:
            cmd = ["mpirun", "-np", str(numProc), mplrsPath,tmpFilePath+str(i), outFilePath+str(i)]
        else:
            cmd = ["mpirun", "-np", str(numProc), mplrsPath, tmpFilePath, outFilePath]
        subprocess.run(cmd)

        if logEliminations:
            subprocess.run(["cp", outFilePath+str(i), tmpFilePath+str(i+1)])
        else:
            subprocess.run(["cp", outFilePath, tmpFilePath])

    if verbose:
        print("Running elimination step " + str(len(dimensionsToEliminate)) + "/" + str(len(dimensionsToEliminate)))
    if logEliminations:
        cmd = ["mpirun", "-np", str(numProc), mplrsPath,tmpFilePath+str(len(dimensionsToEliminate)), outFilePath+str(len(dimensionsToEliminate))]
    else:
        cmd = ["mpirun", "-np", str(numProc), mplrsPath, tmpFilePath, outFilePath]
    subprocess.run(cmd)
    
    subprocess.run(["cp",outFilePath+str(len(dimensionsToEliminate)), outFilePath])


def runMarashiWithPolcoV2(stoicMatrix, projectOntoDimensions, outputDir, verbose = True):
    printDatetime("Start", datetime.now())
    sortedStoicMatrix = sortStoicMatrix(stoicMatrix, projectOntoDimensions)
    
    allIndices = list(range(stoicMatrix.shape[1]))
    #indices = list(range(len(projectOntoDimensions)))
    boolStartLoop = True
    counter = 0
    while boolStartLoop:
        counter += 1
        if(len(allIndices) - 5 < len(projectOntoDimensions)):
            allIndices = list(range(len(projectOntoDimensions)))
            boolStopLoop = False
        else:
            allIndices.pop(-1)
            allIndices.pop(-1)
            allIndices.pop(-1)
            allIndices.pop(-1)
            allIndices.pop(-1)
        subOutputFolder = outputDir + "/polco" + str(counter) + "/"
        if not os.path.exists(subOutputFolder):
            os.mkdir(subOutputFolder)
        print(f'Project onto dimensions: {allIndices}')
            
        eliminationMatrix = buildEliminationMatrix(sortedStoicMatrix, allIndices)

        eqFile = 'file.eq'
        iqFileIdentity = 'identity.iq'
        stage = 'projectionCone'
        outputFile = stage + '_out.zip' 
        if verbose:
            print("\nStart Double Description:")
        rayMatrix = doubleDescription(np.identity(eliminationMatrix.shape[0]), iqFileIdentity, subOutputFolder, outputFile, stage, eliminationMatrix.transpose(), eqFile)
        print("Raymatrix - Shape:", rayMatrix.shape)
        printDatetime("Finished first double description step:", datetime.now())
        # rayMatrix = rayMatrix.T
        if verbose:
            print("RayMatrix:\n", rayMatrix)
        # rayMatrix[:,[1,3]] = rayMatrix[:,[3,1]]
        # print("RayMatrix:\n", rayMatrix)
        interestedMatrix = buildInterestedMatrix(sortedStoicMatrix, allIndices)
        if verbose:
            print("InterestedMatrix:\n", interestedMatrix)
        # projectedConeMatrix = buildProjectedConeMatrix(rayMatrix.transpose(), interestedMatrix)
        projectedConeMatrix = -buildProjectedConeMatrix(rayMatrix, interestedMatrix)
        print("projectedConeMatrix - Shape:", projectedConeMatrix.shape)
        # projectedConeMatrix = removeRowsWithZeros(projectedConeMatrix)
        ineFile = os.path.join(subOutputFolder, "projectedCone_H.ine")
        convertMatrixToHRepresentation(projectedConeMatrix, ineFile)
        printDatetime("Finished converting projected Cone Matrix to H-Rep:", datetime.now())

        #projectedConeMatrix[[2,3],:]= projectedConeMatrix[[3,2],:]
        #projectedConeMatrix = removeRedundancy(projectedConeMatrix)

        redund(16,mplrsPath,ineFile, os.path.join(subOutputFolder, "projectedConeMatrix_final_H.ine"))
        printDatetime("Finished redund step:", datetime.now())

        projectedConeMatrix = getMatrixFromHrepresentation(os.path.join(subOutputFolder, "projectedConeMatrix_final_H.ine"))
        if verbose:
            print("\nprojectedConeMatrix:\n", projectedConeMatrix)
            print("\\Shape:\n", projectedConeMatrix.shape)

    eqFileIdentity = 'identity.eq'
    iqFile = 'file.iq'
    stage = 'proCEMs'
    outputFile = stage + '_out.zip'
    if verbose:
        print("\nStart Double Description:")
    proCEMs = doubleDescription(projectedConeMatrix, iqFile, outputDir, outputFile, stage)
    printDatetime("Finished second double description step:", datetime.now())
    # proCEMs= proCEMs.T
    if verbose:
        print("\n#########\nResulting proCEMs:\n", proCEMs)
    efps = calculateEFPsFromProCEMs(proCEMs)
    printDatetime("Finished calculating efps:", datetime.now())
    if verbose:
        print("EFPs:", efps)
    return proCEMs, efps

def runMarashiWithMPLRS(stoicMatrix, projectOntoDimensions, outputDir, mplrsPath, mplrsProcesses = 4, verbose = True):
    printDatetime("Start:", datetime.now())
    sortedStoicMatrix = sortStoicMatrix(stoicMatrix, projectOntoDimensions)
    eliminationMatrix = buildEliminationMatrix(sortedStoicMatrix, projectOntoDimensions)
    
    hEliminationInePath = os.path.join(outputDir, "elimination_H.ine")
    vEliminationInePath = os.path.join(outputDir, "elimination_V.ine")
    convertEqualitiesAndInequalities2hRep(eliminationMatrix.transpose(), np.identity(eliminationMatrix.shape[0]), hEliminationInePath)
    printDatetime("Finished creating eliminiation_H.ine:", datetime.now())
    mplrs_conversion(mplrsProcesses, mplrsPath, hEliminationInePath, vEliminationInePath)
    printDatetime("Finished mplrs conversion:", datetime.now())
    rayMatrix = getRaysFromVrepresentation(vEliminationInePath)
    if verbose:
        print("RayMatrix from first conversion:\n", rayMatrix)

    interestedMatrix = buildInterestedMatrix(sortedStoicMatrix, projectOntoDimensions)
    if verbose:
        print("InterestedMatrix:\n", interestedMatrix)
    projectedConeMatrix = buildProjectedConeMatrix(rayMatrix, interestedMatrix)
    
    hProjectedConeInePath = os.path.join(outputDir, "projectedCone_H.ine")
    convertMatrixToHRepresentation(-projectedConeMatrix, hProjectedConeInePath)
    printDatetime("Finished converting projected Cone Matrix to H-Rep:", datetime.now())

    hRedundProjectedConeInePath = os.path.join(outputDir, "projectedConeMatrix_H_redund.ine")
    redund(mplrsProcesses, mplrsPath, hProjectedConeInePath, hRedundProjectedConeInePath)
    printDatetime("Finished redund step:", datetime.now())

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
