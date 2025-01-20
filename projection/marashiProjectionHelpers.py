import numpy as np
import copy
from projection.logging_config import logger

def sortStoicMatrix(stoicMatrix: np.matrix, reactions: list, verbose=True):
    '''
    sorts the stoichometric matrix so that the columns specified with indices in ${reactions} are the first columns
    '''
    logger.info("Sorting Stoic Matrix")
    j = 0
    sortedReactions = copy.deepcopy(reactions)
    sortedReactions.sort()
    logger.info(sortedReactions)
    while len(sortedReactions) > 0:
        i = sortedReactions.pop(0)
        if(not (j in sortedReactions)):
            if verbose:
                logger.info(f"Swapping col {str(i)} with {str(j)}")
                # print(f"Swapping col {str(i)} with {str(j)}")
            stoicMatrix[:, [j, i]] = stoicMatrix[:, [i, j]]
        j += 1

    if verbose:
        logger.info("Sorted StoicMatrix:")
        logger.info(stoicMatrix)
    return stoicMatrix

def buildInterestedMatrix(stoicMatrix: np.matrix, reactions: list, verbose=True):
    '''
    returns matrix G according to Marashi paper
    '''
    p = len(reactions) 
    q = stoicMatrix.shape[1]-p
    if( p + q != stoicMatrix.shape[1]):
        raise Exception("Dimensions of columns are incorrect")
    A = stoicMatrix[:,:p]
    if verbose:
        logger.info("A:")
        logger.info(A)
    interestedMatrix = np.block([[-A],[A],[-np.identity(p)],[np.zeros((q,p))]])
    if verbose:
        logger.info("InterestedMatrix:")
        logger.info(interestedMatrix)
        logger.info(f"Dimensions: {interestedMatrix.shape}")
    return interestedMatrix

def buildEliminationMatrix(stoicMatrix: np.matrix, reactions: list, verbose=True):
    '''
    returns matrix H according to Marashi paper
    '''

    p = len(reactions)
    q = stoicMatrix.shape[1]-p
    print(f"elMatrix: p={p}, q={q}")
    if( p + q != stoicMatrix.shape[1]):
        raise Exception("Dimensions of columns are incorrect")
    B = stoicMatrix[:,p:]
    if verbose:
        logger.info("B:")
        logger.info(B)
    eliminationMatrix = np.block([[-B],[B],[np.zeros((p,q))],[-np.identity(q)]])
    if verbose:
        logger.info("EliminationMatrix:")
        logger.info(eliminationMatrix)
        logger.info(f"Dimensions: {eliminationMatrix.shape}")
    return eliminationMatrix
        

def buildProjectedConeMatrix(rayMatrix: np.matrix, interestedMatrix: np.matrix, verbose = True):
    '''
    multiply R and G (Marashi paper)
    '''
    projectedConeMatrix = np.matmul(rayMatrix, interestedMatrix)
    if(verbose):
        logger.info("ProjectedConeMatrix:")
        logger.info(projectedConeMatrix)
    return projectedConeMatrix

def getSupport(vector):
    '''
    returns the support of a vector
    '''
    support = []
    index = 0
    for entry in vector:
        if abs(entry) > 10**-8:
            support.append(index)
        index += 1
    return support    

def calculateEFPsFromProCEMs(proCEMs: np.matrix, verbose=True):
    '''
    calculates the EFPs from the ProCEMs according to Marashi paper
    (not really tested/verified by me yet)
    '''
    # print(proCEMs)
    efps = []
    for colIndex in range(proCEMs.shape[1]):
        proCEM = proCEMs[:,colIndex]
        # print(f"{colIndex}: {proCEM}")
        supportProCEM = getSupport(proCEM)
        # print(supportProCEM)
        Z = copy.deepcopy(supportProCEM)
        for compareColIndex in range(proCEMs.shape[1]):
            if compareColIndex == colIndex:
                continue
            compareProCEM = proCEMs[:,compareColIndex]
            compareSupport = getSupport(compareProCEM)
            # print(Z)
            # print(compareSupport)
            # print("")
            if set(compareSupport).issubset(supportProCEM) and compareSupport != supportProCEM:
                for entry in compareSupport:
                    if entry in Z:
                        Z.remove(entry)
                
        if len(Z) > 0:
            logger.info(supportProCEM)
            efps.append([supportProCEM])
            
    return efps