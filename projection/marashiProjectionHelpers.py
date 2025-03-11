import numpy as np
import copy
from projection.logging_config import logger
import datetime
from fractions import Fraction

def logTimeToCSV(filepath, step, substep, time, date=None):
    if not date:
        date = datetime.datetime.now()
    with open(filepath, "a") as file:
        file.write(f"{step},{substep},{time},{str(date)}\n")

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

def fraction_matmul(A, B):
    """Perform matrix multiplication while keeping Fraction precision."""
    return np.array([[sum(a * b for a, b in zip(row, col)) 
                      for col in zip(*B)] 
                     for row in A], dtype=object)
    
def buildProjectedConeMatrix(rayMatrix, interestedMatrix, verbose = True):
    '''
    multiply R and G (Marashi paper)
    '''
    # projectedConeMatrix = np.matmul(rayMatrix.astype(float), interestedMatrix.astype(float))
    # projectedConeMatrix = float_to_fraction_matrix(projectedConeMatrix)
    projectedConeMatrix = fraction_matmul(rayMatrix, interestedMatrix)
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
    '''
    efps = []
    for rowIndex in range(proCEMs.shape[0]):
        proCEM = proCEMs[rowIndex, :]
        supportProCEM = getSupport(proCEM)
        Z = copy.deepcopy(supportProCEM)
        for compareRowIndex in range(proCEMs.shape[0]):
            if compareRowIndex == rowIndex:
                continue
            compareProCEM = proCEMs[compareRowIndex,:]
            compareSupport = getSupport(compareProCEM)
            if set(compareSupport).issubset(supportProCEM) and compareSupport != supportProCEM:
                for entry in compareSupport:
                    if entry in Z:
                        Z.remove(entry)
                
        if len(Z) > 0:
            logger.info(supportProCEM)
            efps.append([supportProCEM])
      
    return efps