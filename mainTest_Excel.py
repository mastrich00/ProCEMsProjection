import numpy as np
import copy
import subprocess
from zipfile import ZipFile 
import re
import os
from datetime import datetime
import time
#import efmtool
#from projection.marianneBianca import get_blocked_reactions, rm_reactions, split_all_reversible_reactions, indicate_exchange
from projection.projection import runMarashiWithPolco, runMarashiWithMPLRS, runFELWithPolco, runMarashiWithPolcoIterative, runMarashiWithMPLRSSubsets, runMarashiWithPolcoSubsets
#from gmpy2 import mpq, mpfr
from argparse import ArgumentParser, ArgumentTypeError
from projection.logging_config import logger

polcoPath = "../secondDraft/polco.jar"
mplrsPath = "mplrsV7_2"

import pandas as pd
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
    
parentMatrix = pd.read_excel("13015_2011_155_MOESM1_ESM.xls", index_col=0, sheet_name=1)
for i in range(parentMatrix.shape[1]):
    if parentMatrix.iloc[-1][i] == 1.0:
        parentMatrix.insert(parentMatrix.shape[1], list(parentMatrix.columns)[i].strip()+"_inv", -parentMatrix.iloc[:,i])

inpMatrix = parentMatrix.iloc[:-2,:]

wantedColumns = pd.read_excel("13015_2011_155_MOESM1_ESM.xls", sheet_name=2)
wantedColumnsIndices = []
liststoicMatrix = inpMatrix.columns.tolist()
for i in range(len(liststoicMatrix)):
    saveIndex = False
    for name in wantedColumns.values:
        # print(name[0])
        # print(listSubStoicMatrix[i] + "\n")
        if (name[0].strip() == liststoicMatrix[i].strip()) or ((name[0].strip()+"_inv") == liststoicMatrix[i].strip()):
            saveIndex = True
            break
    if saveIndex:
        wantedColumnsIndices.append(i)
logger.info(f"Shape: {parentMatrix.shape}")        
logger.info(f"NumsProjectionDim: {len(wantedColumnsIndices)}")
parentMatrix.to_excel("out.xlsx", index=False)
#import sys
#sys.exit()
inpDims = wantedColumnsIndices
inpMatrix.to_csv("out_excel.csv", index=False)
inpMatrix = inpMatrix.to_numpy(dtype=object)
print(inpDims)

reactions = inpDims
originalReactions = reactions
lenOriginalReactions = len(reactions)
remaining_reactions = [i for i in range(inpMatrix.shape[1]) if i not in reactions]
reactions = reactions + remaining_reactions
inpMatrix = -inpMatrix[:, reactions]
lenCurrentReactions = len(reactions)
stepSize = 5
iteration = 0
while lenCurrentReactions - stepSize > lenOriginalReactions:
    print(f"Iter: {lenCurrentReactions}")
    lenCurrentReactions -= stepSize
    if lenCurrentReactions < lenOriginalReactions:
        break
    tempFolder = f"testResults/polco_iterative_small/iter_{iteration}/"
    if not os.path.exists(tempFolder):
        os.makedirs(tempFolder)
    sortReactions = get_sorted_column_indices_from_array(inpMatrix, lenOriginalReactions)
    sortReactions = list(range(lenOriginalReactions)) + sortReactions
    print(sortReactions)
    inpMatrix = inpMatrix[:,sortReactions]
    inpMatrix, efps = runMarashiWithMPLRSSubsets(inpMatrix, reactions[:lenCurrentReactions], tempFolder, 20, True, True, iteration=iteration,originalProjectionReactions=originalReactions)
    #inpMatrix = -inpMatrix
    iteration += 1
#exit(0)
tempFolder = f"testResults/polco_iterative_small/iter_{iteration}/"
if not os.path.exists(tempFolder):
    os.makedirs(tempFolder)
proCEMs, efps = runMarashiWithMPLRSSubsets(inpMatrix, reactions[:lenOriginalReactions], tempFolder, 20, False, True, iteration=iteration,originalProjectionReactions=originalReactions)

#proCEMs, efps = runMarashiWithPolcoIterative(inpMatrix, inpDims, "testResults/polco_iterative_excel/", 100, False, True)
# proCEMs, efps = runMarashiWithPolco(inpMatrix, inpDims, "testResults/polco_excel/", 28, False, True)
# proCEMs, efps = runMarashiWithMPLRS(inpMatrix, inpDims, "testResults/mplrs_excel/", mplrsPath, 100, False, True)
# proCEMs, efps = runMarashiWithMPLRS(inpMatrix, inpDims, "testResults/mplrs/", mplrsPath)
# proCEMs, efps = runFELWithPolco(inpMatrix, inpDims, "~/secondDraft/testResults/fel/", mplrsPath)
print(f"proCEMS: {len(proCEMs)}")
print(f"efps: {len(efps)}")