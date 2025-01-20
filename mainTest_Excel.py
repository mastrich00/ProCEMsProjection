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

polcoPath = "polco.jar"
mplrsPath = "mplrsV7_2"

import pandas as pd
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
inpMatrix = inpMatrix.to_numpy()
print(inpDims)

reactions = inpDims
originalReactions = reactions
lenOriginalReactions = len(reactions)
remaining_reactions = [i for i in range(inpMatrix.shape[1]) if i not in reactions]
reactions = reactions + remaining_reactions
inpMatrix = inpMatrix[:, reactions]
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
    inpMatrix, efps = runMarashiWithPolcoSubsets(inpMatrix, reactions[:lenCurrentReactions], tempFolder, 4, True, True, iteration=iteration,originalProjectionReactions=originalReactions)
    #inpMatrix = -inpMatrix
    iteration += 1
#exit(0)
tempFolder = f"testResults/polco_iterative_small/iter_{iteration}/"
if not os.path.exists(tempFolder):
    os.makedirs(tempFolder)
proCEMs, efps = runMarashiWithPolcoSubsets(inpMatrix, reactions[:lenOriginalReactions], tempFolder, 4, False, True, iteration=iteration,originalProjectionReactions=originalReactions)

#proCEMs, efps = runMarashiWithPolcoIterative(inpMatrix, inpDims, "testResults/polco_iterative_excel/", 100, False, True)
# proCEMs, efps = runMarashiWithPolco(inpMatrix, inpDims, "testResults/polco_excel/", 28, False, True)
# proCEMs, efps = runMarashiWithMPLRS(inpMatrix, inpDims, "testResults/mplrs_excel/", mplrsPath, 100, False, True)
# proCEMs, efps = runMarashiWithMPLRS(inpMatrix, inpDims, "testResults/mplrs/", mplrsPath)
# proCEMs, efps = runFELWithPolco(inpMatrix, inpDims, "~/secondDraft/testResults/fel/", mplrsPath)
logger.info(f"proCEMS: {len(proCEMs)}")
logger.info(f"efps: {len(efps)}")