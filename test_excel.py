import numpy as np
import copy
import subprocess
from zipfile import ZipFile 
import re
import os
from datetime import datetime
import time
from projection.projection import runMarashiWithMPLRSSubsets, runMarashiWithPolcoSubsets
from argparse import ArgumentParser, ArgumentTypeError
from projection.logging_config import logger
from projection.marashiProjectionHelpers import get_sorted_column_indices_from_array, logTimeToCSV
import pandas as pd

###
# ----------------
# choose inputs
###

parser = ArgumentParser(description='Generate measurement results with specified parameters.')
parser.add_argument('--tool', choices=['mplrs', 'polco'], default='polco',
                    help='Tool to use for computation (default: polco)')
parser.add_argument('--stepsize', type=int, default=5,
                    help='Step size parameter (default: 5)')
parser.add_argument('-n', type=int, default=120,
                    help='Number of threads to use (default: 120)')
parser.add_argument('--core', default='excel',
                    help='Core name identifier for output files (default: excel)')
parser.add_argument('--outdir', default='testResults/measurements',
                    help='Base output directory path (default: testResults/measurements)')
parser.add_argument('--input', default='input/13015_2011_155_MOESM1_ESM.xls',
                    help='Input excel file used in Marashis Paper')
parser.add_argument('--mplrs', default='mplrs',
                    help='Path to mplrs')
parser.add_argument('--polco', default='polco.jar',
                    help='Path to polco')

# Parse arguments
args = parser.parse_args()

tool = args.tool
stepSize = args.stepsize
numThreads = args.n
start = time.time()
core_name = "/" + args.core
outPath = f"{args.outdir}/{tool}_{args.core}_{numThreads}t_{stepSize}ss"
logTimesFile = outPath + core_name + '_times.csv'
inpExcel = args.input
mplrsPath = args.mplrs
polcoPath = args.polco

if not os.path.exists(inpExcel):
    raise f"Input file {inpExcel} does not exist."

os.makedirs(outPath, exist_ok=True)
time_start = time.time()
logTimeToCSV(logTimesFile, "Start", "Start", 0)
time_initial_setup_start = time.time()

###
# ----------------
# parse excel file
###

parentMatrix = pd.read_excel(inpExcel, index_col=0, sheet_name=1)
for i in range(parentMatrix.shape[1]):  # add backward reactions
    if parentMatrix.iloc[-1][i] == 1.0:
        parentMatrix.insert(
            parentMatrix.shape[1],
            list(parentMatrix.columns)[i].strip()+"_inv",
            -parentMatrix.iloc[:, i]
        )

inpMatrix = parentMatrix.iloc[:-2, :]  # remove footer

# parse reactions of interest
wantedColumns = pd.read_excel(inpExcel, sheet_name=2)
wantedColumnsIndices = []
liststoicMatrix = inpMatrix.columns.tolist()
for i in range(len(liststoicMatrix)):
    saveIndex = False
    for name in wantedColumns.values:
        if (name[0].strip() == liststoicMatrix[i].strip()) or ((name[0].strip()+"_inv") == liststoicMatrix[i].strip()):
            saveIndex = True
            break
    if saveIndex:
        wantedColumnsIndices.append(i)

print(f"Shape Input-Stoichiometric-Matrix: {parentMatrix.shape}")
print(f"Project onto Subnetwork of size: {len(wantedColumnsIndices)} reactions")

inpDims = wantedColumnsIndices
inpMatrix = inpMatrix.to_numpy(dtype=object)

reactions = inpDims
originalReactions = reactions
lenOriginalReactions = len(reactions)
remaining_reactions = [i for i in range(inpMatrix.shape[1]) if i not in reactions]
reactions = reactions + remaining_reactions
inpMatrix = -inpMatrix[:, reactions]
lenCurrentReactions = len(reactions)
time_initial_setup_end = time.time()
logTimeToCSV(logTimesFile, "Initial Setup", "Initial Setup", time_initial_setup_end - time_initial_setup_start)

# split all reversible reactions to match Marashi's approach
reversibleList = np.array([True]*lenCurrentReactions, dtype=bool)

###
# ----------------
# projection
###

iteration = 0
while lenCurrentReactions - stepSize > lenOriginalReactions:
    print(f"Iteration: {lenCurrentReactions}")
    lenCurrentReactions -= stepSize
    if lenCurrentReactions < lenOriginalReactions:
        break
    tempFolder = f"{outPath}/iter_{iteration}/"
    if not os.path.exists(tempFolder):
        os.makedirs(tempFolder)
    sortReactions = get_sorted_column_indices_from_array(inpMatrix, lenOriginalReactions)
    sortReactions = list(range(lenOriginalReactions)) + sortReactions
    print(sortReactions)
    inpMatrix = inpMatrix[:, sortReactions]
    if tool == "mplrs":
        inpMatrix, efps = runMarashiWithMPLRSSubsets(
            inpMatrix, reactions[:lenCurrentReactions], tempFolder, numThreads, True, True, iteration=iteration, 
            originalProjectionReactions=originalReactions, logTimesFile=logTimesFile, reversibleList=reversibleList,
            MPLRS_PATH=mplrsPath
        )
    else:
        inpMatrix, efps = runMarashiWithPolcoSubsets(
            inpMatrix, reactions[:lenCurrentReactions], tempFolder, numThreads, True, True, iteration=iteration, 
            originalProjectionReactions=originalReactions, logTimesFile=logTimesFile, reversibleList=reversibleList,
            MPLRS_PATH=mplrsPath, POLCO_PATH=polcoPath
        )
    iteration += 1
    reversibleList = reversibleList[list(range(lenCurrentReactions))]

# final projection:
tempFolder = f"{outPath}/iter_{iteration}/"
if not os.path.exists(tempFolder):
    os.makedirs(tempFolder)
if tool == "mplrs":
    proCEMs, efps = runMarashiWithMPLRSSubsets(
        inpMatrix, reactions[:lenOriginalReactions], tempFolder, numThreads, False, True,
        iteration=iteration, originalProjectionReactions=originalReactions, logTimesFile=logTimesFile,
        reversibleList=reversibleList, MPLRS_PATH=mplrsPath
    )
else:
    proCEMs, efps = runMarashiWithPolcoSubsets(
        inpMatrix, reactions[:lenOriginalReactions], tempFolder, numThreads, False, True,
        iteration=iteration, originalProjectionReactions=originalReactions, logTimesFile=logTimesFile,
        reversibleList=reversibleList, MPLRS_PATH=mplrsPath, POLCO_PATH=polcoPath
    )

print(f"Number proCEMS: {len(proCEMs)}")
print(f"Number EFPs: {len(efps)}")
time_end = time.time()
print(f"Runtime: {time_end - time_start} seconds")
logTimeToCSV(logTimesFile, "End", "End", time_end)