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
from projection.projection import runMarashiWithPolco, runMarashiWithMPLRS, runFELWithPolco, runMarashiWithPolcoIterative, runMarashiWithPolcoSubsets, runMarashiWithMPLRSSubsets
#from gmpy2 import mpq, mpfr
from argparse import ArgumentParser, ArgumentTypeError
from projection.logging_config import logger

polcoPath = "polco.jar"
mplrsPath = "mplrsV7_2"

inpMatrix = np.matrix([[1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0],
                         [0,1,0,0,-1,1,-1,0,0,0,0,0,0,0],
                         [0,0,1,0,1,0,0,-1,0,0,0,0,-4,0],
                         [0,0,0,1,0,0,0,-1,-3,0,0,0,-1,0],
                         [0,0,0,0,0,-1,1,1,0,0,1,0,-4,0],
                         [0,0,0,0,0,0,0,0,2,-1,0,0,0,0],
                         [0,0,0,0,0,0,0,0,0,0,0,1,-1,0],
                         [0,0,0,0,0,0,0,1,1,0,-1,-1,0,-1]])
reactions = [0,9,12,13] # interested in R1, R10, Rbio, Pex
originalReactions = reactions
lenOriginalReactions = len(reactions)
remaining_reactions = [i for i in range(inpMatrix.shape[1]) if i not in reactions]
reactions = reactions + remaining_reactions
inpMatrix = -inpMatrix[:, reactions]
lenCurrentReactions = len(reactions)
stepSize = 2
iteration = 0
while lenCurrentReactions - stepSize > lenOriginalReactions:
    print(f"Iter: {lenCurrentReactions}")
    lenCurrentReactions -= stepSize
    if lenCurrentReactions < lenOriginalReactions:
        break
    inpMatrix, efps = runMarashiWithPolcoSubsets(inpMatrix, reactions[:lenCurrentReactions], "testResults/polco_iterative_small/", 100, True, True, iteration=iteration)
    iteration += 1
#exit(0)
#proCEMs, efps = runMarashiWithMPLRSSubsets(inpMatrix, reactions[:lenOriginalReactions], "testResults/polco_iterative_small/", 100, False, True, iteration=iteration)
proCEMs, efps = runMarashiWithPolcoSubsets(inpMatrix, reactions[:lenOriginalReactions], "testResults/polco_iterative_small/", 100, False, True, iteration=iteration)
#proCEMs, efps = runMarashiWithPolcoIterative(inpMatrix, originalReactions, "testResults/polco_iterative_small/", 100, False, True)
#proCEMs, efps = runMarashiWithPolco(inpMatrix, reactions, "testResults/polco_small/", 17, False, True)
#proCEMs, efps = runMarashiWithMPLRS(inpMatrix, reactions, "testResults/mplrs/", mplrsPath, 20, False, True)
# proCEMs, efps = runFELWithPolco(inpMatrix, inpDims, "~/secondDraft/testResults/fel/", mplrsPath)
logger.info(f"proCEMS: {len(proCEMs)}")
logger.info(f"efps: {len(efps)}")
