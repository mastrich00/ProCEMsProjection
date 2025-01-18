import cobra
import os

# Load the model from an SBML file
model = cobra.io.read_sbml_model(os.path.join("testResults", "polco_mmsyn_sm00", "mmsyn_sm00.xml"))

# Extract stoichiometric matrix
metabolites = [met.id for met in model.metabolites]
reactions = [rxn.id for rxn in model.reactions]

inpReactions = []
# Create the stoichiometric matrix as a DataFrame
inpMatrix = cobra.util.create_stoichiometric_matrix(model)
for i in range(inpMatrix.shape[1]):
    if "EX_" in model.reactions[i].name[0:4]:
        inpReactions.append(i)

# Convert to DataFrame for better visualization
#import pandas as pd
#stoich_df = pd.DataFrame(stoichiometric_matrix, index=metabolites, columns=reactions)

# Display the matrix
#print(stoich_df)
print(inpReactions)



import numpy as np
import copy
import subprocess
import re
import os
from datetime import datetime
import time
#import efmtool
#from projection.marianneBianca import get_blocked_reactions, rm_reactions, split_all_reversible_reactions, indicate_exchange
from projection.projection import runMarashiWithPolco, runMarashiWithMPLRS, runFELWithPolco
#from gmpy2 import mpq, mpfr
from projection.logging_config import logger

polcoPath = "polco.jar"
mplrsPath = "mplrs"


proCEMs, efps = runMarashiWithPolco(inpMatrix, inpReactions, "testResults/polco_mmsyn_sm00/", 20, False, True)
#proCEMs, efps = runMarashiWithMPLRS(inpMatrix, reactions, "testResults/mplrs/", mplrsPath, 20, False, True)
# proCEMs, efps = runFELWithPolco(inpMatrix, inpDims, "~/secondDraft/testResults/fel/", mplrsPath)
logger.info(f"proCEMS: {len(proCEMs)}")
logger.info(f"efps: {len(efps)}")
