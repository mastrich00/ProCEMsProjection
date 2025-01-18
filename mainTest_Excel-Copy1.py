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
from projection.projection import runMarashiWithPolco, runMarashiWithMPLRS, runFELWithPolco
#from gmpy2 import mpq, mpfr
from argparse import ArgumentParser, ArgumentTypeError
from projection.logging_config import logger
from projection.redund import redundIterative, redund

polcoPath = "polco.jar"
mplrsPath = "mplrs"

import pandas as pd

chunkSize = 10000
logger.info(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} with chunkSize {chunkSize}")
# redundIterative(17, mplrsPath, "testResults/polco_excel/projectedCone_H.ine", "testResults/polco_excel/testOutprojectedCone_H.ine", "testResults/polco_excel/temp", chunkSize=chunkSize, verbose=True)
redund(17, mplrsPath, "testResults/polco_excel/testOutprojectedCone_H.ine", "testResults/polco_excel/out_testOutprojectedCone_H.ine", verbose=True)
logger.info(f"End: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
