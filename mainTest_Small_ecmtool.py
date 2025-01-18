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

polcoPath = "polco.jar"
mplrsPath = "mplrs"

inpMatrix = np.matrix([[1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0],
                         [0,1,0,0,-1,1,-1,0,0,0,0,0,0,0],
                         [0,0,1,0,1,0,0,-1,0,0,0,0,-4,0],
                         [0,0,0,1,0,0,0,-1,-3,0,0,0,-1,0],
                         [0,0,0,0,0,-1,1,1,0,0,1,0,-4,0],
                         [0,0,0,0,0,0,0,0,2,-1,0,0,0,0],
                         [0,0,0,0,0,0,0,0,0,0,0,1,-1,0],
                         [0,0,0,0,0,0,0,1,1,0,-1,-1,0,-1]])
reactions = [0,9,12,13] # interested in R1, R10, Rbio, Pex

import libsbml

def create_sbml_from_matrix(stoichiometric_matrix, species_names, reaction_names, output_file):
    # Create an SBML Document
    document = libsbml.SBMLDocument(3, 1)  # SBML Level 3 Version 1
    model = document.createModel()
    model.setId("StoichiometricModel")
    # Add a default compartment
    compartment = model.createCompartment()
    compartment.setId("default")
    compartment.setSize(1.0)
    compartment.setConstant(True)


    # Add Species
    for species in species_names:
        s = model.createSpecies()
        s.setId(species)
        s.setCompartment("default")
        s.setBoundaryCondition(False)
        s.setInitialAmount(0.0)
        s.setHasOnlySubstanceUnits(True)
        s.setConstant(False)

    # Add Reactions and Flux Bounds as Parameters
    for j, reaction in enumerate(reaction_names):
        r = model.createReaction()
        r.setId(reaction)
        r.setReversible(False)
        r.setFast(False)

        # Add flux bounds as parameters directly to the model
        lower_bound = model.createParameter()
        lower_bound.setId(f"LOWER_BOUND_{reaction}")
        lower_bound.setValue(0.0)
        lower_bound.setUnits("mmol_per_gDW_per_hr")
        lower_bound.setConstant(True)

        upper_bound = model.createParameter()
        upper_bound.setId(f"UPPER_BOUND_{reaction}")
        upper_bound.setValue(1000.0)
        upper_bound.setUnits("mmol_per_gDW_per_hr")
        upper_bound.setConstant(True)

        # Add Reactants and Products
        for i, species in enumerate(species_names):  # Loop through species
            stoich = stoichiometric_matrix[i,j]  # Extract the single stoichiometry value
            if stoich < 0:  # Reactant
                sr = r.createReactant()
                sr.setSpecies(species)
                sr.setStoichiometry(float(-stoich))  # Convert to float
                sr.setConstant(True)
            elif stoich > 0:  # Product
                sp = r.createProduct()
                sp.setSpecies(species)
                sp.setStoichiometry(float(stoich))  # Convert to float
                sp.setConstant(True)
                

    # Add Reactions
    for j, reaction in enumerate(reaction_names):  # Loop through reactions
        r = model.createReaction()
        r.setId(reaction)
        r.setReversible(True)
        r.setFast(False)

        # Add Reactants and Products
        for i, species in enumerate(species_names):  # Loop through species
            stoich = stoichiometric_matrix[i,j]  # Extract the single stoichiometry value
            if stoich < 0:  # Reactant
                sr = r.createReactant()
                sr.setSpecies(species)
                sr.setStoichiometry(float(stoich))  # Convert to float
                sr.setConstant(True)
            elif stoich > 0:  # Product
                sp = r.createProduct()
                sp.setSpecies(species)
                sp.setStoichiometry(float(stoich))  # Convert to float
                sp.setConstant(True)


    # Write SBML to File
    libsbml.writeSBMLToFile(document, output_file)
    print(f"SBML file saved to {output_file}")



# Define species and reaction names
species_names = ["A", "B", "C", "D", "E", "F", "G", "H"]
reaction_names = ["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13", "R14"]

# Convert to SBML
create_sbml_from_matrix(inpMatrix, species_names, reaction_names, "model.xml")

# Validate the SBML file 
doc = libsbml.readSBMLFromFile("model.xml") 
if doc.getNumErrors() > 0: 
    print("SBML validation errors:") 
    print(doc.getErrorLog().toString())
else: 
    print("SBML file is valid.")