# Marianne + Biancas Code
#import cobra.io
import cobra
import io
from contextlib import redirect_stderr

def read_model(input_filename):
    """
    Reads metabolic model from sbml file using the 'cobra.io.read_sbml_model' functions. Reads io string during reading
    model and catches exchange reactions added by cobrapy. These reactions are going to be removed from the model again.
    :param input_filename: sbml file with cobrapy compatible metabolic model
    :return:
        - model - cobrapy model
        - reas_added - list of reactions added by cobrapy
    """
    
    stream = io.StringIO()
    with redirect_stderr(stream):
        model = cobra.io.read_sbml_model(input_filename)
    console_output = stream.getvalue()

    reas_added = []
    for text in console_output.split('\n'):
        if text.startswith("Adding exchange reaction"):
            tmp = text.split()
            reas_added.append(tmp[3])

    return model, reas_added


def rm_reactions(model, rea_list):
    """
    Removes exchange reaction that were added by cobrapy.
    :param model: cobrapy model
    :param list rea_list: list of reaction names that will be removed
    :return:
        - model - altered cobrapy model
    """
    model.remove_reactions(rea_list)
    return model

def get_blocked_reactions(model):
    rea_list = []
    counter = 0
    for reaction in model.reactions:
        if (reaction.lower_bound == 0 and reaction.upper_bound == 0) or (reaction.lower_bound > reaction.upper_bound):
            rea_list.append(counter)

        counter += 1
    return rea_list

def write_reaction_direction(model):
    '''Loops through all reactions of a given model and writes an attribute "direction" depending on the reaction borders.'''
    for reaction in model.reactions:
        if reaction.lower_bound < 0 and reaction.upper_bound > 0:
            reaction.direction = "both"
        elif reaction.upper_bound > 0 and reaction.lower_bound >= 0:
            reaction.direction = "forward"
        elif reaction.lower_bound < 0 and reaction.upper_bound <= 0:
            reaction.direction = "backward"
        else:
            reaction.direction = False


def correct_stochiometric_matrix(model, smatrix):
    '''Changes algebraic sign of each column in stochiometric matrix with backward reaction.'''
    for i, reaction in enumerate(model.reactions):
        if reaction.direction == "backward":
            smatrix[:,i] = smatrix[:,i] * -1
    return smatrix

def split_all_reversible_reactions(model, projectOntoDimensions: list = None, verbose: bool = True):
    """
    Splitting all reversible reactions into two cobrapy reaction objects.
    """
    counter = 0
    for reaction in model.reactions:
        if reaction.reversibility: # since my splitted reactions are not reversible, they dont get splitted again
            # create backward irreversible reaction from reversible reaction
            backward_reaction = cobra.Reaction(reaction.id + "_b")
            backward_reaction.name = reaction.name # reaction name is the same by purpose (remerging if name is the same)
            backward_reaction.subsystem = reaction.subsystem
            backward_reaction.lower_bound = 0.  # make it irreversible
            backward_reaction.exchange = reaction.exchange
            backward_reaction.notes = reaction.notes
            if 'sbo' in reaction.annotation:
                backward_reaction.annotation['sbo'] = reaction.annotation['sbo']
            
            # add reaction to model
            model.add_reactions([backward_reaction])
            # add metabolites to reaction
            metabolite_dict = reaction.metabolites
            for object in metabolite_dict:
                backward_reaction.add_metabolites({object: (metabolite_dict[object] * -1)})
            # alter forward reaction to split
            reaction.id = reaction.id + "_f"
            #reaction.name = reaction.name
            reaction.lower_bound = 0
            if verbose:
                print("Splitting reaction '" + reaction.id + "'")
            if counter in projectOntoDimensions:
                projectOntoDimensions.append(len(model.reactions)-1)
                if verbose:
                    print("Added reaction '" + reaction.id + "' to dimension to project to")
            
        counter += 1
    return projectOntoDimensions
            
def split_extern_reversible_reactions(model):
    """
    Splitting all reversible extern reactions into two cobrapy reaction objects.
    """
    for reaction in model.reactions:
        if reaction.reversibility and reaction.exchange: # since my splitted reactions are not reversible, they dont get splitted again
            # create backward irreversible reaction from reversible reaction
            backward_reaction = cobra.Reaction(reaction.id + "_b")
            backward_reaction.name = reaction.name # reaction name is the same by purpose (remerging if name is the same)
            backward_reaction.subsystem = reaction.subsystem
            backward_reaction.lower_bound = 0.  # make it irreversible
            backward_reaction.exchange = True
            backward_reaction.notes = reaction.notes
            if 'sbo' in reaction.annotation:
                backward_reaction.annotation['sbo'] = reaction.annotation['sbo']
            
            # add reaction to model
            model.add_reactions([backward_reaction])
            # add metabolites to reaction
            metabolite_dict = reaction.metabolites
            for object in metabolite_dict:
                backward_reaction.add_metabolites({object: (metabolite_dict[object] * -1)})
            # alter forward reaction to split
            reaction.id = reaction.id + "_f"
            #reaction.name = reaction.name
            reaction.lower_bound = 0

def merge_reactions(model, reaction_f, reaction_b):
    """
    Merges two reaction objects together. Rewrites the forward reaction to original state and deletes backward reaction.
    """
    # rewrite forward reaction
    reaction_f.id = reaction_f.id[:-2]
    reaction_f.lower_bound = -1000
    # delete backward reaction
    model.remove_reactions([reaction_b])

def indicate_exchange(model, extern_compartments):
    """Iterates over each reaction in the cobra model and adds an reaction.exchange attribute (value: True or False) to each reaction.
    This attribute is used later in the code to identify exchange reactions."""
    # if no SBO term and products and reactands on reaction --> cobrapy cant detect it as extern
    # we take each reaction as external if:
    # SBO:0000627
    # OR
    # reaction id starts with 'R_EX_' or 'EX_'
    # OR
    # reaction is one sided (reaction.boundary = True) and no demand or sink (we proof for demand and sink before we proof for boundary)
    # if exchange reaction, than reaction.exchange = True
    for reaction in model.reactions:
        # reaction id starts with 'R_EX_' or 'EX_'
        if 'EX_' in reaction.id[:5]:
            print(f'exchange: {reaction.id}')
            reaction.exchange = True
            continue

        # unfortunately, demand reactions do have the same SBO Term as exchange reactions
        # demand reactions differ from exchange reactions in case of exchange reactions work with metabolites of external compartment
        # we could use that (if exchange sbo term --> proof if one metabolite is in external compartment):
        if 'sbo' in reaction.annotation:
            if reaction.annotation['sbo'] == 'SBO:0000627':
                extern = False
                for metabolite in reaction.metabolites:
                    # compartment name needs to be given and contain 'ex' or compartment id needs to contain 'e'
                    if (metabolite.compartment in extern_compartments and model.compartments[metabolite.compartment]) or 'e' in metabolite.compartment:
                        extern = True
                if extern:
                    # add reaction to exchange reactions
                    print(f'exchange: {reaction.id}')
                    reaction.exchange = True
                    continue

        # kill all demand reactions by id prefix 'R_DM_'
        if 'DM_' in reaction.id[:5]:
            print(f'demand: {reaction.id}')
            reaction.exchange = False
            continue
        # kill all sink reactions by SBO Term
        if 'sbo' in reaction.annotation:
            if reaction.annotation['sbo'] == 'SBO:0000632':
                print(f'sink: {reaction.id}')
                reaction.exchange = False
                continue
        # kill all sink reactions by id prefix 'R_SINK_'
        if 'SINK_' in reaction.id[:6]:
            print(f'sink: {reaction.id}')
            reaction.exchange = False
            continue

        # proof if one sided
        if reaction.boundary:
            # since we already ended iteration if criterium was for sink or demand, all remaining reaction.boundary = True reactions hopefully are exchange reactions
            print(f'exchange: {reaction.id}')
            reaction.exchange = True
            continue
        # if no exchange reaction (iteration reaches this point)
        reaction.exchange = False
    
    # create an input reaction counter, which is needed for the column names of the flux vertices dataframe
    input_reaction_counter = {}
    counter = 0
    for reaction in model.reactions:
        if reaction.exchange == True:
            counter += 1
            input_reaction_counter[reaction.id] = counter

    return input_reaction_counter