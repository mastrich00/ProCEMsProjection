import numpy as np
import sympy as sp
import copy
import subprocess
from zipfile import ZipFile 
import re
import os
from datetime import datetime
import time
#import efmtool
from projection.marianneBianca import get_blocked_reactions, rm_reactions, split_all_reversible_reactions, indicate_exchange
from projection.projection import runMarashiWithPolco, runMarashiWithMPLRS, runFELWithPolco
from gmpy2 import mpq, mpfr
from argparse import ArgumentParser, ArgumentTypeError

from pathlib import Path
home = str(Path('~').expanduser())
projectFolder = os.path.join(home, "secondDraft")

polcoPath = os.path.join(projectFolder,"polco.jar")
mplrsPath = "mplrs"
#polcoPath = "/home/mstasek/Downloads/polco.jar"

import cobra.io
import cobra
#import ecmtool
# model = cobra.io.sbml.read_sbml_model("/home/martin/projects/python-projects/secondDraft/ecmtool/models/e_coli_core.xml")
#model = cobra.io.sbml.read_sbml_model("~/secondDraft/e_coli_core_modified.xml")
# model = cobra.io.sbml.read_sbml_model("/home/martin/projects/python-projects/testmodel_v3.sbml")
# model = cobra.io.sbml.read_sbml_model("/home/martin/projects/python-projects/secondDraft/ecolicore_ucsd.xml")
modelPath = os.path.join(projectFolder, "e_coli_core_modified.xml")
# modelPath = os.path.join(projectFolder, "testmodel_v3.sbml")
model = cobra.io.sbml.read_sbml_model(modelPath)

blocked_reactions = get_blocked_reactions(model)
rm_reactions(model, blocked_reactions)

input_reaction_counter = indicate_exchange(model, "e")
projectOntoDimensions = []
counter = 0
for reaction in model.reactions:
    # print(vars(reaction))
    #if "EX_" in reaction.id[:5]:
    if reaction.exchange or ("RBIOMASS" == reaction.id) or ("BIOMASS" in reaction.id):
        if reaction.id != "EX_RBIOMASS":
            projectOntoDimensions.append(counter)
    counter += 1
# for reaction in model.reactions:
#     print(vars(reaction))
#     #if "EX_" in reaction.id[:5]:
#     if reaction.exchange:
#         projectOntoDimensions.append(counter)
#     counter += 1

print(projectOntoDimensions)
for i in projectOntoDimensions:
    print(f"{i}: {model.reactions[i].id}")
print(len(projectOntoDimensions))
# exit()
#projectOntoDimensions.pop(0)
inpDims = split_all_reversible_reactions(model, projectOntoDimensions)
inpMatrix = cobra.util.array.create_stoichiometric_matrix(model, 'dense')
print(inpMatrix.shape)
print(len(inpDims))

# import sys
# sys.path.append(os.path.abspath("/home/martin/projects/python-projects/ecmtool"))

# from ecmtool import mpi_wrapper
# from ecmtool.conversion_cone import calculate_linearities, calc_C0_dual_extreme_rays, calc_H, calc_C_extreme_rays, post_process_rays
# from ecmtool.helpers import redund, to_fractions, prep_mplrs_input, execute_mplrs, process_mplrs_output
# from ecmtool.helpers import mp_print, unsplit_metabolites, process_all_from_mplrs, print_ecms_direct, normalize_columns, uniqueReadWrite, save_and_print_ecms
# from ecmtool.intersect_directly_mpi import intersect_directly
# from ecmtool.network import add_reaction_tags, Network#, extract_sbml_stoichiometry
# from ecmtool.ecmtool.network import extract_sbml_stoichiometry
# # def preprocess_sbml(model_path, args_hide='', args_add_objective_metabolite=True, args_auto_direction=True, args_use_external_compartment=None, args_prohibit='', args_tag='', 
# #                     args_print_reactions=False, args_print_metabolites=True, args_compress=True, args_only_rays=False, args_hide_all_in_or_outputs='', args_verbose=True, args_scei=True,
# #                     args_remove_infeasible=True):

# #     network = extract_sbml_stoichiometry(model_path, add_objective=args_add_objective_metabolite,
# #                                          determine_inputs_outputs=args_auto_direction,
# #                                          skip_external_reactions=True,
# #                                          use_external_compartment=args_use_external_compartment)

# #     # tagged_reaction_indices = []
# #     external_cycles = None

# #     # adj = get_metabolite_adjacency(network.N)

# #     #if not args.auto_direction:
# #     #    set_inoutputs(args.inputs, args.outputs, network)

# #     if args_hide:
# #         hide_indices = [int(index) for index in args_hide.split(',') if len(index)]
# #         network.hide(hide_indices)

# #     if args_prohibit:
# #         prohibit_indices = [int(index) for index in args_prohibit.split(',') if len(index)]
# #         network.prohibit(prohibit_indices)

# #     tag_ids = []
# #     if args_tag:
# #         tagged_reaction_indices = [int(index) for index in args_tag.split(',') if len(index)]
# #         tag_ids = add_reaction_tags(network, tagged_reaction_indices)

# #     if args_print_reactions:
# #         mp_print('Reactions%s:' % (' before compression' if args_compress else ''))
# #         for index, item in enumerate(network.reactions):
# #             mp_print(index, item.id, item.name, 'reversible' if item.reversible else 'irreversible')

# #     if args_print_metabolites:
# #         mp_print('Metabolites%s:' % (' before compression' if args_compress else ''))
# #         for index, item in enumerate(network.metabolites):
# #             mp_print(index, item.id, item.name, 'external' if item.is_external else 'internal', item.direction)

# #     orig_ids = [m.id for m in network.metabolites]
# #     orig_N = network.N

# #     # if args.direct:
# #     #     # Initialise mpi4py only here, because it can not be started when using mplrs due to
# #     #     # only being able to run one instance at a time, and mpi4py creates an instance on import.
# #     #     mpi_wrapper.mpi_init(mplrs_present=mplrs_present)

# #     #     from ecmtool.intersect_directly_mpi import intersect_directly, remove_cycles, \
# #     #         compress_after_cycle_removing, check_if_intermediate_cone_exists

# #     #     # Check if intermediate cone exists at the given location
# #     #     if args.intermediate_cone_path:
# #     #         check_if_intermediate_cone_exists(args.intermediate_cone_path)

# #     # Split metabolites in input and output
# #     network.split_in_out(args_only_rays)
# #     cycle_removal_boolean = True if not args_only_rays else False

# #     if args_hide_all_in_or_outputs:
# #         hide_indices = [ind for ind, metab in enumerate(network.metabolites) if
# #                         (metab.is_external) & (metab.direction == args_hide_all_in_or_outputs) & (
# #                             not metab.id == 'objective_virtout') & (
# #                                 metab.id.replace("_virtin", "").replace("_virtout", "") not in tag_ids)]
# #         network.hide(hide_indices)

# #     if args_compress:
# #         network.compress(verbose=args_verbose, SCEI=args_scei, cycle_removal=cycle_removal_boolean,
# #                          remove_infeasible=args_remove_infeasible)

# #     # if args.direct:
# #     #     network.split_reversible()
# #     #     network.N = np.transpose(redund(np.transpose(network.N)))

# #     #     R, network, external_cycles = remove_cycles(network.N, network)
# #     #     n_reac_according_to_N = network.N.shape[1]
# #     #     removable_reacs = np.arange(n_reac_according_to_N, len(network.reactions))
# #     #     network.drop_reactions(removable_reacs)
# #     #     network = compress_after_cycle_removing(network)

# #     if args_print_reactions and args_compress:
# #         mp_print('Reactions (after compression):')
# #         for index, item in enumerate(network.reactions):
# #             mp_print(index, item.id, item.name, 'reversible' if item.reversible else 'irreversible')

# #     if args_print_metabolites and args_compress:
# #         mp_print('Metabolites (after compression):')
# #         for index, item in enumerate(network.metabolites):
# #             mp_print(index, item.id, item.name, 'external' if item.is_external else 'internal', item.direction)

# #     return network, external_cycles

# def preprocess_sbml(args):
#     model_path = args.model_path

#     network = extract_sbml_stoichiometry(model_path, add_objective=args.add_objective_metabolite,
#                                          determine_inputs_outputs=args.auto_direction,
#                                          skip_external_reactions=True,
#                                          use_external_compartment=args.use_external_compartment)

#     # tagged_reaction_indices = []
#     external_cycles = None

#     # adj = get_metabolite_adjacency(network.N)

#     if not args.auto_direction:
#         set_inoutputs(args.inputs, args.outputs, network)

#     if args.hide:
#         hide_indices = [int(index) for index in args.hide.split(',') if len(index)]
#         network.hide(hide_indices)

#     if args.prohibit:
#         prohibit_indices = [int(index) for index in args.prohibit.split(',') if len(index)]
#         network.prohibit(prohibit_indices)

#     tag_ids = []
#     if args.tag:
#         tagged_reaction_indices = [int(index) for index in args.tag.split(',') if len(index)]
#         tag_ids = add_reaction_tags(network, tagged_reaction_indices)

#     if args.print_reactions:
#         mp_print('Reactions%s:' % (' before compression' if args.compress else ''))
#         for index, item in enumerate(network.reactions):
#             mp_print(index, item.id, item.name, 'reversible' if item.reversible else 'irreversible')

#     if args.print_metabolites:
#         mp_print('Metabolites%s:' % (' before compression' if args.compress else ''))
#         for index, item in enumerate(network.metabolites):
#             mp_print(index, item.id, item.name, 'external' if item.is_external else 'internal', item.direction)

#     orig_ids = [m.id for m in network.metabolites]
#     orig_N = network.N

#     if args.direct:
#         # Initialise mpi4py only here, because it can not be started when using mplrs due to
#         # only being able to run one instance at a time, and mpi4py creates an instance on import.
#         mpi_wrapper.mpi_init(mplrs_present=mplrs_present)

#         from ecmtool.intersect_directly_mpi import intersect_directly, remove_cycles, \
#             compress_after_cycle_removing, check_if_intermediate_cone_exists

#         # Check if intermediate cone exists at the given location
#         if args.intermediate_cone_path:
#             check_if_intermediate_cone_exists(args.intermediate_cone_path)

#     # Split metabolites in input and output
#     network.split_in_out(args.only_rays)
#     cycle_removal_boolean = True if not args.only_rays else False

#     if args.hide_all_in_or_outputs:
#         hide_indices = [ind for ind, metab in enumerate(network.metabolites) if
#                         (metab.is_external) & (metab.direction == args.hide_all_in_or_outputs) & (
#                             not metab.id == 'objective_virtout') & (
#                                 metab.id.replace("_virtin", "").replace("_virtout", "") not in tag_ids)]
#         network.hide(hide_indices)

#     if args.compress:
#         network.compress(verbose=args.verbose, SCEI=args.scei, cycle_removal=cycle_removal_boolean,
#                          remove_infeasible=args.remove_infeasible)

#     if args.direct:
#         network.split_reversible()
#         network.N = np.transpose(redund(np.transpose(network.N)))

#         R, network, external_cycles = remove_cycles(network.N, network)
#         n_reac_according_to_N = network.N.shape[1]
#         removable_reacs = np.arange(n_reac_according_to_N, len(network.reactions))
#         network.drop_reactions(removable_reacs)
#         network = compress_after_cycle_removing(network)

#     if args.print_reactions and args.compress:
#         mp_print('Reactions (after compression):')
#         for index, item in enumerate(network.reactions):
#             mp_print(index, item.id, item.name, 'reversible' if item.reversible else 'irreversible')

#     if args.print_metabolites and args.compress:
#         mp_print('Metabolites (after compression):')
#         for index, item in enumerate(network.metabolites):
#             mp_print(index, item.id, item.name, 'external' if item.is_external else 'internal', item.direction)

#     return network, external_cycles

# parser = ArgumentParser(
#         description='Calculate Elementary Conversion Modes from an SBML model. For medium-to large networks, be sure to define --inputs and --outputs. This reduces the enumeration problem complexity considerably.')
#     # Choosing what part to run. Not putting this first argument will just run the whole program.
# parser.add_argument('command', nargs='?', default='all',
#                         help='Optional: run only a single step of ecmtool, continuing from the state of the previous step. \n'
#                              'Allowed values (in order of execution): preprocess, direct_intersect (only when --direct true),\n'
#                              'calc_linearities, prep_C0_rays, all_until_mplrs, calc_C0_rays, all_between_mplrs, process_C0_rays, calc_H, prep_C_rays, calc_C_rays,\n'
#                              'process_C_rays, all_from_mplrs, postprocess, save_ecms, make_unique. Omit to run all steps.')

# # Choices for the model on which ECM-calculation is performed
# parser.add_argument('--model_path', type=str, default='ecmtool/models/e_coli_core.xml',
#                     help='Relative or absolute path to an SBML model .xml file')
# parser.add_argument('--add_objective_metabolite', type=bool, default=True,
#                     help='Add a virtual metabolite containing the stoichiometry of the objective function of the model (default: true)')
# parser.add_argument('--use_external_compartment', type=str, default=None,
#                     help='If a string is given, this string indicates how the external compartment in metabolite_ids of SBML-file is marked. By default, dead-end reaction-detection is used to find external metabolites, and no compartment-information. Please check if external compartment detection works by checking metabolite information before compression and with --primt metabolites true')
# parser.add_argument('--auto_direction', type=bool, default=True,
#                     help='Automatically determine external metabolites that can only be consumed or produced (default: true)')
# parser.add_argument('--inputs', type=str, default='',
#                     help='Comma-separated list of external metabolite indices, as given by --print_metabolites true (before compression), that can only be consumed')
# parser.add_argument('--outputs', type=str, default='',
#                     help='Comma-separated list of external metabolite indices, as given by --print_metabolites true (before compression), that can only be produced. '
#                             'If inputs are given, but no outputs, then everything not marked as input is marked as output.'
#                             'If inputs and outputs are given, the possible remainder of external metabolites is marked as both')
# parser.add_argument('--hide', type=str, default='',
#                     help='Comma-separated list of external metabolite indices, as given by --print_metabolites true (before compression), that are transformed into internal metabolites by adding bidirectional exchange reactions')
# parser.add_argument('--prohibit', type=str, default='',
#                     help='Comma-separated list of external metabolite indices, as given by --print_metabolites true (before compression), that are transformed into internal metabolites without adding bidirectional exchange reactions.'
#                             'This metabolite can therefore not be used as input nor output.')
# parser.add_argument('--tag', type=str, default='',
#                     help='Comma-separated list of reaction indices, as given by --print_reactions true (before compression), that will be tagged with new virtual metabolites, such that the reaction flux appears in ECMs.')
# parser.add_argument('--hide_all_in_or_outputs', type=str, default='',
#                     help='String that is either empty, input, or output. If it is inputs or outputs, after splitting metabolites, all inputs or outputs are hidden (objective is always excluded)')
# parser.add_argument('--only_rays', type=bool, default=False,
#                     help='Enable to only return extreme rays, and not elementary modes. This describes the full conversion space, but not all biologically relevant minimal conversions. See (Urbanczik, 2005) (default: false)')

# # Choices for the output
# parser.add_argument('--output_fractions', type=bool, default=False,
#                     help='Determines whether fractions or their approximating floats are stored as outputs.')
# parser.add_argument('--out_path', default='conversion_cone.csv',
#                     help='Relative or absolute path to the .csv file you want to save the calculated conversions to (default: conversion_cone.csv)')
# parser.add_argument('--make_unique', type=bool, default=False,
#                     help='Make sure set of ECMs is unique at the end  (default: False). Setting this to false '
#                             'drastically reduces memory requirements of '
#                             'the postprocessing. When running with direct-method or polco, uniqueness is already '
#                             'guaranteed, but mplrs is known to sometimes create duplicates. One can always make a'
#                             'result unique by running main.py make_unique --out_path <path to calculated conversions>')

# # Choices for the way in which ecm-calculation is performed. This can have large consequences for computational
# # time. See publications for information.
# parser.add_argument('--direct', type=bool, default=False,
#                     help='Enable to intersect with equalities directly. Direct intersection works better than indirect when many metabolites are hidden, and on large networks (default: False)')
# parser.add_argument('--compress', type=bool, default=True,
#                     help='Perform compression to which the conversions are invariant, and reduce the network size considerably (default: True)')
# parser.add_argument('--remove_infeasible', type=bool, default=True,
#                     help='Remove reactions that cannot carry flux during compression. Switch off when this gives rise to numerical linear algebra problems. (default: True)')
# parser.add_argument('--redund_after_polco', type=bool, default=True,
#                     help='(Indirect intersection only) Enables redundant row removal from inequality description of dual cone. Works well with models with relatively many internal metabolites, and when running parrallelized computation using MPI (default: true)')
# parser.add_argument('--polco', type=bool, default=False,
#                     help='Uses polco instead of mplrs for extreme ray enumeration (default: false)')
# parser.add_argument('--processes', type=int, default=3,
#                     help='Numer of processes for calculations (default: 3 - minimum required for mplrs)')
# parser.add_argument('--jvm_mem', type=int, default=None, nargs='*', action='store',
#                     help='Two values given the minimum and maximum memeory for java machine in GB e.g. 50 300 (default: maximum memory available)')
# parser.add_argument('--path2mplrs', type=str, default=None,
#                     help='if mplrs binary is not accessable via PATH variable "mplrs", the absolute path to the binary can be provided with "--path2mplrs" e.g. "--path2mplrs /home/user/mplrs/lrslib-071b/mplrs" ')
# parser.add_argument('--scei', type=bool, default=True, help='Enable to use SCEI compression (default: true)')
# parser.add_argument('--sort_order', type=str, default='min_adj',
#                     help='Order in which internal metabolites should be set to zero (in direct enumeration). Default is to minimize the added adjacencies, other options are: min_lp, max_lp_per_adj, min_connections')
# parser.add_argument('--intermediate_cone_path', type=str, default='',
#                     help='Filename where intermediate cone result can be found. If an empty string is given (default), then no intermediate result is picked up and the calculation is done in full')
# parser.add_argument('--manual_override', type=str, default='',
#                     help='Index indicating which metabolite should be intersected in first step. Advanced option, can be used in combination with --intermediate_cone_path, to pick a specific intersection at a specific time.')

# # Choices on what is printed in the console
# parser.add_argument('--print_metabolites', type=bool, default=True,
#                     help='Print the names and IDs of metabolites in the (compressed) metabolic network (default: true)')
# parser.add_argument('--print_reactions', type=bool, default=False,
#                     help='Print the names and IDs of reactions in the (compressed) metabolic network (default: true)')
# parser.add_argument('--print_conversions', type=bool, default=False,
#                     help='Print the calculated conversion modes (default: false)')
# parser.add_argument('--verbose', type=bool, default=True,
#                     help='Enable to show detailed console output (default: true)')
# parser.add_argument('--timestamp', type=bool, default=True,
#                     help='Determines whether we print timestamps for several steps in the program.')
# args = parser.parse_args()
# #model_path = "/home/martin/projects/python-projects/secondDraft/ecmtool/models/e_coli_core.xml"

os.chdir(projectFolder)

# network, external_cycles = preprocess_sbml(args)
# print("Reactions: ", network.reactions)
# print("Metabolites: ", network.metabolites)
# print("N: ", network.N)
# print(network.N[0][2])
# print(network.external_metabolite_indices())
# for i in range(len(network.metabolites)):
#    print("i:", i,"; name:", network.metabolites[i].id, "; external:", network.metabolites[i].is_external)
# runMarashiWithPolco(network.N, projectOntoDimensions, "/home/martin/projects/python-projects/secondDraft/testResults/polco/", True)
# runMarashiWithMPLRS(network.N, projectOntoDimensions, "/home/martin/projects/python-projects/secondDraft/testResults/mplrs/", mplrsPath, 4, True)
# runFELWithPolco(network.N, projectOntoDimensions, "/home/martin/projects/python-projects/secondDraft/testResults/fel/", mplrsPath, 4, True)
# print(len(network.reactions))

# inpMatrix = np.zeros(network.N.shape)
# for i in range(network.N.shape[0]):
#    for j in range(network.N.shape[1]):
#        val = int(mpfr(network.N[i][j], 53))
#        inpMatrix[i][j] = val
#print(inpMatrix.shape)
#print(network.N)
#file = open("/home/martin/projects/python-projects/secondDraft/test.matrix","w")
#np.savetxt(file, network.N)
#print(network.N.shape)
#inpDims = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46]
# inpDims = [i for i in range(inpMatrix.shape[1])]
#print(inpDims)
# proCEMs, efps = runMarashiWithPolco(inpMatrix, inpDims, os.path.join(projectFolder, "testResults", "polco"), True)
proCEMs, efps = runMarashiWithMPLRS(inpMatrix, inpDims, os.path.join(projectFolder, "testResults", "mplrs"), mplrsPath, mplrsProcesses=32)
#proCEMs, efps = runFELWithPolco(inpMatrix, inpDims, os.path.join(projectFolder, "testResults", "fel"), mplrsPath, mplrsProcesses=16)
print(f"proCEMS: {len(proCEMs)}")
print(f"efps: {len(efps)}")
