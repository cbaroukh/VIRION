# -*- coding:Utf8 -*-

#############################################################################
# Program Python type
# authors: Gerlin et al., 2025
#############################################################################

#############################################################################
# External functions
import pandas as pd
from lxml import etree


#############################################################################
# Local functions

def save_data(data, name):
    with open('{0}.txt'.format(name), 'w') as filehandle:
        for listitem in data:
            filehandle.write('%s\n' % listitem)


#############################################################################
# SBML file
filename = 'input/Sl2183.xml'
output_filename = "output/parser/parser_VYTOP/VYTOP"


# open calib file and assign values
calib = pd.read_csv("input/calibration.csv", sep=";", na_values="na", index_col=0)

print(calib)
print()

nh4_no3 = float(calib.loc['nh4_to_no3_ratio'].value)
leaf_stem_photosynthesis = float(calib.loc['photosynthesis'].value)
cost = float(calib.loc['cost'].value)
xyl = int(calib.loc['xylem_constraint'].value)
atpm_l = float(calib.loc['ATP_maintenance_leaf'].value)
atpm_r = float(calib.loc['ATP_maintenance_root'].value)
atpm_s = float(calib.loc['ATP_maintenance_stem'].value)
leaf_biomass_growth = float(calib.loc['Leaf_biomass_growth'].value)
unit = 'mM'
ratio_Rubisco_carboxylase_oxygenase = float(calib.loc['Rubisco_carboxylase_oxygenase_ratio'].value)
transpiration_limit = float(calib.loc['transpiration_limit'].value)

# organs relative weight from calib file
LEAF_WEIGHT = float(calib.loc['leaf_relative_weight'].value)
STEM_WEIGHT = float(calib.loc['stem_relative_weight'].value)
ROOT_WEIGHT = float(calib.loc['root_relative_weight'].value)


#############################################################################

#   PARSE THE SBML FILE AND GENERATE 3 ORGAN MODEL

#############################################################################


# Exchange compartments
organs = ['root', 'leaf', 'stem']
organs_id = ['r', 'l', 's']
nb_organs = len(organs_id)
exchange_compartments_names = ['xylem', 'phloem']
exchange_compartments_id = ['xyl', 'phl']
exchange_compartments_metabolites = {
    'xyl': ['M_phe__L_c', 'M_tyr__L_c', 'M_asn__L_c', 'M_asp__L_c',
            'M_gln__L_c', 'M_arg__L_c', 'M_sucr_c', 'M_glc__D_c',
            'M_thr__L_c', 'M_pro__L_c', 'M_val__L_c', 'M_so4_c',
            'M_ile__L_c', 'M_leu__L_c', 'M_lys__L_c', 'M_etoh_c',
            'M_no3_c', 'M_nh4_c', 'M_pi_c', 'M_ca2_c', 'M_mg2_c',
            'M_ala__L_c', 'M_fum_c', 'M_h2o_c', 'M_k_c', 'M_no2_c',
            'M_na1_c', 'M_cl_c', 'M_h_c', 'M_fe2_c', 'M_o2_c', 'M_co2_c',
            'M_mobd_c', 'M_mn2_c', 'M_cobalt2_c', 'M_ptrc_c'],
    'phl': ['M_glc__D_c', 'M_sucr_c', 'M_pro__L_c', 'M_mal__L_c',
            'M_4abut_c', 'M_cit_c', 'M_quin_c', 'M_o2_c', 'M_co2_c']}

exchange_compartments_direction_RL = {'xyl': ['r', 'l'],
                                      'phl': ['l', 'r']
                                      }

exchange_compartments_direction_S = {'xyl': ['s'],
                                     'phl': ['s']
                                     }

organs_weights = {'r': ROOT_WEIGHT, 'l': LEAF_WEIGHT}

PLANT_WEIGHT = LEAF_WEIGHT + STEM_WEIGHT + ROOT_WEIGHT


# create a tree by parsing the xml file
tree = etree.parse(filename)

# get the root of the tree
root = tree.getroot()

# get and parse compartments
compartments_name = []
compartments_id = []

listOfCompartments = root[0][3]
print('Parsing VYTOP metabolic network...')
print('General informations on VYTOP 1 organ model')
print('Nb Compartments: ', len(listOfCompartments))

for compartment in listOfCompartments:
    compartments_name.append(compartment.get('name'))
    compartments_id.append(compartment.get('id'))

# get genes nodes and print statistics
listOfGenes = root[0][1]
print('Nb Genes: ', len(listOfGenes))

# get metabolites nodes and print statistics
listOfSpecies = root[0][4]
print('Nb Metabolites: ', len(listOfSpecies))

# get reactions nodes and print statistics
listOfReactions = root[0][6]
print('Nb Reactions: ', len(listOfReactions))

# get the metabolites and reactions of the 1 organ network
metab_id_1comp = []
react_id_1comp = []

for metabolite in listOfSpecies:
    metab_id_1comp.append(metabolite.get('id'))

for reaction in listOfReactions:
    react_id_1comp.append(reaction.get('id'))

# generation of the 3 comp network

# generation of metabolites
metabolites_name = []
metabolites_id = []
metabolites_boundaryCondition = []
boundary_metabolites = []

for metabolite in listOfSpecies:
    for i, organ in enumerate(organs):
        metab_name = metabolite.get('name') + '_' + organ
        metab_id = metabolite.get('id') + '_' + organs_id[i]
        metabolites_name.append(metab_name)
        metabolites_id.append(metab_id)

        if metabolite.get('boundaryCondition') == 'true':
            metabolites_boundaryCondition.append(1)
            boundary_metabolites.append(metab_id)
        else:
            metabolites_boundaryCondition.append(0)


# metabolites in exchange compartment (xylem, phloeme)
for i, comp in enumerate(exchange_compartments_id):
    for metab in exchange_compartments_metabolites[comp]:
        metabolites_id.append(metab[0:-1] + comp)
        metabolites_name.append(metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                + exchange_compartments_names[i])



# generation of reactions
reactions_name = []
reactions_id = []
reactions_reversible = []
reactions_formulas = []

stoichMat_rows = []
stoichMat_columns = []
stoichMat_values = []


for ind, reaction in enumerate(listOfReactions):

    for i, organ in enumerate(organs):
        reaction_formula = ''
        reactions_name.append(reaction.get('name') + '_' + organ)
        reacid = reaction.get('id') + '_' + organs_id[i]
        reactions_id.append(reacid)


        listOfReactants = reaction[- 2]

        for substrate in listOfReactants:
            stoichMat_rows.append('QSSA_' + substrate.get('species') + '_' + organs_id[i])
            stoichMat_columns.append(reacid)
            if 'R_BIOMASS' in reacid and unit == 'mM':
                stoichMat_values.append(- float(substrate.get('stoichiometry'))/1000)
            else:
                stoichMat_values.append(- float(substrate.get('stoichiometry')))

            reaction_formula = reaction_formula + \
                               substrate.get('stoichiometry') + \
                               ' ' + substrate.get('species') + '_' + \
                               organs_id[i] + ' + '


        reaction_formula = reaction_formula[0:-2]

        if reaction.get('reversible') == 'true':

            reactions_reversible.append(1)
            reaction_formula = reaction_formula + '<--> '

        else:
            reactions_reversible.append(0)
            reaction_formula = reaction_formula + '--> '

        listOfProducts = reaction[- 1]

        for product in listOfProducts:
            stoichMat_rows.append('QSSA_' + product.get('species') + '_' + organs_id[i])

            stoichMat_columns.append(reacid)

            stoichMat_values.append(float(product.get('stoichiometry')))

            reaction_formula = reaction_formula + ' ' + \
                               product.get('stoichiometry') + ' ' + \
                               product.get('species') + '_' + \
                               organs_id[i] + ' + '

            reaction_formula = reaction_formula[0:-3]

        reactions_formulas.append(reaction_formula)

# Exchanges between compartments
nb_reac = len(reactions_id)

atp_cost_list = []
indexes_constraint_stem_weight = []
indexes_constraint_leaf_root_weight = []

for i, id in enumerate(exchange_compartments_id):

    for j, metab in enumerate(exchange_compartments_metabolites[id]):

        # organs -> exchange compartment (xylem or phloem)
        # Reaction names and reversibility
        reactions_name.append('Exchange_' + organs[organs_id.index(exchange_compartments_direction_RL[id][0])]
                              + '_' + exchange_compartments_names[i] + '_' + metab)

        reacid = 'Exch_' + exchange_compartments_direction_RL[id][0] + '_' + id + '_' + metab
        reactions_id.append(reacid)
        reactions_reversible.append(0)

        # Stoichiometric matrix
        if metab != 'M_h2o_c':
            stoichMat_rows.append('QSSA_' + 'M_atp_cost_b_' + exchange_compartments_direction_RL[id][0])
            stoichMat_columns.append(reacid)
            stoichMat_values.append(-cost)
            atp_cost_list.append(len(stoichMat_values) - 1)

        stoichMat_rows.append('QSSA_' + metab[0:-1] + id)
        stoichMat_columns.append(reacid)
        stoichMat_values.append(organs_weights[exchange_compartments_direction_RL[id][0]] / PLANT_WEIGHT)
        indexes_constraint_leaf_root_weight.append(len(stoichMat_values) - 1)

        stoichMat_rows.append('QSSA_' + metab[0:-1] + 'c_' + exchange_compartments_direction_RL[id][0])
        stoichMat_columns.append(reacid)
        stoichMat_values.append(-1)

        # Reaction Formula
        if metab != 'M_h2o_c':
            reactions_formulas.append('1 ' + metab + '_' + exchange_compartments_direction_RL[id][0] + ' + '
                                      + str(cost) + ' M_atp_cost_b_' + exchange_compartments_direction_RL[id][0]
                                      + ' --> '
                                      + str(organs_weights[exchange_compartments_direction_RL[id][0]] / PLANT_WEIGHT)
                                      + ' ' + metab[0:-1] + id)

        else:
            reactions_formulas.append('1 ' + metab + '_' + exchange_compartments_direction_RL[id][0]
                                      + ' --> '
                                      + str(organs_weights[exchange_compartments_direction_RL[id][0]] / PLANT_WEIGHT)
                                      + ' ' + metab[0:-1] + id)


        # exchange compartement (xylem or phloem) -> organs
        # Reaction names and reversibility
        reactions_name.append('Exchange_' + exchange_compartments_names[i]
                              + '_' + organs[organs_id.index(exchange_compartments_direction_RL[id][1])]
                              + '_' + metab)
        reacid = 'Exch_' + id + '_' + exchange_compartments_direction_RL[id][1] + '_' + metab
        reactions_id.append(reacid)
        reactions_reversible.append(0)

        # Stoichiometric matrix
        if metab != 'M_h2o_c':
            stoichMat_rows.append('QSSA_' + 'M_atp_cost_b_' + exchange_compartments_direction_RL[id][1])
            stoichMat_columns.append(reacid)
            stoichMat_values.append(-cost)
            atp_cost_list.append(len(stoichMat_values) - 1)

        stoichMat_rows.append('QSSA_' + metab[0:-1] + id)
        stoichMat_columns.append(reacid)
        stoichMat_values.append(-organs_weights[exchange_compartments_direction_RL[id][1]] / PLANT_WEIGHT)
        indexes_constraint_leaf_root_weight.append(len(stoichMat_values) - 1)

        stoichMat_rows.append('QSSA_'+ metab[0:-1] + 'c_' + exchange_compartments_direction_RL[id][1])
        stoichMat_columns.append(reacid)
        stoichMat_values.append(1)

        # Reaction Formula
        if metab != 'M_h2o_c':
            reactions_formulas.append(str(organs_weights[exchange_compartments_direction_RL[id][1]] /PLANT_WEIGHT)
                                      + ' ' + metab[0:-1] + id + ' + '
                                      + str(cost) + ' M_atp_cost_b_' + exchange_compartments_direction_RL[id][1] +
                                      ' --> 1 ' + metab + '_' + exchange_compartments_direction_RL[id][1])

        else:
            reactions_formulas.append(str(organs_weights[exchange_compartments_direction_RL[id][1]] / PLANT_WEIGHT)
                                      + ' ' + metab[0:-1] + id +
                                      ' --> 1 ' + metab + '_' + exchange_compartments_direction_RL[id][1])

        nb_reac = nb_reac + 2

# Additional exchanges with stem (uptake and export can both happen for both phloem and xylem)
for i, id in enumerate(exchange_compartments_id):
    for j, metab in enumerate(exchange_compartments_metabolites[id]):


        # stem -> exchange compartement (xylem or phloem)
        # Reaction names and reversibility
        reactions_name.append('Exchange_' + organs[organs_id.index(exchange_compartments_direction_S[id][0])]
                              + '_' + exchange_compartments_names[i] + '_' + metab)
        reacid = 'Exch_' + exchange_compartments_direction_S[id][0] + '_' + id + '_' + metab
        reactions_id.append(reacid)
        reactions_reversible.append(0)

        # Stoichiometric matrix
        if metab != 'M_h2o_c':
            stoichMat_rows.append('QSSA_' + 'M_atp_cost_b_' + exchange_compartments_direction_S[id][0])
            stoichMat_columns.append(reacid)
            stoichMat_values.append(-cost)
            atp_cost_list.append(len(stoichMat_values) - 1)

        stoichMat_rows.append('QSSA_' + metab[0:-1] + id)
        stoichMat_columns.append(reacid)
        stoichMat_values.append(STEM_WEIGHT/PLANT_WEIGHT)
        indexes_constraint_stem_weight.append(len(stoichMat_values) - 1)

        stoichMat_rows.append('QSSA_' + metab + '_' + exchange_compartments_direction_S[id][0])
        stoichMat_columns.append(reacid)
        stoichMat_values.append(-1.0)

        # Reaction Formula
        if metab != 'M_h2o_c':
            reactions_formulas.append('1 ' + metab + '_' +
                exchange_compartments_direction_S[id][0] + ' + ' +
                str(cost) + ' M_atp_cost_b_' + exchange_compartments_direction_S[id][0] +
                ' --> ' + str(STEM_WEIGHT / PLANT_WEIGHT) + ' ' + metab[0:-1] + id)

        else:
            reactions_formulas.append('1 ' + metab + '_' + exchange_compartments_direction_S[id][0] +
                                      ' --> ' + str(STEM_WEIGHT / PLANT_WEIGHT) + ' ' + metab[0:-1] + id)

        # exchange compartement (xylem or phloem) -> stem
        # Reaction names and reversibility
        reactions_name.append('Exchange_' + exchange_compartments_names[i] + '_'
                              + organs[organs_id.index(exchange_compartments_direction_S[id][0])] + '_' + metab)
        reacid = 'Exch_' + id + '_' + exchange_compartments_direction_S[id][0] + '_' + metab
        reactions_id.append(reacid)
        reactions_reversible.append(0)

        # Stoichiometric matrix
        if metab != 'M_h2o_c':
            stoichMat_rows.append('QSSA_' + 'M_atp_cost_b_' + exchange_compartments_direction_S[id][0])
            stoichMat_columns.append(reacid)
            stoichMat_values.append(-cost)
            atp_cost_list.append(len(stoichMat_values) - 1)

        stoichMat_rows.append('QSSA_' + metab[0:-1] + id)
        stoichMat_columns.append(reacid)
        stoichMat_values.append(-STEM_WEIGHT/PLANT_WEIGHT)
        indexes_constraint_stem_weight.append(len(stoichMat_values) - 1)


        stoichMat_rows.append('QSSA_' + metab + '_' + exchange_compartments_direction_S[id][0])
        stoichMat_columns.append(reacid)
        stoichMat_values.append(1.0)

        nb_reac = nb_reac + 2

        # Reaction Formula
        if metab != 'M_h2o_c':
            reactions_formulas.append(str(STEM_WEIGHT / PLANT_WEIGHT) + ' ' + metab[0:-1] + id + ' + ' +
                                      str(cost) + ' M_atp_cost_b_' + exchange_compartments_direction_S[id][0] +
                                      ' --> 1 ' +
                                      metab + '_' + exchange_compartments_direction_S[id][0])

        else:
            reactions_formulas.append(str(STEM_WEIGHT / PLANT_WEIGHT) + ' ' + metab[0:-1] + id + ' --> 1 ' +
                                      metab + '_' + exchange_compartments_direction_S[id][0])


# reactions for boundary metabolites
exchangereactions_id = []
exchangereactions_name = []
exchangereactions_formulas = []
exchangestoichMat_rows = []
exchangestoichMat_columns = []
exchangestoichMat_values = []


for index, mid in enumerate(boundary_metabolites):
    exchangestoichMat_rows.append('QSSA_' + mid)
    exchangestoichMat_columns.append('Exch_' + mid)
    exchangestoichMat_values.append(-1)
    exchangereactions_id.append('Exch_' + mid)
    exchangereactions_name.append('Exchange of ' + mid)
    exchangereactions_formulas.append(mid + ' <--> ')


#get some statistics
nb_reac = len(reactions_id)
nb_metab = len(metabolites_id)
nb_exchangereac = len(exchangereactions_id)

exchangereactions_reversible = [1]*nb_exchangereac

print('Number of reactions in the VYTOP three organ model: ', nb_reac)
print('Number of metabolites in the VYTOP three organ model: ', nb_metab)
print()

# saved the parsed network file
save_data(metabolites_id, output_filename + "_ID_metabolites")
save_data(boundary_metabolites, output_filename + "_boundary_metabolites")
save_data(metabolites_boundaryCondition, output_filename + "_boundaryCondition_metabolites")

save_data(reactions_id, output_filename + "_ID_reactions")
save_data(reactions_name, output_filename + "_Name_reactions")
save_data(reactions_reversible, output_filename + "_Reversibility_reactions")
save_data(reactions_formulas, output_filename + "_Formulas_reactions")

save_data(stoichMat_rows, output_filename + "_Rows_reactions")
save_data(stoichMat_columns, output_filename + "_Columns_reactions")
save_data(stoichMat_values, output_filename + "_Coefficients_reactions")

save_data(exchangereactions_id, output_filename + "_ID_exchanged_reactions")
save_data(exchangereactions_name, output_filename + "_Name_exchanged_reactions")
save_data(exchangereactions_reversible, output_filename + "_Reversibility_exchanged_reactions")
save_data(exchangereactions_formulas, output_filename + "_Formulas_exchanged_reactions")

save_data(exchangestoichMat_rows, output_filename + "_Rows_exchanged_reactions")
save_data(exchangestoichMat_columns, output_filename + "_Columns_exchanged_reactions")
save_data(exchangestoichMat_values, output_filename + "_Coefficients_exchanged_reactions")

save_data(react_id_1comp, "output/parser/parser_VYTOP/1comp_ID_reactions")
save_data(metab_id_1comp, "output/parser/parser_VYTOP/1comp_ID_metabolites")
