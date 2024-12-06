# -*- coding:Utf8 -*-

#############################################################################
# Program Python type
# authors: Gerlin et al., 2025
#############################################################################

#############################################################################
# External functions
import cplex as cplex
from cplex.exceptions import CplexError
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import FBA_VYTOP as fbaVYTOP
import parser_merge_VIRION as mnetVYTOPRalsto  # open the parsed network
from pickle import *

#############################################################################

#                                MAIN PROGRAM

#############################################################################

option_print = False
option_save_FBA = True

''' choose the simulation you want to perform for the sequential FBAs
1-PhotonsOnly
has no additional constraints: only photon uptake is taken from the fluxes of a healthy plant
2-PhotonsNitrogen
constrains photons + N uptake using values from the fluxes of a healthy plant
3-PhotonsIron
constrains photons + iron uptake using values from the fluxes of a healthy plant
4-PhotonsTranspiration
constrains photons + reduces xylem fluxes according to transpiration reduction
5-All
constrains photons + reduces xylem fluxes according to transpiration reduction
+ iron/nitrogen uptake using values from the fluxes of a healthy plant
6-All-Stem
authorizes nutrient exchanges from stem to xylem with the constrainst defined in 5
7-PhotonsOnly-Stem
authorizes nutrient exchanges from stem to xylem with the constraints defined in 1
8-PhotonsXylem
limits xylem elements from the amount in healthy plant, with the constraints defined for 1
'''

simulation_name = '1-PhotonsOnly'
# simulation_name = '2-PhotonsNitrogen'
# simulation_name = '3-PhotonsIron'
# simulation_name = '4-PhotonsTranspiration'
# simulation_name = '5-All'
# simulation_name = '6-All-Stem'
# simulation_name = '7-PhotonsOnly-Stem'
# simulation_name = '8-PhotonsXylem'

print(simulation_name)
cfu_start = 4
cfu_end = 12
nb_points = 100

# compute different cfus list
step = (cfu_end - cfu_start) / nb_points
cfus = np.arange(cfu_start, cfu_end, step)

# define substrates parameters for Ralstonia
NGAM = 0.0  # 0 if GAM modified, 9 otherwise

substrates = ["R_GLNtex_rs", "R_GLCtex_rs"]

params_growth_WT = np.array([(1.10898478e-02, 4.21294241e-04, 2.63210539e-01, 8.49750564e-01),
                             (4.13159097e-03, 2.42012991e-04, 1.14597371e-01, 1.20941856e-04)])
concentrations_substrates = [3.29, 0.010]

other_substrates = ['R_PHEtex_rs', 'R_VALtex_rs', 'R_LEUtex_rs', 'R_ARGtex_rs', 'R_SUCRtex_rs', 'R_TYRtex_rs',
                    'R_PROtex_rs', 'R_LYStex_rs', 'R_THRtex_rs', 'R_ILEtex_rs', 'R_ASNtex_rs', 'R_ASPtex_rs',
                    'R_ALAtex_rs', 'R_FUMtex_rs']
conso_coeff_other_substrate = [2.85832127e-04, 2.05239626e-04, 3.11072406e-04, 1.47043323e-04,
                               7.42787597e-04, 1.39022128e-04, 9.60740956e-04, 1.86695645e-04,
                               1.12618646e-04, 1.58110980e-04, 7.01713239e-04, 3.68270000e-04, 0, 0]

interesting_reactions_plant = ['R_EX_nh4_c_r', 'R_EX_no3_c_r', 'R_EX_co2_e_l', 'R_EX_co2_e_s', 'R_EX_co2_e_r',
                               'R_EX_photon_h_s', 'R_EX_photon_h_l']

interesting_reactions_RS = substrates + other_substrates + ['R_PTRCtex_rs']

interesting_reactions_xylem = ['Exch_r_xyl_M_ala__L_c', 'Exch_r_xyl_M_arg__L_c', 'Exch_r_xyl_M_asn__L_c',
                               'Exch_r_xyl_M_asp__L_c', 'Exch_r_xyl_M_etoh_c', 'Exch_r_xyl_M_fum_c',
                               'Exch_r_xyl_M_glc__D_c', 'Exch_r_xyl_M_ile__L_c',
                               'Exch_r_xyl_M_leu__L_c', 'Exch_r_xyl_M_lys__L_c', 'Exch_r_xyl_M_phe__L_c',
                               'Exch_r_xyl_M_pro__L_c', 'Exch_r_xyl_M_sucr_c', 'Exch_r_xyl_M_thr__L_c',
                               'Exch_r_xyl_M_tyr__L_c', 'Exch_r_xyl_M_val__L_c']

interesting_reactions_xylem_n = ['Exch_r_xyl_M_gln__L_c', 'Exch_r_xyl_M_nh4_c', 'Exch_r_xyl_M_no3_c',
                                 'Exch_r_xyl_M_no2_c']

interesting_reactions_stem = ['Exch_s_xyl_M_ala__L_c', 'Exch_s_xyl_M_arg__L_c', 'Exch_s_xyl_M_asn__L_c',
                              'Exch_s_xyl_M_asp__L_c', 'Exch_s_xyl_M_etoh_c', 'Exch_s_xyl_M_fum_c',
                              'Exch_s_xyl_M_ile__L_c', 'Exch_s_xyl_M_gln__L_c', 'Exch_s_xyl_M_glc__D_c',
                              'Exch_s_xyl_M_leu__L_c', 'Exch_s_xyl_M_lys__L_c', 'Exch_s_xyl_M_phe__L_c',
                              'Exch_s_xyl_M_pro__L_c', 'Exch_s_xyl_M_sucr_c', 'Exch_s_xyl_M_thr__L_c',
                              'Exch_s_xyl_M_tyr__L_c', 'Exch_s_xyl_M_val__L_c', 'Exch_s_xyl_M_nh4_c',
                              'Exch_s_xyl_M_no3_c', 'Exch_s_xyl_M_no2_c']

Putr_Metabolism = ['R_ORNDC_rs',
                   'R_ACODA_rs',
                   'R_ARGN_rs',
                   'R_ACOTA_rs',
                   'R_AGPR_rs',
                   'R_ACGK_rs',
                   'R_ACGS_rs',
                   'R_PTRCt2pp_rs',
                   'R_PTRCabcpp_rs',
                   'Exch_xyl_s_M_ptrc_c',
                   'Exch_xyl_l_M_ptrc_c',
                   'R_RE1537C_l',
                   'R_PMT_l',
                   'R_NMPTRCOX2_l',
                   'R_1MPYRS_l',
                   'R_PTRCOX1_l',
                   'R_ABUTD_l',
                   'R_ABTArm_l',
                   'R_4abut_tx_m__l',
                   'R_CATp_l',
                   'R_CATp_s',
                   'R_URIC_l',
                   'R_URIC_s',
                   'R_RE1537C_s',
                   'R_PMT_s',
                   'R_NMPTRCOX2_s',
                   'R_1MPYRS_s',
                   'R_PTRCOX1_s',
                   'R_ABUTD_s',
                   'R_ABTArm_s',
                   'R_4abut_tx_m__s',
                   'R_ORNabcpp_rs',
                   'R_OCBT_rs',
                   'R_CITRabcpp_rs',
                   'R_ARGSS_rs',
                   'R_ARGSL_rs',
                   'R_ARGDC_rs',
                   'R_AGMT_rs']

interesting_reactions = interesting_reactions_plant + interesting_reactions_RS + interesting_reactions_xylem \
                        + interesting_reactions_xylem_n + interesting_reactions_stem + Putr_Metabolism
val_atp = str(mnetVYTOPRalsto.mVYTOP.cost)

#############################################################################

#                DEFINE THE CONSTRAINTS ON THE 3 COMP MODEL

#############################################################################
# Upper bounds
ub = [cplex.infinity] * (len(mnetVYTOPRalsto.reactions_id) + len(mnetVYTOPRalsto.exchangereactions_id))

# Lower bounds. Set to 0 for exchange reactions to avoid infinite uptake of substrates other than those wished
lb = [-cplex.infinity if rev == 1 else 0 for rev in mnetVYTOPRalsto.reactions_reversible] \
     + [0] * mnetVYTOPRalsto.nb_exchangereac * 1

# Additional constraints such as ATP maintenance ou leaf biomass growth

# ATP maintenance
lb[mnetVYTOPRalsto.reactions_id.index('R_ATPS_l')] = mnetVYTOPRalsto.mVYTOP.atpm_l
lb[mnetVYTOPRalsto.reactions_id.index('R_ATPS_s')] = mnetVYTOPRalsto.mVYTOP.atpm_s
lb[mnetVYTOPRalsto.reactions_id.index('R_ATPS_r')] = mnetVYTOPRalsto.mVYTOP.atpm_r

# Change of bounds for boundary metabolites allowed to be exchanged
infinite_exchange = ['M_ca2_b_r',
                     'M_cl_b_r',
                     'M_so4_b_r',
                     'M_pi_b_r',
                     'M_fe2_b_r',
                     'M_k_b_r',
                     'M_mg2_b_r',
                     'M_na1_b_r',
                     'M_nh4_b_r',
                     'M_no3_b_r',
                     'M_h2o_b_r',
                     'M_co2_b_r', 'M_co2_b_s', 'M_co2_b_l',
                     'M_o2_b_r', 'M_o2_b_s', 'M_o2_b_l',
                     'M_h_b_r', 'M_h_b_s', 'M_h_b_l',
                     'M_photon_b_s', 'M_photon_b_l',
                     'M_nadh_work_b_r', 'M_nadh_work_b_s', 'M_nadh_work_b_l',
                     'M_nadph_work_b_r', 'M_nadph_work_b_s', 'M_nadph_work_b_l',
                     'M_atp_work_b_r', 'M_atp_work_b_s', 'M_atp_work_b_l',
                     'M_biomass_leaf_b_l',
                     'M_biomass_root_b_r',
                     'M_biomass_stem_b_s',
                     'M_mobd_b_r', 'M_mn2_b_r', 'M_cobalt2_b_r']

for id in infinite_exchange:
    lb[mnetVYTOPRalsto.nb_reac + mnetVYTOPRalsto.exchangereactions_id.index('Exch_' + id)] = -cplex.infinity

fordidden_exchange = ['M_nh4_b_l', 'M_nh4_b_s']

for id in fordidden_exchange:
    lb[mnetVYTOPRalsto.nb_reac + mnetVYTOPRalsto.exchangereactions_id.index('Exch_' + id)] = 0
    ub[mnetVYTOPRalsto.nb_reac + mnetVYTOPRalsto.exchangereactions_id.index('Exch_' + id)] = 0

# [OPTIONAL] Constraint on xylem composition (taken from xylem.csv file)
if mnetVYTOPRalsto.mVYTOP.xyl == 1 or simulation_name=='8-PhotonsXylem':
    xylem = pd.read_csv("input/xylem.csv", sep=";")
    for num, id in enumerate(xylem.id):
        ub[mnetVYTOPRalsto.reactions_id.index(id)] = xylem.ub[num]

# leaf
ub[mnetVYTOPRalsto.reactions_id.index('R_EX_photon_h_l')] = fbaVYTOP.x[
    mnetVYTOPRalsto.reactions_id.index('R_EX_photon_h_l')]

# stem
ub[mnetVYTOPRalsto.reactions_id.index('R_EX_photon_h_s')] = fbaVYTOP.x[
    mnetVYTOPRalsto.reactions_id.index('R_EX_photon_h_s')]

# Variables' initialisation for additional constraints of the LP problem, not on lb or ub
addnames = []
addstoichMat_rows = []
addstoichMat_values = []
addstoichMat_columns = []
addrhs = []
addsense = ""

# Constraint for relative growth rates of each organs (based on experimental data)
addnames.append("Stem to Leaf relative growth rate")
addstoichMat_values.append(0.26 / 0.21)
addstoichMat_values.append(-1)

addstoichMat_rows.append("Stem to Leaf relative growth rate")
addstoichMat_rows.append("Stem to Leaf relative growth rate")

addstoichMat_columns.append('R_BIOMASS_LEAF_l')
addstoichMat_columns.append('R_BIOMASS_STEM_s')

addrhs.append(0)

addsense = addsense + "E"

# leaf stem photosynthesis constraint

addnames.append('leaf stem photosynthesis constraint')

addstoichMat_values.append(mnetVYTOPRalsto.mVYTOP.leaf_stem_photosynthesis)
addstoichMat_values.append(-1)

addstoichMat_rows.append('leaf stem photosynthesis constraint')
addstoichMat_rows.append('leaf stem photosynthesis constraint')

addstoichMat_columns.append('R_EX_photon_h_s')
addstoichMat_columns.append('R_EX_photon_h_l')

addrhs.append(0)

addsense = addsense + "L"

# [OPTIONAL] nitrogen ratio constraint

if not np.isnan(mnetVYTOPRalsto.mVYTOP.nh4_no3):
    addnames.append("NH4 to NO3 uptake ratio")

    addstoichMat_values.append(mnetVYTOPRalsto.mVYTOP.nh4_no3)
    addstoichMat_values.append(-1)

    addstoichMat_rows.append("NH4 to NO3 uptake ratio")
    addstoichMat_rows.append("NH4 to NO3 uptake ratio")

    addstoichMat_columns.append('R_EX_no3_c_r')
    addstoichMat_columns.append('R_EX_nh4_c_r')

    addrhs.append(0)

    addsense = addsense + "E"

# [OPTIONAL] carboxylase oxygenase Rubisco activity ratio

if not np.isnan(mnetVYTOPRalsto.mVYTOP.ratio_Rubisco_carboxylase_oxygenase):
    addnames.append("Rubisco CO2 vs O2 ratio activity leaf")

    addstoichMat_values.append(1)
    addstoichMat_values.append(-mnetVYTOPRalsto.mVYTOP.ratio_Rubisco_carboxylase_oxygenase)

    addstoichMat_rows.append("Rubisco CO2 vs O2 ratio activity leaf")
    addstoichMat_rows.append("Rubisco CO2 vs O2 ratio activity leaf")

    addstoichMat_columns.append('R_RBPCh_l')
    addstoichMat_columns.append('R_RBCh_1_l')

    addrhs.append(0)

    addsense = addsense + "E"

    addnames.append("Rubisco CO2 vs O2 ratio activity stem")

    addstoichMat_values.append(1)
    addstoichMat_values.append(-mnetVYTOPRalsto.mVYTOP.ratio_Rubisco_carboxylase_oxygenase)

    addstoichMat_rows.append("Rubisco CO2 vs O2 ratio activity stem")
    addstoichMat_rows.append("Rubisco CO2 vs O2 ratio activity stem")

    addstoichMat_columns.append('R_RBPCh_s')
    addstoichMat_columns.append('R_RBCh_1_s')

    addrhs.append(0)

    addsense = addsense + "E"

# Additional constraints on Ralstonia

# set bound on import fluxes of carbon substrates
for i, substrate in enumerate(substrates):
    ub[mnetVYTOPRalsto.reactions_id.index(substrate)] = params_growth_WT[i, 0] * params_growth_WT[i, 2] \
                                                        * concentrations_substrates[i] / (
                                                                    concentrations_substrates[i] + params_growth_WT[
                                                                i, 3]) \
                                                        * 1000 * 24

# set constraints for other substrates:
rhs_os = [0] * len(other_substrates)
sense_os = 'L' * len(other_substrates)
constraint_names_os = []
rows_os = []
columns_os = []
values_os = []

for i, other_substrate in enumerate(other_substrates):
    ub[mnetVYTOPRalsto.reactions_id.index(other_substrate)] = cplex.infinity
    constraint_names_os = constraint_names_os + ['SynchronizedAssimilationOf_' + other_substrate]
    rows_os = rows_os + ['SynchronizedAssimilationOf_' + other_substrate] * 3
    columns_os = columns_os + [other_substrate, substrates[0], substrates[1]]
    values_os = values_os + [1, - conso_coeff_other_substrate[i] / params_growth_WT[0, 0],
                             - conso_coeff_other_substrate[i] / params_growth_WT[1, 0]]

# set additional constraints for putrescine excretion
rhs_ptrc = [0]
sense_ptrc = 'E'
constraint_names_ptrc = ['Putrescine excretion']
rows_ptrc = ['Putrescine excretion'] * (len(substrates) + 1)
columns_ptrc = ['R_PTRCtex_rs']
values_ptrc = [1]

for i, substrate in enumerate(substrates):
    columns_ptrc = columns_ptrc + [substrate]
    values_ptrc = values_ptrc + [params_growth_WT[i, 1] / params_growth_WT[i, 0]]

# set constraints on some specific reactions
lb[mnetVYTOPRalsto.reactions_id.index('R_FEOXpp_rs')] = 0
ub[mnetVYTOPRalsto.reactions_id.index('R_FEOXpp_rs')] = 0

lb[mnetVYTOPRalsto.reactions_id.index('R_NGAME_rs')] = NGAM * 24

# Constraint on transpiration, set on root to xylem flow
no_transpiration_effect = ['Exch_r_xyl_M_h2o_c', 'Exch_r_xyl_M_h_c', 'Exch_r_xyl_M_o2_c', 'Exch_r_xyl_M_co2_c',
                           'Exch_r_xyl_M_cobalt2_c', 'Exch_r_xyl_M_mn2_c', 'Exch_r_xyl_M_mobd_c']

if simulation_name in ('4-PhotonsTranspiration', '5-All', '6-All-Stem'):
    print("Transpiration limit on")
    reduction_flow_coeff = -0.2746 * mnetVYTOPRalsto.patho_density + 2.8764
    if reduction_flow_coeff < 0.0:
        reduction_flow_coeff = 0.0

    if mnetVYTOPRalsto.patho_density > 7:
        for num, reac in enumerate(mnetVYTOPRalsto.reactions_id):
            if 'Exch_r_xyl_' in reac and reac not in no_transpiration_effect:
                ub[num] = fbaVYTOP.x[mnetVYTOPRalsto.reactions_id.index(reac)] * reduction_flow_coeff
else:
    print("Transpiration limit off")

# set root growth which continues according to experimental data (facultative)

lb[mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_ROOT_r')] = fbaVYTOP.x[
    mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_ROOT_r')]

# forbid stem to contribute to xylem
if not (simulation_name in ('6-All-Stem', '7-PhotonsOnly-Stem')):
    for num, reac in enumerate(mnetVYTOPRalsto.reactions_id):
        if 'Exch_s_xyl_' in reac:
            ub[num] = 0
            lb[num] = 0

# remove the possiblity for ptrc to go to sink

lb[mnetVYTOPRalsto.nb_reac + mnetVYTOPRalsto.exchangereactions_id.index('Exch_sink_M_ptrc_xyl')] = 0
ub[mnetVYTOPRalsto.nb_reac + mnetVYTOPRalsto.exchangereactions_id.index('Exch_sink_M_ptrc_xyl')] = 0

# Right-hand side of the LP problem
rhs = [0] * mnetVYTOPRalsto.nb_metab + addrhs + rhs_os + rhs_ptrc

# Senses of the LP problem
sense = "E" * mnetVYTOPRalsto.nb_metab + addsense + sense_os + sense_ptrc

# Objective function of the LP problem (minimizations of photons import, stem and leaf)
obj = [0.0] * (len(mnetVYTOPRalsto.reactions_id) + len(mnetVYTOPRalsto.exchangereactions_id))
obj[mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_rs')] = 1

# Constraint Names
constraint_names = ["QSSA_" + id for id in
                    mnetVYTOPRalsto.metabolites_id] + addnames + constraint_names_os + constraint_names_ptrc

# variable names
variable_names = mnetVYTOPRalsto.reactions_id + mnetVYTOPRalsto.exchangereactions_id

old_patho_weight = mnetVYTOPRalsto.patho_weight

# initiate result variables
objectives1 = ['na'] * len(cfus)
objectives2 = ['na'] * len(cfus)
objectives3 = ['na'] * len(cfus)

status1 = ['na'] * len(cfus)
status2 = ['na'] * len(cfus)
status3 = ['na'] * len(cfus)

interesting_reactions_fluxes = {reac_id: ['na'] * len(cfus) for reac_id in interesting_reactions}

if simulation_name in ('3-PhotonsIron', '5-All', '6-All-Stem'):
    print('IRON CONSTRAINED')
    ub[mnetVYTOPRalsto.reactions_id.index('R_EX_fe2_c_r')] = fbaVYTOP.x[
        mnetVYTOPRalsto.reactions_id.index('R_EX_fe2_c_r')]

if simulation_name in ('2-PhotonsNitrogen', '5-All', '6-All-Stem'):
    print('NITROGEN CONSTRAINED')
    ub[mnetVYTOPRalsto.reactions_id.index('R_EX_nh4_c_r')] = fbaVYTOP.x[
        mnetVYTOPRalsto.reactions_id.index('R_EX_nh4_c_r')]
    ub[mnetVYTOPRalsto.reactions_id.index('R_EX_no3_c_r')] = fbaVYTOP.x[
        mnetVYTOPRalsto.reactions_id.index('R_EX_no3_c_r')]

# solve VIRION for different cfus
for i, patho_density in enumerate(cfus):
    patho_weight = 1.11E-11 * 10 ** patho_density
    print("Ralstonia density: ", patho_density)

    # change stoichmat_values involving pathogen weight
    for index in mnetVYTOPRalsto.stoichMat_indices_for_patho_weight:
        mnetVYTOPRalsto.stoichMat_values[index] = mnetVYTOPRalsto.stoichMat_values[
                                                      index] * patho_weight / old_patho_weight

    if simulation_name in ('4-PhotonsTranspiration', '5-All', '6-All-Stem'):

        # change transpiration according to pathogen weight
        reduction_flow_coeff = -0.2746 * patho_density + 2.8764
        if reduction_flow_coeff < 0.0:
            reduction_flow_coeff = 0.0

        if patho_density > 7:
            for num, reac in enumerate(mnetVYTOPRalsto.reactions_id):
                if 'Exch_r_xyl_' in reac and reac not in no_transpiration_effect:
                    ub[num] = fbaVYTOP.x[mnetVYTOPRalsto.reactions_id.index(reac)] * reduction_flow_coeff
                    print(ub[num])
                    print(fbaVYTOP.x[mnetVYTOPRalsto.reactions_id.index(reac)])

    old_patho_weight = patho_weight

    # constraints
    constraints = zip(
        mnetVYTOPRalsto.stoichMat_rows +
        mnetVYTOPRalsto.exchangestoichMat_rows + addstoichMat_rows + rows_os + rows_ptrc,
        mnetVYTOPRalsto.stoichMat_columns + mnetVYTOPRalsto.exchangestoichMat_columns +
        addstoichMat_columns + columns_os + columns_ptrc,
        mnetVYTOPRalsto.stoichMat_values + mnetVYTOPRalsto.exchangestoichMat_values +
        addstoichMat_values + values_os + values_ptrc)

    #############################################################################

    #                      FLUX BALANCE ANALYSIS
    #                     1) Max Biom Ralsto
    #                     2) Max Biom Tomato
    #                     3) Min Sum(abs(flux))

    #############################################################################

    try:
        #####################################
        # FBA 1) Max Biom Ralsto
        #####################################

        print()
        print("Solving  Max Biom Ralsto FBA for VIRION model")

        probFBA = cplex.Cplex()

        probFBA.set_problem_name("FBATomatoMaxBiomRalsto")

        probFBA.objective.set_sense(probFBA.objective.sense.maximize)

        probFBA.linear_constraints.add(rhs=rhs, senses=sense,
                                       names=constraint_names)

        probFBA.variables.add(obj=obj, lb=lb, ub=ub, names=variable_names)

        probFBA.linear_constraints.set_coefficients(constraints)

        probFBA.write("output/LP/" + "FBA2-MaxBiomRalsto_" + simulation_name + "_probFBA.lp")

        probFBA.solve()

        x = probFBA.solution.get_values()
        growthRalsto = probFBA.solution.get_objective_value()
        status1[i] = probFBA.solution.get_status()
        if probFBA.solution.get_status() == 1:
            objectives1[i] = growthRalsto

        # Print results
        print("Solution status = ", probFBA.solution.get_status(), ":", )
        print(probFBA.solution.status[probFBA.solution.get_status()])
        print("Objective value  = ", probFBA.solution.get_objective_value())
        print("Iteration count = ", probFBA.solution.progress.get_num_iterations())

        if option_print:

            print("NH4", x[mnetVYTOPRalsto.reactions_id.index('R_EX_nh4_c_r')])
            print("NO3", x[mnetVYTOPRalsto.reactions_id.index('R_EX_no3_c_r')])
            print("co2_l", x[mnetVYTOPRalsto.reactions_id.index('R_EX_co2_e_l')])
            print("co2_s", x[mnetVYTOPRalsto.reactions_id.index('R_EX_co2_e_s')])
            print("co2_r", x[mnetVYTOPRalsto.reactions_id.index('R_EX_co2_e_r')])
            print('PS stem', x[mnetVYTOPRalsto.reactions_id.index('R_EX_photon_h_s')])
            print('PS leaf', x[mnetVYTOPRalsto.reactions_id.index('R_EX_photon_h_l')])

            print('Leaf growth', x[mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_LEAF_l')])
            print('Stem growth', x[mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_STEM_s')])
            print('Root growth', x[mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_ROOT_r')])

            for substrate in substrates:
                print('RS import ' + substrate, x[mnetVYTOPRalsto.reactions_id.index(substrate)])

            for other_substrate in other_substrates:
                print('RS import ' + other_substrate, x[mnetVYTOPRalsto.reactions_id.index(other_substrate)])

            print('RS putrescine', x[mnetVYTOPRalsto.reactions_id.index('R_PTRCtex_rs')])

        #####################################
        # FBA 2) Max Biom Tomato
        #####################################

        print()
        print("Solving  Max Biom Tomato FBA for VIRION model")

        probFBA.set_problem_name("FBATomatoMaxBiomTomato")

        probFBA.objective.set_linear('R_BIOMASS_LEAF_l', 1.0)
        probFBA.objective.set_linear('R_BIOMASS_rs', 0.0)

        probFBA.variables.set_lower_bounds('R_BIOMASS_rs', growthRalsto)

        probFBA.variables.set_lower_bounds('R_ATPS_l', mnetVYTOPRalsto.mVYTOP.atpm_l)
        probFBA.variables.set_lower_bounds('R_ATPS_s', mnetVYTOPRalsto.mVYTOP.atpm_s)

        probFBA.write("output/LP/" + "FBA3-MaxBiomTomato_" + simulation_name + "_probFBA.lp")

        probFBA.solve()

        x = probFBA.solution.get_values()
        growthLeaf = probFBA.solution.get_objective_value()
        status2[i] = probFBA.solution.get_status()
        if probFBA.solution.get_status() == 1:
            objectives2[i] = growthLeaf

        # Print results
        print("Solution status = ", probFBA.solution.get_status(), ":", )
        print(probFBA.solution.status[probFBA.solution.get_status()])
        print("Objective value  = ", probFBA.solution.get_objective_value())
        print("Iteration count = ", probFBA.solution.progress.get_num_iterations())

        if option_print:

            print("NH4", x[mnetVYTOPRalsto.reactions_id.index('R_EX_nh4_c_r')])
            print("NO3", x[mnetVYTOPRalsto.reactions_id.index('R_EX_no3_c_r')])
            print("co2_l", x[mnetVYTOPRalsto.reactions_id.index('R_EX_co2_e_l')])
            print("co2_s", x[mnetVYTOPRalsto.reactions_id.index('R_EX_co2_e_s')])
            print("co2_r", x[mnetVYTOPRalsto.reactions_id.index('R_EX_co2_e_r')])
            print('PS stem', x[mnetVYTOPRalsto.reactions_id.index('R_EX_photon_h_s')])
            print('PS leaf', x[mnetVYTOPRalsto.reactions_id.index('R_EX_photon_h_l')])

            print('Leaf growth', x[mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_LEAF_l')])
            print('Stem growth', x[mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_STEM_s')])
            print('Root growth', x[mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_ROOT_r')])

            for substrate in substrates:
                print('RS import ' + substrate, x[mnetVYTOPRalsto.reactions_id.index(substrate)])

            for other_substrate in other_substrates:
                print('RS import ' + other_substrate, x[mnetVYTOPRalsto.reactions_id.index(other_substrate)])

            print('RS putrescine', x[mnetVYTOPRalsto.reactions_id.index('R_PTRCtex_rs')])

        #####################################
        # FBA 3) Min Sum(abs(flux))
        #####################################

        print()
        print("Solving  Min Sum Flux FBA for VIRION model")

        # define new name of optimization problem
        probFBA.set_problem_name("FBATomatoRSMinFluxSum")

        # add constraints for biomass growth from previous FBA
        probFBA.variables.set_lower_bounds('R_BIOMASS_LEAF_l', growthLeaf)

        # define additional variables names
        variable_names2 = ['MinVar_' + id for id in mnetVYTOPRalsto.reactions_id]

        # define new obj
        probFBA.objective.set_sense(probFBA.objective.sense.minimize)
        probFBA.objective.set_linear('R_BIOMASS_LEAF_l', 0.0)
        obj2 = [1] * len(mnetVYTOPRalsto.reactions_id)

        # new bounds for additional variables
        lb2 = [-cplex.infinity] * len(mnetVYTOPRalsto.reactions_id)

        ub2 = [cplex.infinity] * len(mnetVYTOPRalsto.reactions_id)

        # add new variables
        probFBA.variables.add(obj=obj2, lb=lb2, ub=ub2, names=variable_names2)

        # define additional rhs
        rhs2 = [0] * len(mnetVYTOPRalsto.reactions_id) * 2

        # define additional sens
        sense2 = "L" * len(mnetVYTOPRalsto.reactions_id) * 2

        # define additional constraint names
        constraint_names2 = ["MinFluxConstraint1_" + id for id in mnetVYTOPRalsto.reactions_id] + \
                            ["MinFluxConstraint2_" + id for id in mnetVYTOPRalsto.reactions_id]

        # define additional constraints for the LP problem
        minstoichMat_rows = ["MinFluxConstraint1_" + id for id in mnetVYTOPRalsto.reactions_id] \
                            + ["MinFluxConstraint1_" + id for id in mnetVYTOPRalsto.reactions_id] \
                            + ["MinFluxConstraint2_" + id for id in mnetVYTOPRalsto.reactions_id] \
                            + ["MinFluxConstraint2_" + id for id in mnetVYTOPRalsto.reactions_id]

        minstoichMat_columns = mnetVYTOPRalsto.reactions_id \
                               + ['MinVar_' + id for id in mnetVYTOPRalsto.reactions_id] \
                               + mnetVYTOPRalsto.reactions_id \
                               + ['MinVar_' + id for id in mnetVYTOPRalsto.reactions_id]

        minstoichMat_values = [1] * len(mnetVYTOPRalsto.reactions_id) + [-1] * len(mnetVYTOPRalsto.reactions_id) * 3

        constraints2 = zip(minstoichMat_rows, minstoichMat_columns, minstoichMat_values)

        probFBA.linear_constraints.add(rhs=rhs2, senses=sense2,
                                       names=constraint_names2)

        probFBA.linear_constraints.set_coefficients(constraints2)

        probFBA.write("output/LP/" + "FBA4-MinFluxSum_" + simulation_name + "_probFBA.lp")

        probFBA.solve()

        x = probFBA.solution.get_values()
        status3[i] = probFBA.solution.get_status()
        if probFBA.solution.get_status() == 1:
            objectives3[i] = probFBA.solution.get_objective_value()
            for reac_id in interesting_reactions:
                interesting_reactions_fluxes[reac_id][i] = x[mnetVYTOPRalsto.reactions_id.index(reac_id)]

        # Print results
        print("Solution status = ", probFBA.solution.get_status(), ":", )
        print(probFBA.solution.status[probFBA.solution.get_status()])
        print("Objective value  = ", probFBA.solution.get_objective_value())
        print("Iteration count = ", probFBA.solution.progress.get_num_iterations())

        if option_print:

            print("NH4", x[mnetVYTOPRalsto.reactions_id.index('R_EX_nh4_c_r')])
            print("NO3", x[mnetVYTOPRalsto.reactions_id.index('R_EX_no3_c_r')])
            print("co2_l", x[mnetVYTOPRalsto.reactions_id.index('R_EX_co2_e_l')])
            print("co2_s", x[mnetVYTOPRalsto.reactions_id.index('R_EX_co2_e_s')])
            print("co2_r", x[mnetVYTOPRalsto.reactions_id.index('R_EX_co2_e_r')])
            print('PS stem', x[mnetVYTOPRalsto.reactions_id.index('R_EX_photon_h_s')])
            print('PS leaf', x[mnetVYTOPRalsto.reactions_id.index('R_EX_photon_h_l')])

            print('Leaf growth', x[mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_LEAF_l')])
            print('Stem growth', x[mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_STEM_s')])
            print('Root growth', x[mnetVYTOPRalsto.reactions_id.index('R_BIOMASS_ROOT_r')])

            for substrate in substrates:
                print('RS import ' + substrate, x[mnetVYTOPRalsto.reactions_id.index(substrate)])

            for other_substrate in other_substrates:
                print('RS import ' + other_substrate, x[mnetVYTOPRalsto.reactions_id.index(other_substrate)])

            print('RS putrescine', x[mnetVYTOPRalsto.reactions_id.index('R_PTRCtex_rs')])

        if option_save_FBA:
            # save results
            df1 = pd.DataFrame({
                'Reaction Name': mnetVYTOPRalsto.reactions_name + mnetVYTOPRalsto.exchangereactions_name,
                'Reaction Number': range(0, len(mnetVYTOPRalsto.reactions_id + mnetVYTOPRalsto.exchangereactions_id)),
                'Reaction Id': mnetVYTOPRalsto.reactions_id + mnetVYTOPRalsto.exchangereactions_id,
                'Reaction Formula Id': mnetVYTOPRalsto.reactions_formulas + mnetVYTOPRalsto.exchangereactions_formulas,
                'FBA value': x[0:len(mnetVYTOPRalsto.reactions_name + mnetVYTOPRalsto.exchangereactions_name)],
                'Absolute FBA value:': [abs(v) for v in x[0:len(
                    mnetVYTOPRalsto.reactions_name + mnetVYTOPRalsto.exchangereactions_name)]],
                'Lower bound:': probFBA.variables.get_lower_bounds()[
                                0:len(mnetVYTOPRalsto.reactions_name + mnetVYTOPRalsto.exchangereactions_name)],
                'Upper bound:': probFBA.variables.get_upper_bounds()[
                                0:len(mnetVYTOPRalsto.reactions_name + mnetVYTOPRalsto.exchangereactions_name)]
            })

            df1 = df1.sort_values(by='Absolute FBA value:', ascending=False)
            df1.to_excel('output/FBA/FBA_VIRION_'
                         + simulation_name + '_CFU' + str(patho_density) + '.xlsx',
                         sheet_name='FBAoutput', index=False)

    except CplexError as exc:
        print(exc)

df = pd.DataFrame({
    'CFUs': cfus,
    "Growth rate Rs": objectives1,
    "Growth rate Tomato": objectives2,
    "Min Flux Sum": objectives3,
    "Status Max Biom Rs": status1,
    "Status Max Biom Tomato": status2,
    "Status Min Flux Sum": status3,
})

df.to_excel('output/FBA/FBA_VIRION_' + simulation_name + '_allCFUs.xlsx',
            sheet_name='FBAoutput', index=False)

# plot results
index_to_plot = [i for i, status in enumerate(status1) if status1[i] == 1 and status2[i] == 1 and status3[i] == 1]
cfus_to_plot = [cfus[i] for i in index_to_plot]
growth_ralsto_to_plot = [objectives1[i] / 24 for i in index_to_plot]
growth_leaf_to_plot = [objectives2[i] for i in index_to_plot]
interesting_reactions_flux_to_plot = {}
for reac_id in interesting_reactions:
    interesting_reactions_flux_to_plot[reac_id] = [interesting_reactions_fluxes[reac_id][i] for i in index_to_plot]

filename = simulation_name + '_results_' + str(nb_points) + '_points'
f = open('output/values/'+filename, "wb")
dump(index_to_plot, f)
dump(cfus_to_plot, f)
dump(growth_leaf_to_plot, f)
dump(growth_ralsto_to_plot, f)
dump(interesting_reactions_flux_to_plot, f)
dump(interesting_reactions, f)
f.close()

fig, axs = plt.subplots(5, 1, figsize=(25, 20))

# plot growth rates of plant and tomato
ax = axs[0]
ax.plot(cfus_to_plot, [objectives1[i] / 24 for i in index_to_plot], 'r-', label="GR Ralstonia (h-1)")
ax.plot(cfus_to_plot, [objectives2[i] for i in index_to_plot], 'g-', label="GR Tomato (d-1)")
ax.set_ylabel('Growth rate (d-1 for Tomato, h-1 for Ralstonia)')
ax.set_xlabel('cfus (log10(cfu/gFW Stem)')
ax.set_title("Growth rates")
ax.legend()
ax.grid(True)
ax.set_ybound(lower=0)

# plot interesting fluxes RS
ax = axs[1]
for reac_id in interesting_reactions_RS:
    ax.plot(cfus_to_plot, [interesting_reactions_fluxes[reac_id][i] for i in index_to_plot], '-', label=reac_id)

ax.set_ylabel('Flux rate (mmol.d-1.gplant-1)')
ax.set_xlabel('cfus (log10(cfu/gFW Stem)')
ax.set_title("Ralstonia fluxes")
ax.legend()
ax.grid(True)
ax.set_ybound(lower=0)

# plot interesting fluxes tomato
ax = axs[2]
for reac_id in interesting_reactions_plant:
    ax.plot(cfus_to_plot, [interesting_reactions_fluxes[reac_id][i] for i in index_to_plot], '-', label=reac_id)

ax.set_ylabel('Flux rate (mmol.d-1.gplant-1)')
ax.set_xlabel('cfus (log10(cfu/gFW Stem)')
ax.set_title("Tomato fluxes")
ax.legend()
ax.grid(True)
ax.set_ybound(lower=0)

# plot interesting fluxes xylem
ax = axs[3]
for reac_id in interesting_reactions_xylem:
    ax.plot(cfus_to_plot, [interesting_reactions_fluxes[reac_id][i] for i in index_to_plot], '-', label=reac_id)

ax.set_ylabel('Flux rate (mmol.d-1.gplant-1)')
ax.set_xlabel('cfus (log10(cfu/gFW Stem)')
ax.set_title("Xylem fluxes")
ax.legend()
ax.grid(True)
ax.set_ybound(lower=0)

# plot interesting fluxes xylem
ax = axs[4]
for reac_id in interesting_reactions_xylem_n:
    ax.plot(cfus_to_plot, [interesting_reactions_fluxes[reac_id][i] for i in index_to_plot], '-', label=reac_id)

ax.set_ylabel('Flux rate (mmol.d-1.gplant-1)')
ax.set_xlabel('cfus (log10(cfu/gFW Stem)')
ax.set_title("nitrogen-related xylem fluxes")
ax.legend()
ax.grid(True)
ax.set_ybound(lower=0)

plt.tight_layout()
plt.savefig('output/png/VIRION_' + simulation_name + '.png')
# plt.show()

plt.clf()
