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
import parser_sbml_VYTOP as mnetVYTOP  # open the parsed network

#############################################################################

#                                MAIN PROGRAM

#############################################################################

#############################################################################

#                DEFINE THE CONSTRAINTS ON THE 3 COMP MODEL

#############################################################################
# Upper bounds
ub = [cplex.infinity] * len(mnetVYTOP.reactions_id) + [cplex.infinity] * mnetVYTOP.nb_exchangereac

# Lower bounds. Set to 0 for exchange reactions to avoid infinite uptake of substrates other than those wished
lb = [-cplex.infinity if rev == 1 else 0 for rev in mnetVYTOP.reactions_reversible] \
     + [0] * mnetVYTOP.nb_exchangereac


# Additional constraints such as ATP maintenance ou leaf biomass growth
# leaf biomass growth
lb[mnetVYTOP.reactions_id.index('R_BIOMASS_LEAF_l')] = mnetVYTOP.leaf_biomass_growth

# ATP maintenance
lb[mnetVYTOP.reactions_id.index('R_ATPS_l')] = mnetVYTOP.atpm_l
lb[mnetVYTOP.reactions_id.index('R_ATPS_s')] = mnetVYTOP.atpm_s
lb[mnetVYTOP.reactions_id.index('R_ATPS_r')] = mnetVYTOP.atpm_r

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
                     'M_mobd_b_r',
                     'M_mn2_b_r',
                     'M_cobalt2_b_r',
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
                     'M_biomass_stem_b_s']

for id in infinite_exchange:
    lb[mnetVYTOP.nb_reac + mnetVYTOP.exchangereactions_id.index('Exch_'+id)] = -cplex.infinity
    ub[mnetVYTOP.nb_reac + mnetVYTOP.exchangereactions_id.index('Exch_'+id)] = cplex.infinity

fordidden_exchange = ['M_nh4_b_l', 'M_nh4_b_s']

for id in fordidden_exchange:
    lb[mnetVYTOP.nb_reac + mnetVYTOP.exchangereactions_id.index('Exch_' + id)] = 0
    ub[mnetVYTOP.nb_reac + mnetVYTOP.exchangereactions_id.index('Exch_'+id)] = 0

# [OPTIONAL] Constraint on xylem composition (taken from xylem.csv file)
if mnetVYTOP.xyl == 1:
    xylem = pd.read_csv("input/xylem.csv", sep=";")
    for num, id in enumerate(xylem.id):
        lb[mnetVYTOP.reactions_id.index(id)] = xylem.lb[num]
        ub[mnetVYTOP.reactions_id.index(id)] = xylem.ub[num]

# Right-hand side of the LP problem
rhs = [0] * mnetVYTOP.nb_metab


# Senses of the LP problem
sense = "E" * mnetVYTOP.nb_metab


# Objective function of the LP problem (minimizations of photons import, stem and leaf)
obj = [0.0] * (len(mnetVYTOP.reactions_id) + len(mnetVYTOP.exchangereactions_id))

# leaf
num_R_photons_leaf = len(mnetVYTOP.reactions_id) + mnetVYTOP.exchangereactions_id.index('Exch_M_photon_b_l')
lb[num_R_photons_leaf] = -cplex.infinity
ub[num_R_photons_leaf] = cplex.infinity

# stem
num_R_photons_stem = len(mnetVYTOP.reactions_id) + mnetVYTOP.exchangereactions_id.index('Exch_M_photon_b_s')
lb[num_R_photons_stem] = -cplex.infinity
ub[num_R_photons_stem] = cplex.infinity

obj[num_R_photons_leaf] = float(-1.0 * mnetVYTOP.LEAF_WEIGHT)
obj[num_R_photons_stem] = float(-1.0 * mnetVYTOP.STEM_WEIGHT)

# Variables' initialisation for additional constraints of the LP problem, not on lb or ub
addnames = []
addstoichMat_rows = []
addstoichMat_values = []
addstoichMat_columns = []
addrhs = []
addsense = ""

# Constraint for relative growth rates of each organ (based on experimental data)
addnames.append("Root to Leaf relative growth rate")

addstoichMat_values.append(0.17/0.21)
addstoichMat_values.append(-1)

addstoichMat_rows.append("Root to Leaf relative growth rate")
addstoichMat_rows.append("Root to Leaf relative growth rate")

addstoichMat_columns.append('R_BIOMASS_LEAF_l')
addstoichMat_columns.append('R_BIOMASS_ROOT_r')

addrhs.append(0)

addsense = addsense + "E"


addnames.append("Stem to Leaf relative growth rate")
addstoichMat_values.append(0.26/0.21)
addstoichMat_values.append(-1)

addstoichMat_rows.append("Stem to Leaf relative growth rate")
addstoichMat_rows.append("Stem to Leaf relative growth rate")

addstoichMat_columns.append('R_BIOMASS_LEAF_l')
addstoichMat_columns.append('R_BIOMASS_STEM_s')

addrhs.append(0)

addsense = addsense + "E"

# leaf stem photosynthesis constraint
addnames.append('leaf stem photosynthesis constraint')

addstoichMat_values.append(mnetVYTOP.leaf_stem_photosynthesis)
addstoichMat_values.append(-1)

addstoichMat_rows.append('leaf stem photosynthesis constraint')
addstoichMat_rows.append('leaf stem photosynthesis constraint')

addstoichMat_columns.append('R_EX_photon_h_s')
addstoichMat_columns.append('R_EX_photon_h_l')

addrhs.append(0)

addsense = addsense + "L"

# [OPTIONAL] nitrogen ratio constraint
if not np.isnan(mnetVYTOP.nh4_no3):
    addnames.append("NH4 to NO3 uptake ratio")

    addstoichMat_values.append(mnetVYTOP.nh4_no3)
    addstoichMat_values.append(-1)

    addstoichMat_rows.append("NH4 to NO3 uptake ratio")
    addstoichMat_rows.append("NH4 to NO3 uptake ratio")

    addstoichMat_columns.append('R_EX_no3_c_r')
    addstoichMat_columns.append('R_EX_nh4_c_r')

    addrhs.append(0)

    addsense = addsense + "E"

# [OPTIONAL] carboxylase oxygenase Rubisco activity ratio
if not np.isnan(mnetVYTOP.ratio_Rubisco_carboxylase_oxygenase):
    addnames.append("Rubisco CO2 vs O2 ratio activity leaf")

    addstoichMat_values.append(1)
    addstoichMat_values.append(-mnetVYTOP.ratio_Rubisco_carboxylase_oxygenase)

    addstoichMat_rows.append("Rubisco CO2 vs O2 ratio activity leaf")
    addstoichMat_rows.append("Rubisco CO2 vs O2 ratio activity leaf")

    addstoichMat_columns.append('R_RBPCh_l')
    addstoichMat_columns.append('R_RBCh_1_l')

    addrhs.append(0)

    addsense = addsense + "E"

    addnames.append("Rubisco CO2 vs O2 ratio activity stem")

    addstoichMat_values.append(1)
    addstoichMat_values.append(-mnetVYTOP.ratio_Rubisco_carboxylase_oxygenase)

    addstoichMat_rows.append("Rubisco CO2 vs O2 ratio activity stem")
    addstoichMat_rows.append("Rubisco CO2 vs O2 ratio activity stem")

    addstoichMat_columns.append('R_RBPCh_s')
    addstoichMat_columns.append('R_RBCh_1_s')

    addrhs.append(0)

    addsense = addsense + "E"


constraints = zip(mnetVYTOP.stoichMat_rows + mnetVYTOP.exchangestoichMat_rows + addstoichMat_rows,
            mnetVYTOP.stoichMat_columns + mnetVYTOP.exchangestoichMat_columns + addstoichMat_columns,
            mnetVYTOP.stoichMat_values + mnetVYTOP.exchangestoichMat_values + addstoichMat_values)

constraint_names = ["QSSA_" + id for id in mnetVYTOP.metabolites_id] + addnames

variable_names = mnetVYTOP.reactions_id + mnetVYTOP.exchangereactions_id

rhs = rhs + addrhs

sense = sense + addsense
#############################################################################

#                      FLUX BALANCE ANALYSIS
#                     1) Photon uptake minimization
#                     2) ABS(Flux) sum minimization

#############################################################################

print()
print("Solving VYTOP Min Photon FBA...")

try:
    #####################################
    # FBA 1) Photon uptake minimization
    #####################################

    probFBA = cplex.Cplex()

    probFBA.set_problem_name("FBATomatoMinPhotons")

    probFBA.objective.set_sense(probFBA.objective.sense.minimize)

    probFBA.linear_constraints.add(rhs=rhs, senses=sense,
                                   names=constraint_names)

    probFBA.variables.add(obj=obj, lb=lb, ub=ub, names=variable_names)

    probFBA.linear_constraints.set_coefficients(constraints)

    probFBA.write("output/LP/" + "FBA1-FBATomatoMinPhotons" + "_probFBA.lp")

    probFBA.solve()

    x = probFBA.solution.get_values()
    photon_stem = x[mnetVYTOP.reactions_id.index('R_EX_photon_h_s')]
    photon_leaf = x[mnetVYTOP.reactions_id.index('R_EX_photon_h_l')]

    # Print results

    print("Solution status = ", probFBA.solution.get_status(), ":",)
    print(probFBA.solution.status[probFBA.solution.get_status()])
    print("Objective value  = ", probFBA.solution.get_objective_value())
    print("Iteration count = ", probFBA.solution.progress.get_num_iterations())

    print("NH4", x[mnetVYTOP.reactions_id.index('R_EX_nh4_c_r')])
    print("NO3", x[mnetVYTOP.reactions_id.index('R_EX_no3_c_r')])
    print("co2_l", x[mnetVYTOP.reactions_id.index('R_EX_co2_e_l')])
    print("co2_s", x[mnetVYTOP.reactions_id.index('R_EX_co2_e_s')])
    print("co2_r", x[mnetVYTOP.reactions_id.index('R_EX_co2_e_r')])
    print('PS stem', x[mnetVYTOP.reactions_id.index('R_EX_photon_h_s')])
    print('PS leaf', x[mnetVYTOP.reactions_id.index('R_EX_photon_h_l')])

    print('Leaf growth', x[mnetVYTOP.reactions_id.index('R_BIOMASS_LEAF_l')])
    print('Stem growth', x[mnetVYTOP.reactions_id.index('R_BIOMASS_STEM_s')])
    print('Root growth', x[mnetVYTOP.reactions_id.index('R_BIOMASS_ROOT_r')])

    #####################################
    # FBA 2) Min Sum(abs(flux))
    #####################################

    # define new name of optimization problem
    probFBA.set_problem_name("FBATomatoRSMinFluxSum")

    # define new name of optimization problem
    probFBA.set_problem_name("FBATomatoRSMinFluxSum")

    # add constraints for photon intake from previous FBA
    probFBA.variables.set_upper_bounds('R_EX_photon_h_l', photon_leaf)
    probFBA.variables.set_upper_bounds('R_EX_photon_h_s', photon_stem)

    # define additional variables names
    variable_names2 = ['MinVar_' + id for id in mnetVYTOP.reactions_id]

    # define new obj
    probFBA.objective.set_sense(probFBA.objective.sense.minimize)
    probFBA.objective.set_linear('Exch_M_photon_b_l', 0.0)
    probFBA.objective.set_linear('Exch_M_photon_b_s', 0.0)

    obj2 = [1] * len(mnetVYTOP.reactions_id)

    # new bounds for additional variables
    lb2 = [-cplex.infinity] * len(mnetVYTOP.reactions_id)

    ub2 = [cplex.infinity] * len(mnetVYTOP.reactions_id)

    # add new variables
    probFBA.variables.add(obj=obj2, lb=lb2, ub=ub2, names=variable_names2)

    # define additional rhs
    rhs2 = [0] * len(mnetVYTOP.reactions_id) * 2

    # define additional sens
    sense2 = "L" * len(mnetVYTOP.reactions_id) * 2

    # define additional constraint names
    constraint_names2 = ["MinFluxConstraint1_" + id for id in mnetVYTOP.reactions_id] + \
                        ["MinFluxConstraint2_" + id for id in mnetVYTOP.reactions_id]

    # define additional constraints for the LP problem
    minstoichMat_rows = ["MinFluxConstraint1_" + id for id in mnetVYTOP.reactions_id] \
                        + ["MinFluxConstraint1_" + id for id in mnetVYTOP.reactions_id] \
                        + ["MinFluxConstraint2_" + id for id in mnetVYTOP.reactions_id] \
                        + ["MinFluxConstraint2_" + id for id in mnetVYTOP.reactions_id]

    minstoichMat_columns = mnetVYTOP.reactions_id \
                           + ['MinVar_' + id for id in mnetVYTOP.reactions_id] \
                           + mnetVYTOP.reactions_id \
                           + ['MinVar_' + id for id in mnetVYTOP.reactions_id]

    minstoichMat_values = [1] * len(mnetVYTOP.reactions_id) + [-1] * len(mnetVYTOP.reactions_id) * 3

    constraints2 = zip(minstoichMat_rows, minstoichMat_columns, minstoichMat_values)

    probFBA.linear_constraints.add(rhs=rhs2, senses=sense2,
                                   names=constraint_names2)

    probFBA.linear_constraints.set_coefficients(constraints2)

    probFBA.write("output/LP/" + "FBA1-TomatoMinFluxSum" + "_probFBA.lp")

    probFBA.solve()

    x = probFBA.solution.get_values()

    # Print results
    print("Solution status = ", probFBA.solution.get_status(), ":", )
    print(probFBA.solution.status[probFBA.solution.get_status()])
    print("Objective value  = ", probFBA.solution.get_objective_value())
    print("Iteration count = ", probFBA.solution.progress.get_num_iterations())

    print("NH4", x[mnetVYTOP.reactions_id.index('R_EX_nh4_c_r')])
    print("NO3", x[mnetVYTOP.reactions_id.index('R_EX_no3_c_r')])
    print("co2_l", x[mnetVYTOP.reactions_id.index('R_EX_co2_e_l')])
    print("co2_s", x[mnetVYTOP.reactions_id.index('R_EX_co2_e_s')])
    print("co2_r", x[mnetVYTOP.reactions_id.index('R_EX_co2_e_r')])
    print('PS stem', x[mnetVYTOP.reactions_id.index('R_EX_photon_h_s')])
    print('PS leaf', x[mnetVYTOP.reactions_id.index('R_EX_photon_h_l')])

    print('Leaf growth', x[mnetVYTOP.reactions_id.index('R_BIOMASS_LEAF_l')])
    print('Stem growth', x[mnetVYTOP.reactions_id.index('R_BIOMASS_STEM_s')])
    print('Root growth', x[mnetVYTOP.reactions_id.index('R_BIOMASS_ROOT_r')])

    # save results

    df1 = pd.DataFrame({
        'Reaction Name': mnetVYTOP.reactions_name + mnetVYTOP.exchangereactions_name,
        'Reaction Number': range(0, len(mnetVYTOP.reactions_id + mnetVYTOP.exchangereactions_id)),
        'Reaction Id': mnetVYTOP.reactions_id + mnetVYTOP.exchangereactions_id,
        'Reaction Formula Id': mnetVYTOP.reactions_formulas + mnetVYTOP.exchangereactions_formulas,
        'FBA value': x[0:len(mnetVYTOP.reactions_name + mnetVYTOP.exchangereactions_name)],
        'Absolute FBA value:': [abs(v) for v in x[
                                                0:len(mnetVYTOP.reactions_name + mnetVYTOP.exchangereactions_name)]],
        'Lower Bound:': probFBA.variables.get_lower_bounds()[
                        0:len(mnetVYTOP.reactions_name + mnetVYTOP.exchangereactions_name)],
        'Upper Bound:': probFBA.variables.get_upper_bounds()[
                        0:len(mnetVYTOP.reactions_name + mnetVYTOP.exchangereactions_name)]
    })

    val_atp = str(mnetVYTOP.cost)

    df1 = df1.sort_values(by='Absolute FBA value:', ascending=False)
    df1.to_excel('output/FBA/FBA_VYTOP.xlsx',
                 sheet_name='FBAoutput', index=False)

except CplexError as exc:
    print(exc)
