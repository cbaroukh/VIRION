# -*- coding:Utf8 -*-

#############################################################################
# Program Python type
# authors: Gerlin et al., 2025
#############################################################################

#############################################################################
# External functions
import parser_sbml_Ralsto as mRalsto
import parser_sbml_VYTOP as mVYTOP
import math

output_filename = "output/parser/parser_VIRION/VIRION"
#########################################################################
# patho weight
patho_density = 7


patho_weight = 1.11E-11 * 10**patho_density

sink_metabolites = ['M_ptrc_xyl']

########################################################################
# Addition of Ralstonia network to plant network

print()
print("Merging VYTOP and RS model...")
print("Ralstonia density: ", patho_density)

#Metabolites

metabolites_name = mVYTOP.metabolites_name
metabolites_id = mVYTOP.metabolites_id
metabolites_boundaryCondition = mVYTOP.metabolites_boundaryCondition
boundary_metabolites = mVYTOP.boundary_metabolites

for (i, mid) in enumerate(mRalsto.metabolites_id):
    if mid[-2:] == "_e":
        if mid[0:-2] + "_xyl" in mVYTOP.metabolites_id:
            print(mid, ": metabolite already in tomato model")
        else:
            metabolites_name.append(mid[0:-2] + " Xylem")
            metabolites_id.append(mid[0:-2] + "_xyl")
            metabolites_boundaryCondition.append(mRalsto.metabolites_boundaryCondition[i])
            boundary_metabolites.append(mid[0:-2] + "_xyl")
    else:
        metabolites_name.append(mid + " Ralstonia")
        metabolites_id.append(mid + "_rs")
        metabolites_boundaryCondition.append(mRalsto.metabolites_boundaryCondition[i])


reactions_name = mVYTOP.reactions_name + [
    reaction_name + " Ralstonia" for reaction_name in mRalsto.reactions_name]
reactions_id = mVYTOP.reactions_id + [reaction_id + "_rs" for reaction_id in mRalsto.reactions_id]
reactions_reversible = mVYTOP.reactions_reversible + mRalsto.reactions_reversible

reactions_formulas = mVYTOP.reactions_formulas

for reac_formula in mRalsto.reactions_formulas:
    new_formula = ''
    if ' <--> ' in reac_formula:
        arrow = ' <--> '
    elif ' <-->' in reac_formula:
        arrow = ' <-->'
    elif ' --> ' in reac_formula:
        arrow = ' --> '
    elif ' -->' in reac_formula:
        arrow = ' -->'
    else:
        print("Problem finding arrow")

    substrates = reac_formula.split(arrow)[0].split(' + ')
    products = reac_formula.split(arrow)[1].split(' + ')

    if substrates != ['']:
        for stoich_substrate in substrates:

            if " " in stoich_substrate:
                stoich = float(stoich_substrate.split(' ')[0])
                substrate = stoich_substrate.split(' ')[1]
            else:
                stoich = 1
                substrate = stoich_substrate

            if substrate[-2:] == "_e":
                substrate = substrate[0:-2] + "_xyl"
                stoich = stoich * patho_weight / mVYTOP.PLANT_WEIGHT
            else:
                substrate = substrate + "_rs"

            new_formula = new_formula + str(stoich) + " " + substrate + " + "

    new_formula = new_formula[0:-3] + arrow

    if products != ['']:
        for stoich_product in products:

            if " " in stoich_product:
                stoich = float(stoich_product.split(' ')[0])
                product = stoich_product.split(' ')[1]
            else:
                stoich = 1
                product = stoich_product

            if product[-2:] == "_e":
                product = product[0:-2] + "_xyl"
                stoich * patho_weight / mVYTOP.PLANT_WEIGHT
            else:
                product = product + "_rs"

            new_formula = new_formula + str(stoich) + " " + product + ' + '

    new_formula = new_formula[0:-3]

    reactions_formulas.append(new_formula)


stoichMat_rows = mVYTOP.stoichMat_rows + \
                 [row[0:-2] + "_xyl" if row[-2:] == "_e" else row + "_rs" for row in mRalsto.stoichMat_rows]
stoichMat_columns = mVYTOP.stoichMat_columns + [column + "_rs" for column in mRalsto.stoichMat_columns]
stoichMat_values = mVYTOP.stoichMat_values

temp_len_stoichMat_values = len(stoichMat_values)
stoichMat_indices_for_patho_weight = []

for i, row in enumerate(mRalsto.stoichMat_rows):
    val = mRalsto.stoichMat_values[i]
    if row[-2:] == "_e":
        new_val = val * patho_weight / mVYTOP.PLANT_WEIGHT
        stoichMat_values.append(new_val)
        stoichMat_indices_for_patho_weight.append(i + temp_len_stoichMat_values)
    else:
        stoichMat_values.append(val)

exchangereactions_id = ['Exch_' + mid for mid in boundary_metabolites]
exchangereactions_name = ['Exchange of ' + mid for mid in boundary_metabolites]
exchangereactions_reversible = [1]*len(boundary_metabolites)
exchangereactions_formulas = [mid + ' <--> ' for mid in boundary_metabolites]

exchangestoichMat_rows = ['QSSA_' + mid for mid in boundary_metabolites]
exchangestoichMat_columns = exchangereactions_id
exchangestoichMat_values = [-1]*len(boundary_metabolites)

# adding sink reactions for very specific metabolites present in VYTOP xylem
# but possibly not uptaken by the plant

exchangereactions_id = exchangereactions_id + \
                       ['Exch_sink_' + sink_metabolite for sink_metabolite in sink_metabolites]
exchangereactions_name = exchangereactions_name + \
                         ['Exchange of sink ' + sink_metabolite for sink_metabolite in sink_metabolites]
exchangereactions_reversible = exchangereactions_reversible + [1]*len(sink_metabolites)
exchangereactions_formulas = exchangereactions_formulas + \
                             [sink_metabolite + ' <--> ' for sink_metabolite in sink_metabolites]

exchangestoichMat_rows = exchangestoichMat_rows + [
    'QSSA_' + sink_metabolite for sink_metabolite in sink_metabolites]
exchangestoichMat_columns = exchangereactions_id
exchangestoichMat_values = exchangestoichMat_values + [-1]*len(sink_metabolites)


nb_exchangereac = len(exchangereactions_id)

nb_reac = len(reactions_id)
nb_metab = len(metabolites_id)

#########################################
# Saving data into .txt files

mVYTOP.save_data(metabolites_id, output_filename + "_ID_metabolites")
mVYTOP.save_data(boundary_metabolites, output_filename + "_boundary_metabolites")
mVYTOP.save_data(metabolites_boundaryCondition, output_filename + "_boundaryCondition_metabolites")

mVYTOP.save_data(reactions_id, output_filename + "_ID_reactions")
mVYTOP.save_data(reactions_name, output_filename + "_Name_reactions")
mVYTOP.save_data(reactions_reversible, output_filename + "_Reversibility_reactions")
mVYTOP.save_data(reactions_formulas, output_filename + "_Formulas_reactions")

mVYTOP.save_data(stoichMat_rows, output_filename + "_Rows_reactions")
mVYTOP.save_data(stoichMat_columns, output_filename + "_Columns_reactions")
mVYTOP.save_data(stoichMat_values, output_filename + "_Coefficients_reactions")


mVYTOP.save_data(exchangereactions_id, output_filename + "_ID_exchanged_reactions")
mVYTOP.save_data(exchangereactions_name, output_filename + "_Name_exchanged_reactions")
mVYTOP.save_data(exchangereactions_reversible, output_filename + "_Reversibility_exchanged_reactions")
mVYTOP.save_data(exchangereactions_formulas, output_filename + "_Formulas_exchanged_reactions")

mVYTOP.save_data(exchangestoichMat_rows, output_filename + "_Rows_exchanged_reactions")
mVYTOP.save_data(exchangestoichMat_columns, output_filename + "_Columns_exchanged_reactions")
mVYTOP.save_data(exchangestoichMat_values, output_filename + "_Coefficients_exchanged_reactions")

mVYTOP.save_data(stoichMat_indices_for_patho_weight, output_filename + "_stoichMat_indices_for_patho_weight")
