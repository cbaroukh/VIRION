from lxml import etree

##########################
# input data

input_folder = 'input/'
strain = 'GMI1000'

filename = input_folder + 'RS' + strain + '.xml'

output_folder = 'output/parser/parser_RS/'
output_filename = output_folder + 'RS' + strain


##############################
# functions to save data in .txt files
def save_data(data, name):

    with open('{0}.txt'.format(name), 'w') as filehandle:
        for listitem in data:
            filehandle.write('%s\n' % listitem)

##################################
# reading input .xml network file

print()
print("Parsing RS metabolic network...")
# create a tree by parsing the xml file
tree = etree.parse(filename)

# get the root of the tree
root = tree.getroot()

# Initialize vectors
reactions_name = []
reactions_id = []
reactions_reversible = []
reactions_formulas = []

metabolites_id = []
metabolites_name = []
metabolites_compartment = []
metabolites_boundaryCondition = []
boundary_metabolites = []

stoichMat_rows = []
stoichMat_columns = []
stoichMat_values = []

compartments_name = []
compartments_id = []

# get and parse compartments
listOfCompartments = root[0][4]

for compartment in listOfCompartments:
    compartments_name.append(compartment.get('name'))
    compartments_id.append(compartment.get('id'))

# get genes nodes
listOfGenes = root[0][1]

# get metabolites nodes
listOfSpecies = root[0][5]

# get reactions nodes
listOfReactions = root[0][7]

# parse metabolites
for metabolite in listOfSpecies:
    metabolites_id.append(metabolite.get('id'))
    metabolites_name.append(metabolite.get('name'))

    if metabolite.get('boundaryCondition') == 'true':
        metabolites_boundaryCondition.append(1)
        boundary_metabolites.append(metabolite.get('id'))
    else:
        metabolites_boundaryCondition.append(0)

# parse reactions
for ind, reaction in enumerate(listOfReactions):
    reactions_name.append(reaction.get('name'))
    reactions_id.append(reaction.get('id'))

    if reaction.get('reversible') == "true":
        reactions_reversible.append(1)
    else:
        reactions_reversible.append(0)

    reaction_formula = ''


    listOfReactants = reaction[len(reaction)-2]
    listOfProducts = reaction[len(reaction)-1]

    for substrate in listOfReactants:
        stoichMat_rows.append('QSSA_' + substrate.get("species"))
        stoichMat_columns.append(reaction.get('id'))
        stoichMat_values.append(- float(substrate.get('stoichiometry')))
        reaction_formula = reaction_formula + substrate.get('stoichiometry') + ' ' + substrate.get('species') + ' + '

    reaction_formula = reaction_formula[0:-2]

    if reaction.get('reversible') == 'true':
        reaction_formula = reaction_formula + '<--> '
    else:
        reaction_formula = reaction_formula + '--> '

    if len(listOfProducts) == 0:
        reaction_formula = reaction_formula[0:-1]
    else:
        for product in listOfProducts:
            stoichMat_rows.append('QSSA_' + product.get("species"))
            stoichMat_columns.append(reaction.get('id'))
            stoichMat_values.append(float(product.get('stoichiometry')))
            reaction_formula = reaction_formula + product.get('stoichiometry') + ' ' + product.get('species') + ' + '

        reaction_formula = reaction_formula[0:-3]

    reactions_formulas.append(reaction_formula)

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

# get some statistics
nb_reac = len(reactions_id)
nb_metab = len(metabolites_id)
nb_exchangereac = len(exchangereactions_id)

exchangereactions_reversible = [1]*nb_exchangereac


################################################
# Printing informations on the parsed network

print('General informations on RS metabolic network')
print('Nb Metabolites: ', nb_metab)
print('Nb Reactions: ', nb_reac)
print('Nb Compartments:', len(listOfCompartments))
print('Nb Genes:', len(listOfGenes))
print('Number of boundary metabolites', len(boundary_metabolites))
print('Number of exchange reactions', nb_exchangereac)

print()


#########################################
# Saving data into .txt files

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


