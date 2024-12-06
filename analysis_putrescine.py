from pickle import *
import matplotlib.pyplot as plt

def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


def save_data(data, name):
    with open('{0}.txt'.format(name), 'w') as filehandle:
        for listitem in data:
            filehandle.write('%s\n' % listitem)

Putr_Metabolism = ['R_ORNDC_rs',
                   'R_ACODA_rs',
                   'R_ARGN_rs',
                   'R_ACOTA_rs',
                    'R_AGPR_rs',
                    'R_ACGK_rs',
                    'R_ACGS_rs',
                    'R_PTRCtex_rs',
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
                    'R_RE1537C_s',
                    'R_PMT_s',
                    'R_NMPTRCOX2_s',
                    'R_1MPYRS_s',
                    'R_PTRCOX1_s',
                    'R_ABUTD_s',
                    'R_ABTArm_s',
                    'R_4abut_tx_m__s',
                    'R_CATp_l',
                   'R_CATp_s',
                   'R_URIC_l',
                   'R_URIC_s',
                   'R_ORNabcpp_rs',
                    'R_OCBT_rs',
                    'R_CITRabcpp_rs',
                    'R_ARGSS_rs',
                    'R_ARGSL_rs',
                    'R_ARGDC_rs',
                    'R_AGMT_rs']

files = ['5-All']

index = []
cfus = []
growth_leaf = []
growth_ralsto = []
interesting_reactions_flux = []
interesting_reactions = []
index_cfus_6_5 = []
index_cfus_8_5 = []

for i in files:
    f = open('output/values/' + i +'_results_100_points', "rb")
    index.append(load(f))
    cfus.append(load(f))
    growth_leaf.append(load(f))
    growth_ralsto.append(load(f))
    interesting_reactions_flux.append(load(f))
    interesting_reactions.append(load(f))

for c, i in enumerate(files):
    index_cfus_6_5.append(cfus[c].index(closest(cfus[c], 6.5)))
    index_cfus_8_5.append(cfus[c].index(closest(cfus[c], 8.5)))

list_to_save_6_5 = []
list_to_save_6_5.append(files[0])
list_to_save_6_5.append('CFU 6.5')
list_to_save_8_5 = []
list_to_save_8_5.append(files[0])
list_to_save_8_5.append('CFU 8.5')

for i in Putr_Metabolism:
    list_to_save_6_5.append(i)
    list_to_save_8_5.append(i)
    list_to_save_6_5.append(interesting_reactions_flux[0][i][index_cfus_6_5[0]])
    list_to_save_8_5.append(interesting_reactions_flux[0][i][index_cfus_8_5[0]])

save_data(list_to_save_6_5,'output/FBA/putrescine_metabolism_6_5')
save_data(list_to_save_8_5,'output/FBA/putrescine_metabolism_8_5')

labelsx = 'log CFU/gFW'
linewidth = [3, 3]
markersize= [0, 0]
lines = ['-', ':']

plt.figure()
plt.plot(cfus[0], interesting_reactions_flux[0]['R_PTRCOX1_s'], lines[0], linewidth=linewidth[0],
         markersize=markersize[0], color='black')
plt.plot(cfus[0], interesting_reactions_flux[0]['R_PTRCOX1_l'], lines[1], linewidth=linewidth[1],
         markersize=markersize[0], color='black')
plt.title('Putrescine-related H$_{2}$O$_{2}$ production')
plt.ylabel('$\mu$mol/g/day')
plt.legend(['stem', 'leaf'], frameon=False)
plt.xlabel(labelsx)

plt.show()