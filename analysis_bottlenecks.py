from pickle import *
import matplotlib.pyplot as plt

files = ['1-PhotonsOnly', '2-PhotonsNitrogen', '3-PhotonsIron',
              '4-PhotonsTranspiration', '5-All']

Stem_Only = ['Exch_s_xyl_M_ala__L_c', 'Exch_s_xyl_M_arg__L_c', 'Exch_s_xyl_M_asn__L_c',
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
                   'R_AGMT_rs', 'R_CATp_l', 'R_CATp_s', 'R_URIC_l', 'R_URIC_s']

index = []
cfus = []
growth_leaf = []
growth_ralsto = []
interesting_reactions_flux = []
interesting_reactions = []

for i in files:
    f = open('output/values/'+i+'_results_100_points', "rb")
    index.append(load(f))
    cfus.append(load(f))
    growth_leaf.append(load(f))
    growth_ralsto.append(load(f))
    interesting_reactions_flux.append(load(f))
    interesting_reactions.append(load(f))

Main_Figure = ['R_GLNtex_rs', 'R_GLCtex_rs', 'Exch_r_xyl_M_glc__D_c', 'Exch_r_xyl_M_gln__L_c']
Supp_Figure_Plant = ['R_EX_nh4_c_r', 'R_EX_no3_c_r', 'R_EX_co2_e_l', 'R_EX_co2_e_s', 'R_EX_co2_e_r',
                         'R_EX_photon_h_s', 'R_EX_photon_h_l']
Supp_Figure_xylem = ['Exch_r_xyl_M_ala__L_c', 'Exch_r_xyl_M_arg__L_c', 'Exch_r_xyl_M_asn__L_c',
                               'Exch_r_xyl_M_asp__L_c', 'Exch_r_xyl_M_etoh_c', 'Exch_r_xyl_M_fum_c',
                               'Exch_r_xyl_M_ile__L_c',
                               'Exch_r_xyl_M_leu__L_c', 'Exch_r_xyl_M_lys__L_c', 'Exch_r_xyl_M_phe__L_c',
                               'Exch_r_xyl_M_pro__L_c', 'Exch_r_xyl_M_sucr_c', 'Exch_r_xyl_M_thr__L_c',
                               'Exch_r_xyl_M_tyr__L_c', 'Exch_r_xyl_M_val__L_c',
                     'Exch_r_xyl_M_nh4_c', 'Exch_r_xyl_M_no3_c', 'Exch_r_xyl_M_no2_c']

Supp_Figure_Patho = interesting_reactions[0]

for i in Main_Figure:
    Supp_Figure_Patho.remove(i)
for i in Supp_Figure_Plant:
    Supp_Figure_Patho.remove(i)
for i in Supp_Figure_xylem:
    Supp_Figure_Patho.remove(i)
for i in Stem_Only:
    Supp_Figure_Patho.remove(i)
for i in Putr_Metabolism:
    Supp_Figure_Patho.remove(i)

linewidth = [3, 3, 0, 0, 3, 0]
markersize= [0, 0, 2.5, 2.5, 0, 2.5]
lines = ['-', '--', '^', '.', ':', '*']
labelsy = '$\mu$mol/g/day'
labelsx = 'log CFU/gFW'

Supp_Figure_Patho_Txt = []
Supp_Figure_xylem_Txt = []

for c, i in enumerate(Supp_Figure_Patho):
    Supp_Figure_Patho_Txt.append(i[2:len(i)-6])

for c, i in enumerate(Supp_Figure_xylem):
    Supp_Figure_xylem_Txt.append(i[13:len(i)-2])

fig, axs = plt.subplots(3, 6, facecolor='w', edgecolor='k')
axs = axs.ravel()
for c, i in enumerate(Supp_Figure_xylem):
    for count, j in enumerate(files):
        axs[c].plot(cfus[count], interesting_reactions_flux[count][i], lines[count],
                    linewidth=linewidth[count],
                    markersize=markersize[count])
    axs[c].set_title(Supp_Figure_xylem_Txt[c])
axs[0].set_ylabel(labelsy)
axs[12].set_xlabel(labelsx)

fig, bxs = plt.subplots(4, 2, facecolor='w', edgecolor='k')
bxs = bxs.ravel()
for c, i in enumerate(Supp_Figure_Plant):
    for count, j in enumerate(files):
        bxs[c].plot(cfus[count], interesting_reactions_flux[count][i], lines[count],
                    linewidth=linewidth[count],
                    markersize=markersize[count])
    bxs[c].set_title(i)
bxs[0].set_ylabel(labelsy)
bxs[6].set_xlabel(labelsx)

fig, cxs = plt.subplots(3, 5, facecolor='w', edgecolor='k')
cxs = cxs.ravel()

for c, i in enumerate(Supp_Figure_Patho):
    for count, j in enumerate(files):
        cxs[c].plot(cfus[count], interesting_reactions_flux[count][i], lines[count],
                    linewidth=linewidth[count], markersize=markersize[count])
    cxs[c].set_title(Supp_Figure_Patho_Txt[c])
cxs[0].set_ylabel(labelsy)
cxs[10].set_xlabel(labelsx)

plt.figure()

ax1 = plt.subplot(611)
for c, i in enumerate(files):
    plt.plot(cfus[c], growth_ralsto[c], lines[c], linewidth=linewidth[c], markersize=markersize[c])
plt.tick_params('x', labelbottom=False)
plt.title('Growth rate - $\it{R. solanacearum}$')
plt.ylabel('hour$^{-1}$')
plt.legend(files, frameon=False)

# share x only
ax2 = plt.subplot(612, sharex=ax1)
for c, i in enumerate(files):
    plt.plot(cfus[c], growth_leaf[c], lines[c], linewidth=linewidth[c], markersize=markersize[c])
# make these tick labels invisible
plt.tick_params('x', labelbottom=False)
plt.title('Growth rate - tomato leaf')
plt.ylabel('day$^{-1}$')

# share x only
ax3 = plt.subplot(613)
for c, i in enumerate(files):
    plt.plot(cfus[c], interesting_reactions_flux[c]['R_GLNtex_rs'], lines[c],
             linewidth=linewidth[c], markersize=markersize[c])
# make these tick labels invisible
plt.tick_params('x', labelbottom=False)
plt.title('L-Gln uptake rate - $\it{R. solanacearum}$')
plt.ylabel(labelsy)

# share x only
ax4 = plt.subplot(614, sharex=ax3)
for c, i in enumerate(files):
    plt.plot(cfus[c], interesting_reactions_flux[c]['Exch_r_xyl_M_gln__L_c'], lines[c],
             linewidth=linewidth[c], markersize=markersize[c])
# make these tick labels invisible
plt.tick_params('x', labelbottom=False)
plt.title('L-Gln uptake rate - xylem')
plt.ylabel(labelsy)

# share x only
ax5 = plt.subplot(615)
for c, i in enumerate(files):
    plt.plot(cfus[c], interesting_reactions_flux[c]['R_GLCtex_rs'], lines[c],
             linewidth=linewidth[c], markersize=markersize[c])
plt.title('Glc uptake rate - $\it{R. solanacearum}$')
plt.ylabel(labelsy)
plt.tick_params('x', labelbottom=False)

ax6 = plt.subplot(616, sharex=ax5)
for c, i in enumerate(files):
    plt.plot(cfus[c], interesting_reactions_flux[c]['Exch_r_xyl_M_glc__D_c'], lines[c],
             linewidth=linewidth[c], markersize=markersize[c])
plt.title('Glc uptake rate - xylem')
plt.ylabel(labelsy)
plt.xlabel(labelsx)


plt.show()

