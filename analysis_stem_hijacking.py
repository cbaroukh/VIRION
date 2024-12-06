from pickle import *
import matplotlib.pyplot as plt
import pandas as pd

def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


files = ['5-All', '6-All-Stem']

index = []
cfus = []
growth_leaf = []
growth_ralsto = []
interesting_reactions_flux = []
interesting_reactions = []
index_cfus_8_5 = []

for i in files:
    f = open('output/values/' + i + '_results_100_points', "rb")
    index.append(load(f))
    cfus.append(load(f))
    growth_leaf.append(load(f))
    growth_ralsto.append(load(f))
    interesting_reactions_flux.append(load(f))
    interesting_reactions.append(load(f))

for c, i in enumerate(files):
    index_cfus_8_5.append(cfus[c].index(closest(cfus[c], 8.5)))

Supp_Figure_xylem = ['Exch_r_xyl_M_ala__L_c', 'Exch_r_xyl_M_arg__L_c', 'Exch_r_xyl_M_asn__L_c',
                               'Exch_r_xyl_M_asp__L_c', 'Exch_r_xyl_M_etoh_c', 'Exch_r_xyl_M_fum_c',
                               'Exch_r_xyl_M_ile__L_c', 'Exch_r_xyl_M_gln__L_c', 'Exch_r_xyl_M_glc__Dc'
                               'Exch_r_xyl_M_leu__L_c', 'Exch_r_xyl_M_lys__L_c', 'Exch_r_xyl_M_phe__L_c',
                               'Exch_r_xyl_M_pro__L_c', 'Exch_r_xyl_M_sucr_c', 'Exch_r_xyl_M_thr__L_c',
                               'Exch_r_xyl_M_tyr__L_c', 'Exch_r_xyl_M_val__L_c']

Stem_Only = [
    'Exch_s_xyl_M_ala__L_c', 'Exch_s_xyl_M_arg__L_c', 'Exch_s_xyl_M_asn__L_c',
    'Exch_s_xyl_M_asp__L_c', 'Exch_s_xyl_M_etoh_c', 'Exch_s_xyl_M_fum_c',
    'Exch_s_xyl_M_ile__L_c', 'Exch_s_xyl_M_gln__L_c', 'Exch_s_xyl_M_glc__D_c',
    'Exch_s_xyl_M_leu__L_c', 'Exch_s_xyl_M_lys__L_c', 'Exch_s_xyl_M_phe__L_c',
    'Exch_s_xyl_M_pro__L_c', 'Exch_s_xyl_M_sucr_c', 'Exch_s_xyl_M_thr__L_c',
    'Exch_s_xyl_M_tyr__L_c', 'Exch_s_xyl_M_val__L_c', 'Exch_s_xyl_M_nh4_c',
    'Exch_s_xyl_M_no3_c', 'Exch_s_xyl_M_no2_c']

Stem_To_Show = []
Stem_Value_at_8_5 = []
for i in Stem_Only:
    if interesting_reactions_flux[1][i][index_cfus_8_5[1]] != 0:
        Stem_To_Show.append(i)
        Stem_Value_at_8_5.append(interesting_reactions_flux[1][i][index_cfus_8_5[1]])

for c, i in enumerate(Stem_To_Show):
    Stem_To_Show[c] = i[13:len(i)-2]
    if i[len(i)-5:] == '__L_c':
        Stem_To_Show[c] = 'L-'+i[13:len(i)-5].title()
    else:
        if i[len(i)-5:] == '__D_c':
            Stem_To_Show[c].title()
            Stem_To_Show[c] = 'D-' + i[13:len(i) - 5].title()
        else:
            Stem_To_Show[c] = i[13:len(i)-2].title()

df = pd.DataFrame(
   dict(
      names=Stem_To_Show, values=Stem_Value_at_8_5
   )
)
plt.figure()

plt.plot()
df_sorted = df.sort_values('values', ascending=[False])
plt.bar('names', 'values', data=df_sorted, color='black')
plt.xticks(rotation=90)
plt.ylabel('$\mu$mol/g/day')
plt.title('Stem to xylem fluxes at 10$^{8.5}$ CFU/gFW')

labelsx = 'log CFU/gFW'

linewidth = [3, 3]
markersize= [0, 0]
lines = [':', '-']
labelsy = '$\mu$mol/g/day'
labelsx = 'log CFU/gFW'
plt.figure()
ax1 = plt.subplot(211)
for c, i in enumerate(files):
    plt.plot(cfus[c], growth_ralsto[c], lines[c], linewidth=linewidth[c], markersize=markersize[c], color='black')
plt.tick_params('x', labelbottom=False)
plt.title('Growth rate - $\it{R. solanacearum}$')
plt.ylabel('hour$^{-1}$')
plt.legend(['Fluxes deactivated', 'Fluxes activated'], frameon=False)

# share x only
ax2 = plt.subplot(212, sharex=ax1)
for c, i in enumerate(files):
    plt.plot(cfus[c], growth_leaf[c], lines[c], linewidth=linewidth[c], markersize=markersize[c], color='black')
plt.title('Growth rate - tomato leaf')
plt.ylabel('day$^{-1}$')
plt.xlabel(labelsx)

plt.show()
