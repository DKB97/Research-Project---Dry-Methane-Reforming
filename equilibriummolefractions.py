# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 09:00:30 2023

@author: Gebruiker
"""


import kinetics . tgm as TGM
import kinetics . thermp as TP
from kinetics . physconst import RG
import numpy as np
import matplotlib . pyplot as plt
import matplotlib


# =============================================================================
# Goal of this code is to plot the mole fractions of the main species for differen reactant ratios vs temp
# =============================================================================
tdb = TP.substfromdb() # ThermoSubstance instances from data base
ch4 = tdb['CH4']
co2 = tdb['CO2']
co = tdb['CO']
h2 = tdb['H2']
c = tdb['C']
h2o = tdb['H2O']
c_s = tdb['K*C'] 
h = tdb['H']
oh = tdb['OH']
o = tdb['O']

fuel = TGM.ThermoGasMix ({ch4:1})
oxidizer = TGM.ThermoGasMix ({co2:1})
alph = np.linspace(0.1 ,.9 ,5 )

Pres=1
# Pressures = np.linspace(0.1,10 ,200)

Temprange= np.linspace(400,1600,200)
data = np.zeros((len(Temprange),len( alph)))
data_Y = np.zeros(( len(Temprange) ,len( alph ) ,10) )
Temp=1200
for j , alpha in enumerate ( alph ):
    # defining the reactant ratio using alpha
    X_CH4 = alpha
    X_CO2 = 1- alpha
    mix = X_CH4 * fuel + X_CO2 * oxidizer
    mix.set_statetp (300 , Pres )
    # storing the equilibrium values for the species considered
    for i , Temp in enumerate ( Temprange ):
        try:
            eqmix = mix.expand_elemsubset( tdb . values () )
            eqmix, _ = TGM.set_equilstate(eqmix, pres=Pres, temp=Temp)
            data [i ][ j] = alpha
            data_Y [ i ][ j ][0] = eqmix . xfr ( ch4 ) 
            data_Y [ i ][ j ][1] = eqmix . xfr ( co2 )
            data_Y [ i ][ j ][2] = eqmix . xfr ( h2 )
            data_Y [ i ][ j ][3] = eqmix . xfr ( co ) 
            data_Y [ i ][ j ][4] = eqmix . xfr ( h2o ) 
            data_Y [ i ][ j ][5] = eqmix . xfr ( c_s ) 
            # data_Y [ i ][ j ][6] = eqmix . xfr ( c) # c
            # data_Y [ i ][ j ][7] = eqmix . xfr ( h) # h
            # data_Y [ i ][ j ][8] = eqmix . xfr ( oh ) # oh
            # data_Y [ i ][ j ][9] = eqmix . xfr ( o) # o
        except Exception as e:
            continue    
cm = 1/2.54 

# list of components to be considered
components=['CH_4', 'CO_2','H_2','CO','H_2O', 'C_s']
# componens names in latex
componentsr=[r'CH$_4$', r'CO$_2$',r'H$_2$',r'CO',r'H$_2$O', r'C$_s$']
fig, axs = plt.subplots(3, 2, figsize=(25*cm, 23*cm))  # Create a figure with 3 rows and 2 columns of subplots

lines = []  # To store the Line2D objects for the legend
labels = []  # To store the labels for the legend

for i in range(6):  # Only loop over the first 6 components
    ax = axs[i//2, i%2]  # Select the correct subplot
    for j, alpha in enumerate(alph):
        line, = ax.plot(Temprange, data_Y[:, j, i], linewidth=1)
        if i == 0:  # Only add the lines/labels to the legend for the first subplot
            lines.append(line)
            labels.append(r'$\alpha$ =' + "{:.1f}".format(alpha))
    ax.set_xlabel('Temperature (K)')
    ax.set_ylim([0,0.7])
    ax.set_xlim([400,1600])
    ax.set_ylabel(r'$X_{{{}}}$'.format(components[i]), fontsize=12)
    ax.text(0.9, 0.95, componentsr[i], ha='right', va='top', fontsize=12, transform=ax.transAxes)

  # Adjust the layout so everything fits, leaving space at the top for the legend
fig.legend(lines, labels, loc='lower center', ncol=len(alph), bbox_to_anchor=(0.5, 1),fontsize=12)  # create a single horizontal legend for the figure
plt.subplots_adjust(top=0.98)  # adjust the top margin to make room for the legend
plt.savefig('C:/Users/Gebruiker/Documents/Research Project/Figures/molefracallnew.pdf', bbox_inches = 'tight')
plt.show()