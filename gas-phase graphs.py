# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 11:33:47 2024

@author: Gebruiker
"""

import kinetics . tgm as TGM
import kinetics . thermp as TP
from kinetics . physconst import RG

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# goal of this code is to make kinetic plots for different SOFC operating conditions of DRM
# of the main species and also include C5+ 
# =============================================================================

import kinetics . tgm as TGM
import kinetics . thermp as TP
from kinetics . physconst import RG

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines



# Time vector
time = np.logspace(-1, 3, num=1000)  # s


# define the gas mixtures and kinetics
gasDRMP = ct.Solution("C:\\Users\\Gebruiker\\Documents\\Research Project\\Mechs\\PolyMech2019.yaml")

# Data vectors for storing the product concentrations
concentrationsDRM = np.zeros((len(time), gasDRMP.n_species))

# alp = np.linspace(0.4 ,0.5 ,5 , endpoint='on')


# for equilibrium calculations
tdb = TP.substfromdb() # ThermoSubstance instances from data base
ch4 = tdb['CH4']
co2 = tdb['CO2']
co = tdb['CO']
h2 = tdb['H2']
c = tdb['C']
h2o = tdb['H2O']
c_s = tdb['K*C'] #k is the first letter of the word "condensed" in Russia
h = tdb['H']
oh = tdb['OH']
o = tdb['O']
o2 = tdb['O2']

# alp = np.linspace(0.4 ,0.5 ,5 , endpoint='on')


fuel = TGM.ThermoGasMix ({ch4:1})
oxidizerDRM = TGM.ThermoGasMix ({co2:1})
oxidizerSMR= TGM.ThermoGasMix ({h2o:1})
alpha=0.5

Pres=1 # in atm

# Tempvalues=[1100,1300,1500]
Tempvalues=[1173,1373,1573]
alphas=[0.25,0.5,0.75]
# a list of elements to be plotted
elements=['CH4','CO2', 'CO', 'H2', 'H2O', 'C5+'] # C5+ stands for the hydrocarbons with +5C 
elementslatex=[r'CH$_4$',r'CO$_2$', 'CO', r'H$_2$', r'H$_2$O', 'C5+'] # C5+ stands for the hydrocarbons with +5C 


fig, axs = plt.subplots(len(elements), len(Tempvalues), figsize=(15, 20))

lines = []  # list to store the lines for the legend
labels = []  # list to store the labels for the legend
indices_C5plus = [i for i, sp in enumerate(gasDRMP.species()) if sp.composition.get('C', 0) > 4]
min_mole_fractions = [np.inf] * len(elements)
max_mole_fractions = [-np.inf] * len(elements)

horizontal_line = mlines.Line2D([], [], color='black', linestyle='dotted', label='Carbon deposition threshold')

for j, T in enumerate(Tempvalues):
    for i, alpha in enumerate(alphas):
        X_CH4 = alpha
        X_CO2 = 1- alpha
        X_H2O = 1- alpha
        mixDRM = X_CH4 * fuel + X_CO2 * oxidizerDRM
        eqmixDRM = mixDRM.expand_elemsubset(tdb.values())
        eqmixDRM, _ = TGM.set_equilstate(eqmixDRM, pres=Pres, temp=T)
        
        gasDRMP.TPX = T, ct.one_atm, f'CH4:{alpha}, CO2:{1-alpha}'
        
        rDRM = ct.IdealGasConstPressureReactor(gasDRMP,energy='off')
        simDRM = ct.ReactorNet([rDRM])
    
        for i,t in enumerate(time):
            simDRM.advance(time[i])
            concentrationsDRM[i,:] = gasDRMP.X
            
            
        for idx, element in enumerate(elements):
            
            
            if element == 'C5+':
                
                 
                # Sum the mole fractions of species with more than 5 carbon atoms
                mole_fractions = np.sum(concentrationsDRM[:, indices_C5plus], axis=1)
                line, = axs[idx, j].plot(time, mole_fractions)
                
                # Add a black dotted line at mole fraction 0.001
                axs[idx, j].axhline(y=0.001, color='black', linestyle='dotted', label='Mole fraction = 0.001')
                
                # Update the minimum and maximum mole fractions for this row
                min_mole_fractions[idx] = min(min_mole_fractions[idx], np.min(mole_fractions))
                max_mole_fractions[idx] = max(max_mole_fractions[idx], np.max(mole_fractions))
            else:
                mole_fractions = concentrationsDRM[:,gasDRMP.species_index(element)]
                line, = axs[idx, j].plot(time, mole_fractions)
                
                # Update the minimum and maximum mole fractions for this row
                min_mole_fractions[idx] = min(min_mole_fractions[idx], np.min(mole_fractions))
                max_mole_fractions[idx] = max(max_mole_fractions[idx], np.max(mole_fractions))

        # create the labels for alpha so that one legend can be made
        if j == 0:    
            
            lines.append(line)
            labels.append(fr'$\alpha$ = {alpha}')
            
            # Set the y-axis range for each row
    for idx in range(len(elements)):
        for j,T in enumerate(Tempvalues):
            axs[idx, j].set_ylim([min_mole_fractions[idx], max_mole_fractions[idx]*1.1])

     
            axs[idx, j].set_title(f'{elementslatex[idx]} at T={T} K', fontsize=15)
            axs[idx, j].set_xlabel("Time (s)", fontsize=15)
            axs[idx, j].set_ylabel(r'Mole fraction', fontsize=15)
            axs[idx, j].set_xscale('log')
            axs[idx, j].set_xlim([min(time), max(time)])
        
            
lines.append(horizontal_line)
labels.append(horizontal_line.get_label())
plt.tight_layout()
fig.legend(lines, labels, loc='lower center', ncol=len(lines), bbox_to_anchor=(0.5, 1),fontsize=18)  # create a single horizontal legend for the figure
plt.subplots_adjust(top=0.98)  # adjust the top margin to make room for the legend
plt.savefig('C:/Users/Gebruiker/Documents/Research Project/Figures/gasphasepolymech.pdf', bbox_inches = 'tight')
plt.show()
