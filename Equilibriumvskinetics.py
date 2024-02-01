# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:18:44 2024

@author: Gebruiker
"""


import kinetics . tgm as TGM
import kinetics . thermp as TP
from kinetics . physconst import RG

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# this code has to goal of plotting equilibrium value vs gas-phase kinetic values at res time of 1000s


# Time vector
time = np.logspace(-1, 8, num=1000)  # s


# define the gas mixtures and kinetics
gasDRMP = ct.Solution("C:\\Users\\Gebruiker\\Documents\\Research Project\\Mechs\\PolyMech2019.yaml")





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

# array of all species in thermodb containing more than 5 carbon atoms
speciesc5 = ['C5', 'C5H6', 'C5H8', 'C5H10', 'C5H12', 'C6H6', 'C6H10', 'C6H12', 'C6H14', 'C7H8', 'C7H10', 'C7H12', 'C7H14', 'C7H16', 'C8H8', 'C8H14', 'C8H16', 'C8H18', 'C9H20', 'C10H8', 'C10H10', 'C10H12', 'C10H18', 'C10H22', 'C14H10', 'C14H12', 'C14H14', 'C14H18', 'C14H24', 'C16H10', 'C18H12', 'C20H12', 'C22H12', 'C22H14', 'C24H12', 'C24H14', 'C26H14', 'C28H14', 'C30H14', 'C30H16', 'C32H14']






# operating conditions
alpha=0.5
time=1000
Tempvalues= np.linspace(773,1473,100, endpoint='on')
Pres=1 # in atm

indices_C5plus = [i for i, sp in enumerate(gasDRMP.species()) if sp.composition.get('C', 0) > 4]

 # Sum the mole fractions of species with more than 5 carbon atoms


# Data vectors for storing the product concentrations
concentrationsDRM = np.zeros((len(Tempvalues),gasDRMP.n_species))
concentrationsC5=[]
equilibriumCH4=[]
equilibriumC5=[]
for i, T in enumerate(Tempvalues):
    # equilibirum
    fuel = TGM.ThermoGasMix ({ch4:1})
    oxidizerDRM = TGM.ThermoGasMix ({co2:1})
    oxidizerSMR= TGM.ThermoGasMix ({h2o:1})
    X_CH4 = alpha
    X_CO2 = 1- alpha
    X_H2O = 1- alpha
    mixDRM = X_CH4 * fuel + X_CO2 * oxidizerDRM
    eqmixDRM = mixDRM.expand_elemsubset(tdb.values())
    eqmixDRM, _ = TGM.set_equilstate(eqmixDRM, pres=Pres, temp=T)
    equilibriumCH4.append(eqmixDRM . xfr ( ch4 ))
    
    sum_mole_fractions = 0
    for specie in speciesc5:
        sum_mole_fractions += eqmixDRM.xfr(tdb[specie])
    equilibriumC5.append(sum_mole_fractions)
    
    # kinetics
    gasDRMP.TPX = T, ct.one_atm, f'CH4:{alpha}, CO2:{1-alpha}'        
    rDRM = ct.IdealGasConstPressureReactor(gasDRMP,energy='off')
    simDRM = ct.ReactorNet([rDRM])
    simDRM.advance(time)
    concentrationsDRM[i,:] = gasDRMP.X
    mole_fractions = np.sum(concentrationsDRM[i, indices_C5plus])
    concentrationsC5.append(mole_fractions)

plt.plot(Tempvalues, concentrationsDRM[:,gasDRMP.species_index('CH4')], 'r', label=r'CH$_4$')
plt.plot(Tempvalues, equilibriumCH4, 'r', linestyle='--', label=r'Equilibrium CH$_4$')
plt.plot(Tempvalues, concentrationsC5, 'b', label=r'C$^{5+}$')
plt.plot(Tempvalues, equilibriumC5, 'b', linestyle='--',label=r'Equilibrium C$^{5+}$' )
plt.ylabel('Mole fraction')
plt.xlabel('Temperature (K)')
plt.legend()
plt.savefig('C:/Users/Gebruiker/Documents/Research Project/Figures/equilibriumvskinetics.pdf', bbox_inches = 'tight')

plt.show()
