# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 16:44:56 2023

@author: P306286
"""

import kinetics . tgm as TGM
import kinetics . thermp as TP
from kinetics . physconst import RG
import numpy as np
import matplotlib . pyplot as plt
import matplotlib

# =============================================================================
# Goal of this code is to plot the equilibrium constant versus temperature 
# for the most important reactions in the system
# =============================================================================
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




fuel = TGM.ThermoGasMix({ch4:1})
oxidizer = TGM.ThermoGasMix({co2:1})

alpha = 0.5 # mole fraction of CH4, you can change it to 0.4, see the difference of figures
Pres = 1 # unit: atm
# Temp = 600 # unit: K
X_CH4 = alpha
X_CO2 = 1-alpha
mix = X_CH4*fuel + X_CO2*oxidizer # the class allows us to use "+" and "*" to get
# mixture
mix.set_statetp(300, 1)
eqmix = mix.expand_elemsubset(tdb.values()) # expand more species into the system,


Temprange= np.linspace(800,1600,36)


eqconstantDRM=np.zeros(len(Temprange))
# methane cracking CH4-> 2H2 +C
eqconstant_methane_cracking=np.zeros(len(Temprange))
 # Boudouard reaction: 2CO <=> CO2 + C
eqconstant_boudouard= np.zeros(len(Temprange))
# Reverse carbon gasification:  CO + H2<=> H2O + C


eqconstant_rcg=np.zeros(len(Temprange))
# Reverse water gas shift: CO2 + H2 <=> CO + H2O
eqconstant_rwgs=np.zeros(len(Temprange))
# Steam methane reforming: CH4 + H2O <=> CO + 3H2
eqconstant_smr=np.zeros(len(Temprange))


for i, Temp in enumerate(Temprange):
    eqmix, _ = TGM.set_equilstate(eqmix, pres=Pres, temp=Temp)
    P=Pres
    # eqconstantDRM[i]= ((((eqmix.xfr(co))/eqmix.xfr_g*Pres)**2) *((eqmix.xfr(h2)/eqmix.xfr_g*Pres)**2))/((eqmix.xfr(ch4)/eqmix.xfr_g*Pres)*(eqmix.xfr(co2)/eqmix.xfr_g*Pres))
    eqconstantDRM[i] = np.exp(-(2*co.gort(Temp) + 2*h2.gort(Temp) - ch4.gort(Temp) - co2.gort(Temp)))

    # eqconstant_methane_cracking[i] = ((eqmix.xfr(c)/eqmix.xfr_g*Pres)*(eqmix.xfr(h2)/eqmix.xfr_g*Pres)**2) / (eqmix.xfr(ch4)/eqmix.xfr_g*Pres)
    eqconstant_methane_cracking[i]=np.exp( -(c_s.gort(Temp) + 2*h2.gort(Temp) - ch4.gort(Temp) ))
    # eqconstant_boudouard[i] = ((eqmix.xfr(co)/eqmix.xfr_g*Pres)**2) / ((eqmix.xfr(co2)/eqmix.xfr_g*Pres)*(eqmix.xfr(c)/eqmix.xfr_g*Pres))
    eqconstant_boudouard[i]= np.exp( -(c_s.gort(Temp) + co2.gort(Temp) - 2*co.gort(Temp) ))
    # eqconstant_rcg[i] = ((eqmix.xfr(co)/eqmix.xfr_g*Pres)**2) / ((eqmix.xfr(co2)/eqmix.xfr_g*Pres)*(eqmix.xfr(c)/eqmix.xfr_g*Pres))
    # eqconstant_rcg[i] = np.exp(-(2*co.gort(Temp) - co2.gort(Temp) - c_s.gort(Temp))) # other reaction, maybe also interesting
    eqconstant_rcg[i] = np.exp(-(h2o.gort(Temp) + c_s.gort(Temp)- co.gort(Temp) - h2.gort(Temp)))
    # eqconstant_rwgs[i] = ((eqmix.xfr(co)/eqmix.xfr_g*Pres)*(eqmix.xfr(h2o)/eqmix.xfr_g*Pres)) / ((eqmix.xfr(co2)/eqmix.xfr_g*Pres)*(eqmix.xfr(h2)/eqmix.xfr_g*Pres))
    eqconstant_rwgs[i] = np.exp(-(co.gort(Temp) + h2o.gort(Temp) - co2.gort(Temp) - h2.gort(Temp)))
    eqconstant_smr[i] = np.exp(-(co.gort(Temp) + 3*h2.gort(Temp) - ch4.gort(Temp) - h2o.gort(Temp)))
plt.figure(figsize=(10,5))
plt.plot(Temprange, np.log(eqconstantDRM), label='CH4 + CO2 ⇌ 2CO + 2H2')
plt.plot(Temprange, np.log(eqconstant_rwgs), label='CO2 + H2 ⇌ CO + H2O')
plt.plot(Temprange, np.log(eqconstant_methane_cracking), label='CH4 ⇌ 2H2 +C')
plt.plot(Temprange, np.log(eqconstant_rcg), label='CO + H2 ⇌ H2O + C')
plt.plot(Temprange, np.log(eqconstant_boudouard), label='2CO ⇌ CO2 + C')


# plt.plot(Temprange, np.log(eqconstant_smr), label='Steam Methane Reforming')
plt.xlabel('Temperature (K)', fontsize=12)
plt.ylabel('ln(K)', fontsize=12)
plt.legend()
plt.xlim([800,1600])
plt.savefig('C:/Users/Gebruiker/Documents/Research Project/Figures/eqconstantall.pdf', bbox_inches = 'tight')
plt.show()

# make a separate version for the presenation with chemical formulas instead of name