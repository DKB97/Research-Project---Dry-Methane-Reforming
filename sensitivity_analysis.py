#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 14:55:30 2024

@author: p306286
"""
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
import time
# sensitivity analysis for poly Mech at residence time of 10 s
time_start = time.time()

gas = ct.Solution("C:\\Users\\Gebruiker\\Documents\\Research Project\\Mechs\\PolyMech2019.yaml")

temp = 1373
pres = ct.one_atm

# gas.TPX = temp, pres, 'CH4:1, H2O:1'
gas.TPX = temp, pres, 'CH4:1, CO2:1'
r = ct.IdealGasConstPressureReactor(gas, name = '', energy='off')
# r = ct.IdealGasConstPressureReactor(gas, name = 'R1', energy='off')
sim = ct.ReactorNet([r])

no_rxn = gas.n_reactions

for i in range(no_rxn):
    r.add_sensitivity_reaction(i)


sim.advance(10) # change the reaction time

maxi = [0]*no_rxn

for j in range (no_rxn):
    s = sim.sensitivity('CH4',j)	#sensitivity of ... to reaction 2
    maxi[j] = s
 

# Sorting
maxi_1 = np.copy(maxi)
maxi_1.sort() # sorting by value
left_ptr = 0
right_ptr = len(maxi_1) - 1
ctr = 0
final = []
while (ctr < no_rxn):
    if np.abs(maxi_1[left_ptr]) > np.abs(maxi_1[right_ptr]): # Sorting by absolute value
        final.append(maxi_1[left_ptr])
        left_ptr += 1
        ctr  = 	ctr + 1
    else:
        final.append(maxi_1[right_ptr])
        right_ptr = right_ptr - 1
        ctr  = ctr + 1

# Making a plot
ctr = 0
final1 = final[0:10]
final1.reverse()
rxn = []

while (ctr < 10):
    index = maxi.index(final1[ctr])
    rxn.append(sim.sensitivity_parameter_name(index).replace(": ",""))
    ctr += 1

plt.figure(1)
plt.barh(rxn,final1)
plt.xlabel(r'CH$_4$ sensitivity')
plt.ylabel('Reactions')
plt.tick_params(left = False)
plt.tight_layout()
plt.savefig('C:/Users/Gebruiker/Documents/Research Project/Figures/sensitivity_H2_1600.pdf',bbox_inches = 'tight')
plt.show()
time_end = time.time()
print('time consuming = ', time_end-time_start)