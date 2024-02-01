# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 12:15:53 2023

@author: Gebruiker
"""

import kinetics . tgm as TGM
import kinetics . thermp as TP
from kinetics . physconst import RG

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


# Define the gas mixture and initial conditions


# set initial conditions
Temp=1323
pres = ct.one_atm
alpha=0.5
comp= f'CH4:{alpha}, CO2:{1-alpha}'
# set kinetic paratemers
n1=1
n2=1
b1=0


# Time vector
time = np.logspace(-3,3, num=1000)  # s

# temperature values
Tempvalues = np.linspace(1200 ,1600 ,5 , endpoint = True ) # this range is based on the carbon free regime graph

# now we set the 'base' case according to the fitted variables
# we need this to determine the % change compared to the change with different input variable

# Define the gas mixture and kinetics
gas0 = ct.Solution("C:\\Users\\Gebruiker\\Documents\\Research Project\\mechsfromfit\\FittingMechGuo2012_set2.yaml")

# Create a batch reactor with the gas mixture
r0 = ct.IdealGasConstPressureReactor(gas0, energy ='off')

species = gas0.species()

# A1=74.4358e+04 # from Guo fit set 1
# Ea1= 9.63403e+07 # from Guo fit set 1

A1=372.038e+04 # from Guo fit set 2
Ea1=11.3255e+07 # from Guo fit set 2

# A1=1948.4e+04 # from Guo fit set 3
# Ea1=13.0632e+07 # from Guo fit set 3

# change of the pre-exponential factor A with resprect to fitted value
changeEa=0.1

# new kinetic parameter set with 10% change in pre-exponential factor A
R1 = ct.Reaction({"CH4": n1, "CO2": n2},{"CH4": n1-n2, "CO2": 0, "H2": 2*n2, "CO": n2*2},
                rate={"A": A1, "b": b1, "Ea": Ea1+Ea1*changeEa})

R1.allow_negative_orders = True

# now make a custom reaction for the sensitity analysis 
custom_reactions = []
custom_reactions.append(R1)
# custom_reactions.append(R2)
gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                       species=species, reactions=custom_reactions)
r = ct.IdealGasConstPressureReactor(gas,energy='off')

# a third one for comparison
R2 = ct.Reaction({"CH4": n1, "CO2": n2},{"CH4": n1-n2, "CO2": 0, "H2": 2*n2, "CO": n2*2},
                rate={"A": A1, "b": b1, "Ea": Ea1-Ea1*changeEa})

R2.allow_negative_orders = True

# now make a custom reaction for the sensitity analysis 
custom_reactions = []
custom_reactions.append(R2)
# custom_reactions.append(R2)
gas2 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                       species=species, reactions=custom_reactions)
r2 = ct.IdealGasConstPressureReactor(gas2,energy='off')
ch4conversion0=[]
ch4conversion=[]
ch4conversion2=[]

ch4concentration0=[]
ch4concentration=[]
ch4concentration2=[]


sensitivityEa=[]
sensitivity2=[]
# now we run the two equal gas mixtures over time and see what is the effect of difference in kinetic parameter
for i, t in enumerate(time):
    gas0.TPX = Temp, pres, comp
    gas.TPX = Temp, pres, comp
    gas2.TPX = Temp, pres, comp
    
    r0.syncState() #sync the state of gas
    r.syncState()
    r2.syncState()
    
    Y_in0 = gas0['CH4'].Y[0]
    Y_in = gas['CH4'].Y[0]
    Y_in2 = gas['CH4'].Y[0]
    
    sim0 = ct.ReactorNet([r0])
    sim = ct.ReactorNet([r])
    sim2 = ct.ReactorNet([r2])
    
    sim0.advance(time[i])
    sim.advance(time[i])
    sim2.advance(time[i])
    
    Y_out0 = gas0['CH4'].X[0] # CH4 mass fraction can be considered output of the model
    Y_out = gas['CH4'].X[0]
    Y_out2 = gas2['CH4'].X[0]
    
    y_bar0 = 1-Y_out0/Y_in0  #CH4 conversion can be considered output of the model
    ch4conversion0.append(y_bar0)
    
    y_bar = 1-Y_out/Y_in  #CH4 conversion can be considered output of the model
    ch4conversion.append(y_bar)
    
    y_bar2 = 1-Y_out2/Y_in2  #CH4 conversion can be considered output of the model
    ch4conversion.append(y_bar2)
    
    ch4concentration0.append(Y_out0)
    ch4concentration.append(Y_out)
    ch4concentration2.append(Y_out2)

    sensitivityEa.append(abs(Y_out-Y_out0)/abs(Ea1+Ea1*changeEa)*Ea1/ch4concentration0[0])
    # sensitivity2.append(abs(Y_out-Y_out0)/abs(Ea1+Ea1*changeEa)*Ea1/ch4concentration0[0])
    
 

fig, ax1 = plt.subplots()

color = 'tab:red'
# The first plot, data1, using the first y-axis
ax1.set_xlabel("Time (s)")
ax1.set_ylabel(r'S$^{rel}_{Ea}$',fontsize=12,color=color, rotation=0)

ax1.plot(time,sensitivityEa, color=color, label=r'S$_{Ea}$')
ax1.tick_params(axis='y', labelcolor=color)
handles1, labels1 = ax1.get_legend_handles_labels()


ax1.set_xscale('log')
ax1.tick_params(axis='y')

# Creating the second y-axis that shares the same x-axis
ax2 = ax1.twinx()

color = 'tab:green'
# The second plot, data2, using the second y-axis

ax2.plot(time,ch4concentration, label=f'Ea ={(Ea1+Ea1*changeEa):.2e} J/kmol', linestyle='--',color=color)
ax2.plot(time,ch4concentration0, label=f'Ea={Ea1:.2e} J/kmol',linestyle='--')
ax2.plot(time,ch4concentration2, label=f'Ea ={(Ea1-Ea1*changeEa):.2e} J/kmol',linestyle='--')
ax2.set_ylabel(r'Mole fraction CH$_4$')
ax2.set_xscale('log')

handles2, labels2 = ax2.get_legend_handles_labels()
fig.tight_layout()  # To ensure the right y-label is not slightly clipped
handles = handles1 + handles2
labels = labels1 + labels2
fig.legend(handles, labels,bbox_to_anchor=(0.88, .95), fontsize=8.3)
plt.savefig('C:\\Users\\Gebruiker\\Documents\Research Project\\Figures\\sensitivityGuoEa.pdf')
plt.show()


