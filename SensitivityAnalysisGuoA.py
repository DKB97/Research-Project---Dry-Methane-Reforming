# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 13:28:59 2023

@author: Gebruiker
"""

import kinetics . tgm as TGM
import kinetics . thermp as TP
from kinetics . physconst import RG

import pandas as pd
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

k = A1*np.exp(-Ea1/((ct.gas_constant*Temp)))
print(k)

# change of the pre-exponential factor A with resprect to fitted value
changeA=0.1

# new kinetic parameter set with changeA% change in pre-exponential factor A
R1 = ct.Reaction({"CH4": n1, "CO2": n2},{"CH4": n1-n2, "CO2": 0, "H2": 2*n2, "CO": n2*2},
                rate={"A":A1+ A1*changeA, "b": b1, "Ea": Ea1})

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
                rate={"A": A1-A1*changeA, "b": b1, "Ea": Ea1})

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

sensitivityA=[]
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

    sensitivityA.append((abs(Y_out-Y_out0))/abs(A1*changeA-A1)*A1/ch4concentration0[0])
    # sensitivity2.append(abs(y_bar-y_bar0)/abs(Ea1*changeA-Ea1))
    # sensitivity2.append(((abs(Y_out-Y_out0)/abs(A1+A1*changeA)*A1/ch4concentration0[0] + abs(Y_out-Y_out2)/abs(A1-A1*changeA)*A1/ch4concentration2[0]))/2)

fig, ax1 = plt.subplots()
fig, ax1 = plt.subplots()

color = 'tab:red'
# The first plot, data1, using the first y-axis
ax1.set_xlabel("Time (s)")
ax1.set_ylabel(r'S$^{rel}_A$',fontsize=12,rotation=0,color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.plot(time,sensitivityA, color=color, label=r'S$_{A}$')
handles1, labels1 = ax1.get_legend_handles_labels()


ax1.set_xscale('log')
ax1.tick_params(axis='y')

# Creating the second y-axis that shares the same x-axis
ax2 = ax1.twinx()

color = 'tab:green'
# The second plot, data2, using the second y-axis
ax2.set_ylabel(r'Mole fraction CH$_4$')
ax2.plot(time,ch4concentration, label=f'A ={(A1*changeA+A1):.2e} m$^3$/kmol/s', linestyle='--',color=color)
ax2.plot(time,ch4concentration0, label=f'A={A1:.2e} m$^3$/kmol/s',linestyle='--')
ax2.plot(time,ch4concentration2, label=f'A ={(A1-A1*changeA):.2e} m$^3$/kmol/s',linestyle='--')
ax2.set_xscale('log')

handles2, labels2 = ax2.get_legend_handles_labels()
fig.tight_layout()  # To ensure the right y-label is not slightly clipped
handles = handles1 + handles2
labels = labels1 + labels2
fig.legend(handles, labels,bbox_to_anchor=(0.89, .95), fontsize=8.3)
plt.savefig('C:\\Users\\Gebruiker\\Documents\Research Project\\Figures\\sensitivityA.pdf')
plt.show()


    
# plt.plot(time,ch4concentration, label=f'A ={(A1*changeA)/10**3:.2f} cm$^3$/mol/s')
# plt.plot(time,ch4concentration0, label=f'A={A1/10**3:.2f} cm$^3$/mol/s')
# plt.plot(time,ch4concentration2, label=f'A ={(A1/changeA)/10**3:.2f} cm$^3$/mol/s')
# plt.xlabel("Time (s)", fontsize=12)
# plt.ylabel('Mole fraction CH4', fontsize=12)
# plt.xscale('log')
# plt.legend()
# plt.show()    
#%%
# now we plot the sensitivity and the experimental datapoints
expdata = pd.ExcelFile(r"C:\Users\Gebruiker\Documents\Research Project\DataGraphs\Expdatarightformat.xlsx")
data_Guo = expdata.parse('Guo')
timeGuo = data_Guo.loc[:,"Residence_time[s]"]

zeros=np.zeros(len(data_Guo.loc[20:29,"Residence_time[s]"]))

fig, ax1 = plt.subplots()

color = 'tab:red'
# The first plot, data1, using the first y-axis
ax1.set_xlabel("Time (s)")
ax1.set_ylabel(r'S$^{\text{rel}}$',fontsize=12, rotation=0)
ax1.plot(time,sensitivityA,'b', label=r'A')
ax1.plot(time,sensitivityEa,'r', label=r'E$_a$')


heights = np.interp(data_Guo.loc[20:29,"Residence_time[s]"], time,sensitivityEa)
ax1.vlines(data_Guo.loc[20:29,"Residence_time[s]"],0, heights,'g', linewidth=0.8,label='Exp. res. time')
# for xc in data_Guo.loc[20:29,"Residence_time[s]"]:
#     plt.axvlines(x=xc, color='g', clip_on=True, label='Experiment res. time')
# ax1.plot(data_Guo.loc[20:29,"Residence_time[s]"],zeros, 'g+', label='Experiment')
handles1, labels1 = ax1.get_legend_handles_labels()


ax1.set_xscale('log')
ax1.tick_params(axis='y')
plt.ylim([0,0.23])
plt.xlim([10**-3,10**3])
plt.legend()
plt.savefig('C:\\Users\\Gebruiker\\Documents\Research Project\\Figures\\sensitivityvsfittedGuoA.pdf')
plt.show()






