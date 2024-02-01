# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 10:16:22 2023

@author: Gebruiker
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 11:05:13 2023

@author: Gebruiker
"""
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import cantera as ct
import matplotlib.pyplot as plt


import kinetics . tgm as TGM
import kinetics . thermp as TP
from kinetics . physconst import RG

# =============================================================================
# This code compares the fitted values with the experimental ones to look at the agreement between simulation and experiment
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
o2 = tdb['O2']


fuel = TGM.ThermoGasMix ({ch4:1})
oxidizer = TGM.ThermoGasMix ({co2:1})

Pres=1 # in atm
alpha=0.5
T=1223
 #K

time = np.linspace(0, 30, num=100)  # s
# import experiment data
expdata = pd.ExcelFile(r"C:\Users\Gebruiker\Documents\Research Project\DataGraphs\Expdatarightformat.xlsx")
data_Guo = expdata.parse('Guo')

plt.plot(data_Guo.loc[:9,"Residence_time[s]"],data_Guo.loc[:9,"Conversion_CH4"], 'r+', label='Experiment')



# convert conversion ration to mole fraction

# define the number of mechanisms
nrmechs = 3

# Define the gas mixture and kinetics from the obtained mechanism
gasmechs = []

# Loop over the file paths and create the ct.Solution objects
for i in range(nrmechs):
    gasmechs.append(ct.Solution(f'C:\\Users\\Gebruiker\\Documents\\Research Project\\mechsfromfit\\FittingMechGuo2012_set{i+1}.yaml'))

# gasmechs.append(ct.Solution('C:\\Users\\Gebruiker\\Documents\\Research Project\\Mechanisms\\FittingMechGuo2012_8Parameters_reversible_RWGS_1.yaml'))
# nrmechs=9
rmechs=[]
simmechs=[]
Y_in=[]    
for k in range(len(gasmechs)):
    gasmechs[k].TPX= T, ct.one_atm, f'CH4:{alpha}, CO2:{1-alpha}'
    gasmechs[k]
    rmechs.append(ct.IdealGasConstPressureReactor(gasmechs[k],energy='off'))
    simmechs.append(ct.ReactorNet([rmechs[k]]))
    Y_in.append(gasmechs[k]['CH4'].Y)	

    

# make empty array to fill with concentrations
concentrationsmech = np.zeros((len(time), gasmechs[1].n_species,nrmechs))

# add the equilibrium value
X_CH4 = alpha
X_CO2 = 1- alpha
mix = X_CH4 * fuel + X_CO2 * oxidizer
eqmix = mix.expand_elemsubset( tdb . values () )
eqmix, _ = TGM.set_equilstate(eqmix, pres=Pres, temp=T)



# also add the GRI 3.0 mechanism for comparison
gasgrimech3 = ct.Solution('gri30.yaml')
gasgrimech3.TPX = T, ct.one_atm, f'CH4:{alpha}, CO2:{1-alpha}'
rgrimech3 = ct.IdealGasConstPressureReactor(gasgrimech3,energy='off')
simgrimech3 = ct.ReactorNet([rgrimech3])
concentrationsgrimech3 = np.zeros((len(time), gasgrimech3.n_species))

  
# do the simulation over time for the 5 mechanisms and GRI Mech 

CH4conv=  np.zeros((len(time), nrmechs))                
for i,t in enumerate(time):
    # advance the solution in time
    simgrimech3.advance(time[i])
     
    concentrationsgrimech3[i, :]= gasgrimech3.X
    
    # fitted mechanisms
    # print(i)
        
    for j in range(nrmechs):
        
        
        
        simmechs[j].advance(time[i])
        concentrationsmech[i, :,j] = gasmechs[j].X
        
        Y_out = gasmechs[j]['CH4'].Y
        y_bar = 1-Y_out/Y_in[j]
        CH4conv[i,j]=y_bar
# plot the conversion ratio vs time
for m in range(nrmechs):
    plt.plot(time,CH4conv[:,m], label =f'Set {m+1}')


plt.xlabel("Residence time (s)")
plt.ylabel(r'CH$_4$ conversion ratio')
plt.title(f'T={T}K')
plt.legend()
plt.savefig (f'C:/Users/Gebruiker/Documents/Research Project/Figures/fittedvsexpGuoT={T}.pdf', bbox_inches = 'tight') # better for latex
plt.show()
#%%
# plot also for K=1273

Pres=1 # in atm
alpha=0.5
T=1273
 #K

time = np.linspace(0, 30, num=100)  # s

# import experiment data
expdata = pd.ExcelFile(r"C:\Users\Gebruiker\Documents\Research Project\DataGraphs\Expdatarightformat.xlsx")
data_Guo = expdata.parse('Guo')

plt.plot(data_Guo.loc[10:19,"Residence_time[s]"],data_Guo.loc[10:19,"Conversion_CH4"], 'r+', label='Experiment')


# convert conversion ration to mole fraction

# define the number of mechanisms
nrmechs = 3

# Loop over the file paths and create the ct.Solution objects
for i in range(nrmechs):
    gasmechs.append(ct.Solution(f'C:\\Users\\Gebruiker\\Documents\\Research Project\\mechsfromfit\\FittingMechGuo2012_set{i+1}.yaml'))


rmechs=[]
simmechs=[]
Y_in=[]    
for k in range(len(gasmechs)):
    gasmechs[k].TPX= T, ct.one_atm, f'CH4:{alpha}, CO2:{1-alpha}'
    gasmechs[k]
    rmechs.append(ct.IdealGasConstPressureReactor(gasmechs[k],energy='off'))
    simmechs.append(ct.ReactorNet([rmechs[k]]))
    Y_in.append(gasmechs[k]['CH4'].Y)	

    

# make empty array to fill with concentrations
concentrationsmech = np.zeros((len(time), gasmechs[1].n_species,nrmechs))

# add the equilibrium value
X_CH4 = alpha
X_CO2 = 1- alpha
mix = X_CH4 * fuel + X_CO2 * oxidizer
eqmix = mix.expand_elemsubset( tdb . values () )
eqmix, _ = TGM.set_equilstate(eqmix, pres=Pres, temp=T)



# also add the GRI 3.0 mechanism for comparison
gasgrimech3 = ct.Solution('gri30.yaml')
gasgrimech3.TPX = T, ct.one_atm, f'CH4:{alpha}, CO2:{1-alpha}'
rgrimech3 = ct.IdealGasConstPressureReactor(gasgrimech3,energy='off')
simgrimech3 = ct.ReactorNet([rgrimech3])
concentrationsgrimech3 = np.zeros((len(time), gasgrimech3.n_species))

  
# do the simulation over time for the 5 mechanisms and GRI Mech 

CH4conv=  np.zeros((len(time), nrmechs))                
for i,t in enumerate(time):
    # advance the solution in time
    simgrimech3.advance(time[i])
     
    concentrationsgrimech3[i, :]= gasgrimech3.X
    
    # fitted mechanisms
    # print(i)
        
    for j in range(nrmechs):
        
        
        
        simmechs[j].advance(time[i])
        concentrationsmech[i, :,j] = gasmechs[j].X
        
        Y_out = gasmechs[j]['CH4'].Y
        y_bar = 1-Y_out/Y_in[j]
        CH4conv[i,j]=y_bar
# plot the conversion ratio vs time
for m in range(nrmechs):
    # plt.plot(time, (0.5-concentrationsmech[:, gasmechs[m].species_index('CH4'),m])/0.5, label =f'Mechanism {m+1}')
    plt.plot(time,CH4conv[:,m], label =f'Set {m+1}')
# Define the gas mixture and kinetics from the obtained mechanism
gasmechs = []



plt.xlabel("Residence time (s)")
plt.ylabel(r'CH$_4$ conversion ratio')
plt.title(f'T={T}K')
plt.legend()
plt.savefig (f'C:/Users/Gebruiker/Documents/Research Project/Figures/fittedvsexpGuoT={T}.pdf', bbox_inches = 'tight') # better for latex


plt.show()
#%%

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import cantera as ct
import matplotlib.pyplot as plt


import kinetics . tgm as TGM
import kinetics . thermp as TP
from kinetics . physconst import RG
# also for K=1323
Pres=1 # in atm
alpha=0.5
T=1323
 #K

time = np.linspace(0, 32, num=100)  # s

# import experiment data
expdata = pd.ExcelFile(r"C:\Users\Gebruiker\Documents\Research Project\DataGraphs\Expdatarightformat.xlsx")
data_Guo = expdata.parse('Guo')

plt.plot(data_Guo.loc[20:29,"Residence_time[s]"],data_Guo.loc[20:29,"Conversion_CH4"], 'r+', label='Experiment')


# convert conversion ration to mole fraction

# define the number of mechanisms
nrmechs = 4

# Define the gas mixture and kinetics from the obtained mechanism
gasmechs = []


# Loop over the file paths and create the ct.Solution objects
for i in range(nrmechs):
    gasmechs.append(ct.Solution(f'C:\\Users\\Gebruiker\\Documents\\Research Project\\mechsfromfit\\FittingMechGuo2012_set{i+1}.yaml'))


rmechs=[]
simmechs=[]
Y_in=[]    
X_in=[]
for k in range(len(gasmechs)):
    gasmechs[k].TPX= T, ct.one_atm, f'CH4:{alpha}, CO2:{1-alpha}'
    gasmechs[k]
    rmechs.append(ct.IdealGasConstPressureReactor(gasmechs[k],energy='off'))
    simmechs.append(ct.ReactorNet([rmechs[k]]))
    Y_in.append(gasmechs[k]['CH4'].Y)	
    X_in.append(gasmechs[k]['CH4'].X)
    

# make empty array to fill with concentrations
concentrationsmech = np.zeros((len(time), gasmechs[1].n_species,nrmechs))





# also add the GRI 3.0 mechanism for comparison
gasgrimech3 = ct.Solution('gri30.yaml')
gasgrimech3.TPX = T, ct.one_atm, f'CH4:{alpha}, CO2:{1-alpha}'
rgrimech3 = ct.IdealGasConstPressureReactor(gasgrimech3,energy='off')
simgrimech3 = ct.ReactorNet([rgrimech3])
concentrationsgrimech3 = np.zeros((len(time), gasgrimech3.n_species))

  
# do the simulation over time for the 5 mechanisms and GRI Mech 

CH4convY=  np.zeros((len(time), nrmechs))   
CH4convX=  np.zeros((len(time), nrmechs))              
for i,t in enumerate(time):
    # advance the solution in time
    simgrimech3.advance(time[i])
     
    concentrationsgrimech3[i, :]= gasgrimech3.X
    
    # fitted mechanisms
    # print(i)
        
    for j in range(nrmechs):
        
        
        
        simmechs[j].advance(time[i])
        concentrationsmech[i, :,j] = gasmechs[j].X
        
        Y_out = gasmechs[j]['CH4'].Y
        y_bar = 1-Y_out/Y_in[j]
        
        X_out = gasmechs[j]['CH4'].X
        x_bar = 1-X_out/X_in[j]
        CH4convY[i,j]=y_bar
        CH4convX[i,j]=x_bar
# plot the conversion ratio vs time
for m in range(3):
    # plt.plot(time, (0.5-concentrationsmech[:, gasmechs[m].species_index('CH4'),m])/0.5, label =f'Mechanism {m+1}')
    plt.plot(time,CH4convY[:,m], label =f'Set {m+1}')
  

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


fuel = TGM.ThermoGasMix ({ch4:1})
oxidizer = TGM.ThermoGasMix ({co2:1})
# add the equilibrium value
X_CH4 = alpha
X_CO2 = 1- alpha
mix = X_CH4 * fuel + X_CO2 * oxidizer

# eqmix = mix.expand_elemsubset( tdb . values () )
eqmix = mix.expand_elemsubset([ch4,co2,co,h2])
eqmix, _ = TGM.set_equilstate(eqmix, pres=Pres, temp=T)   
 
plt.axhline((Y_in[0]-eqmix.mfr(ch4))/Y_in[0], color='k',linestyle=':', label='equilibrium')



plt.xlabel("Residence time (s)")
plt.ylabel(r'CH$_4$ conversion ratio')
plt.title(f'T={T}K')
plt.ylim([0,1])
plt.xlim([0,32])
plt.legend()
plt.savefig (f'C:/Users/Gebruiker/Documents/Research Project/Figures/fittedvsexpGuoT={T}.pdf', bbox_inches = 'tight') # better for latex
plt.show()

