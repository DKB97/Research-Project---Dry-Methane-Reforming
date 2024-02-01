# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 11:02:27 2024

@author: Gebruiker
"""
# this code has the goal to compare gas-phase reactio with catalytic

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import cantera as ct
import matplotlib.pyplot as plt


import kinetics . tgm as TGM
import kinetics . thermp as TP
from kinetics . physconst import RG

gasDRMP = ct.Solution("C:\\Users\\Gebruiker\\Documents\\Research Project\\Mechs\\PolyMech2019.yaml")

Pres=1 # in atm
alpha=0.5
T=1323
 #K

time = np.linspace(0, 30, num=100)  # s




# define the number of sets of experimental data
nrmechs = 3

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



# also add the polymech mechanism for comparison (name refers to gri which was used first)
gasgrimech3 = ct.Solution('C:\\Users\\Gebruiker\\Documents\\Research Project\\Mechs\\PolyMech2019.yaml')
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
for m in range(1):
    plt.plot(time, concentrationsmech[:, gasmechs[m].species_index('CH4'),m], label ='Catalytic DRM')
    plt.plot(time, concentrationsgrimech3[:, gasgrimech3.species_index('CH4')], label ='Gas-phase DRM')

    # plt.plot(time,CH4convY[:,m], label =f'Set {m+1}')
    # plt.plot(X_out = gasmechs[j]['CH4'].X, label =f'Set {m+1}')
  

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
 
plt.axhline(eqmix.mfr(ch4), color='k',linestyle=':', label='Equilibrium')


plt.xlabel("Residence time (s)")
plt.ylabel(r'CH$_4$ mole fraction')
plt.xlim([0,30])
plt.title(f'T={T}K')
plt.legend()
plt.savefig ('C:/Users/Gebruiker/Documents/Research Project/Figures/gasvscat.pdf', bbox_inches = 'tight') # better for latex
plt.show()
