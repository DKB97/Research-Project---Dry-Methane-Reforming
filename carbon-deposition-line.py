# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 14:47:12 2023

@author: Gebruiker
"""

import kinetics . tgm as TGM
import kinetics . thermp as TP
from kinetics . physconst import RG
import numpy as np
import matplotlib . pyplot as plt
import matplotlib

# Using the code of Li Yang
# =============================================================================
# Goal of this code is to plot the carbon deposition line for different pressures vs temperature
# for the dry reforming of methane reaction
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
alph = np.linspace(0 ,1 ,301 , endpoint='on')

count = 0

Pres = [1]
first = np.linspace (600 ,1200 ,20 , endpoint = False )
second = np.linspace (1200 ,1600 ,11 , endpoint = True )
Temp = np.append( first , second )
data = np.zeros((len(Temp),len( Pres)))
data_Y = np.zeros(( len(Temp) ,len( Pres ) ,10) )
critica = 0
for j , pressure in enumerate ( Pres ):
 for i , Temperature in enumerate ( Temp ):
     for alpha in alph :
         
         # if alpha < critica : # this is for accelarating computation speed .
         #     print ( alpha )
         #     continue
         # else :
         X_CH4 = alpha
         X_CO2 = 1- alpha
         mix = X_CH4 * fuel + X_CO2 * oxidizer
         mix.set_statetp (300 , pressure )
         try :
             eqmix = mix.expand_elemsubset( tdb . values () )
             eqmix , _ = TGM . set_equilstate ( eqmix , pres = pressure , temp = Temperature )
         except :
             count +=1
             print ('without convergence', str ( count ))
         if eqmix . xfr ( c_s ) == 0:
             data [i ][ j] = 0
         else :
             data [i ][ j] = alpha
             data_Y [ i ][ j ][0] = eqmix . xfr ( ch4 ) # ch4
             data_Y [ i ][ j ][1] = eqmix . xfr ( h2o ) # h2o
             data_Y [ i ][ j ][2] = eqmix . xfr ( h2 ) # h2
             data_Y [ i ][ j ][3] = eqmix . xfr ( co2 ) # co2
             data_Y [ i ][ j ][4] = eqmix . xfr ( co ) # co
             data_Y [ i ][ j ][5] = eqmix . xfr ( c_s ) # c_s
             data_Y [ i ][ j ][6] = eqmix . xfr ( c) # c
             data_Y [ i ][ j ][7] = eqmix . xfr ( h) # h
             data_Y [ i ][ j ][8] = eqmix . xfr ( oh ) # oh
             data_Y [ i ][ j ][9] = eqmix . xfr ( o) # o
             break
     critica = alpha -0.1
     print ( 'critica = ', critica )

 # Carbon deposition line
cm = 1/2.54 # The conversion between the default units ( inches ) and cm
-figure , ax = plt . subplots ( figsize = [13* cm , (0 + 8) * cm ])
Pres = [0.1 ,0.5 ,1 ,2 ,5 ,10]


i =0
# for j, P in enumerate ( Pres ):
    # ax.plot( Temp , data [: , i], linewidth =1 , label ="{:.1f}".format(P) +' bar')
    
ax.plot( Temp , data [: , i], linewidth =1)
ax . set_ylim ([0 ,1])
ax.set_xlim([600,1600])
ax . set_xlabel (r'T$ \,\mathrm{(K)}$')
ax . set_ylabel (r'$\alpha$')
ax . text (1200 , 0.2 , 'carbon free regime')

plt . savefig ('C:/Users/Gebruiker/Documents/Research Project/Figures/carbon_deposition_line_atmpressure.pdf', bbox_inches = 'tight')
plt.show()
