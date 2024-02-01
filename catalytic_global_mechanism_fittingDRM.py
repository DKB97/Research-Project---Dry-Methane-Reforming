import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import cantera as ct
import matplotlib.pyplot as plt

# code used for fitting the kinetic parameters from experimental data

expdata = pd.ExcelFile(r"C:\Users\Gebruiker\Documents\Research Project\DataGraphs\Expdatarightformat.xlsx")

data_Guo = expdata.parse('Guo')

pres = ct.one_atm  # constant pressure [Pa]

# putting all the datasets in one array

data_array = []


# most optimal dataset:

data_array = [('Guo', data_Guo)]


counter = 0
def NewMech(data, A11, Ea11):
# def NewMech(x, A1, b1, A2, b2, A3, b3, A4, b4, n1, n2, n3, n4, n5, n6, n7, n8):
    # Make variables in a same Magnitude ??
    A1 = A11*1e4 # So that the fitted value is in 10 cm^3/mol/s, cantera input is in m^3/kmol/s, have a look at the unit m³/kmol/s
    # flowchart, how to exclude useless data points. that have measurement error that contribute negatively
    # b1 = b11
    b1 = 0
    Ea1 = Ea11*1e7 # So that the fitted value is in  10 kJ/mol whereas the standard cantera unit is J/kmol
    
    # since the optimal composition is 1:1, we assume reactant orders are 1 and 1
    n1 = 1
    n2 = 1
    global counter # counte the number of iterations
    counter += 1
    if counter%10 == 0:
        print('iteration counter: ',counter)
    mech = "gri30.yaml"
    gas0 = ct.Solution(mech)
    species = gas0.species()
    custom_reactions = []
    # R1 = ct.Reaction(equation=f"{}CH4 + CO2 <=> 2 CO + 2 H2",
    #                 rate={"A": A1, "b": b1, "Ea": Ea1}) # Ea J/kmol)
    R1 = ct.Reaction({"CH4": n1, "CO2": n2},{"CH4": n1-n2, "CO2": 0, "H2": 2*n2, "CO": n2*2},
                    rate={"A": A1, "b": b1, "Ea": Ea1})
    # R1.orders={'CH4': n1, 'CO2': n2}
    R1.allow_negative_orders = True

    custom_reactions.append(R1)
    # custom_reactions.append(R2)
    gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                       species=species, reactions=custom_reactions)
    r = ct.IdealGasConstPressureReactor(gas,energy='off')
    # sim = ct.ReactorNet([r])
    # define time, space, and other information vectors
    Temp = data.loc[:,"Temperature[K]"] # react temperature [K]
    comp = data.loc[:,"Composition"] # composition
    time = data.loc[:,"Residence_time[s]"] # residence time
    
    value = []
    for i in data.index:
        # print(i)
        gas.TPX = Temp[i], pres, comp[i]
        r.syncState() #sync the state of gas
        Y_in = gas['CH4'].Y[0]
        sim = ct.ReactorNet([r])
        # print(r.thermo.X)
        # print("time: ", sim.time)
        sim.advance(time[i])
        Y_out = gas['CH4'].Y[0]
        y_bar = 1-Y_out/Y_in
        # print("time: ", sim.time)
        value.append(y_bar) #[0],this is an 1x1 array, should use [0] to get number inside this array, otherwise there will be error due to data format problem.
        # print("value: ", value)
    return value


#%%
# create arrays containing fitting parameters for the mechanisms only for R^2>0.9
fitparsGuo=[]
for j, l in enumerate(data_array):
    success = False 
    for i in range(5):
        success = False 
        for k in range (2):
            success = False
            try:
                popt, pcov = curve_fit(NewMech, data_array[j][1], data_array[j][1].loc[:,"Conversion_CH4"], p0=np.asarray([(k*3+2)*10**i,10]), bounds=((10**i,5), (10**(i+1),15)), method='dogbox')
                
                success = True  # If fitting is successful, set the flag to True
            except Exception as e:
                print(f"An error occurred: {e}")
                continue
            if success:
                fittingData = NewMech(data_array[j][1], *popt)
                compare = np.stack((data_array[j][1].loc[:,"Temperature[K]"], data_array[j][1].loc[:,"Conversion_CH4"], fittingData), axis=1)
                column_values = ['Temperature[K]', 'X_CH4', 'X_CH4_Fit']
                df = pd.DataFrame (data = compare, 
                                    # index = index_values, 
                                  columns = column_values)
                df['error'] = 1 - df.loc[:,"X_CH4_Fit"]/df.loc[:,"X_CH4"]
                y_fit = df.loc[:,"X_CH4_Fit"]
                # residual sum of squares
                ss_res = np.sum((data_array[j][1].loc[:,"Conversion_CH4"] - y_fit) ** 2)
                
                # total sum of squares
                ss_tot = np.sum((data_array[j][1].loc[:,"Conversion_CH4"] - np.mean(data_array[j][1].loc[:,"Conversion_CH4"])) ** 2)
                
                # r-squared
                r2 = 1 - (ss_res / ss_tot)
                # only selected fitted values with R squared higher than 0.9
                if r2>0.9:
                    fitparsGuo.append((data_array[j][0], popt, 'R^2 = ', r2,'error^2 = ',sum(df['error']**2)))




