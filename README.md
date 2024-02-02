The codes and data provided in this repository are used for a numerical analysis of the Dry Reforming of Methane (DRM) reaction.
- Thermodynamic equilibrium calculations use the kinetics package by Mokhov (2017). 
- For simulation kinetics, the ideal gas constant pressure reactor from Cantera is used
- The reaction mechanism (PolyMech) for simulating gas-phase DRM is developed by Porras et al. (2019)
- Three versions of a simplified catalytic reaction mechanism are included, all having a different set of kinetic paramters that followed from fitting.
  
The codes incudes path to datasets and reaction mechanism related to my own directories. Do not forget to change these paths to your own directories

Research abstract:
Dry reforming of methane is a promising reaction from a sustainability point of view
since two of the most abundant greenhouse gases, CH4 and CO2, are converted into
valuable syngas. Moreover, biogas which is already a mixture of CH4 and CO2 can be
used as fuel or CO2 can be used from carbon capture. The product of the reaction,
syngas, can be used as fuel for a SOFC or for the production of a variety of chemicals.
However, the development of widespread use of DRM has slowed down due to catalysts
issues as a consequence of carbon deposition. This research aims at finding the optimal
reaction circumstances (reactants composition and temperature) to have maximum
hydrogen yield without carbon deposition. In addition, it seeks to provide numerical
methods that can be used to simulate DRM kinetics, both gas-phase and catalytic.
By doing thermodynamic equilibrium calculations, it was shown that a temperature
of at least 1300 K and a reactant composition of 1:1 results in maximum hydrogen
production without carbon deposition. The remainder of the research focuses on the
chemical kinetics of the reaction. The chemical reaction mechanism PolyMech2019 is
used for making the simulations. Gas-phase simulations show that there is negligible
conversion of CH4 for SOFC operating conditions (max. 1200 K). The construction
of a simplified reaction mechanism consisting of only one reaction, the DRM reaction
is discussed in the chapter on heterogeneous kinetics. Experimental data is used for
fitting kinetic parameters of which the simplified mechanism consists. The resulting
simulations showed good agreement with experimental values for high temperatures
(>1300 K). Lastly, a general method for making simplified reaction mechanisms is
suggested. The results of a sensitivity analysis indicate the range of residence times
for which the system is most sensitive for changes in kinetic parameters. It is therefore
this range of residence times in which the experiments should be performed in order
to obtain the most accurate kinetic parameters from fitting.

