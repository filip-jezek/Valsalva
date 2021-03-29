# How to run the model in OpenModelica

## Installation prerequisites
Follow the prerequisites in https://github.com/filip-jezek/Valsalva/blob/master/README.md

## Loading the model
IMPORTANT!! The modelica source file formatting is overwritten just by loading into OpenModelica without prompt. Might cause some troubles when opened in two editors at once.

- Load the Physiolibrary 2.4
  - File - load - Physiolibrary/Physiolibrary/Package.mo
- Open the ADAN86.mo model
  - File - open - Valsalva/ADAN86.mo

## running the baseline simulation
- navigate to ADAN_main.SystemicTree.CardiovascularSystem and doubleclick to open. This is a base class, all changes applied here are inherited in other models.
- Select Simulation - Simulate from menu. Set simulation time (here 30s) and some reasonable interval (e.g. 0.02s). Make sure the CVODE solver is selected, as the simulation might be unstable with other solvers. Click OK. The OM editor might hang for a while, please be patient.
- Select variables to plot in right panel "Variables Browser"
- To draw X-Y plot (e.g. PV loop) select the New parametric plot button and click on the variable for X axis (e.g. V_LV and then on the variable intended for Y axis (e.g. P_LV)
- To adjust a parameter, find it in Variables browser and modify it. Click the Re-simulate button. The simulation will rerun with updated parameter value. Note, tha this change is not saved in the model.

## Running advanced use cases

## Modyfiyng the model


