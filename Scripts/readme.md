# Tools
These scripts are powerful tool to build / modify plaintext modelica source codes.

## Initialization at steady state
- [buildInitials.py](buildInitials.py) Generates new .mo model, that extends from the base one and is parametrized to start at exact moment from the .mat of the base model. The states being set are read from the [states.csv](states.csv)

## Building parameter files
Builds list of Modelica parameters as a separate model, using csv input (merging objects and values)
- BuildParameterFileFromCsv.py
- [SystemicTissueParameters.mo](SystemicTissueParameters.mo)

## Snapshot of the model results at given timepoint
- [Read_result2.py](Read_result2.py) Saves the flow values from given timepoint and recalculates them for different components. Outputs as TerminalsVenousParameters.csv .
- [readTerminalParams.py](readTerminalParams.py) Reads parameters of venous terminals for reparametrization

# Common
- [ModelicaClass.py](ModelicaClass.py) Data structure for traversing nodes for parametrization of its cubcomponents.
- TerminalDS: Data structure for the tissue (arterial terminal)

## Alter component parametrization
- [modifyModel.py](modifyModel.py) Includes new parameters into the systemic components. Not used anymore.

## Matlab scripts
- [import_txt_data.m](import_txt_data.m) - imports

## Modelica scripts
Offers direct recalculation of parameters and reloading them in the Dymola. However, Dymola does not clear cache (as of summer 2019) and is thus translating previously loaded file.
Still, usable for calling python scripts from the Modelica.
- [Recalculate_parameters.mos](Recalculate_parameters.mos)

## Not usable anymore
- ? [GetEndVals.m](GetEndVals.m) the same as already imlpemented in python