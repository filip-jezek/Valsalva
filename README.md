# About
Model of cardiovascular model is able to describe  supine normal, tilt, valsalva maneuver and exercise of a healthy subject.

# Required

Implemented in Modelica language. Requires:

- [Modelica Standard library 3.2.3](https://github.com/modelica/ModelicaStandardLibrary/releases/tag/v3.2.3%2Bbuild.4) (Usually already a part of the Modelica environment)
- [Physiolibrary 2.4.0](https://github.com/filip-jezek/Physiolibrary/releases/tag/v2.4)

The model has not been upgraded to MSL 4.0 just to support OpenModelica. Developed in Dymola 2021, tested in OpenModelica 1.18-dev66.

# Installation
- Install a Modelica tool, e.g. a latest version (>= v1.18-dev66) of the OpenModelica (e.g. https://build.openmodelica.org/omc/builds/windows/nightly-builds/64bit/OpenModelica-latest.exe)
  - The MSL (Modelica Standard Library) 3.2.3 should be already loaded or load it manually
  - Current version is NOT compatible with MSL 4.0
- Download Physiolibrary from https://github.com/filip-jezek/Physiolibrary/releases/tag/v2.4
- Load the Physiolibrary into the Modelica tool
- Load the ADAN-86.mo model into the Modelica tool
  - If MSL 4.0 is used instead, the Modelica tool (e.g. Dymola) should prompt for an automatic model update. If this does not work, grab a MSL 3.2.3 instead from e.g. https://github.com/modelica/ModelicaStandardLibrary/releases/tag/v3.2.3%2Bbuild.4
  
# Model Simulation
To simulate the following main use-cases, run:
- Supine baseline - ADAN_main.SystemicTree.Baseline.CVS_baseline
- Valsalva maneuver - ADAN_main.SystemicTree.Valsalva.CVS_valsalva
- 60° HUT - ADAN_main.SystemicTree.Tilt.CVS_tiltable
- Exercise (90% of maximal) - ADAN_main.SystemicTree.Exercise.CVS_exercise
- Exercise (Stepping from 0 -- 100%) - ADAN_main.SystemicTree.Exercise.CVS_Exercise_stepping
- Base model for extensions - ADAN_main.SystemicTree.CardiovascularSystem

In OpenModelica, make sure the CVODE algorithm is used. Model seems unstable using other integrators. New to OpenModelica? See out [OpenModelica quick-start guide](Doc/OpenModelica.md)

