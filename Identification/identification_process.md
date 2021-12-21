Notes: All python calculations are done from within the folder Identification, with launch.json including working directory, e.g.:

```
"cwd" : "..\\Identification\\SimulationEvaluation\\"
```

1. Run the baseline optimization
  a. Prepare the optim with baroreflex OFF
    1. Select the folder for identification, e.g. optimizeBaselinTriSegLumens for baseline. The folder must have opt.ini in it.
    2. Make sure your cost_function.py includes all the targets
    3. Prepare the params_for_optim.txt by listing all the optimized parameters
    4. Build the model in Dymola - creates the dymosim.exe and dsin.txt
    5. Run the manipulate_dsin.py to generate the opt_command and dsin_template.txt
  b. Run the GenOpt optim 
  c. Once finished, run the post_process_optim to generate the Modelica model with optimized params
  d. Sweep the baroreflex rising rate to get a steady state
  e. Initialize in steady state
    1. Get all the model states into states.txt
    2. Run the build_initials.py to get the steady state init
2. Run the combined optimization
  a. Based on the steady state, regenerate all the FMIs for the combined optimization
  b. Copy the settings from baseline into the combined model
  c. prepare the optim as in the baseline, yet for combined cost_function
  d. run the post_process_optim
