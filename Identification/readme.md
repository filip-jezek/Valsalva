# Install #

For running the identification one needs: 

1. Dassault Dymola for running simulations. Could be replaced by any Modelica tool, eg. OpenModelica. However the templates have to be changed then.

1. Java SDK
1. GenOpt - https://simulationresearch.lbl.gov/GO/download.html
  - run isntall file
2. Anaconda with python3
  - install anaconda
  - run anaconda prompt
  - type 
  
    conda env create -f env.yml
    conda activate CVS-identification

  

#Process of identification#


## Structure ##
1. Invoking GenOpt with opt.ini
  2. opt.ini specifies all the files that needs to be copied into the generated folders in the both of Template and Input sections.
  3. use the template to specify the varied params e.g. %MyParam01% in dsinTemplate.txt
2. opt.ini reads DymolaWinXPplus.cfg
3. DymolaWinXPplus.cfg describes the command line called (i.e. dymosim commandline through runSim.bat). The working directory is child of the ini location (temp folders to run in parallel)
4. command.txt specifies optimized parameters (e.g. MyParam01), details of optimization algorithm and its parameters
5. After the dymosim, the post-process.py (two dirs above the working dir) is run, which runs the objective-function.py (one dir above the working dir) to calculate the costs

##Process##
1. In Identification folder, create a new folder - e.g. TestRun
1. In TestRun/ create files opt.ini, command.txt
1. Simulate your model in Dymola, locate the dymosim file. For larger result files define top-level output aliases and translate the models with only output variables as an output.
1. Copy the dymosim and dsin.txt of your model into your TestRun folder
1. Prepare the dsinTemplate.txt from dsin.txt, or preferably from dsfinal to start in steadystate and copy it 
2. Describe the params in command.txt and add the placeholders in dsinTemplate
3. Edit cost-function.py and place it into your folder
5. Run GenOpt and open the TestRun/ini.txt file


##Calling Dymosim ##

* use dsfinal.txt instead of dsin.txt - initialization using last values of previous simulation
* dymosim -i generates dsin.txt
* 2nd column of the data struct means the value

## Objective function ##
The objective function (`cost-function.py` by default) and takes dict of variables names and numpy.ndarray

The objective function is unique to each identification folder.

## Debug ##
For visual studio code debugging, use the TestRun/Debug as a working directory and move dsres.mat file into it. It enables to debug the cost function using the same model.