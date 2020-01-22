*Process of identification*

** Structure **
1. Invoking GenOpt with opt.ini
  2. opt.ini specifies all the files that needs to be copied into the generated folders in the both of Template and Input sections.
  3. use the template to specify the varied params e.g. %Moje01% 
2. opt.ini reads DymolaWinXPplus.cfg
3. DymolaWinXPplus.cfg describes the command line called (i.e. dymosim commandline)
4. command.txt specifies optimized parameters, details of optimization algorithm and its parameters
5. After the dymosim, the post-process.py is run, which runs the objective-function.py to calculate the costs

**Process**
1. Prepare the template, preferably from dsfinal to start in steadystate
2. describe the params in command.txt
3. prepare / edit cost-function.py (or whatever the File4 is pointing to)



**Calling Dymosim**

* use dsfinal.txt instead of dsin.txt - initialization using last values of previous simulation
* dymosim -i generates dsin.txt
* 2nd column of the data struct means the value

** Objective function **
The objective function (`cost-function.py` by default) and takes the results filename as an argument

** Varying the objective function **
The file in the Template section in the command might be changed and it still should be copied with rename automatically, thus we can easily switch between different cost objectives.
E.g.

```
Simulation {
  Files {
	Template {
//    File1 = modelicaScheduleTemplate.txt;
	  File1 = dsinTemplate.txt;
	  File2 = modelicaSchedule.txt;
	  File3 = post-process.py;
	  File4 = cost-function-valsalva.py;

	}
	Input {
	  // Uncomment the line below to save the control sequence
	  // SavePath1 = "Simulation.Files.Template.Path1\\schedules";
	  File1 = dsin.txt;
	  File2 = modelicaSchedule.txt;
	  File3 = post-process.py; 
	  File4 = cost-function.py; 
```