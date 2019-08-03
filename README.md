# Overview
This libraru contains code for simulating the motion of a system of N free-bodies interacting with each-other in accordance to some power-law. Specifically the scripts BasicErrorPlot.py, and MonteCarlo.py, were used to simulate the motion of the free-bodies in our paper titled: "Ideas from bounded confidence theory applied to dynamical networks of interacting free-bodies" which was submitted for publication to the academic journal "Frontiers in Physics" (by Gary O'Keeffe, and Ioannis Dassios).

### BasicErrorPlot.py
Contains the code used to generate the plots shown in Figure 3 in our paper.

### MonteCarlo.py
Contains the code used to generate the plots shown in Figure 5 in our paper.

# Refactoring
I recently refactored this code to make it more readable and usable in anticipation of publishing it here. I created two scripts with the refactored code:

### nBodySimulator.py
Contains a library of functions which can be used to simulate the motion of N free-bodies. 

### SimulationScripts.py
Demonstrates an example of how these refactored functions can be used. Thus, we have provided the tools for anyone interested to explore our method beyond the scope of the case-study in our paper.

