The present SPH Code reads three important input files including:

1- "data.plt," which includes the initial data of the collapse problem.

2- "materials.mat,” which contains material properties of collapsing column, walls, and bed.

3- "options.opt," which sets the input parameters for simulation including final time, CFL number, particle shifting coefficient, etc.


Notes:
	For simulation of additional test cases, please note that "data.plt" file should be created using a mesh generation software with the following format:
Variables= "Num" "PType" "X" "Y" "VX" "VY" "P" "sxx" "sxy" "syy" "m" "h" "rho" "GammaDot" "Vorticity" "Okubo_Weiss" 
Zone

0	0	2.1e-01	5.22e-01	0	0	0	0	0	0	0	0	0	0	0	0

    ...
    
    ...
    
    ... 
    

Where "Num" denotes the index of each SPH particle and starts from 0. "PType" is the type of each SPH particle (0= wall, 11= overlying fluid, and 15= bed material). "X" and "Y" denote the initial position of each SPH particle. "m" and "rho" are the mass and density of each SPH particle, respectively. "h" is the smoothing radius and equals to 2.6 times initial particle spacing. The other parameters are initially set to zero and after solving the problem, their values will be updated at different instances.

-	Runout of the flow is written in “Runout.dat” file.

-	Data of flow in different monitoring sections along the path including momentum, moving mass, etc. are written in “type1_data_X_*.dat” file, where * denotes the index of the monitoring section.  

-	Tecplot file "lay500-vector--drop.lay" reads the output file (data.plt) which is in the aforementioned format.  

-	For adjustment of slip velocity of granular materials along the walls, change the "Ls" parameter in “projection.cpp” C++ source file.


Note:
A "Solver" folder should be created with a "boost" folder located within it. Please make sure that "boost" folder of the MTL4 solver (http://old.simunova.com/en/node/189) replaces the existing empty "boost" folder. 


