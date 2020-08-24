# BELLO
BELLO code, a post-processing script-tool created for the automatic analysis and extraction of structural characteristics of disordered and amorphous systems. BELLO is agnostic to the code that generated single configurations or trajectories. Its capabilities include calculation of order parameter q, folded structure identification and statistics,  and detailed atomic coordination number and pair/angle-distribution functions.

to run it:

1- in the trajectory file, remove all atom numbers at the start of each frame. also change all the names "Se" to "Te"

2- rename the trajectory file to "initial.txt"

3- run "Tetrahedral finder Boundry condition" code (for applying boundry condition to the system)

4- run "Tetrahedral finder with angles" code 

5/6- run "angle-sorter" & "locality stats (detailed coordination number)" 
