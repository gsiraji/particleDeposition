# particleDeposition
Simulation of deposition of particles in a fluid network

Read the preprint: 
https://doi.org/10.48550/arXiv.2210.01780

Comments and feedback: gessiraji@brandeis.edu



particlesim.m: the main program. Requires particle.m and G.mat to run

particle.m: the particle object.

graphcalculator.m: creates G.mat from the image of the porous medium. Requires the package Skel2Graph3D (https://www.mathworks.com/matlabcentral/fileexchange/43527-skel2graph-3d) and all its dependencies. Requires FindEdgeFlows.m, FindEdgeProps.m, potSolver.m, Im2Graph.m to run.

expFx.m: creates the F(x) plots from particle set.


Functions:
FindEdgeFlows.m: calculates the edge flow rate and shear stress for the network G. 

FindEdgeProps.m: calculates the edge geometrical properties (diameter).

potSolver.m: calculates the nodal pressures for the network G. 

Im2Graph.m: creates the graph G from an image.


![alt text](https://github.com/gsiraji/particleDeposition/blob/main/tiles1to655to195edit.pdf)
