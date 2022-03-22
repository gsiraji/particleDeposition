# particleDeposition
Simulation of deposition of particles in a fluid network


particlesim.m: the main program. Requires particle.m and G.mat to run
particle.m: the particle object.
graphcalculator.m: creates G.mat from the image of the porous medium. Requires the package Skel2Graph3D (https://www.mathworks.com/matlabcentral/fileexchange/43527-skel2graph-3d) and all its dependencies. Requires FindEdgeFlows.m, FindEdgeProps.m, potSolver.m, Im2Graph.m to run.
expFx.m: creates the F(x) plots from particle set.
