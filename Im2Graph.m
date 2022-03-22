function [G0, Nodes0, Edges0, dist] = Im2Graph(Im)
% Calculate the network given the image of the porous medium
%  Input:
%       Im    :=    image of the medium/device
%       
%  Output:
%       G0    :=    the network as calculated by Skel2Graph3D 
%       Nodes0:=    the nodes/junctions in the network as calculated by...
%                   Skel2Graph3D
%       Edges0:=    the edges/channels in the network as calculated by...
%                   Skel2Graph3D
%       dist  :=    the closest nonzero pixel ...
%                   (used for calculating the diameter of channels later)
%
% Mar 18 2022 - gessiraji@brandeis.edu 


ThR=graythresh(Im)/1.5; % find the gray threshold then binarize the image 
BW = imbinarize(Im,ThR); % binarize the image

C = bwskel(BW(:,:,1)); % create the image skeleton from the binary image

dist = bwdist(imcomplement(BW(:,:,1))); % calculate the distance to the...
% nearest bead from each pixel

% use the package 'Skel2Graph' to calculate the graph corresponding to
% the skeleton
[G0,Nodes0, Edges0] = Skel2Graph3D(C,100);

