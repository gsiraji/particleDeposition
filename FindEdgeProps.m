function [edge_thic] = FindEdgeProps(G,dists)
% Calculate the channel diameter (pixels) for all channels in the network
%  Input:
%       G    :=     the network 
%       dist  :=    the closest nonzero pixel 
%                   
%       
%  Output:
%       edge_thic    :=    array of channel diameters in pixels
%
% Mar 18 2022 - gessiraji@brandeis.edu 

% approximate the diameter as the width at the center of the channel
twonodes = (G.Edges.EndNodes(:,:)); % find the endnodes
% find the cartesian coordinates in pixels for the center
% these are the indices for array of shortest distances to nonzero pixels
cntr_x = round((G.Nodes.comx(twonodes(:,1)) + G.Nodes.comx(twonodes(:,2)))/(2));
cntr_y = round((G.Nodes.comy(twonodes(:,1)) + G.Nodes.comy(twonodes(:,2)))/(2));
ind = sub2ind(size(dists),cntr_y,cntr_x); % use linear indices 
edge_thic = 2.*dists(ind); % calculate the diameter from this radius

    %     if edge_thic(ed) == 0
%         edge_thic(ed) = 11.967*ptm;  %%%% 'eliminating' these edges
%     end

