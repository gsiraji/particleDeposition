function [edgeQ,shear] = FindEdgeFlows(G,mu)
% Calculate the shear stress and flow rate for all edges in a network
%  Input:
%       G    = the graph/network
%       mu   = dynamic viscosity
%  Output:
%       edgeQ  = the flow rate array for all edges in G
%       shear  = the shear stress for all edges in G 
%
% Mar 18 2022 - gessiraji@brandeis.edu 

edgeQ = ones(numedges(G),1);    % initialize edge flow rate array
shear = ones(numedges(G),1);    % initialize edge shear stress array

for ed =1:1:numedges(G)         % loop over all edges  
    twonodes = (G.Edges.EndNodes(ed,:)); % get the endnodes for each edge
    twonodes_x = G.Nodes(twonodes,:).comx; % get the x-position for each endnode
    [~, max_x_idx] = max(twonodes_x); % find the endnode further down the stream
    
    edq  = (G.Nodes.Potentials(twonodes([1,2] ~= max_x_idx))...
        - G.Nodes.Potentials(twonodes(max_x_idx)))...
        ./G.Edges.Resistances(ed); % assign the correct sign to the flow
    
    edgeQ(ed) = edq; % assign the value to edge flow rate array
    shear(ed) = abs(edq*32*mu/(pi*(G.Edges.Widths(ed,1))^3)); % calculate the shear sress: ...
    % tau_w = 32*mu*q/(pi d^3) d is diameter

end
