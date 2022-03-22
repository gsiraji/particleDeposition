function G2 = potSolver(G2,V,bdry_lm,G01)
%  Calculate the node potentials for all nodes in the network
%  Input:
%       G2    :=     the network (before potentials update)
%       V     :=    applied pressure gradient (potential)
%       bdry_lm:=   inlet and outlet boundaries - constant potential nodes            
%       G01   := reference network (in case nodes are disonnected in G2)
%  Output:
%       G2    :=    the network (with updated potentials)
%
% Mar 18 2022 - gessiraji@brandeis.edu 

% separate the biggest connected component
[bin,binsize] = conncomp(G2); 
idx = binsize(bin) == max(binsize);
SG00 = subgraph(G2, idx);

idx_B100 =[];  %%%%%%%%%%% separate the boundary and center nodes %%%%%%%%%%
idx_B200 = [];
idx_C00 = [];
    for iii=1:1:numnodes(SG00)
        if ismember(SG00.Nodes.ID(iii,1),bdry_lm.left) 
            idx_B100 = [idx_B100; iii];
        elseif ismember(SG00.Nodes.ID(iii,1),bdry_lm.right) 
            idx_B200 = [idx_B200; iii];
        else
            idx_C00 = [idx_C00; iii];
        end
    end
% rearrange the nodes
G2 = subgraph(SG00, [idx_B100; idx_B200; idx_C00]);   
limit1 = length(idx_B100);
limit2 = length(idx_B200);
 
I = incidence(G2);
W = double(diag(1./G2.Edges.Resistances)); % conductances
% the weighted graph Laplacian
L = I*W*I';
% set up the schur's complement of the Laplacian 
lb = limit1+limit2;
lL = length(L);
LBB = L(1:lb,1:lb);
LBC = L(1:lb,lb+1:lL);
LCB = L(lb+1:lL,1:lb);
LCC = L(lb+1:lL,lb+1:lL);
LCCinv = inv(LCC);
LS = LBB - (LBC*(LCCinv*LCB));
%%% set up the voltage vector for boundary nodes
psi_B100 = V.*ones(limit1,1); % high pressure nodes
psi_B200 = zeros(limit2,1); % zero pressure nodes
psi_B00 = [psi_B100; psi_B200];

% boundary current
JB = LS*psi_B00;
% bluk current is zero
J = zeros(lL,1);
J(1:lb) = JB;
% solve for the potentials
pot_vec = L\J;
%psi_C = -LCCinv*(LCB*psi_B); %the solution for the special case of two b%nodes
CV1 = -pot_vec(1) + V;
G2.Nodes.Potentials = CV1+pot_vec; % getting the unique solution by a shift
% use the ref graph for nodes that were not in the connected component
    for iii = 1:1:numnodes(G01)
        if ismember(G01.Nodes.ID(iii),G2.Nodes.ID)
            continue
        else
            new_node  = G01.Nodes(iii,:);
            new_node.Potentials = 0;
            G2 = addnode(G2,new_node);
        end
    end


end