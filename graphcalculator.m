%% particlesim_firstpart.m 
% needs Im2Graph.m, Skel2Graph3D.m, potSolver.m, FindEdgeFlows.m, and more* to run
% * the Skel2Graph3D.m has other dependencies
Im = imread('Beads.tif');Im = Im(:,1903:end,:);

V1 = 160000;

ptm = (1e-6)/1.61; % pixels to meters

w = size(Im,1);l = size(Im,2);boundary_limits = [50, l-54];
mu = 60e-3;depth = 1e-3;p_dmeter = 1e-6;
Nstep = 499;tot_time = 250;dt = tot_time/Nstep; % dt has to be < 0.002

% make a graph from the bead image
[G1, Nodes1, Links1, dists] = Im2Graph(Im);
%%
G0 = graph(G1);
G0.Nodes.comx = [Nodes1.comy]';
G0.Nodes.comy = [Nodes1.comx]';
G0.Edges.Widths = ptm.*FindEdgeProps(G0,dists);
%%
G0 = rmedge(G0, find(G0.Edges.Widths(:,1) < p_dmeter));
%%

% create the connected subgraph
[bin,binsize] = conncomp(G0);
idx = binsize(bin) == max(binsize);
SG = subgraph(G0, idx);
SG.Nodes.comx = ptm.*SG.Nodes.comx;
SG.Nodes.comy = ptm.*SG.Nodes.comy;

%%%%%%%%%% separate the boundary and center nodes %%%%%%%%%%

idx_C = [];
bk = boundary(SG.Nodes.comx,SG.Nodes.comy);
bk_left = [bk(end-3:end-1);bk(1:13)];
% bk_left = bk(1:10);
bk_right = bk(83:97);
idx_B1 = bk_left';
idx_B2 = bk_right';
%%
nodeset = 1:1:numnodes(SG);
% nodeset = nodeset(nodeset ~= bk_left');
for nodenum = nodeset
        if ismember(nodenum,bk_left) || ismember(nodenum,bk_right)
            continue
%         elseif nodez_SG(nodenum).comy > boundary_limits(2)
%                idx_B2 = [idx_B2 nodenum];
        else       
            if degree(SG,nodenum) > 1
                idx_C = [idx_C nodenum];
            end
        end
end
G = subgraph(SG, [idx_B1 idx_B2 idx_C]);   

%%
G.Nodes.ID = (1:1:numnodes(G))';
G.Edges.ID = (1:1:numedges(G))';
limit1 = length(idx_B1);
limit2 = length(idx_B2);
Boundary2.left = G.Nodes.ID(1:limit1);
Boundary2.right = G.Nodes.ID(1+limit1:limit1+limit2);
G.Edges.Lengths = ptm.*(G.Edges.Weight);
G.Edges.Resistances = (16*8*mu/pi).*G.Edges.Lengths./(G.Edges.Widths).^4; %resistance = 8*mu*L/(pi*R^4)
% G.Edges.Resistances =12*mu.*G.Edges.Lengths./((G.Edges.Widths.^3).*1e-3.*(1-0.063.*G.Edges.Widths./1e-3));
%%
G = potSolver(G,V1,Boundary2,G);
[G.Edges.Flows,G.Edges.Shear] = FindEdgeFlows(G,mu);
G.Edges.Open = ones(numedges(G),1);