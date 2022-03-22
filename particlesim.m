%% particlesim.m 
% main file - shortcut version -- line numbers 7 and 17 get modified

mu = 60e-3;p_dmeter = 1e-6;ptm = (1e-6)/1.61; % dynamic vis, particle diameter, pixels to meters
Nstep = 800;tot_time = 50;dt = tot_time/Nstep; % dt has to be < 0.5

V1 = 241000; % high pressure
Vref = 1.6e+05; % pressure in imported graph G
G = load('G.mat').G; %load pre-solved network
G.Edges.Flows = G.Edges.Flows.*V1./Vref; %adjust the flow
G.Edges.Shear = G.Edges.Shear.*V1./Vref; % adjust the shear stress
G.Nodes.Potentials = G.Nodes.Potentials.*V1./Vref; %adjust the pressure at nodes 
G2 = G; % working network
Boundary2.left = G.Nodes.ID(1:16);Boundary2.right = G.Nodes.ID(1+16:16+15);
beginning_set = Boundary2.left; %inlet nodes
target1 = Boundary2.right;%outlet nodes
tau_d = 180; %shear stress threshold for deposition
tau_d = tau_d/10;
tau_er = tau_d; %shear stress threshold for erosion
dep_c = 50/max(G.Edges.Shear); %dep_probability = dep_c.*(tau_d - G.Edges.Shear).*dt; %max is less than 1 for max(tau_d)=2200
er_c = 8/max(G.Edges.Shear); %er_probability = er_c.*(G.Edges.Shear-tau_er).*dt; %max is less than 1 when tau_er = 0

% local_pot0 = abs(G.Edges.Flows).*(G.Edges.Resistances); %local pressure
% across an edge

%% 
tlcf = mean(abs(G.Edges.Flows))./(1e5); %break if the flow falls below this threshold

t_transit = max(G.Nodes.comx)./(mean(G.Edges.Flows)./(mean(G.Edges.Widths).^2)); %transit time for particles
passed_time = 0; % passed time in simulation
inj_step = 2*floor(t_transit/dt); %particle injection step

Nparticle = 230; %number of particles injected at each step
% inedges = find(Inc(Boundary2.left,:));
% avg_capacity = sum(floor(G.Edges.Widths(inedges,:)./p_dmeter),'all');
% %running at about 85% capacity
% inj_rate = Nparticle/inj_step;

tot_Nparticle = (floor(Nstep./inj_step)+1)*Nparticle; % total # of particles
if inj_step > Nstep
    tot_Nparticle = Nparticle; 
end
Nparticle_active = Nparticle; %current # of injected and active particles
particle_set = particle.empty(tot_Nparticle,0); %initialized particle set

%initilize the particles (randomly assign them at inlet nodes)
for p_idx = 1:1:tot_Nparticle
    particle_set(p_idx) = particle();
    particle_set(p_idx) = init_particle(particle_set(p_idx),G.Nodes(beginning_set,:));
end

%%
p_out = 0; % zero particles have exited the medium 
k1 = 1; % counter
for step= 1:1:Nstep
    GID = G2.Edges.ID; %current node IDs of G2
    
    if p_out > tot_Nparticle-2 %if all particles exited, stop
        break
    end

    % inject particles at injection step
    if passed_time > 0
        if mod(step, inj_step) == 0 
            Nparticle_active = Nparticle_active + Nparticle;
        end
    end

    %swap entries for particles that exited
    %choose next action for other active particles
  for p_idx = 1:1:Nparticle_active
%       particle_set(p_idx) = checkOut(particle_set(p_idx),target1);
     [Nparticle_active,p_out,particle_set] = pSwap(particle_set, p_idx, Nparticle_active, tot_Nparticle,p_out,target1);
     [particle_set(p_idx),G2,G] = next_action(particle_set(p_idx),...
         dt,mu,er_c,dep_c,tau_er,tau_d,p_dmeter,G2,G); 
  end
  %update edge resistances
    G2.Edges.Resistances = (16*8*mu/pi).*G2.Edges.Lengths./(G2.Edges.Widths).^4;
   %redo edge flow calculations if edges are clogged or unclogged
    if size(unique(G2.Edges.ID),1) ~= size(unique(GID),1) | unique(G2.Edges.ID) ~= unique(GID)
        G2 = potSolver(G2,V1,Boundary2,G);
        [G2.Edges.Flows,G2.Edges.Shear] = FindEdgeFlows(G2,mu);
    end
    %save the graph and the particle set 
    if mod(step,666) == 0  %666
         step
        save([pwd strcat('/p',num2str(V1/1000),'k',...
num2str(tau_er),num2str(k1),'t2.mat')],'particle_set');
        save([pwd strcat('/G',num2str(V1/1000),'k',...
 num2str(tau_er),num2str(k1),'t2.mat')],'G2');
        
        k1 = k1 +1;
    end
%stop if the mean flow rate is below the threshold 
    if abs(mean(G2.Edges.Flows)) < tlcf
        'broken'
        break
    end
    passed_time = dt*step; %update the current time
end  
%%
% p = Gplot(G2,ptm,G2.Edges.Shear,winter,particle_set, 1); %plot the
% particles
% set(gca,  'CLim', [0 4]); %set the limits on the colorbar
%%

function p = Gplot(G,ptm,edgecdata,cmp,particle_set, pplot)

figure(1);

p = plot(G);
p.XData = (G.Nodes.comx)./ptm;
p.YData =(G.Nodes.comy)./ptm;
p.EdgeCData = edgecdata; 
p.Marker = 'none';
colormap(cmp)
p.LineWidth = (G.Edges.Widths).*5e5;
colorbar
if pplot == 1
    Nparticle_active = length(particle_set);
    hold on;
        for iii = 1:1:Nparticle_active
            if particle_set(1,iii).deposited == 1
                p_color = 'r';
                e_color = 'r';
            else
                p_color = 'w';
                e_color = 'k';
            end
            if particle_set(1,iii).out == 0
                plot((particle_set(1,iii).comx)./ptm,(particle_set(1,iii).comy)./ptm,'o',...
                    'lineWidth',1, 'MarkerFaceColor', p_color,'MarkerEdgeColor',e_color,...
                    'markerSize', 3.2);hold on;
            end
        end
end



end


function [edgeQ,shear] = FindEdgeFlows(G,mu)
nodez_HSG = G.Nodes;
Egraph = G.Edges;
El = numedges(G);
edgeQ = ones(El,1);
shear = ones(El,1);
for ed =1:1:El
    twonodes = (Egraph.EndNodes(ed,:));
    twonodes_x = nodez_HSG(twonodes,:).comx;
    [~, max_x_idx] = max(twonodes_x);
%     [~, min_x_idx] = min(twonodes_x);
    edq  = (nodez_HSG.Potentials(twonodes([1,2] ~= max_x_idx)) - nodez_HSG.Potentials(twonodes(max_x_idx)))./Egraph.Resistances(ed);
    edgeQ(ed) = edq; % dividing the thickness of the system (1e-3)...
    % by characteristic pore size (1e-5) results in ~100
    shear(ed) = abs(edq*32*mu/(pi*(Egraph.Widths(ed,1))^3)); % tau_w = 4*mu*q/pi r^3
    %     lambda = 24./((1-0.351.*edge_thic(ed)./depth).*(1+edge_thic(ed)./depth)).^2;
%     shear_rate2(ed) = abs(ed_col.^2.*res_vec(ed).*lambda./(8*(edge_thic(ed).*depth)^2));
    

end
end

function [Graph00, Nodes00, Links00, disto] = Im2Graph(Im1, div)
    % first find the gray threshold then binarize the image 
    ThR=graythresh(Im1)/div;
    BW = imbinarize(Im1,ThR);
    % create the image skeleton from the binary image
    C = bwskel(BW(:,:,1));
    % calculate the distance to the nearest bead from each pixel
    disto = bwdist(imcomplement(BW(:,:,1)));
    % use the package 'Skel2Graph' to calculate the graph corresponding to
    % the skeleton
    [Graph00,Nodes00, Links00] = Skel2Graph3D(C,100);
end

function [edge_thic] = FindEdgeProps(G,dists)
El = numedges(G);
edge_thic = ones(El,1);
nodez = G.Nodes;
% edge_leng = ones(El,1);
% cntr_xs = ones(El,1);
for ed = 1:1:El
    
    twonodes = (G.Edges.EndNodes(ed,:));
    cntr_y = round((nodez.comx(twonodes(1)) + nodez.comx(twonodes(2)))/(2));%.*ptm .*ptm
    cntr_x = round((nodez.comy(twonodes(1)) + nodez.comy(twonodes(2)))/(2));
    edge_thic(ed) = 2.*dists(cntr_x,cntr_y); %*log(2)
%     if edge_thic(ed) == 0
%         edge_thic(ed) = 11.967*ptm;  %%%% 'eliminating' these edges
%     end
end
end

function G00 = potSolver(G00,V00,bdry_lm,G01)

[bin,binsize] = conncomp(G00);%,'Type','weak');
idx = binsize(bin) == max(binsize);
SG00 = subgraph(G00, idx);

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

G00 = subgraph(SG00, [idx_B100; idx_B200; idx_C00]);   
limit1 = length(idx_B100);
limit2 = length(idx_B200);

I = incidence(G00);

% D = (abs(I)*ws).^(-1/2);
% D = diag(D);
W = double(diag(1./G00.Edges.Resistances));
L = I*W*I';

lb = limit1+limit2;
lL = length(L);
LBB = L(1:lb,1:lb);
LBC = L(1:lb,lb+1:lL);
LCB = L(lb+1:lL,1:lb);
LCC = L(lb+1:lL,lb+1:lL);

psi_B100 = V00.*ones(limit1,1); %%% set up the voltage vector for boundaries
psi_B200 = zeros(limit2,1);

psi_B00 = [psi_B100; psi_B200];
LCCinv = inv(LCC);
LS = LBB - (LBC*(LCCinv*LCB));
JB = LS*psi_B00;

J = zeros(lL,1);
J(1:lb) = JB;
pot_vec = L\J;
%psi_C = -LCCinv*(LCB*psi_B); %%the solution for the special case of two b%nodes
CV1 = -pot_vec(1) + V00;
% testtest=CV1.*ones(lL,1);
pot_vec = CV1+pot_vec; %% getting the unique solution by a shift

G00.Nodes.Potentials = pot_vec;
    for iii = 1:1:numnodes(G01)
        if ismember(G01.Nodes.ID(iii),G00.Nodes.ID)
            continue
        else
            new_node  = G01.Nodes(iii,:);
            new_node.Potentials = 0;
            G00 = addnode(G00,new_node);
        end
    end


end

% swap out the spots for particles that have exited the medium
function [Nparticle_active,p_out,particle_set] = pSwap(particle_set, p_idx, Nparticle_active, tot_Nparticle,p_out,targetnodes)
    if ismember(particle_set(p_idx).node_num,targetnodes) % check if the particle has reached the end        
            particle_set(p_idx).out = 1;
            particle_set([p_idx tot_Nparticle-p_out]) = particle_set([tot_Nparticle-p_out p_idx]); % swap the out particle with the last particle that is not out
            particle_set([p_idx Nparticle_active]) = particle_set([Nparticle_active p_idx]); % swap the new particle with the first element of next set of injected particles
            Nparticle_active = Nparticle_active - 1; % update the number of active particles
            p_out = p_out + 1; % update the number of out particles
    end

end

% check if the particle has exited the medium
%        function obj = checkOut(obj,targetnodes)
%            if ismember(obj.node_num,targetnodes)
%                obj.out = 1;
%            end
%        end
