classdef particle
   properties
      edge_num {mustBeNumericOrLogical}
      node_num {mustBeNumericOrLogical}
      edge_x {mustBeNumericOrLogical}
      deposited {mustBeNumericOrLogical}
      edge_history {mustBeNumericOrLogical}
      destination_node {mustBeNumericOrLogical}
      comx(1,1) {mustBeNumeric, mustBeFinite}
      comy(1,1) {mustBeNumeric, mustBeFinite}
      out {mustBeNumericOrLogical}
   end
   methods
       % initialize the particle
       function obj = init_particle(obj, Nodez)
           Node_index = randi(size(Nodez,1));
           obj.node_num = Nodez.ID(Node_index);% randomly choose the node that particle begins moving forward from
            % set the particle position
           obj.comx =  Nodez.comx(Node_index);
           obj.comy =  Nodez.comy(Node_index);
           obj.edge_num = 0;            
           obj.edge_x = 0;
           obj.deposited= 0;
           obj.edge_history = 0;
           obj.destination_node = 0;
           obj.out = 0;
       end
       

       % choose the particle's next action
       % returns the update particle and graph, takes in particle, graph,
       % time-step, probbility of depositing, probability of erosion, and
       % particle diameter
       % edit: got rid of depth and varibales calculable from G and G2 

       function [objdt,G2,G1] = next_action(obj,deltaT,mu,er_c,dep_c,tau_er,tau_dep,p_diam,G2,G1)

           if (nargin~=10)
               error('Error:  Invalid number of input arguments.')
           end

               % check if the particle is following along an edge, if not choose
               % an edge for the particle to follow
               if obj.edge_num == 0
                    Inc = incidence(G2);
                    Edgez = G2.Edges;
                    Nodez = G2.Nodes;
                  objdt = choose_edge(obj,Inc,Edgez,Nodez);
              
               else
               % if the particle is already following along an edge, move the particle forward according
               % to the flow rate and time-step size  
                    Edgez = G2.Edges;
                    Nodez = G2.Nodes;
                  [objdt,G2,G1] =  move_p(obj,Edgez,deltaT,mu,er_c,dep_c,tau_er,tau_dep,Nodez,p_diam,G2,G1); 
               end
           
       end
       
      % choose the next edge at each junction
      function obj = choose_edge(obj,Inc,Edgez,Nodez)
         tlcf = mean(abs(Edgez.Flows))./1e5;
         % find the potential edges to follow
         Node_index = find(Nodez.ID == obj.node_num);
         e_set = find(Inc(Node_index,:));
%          end_node_list = zeros(50,1);
%          i = 1;
         for edg = e_set
             endpointz = Edgez.EndNodes(edg,:);
             end_node2 = endpointz(endpointz ~= Node_index);
                %only choose the edge if it has a lower potential
                if Nodez.Potentials(Node_index) < Nodez.Potentials(end_node2)
                    e_set = e_set(e_set~=edg);
%                     end_node_list(i)  = end_node2;
                end
%                 i = i + 1;
         end

         e_set = e_set(Edgez.Flows(e_set) > tlcf);
         if isempty(e_set) == 0
%              obj.node_num = Nodez.ID(end_node_list(randi(length(end_node_list))));
             % get the probability of choosing each edge based on its flow rate
             e_prob = abs(Edgez.Flows(e_set));
             %%normalize e_prob
               prob_norm=[0 e_prob']/sum(e_prob);
                %%create cumlative distribution
                p_dist=cumsum(prob_norm);
            %%calculate which bin the random number falls into (which edge the particle selects)
                [~,~,inds] = histcounts(rand,p_dist); 
                obj.edge_num = Edgez.ID(e_set(inds));
                % set the particle destination node
                endpointz = (Edgez.EndNodes(e_set(inds),:));
                obj.destination_node = Nodez.ID(endpointz(endpointz ~= Node_index));
                obj.comx =  Nodez.comx(Node_index);
                obj.comy =  Nodez.comy(Node_index);
                obj.node_num = 0;
         end
         
      end
      % move the particle along an edge, deposit w probability dep_prob
            function [obj,G2,G1] = move_p(obj,Edgez,deltaT,mu,er_c,dep_c,tau_er,tau_dep,Nodez,p_diam,G2, G1)
                
                if G1.Edges.Open(obj.edge_num) == 1 % check if the channel is open
                    Edge_index = find(Edgez.ID == obj.edge_num,1); % find the edge index corresponding to the edge ID
                    if ~isempty(Edge_index) % why would it be empty?
                      %move the particle if it is not deposited
                      if obj.deposited == 0
                          [obj,G2,G1] = traverse(obj, G1, G2,deltaT,Nodez);
                          if abs(obj.edge_x) >= (Edgez.Lengths(Edge_index))
                                obj.node_num = obj.destination_node;
                                dest_node_idx = find(Nodez.ID == obj.destination_node);
                                obj.comy = Nodez.comy(dest_node_idx);
                                obj.comx = Nodez.comx(dest_node_idx);
                                obj.edge_x = 0;
                                obj.edge_num = 0;
                          end
                        %depositing particles
                        [obj,G2,G1] = deposit(obj, G1, G2,dep_c,tau_dep,p_diam,deltaT,mu);
                      else
%                         % erode particles if deposited w/ some prob.
                        [obj,G2,G1] = erode(obj, G1, G2,er_c,tau_er,p_diam,deltaT,mu); 
                      end
                    end
                else
                    if obj.deposited == 1
                        [obj,G2,G1] = erode(obj, G1, G2,er_c,tau_er,p_diam,deltaT,mu);  
                    end
                end
                    
            end
                
%% erosion function
            function [obj,G2,G1] = erode(obj, G1, G2,er_c,tau_er, p_diam,deltaT,mu)
                if G1.Edges.Open(obj.edge_num) == 0 % check if the channel is clogged
                    Edge_index = obj.edge_num; %find(G1.Edges.ID == obj.edge_num); % check the reference graph for the clogged channel
                    tau =  G1.Edges.Shear(Edge_index); % use the shear of the edge in the reference graph
                    erosion_prob = er_c.*(tau - tau_er).*deltaT; % calculate the probability of erosion
                    if erosion_prob > rand
                    obj.deposited = 0;  % the particle erodes
                        G1.Edges.Open(Edge_index) = 1; % unclog the channel
                        G1.Edges.Widths(Edge_index) = G1.Edges.Widths(Edge_index)+ p_diam; % restore the channel capacity
                        G1.Edges.Resistances(Edge_index) =...
                            (16*8*mu/pi).*G1.Edges.Lengths(Edge_index)./(G1.Edges.Widths(Edge_index)).^4; %update the resistance
                        % restore the edge in the current graph
                        new_edge = G1.Edges(Edge_index,:); 
                        % correct the node number by using node ID's
                        end_node1 = new_edge.EndNodes(:,1); 
                        end_node1ID = G1.Nodes.ID(end_node1);
                        end_node2 = new_edge.EndNodes(:,2);
                        end_node2ID = G1.Nodes.ID(end_node2);
                        new_edge.EndNodes(:,1)  = find(G2.Nodes.ID == end_node1ID,1); % use the correct node index
                        new_edge.EndNodes(:,2) = find(G2.Nodes.ID == end_node2ID,1); % use the correct node index
                        G2 = addedge(G2,new_edge); % add the edge back into the network
                    end
                else % if the channel is open use the current graph
                    Edge_index = find(G2.Edges.ID == obj.edge_num); % use the current network for finding the edge index
                    tau = G2.Edges.Shear(Edge_index); % find the shear stress using the edge index
                    erosion_prob = er_c.*(tau - tau_er).*deltaT; % calculate the erosion probability
                    if erosion_prob > rand
                    obj.deposited = 0;  % the particle erodes
                        G2.Edges.Widths(Edge_index) = G2.Edges.Widths(Edge_index)+p_diam; % restore the channel capacity
                        G2.Edges.Resistances(Edge_index) = ...
                            (16*8*mu/pi).*G2.Edges.Lengths(Edge_index)./(G2.Edges.Widths(Edge_index)).^4; %update the channel resistance
                        Edge_index = obj.edge_num; %find(G1.Edges.ID == obj.edge_num); % find the edge index in reference graph
                        G1.Edges.Widths(Edge_index) = G1.Edges.Widths(Edge_index)+p_diam; % restore channel capacity
                        G1.Edges.Resistances(Edge_index) =...
                            (16*8*mu/pi).*G1.Edges.Lengths(Edge_index)./(G1.Edges.Widths(Edge_index)).^4; %update resistance  
                    end
                end
            end
%%  deposition function  
            function [obj,G2,G1] = deposit(obj, G1, G2,dep_c,tau_dep, p_diam, deltaT,mu)
                Edge_index = find(G2.Edges.ID == obj.edge_num); % find the index correspond to the edge ID
                tau = G2.Edges.Shear(Edge_index); % get the edge shear stress
                dep_prob = dep_c.*(tau_dep - tau).*deltaT; % calculate the deposition probability
                    if dep_prob > rand
                        obj.deposited = 1; %deposit the particle 
                        G2.Edges.Widths(Edge_index) = G2.Edges.Widths(Edge_index)-p_diam; %update the edge capacity
                        G2.Edges.Resistances(Edge_index) =... %update the edge resistance
                            (16*8*mu/pi).*G2.Edges.Lengths(Edge_index)./(G2.Edges.Widths(Edge_index)).^4;
                        % clog the edge if the capacity has reached
                        if G2.Edges.Widths(Edge_index) <= p_diam
                            G2.Edges.Open(Edge_index) = 0; % clog the edge in current graph
                            G2 = rmedge(G2, Edge_index); % remove the cogged edge
%                             Edge_index = find(G1.Edges.ID == obj.edge_num); % find the edge in the ref graph
%                             G1.Edges.Widths(Edge_index) = G1.Edges.Widths(Edge_index)-p_diam; %update the edge capacity
%                             G1.Edges.Resistances(Edge_index) =...
%                             (16*8*mu/pi).*G1.Edges.Lengths(Edge_index)./(G1.Edges.Widths(Edge_index)).^4; %update resistance  
                            G1.Edges.Open(obj.edge_num) = 0;
                        end
                        % update the edge in ref graph
                        Edge_index = obj.edge_num; %find(G1.Edges.ID == obj.edge_num);
                        G1.Edges.Widths(Edge_index) = G1.Edges.Widths(Edge_index)-p_diam; %update the edge capacity
                        G1.Edges.Resistances(Edge_index) =...
                            (16*8*mu/pi).*G1.Edges.Lengths(Edge_index)./(G1.Edges.Widths(Edge_index)).^4;%update resistance  
                    end
            end

%% traverse function
            function [obj,G2,G1] = traverse(obj, G1, G2,deltaT,Nodez)
                Edge_index = find(G2.Edges.ID == obj.edge_num);
                Edge_i = G2.Edges(Edge_index,:);
                crossA = (Edge_i.Widths).^2;
                velo = (Edge_i.Flows)./crossA;
                obj.edge_x  = obj.edge_x + deltaT.*velo ;
                endpointz = (Edge_i.EndNodes(:));
                dest_node_idx = find(Nodez.ID == obj.destination_node);
                node_from = endpointz(endpointz ~= dest_node_idx);
                velx = abs(velo).*(Nodez.comx(dest_node_idx) - Nodez.comx(node_from))/(Edge_i.Lengths);
                vely = abs(velo).*(Nodez.comy(dest_node_idx) - Nodez.comy(node_from))/(Edge_i.Lengths);
                obj.comy = obj.comy + deltaT.*vely;
                obj.comx = obj.comx + deltaT.*velx;
            end
            
   end
end