%% expFx
V1 = 260;

if V1 == 260
    pstr = '.tif'; Im = imread('Beads.tif'); Th=graythresh(Im)/1.5; BW = imbinarize(Im,Th);
    N = 50;
%     x_vec = ptm.*(1:8790)';
    x_vec = (1:8790)'./8790;
else
    N = 30;
    pstr = 'localized.tif'; Im = imread('Beadslocalized2.tif'); Th=graythresh(Im)/1.5; BW = imbinarize(Im,Th);
%     x_vec = ptm.*(1:9928)';
    x_vec = (1:9928)'./9928;
end

% p00 = double(BW);
Fx_vec = zeros(size(Im, 2),N);

for k1 = 1:1:N
 
    particles1 = imread(strcat('Particles',num2str((k1)),pstr));
    BWp  = imbinarize(particles1,graythresh(particles1)/3);   
    p11 = double(particles1);
 
    dep = sum(p11, 'all'); 
    temp_vec = sum(p11, 1);
    Fx_vec(:,k1) = cumsum(temp_vec, 2)./dep; 
end

%% plot
% nexttile
timevec = 1:1:N;
timevec = timevec./N;
ptm = (1e-6)/1.61;

for k1 = 1:1:N
    plot(x_vec, Fx_vec(:,k1), 'LineWidth', 3,'color', [(1-timevec(k1)) (1-timevec(k1)) 0.5]);hold on; 
%     plot(x_vec, Fx_vec(:,k1), 'LineWidth', 3,'color', [1 0.8*(1-timevec(k1)) 0.8*(1-timevec(k1))]);hold on; 
end
text1 = strcat('$',num2str(V1),'kPa$');
% title(strcat(num2str(V1),'kPa'), 'fontsize', 16, 'fontname', 'Times New Roman')
set(gca, 'fontsize', 18, 'fontname', 'Times New Roman');
xlabel('$x$', 'interpreter', 'latex')
ylabel('$F(x)$', 'interpreter', 'latex')
annotation('textbox',...
    [0.71 0.6 0.322 0.081],...
    'String',text1,...
    'Interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');
ylim([0 1])


%%

Im = imread('Beads.tif'); Im = Im(:,1903:end,:); %6887

% psymap1 = [255 255 255
% 55 145 230]./255;
%%
% nexttile

N = 6;  %60 for t1 204,50,20-6 for 50t2  - 20 for t2,204,20 -
k3 = 5;k1 = N-k3;
dep_l = ones(1,N);
taucr = 20;
trial = 3;
Pressure =56; 
% tilez = tiledlayout(1,2);
% for Pressure = P_set
ptm = (1e-6)/1.61;
timevec = 1:N;
timevec = timevec./N;
timevec = 1:k3;
timevec = timevec./k3;
%%%
%% SIM FX PLOTS
PX = zeros(1610,N);l=1;

for k2 = 1:k3
    

if trial == 0
    particle_set3 = load(strcat('p',num2str(Pressure),'k',num2str(taucr),num2str(k1),'.mat'));
else
    particle_set3 = load(strcat('p',num2str(Pressure),'k',num2str(taucr),num2str(k1),'t',num2str(trial),'.mat'));
end  
particle_set = particle_set3.particle_set;

%     G2 = load(strcat('G',Pressure,num2str(k1),'.mat')).G2;
    Nparticle_active = length(particle_set3.particle_set);
%% uncomment for FX plots

Np_active = Nparticle_active;
% Np_active = length(particle_set);
PX240k1 = ones(Np_active,1);
iiib = 1;
Edges_set = [];
    for iii = 1:1:Np_active   
            if particle_set(1,iii).out == 0 %particle_set3.particle_set(1,iii).out == 0
                if particle_set(1,iii).deposited == 1 %particle_set3.particle_set(1,iii).deposited == 1
                    PX240k1(iiib,1) = particle_set(1,iii).comx;
                    Edges_set(iiib,1) = particle_set(1,iii).edge_num;
                    Edges_set(iiib,2) = iii;
%                     PY310k1(iii,1) = particle_set3.particle_set(1,iii).comy;
                    iiib = iiib+1;
                end
            end
    end
%%
% ID_in = G1.Edges.ID;
% ID_fin = G2_f.Edges.ID;
% i = 1;
% for ii = 1:length(Edges_set)
%     if isempty(find(ID_fin(:,1) == Edges_set(i,1),1)) == 0
%         ID_cmpi(i) = find(ID_in(:,1) == Edges_set(i,1),1);
%         ID_cmpf(i) = find(ID_fin(:,1) == Edges_set(i,1),1);
%         flow_in(i)  = G1.Edges.Flows(ID_cmpi(i));
%         flow_fin(i)  = G2_f.Edges.Flows(ID_cmpf(i));
%         flow_cmp(i) = flow_in(i) - flow_fin(i);
%         i =i +1;
%     end
% end

%%
    PX240k1 = sort(PX240k1(PX240k1<1));
    Fx_n = PX240k1./(6888*ptm);
    Fx_n1 = Fx_n(Fx_n > 0.012);
    
    
%     PX(:,k1) = PX240k1; 
    k1 = k1 - N + k3 +1;
    colors = [[1 0.8*(1-timevec(k1)) 0.8*(1-timevec(k1))];...
    [(1-0.8*timevec(k1)) (1-timevec(k1)) 0.5];...
    [(1-0.8*timevec(k1)) 0.2*(1-timevec(k1)) 0.7]];
    COL = colors(3,:);
    dep_l(k1) =  length(Fx_n);
    dep_l1(k1) =  length(Fx_n1);
%     plot(Fx_n1,(1:dep_l1(k1))./dep_l(k1),'LineWidth', 3,'LineWidth', 3,'color', [(1-0.8*timevec(k1)) (1-timevec(k1)) 0.5]);hold on;
    plot(Fx_n,(1:dep_l(k1))./dep_l(k1),'LineWidth', 3,'LineWidth', 3,'color', COL);hold on;
    k1 = k1 + N - k3 ;
end
text2 = strcat('$',num2str(Pressure),'kPa, ', '\tau = ', num2str(taucr),'Pa $');
% title(strcat(num2str(Pressure),'kPa'), 'fontsize', 16, 'fontname', 'Times New Roman')
set(gca, 'fontsize', 18, 'fontname', 'Times New Roman');
xlabel('$x$', 'interpreter', 'latex')
ylabel('$F(x)$', 'interpreter', 'latex')
ylim([0 1])
% xlim([0 8790*ptm])
xlim([0 1])
% annotation('textbox',...
%     [0.71 0.136 0.322 0.081],...
%     'String',text2,...
%     'Interpreter','latex',...
%     'FontWeight','bold',...
%     'FontSize',18,...
%     'FontName','Times New Roman',...
%     'FitBoxToText','off',...
%     'EdgeColor','none');

% jj = jj+1;
% end
%%
% plotparticles(particle_set,Im,ptm)
%%
function plotparticles(particle_set,Im2,ptm)
h= imshow(Im2);hold on;set(h, 'AlphaData', 0.25);
Nparticle_active = length(particle_set);
        for iii = 1:1:Nparticle_active
            if particle_set(1,iii).deposited == 1
                p_color = 'r';
                e_color = 'r';
            else
                p_color = 'w';
                e_color = 'k';
            end
            if particle_set(1,iii).out == 0
                plot((particle_set(1,iii).comx)./ptm,(particle_set(1,iii).comy)./ptm,'o', 'lineWidth',1, 'MarkerFaceColor', p_color,'MarkerEdgeColor',e_color, 'markerSize', 3.2);hold on;
            end
        end


end
 