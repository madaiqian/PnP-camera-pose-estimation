clear; clc;
addpath(genpath('rpnp1.0'));
addpath(genpath('C:\zheng\work0919\PnP+f\UPnP'))

warning off;

% experimental parameters
ntest= 500;
nl= 2;
npts= 6:100:1006;

name= { 'UPnPf', 'UPnPf+GN', 'GPnPf',   'GPnPf+GN'};
f= { @upnp_interface,  @upnp_GN_interface, @GPnP_f, @GPnP_f_GN};
marker= { 'o', 'd', '>', 's'};
color= {'g','b','m','k'};
markerfacecolor=  {'g','n','m','n'};
linestyle= {'-','-','-','-'};

% name= {'EPnP+GN', 'RPnP', 'DLS',    'OPnP'};
% f= { @EPnP_GN,  @RPnP, @dls_pnp,@OPnP};
% marker= { 'o', 'd', '^',  's'};
% color= {'g','b','m','k'};
% markerfacecolor=  {'g','n','m','n'};

method_list= struct('name', name, 'f', f, 't', zeros(size(npts)),...
    'marker', marker, 'color', color, 'markerfacecolor', markerfacecolor);


% camera's parameters
width= 640;
height= 480;
f= 800;

for j= 1:length(npts)
    npt= npts(j);
    fprintf('%d points\n',npt);
    
    % generate experimental data
    para= cell(ntest,2);
    for i= 1:ntest
        % generate 3d coordinates in camera space
        Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
        t= mean(Xc,2);
        R= rodrigues(randn(3,1));
        XXw= inv(R)*(Xc-repmat(t,1,npt));
        
        % projection
        xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
        xxn= xx+randn(2,npt)*nl;
		xxn= xxn/f;
        % save
        para{i,1}= XXw;
        para{i,2}= xxn;
    end

    for k= 1:length(method_list)
        tic;
        for i= 1:ntest
            XXw= para{i,1};
            xxn= para{i,2};
            method_list(k).f(XXw,xxn);
        end
        t= toc; method_list(k).t(j)= t/ntest*1000;
        disp([method_list(k).name ' - ' num2str(t) ' s']);
    end


end

close all;

figure('color','w','position',[100 100 360 320]);
hold all;
box on;
p= zeros(size(method_list));
for k= 1:length(method_list)
    p(k)= plot(npts,method_list(k).t,'color',method_list(k).color,...
        'marker',method_list(k).marker,...
        'markerfacecolor',method_list(k).markerfacecolor,...
        'displayname',method_list(k).name, ...
        'LineWidth',2);
end
legend(p,2);
xlim(npts([1,end]));

xtick= 6:200:1006;
set(gca,'xtick',xtick);

xlabel('Number of Points','FontSize',11);
ylabel('Computational Time (milliseconds)','FontSize',11);

rmpath(genpath('C:\zheng\work0919\PnP+f\UPnP'))
rmpath(genpath('rpnp1.0'));