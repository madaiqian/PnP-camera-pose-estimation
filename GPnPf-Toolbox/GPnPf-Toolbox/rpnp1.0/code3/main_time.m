clear; clc;
addpath epnp;
addpath lhm;
addpath func;

% experimental parameters
ntest= 1000;
nl= 5;
npts= 4:50:1000;

name= {'EPnP','RPnP','DLS','ASPnP'};
f= {@EPnP, @RPnP,@dls_pnp,@ASPnP};
marker= { 'o', '+', '^','s'};
color= {'b','c','m','k'};
markerfacecolor= {'n','r','m','k'};

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
        % generate 3d coordinates in target space
        XYZrng= [-0.5 0.5]*height;
        XXw= xrand(3,npt,XYZrng);
        % generate rotation matrix and translation vector
        R= rodrigues(randn(3,1));
        t= [(rand-0.5)*width; (rand-0.5)*height; f]*rand*2 + [0; 0; f];
        % project to image points
        XX= R*XXw+repmat(t,1,size(XXw,2));
        xx= [XX(1,:)./XX(3,:); XX(2,:)./XX(3,:)]*f;
        xxn= xx+randn(2,npt)*nl;
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
        'displayname',method_list(k).name);
end
legend(p,2);
xlim(npts([1,end]));

xtick= 4:50:1000;
set(gca,'xtick',xtick);

xlabel('Number of Points');
ylabel('Computational Time (milliseconds)');

