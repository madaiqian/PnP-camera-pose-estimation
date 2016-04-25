clear; clc;
addpath epnp;
addpath lhm;
addpath func;
addpath sp;

% experimental parameters
nl= 3;
npts= 4:1:15;
num= 2000;

% compared methods
A= zeros(size(npts));
B= zeros(num,1);
% name= {'Best', '4VB', 'RPnP', 'ASPnP_un',          'ASPnP'};
% f= {    @robust_dls_pnp_my2, @PnP_3Order_4Variable_Floating_Grobner_Balance,  @RPnP, @ASPnP_unnormalized, @ASPnP};
% marker= { 'x', 'o', 'd', '^','s'};
% color= {'r','g','b','m','k'};
% markerfacecolor=  {'r','n','n','m','k'};
name= {'RPnP', 'ASPnP'};
f= {   @RPnP, @ASPnP};
marker= { 'x', 'o'};
color= {'r','g'};
markerfacecolor=  {'r','n'};

method_list= struct('name', name, 'f', f, 'mean_r', A, 'mean_t', A,...
    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, 'r', B, 't', B,...
    'marker', marker, 'color', color, 'markerfacecolor', markerfacecolor);

warning off all
% experiments
for i= 1:length(npts)
    
    npt= npts(i);
    fprintf('npt = %d: ',npt);
    
    for k= 1:length(method_list)
        method_list(k).r = zeros(1,num);
        method_list(k).t = zeros(1,num);
    end
    
    index_fail = [];
    for j= 1:num
        
        % camera's parameters
        width= 640;
        height= 480;
        f= 800;
        
        % generate 3d coordinates in camera space
        XXw= [xrand(2,npt,[-2 2]); zeros(1,npt)];
		R= rodrigues(randn(3,1));
        t= [rand-0.5;rand-0.5;rand*8+4];
        Xc= R*XXw+repmat(t,1,npt);
        
        % projection
        xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
        xxn= xx+randn(2,npt)*nl;

        % pose estimation
        for k= 1:length(method_list)
            [R1,t1]= method_list(k).f(XXw,xxn/f);
            if size(t1,2) ~= 1
                index_fail = [index_fail, j];
                break;
            end
            y= cal_pose_err([R1 t1],[R t]);
            method_list(k).r(j)= y(1);
            method_list(k).t(j)= y(2);
        end

        showpercent(j,num);
    end
    
    fprintf('\n');
	
    % save result
    for k= 1:length(method_list)
        method_list(k).r(index_fail) = [];
        method_list(k).t(index_fail) = [];
        
        method_list(k).mean_r(i)= mean(method_list(k).r);
        method_list(k).mean_t(i)= mean(method_list(k).t);
        method_list(k).med_r(i)= median(method_list(k).r);
        method_list(k).med_t(i)= median(method_list(k).t);
        method_list(k).std_r(i)= std(method_list(k).r);
        method_list(k).std_t(i)= std(method_list(k).t);
    end
    
end

close all;
yrange= [0 20];

i= 0; w= 300; h= 300;

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'mean_r','Mean Rotation Error',...
    'Number of Points','Rotation Error (degrees)');

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'med_r','Median Rotation Error',...
    'Number of Points','Rotation Error (degrees)');

yrange= [0 30];
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'mean_t','Mean Translation Error',...
    'Number of Points','Translation Error (%)');

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'med_t','Median Translation Error',...
    'Number of Points','Translation Error (%)');
