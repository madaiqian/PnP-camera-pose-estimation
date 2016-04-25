clear; clc;
addpath(genpath('rpnp1.0'));
addpath(genpath('DLT'));
addpath(genpath('C:\zheng\work0919\PnP+f\UPnP'))

warning off;
% experimental parameters
nl= 2;
npts= 4:1:15;
num= 500;

% compared methods
A= zeros(size(npts));
B= zeros(num,size(npts));
name= {'DLT', 'UPnPf', 'UPnPf+GN', 'GPnPf',   'GPnPf+GN', 'RPnP'};
f= {   @DLT, @upnp_interface,  @upnp_GN_interface, @GPnP_f, @GPnP_f_GN, @RPnP_interface};
marker= { 'x', 'o', 'd', '>', 's', '+'};
color= {'r','g','b','m','k','c'};
markerfacecolor=  {'r','g','n','m','n','n'};
linestyle= {'-','-','-','-','-','-'};

method_list= struct('name', name, 'f', f, 'mean_r', A, 'mean_t', A, 'mean_foc', A, 'mean_reproj', A, ...
    'med_r', A, 'med_t', A, 'med_foc', A, 'med_reproj', A, 'r', B, 't', B,...
    'foc', B, 'reproj', B, 'marker', marker, 'color', color, ...
    'markerfacecolor', markerfacecolor, 'linestyle', linestyle);

% experiments
for i= 1:length(npts)
    
    npt= npts(i);
    fprintf('npt = %d: ',npt);
    
  
    index_fail = [];
    
    for j= 1:num
        
        % camera's parameters
        width= 640;
        height= 480;
        f= rand(1)*1800+200; %random focal length in [200,2000]
        
        % generate 3d coordinates in camera space
        Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[1 2]); xrand(1,npt,[4 8])];
        t= mean(Xc,2);
        R= rodrigues(randn(3,1));
        XXw= inv(R)*(Xc-repmat(t,1,npt));
        
        % projection
        xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
        xxn= xx+randn(2,npt)*nl;

        if npt >= 6
            % pose estimation
            for k= 1:length(method_list)
                 if strcmp(method_list(k).name, 'RPnP') %using ground-truth focal length
                     [f1,R1,t1]= method_list(k).f(XXw,xxn,diag([f,f,1]));
                 else
                    try
                        [f1,R1,t1]= method_list(k).f(XXw,xxn);
                    catch
                        fprintf(['   The solver - ',method_list(k).name,' - encounters internal errors! \n']);
                        index_fail = [index_fail, j];
                        break;
                    end
                 end

                %no solution
                while size(t1,2) < 1
                    [f1,R1,t1]= method_list(k).f(XXw,xxn);                    
                end
                
                %choose the solution with smallest error 
                index_best = 1;
                error = inf;
                for jjj = 1:size(R1,3)
                    tempy = cal_pose_err([R1(:,:,jjj) t1(:,jjj)],[R t]);
                    if sum(tempy) < error
                        y = tempy;
                        error = sum(tempy);
                        index_best = jjj;
                    end
                end

                method_list(k).r(j,i)= y(1);
                method_list(k).t(j,i)= y(2);
                method_list(k).foc(j,i)= abs(f1(index_best)-f)/f*100;
                reproj = R1(:,:,index_best)*XXw+t1(:,index_best)*ones(1,npt);
                reproj = reproj./repmat(reproj(3,:),3,1);
                err = xxn-f1(index_best)*reproj(1:2,:);
                err = sqrt(sum(sum(err.*err))/npt);
                method_list(k).reproj(j,i) = err;
            end
        else
            for k= 1:length(method_list)
                 if strcmp(method_list(k).name, 'DLT') || strcmp(method_list(k).name, 'UPnPf') || strcmp(method_list(k).name, 'UPnPf+GN') 
                     method_list(k).r(j,i)= nan;
                     method_list(k).t(j,i)= nan;
                     method_list(k).foc(j,i) = nan;
                     method_list(k).reproj(j,i) = nan;
                 elseif strcmp(method_list(k).name, 'RPnP') 
                     [f1,R1,t1]= method_list(k).f(XXw,xxn,diag([f,f,1]));
                     
                    %choose the solution with smallest error 
                    index_best = 1;
                    error = inf;
                    for jjj = 1:size(R1,3)
                        tempy = cal_pose_err([R1(:,:,jjj) t1(:,jjj)],[R t]);
                        if sum(tempy) < error
                            y = tempy;
                            error = sum(tempy);
                            index_best = jjj;
                        end
                    end

                    method_list(k).r(j,i)= y(1);
                    method_list(k).t(j,i)= y(2);
                    method_list(k).foc(j,i)= abs(f1(index_best)-f)/f*100;
                    reproj = R1(:,:,index_best)*XXw+t1(:,index_best)*ones(1,npt);
                    reproj = reproj./repmat(reproj(3,:),3,1);
                    err = xxn-f1(index_best)*reproj(1:2,:);
                    err = sqrt(sum(sum(err.*err))/npt);
                    method_list(k).reproj(j,i) = err;
                    
                 else
                    try
                        [f1,R1,t1]= method_list(k).f(XXw,xxn);
                    catch
                        fprintf(['   The solver - ',method_list(k).name,' - encounters internal errors! \n']);
                        index_fail = [index_fail, j];
                        break;
                    end
                    %no solution
                    while size(t1,2) < 1
                        [f1,R1,t1]= method_list(k).f(XXw,xxn);
                    end
                    
                    %choose the solution with smallest error 
                    index_best = 1;
                    error = inf;
                    for jjj = 1:size(R1,3)
                        tempy = cal_pose_err([R1(:,:,jjj) t1(:,jjj)],[R t]);
                        if sum(tempy) < error
                            y = tempy;
                            error = sum(tempy);
                            index_best = jjj;
                        end
                    end

                    method_list(k).r(j,i)= y(1);
                    method_list(k).t(j,i)= y(2);
                    method_list(k).foc(j,i)= abs(f1(index_best)-f)/f*100;
                    reproj = R1(:,:,index_best)*XXw+t1(:,index_best)*ones(1,npt);
                    reproj = reproj./repmat(reproj(3,:),3,1);
                    err = xxn-f1(index_best)*reproj(1:2,:);
                    err = sqrt(sum(sum(err.*err))/npt);
                    method_list(k).reproj(j,i) = err;
                 end
            end
        end
    
        showpercent(j,num);
    end
    fprintf('\n');
    
    % save result
    for k= 1:length(method_list)
       
        method_list(k).mean_r(i)= mean(method_list(k).r(:,i));
        method_list(k).mean_t(i)= mean(method_list(k).t(:,i));
        method_list(k).mean_foc(i)= mean(method_list(k).foc(:,i));
        method_list(k).mean_reproj(i)= mean(method_list(k).reproj(:,i));
        
        method_list(k).med_r(i)= median(method_list(k).r(:,i));
        method_list(k).med_t(i)= median(method_list(k).t(:,i));
        method_list(k).med_foc(i)= median(method_list(k).foc(:,i));
        method_list(k).med_reproj(i)= median(method_list(k).reproj(:,i));

    end
end

%draw the boxplot of rotation error
close all;
yrange = [0 10];
i= 0; w= 300; h= 300;

figure('color','w','position',[w*i,100,w,h]);i=i+1;
boxplot(method_list(2).r,npts);
ylim(yrange); set(gca,'xtick',npts);
title('UPnPf','FontSize',12,'FontName','Arial');
xlabel('Number of Points','FontSize',11);
ylabel('Rotation Error (degrees)','FontSize',11);

figure('color','w','position',[w*i,100,w,h]);i=i+1;
boxplot(method_list(3).r,npts);
ylim(yrange); set(gca,'xtick',npts);
title('UPnPf+GN','FontSize',12,'FontName','Arial');
xlabel('Number of Points','FontSize',11);
ylabel('Rotation Error (degrees)','FontSize',11);

figure('color','w','position',[w*i,100,w,h]);i=i+1;
boxplot(method_list(4).r,npts);
ylim(yrange); set(gca,'xtick',npts);
title('GPnPf','FontSize',12,'FontName','Arial');
xlabel('Number of Points','FontSize',11);
ylabel('Rotation Error (degrees)','FontSize',11);

figure('color','w','position',[w*i,100,w,h]);i=i+1;
boxplot(method_list(5).r,npts);
ylim(yrange); set(gca,'xtick',npts);
title('GPnPf+GN','FontSize',12,'FontName','Arial');
xlabel('Number of Points','FontSize',11);
ylabel('Rotation Error (degrees)','FontSize',11);

rmpath(genpath('rpnp1.0'));
rmpath(genpath('DLT'));
rmpath(genpath('C:\zheng\work0919\PnP+f\UPnP'))

% close all;
% yrange= [0 15];
% 
% i= 0; w= 300; h= 300;
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'mean_r','Mean Rotation Error',...
%     'Number of Points','Rotation Error (degrees)');
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'med_r','Median Rotation Error',...
%     'Number of Points','Rotation Error (degrees)');
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'mean_t','Mean Translation Error',...
%     'Number of Points','Translation Error (%)');
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'med_t','Median Translation Error',...
%     'Number of Points','Translation Error (%)');
% 
% i= 0; w= 300; h= 300;
% 
% figure('color','w','position',[w*i,400,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'mean_foc','Mean Focal Length Error',...
%     'Number of Points','Focal Length Error (%)');
% 
% figure('color','w','position',[w*i,400,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'med_foc','Median Focal Length Error',...
%     'Number of Points','Focal Length Error (%)');
% 
% figure('color','w','position',[w*i,400,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'mean_reproj','Mean Reproj Error',...
%     'Number of Points','Reproj Error (pixels)');
% 
% figure('color','w','position',[w*i,400,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'med_reproj','Median Reproj Error',...
%     'Number of Points','Reproj Error (pixels)');
% 
