clear; clc;
addpath(genpath('rpnp1.0'));
addpath(genpath('p4pf-top'))

% experimental parameters
npt= 4;
num= 4000;

err_list1 = zeros(1,num);
err_list2 = zeros(1,num);
err_list3 = zeros(1,num);
err_list4 = zeros(1,num);

fail_index=[];
% experiments
for i= 1:num
    i
    % camera's parameters
    width= 640;
    height= 480;
    f= rand(1)*1800+200; %random focal length in [200,2000]
    %f= 800;
    
    % generate 3d coordinates in camera space
    Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
    t= mean(Xc,2);
    R= rodrigues(randn(3,1));
    XXw= inv(R)*(Xc-repmat(t,1,npt));

    % projection
    xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
    
    %best distance solver using GB
    [f_BD] = Best_Distance_GB(XXw,xx);
    
    %best ratio solver using GB
    [f_BR] = Best_Ratio_GB(XXw,xx);
    
    %ratio solver using HVR
    [f_RHVR] = Best_Ratio_HVR(XXw,xx);
    
    %our solver (angle+trans)
    [f_AT] = RPnP_f_Test_AT(XXw,xx);
    
    if length(f_BD)<1 || length(f_BR)<1 || length(f_RHVR)<1 || length(f_AT)<1 
        fail_index = [fail_index i];
    else    
        err_list1(i) = min(abs(f_BD-f)/f);
        err_list2(i) = min(abs(f_BR-f)/f);
        err_list3(i) = min(abs(f_RHVR-f)/f);
        err_list4(i) = min(abs(f_AT-f)/f);
    end
 
end
err_list1(fail_index)=[];
err_list2(fail_index)=[];
err_list3(fail_index)=[];
err_list4(fail_index)=[];

%clamp those zero errors to the smallest nonzero values
temp = err_list1; temp(temp<=0) = [];
minerror = min(temp);
err_list1(err_list1<=0) = minerror;  

temp = err_list2; temp(temp<=0) = [];
minerror = min(temp);
err_list2(err_list2<=0) = minerror;  

temp = err_list3; temp(temp<=0) = [];
minerror = min(temp);
err_list3(err_list3<=0) = minerror;   

temp = err_list4; temp(temp<=0) = [];
minerror = min(temp);
err_list4(err_list4<=0) = minerror;  
    
[hh1,bb1] = hist(log10(err_list1),20);
[hh2,bb2] = hist(log10(err_list2),20);
[hh3,bb3] = hist(log10(err_list3),20);
[hh4,bb4] = hist(log10(err_list4),20);

figure('color','w','position',[100 100 360 320]);
hold all;
box on;

plot(bb1,hh1,'r--','linewidth',3);hold on;
plot(bb2,hh2,'b-.','linewidth',3); hold on;
plot(bb3,hh3,'m-','linewidth',3); hold on;
plot(bb4,hh4,'g-','linewidth',3); hold on;

xlim([-18 0]);
xtick= -18:2:0;
set(gca,'xtick',xtick);

%xlabel('Log_{10} Absolute Unit Quaternion Error','FontSize',11);
%ylabel('Number of Counts','FontSize',11);
legend('Dist-Best','Ratio-Best', 'Ratio-HVR','GP4Pf');

rmpath(genpath('rpnp1.0'));
rmpath(genpath('p4pf-top'))