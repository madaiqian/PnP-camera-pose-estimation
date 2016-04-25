clear; clc;
addpath(genpath('rpnp1.0'));

% experimental parameters
npt= 20;
num= 4000;

err_list1 = zeros(1,num);
err_list2 = zeros(1,num);
err_list3 = zeros(1,num);

% experiments
for i= 1:num
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
    
    %call various solvers
    [f_R f_A f_AT] = RPnP_f_Test(XXw,xx);
    
    err_list1(i) = min(abs(f_R-f)/f);
    err_list2(i) = min(abs(f_A-f)/f);
    err_list3(i) = min(abs(f_AT-f)/f);

end

temp = err_list3; temp(temp<=0) = [];
minerror = min(temp);
err_list3(err_list3<=0) = minerror;   

  
[hh1,bb1] = hist(log10(err_list1),20);
[hh2,bb2] = hist(log10(err_list2),20);
[hh3,bb3] = hist(log10(err_list3),20);

figure('color','w','position',[100 100 360 320]);
hold all;
box on;

plot(bb1,hh1,'r-','linewidth',3);hold on;
plot(bb2,hh2,'b-','linewidth',3); hold on;
plot(bb3,hh3,'g-','linewidth',3); hold on;

xlim([-18 0]);
xtick= -18:2:0;
set(gca,'xtick',xtick);

%xlabel('Log_{10} Absolute Unit Quaternion Error','FontSize',11);
%ylabel('Number of Counts','FontSize',11);
legend('Ratio-Bi','Angle-Bi', 'GPnPf');

rmpath(genpath('rpnp1.0'));
