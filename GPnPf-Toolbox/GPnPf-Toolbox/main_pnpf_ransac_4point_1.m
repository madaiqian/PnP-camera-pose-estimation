clear; clc;
addpath(genpath('rpnp1.0'));
addpath(genpath('p4pf-top'))

% experimental parameters
nl = 0.5; %noise level
inthr = 2; %threshold of inliers
outratio = 0.5; %outlier ratio

npt= 1000; %number of points
numransac= 200; %number of RANSAC samplings
numrept = 200; %number of repetion

count_list1 = zeros(numrept,numransac);
count_list2 = zeros(numrept,numransac);
count_list3 = zeros(numrept,numransac);
count_list4 = zeros(numrept,numransac);
count_list5 = zeros(numrept,numransac);

fail_index=[];
% experiments
for i= 1:numrept
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
    
    % add Gaussian noise on all points
    xxn= xx+randn(2,npt)*nl;
    
    % synthetize outliers
    ind = randsample(npt, npt*outratio);
    xxn(:,ind) = xxn(:,ind)+randn(2,length(ind))*50;
    
    maxnum1 = 0; maxnum2 = 0; maxnum3 = 0; maxnum4 = 0; maxnum5 = 0;
    fail_list = [];
    for j = 1:numransac
        %randomly choose 4 points 
        ind4 = randsample(npt, 4);
        
        try
            %call different solvers - the best distance solver
            [ff RR tt] = Best_Distance_GB(XXw(:,ind4),xxn(:,ind4));

            %count the number of inliers
            for k = 1:length(ff)
                proj = RR(:,:,k)*XXw + tt(:,k)*ones(1, npt);
                proj = proj./repmat(proj(3,:),3,1);
                err = xxn - ff(k)*proj(1:2,:);
                err = sqrt(sum(err.*err,1));
                tempnum = length(find(err<=inthr));
                if  tempnum > maxnum1
                    maxnum1 = tempnum;
                end
            end
            count_list1(i,j) = maxnum1;
        catch
            fprintf('The best distance solver encounters internal errors\n');
        end
        
        try
            %call different solvers - the best ratio solver
            [ff RR tt] = Best_Ratio_GB(XXw(:,ind4),xxn(:,ind4));

            %count the number of inliers
            for k = 1:length(ff)
                proj = RR(:,:,k)*XXw + tt(:,k)*ones(1, npt);
                proj = proj./repmat(proj(3,:),3,1);
                err = xxn - ff(k)*proj(1:2,:);
                err = sqrt(sum(err.*err,1));
                tempnum = length(find(err<=inthr));
                if  tempnum > maxnum2
                    maxnum2 = tempnum;
                end
            end
            count_list2(i,j) = maxnum2;        
        catch
            fprintf('The best ratio solver encounters internal errors\n');
        end
        
        try
             %call different solvers - the GP4Pf solver without GN (1 axis)
            [ff RR tt] = GPnP_f(XXw(:,ind4),xxn(:,ind4));

            %count the number of inliers
            for k = 1:length(ff)
                proj = RR(:,:,k)*XXw + tt(:,k)*ones(1, npt);
                proj = proj./repmat(proj(3,:),3,1);
                err = xxn - ff(k)*proj(1:2,:);
                err = sqrt(sum(err.*err,1));
                tempnum = length(find(err<=inthr));
                if  tempnum > maxnum3
                    maxnum3 = tempnum;
                end
            end
            count_list3(i,j) = maxnum3;
        catch
            fprintf('The GP4Pf solver encounters internal errors\n');
        end
        
        try
             %call different solvers - the GP4Pf solver with GN
            [ff RR tt] = GPnP_f_GN(XXw(:,ind4),xxn(:,ind4));

            %count the number of inliers
            for k = 1:length(ff)
                proj = RR(:,:,k)*XXw + tt(:,k)*ones(1, npt);
                proj = proj./repmat(proj(3,:),3,1);
                err = xxn - ff(k)*proj(1:2,:);
                err = sqrt(sum(err.*err,1));
                tempnum = length(find(err<=inthr));
                if  tempnum > maxnum4
                    maxnum4 = tempnum;
                end
            end
            count_list4(i,j) = maxnum4;
        catch
            fprintf('The GP4Pf+GN solver encounters internal errors\n');
        end
        
        try
            %call different solvers - the RPnP (n=4) solver with known f
            [ff RR tt] = RPnP_interface(XXw(:,ind4),xxn(:,ind4),diag([f,f,1]));

            %count the number of inliers
            for k = 1:length(ff)
                proj = RR(:,:,k)*XXw + tt(:,k)*ones(1, npt);
                proj = proj./repmat(proj(3,:),3,1);
                err = xxn - ff(k)*proj(1:2,:);
                err = sqrt(sum(err.*err,1));
                tempnum = length(find(err<=inthr));
                if  tempnum > maxnum5
                    maxnum5 = tempnum;
                end
            end
            count_list5(i,j) = maxnum5;
        catch
            fprintf('The RPnP(n=4) solver encounters internal errors\n');
        end
        

    end
    %
    fprintf('The %d-th iteration finished\n', i);
end

figure('color','w','position',[100 100 360 320]);
hold all;
box on;

index = [1:10:numransac numransac];
plot(index,mean(count_list5(:,index),1),'m-','linewidth',3); hold on;
plot(index,mean(count_list1(:,index),1),'r--','linewidth',3);hold on;
plot(index,mean(count_list2(:,index),1),'b-.','linewidth',3); hold on;
plot(index,mean(count_list3(:,index),1),'g-','linewidth',3); hold on;
plot(index,mean(count_list4(:,index),1),'c-','linewidth',3); hold on;

xlim([1 numransac]);
xtick= 1:numransac/10:numransac+1;
xtick(2:end) = xtick(2:end)-1;
set(gca,'xtick',xtick);

xlabel('RANSAC Iteration','FontSize',11);
ylabel('Number of Inliers','FontSize',11);
legend('RPnP(n=4)','Dist-Best','Ratio-Best', 'GP4Pf','GP4Pf+GN');

rmpath(genpath('rpnp1.0'));
rmpath(genpath('p4pf-top'))