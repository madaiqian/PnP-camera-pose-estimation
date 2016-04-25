function [mask] = RANSAC_GP4Pf(XXw,xxn,inthr,maxItr)
npt = size(xxn,2);

mask = zeros(1,npt);
maxInlier = 0;

for j = 1:maxItr
    %randomly choose 4 points 
    ind4 = randsample(npt, 4);
    
    %call different solvers - the GP4Pf solver without GN (1 axis)
    [ff RR tt] = GPnP_f(XXw(:,ind4),xxn(:,ind4));

    %count the number of inliers
    for k = 1:length(ff)
        proj = RR(:,:,k)*XXw + tt(:,k)*ones(1, npt);
        proj = proj./repmat(proj(3,:),3,1);
        err = xxn - ff(k)*proj(1:2,:);
        err = sqrt(sum(err.*err,1));
        
        tempnum = length(find(err<=inthr));
        if  tempnum > maxInlier
            maxInlier = tempnum;
            mask(find(err<inthr))=1;
        end
    end
end