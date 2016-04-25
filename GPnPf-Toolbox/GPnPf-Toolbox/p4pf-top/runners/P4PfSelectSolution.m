% select the best solution w.r.t. image reprojection error
%
% by Martin Bujnak (c) may 2012

function [f R t za zb zc zd] = P4PfSelectSolution(za, zb, zc, zd, fs, m2D, M3D)

    solsCnt = length(fs);
    if solsCnt == 0
       
        f = [];
        R = [];
        t = [];
        za = [];
        zb = [];
        zc = [];
        zd = [];        
        return;
    end
    
    maxErrSoFar = 1e+100;
    bestSolution = 1;
    R = zeros(3,3);
    t = zeros(3);
    for i=1:solsCnt
        
        % recover camera rotation and translation
        % create p3d points in a camera coordinate system (using depths)

        p3dc(:, 1) = za(i) * [m2D(:, 1); fs(i)];
        p3dc(:, 2) = zb(i) * [m2D(:, 2); fs(i)];
        p3dc(:, 3) = zc(i) * [m2D(:, 3); fs(i)];
        p3dc(:, 4) = zd(i) * [m2D(:, 4); fs(i)];

        % calc camera
        [Ri ti] = GetRigidTransform2(M3D, p3dc, [], false);
        
        M = Ri*M3D + repmat(ti,1, 4);
        M(1,:) = fs(i)*M(1,:) ./ (M(3,:));
        M(2,:) = fs(i)*M(2,:) ./ (M(3,:));
        
        err = sum(sum((M(1:2,:) - m2D).^2));
        
        if (err < maxErrSoFar)
            
            bestSolution = i;
            maxErrSoFar = err;
            R = Ri;
            t = ti;
        end
    end

    f = fs(bestSolution);
    za = za(bestSolution);
    zb = zb(bestSolution);    
    zc = zc(bestSolution);
    zd = zd(bestSolution);    
end