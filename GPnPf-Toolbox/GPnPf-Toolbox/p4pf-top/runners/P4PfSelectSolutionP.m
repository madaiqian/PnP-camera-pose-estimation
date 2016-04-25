% select the best solution w.r.t. image reprojection error
%
% by Martin Bujnak (c) may 2012

function [f R t] = P4PfSelectSolutionP(fs, Rs, ts, m2D, M3D)

    if nargin < 4
        
        P = fs;
        m2D = Rs;
        M3D = ts;
        R = [];
        t = [];
        
        solsCnt = size(P,1)/3;
        if solsCnt < 2

            f = P;
            return;
        end
        
        M3D(4,:) = 1;
        
        maxErrSoFar = 1e+100;
        bestSolution = 1;
        for i=1:solsCnt

            M = P((i*3-2):(i*3), :)*M3D;
            M(1,:) = fs(i)*M(1,:) ./ (M(3,:));
            M(2,:) = fs(i)*M(2,:) ./ (M(3,:));

            err = sum(sum((M(1:2,:) - m2D).^2));

            if (err < maxErrSoFar)

                bestSolution = i;
                maxErrSoFar = err;
            end
        end

        f = P((bestSolution*3-2):(bestSolution*3), :);
        
    else

        solsCnt = length(fs);
        if solsCnt < 2

            f = fs;
            R = Rs;
            t = ts;
            return;
        end

        maxErrSoFar = 1e+100;
        bestSolution = 1;
        for i=1:solsCnt

            M = Rs(:,:,i)*M3D + repmat(ts(:,i),1, 4);
            M(1,:) = fs(i)*M(1,:) ./ (M(3,:));
            M(2,:) = fs(i)*M(2,:) ./ (M(3,:));

            err = sum(sum((M(1:2,:) - m2D).^2));

            if (err < maxErrSoFar)

                bestSolution = i;
                maxErrSoFar = err;
            end
        end

        f = fs(bestSolution);
        R = Rs(:,:,bestSolution);
        t = ts(:,bestSolution);    
    end
end