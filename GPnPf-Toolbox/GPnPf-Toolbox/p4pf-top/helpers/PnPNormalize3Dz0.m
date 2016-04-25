% transform 3D space so that 
%       X(3,1:3) = 0            ;; first three points are on the z=0 plane
%       X(:,1) = [0;0;0]        
%       X(:,2) = [1;0;0]        
%
% by Martin Bujnak, (c)may2010

function [X M] = PnPNormalize3Dz0(X)

    % make the first point 0,0,0
    X(4,:) = 1;
    T1 = [1 0 0 -X(1,1)
          0 1 0 -X(2,1)
          0 0 1 -X(3,1)
          0 0 0 1];
    X = T1*X;

    % rotate space to z = 0 & make first point unit
    r1 = cross(X(1:3,2), X(1:3,3));
    r1 = r1/norm(r1);
    r2 = cross(r1, X(1:3,2));
    r2 = r2/norm(r2);
    r3 = cross(r2, r1);

    s = 1/norm(X(1:3,2));
    sR = s*[r3 r2 r1];
    X = sR'*X(1:3,:);

    M = [sR' zeros(3,1);0 0 0 1]*T1;
end

