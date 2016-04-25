%
% p4p-f using quaternions, formulation with [wtz;0;wtz] in translation
% vector
%
% by Martin Bujnak (c) 2012
%
% correct solver plugins - p4pf_quat2.m p4pf_quat2_wdcb.m p4pf_quat3_wdcb.m 

function [fs Ps meta Mr] = P4PfQuat(pSolver, m, M)

    % transform point so that:
    %  M(1) = [0 0 0 1];
    %  M(2) = [1 0 0 1];
    [Mz Tz] = PnPNormalize3Dz0(M);

    % rotate and scale 2D points so that 
    %  u(1) = [1 0]

    a = norm(m(1:2,1));
    nx = m(:,1) / a;
    nx(3) = 0;
    s = -nx(2);
    c = nx(1);
    a=1/a;
    R = [a*c a*-s 0
         a*s a*c 0
         0 0 1];
    u = R*m;
    scl=a;

    %Pz = R*Pgt*inv(Tz);
    %Pz / Pz(end)

    [ws bs ds cs meta Mr] = pSolver(u(1:2,1:4), Mz(1:3,1:4));

    ith = 0;
    cnt = length(ws);
    fs = zeros(cnt);
    Ps = zeros(3, 4, cnt);
    for s = 1:cnt

        w = ws(s);
        if (w < 0)
            continue;
        end
        
        a = 1;
        b = bs(s);
        c = cs(s);
        d = ds(s);

        wtz = -(a^2+b^2-c^2-d^2) + u(1,2)/u(2,2)*(2*a*d+2*b*c);

        P = [ a^2+b^2-c^2-d^2    2*b*c-2*a*d     2*a*c+2*b*d          wtz
              2*a*d+2*b*c      a^2-b^2+c^2-d^2   2*c*d-2*a*b           0
              w*(2*b*d-2*a*c)  w*(2*a*b+2*c*d)   w*(a^2-b^2-c^2+d^2)  wtz]; 

        ith = ith + 1;  
        Ps(:,:,ith) = inv(R)*P*Tz;
        fs(ith) = 1/(w*scl);
    end
    
    fs = fs(1:ith);
    Ps = Ps(:,:,1:ith);
end