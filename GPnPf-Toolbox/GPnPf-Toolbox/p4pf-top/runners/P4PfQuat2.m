%
% p4p-f using quaternions, formulation with [wtz;0;wtz] in translation
% vector
%
% by Martin Bujnak (c) 2012
%
% plugins - p4pf_quat101.m

function [fs Ps meta M] = P4PfQuat2(pSolver, m, M)

    % transform point so that:
    %  M(1) = [0 0 0 1];
    %  M(2) = [1 0 0 1];
    [Mz Tz] = PnPNormalize3Dz0(M);

    % rotate and scale 2D points so that 
    %  u(2) = [1 0]

    a = norm(m(1:2,2));
    nx = m(:,2) / a;
    nx(3) = 0;
    s = -nx(2);
    c = nx(1);
    a=1/a;
    R = [a*c a*-s 0
         a*s a*c 0
         0 0 1];
    u = R*m;
    scl=a;

    [ws ds cs bs meta M] = pSolver(u(1:2,1:4), Mz(1:3,1:4));

    ith = 0;
    cnt = length(ws);
    fs = zeros(cnt);
    Ps = zeros(3, 4, cnt);
    for s = 1:length(ws)

        w = ws(s);
        if (w < 0)
            continue;
        end

        b = bs(s);
        c = cs(s);
        d = ds(s);

        P = [1+b^2-c^2-d^2,          2*b*c-2*d,          2*c+2*b*d,  u(1,1)/u(2,1)*(-2*d-2*b*c)
             2*d+2*b*c,          1-b^2+c^2-d^2,          2*c*d-2*b,         -2*d-2*b*c
             w*(2*b*d-2*c),      w*(2*b+2*c*d),  w*(1-b^2-c^2+d^2),    (-2*d-2*b*c)/u(2,1)];

        ith = ith + 1;
        Ps(:,:,ith) = inv(R)*P*Tz;
        fs(ith) = 1/(w*scl);
    end
    
    fs = fs(1:ith);
    Ps = Ps(:,:,1:ith);
end