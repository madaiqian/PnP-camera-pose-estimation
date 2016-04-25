% calculate f+RT given 4 projections of 3D points laying on a plane
%
% by Martin Bujnak, (c) 2010

function [f RT] = CalcP4PfPlanar(X, m)

    [X T] = PnPNormalize3Dz0(X);

    X = X(1:3,:);
    X(3,:)=1;
    m(3,:)=1;
    [H] = CalcH4pt(X, m);

    a = H(:,1);
    b = H(:,2);
    f2 = -(a(1)*b(1) + a(2)*b(2))/(a(3)*b(3));
    f = sqrt(f2);
    a(3) = a(3)*f;
    b(3) = b(3)*f;
    scale = 1/norm(a);
    a = a*scale;
    b = b*scale;
    RT = [a b cross(a, b) scale*[H(1:2, 3); f*H(3,3)]] * T;

end