function [f R t] = HOMO(XX,xx)
%direct linear transformation for projection matrix P
H = HomoNormDLT(XX(1:2,:),xx);

%estimate f using zhang's method
w = abs(-(H(3,1)*H(3,2)/(H(1,1)*H(1,2)+H(2,1)*H(2,2))));
f = sqrt(1/w);

Hw = [H(1,:)/f; H(2,:)/f; H(3,:)];

%decompose H into R, t
scale = mean([1/norm(Hw(:,1)), 1/norm(Hw(:,2))]);
t = Hw(:,3)*scale;
r1 = Hw(:,1)*scale;
r2 = Hw(:,2)*scale;
r3 = cross(r1,r2);
R = [r1 r2 r3];

%orthogonal constraint
[U S V] = svd(R);
R = U*V.';

if det(R) < 0
    R = [R(:,1) R(:,2) -R(:,3)];
end

%account for the sign ambiguity
proj = R*XX+t*ones(1,size(XX,2));
if mean(proj(3,:)) < 0
    t = -t;
    R = [-R(:,1) -R(:,2) R(:,3)];
end

end