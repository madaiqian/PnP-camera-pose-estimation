function [f_r R_r t_r] = GPnP_f_GN(XX,xx,KK)
%the average depth of p0 and p1 is 1, so that t0+t1=2
%translate t1 to t1+1 (14 monomials)
%translate f to f-(u0^2+v0^2+u1^2+v1^2)/2 (f=0 is translated)

if nargin < 3
    KK = zeros(3,3);
end

%translate the principle point
xx(1,:) = xx(1,:) - KK(1,3);
xx(2,:) = xx(2,:) - KK(2,3);

n = length(xx);
XXw = XX;

%normalize 2D coordinate
scale = max(max(abs(xx)));
xxv = [xx/scale; ones(1,n)];

% selecting an edge $P_{i1}P_{i2}$ by n random sampling
ri = ceil(rand(1,n)*n);
rj = ceil(rand(1,n)*n);

%2D distance
%dtemp = sum((xxv(:,ri) - xxv(:,rj)).^2,1);

%3D distance
dtemp = sum((XX(:,ri) - XX(:,rj)).^2,1);
[mindist index] = max(dtemp);
i1 = ri(index); i2 = rj(index);

%avoid the extreme case i1==i2
if i1 == i2
   if i1==1 
       i2 = i1+1;
   elseif i1==n
       i2 = i1-1;
   else
       i2 = i1+1;
   end
end       

% calculating the rotation matrix of $O_aX_aY_aZ_a$.
p1= XX(:,i1);
p2= XX(:,i2);
p0= (p1+p2)/2;
x= p2-p0; x= x/norm(x);
if abs([0 1 0]*x) < abs([0 0 1]*x)
    z= xcross(x,[0; 1; 0]); z= z/norm(z);
    y= xcross(z, x); y= y/norm(y);
else
    y= xcross([0; 0; 1], x); y= y/norm(y);
    z= xcross(x,y); z= z/norm(z);
end
Ro= [x y z];

% transforming the reference points form orignial object space 
% to the new coordinate frame  $O_aX_aY_aZ_a$.
XX= Ro.'*(XX-repmat(p0,1,n));
scale3d = max(max(abs(XX)));
XX = XX/scale3d;

% Dividing the n-point set into (n-2) 3-point subsets
% and setting up the P3P equations
u0 = xxv(1,i1); v0 = xxv(2,i1);
u1 = xxv(1,i2); v1 = xxv(2,i2);

%constant coeffients, depends only on the selection of axis
d01= norm(XX(:,i1)-XX(:,i2));

%constant parts w.r.t. u0, v0, u1, v1, d01
const1 = [16, -16];
const2 = [16, -16, 4, - 16*u0 - 16*u1, 24*u0 + 8*u1, -8*u0, 16, -16, 4, - 16*v0 - 16*v1, 24*v0 + 8*v1, -8*v0, 8*(u0 + u1)^2 + 8*(v0 + v1)^2, - 4*(u0 + u1)^2 - 4*(v0 + v1)^2, - 8*u0*(u0 + u1) - 8*v0*(v0 + v1), 4*u0^2 + 4*v0^2];
const3 = [8*(u0 + u1)^2 + 8*(v0 + v1)^2, -4*(u0 + u1)^2, - 8*v0*(v0 + v1) - 4*(u0 + u1)^2 - 4*(v0 + v1)^2, 4*v0*(v0 + v1) + (u0 + u1)^2, -8*(u0 + u1)*(v0 + v1), 8*v0*(u0 + u1) + 8*u0*(v0 + v1), 2*(u0 + u1)*(v0 + v1) - 4*u0*(v0 + v1) - 4*v0*(u0 + u1), -(4*u0 + 4*u1)*((u0 + u1)^2 + (v0 + v1)^2), 2*(u0 + u1)*((u0 + u1)^2 + (v0 + v1)^2) + 4*u0*((u0 + u1)^2 + (v0 + v1)^2), 4*v0^2*(u0 + u1) - 2*(u0 + u1)*(u0*(u0 + u1) + v0*(v0 + v1)) - 4*u0*v0*(v0 + v1), 8*(u0 + u1)^2 + 8*(v0 + v1)^2, -4*(v0 + v1)^2, - 8*u0*(u0 + u1) - 4*(u0 + u1)^2 - 4*(v0 + v1)^2, 4*u0*(u0 + u1) + (v0 + v1)^2, -(4*v0 + 4*v1)*((u0 + u1)^2 + (v0 + v1)^2), 2*(v0 + v1)*((u0 + u1)^2 + (v0 + v1)^2) + 4*v0*((u0 + u1)^2 + (v0 + v1)^2), 4*u0^2*(v0 + v1) - 2*(v0 + v1)*(u0*(u0 + u1) + v0*(v0 + v1)) - 4*u0*v0*(u0 + u1), ((u0 + u1)^2 + (v0 + v1)^2)^2, -2*(u0*(u0 + u1) + v0*(v0 + v1))*((u0 + u1)^2 + (v0 + v1)^2), (u0*(u0 + u1) + v0*(v0 + v1))^2];
const4 = [((u0 + u1)^2 + (v0 + v1)^2)^2, -(u0 + u1)^2*((u0 + u1)^2 + (v0 + v1)^2), -2*v0*(v0 + v1)*((u0 + u1)^2 + (v0 + v1)^2), v0^2*(u0 + u1)^2 + v0^2*(v0 + v1)^2, -2*(u0 + u1)*(v0 + v1)*((u0 + u1)^2 + (v0 + v1)^2), 2*v0*(u0 + u1)*((u0 + u1)^2 + (v0 + v1)^2) + 2*u0*(v0 + v1)*((u0 + u1)^2 + (v0 + v1)^2), - 2*u0*v0*(u0 + u1)^2 - 2*u0*v0*(v0 + v1)^2, ((u0 + u1)^2 + (v0 + v1)^2)^2, -(v0 + v1)^2*((u0 + u1)^2 + (v0 + v1)^2), -2*u0*(u0 + u1)*((u0 + u1)^2 + (v0 + v1)^2), u0^2*(u0 + u1)^2 + u0^2*(v0 + v1)^2];
const5 = [16, -8, 16*u0 - 16*u1, -32*u0, 16*u0, 16, -8, 16*v0 - 16*v1, -32*v0, 16*v0, 16*u1^2 - 16*u0^2 - 16*v0^2 + 16*v1^2, 8*u0^2 - 8*u1^2 + 8*v0^2 - 8*v1^2, 16*u0^2 + 16*v0^2, - 8*u0^2 - 8*v0^2];
const6 = [- 16*(u0 + u1)*(u0 - u1) - 16*(v0 + v1)*(v0 - v1), 8*(u0 + u1)*(u0 - u1), 8*v0*(v0 + v1) + 8*(u0 + u1)*(u0 - u1) + 8*(v0 + v1)*(v0 - v1) + 8*v0*(v0 - v1) + 4*(u0 + u1)^2 + 4*(v0 + v1)^2, - 4*u0*(u0 + u1) - 8*v0*(v0 + v1) - 4*v0*(v0 - v1), 8*(u0 + u1)*(v0 - v1) + 8*(v0 + v1)*(u0 - u1), - 8*v0*(u0 + u1) - 8*u0*(v0 + v1) - 8*v0*(u0 - u1) - 8*u0*(v0 - v1), 4*v0*(u0 + u1) + 4*u0*(v0 + v1) + 4*v0*(u0 - u1) + 4*u0*(v0 - v1), (4*u0 - 4*u1)*((u0 + u1)^2 + (v0 + v1)^2) + (2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(4*u0 + 4*u1), - 8*u0*((u0 + u1)^2 + (v0 + v1)^2) - 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 + u1) - 4*u0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)), 2*(u0 + u1)*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1)) + 4*u0*(u0*(u0 + u1) + v0*(v0 + v1)) - 4*v0^2*(u0 + u1) - 4*v0*(v0*(u0 + u1) + v0*(u0 - u1)) + 4*u0*(v0*(v0 + v1) + v0*(v0 - v1)) + 4*u0*v0*(v0 + v1), - 16*(u0 + u1)*(u0 - u1) - 16*(v0 + v1)*(v0 - v1), 8*(v0 + v1)*(v0 - v1), 8*u0*(u0 + u1) + 8*(u0 + u1)*(u0 - u1) + 8*(v0 + v1)*(v0 - v1) + 8*u0*(u0 - u1) + 4*(u0 + u1)^2 + 4*(v0 + v1)^2, - 8*u0*(u0 + u1) - 4*v0*(v0 + v1) - 4*u0*(u0 - u1), (4*v0 - 4*v1)*((u0 + u1)^2 + (v0 + v1)^2) + (2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(4*v0 + 4*v1), - 8*v0*((u0 + u1)^2 + (v0 + v1)^2) - 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(v0 + v1) - 4*v0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)), 2*(v0 + v1)*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1)) + 4*v0*(u0*(u0 + u1) + v0*(v0 + v1)) - 4*u0^2*(v0 + v1) + 4*v0*(u0*(u0 + u1) + u0*(u0 - u1)) - 4*u0*(u0*(v0 + v1) + u0*(v0 - v1)) + 4*u0*v0*(u0 + u1), -2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*((u0 + u1)^2 + (v0 + v1)^2), 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0*(u0 + u1) + v0*(v0 + v1)) + 2*((u0 + u1)^2 + (v0 + v1)^2)*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1)), -2*(u0*(u0 + u1) + v0*(v0 + v1))*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1))];
const7 = [-2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*((u0 + u1)^2 + (v0 + v1)^2), (2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 + u1)^2 + 2*(u0 + u1)*(u0 - u1)*((u0 + u1)^2 + (v0 + v1)^2), 2*(v0*(v0 + v1) + v0*(v0 - v1))*((u0 + u1)^2 + (v0 + v1)^2) + 2*v0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(v0 + v1), - 2*v0*(u0 + u1)*(v0*(u0 + u1) + v0*(u0 - u1)) - 2*v0*(v0*(v0 + v1) + v0*(v0 - v1))*(v0 + v1), 2*((u0 + u1)*(v0 - v1) + (v0 + v1)*(u0 - u1))*((u0 + u1)^2 + (v0 + v1)^2) + 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 + u1)*(v0 + v1), - 2*((u0 + u1)^2 + (v0 + v1)^2)*(v0*(u0 + u1) + v0*(u0 - u1)) - 2*((u0 + u1)^2 + (v0 + v1)^2)*(u0*(v0 + v1) + u0*(v0 - v1)) - 2*v0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 + u1) - 2*u0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(v0 + v1), 2*v0*(u0 + u1)*(u0*(u0 + u1) + u0*(u0 - u1)) + 2*u0*(u0 + u1)*(v0*(u0 + u1) + v0*(u0 - u1)) + 2*v0*(v0 + v1)*(u0*(v0 + v1) + u0*(v0 - v1)) + 2*u0*(v0*(v0 + v1) + v0*(v0 - v1))*(v0 + v1), -2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*((u0 + u1)^2 + (v0 + v1)^2), (2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(v0 + v1)^2 + 2*(v0 + v1)*(v0 - v1)*((u0 + u1)^2 + (v0 + v1)^2), 2*((u0 + u1)^2 + (v0 + v1)^2)*(u0*(u0 + u1) + u0*(u0 - u1)) + 2*u0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 + u1), - 2*u0*(u0 + u1)*(u0*(u0 + u1) + u0*(u0 - u1)) - 2*u0*(v0 + v1)*(u0*(v0 + v1) + u0*(v0 - v1))];
const8 = [4, 8*u0 - 8*u1, -8*u0, 4, 8*v0 - 8*v1, -8*v0, 8*(u0 - u1)^2 + 8*(v0 - v1)^2, - 4*(u0 - u1)^2 - 4*(v0 - v1)^2, - 8*u0*(u0 - u1) - 8*v0*(v0 - v1), 4*u0^2 + 4*v0^2];
const9 = [8*(u0 - u1)^2 + 8*(v0 - v1)^2, -4*(u0 - u1)^2, - 4*(u0 - u1)^2 - 4*(v0 - v1)^2 - 8*(u0 + u1)*(u0 - u1) - 8*(v0 + v1)*(v0 - v1) - 8*v0*(v0 - v1), 4*v0*(v0 + v1) + 2*(u0 + u1)*(u0 - u1) + 8*v0*(v0 - v1) + 4*u0^2, -8*(u0 - u1)*(v0 - v1), 8*v0*(u0 - u1) + 8*u0*(v0 - v1), 8*u0*v0 - 4*u0*(v0 + v1) - 4*v0*(u0 + u1) + 2*(u0 + u1)*(v0 - v1) + 2*(v0 + v1)*(u0 - u1) - 8*v0*(u0 - u1) - 8*u0*(v0 - v1), - (2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(4*u0 - 4*u1) - (4*u0 + 4*u1)*((u0 - u1)^2 + (v0 - v1)^2), 2*(u0 + u1)*((u0 - u1)^2 + (v0 - v1)^2) + 4*u0*((u0 - u1)^2 + (v0 - v1)^2) + 2*(u0 - u1)*((u0 + u1)^2 + (v0 + v1)^2) + 8*u0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)), 4*v0*(v0*(u0 + u1) + v0*(u0 - u1)) - 4*u0*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1)) - 2*(u0*(u0 + u1) + v0*(v0 + v1))*(u0 - u1) - 2*(u0 + u1)*(u0*(u0 - u1) + v0*(v0 - v1)) - 4*u0*(v0*(v0 + v1) + v0*(v0 - v1)) + 4*v0^2*(u0 - u1) - 4*u0*v0*(v0 - v1), 8*(u0 - u1)^2 + 8*(v0 - v1)^2, -4*(v0 - v1)^2, - 4*(u0 - u1)^2 - 4*(v0 - v1)^2 - 8*(u0 + u1)*(u0 - u1) - 8*(v0 + v1)*(v0 - v1) - 8*u0*(u0 - u1), 4*u0*(u0 + u1) + 2*(v0 + v1)*(v0 - v1) + 8*u0*(u0 - u1) + 4*v0^2, - (2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(4*v0 - 4*v1) - (4*v0 + 4*v1)*((u0 - u1)^2 + (v0 - v1)^2), 2*(v0 + v1)*((u0 - u1)^2 + (v0 - v1)^2) + 4*v0*((u0 - u1)^2 + (v0 - v1)^2) + 2*(v0 - v1)*((u0 + u1)^2 + (v0 + v1)^2) + 8*v0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)), 4*u0*(u0*(v0 + v1) + u0*(v0 - v1)) - 4*v0*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1)) - 2*(u0*(u0 + u1) + v0*(v0 + v1))*(v0 - v1) - 4*v0*(u0*(u0 + u1) + u0*(u0 - u1)) - 2*(v0 + v1)*(u0*(u0 - u1) + v0*(v0 - v1)) + 4*u0^2*(v0 - v1) - 4*u0*v0*(u0 - u1), (2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))^2 + 2*((u0 - u1)^2 + (v0 - v1)^2)*((u0 + u1)^2 + (v0 + v1)^2), - 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1)) - 2*(u0*(u0 - u1) + v0*(v0 - v1))*((u0 + u1)^2 + (v0 + v1)^2) - 2*((u0 - u1)^2 + (v0 - v1)^2)*(u0*(u0 + u1) + v0*(v0 + v1)), 2*(u0*(u0 - u1) + v0*(v0 - v1))*(u0*(u0 + u1) + v0*(v0 + v1)) + (u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1))^2];
const10 = [(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))^2 + 2*((u0 - u1)^2 + (v0 - v1)^2)*((u0 + u1)^2 + (v0 + v1)^2), - (u0 - u1)^2*((u0 + u1)^2 + (v0 + v1)^2) - (u0 + u1)^2*((u0 - u1)^2 + (v0 - v1)^2) - 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 + u1)*(u0 - u1), - 2*(v0*(v0 + v1) + v0*(v0 - v1))*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)) - 2*v0*(v0 - v1)*((u0 + u1)^2 + (v0 + v1)^2) - 2*v0*(v0 + v1)*((u0 - u1)^2 + (v0 - v1)^2), (v0*(u0 + u1) + v0*(u0 - u1))^2 + (v0*(v0 + v1) + v0*(v0 - v1))^2 + 2*v0^2*(u0 + u1)*(u0 - u1) + 2*v0^2*(v0 + v1)*(v0 - v1), - 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 + u1)*(v0 - v1) - 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(v0 + v1)*(u0 - u1) - 2*(u0 + u1)*(v0 + v1)*((u0 - u1)^2 + (v0 - v1)^2) - 2*(u0 - u1)*(v0 - v1)*((u0 + u1)^2 + (v0 + v1)^2), 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(v0*(u0 + u1) + v0*(u0 - u1)) + 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0*(v0 + v1) + u0*(v0 - v1)) + 2*v0*(u0 - u1)*((u0 + u1)^2 + (v0 + v1)^2) + 2*u0*(v0 - v1)*((u0 + u1)^2 + (v0 + v1)^2) + 2*v0*(u0 + u1)*((u0 - u1)^2 + (v0 - v1)^2) + 2*u0*(v0 + v1)*((u0 - u1)^2 + (v0 - v1)^2), - 2*(u0*(u0 + u1) + u0*(u0 - u1))*(v0*(u0 + u1) + v0*(u0 - u1)) - 2*(v0*(v0 + v1) + v0*(v0 - v1))*(u0*(v0 + v1) + u0*(v0 - v1)) - 4*u0*v0*(u0 + u1)*(u0 - u1) - 4*u0*v0*(v0 + v1)*(v0 - v1), (2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))^2 + 2*((u0 - u1)^2 + (v0 - v1)^2)*((u0 + u1)^2 + (v0 + v1)^2), - (v0 - v1)^2*((u0 + u1)^2 + (v0 + v1)^2) - (v0 + v1)^2*((u0 - u1)^2 + (v0 - v1)^2) - 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(v0 + v1)*(v0 - v1), - 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0*(u0 + u1) + u0*(u0 - u1)) - 2*u0*(u0 - u1)*((u0 + u1)^2 + (v0 + v1)^2) - 2*u0*(u0 + u1)*((u0 - u1)^2 + (v0 - v1)^2), (u0*(u0 + u1) + u0*(u0 - u1))^2 + (u0*(v0 + v1) + u0*(v0 - v1))^2 + 2*u0^2*(u0 + u1)*(u0 - u1) + 2*u0^2*(v0 + v1)*(v0 - v1)];
const11 = [4*(u0 - u1)^2 + 4*(v0 - v1)^2, - 4*u0*(u0 - u1) - 4*v0*(v0 - v1), (4*u0 - 4*u1)*((u0 - u1)^2 + (v0 - v1)^2), - 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 - u1) - 8*u0*((u0 - u1)^2 + (v0 - v1)^2), 4*u0*(u0*(u0 - u1) + v0*(v0 - v1)) + 2*(u0 - u1)*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1)) - 4*v0^2*(u0 - u1) + 4*u0*v0*(v0 - v1), 4*(u0 - u1)^2 + 4*(v0 - v1)^2, - 4*u0*(u0 - u1) - 4*v0*(v0 - v1), (4*v0 - 4*v1)*((u0 - u1)^2 + (v0 - v1)^2), - 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(v0 - v1) - 8*v0*((u0 - u1)^2 + (v0 - v1)^2), 4*v0*(u0*(u0 - u1) + v0*(v0 - v1)) + 2*(v0 - v1)*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1)) - 4*u0^2*(v0 - v1) + 4*u0*v0*(u0 - u1), -2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*((u0 - u1)^2 + (v0 - v1)^2), 2*((u0 - u1)^2 + (v0 - v1)^2)*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1)) + 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0*(u0 - u1) + v0*(v0 - v1)), -2*(u0*(u0 - u1) + v0*(v0 - v1))*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1))];
const12 = [-2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*((u0 - u1)^2 + (v0 - v1)^2), (2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 - u1)^2 + 2*(u0 + u1)*((u0 - u1)^2 + (v0 - v1)^2)*(u0 - u1), 2*(v0*(v0 + v1) + v0*(v0 - v1))*((u0 - u1)^2 + (v0 - v1)^2) + 2*v0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(v0 - v1), - 2*v0*(u0 - u1)*(v0*(u0 + u1) + v0*(u0 - u1)) - 2*v0*(v0*(v0 + v1) + v0*(v0 - v1))*(v0 - v1), 2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 - u1)*(v0 - v1) + 2*(u0 + u1)*((u0 - u1)^2 + (v0 - v1)^2)*(v0 - v1) + 2*(v0 + v1)*((u0 - u1)^2 + (v0 - v1)^2)*(u0 - u1), - 2*((u0 - u1)^2 + (v0 - v1)^2)*(v0*(u0 + u1) + v0*(u0 - u1)) - 2*((u0 - u1)^2 + (v0 - v1)^2)*(u0*(v0 + v1) + u0*(v0 - v1)) - 2*v0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 - u1) - 2*u0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(v0 - v1), 2*v0*(u0 - u1)*(u0*(u0 + u1) + u0*(u0 - u1)) + 2*u0*(u0 - u1)*(v0*(u0 + u1) + v0*(u0 - u1)) + 2*v0*(v0 - v1)*(u0*(v0 + v1) + u0*(v0 - v1)) + 2*u0*(v0*(v0 + v1) + v0*(v0 - v1))*(v0 - v1), -2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*((u0 - u1)^2 + (v0 - v1)^2), (2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(v0 - v1)^2 + 2*(v0 + v1)*((u0 - u1)^2 + (v0 - v1)^2)*(v0 - v1), 2*((u0 - u1)^2 + (v0 - v1)^2)*(u0*(u0 + u1) + u0*(u0 - u1)) + 2*u0*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u0 - u1), - 2*u0*(u0 - u1)*(u0*(u0 + u1) + u0*(u0 - u1)) - 2*u0*(v0 - v1)*(u0*(v0 + v1) + u0*(v0 - v1))];
const13 = [(u0 - u1)^2, 2*(u0 - u1)*(v0 - v1), 2*((u0 - u1)^2 + (v0 - v1)^2)*(u0 - u1), -2*(u0*(u0 - u1) + v0*(v0 - v1))*(u0 - u1), (v0 - v1)^2, 2*((u0 - u1)^2 + (v0 - v1)^2)*(v0 - v1), -2*(u0*(u0 - u1) + v0*(v0 - v1))*(v0 - v1), ((u0 - u1)^2 + (v0 - v1)^2)^2, -2*(u0*(u0 - u1) + v0*(v0 - v1))*((u0 - u1)^2 + (v0 - v1)^2), (u0*(u0 - u1) + v0*(v0 - v1))^2];
const14 = [((u0 - u1)^2 + (v0 - v1)^2)^2, -((u0 - u1)^2 + (v0 - v1)^2)*(u0 - u1)^2, -2*v0*((u0 - u1)^2 + (v0 - v1)^2)*(v0 - v1), v0^2*(u0 - u1)^2 + v0^2*(v0 - v1)^2, -2*((u0 - u1)^2 + (v0 - v1)^2)*(u0 - u1)*(v0 - v1), 2*v0*((u0 - u1)^2 + (v0 - v1)^2)*(u0 - u1) + 2*u0*((u0 - u1)^2 + (v0 - v1)^2)*(v0 - v1), - 2*u0*v0*(u0 - u1)^2 - 2*u0*v0*(v0 - v1)^2, ((u0 - u1)^2 + (v0 - v1)^2)^2, -((u0 - u1)^2 + (v0 - v1)^2)*(v0 - v1)^2, -2*u0*((u0 - u1)^2 + (v0 - v1)^2)*(u0 - u1), u0^2*(u0 - u1)^2 + u0^2*(v0 - v1)^2];

%varying parts w.r.t. k, cs, u2, v2
tempxx = xxv; tempxx(:,[i1,i2]) = [];
tempXX = XX; tempXX(:,[i1,i2]) = [];
u2= tempxx(1,:)'; v2 = tempxx(2,:)';
d02= sqrt(sum((repmat(XX(:,i1),1,n-2) - tempXX).^2,1)');
d12= sqrt(sum((repmat(XX(:,i2),1,n-2) - tempXX).^2,1)');

k = d02/d01; %the distance ratio
cs = (d01^2+d02.^2 - d12.^2)./(2*d01*d02); 
        
var1 = [cs.^2.*k.^2, k.^2];
var2 = [cs.^2.*k.^2.*u2.^2, cs.*k.*u2.^2, u2.^2, k.^2.*u2, cs.*k.*u2, u2, cs.^2.*k.^2.*v2.^2, cs.*k.*v2.^2, v2.^2, k.^2.*v2, cs.*k.*v2, v2, cs.^2.*k.^2, k.^2, cs.*k, ones(n-2,1)];
var3 = [cs.^2.*k.^2.*u2.^2, k.^2.*u2.^2, cs.*k.*u2.^2, u2.^2, k.^2.*u2.*v2, cs.*k.*u2.*v2, u2.*v2, k.^2.*u2, cs.*k.*u2, u2, cs.^2.*k.^2.*v2.^2, k.^2.*v2.^2, cs.*k.*v2.^2, v2.^2, k.^2.*v2, cs.*k.*v2, v2, cs.^2.*k.^2, cs.*k, ones(n-2,1)];
var4 = [cs.^2.*k.^2.*u2.^2, k.^2.*u2.^2, cs.*k.*u2.^2, u2.^2, k.^2.*u2.*v2, cs.*k.*u2.*v2, u2.*v2, cs.^2.*k.^2.*v2.^2, k.^2.*v2.^2, cs.*k.*v2.^2, v2.^2];
var5 = [cs.*k.*u2.^2, u2.^2, k.^2.*u2, cs.*k.*u2, u2, cs.*k.*v2.^2, v2.^2, k.^2.*v2, cs.*k.*v2, v2, cs.^2.*k.^2, k.^2, cs.*k, ones(n-2,1)];
var6 = [cs.^2.*k.^2.*u2.^2, k.^2.*u2.^2, cs.*k.*u2.^2, u2.^2, k.^2.*u2.*v2, cs.*k.*u2.*v2, u2.*v2, k.^2.*u2, cs.*k.*u2, u2, cs.^2.*k.^2.*v2.^2, k.^2.*v2.^2, cs.*k.*v2.^2, v2.^2, k.^2.*v2, cs.*k.*v2, v2, cs.^2.*k.^2, cs.*k, ones(n-2,1)];
var7 = [cs.^2.*k.^2.*u2.^2, k.^2.*u2.^2, cs.*k.*u2.^2, u2.^2, k.^2.*u2.*v2, cs.*k.*u2.*v2, u2.*v2, cs.^2.*k.^2.*v2.^2, k.^2.*v2.^2, cs.*k.*v2.^2, v2.^2];
var8 = [u2.^2, cs.*k.*u2, u2, v2.^2, cs.*k.*v2, v2, cs.^2.*k.^2, k.^2, cs.*k, ones(n-2,1)];
var9 = [cs.^2.*k.^2.*u2.^2, k.^2.*u2.^2, cs.*k.*u2.^2, u2.^2, k.^2.*u2.*v2, cs.*k.*u2.*v2, u2.*v2, k.^2.*u2, cs.*k.*u2, u2, cs.^2.*k.^2.*v2.^2, k.^2.*v2.^2, cs.*k.*v2.^2, v2.^2, k.^2.*v2, cs.*k.*v2, v2, cs.^2.*k.^2, cs.*k, ones(n-2,1)];
var10 = [cs.^2.*k.^2.*u2.^2, k.^2.*u2.^2, cs.*k.*u2.^2, u2.^2, k.^2.*u2.*v2, cs.*k.*u2.*v2, u2.*v2, cs.^2.*k.^2.*v2.^2, k.^2.*v2.^2, cs.*k.*v2.^2, v2.^2];
var11 = [cs.*k.*u2.^2, u2.^2, k.^2.*u2, cs.*k.*u2, u2, cs.*k.*v2.^2, v2.^2, k.^2.*v2, cs.*k.*v2, v2, cs.^2.*k.^2, cs.*k, ones(n-2,1)];
var12 = [cs.^2.*k.^2.*u2.^2, k.^2.*u2.^2, cs.*k.*u2.^2, u2.^2, k.^2.*u2.*v2, cs.*k.*u2.*v2, u2.*v2, cs.^2.*k.^2.*v2.^2, k.^2.*v2.^2, cs.*k.*v2.^2, v2.^2];
var13 = [u2.^2, u2.*v2, cs.*k.*u2, u2, v2.^2, cs.*k.*v2, v2, cs.^2.*k.^2, cs.*k, ones(n-2,1)];
var14 = [cs.^2.*k.^2.*u2.^2, k.^2.*u2.^2, cs.*k.*u2.^2, u2.^2, k.^2.*u2.*v2, cs.*k.*u2.*v2, u2.*v2, cs.^2.*k.^2.*v2.^2, k.^2.*v2.^2, cs.*k.*v2.^2, v2.^2];

%coefficient matrix
coeffMat = [var1*const1' var2*const2' var3*const3' var4*const4' var5*const5' var6*const6' var7*const7' var8*const8' var9*const9' var10*const10' var11*const11' var12*const12' var13*const13' var14*const14'];
            
if n == 4 %4 points
    [ftemp ttemp] = PnP_f_Solver_4order_f3t4_Transt_1_t(coeffMat(1,:),coeffMat(2,:));
elseif n > 4 %5 points or more
    %vec1 = [ f^5*t1^8, f^4*t1^8, f^4*t1^7, f^4*t1^6, f^3*t1^8, f^3*t1^7, f^3*t1^6, f^3*t1^5, f^3*t1^4, f^2*t1^8, f^2*t1^7, f^2*t1^6, f^2*t1^5, f^2*t1^4, f^2*t1^3, f^2*t1^2, f*t1^8, f*t1^7, f*t1^6, f*t1^5, f*t1^4, f*t1^3, f*t1^2, f*t1, f, t1^8, t1^7, t1^6, t1^5, t1^4, t1^3, t1^2, t1, 1]
    %vec2 = [ f^6*t1^7, f^5*t1^7, f^5*t1^6, f^5*t1^5, f^4*t1^7, f^4*t1^6, f^4*t1^5, f^4*t1^4, f^4*t1^3, f^3*t1^7, f^3*t1^6, f^3*t1^5, f^3*t1^4, f^3*t1^3, f^3*t1^2, f^3*t1, f^2*t1^7, f^2*t1^6, f^2*t1^5, f^2*t1^4, f^2*t1^3, f^2*t1^2, f^2*t1, f^2, f*t1^7, f*t1^6, f*t1^5, f*t1^4, f*t1^3, f*t1^2, f*t1, f, t1^7, t1^6, t1^5, t1^4, t1^3, t1^2, t1, 1]
    e1 = coeffMat(:,1);   e2 = coeffMat(:,2);   e3 = coeffMat(:,3);   e4 = coeffMat(:,4);   e5 = coeffMat(:,5);   e6 = coeffMat(:,6);   e7 = coeffMat(:,7);   e8 = coeffMat(:,8);   e9 = coeffMat(:,9);   e10 = coeffMat(:,10);
    e11 = coeffMat(:,11); e12 = coeffMat(:,12); e13 = coeffMat(:,13); e14 = coeffMat(:,14); 

    %derivative w.r.t. f (34 terms)
    temp = [ 6*e1.^2, 10*e1.*e2, 10*e1.*e5, 10*e1.*e8, 4*e2.^2 + 8*e1.*e3, 8*e1.*e6 + 8*e2.*e5, 4*e5.^2 + 8*e1.*e9 + 8*e2.*e8, 8*e1.*e11 + 8*e5.*e8, 4*e8.^2 + 8*e1.*e13, 6*e1.*e4 + 6*e2.*e3, 6*e1.*e7 + 6*e2.*e6 + 6*e3.*e5, 6*e1.*e10 + 6*e2.*e9 + 6*e3.*e8 + 6*e5.*e6, 6*e1.*e12 + 6*e2.*e11 + 6*e5.*e9 + 6*e6.*e8, 6*e1.*e14 + 6*e2.*e13 + 6*e5.*e11 + 6*e8.*e9, 6*e5.*e13 + 6*e8.*e11, 6*e8.*e13, 2*e3.^2 + 4*e2.*e4, 4*e2.*e7 + 4*e3.*e6 + 4*e4.*e5, 2*e6.^2 + 4*e2.*e10 + 4*e3.*e9 + 4*e4.*e8 + 4*e5.*e7, 4*e2.*e12 + 4*e3.*e11 + 4*e5.*e10 + 4*e6.*e9 + 4*e7.*e8, 2*e9.^2 + 4*e2.*e14 + 4*e3.*e13 + 4*e5.*e12 + 4*e6.*e11 + 4*e8.*e10, 4*e5.*e14 + 4*e6.*e13 + 4*e8.*e12 + 4*e9.*e11, 2*e11.^2 + 4*e8.*e14 + 4*e9.*e13, 4*e11.*e13, 2*e13.^2, 2*e3.*e4, 2*e3.*e7 + 2*e4.*e6, 2*e3.*e10 + 2*e4.*e9 + 2*e6.*e7, 2*e3.*e12 + 2*e4.*e11 + 2*e6.*e10 + 2*e7.*e9, 2*e3.*e14 + 2*e4.*e13 + 2*e6.*e12 + 2*e7.*e11 + 2*e9.*e10, 2*e6.*e14 + 2*e7.*e13 + 2*e9.*e12 + 2*e10.*e11, 2*e9.*e14 + 2*e10.*e13 + 2*e11.*e12, 2*e11.*e14 + 2*e12.*e13, 2*e13.*e14] ;
    coeff34_f = sum(temp,1);

    %derivative w.r.t. t (40 terms)
    temp = [ 8*e1.^2, 16*e1.*e2, 14*e1.*e5, 12*e1.*e8, 8*e2.^2 + 16*e1.*e3, 14*e1.*e6 + 14*e2.*e5, 6*e5.^2 + 12*e1.*e9 + 12*e2.*e8, 10*e1.*e11 + 10*e5.*e8, 4*e8.^2 + 8*e1.*e13, 16*e1.*e4 + 16*e2.*e3, 14*e1.*e7 + 14*e2.*e6 + 14*e3.*e5, 12*e1.*e10 + 12*e2.*e9 + 12*e3.*e8 + 12*e5.*e6, 10*e1.*e12 + 10*e2.*e11 + 10*e5.*e9 + 10*e6.*e8, 8*e1.*e14 + 8*e2.*e13 + 8*e5.*e11 + 8*e8.*e9, 6*e5.*e13 + 6*e8.*e11, 4*e8.*e13, 8*e3.^2 + 16*e2.*e4, 14*e2.*e7 + 14*e3.*e6 + 14*e4.*e5, 6*e6.^2 + 12*e2.*e10 + 12*e3.*e9 + 12*e4.*e8 + 12*e5.*e7, 10*e2.*e12 + 10*e3.*e11 + 10*e5.*e10 + 10*e6.*e9 + 10*e7.*e8, 4*e9.^2 + 8*e2.*e14 + 8*e3.*e13 + 8*e5.*e12 + 8*e6.*e11 + 8*e8.*e10, 6*e5.*e14 + 6*e6.*e13 + 6*e8.*e12 + 6*e9.*e11, 2*e11.^2 + 4*e8.*e14 + 4*e9.*e13, 2*e11.*e13, 16*e3.*e4, 14*e3.*e7 + 14*e4.*e6, 12*e3.*e10 + 12*e4.*e9 + 12*e6.*e7, 10*e3.*e12 + 10*e4.*e11 + 10*e6.*e10 + 10*e7.*e9, 8*e3.*e14 + 8*e4.*e13 + 8*e6.*e12 + 8*e7.*e11 + 8*e9.*e10, 6*e6.*e14 + 6*e7.*e13 + 6*e9.*e12 + 6*e10.*e11, 4*e9.*e14 + 4*e10.*e13 + 4*e11.*e12, 2*e11.*e14 + 2*e12.*e13, 8*e4.^2, 14*e4.*e7, 6*e7.^2 + 12*e4.*e10, 10*e4.*e12 + 10*e7.*e10, 4*e10.^2 + 8*e4.*e14 + 8*e7.*e12, 6*e7.*e14 + 6*e10.*e12, 2*e12.^2 + 4*e10.*e14, 2*e12.*e14];
    coeff40_t = sum(temp,1);  
   
    %normalize scale
    coeff34_f = coeff34_f/median(abs(coeff34_f));
    coeff40_t = coeff40_t/median(abs(coeff40_t));
    
    [ftemp ttemp] = PnP_f_Solver_NEW_Transt_1_t(coeff34_f,coeff40_t);
end

%de-translation
ttemp = ttemp + 1;

%keep solutions in a vaild interval
index = find((ttemp>0.1).*(ttemp<10).*(ftemp>0).*(ftemp<1000));
ttemp = ttemp(index);
ftemp = ftemp(index);

%in case of no solution
if isempty(ttemp)||isempty(ftemp)
    f_r = []; R_r=[]; t_r=[];
    return;
end

%denormalize f
ft = [sqrt(ftemp)*scale ttemp];

num = size(ft,1);
%determine the rotation axis 
U = (2-repmat(ft(:,2),1,3)).*[repmat(xx(1,i1),num,1) repmat(xx(2,i1),num,1) ft(:,1)];
V = repmat(ft(:,2),1,3).*[repmat(xx(1,i2),num,1) repmat(xx(2,i2),num,1) ft(:,1)];
axis = V - U; %the axis directing from p1 to p2
axis = axis./repmat(sqrt(sum(axis.^2,2)),1,3);

%calculate the rotation angle and translation
index = 1;
for i = 1:num
    [Rttemp] = CalculateAngleTranslation(axis(i,:),ft(i,1),XX*scale3d,[xx;ones(1,size(xx,2))]);
    if ~isempty(Rttemp)
        RR(:,:,index) = Rttemp(:,1:3)*Ro.';
        tt(:,index) = Rttemp(:,4)-Rttemp(:,1:3)*Ro.'*p0;
        ft_new(index,:) = ft(i,:);
        index = index+1;
    end
end

if index == 1
    f_r=[]; R_r=[]; t_r=[];
    return;
end

%the data matrix Q for refinement
%3D points after translation to centroid
Ucent = mean(XXw,2);
Um = XXw - repmat(Ucent,1,n);
xm = Um(1,:)'; ym = Um(2,:)'; zm = Um(3,:)';

xx = xx/scale;
ucent = mean(xx,2);
um = xx - repmat(ucent,1,n);
uxm = um(1,:).'; 
uym = um(2,:).';

uUx = repmat(xx(1,:),3,1).*Um;
uUy = repmat(xx(2,:),3,1).*Um;

wx = repmat(mean(uUx,2),1,n) - uUx;
wx1 = wx(1,:).'; wx2 = wx(2,:).'; wx3 = wx(3,:).'; 
wy = repmat(mean(uUy,2),1,n) - uUy;
wy1 = wy(1,:).'; wy2 = wy(2,:).'; wy3 = wy(3,:).'; 

%construct matrix N: 2n*19
N = zeros(2*n,19);
N(1:2:end,:) = [xm zeros(n,1) 2*zm -2*ym xm 2*ym 2*zm -xm zeros(n,1) -xm wx3 2*wx2 -2*wx1 -wx3 2*wx1 -wx3 2*wx2 wx3 -uxm];
N(2:2:end,:) = [ym -2*zm zeros(n,1) 2*xm -ym 2*xm zeros(n,1) ym 2*zm -ym wy3 2*wy2 -2*wy1 -wy3 2*wy1 -wy3 2*wy2 wy3 -uym];

MTM = N(:,1:18).'*N(:,1:18);
MTb = -N(:,1:18).'*N(:,19);
Q = N.'*N;
%Q = Q/max(max(abs(Q)));

%rotation matrix to quaternion
num = size(tt,2);
q_r = zeros(4,num); ff_r = zeros(1,num); obj_r = zeros(1,num);

for i  = 1:num
    depth = RR(3,:,i)*XXw+tt(3,i)*ones(1,n);
    q0 = matrix2quaternion(RR(:,:,i));
    q0 = q0*1/sqrt(mean(depth));

    [ff_r(i), q_r(:,i), obj_r(i)] = RefineDampedNewton(ft_new(i,1)/scale,q0,MTM,MTb,Q);
    
    %denormalize f
    ff_r(i) = ff_r(i)*scale;
    
    %recover rotation and translation
    a = q_r(1,i); b = q_r(2,i); c = q_r(3,i); d = q_r(4,i); 

    lambda = 1/(a^2+b^2+c^2+d^2);

    a = a*sqrt(lambda); b = b*sqrt(lambda); c = c*sqrt(lambda); d = d*sqrt(lambda);

    R_r = [a^2+b^2-c^2-d^2     2*b*c-2*a*d     2*b*d+2*a*c
         2*b*c+2*a*d         a^2-b^2+c^2-d^2 2*c*d-2*a*b
         2*b*d-2*a*c         2*c*d+2*a*b   a^2-b^2-c^2+d^2];    

    RR(:,:,i) = R_r;
    tt(:,i) = [scale*(lambda*ucent(1)+R_r(3,:)*mean(uUx,2))/ff_r(i)-R_r(1,:)*Ucent;
           scale*(lambda*ucent(2)+R_r(3,:)*mean(uUy,2))/ff_r(i)-R_r(2,:)*Ucent;
           lambda-R_r(3,:)*Ucent];
  
end

if n <= 6  %return all solutions
    f_r = ff_r;
    R_r = RR;
    t_r = tt;
else %return the one with smallest reprojection error
    %calculate the reprojection error 
    for i = 1:size(tt,2)
        proj = RR(:,:,i)*XXw + tt(:,i)*ones(1, n);
        proj = proj./repmat(proj(3,:),3,1);
        err(i) = norm(xx*scale - ff_r(i)*proj(1:2,:),'fro');
    end
    %choose the one with smallest error 
    [minerr index] = min(err);
    f_r = ff_r(index);
    R_r = RR(:,:,index);
    t_r = tt(:,index);    
end

end

function c = xcross(a,b)

c = [a(2)*b(3)-a(3)*b(2);
     a(3)*b(1)-a(1)*b(3);
     a(1)*b(2)-a(2)*b(1)];
 
end
 