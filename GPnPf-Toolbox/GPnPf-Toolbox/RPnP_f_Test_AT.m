function [f_AT] = RPnP_f_Test_AT(XX,xx)
%f_AT: angle + translation

n= length(xx);
XXw= XX;

scale = max(max(abs(xx)));
xxv= [xx/scale; ones(1,n)];

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
scale2 = max(max(abs(XX)));
XX = XX/scale2;

% Dividing the n-point set into (n-2) 3-point subsets
% and setting up the P3P equations

u0 = xxv(1,i1); v0 = xxv(2,i1);
u1 = xxv(1,i2); v1 = xxv(2,i2);

%constant coeffients, depends only on the selection of axis
d01= norm(XX(:,i1)-XX(:,i2));

%using angle constraints+translation, 14 monomials, t1, f two variables (average version)
%nn=[ f^3*t1^4, f^2*t1^4, f*t1^4, t1^4, f^2*t1^3, f*t1^3, t1^3, f^2*t1^2, f*t1^2, t1^2, f*t1, t1, f, 1]

coeffMat = zeros(n-2,14);

j= 0;
for i= 1:n
    if i == i1 || i == i2
        continue;
    end
    j= j+1;

    u2= xxv(1,i); v2 = xxv(2,i);
    d02= norm(XX(:,i1)-XX(:,i));
    d12= norm(XX(:,i)-XX(:,i2));

    k = d02/d01; %the distance ratio
    cs = (d01^2+d02^2 - d12^2)/(2*d01*d02); %the cosine of the angle between do1 and do2

    mm1 = 16*k^2*(cs^2 - 1);
    mm2 = (2*u0 + u2*(4*cs*k - 2))^2 - 4*k^2*((u0 + u1)^2 + (v0 + v1)^2) - 4*k^2*(4*u2*(u0 + u1) + 4*v2*(v0 + v1)) + (2*v0 + v2*(4*cs*k - 2))^2 + 8*cs*k*(u2*(u0 + u1) - u0*(u0 + u1) - v0*(v0 + v1) + v2*(v0 + v1) + cs*k*((u0 + u1)^2 + (v0 + v1)^2));
    mm3 = (u2*(u0 + u1) - u0*(u0 + u1) - v0*(v0 + v1) + v2*(v0 + v1) + cs*k*((u0 + u1)^2 + (v0 + v1)^2))^2 - 2*(u2*(u0*(u0 + u1) + v0*(v0 + v1) - cs*k*((u0 + u1)^2 + (v0 + v1)^2)) - u0*(u2*(u0 + u1) + v2*(v0 + v1)))*(2*u0 + u2*(4*cs*k - 2)) - 2*(v2*(u0*(u0 + u1) + v0*(v0 + v1) - cs*k*((u0 + u1)^2 + (v0 + v1)^2)) - v0*(u2*(u0 + u1) + v2*(v0 + v1)))*(2*v0 + v2*(4*cs*k - 2)) - 4*k^2*(u2*(u0 + u1) + v2*(v0 + v1))^2 - k^2*(4*u2*(u0 + u1) + 4*v2*(v0 + v1))*((u0 + u1)^2 + (v0 + v1)^2);
    mm4 = (u2*(u0*(u0 + u1) + v0*(v0 + v1) - cs*k*((u0 + u1)^2 + (v0 + v1)^2)) - u0*(u2*(u0 + u1) + v2*(v0 + v1)))^2 + (v2*(u0*(u0 + u1) + v0*(v0 + v1) - cs*k*((u0 + u1)^2 + (v0 + v1)^2)) - v0*(u2*(u0 + u1) + v2*(v0 + v1)))^2 - k^2*(u2*(u0 + u1) + v2*(v0 + v1))^2*((u0 + u1)^2 + (v0 + v1)^2);
    mm5 = 16*u0*u2 + 16*v0*v2 - 8*u0^2 - 8*u2^2 - 8*v0^2 - 8*v2^2 + 8*k^2*u0^2 - 8*k^2*u1^2 + 8*k^2*v0^2 - 8*k^2*v1^2 + 16*cs*k*u0^2 + 16*cs*k*u2^2 + 16*cs*k*v0^2 + 16*cs*k*v2^2 + 16*k^2*u0*u2 - 16*k^2*u1*u2 + 16*k^2*v0*v2 - 16*k^2*v1*v2 - 16*cs^2*k^2*u0^2 + 16*cs^2*k^2*u1^2 - 16*cs^2*k^2*v0^2 + 16*cs^2*k^2*v1^2 - 32*cs*k*u0*u2 - 32*cs*k*v0*v2;
    mm6 = 2*(2*u0 - 2*u2)*(u2*(u0*(u0 + u1) + v0*(v0 + v1) - cs*k*((u0 + u1)^2 + (v0 + v1)^2)) - u0*(u2*(u0 + u1) + v2*(v0 + v1))) - 2*(u2*(u0 + u1) - u0*(u0 + u1) - v0*(v0 + v1) + v2*(v0 + v1) + cs*k*((u0 + u1)^2 + (v0 + v1)^2))*(u2*(u0 + u1) - u0*(u0 + u1) - v0*(v0 + v1) + v2*(v0 + v1) - u0*(u0 - u1) + u2*(u0 - u1) - v0*(v0 - v1) + v2*(v0 - v1) + cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))) + 2*(2*v0 - 2*v2)*(v2*(u0*(u0 + u1) + v0*(v0 + v1) - cs*k*((u0 + u1)^2 + (v0 + v1)^2)) - v0*(u2*(u0 + u1) + v2*(v0 + v1))) - 2*(2*u0 + u2*(4*cs*k - 2))*(u0*(u2*(u0 - u1) + v2*(v0 - v1)) + u0*(u2*(u0 + u1) + v2*(v0 + v1)) - u2*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1) - cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)))) - 2*(2*v0 + v2*(4*cs*k - 2))*(v0*(u2*(u0 - u1) + v2*(v0 - v1)) + v0*(u2*(u0 + u1) + v2*(v0 + v1)) - v2*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1) - cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)))) + 8*k^2*(u2*(u0 - u1) + v2*(v0 - v1))*(u2*(u0 + u1) + v2*(v0 + v1)) + k^2*(4*u2*(u0 - u1) + 4*v2*(v0 - v1))*((u0 + u1)^2 + (v0 + v1)^2) + k^2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(4*u2*(u0 + u1) + 4*v2*(v0 + v1));
    mm7 = 2*(u2*(u0*(u0 + u1) + v0*(v0 + v1) - cs*k*((u0 + u1)^2 + (v0 + v1)^2)) - u0*(u2*(u0 + u1) + v2*(v0 + v1)))*(u0*(u2*(u0 - u1) + v2*(v0 - v1)) + u0*(u2*(u0 + u1) + v2*(v0 + v1)) - u2*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1) - cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)))) + 2*(v2*(u0*(u0 + u1) + v0*(v0 + v1) - cs*k*((u0 + u1)^2 + (v0 + v1)^2)) - v0*(u2*(u0 + u1) + v2*(v0 + v1)))*(v0*(u2*(u0 - u1) + v2*(v0 - v1)) + v0*(u2*(u0 + u1) + v2*(v0 + v1)) - v2*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1) - cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)))) + k^2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u2*(u0 + u1) + v2*(v0 + v1))^2 + 2*k^2*(u2*(u0 - u1) + v2*(v0 - v1))*(u2*(u0 + u1) + v2*(v0 + v1))*((u0 + u1)^2 + (v0 + v1)^2);
    mm8 = (2*u0 - 2*u2)^2 + (2*v0 - 2*v2)^2 - 4*k^2*((u0 - u1)^2 + (v0 - v1)^2) + 8*cs*k*(u2*(u0 - u1) - u0*(u0 - u1) - v0*(v0 - v1) + v2*(v0 - v1) + cs*k*((u0 - u1)^2 + (v0 - v1)^2));
    mm9 = 2*(2*u0 - 2*u2)*(u0*(u2*(u0 - u1) + v2*(v0 - v1)) + u0*(u2*(u0 + u1) + v2*(v0 + v1)) - u2*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1) - cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)))) + 2*(2*v0 - 2*v2)*(v0*(u2*(u0 - u1) + v2*(v0 - v1)) + v0*(u2*(u0 + u1) + v2*(v0 + v1)) - v2*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1) - cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)))) + 2*(u2*(u0 + u1) - u0*(u0 + u1) - v0*(v0 + v1) + v2*(v0 + v1) + cs*k*((u0 + u1)^2 + (v0 + v1)^2))*(u2*(u0 - u1) - u0*(u0 - u1) - v0*(v0 - v1) + v2*(v0 - v1) + cs*k*((u0 - u1)^2 + (v0 - v1)^2)) + 2*(u0*(u2*(u0 - u1) + v2*(v0 - v1)) - u2*(u0*(u0 - u1) + v0*(v0 - v1) - cs*k*((u0 - u1)^2 + (v0 - v1)^2)))*(2*u0 + u2*(4*cs*k - 2)) + 2*(v0*(u2*(u0 - u1) + v2*(v0 - v1)) - v2*(u0*(u0 - u1) + v0*(v0 - v1) - cs*k*((u0 - u1)^2 + (v0 - v1)^2)))*(2*v0 + v2*(4*cs*k - 2)) + (u2*(u0 + u1) - u0*(u0 + u1) - v0*(v0 + v1) + v2*(v0 + v1) - u0*(u0 - u1) + u2*(u0 - u1) - v0*(v0 - v1) + v2*(v0 - v1) + cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)))^2 - 4*k^2*(u2*(u0 - u1) + v2*(v0 - v1))^2 - k^2*((u0 - u1)^2 + (v0 - v1)^2)*(4*u2*(u0 + u1) + 4*v2*(v0 + v1)) - k^2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(4*u2*(u0 - u1) + 4*v2*(v0 - v1));
    mm10 = (u0*(u2*(u0 - u1) + v2*(v0 - v1)) + u0*(u2*(u0 + u1) + v2*(v0 + v1)) - u2*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1) - cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))))^2 + (v0*(u2*(u0 - u1) + v2*(v0 - v1)) + v0*(u2*(u0 + u1) + v2*(v0 + v1)) - v2*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1) - cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))))^2 - 2*(u2*(u0*(u0 + u1) + v0*(v0 + v1) - cs*k*((u0 + u1)^2 + (v0 + v1)^2)) - u0*(u2*(u0 + u1) + v2*(v0 + v1)))*(u0*(u2*(u0 - u1) + v2*(v0 - v1)) - u2*(u0*(u0 - u1) + v0*(v0 - v1) - cs*k*((u0 - u1)^2 + (v0 - v1)^2))) - 2*(v2*(u0*(u0 + u1) + v0*(v0 + v1) - cs*k*((u0 + u1)^2 + (v0 + v1)^2)) - v0*(u2*(u0 + u1) + v2*(v0 + v1)))*(v0*(u2*(u0 - u1) + v2*(v0 - v1)) - v2*(u0*(u0 - u1) + v0*(v0 - v1) - cs*k*((u0 - u1)^2 + (v0 - v1)^2))) - k^2*((u0 - u1)^2 + (v0 - v1)^2)*(u2*(u0 + u1) + v2*(v0 + v1))^2 - k^2*(u2*(u0 - u1) + v2*(v0 - v1))^2*((u0 + u1)^2 + (v0 + v1)^2) - 2*k^2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u2*(u0 - u1) + v2*(v0 - v1))*(u2*(u0 + u1) + v2*(v0 + v1));
    mm11 = k^2*(4*u2*(u0 - u1) + 4*v2*(v0 - v1))*((u0 - u1)^2 + (v0 - v1)^2) - 2*(2*u0 - 2*u2)*(u0*(u2*(u0 - u1) + v2*(v0 - v1)) - u2*(u0*(u0 - u1) + v0*(v0 - v1) - cs*k*((u0 - u1)^2 + (v0 - v1)^2))) - 2*(2*v0 - 2*v2)*(v0*(u2*(u0 - u1) + v2*(v0 - v1)) - v2*(u0*(u0 - u1) + v0*(v0 - v1) - cs*k*((u0 - u1)^2 + (v0 - v1)^2))) - 2*(u2*(u0 - u1) - u0*(u0 - u1) - v0*(v0 - v1) + v2*(v0 - v1) + cs*k*((u0 - u1)^2 + (v0 - v1)^2))*(u2*(u0 + u1) - u0*(u0 + u1) - v0*(v0 + v1) + v2*(v0 + v1) - u0*(u0 - u1) + u2*(u0 - u1) - v0*(v0 - v1) + v2*(v0 - v1) + cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)));
    mm12 = k^2*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1))*(u2*(u0 - u1) + v2*(v0 - v1))^2 - 2*(v0*(u2*(u0 - u1) + v2*(v0 - v1)) - v2*(u0*(u0 - u1) + v0*(v0 - v1) - cs*k*((u0 - u1)^2 + (v0 - v1)^2)))*(v0*(u2*(u0 - u1) + v2*(v0 - v1)) + v0*(u2*(u0 + u1) + v2*(v0 + v1)) - v2*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1) - cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)))) - 2*(u0*(u2*(u0 - u1) + v2*(v0 - v1)) - u2*(u0*(u0 - u1) + v0*(v0 - v1) - cs*k*((u0 - u1)^2 + (v0 - v1)^2)))*(u0*(u2*(u0 - u1) + v2*(v0 - v1)) + u0*(u2*(u0 + u1) + v2*(v0 + v1)) - u2*(u0*(u0 + u1) + v0*(v0 + v1) + u0*(u0 - u1) + v0*(v0 - v1) - cs*k*(2*(u0 + u1)*(u0 - u1) + 2*(v0 + v1)*(v0 - v1)))) + 2*k^2*(u2*(u0 - u1) + v2*(v0 - v1))*((u0 - u1)^2 + (v0 - v1)^2)*(u2*(u0 + u1) + v2*(v0 + v1));
    mm13 = (u2*(u0 - u1) - u0*(u0 - u1) - v0*(v0 - v1) + v2*(v0 - v1) + cs*k*((u0 - u1)^2 + (v0 - v1)^2))^2;
    mm14 = (u0*(u2*(u0 - u1) + v2*(v0 - v1)) - u2*(u0*(u0 - u1) + v0*(v0 - v1) - cs*k*((u0 - u1)^2 + (v0 - v1)^2)))^2 + (v0*(u2*(u0 - u1) + v2*(v0 - v1)) - v2*(u0*(u0 - u1) + v0*(v0 - v1) - cs*k*((u0 - u1)^2 + (v0 - v1)^2)))^2 - k^2*(u2*(u0 - u1) + v2*(v0 - v1))^2*((u0 - u1)^2 + (v0 - v1)^2);

    %calculate coefficients of the polynomial equations
    coeffMat(j,:) = [mm1,mm2,mm3,mm4,mm5,mm6,mm7,mm8,mm9,mm10,mm11,mm12,mm13,mm14];
end

%p4p+f problem
if n == 4 %4 points
    [f_AT t_AT] = PnP_f_Solver_4order_f3t4_Transt_1_t(coeffMat(1,:),coeffMat(2,:));
end

%denormalize f 
f_AT = sqrt(f_AT)*scale;

end

function c = xcross(a,b)

c = [a(2)*b(3)-a(3)*b(2);
     a(3)*b(1)-a(1)*b(3);
     a(1)*b(2)-a(2)*b(1)];
 
end
 