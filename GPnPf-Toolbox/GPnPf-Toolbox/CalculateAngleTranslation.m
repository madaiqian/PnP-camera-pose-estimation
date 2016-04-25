function [Rt obj] = CalculateAngleTranslation(axis,f,XX,xx)
%the third column of Rx is rotation axis
%then R=Rx*[0 c  -s
%           0 s  c
%           1 0  0];

%construct the first part of R
axis = axis';
if abs([0 1 0]*axis) < abs([0 0 1]*axis)
    z= xcross(axis,[0; 1; 0]); z= z/norm(z);
    y= xcross(z, axis); y= y/norm(y);
else
    y= xcross([0; 0; 1], axis); y= y/norm(y);
    z= xcross(axis,y); z= z/norm(z);
end
Rx= [y z axis];

%reshufle XX
%XX = XX([2,3,1],:);

%construct the quadratic eigenvalue problem
npt = size(XX,2);
xx(1,:) = xx(1,:)/f;
xx(2,:) = xx(2,:)/f;
xx = Rx.'*xx;
xx = xx./repmat(xx(3,:),3,1);
AT = zeros(3,2*npt);
AT(1,1:2:end) = 1;
AT(2,2:2:end) = 1;
AT(3,1:2:end) = -xx(1,:);
AT(3,2:2:end) = -xx(2,:);

B = zeros(2*npt,3);
B(1:2:end,1) = -XX(2,:).';
B(2:2:end,1) = -XX(3,:).';
B(1:2:end,2) = XX(3,:).';
B(2:2:end,2) = -XX(2,:).';
B(1:2:end,3) = (xx(1,:).*XX(1,:)).';
B(2:2:end,3) = (xx(2,:).*XX(1,:)).';

ATB = AT*B; 
BTB = B.'*B;
ATA = [npt 0 -sum(xx(1,:)); 0 npt -sum(xx(2,:)); -sum(xx(1,:)) -sum(xx(2,:)) sum(xx(1,:).^2+xx(2,:).^2)];
invATAATB = ATA\ATB;
G = -ATB.'*invATAATB + BTB;

%solve the quadratic eigenvalue problem (correct)
H = G(1:2,1:2);
g = -G(1:2,3);
A2 = speye(2);
A1 = -2*H;
gn = g;
A0 = H*H - gn*gn';
% QEP
lambda = polyeig(A0,A1,A2); 

% smallest real eigenvalue
lambda = lambda(imag(lambda)==0);
if isempty(lambda)
    error('spherelsq: problem with floating point accuracy');
end
lambda = sort(lambda);

for i = 1:length(lambda)
    % solution
    cs = (H-lambda(i)*speye(size(H))) \ g;

    %the translation
    trans = Rx*invATAATB*[cs;1];

    %the rotation matrix
    R = Rx*[0 cs(1) -cs(2); 0 cs(2) cs(1); 1 0 0];

    %chirality conditions
    depth = R*XX + repmat(trans,1,npt);
    depth = depth(3,:); 
    if mean(depth) > 0
        Rt = [R trans];
        obj = [cs;1].'*G*[cs;1];
        break;
    else
        Rt = []; obj = [];
    end
end
end

function c = xcross(a,b)

c = [a(2)*b(3)-a(3)*b(2);
     a(3)*b(1)-a(1)*b(3);
     a(1)*b(2)-a(2)*b(1)];
 
end
 