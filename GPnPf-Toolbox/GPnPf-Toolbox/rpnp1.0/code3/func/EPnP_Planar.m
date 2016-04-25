function [Rr,tr] = EPnP_Planar(XX,xx)

[R,t]= efficient_pnp_planar(XX.',xx.',diag([1 1 1]));

%in the source code of EPnP, the planar case was not carefully addressed. 
%It is known that there are two symmetric (globally optimal) solutions in the
%planar case. One is in front of the camera, while the other is behind the
%camera, thus physically infeasible. 
%EPnP did not consider this situation. 

%If [r1 r2 r3] is a solution, so is [-r1 -r2 r3].
Rr(:,:,1) = R;
Rr(:,:,2) = [-R(:,1) -R(:,2) R(:,3)];

%for each rotation, calculate the translation vector
n = size(XX,2);
M = zeros(2*n,3);
v = zeros(2*n,1);

M(1:2:end,1) = 1;
M(2:2:end,2) = 1;
M(:,3) = -xx(:);

for i = 1:n
    v(2*(i-1)+1) = (xx(1,i)*Rr(3,:,1)-Rr(1,:,1))*XX(:,i);
    v(2*(i-1)+2) = (xx(2,i)*Rr(3,:,1)-Rr(2,:,1))*XX(:,i);
end
tr(:,1) = pinv(M)*v;

tr(:,2) = -tr(:,1);

%remove the one behind the camera
index = [];
for i = 1:2
    temp = Rr(:,:,i)*XX + repmat(tr(:,i),1,n);
    if ~isempty(find(temp(3,:)<0))
        index = [index i];
    end
end
Rr(:,:,index) = []; tr(:,index) = [];
return