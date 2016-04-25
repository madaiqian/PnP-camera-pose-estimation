% COMPUTE_RHO
function rho = UPnP_compute_rho(Cw)

c1=Cw(1,:);
c2=Cw(2,:);
c3=Cw(3,:);
c4=Cw(4,:);

rho(1) = dist2(c1,c2);
rho(2) = dist2(c1,c3);
rho(3) = dist2(c1,c4);
rho(4) = dist2(c2,c3);
rho(5) = dist2(c2,c4);
rho(6) = dist2(c3,c4);

%% Finished
