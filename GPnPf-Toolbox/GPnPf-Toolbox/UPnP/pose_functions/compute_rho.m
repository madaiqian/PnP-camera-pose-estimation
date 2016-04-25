function rho=compute_rho(Cw)

c0=Cw(1,:);
c1=Cw(2,:);
c2=Cw(3,:);
c3=Cw(4,:);

rho(1) = dist2(c0,c1);
rho(2) = dist2(c0,c2);
rho(3) = dist2(c0,c3);
rho(4) = dist2(c1,c2);
rho(5) = dist2(c1,c3);
rho(6) = dist2(c2,c3);

