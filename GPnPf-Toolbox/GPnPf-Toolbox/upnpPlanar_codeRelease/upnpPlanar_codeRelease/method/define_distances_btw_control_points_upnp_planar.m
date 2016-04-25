function dsq=define_distances_btw_control_points_upnp_planar(Cw)

%relative coordinates of the control points
c1=Cw(1,:);
c2=Cw(2,:);
c3=Cw(3,:);

d12=(c1(1)-c2(1))^2 + (c1(2)-c2(2))^2 + (c1(3)-c2(3))^2;
d13=(c1(1)-c3(1))^2 + (c1(2)-c3(2))^2 + (c1(3)-c3(3))^2;
d23=(c2(1)-c3(1))^2 + (c2(2)-c3(2))^2 + (c2(3)-c3(3))^2;

dsq=[d12,d13,d23]';