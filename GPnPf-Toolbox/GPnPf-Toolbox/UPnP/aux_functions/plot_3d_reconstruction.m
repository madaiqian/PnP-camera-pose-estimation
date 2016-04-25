function plot_3d_reconstruction(point,method,A,R,T,h)

figure(h); hold on;
n=size(point,2);
%err=0;
for i=1:n
   plot3(point(i).Xcam(1),point(i).Xcam(2),point(i).Xcam(3),'.','color',[1 0 0]);
   plot3(point(i).Xcam_est(1),point(i).Xcam_est(2),point(i).Xcam_est(3),'o','color',[0 0 1]);
end
txt=strcat('Reconstr. Error');
title(txt,'fontsize',20);
grid on;  

draw_camera(A,R,T,'',1,h);
view(66,50);