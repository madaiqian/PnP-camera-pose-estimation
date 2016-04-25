function draw_noisy_input_data(point)

n=size(point,2);

figure; hold on;
for i=1:n
    plot(point(i).Ximg_pix_true(1),point(i).Ximg_pix_true(2),'.','color',[1 0 0]);
    plot(point(i).Ximg_pix(1),point(i).Ximg_pix(2),'o','color',[0 0 1],'markersize',5);
end
title('Noise in image plane','fontsize',20);
grid on;
legend('True 2D','Observation');
