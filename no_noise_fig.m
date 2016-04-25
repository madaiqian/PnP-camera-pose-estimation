data=load('no_noise.txt');
num=11;
[h1,b1]=hist(log10(data(:,1)),num);
[h2,b2]=hist(log10(data(:,2)),num);
[h3,b3]=hist(log10(data(:,3)),num);
[h4,b4]=hist(log10(data(:,4)),num);
[h5,b5]=hist(log10(data(:,5)),num);
[h6,b6]=hist(log10(data(:,6)),num);
h1=h1./sum(h1);
h2=h2./sum(h2);
h3=h3./sum(h3);
h4=h4./sum(h4);
h5=h5./sum(h5);
h6=h6./sum(h6);


figure('color','w','position',[100 100 360 320]);
hold all;
box on;

wide=2.5;
plot(b1,h1,'r-','linewidth',wide);hold on;
plot(b2,h2,'b-','linewidth',wide);hold on;
plot(b3,h3,'g-','linewidth',wide);hold on;
plot(b4,h4,'c-','linewidth',wide);hold on;
plot(b5,h5,'y-','linewidth',wide);hold on;
plot(b6,h6,'k-','linewidth',wide);hold on;

legend('P3.5P','Bujnak','Ratio','Angle','Zheng-polyeig','Zheng-CP');
title('4000 no-noise experiments');
ylabel('number(%)');
xlabel('accuracy(  log10(err) )');


yticklabel={};
for ii=0:10:100
    yticklabel=[yticklabel,sprintf('%d%%',ii)];
end
set(gca,'yticklabel',yticklabel);