clear;close all;clc;

load('F50 vs Force 1 second output.mat','FO_F50','Nall','F50')
FO_F50_1s=FO_F50;
load('F50 vs Force 1 minute output.mat','FO_F50','Nall','F50')
figure
loglog(6:3000,FO_F50(6:3000),'--k','LineWidth',1)
loglog(6:3000,FO_F50_1s(6:3000),':k','LineWidth',1)
hold on
plot(Nall,F50,'o','LineWidth',1,'MarkerFaceColor','w')
xlabel('N')
ylabel('F50 (pN)')
axis([-inf inf -inf inf])
axis square
grid on
set(gca,'FontSize',14,'FontName','Arial','XTick',[10 100 1000],...
    'XTickLabel',{'10','100','1,000'},'YTick',[1 10 100],'YTickLabel',...
    {'1','10','100'})
hold on