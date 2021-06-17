clear;close all;clc;

%% Plot koff vs Force
figure
F=logspace(-1,3,1000);
loglog(F,6.86*10^7*exp(-14*(2.53-(0.15)*F/4.114)),'LineWidth',1)
hold on
loglog(F,6.86*10^7*exp(-14*(2.53-(0.15)*F/400/4.114)),'LineWidth',1)
grid on
set(gca,'FontSize',14,'FontName','Arial','XTick',[.1 1 10 100 1000])
grid on
axis square
xlabel('Force magnitude (pN)')
ylabel('Rupture Rate (s-1)')
axis([0 inf 10^-8 10^-3])

%% Initialize PDF calculations
F=0.01:.01:60;
kBT=4.114;
Npts=121;
load('FEX15bpAndTethTst.mat','F','xTGT','xTethTst','Gtgt','Gtst','xTeth','Gteth')

FOg=fit(xTeth',Gteth','linearinterp');
FOgOn=fit(xTethTst',Gtst','linearinterp');

Width=31;
x=linspace(-Width,Width,Npts*2+1);
y=x;
z=linspace(0,Width,Npts+1);
[x,y,z]=meshgrid(x,y,z);
r=sqrt(x.^2+y.^2+z.^2);

[xf,yf]=meshgrid(linspace(0,Width*2,Npts*2));
xf=xf(:);
yf=yf(:);
kOnAll=x(:);
A=max(factor(numel(x)))
% A=61*11*11;
for i=1:(numel(x)/A)
    dx=xf-x(i*A-(0:(A-1)));
    dy=yf-y(i*A-(0:(A-1)));
    rAll=sqrt(dx.^2+dy.^2);
    if i==1
        kOn=rAll;
    end
    kOn(1:numel(kOn))=exp(-FOgOn(rAll(:))/kBT);
    kOnAll(i*A-(0:(A-1)))=sum(kOn,1);
        i/(numel(x)/A)
end
clear dx dy rAll kOn

%% Calculate energy landscapes and PDF
ForceField=[0 logspace(-2,2.5,50)];
for j=1:length(ForceField)

    En{1}=zeros(size(x));
    En{1}(1:numel(En{1}))=FOg(r(:))/kBT+ForceField(j)*x(:)/kBT;
    PDF{1}=exp(-En{1});
    PDF{1}=PDF{1}/sum(PDF{1},'all');
    AvgDisp{1}(j)=sum(x.*PDF{1},'all');
    kOnCum(j,1)=sum(PDF{1}(:).*kOnAll);
    
    
    En{2}=En{1};
    En{2}(1:numel(En{2}))=400*FOg(r(:))/kBT+ForceField(j)*x(:)/kBT;
    PDF{2}=exp(-En{2});
    PDF{2}=PDF{2}/sum(PDF{2},'all');
    AvgDisp{2}(j)=sum(x.*PDF{2},'all');
    kOnCum(j,2)=sum(PDF{2}(:).*kOnAll);
    
    j/length(ForceField)
end

figure
subplot(1,2,1)
scale=(10^4/60)/kOnCum(1);
kOnCum=kOnCum*scale;
semilogx(ForceField,kOnCum,'o','LineWidth',1,'MarkerFaceColor','w')
axis([.01 315 -inf inf])
xlabel('Force (pN)')
ylabel('On Rate (s-1)')
set(gca,'FontSize',14,'FontName','Arial','XDir','reverse')
axis([.01 315 0 350])
axis square
grid on

subplot(1,2,2)
semilogx(ForceField,-AvgDisp{1},'o','LineWidth',1,'MarkerFaceColor','w')
hold on
semilogx(ForceField,-AvgDisp{2},'o','LineWidth',1,'MarkerFaceColor','w')
xlabel('Force (pN)')
ylabel('Mean Displacement')
legend('1 Tether','400 Tethers')
grid on
axis square
axis([-inf inf 0 25])
set(gca,'FontSize',14,'FontName','Arial')


% ForceField=[0 -logspace(-2,2.5,50)];
% for j=1:length(ForceField)
% 
%     En{1}=zeros(size(x));
%     En{1}(1:numel(En{1}))=FOg(r(:))/kBT+ForceField(j)*x(:)/kBT;
%     PDF{1}=exp(-En{1});
%     PDF{1}=PDF{1}/sum(PDF{1},'all');
%     kOnCum(j,1)=sum(PDF{1}(:).*kOnAll);
%     
%     
%     En{2}=En{1};
%     En{2}(1:numel(En{2}))=400*FOg(r(:))/kBT+ForceField(j)*x(:)/kBT;
%     PDF{2}=exp(-En{2});
%     PDF{2}=PDF{2}/sum(PDF{2},'all');
%     kOnCum(j,2)=sum(PDF{2}(:).*kOnAll);
%     
%     j/length(ForceField)
% end
% 
% subplot(1,2,2)
% % kOnCum=kOnCum;
% semilogx(-ForceField,kOnCum,'o','LineWidth',1,'MarkerFaceColor','w')
% axis([.01 315 -inf inf])
% xlabel('Force (pN)')
% ylabel('On Rate (s-1)')
% set(gca,'FontSize',14,'FontName','Arial')
% axis([.01 315 0 350])
% axis square
% grid on