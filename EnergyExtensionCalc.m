% clear;close all;clc;

%% Initialize
F=0.01:.01:40;
kBT=4.114;
Npts=60;
SampleInt=4;

%% Create Force-extension curve
                xAnchor = WLCapprox(F,15,.34,53,kBT);
         xSurfaceSpacer = WLCapprox(F,15,.6,1.07,kBT);
        xParticleSpacer = WLCapprox(F,15,.6,1.07,kBT);
         xPairingRegion = WLCapprox(F,15,.33,100,kBT);
      xPairingRegionTst = WLCapprox(F,15,.53,2.7,kBT);
  xInterDuplexSpacerRNA = WLCapprox(F,5,.75,.67,kBT);
  
      xTeth = xSurfaceSpacer + xAnchor + xInterDuplexSpacerRNA + ...
                xPairingRegion + xParticleSpacer;
      xTst = xSurfaceSpacer + xAnchor + xInterDuplexSpacerRNA + ...
                xPairingRegionTst + xParticleSpacer;

F=[0 F];
xTeth=[0 xTeth];
xTst=[0 xTst];
Gteth=cumtrapz(xTeth,F);
Gtst=cumtrapz(xTst,F);


%% Plot force-extension curve
FOg=fit(xTeth',Gteth','linearinterp');
FOgTst=fit(xTeth',Gtst','linearinterp');
delta=linspace(0,max(xTeth),1000);
Gint=FOg(delta)';
figure('OuterPosition',[481.0000  473.0000  839.2000  371.2000])
subplot(1,2,1)
plot(delta(1:20:end),FOg(delta(1:20:end)),'k-','LineWidth',1.5)
axis([0 inf 0 inf])
xlabel('End-to-end extension (nm)')
ylabel('Gtether (pN nm)')
set(gca,'FontSize',14,'FontName','Arial','XTick',0:5:25)
grid on

subplot(1,2,1)
plot(delta(1:20:end),FOgTst(delta(1:20:end)),'k-','LineWidth',1.5)
axis([0 inf 0 inf])
xlabel('Foot-fuel Distance (nm)')
ylabel('Gtst (pN nm)')
set(gca,'FontSize',14,'FontName','Arial','XTick',0:5:25)
grid on