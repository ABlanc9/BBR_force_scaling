clear;close all;clc;
Nall=[15:3:30 40 60 100 200 350 500 800 1100 1550 2000];
cmap=jet(length(Nall));
kdeg=25;
kon=18.4;
tmax=1;
load('Fixed Force Dependant On Rate Output.mat','kOnCumUpper','ForceFieldUpper')
F=[ForceFieldUpper];
kon=kOnCumUpper*kon/max(kOnCumUpper(:));

%% Force against movement direction
figure(1)
for k=1:length(Nall)
    N=Nall(k);
    if k==1
        Find{k}=2:length(F);
    else
        Cutoff=max([1 Find{k-1}(find(Bound{k-1}<.99,1))]);
        Find{k}=Cutoff:length(F);
    end
    Find{k}=[1 Find{k}];
    
    for i=1:length(Find{k})
        n=0:N;
        Init=zeros(size(n));
        Init(end)=1;
        krup=6.86*10^7*exp(-14*(2.53-(0.15)*F(Find{k}(i))./n/4.114));
        koff(1)=0;
        koff=krup+kdeg;
        koff(1)=kdeg;
        dt=.01/max(N*[koff kon(Find{k}(i),:)]);
        MM=eye(N+1);
        for j=2:max(N)
            MM(j,j-1)=dt*koff(j)*(j-1);
            MM(j,j+1)=dt*kon(Find{k}(i),j)*(N-j+1);
            MM(j,j)=1-MM(j,j+1)-MM(j,j-1);
        end
        MM(end,end-1)=dt*kdeg*N;
        MM(end,end)=1-MM(end,end-1);

        NumStep=round(tmax/dt);
        Frac=Init*MM^NumStep;

        Bound{k}(i)=1-Frac(1);
        Bound{k}(isnan(Bound{k}))=0;
        Bound{k}(Bound{k}<0)=0;
        if i==1
            AvgPoly(k)=sum(Frac(2:end).*n(2:end))/sum(Frac(2:end));
        end

        [i/length(Find{k}) k/length(Nall)]
        if Bound{k}(i)<.0001
            Find{k}((i+1):end)=[];
            break
        end
    end
    
    DiscInd=[1 find(Bound{k}>.999999)];
    if length(DiscInd)>1
        DiscInd(end)=[];
    end
    
    Bound{k}(DiscInd)=[];
    Find{k}(DiscInd)=[];
    FO{k}=fit(Bound{k}',F(Find{k}),'linearinterp');
    F50(k)=FO{k}(0.5);
    semilogx(F(Find{k}),Bound{k},'s-','LineWidth',1,'MarkerFaceColor','w',...
        'MarkerEdgeColor',cmap(k,:),'Color',cmap(k,:))
    hold on
    drawnow
end

grid on
xlabel('Force (pN)')
ylabel(['Percent attached' newline 'after one minute'])
h=legend(strsplit(num2str(Nall)),'Location','eastoutside');
title(h,'N')
set(gca,'FontSize',14,'FontName','Arial','YTick',[0 .5 1])
axis([0 100 0 1])

figure(2)
subplot(3,2,1)
plot(Nall,F50,'o')
subplot(3,2,2)
plot(Nall,AvgPoly,'o')
subplot(3,1,[2 3])
FO_F50=fit(Nall(7:end)',F50(7:end)','poly1');
plot(1:(max(AvgPoly)+10),FO_F50(1:(max(AvgPoly)+10)),'--k')
hold on
plot(AvgPoly,F50,'o','LineWidth',1,'MarkerFaceColor','w')

% %% Force in movement direction
% load('Fixed Force Dependant On Rate Output Reverse.mat','kOnCumLower','ForceFieldLower')
% F=-flipud(ForceFieldLower);
% kon=fliplr(kOnCumLower*18.4/max(kOnCumUpper(:)));
% 
% figure(1)
% subplot(2,1,2)
% clear Bound Find
% tmax=1;
% Nall=[3:20];
% 
% for k=1:length(Nall)
%     N=Nall(k);
% %     if k==1
%         Find{k}=1:length(F);
% %     else
% %         Cutoff=max([1 Find{k-1}(find(Bound{k-1}<.99,1))]);
% %         Find{k}=Cutoff:length(F);
% %     end
%     
%     for i=1:length(Find{k})
%         n=0:N;
%         Init=zeros(size(n));
%         Init(end)=1;
%         krup=6.86*10^7*exp(-14*(2.53-(0.15)*F(Find{k}(i))./n/4.114));
%         koff(1)=0;
%         koff=krup+kdeg;
%         koff(1)=kdeg;
%         dt=.01/max(N*[koff kon(Find{k}(i),:)]);
%         MM=eye(N+1);
%         for j=2:max(N)
%             MM(j,j-1)=dt*koff(j)*(j-1);
%             MM(j,j+1)=dt*kon(Find{k}(i),j)*(N-j+1);
%             MM(j,j)=1-MM(j,j+1)-MM(j,j-1);
%         end
%         MM(end,end-1)=dt*kdeg*N;
%         MM(end,end)=1-MM(end,end-1);
% 
%         NumStep=round(tmax/dt);
%         Frac=Init*MM^NumStep;
% 
%         Bound{k}(i)=1-Frac(1);
%         Bound{k}(isnan(Bound{k}))=0;
%         Bound{k}(Bound{k}<0)=0;
%         if i==1
%             AvgPoly(k)=sum(Frac(2:end).*n(2:end))/sum(Frac(2:end));
%         end
% 
%         [i/length(Find{k}) k/length(Nall)]
%         if Bound{k}(i)<.0001
%             Find{k}((i+1):end)=[];
%             break
%         end
%     end
%     
% %     DiscInd=find(Bound{k}>.999999);
% %     if ~isempty(DiscInd)
% %         DiscInd(end)=[];
% %     end
%     
% %     Bound{k}(DiscInd)=[];
% %     Find{k}(DiscInd)=[];
% %     FO{k}=fit(Bound{k}',F(Find{k}),'linearinterp');
% %     F50(k)=FO{k}(0.5);
%     semilogx(F(Find{k}),Bound{k},'s-','LineWidth',1,'MarkerFaceColor','w',...
%         'MarkerEdgeColor',cmap(k,:),'Color',cmap(k,:))
%     hold on
%     drawnow
% end
% 
% 
% grid on
% xlabel('Force (pN)')
% ylabel(['Percent attached' newline 'after one minute'])
% h=legend(strsplit(num2str(Nall)),'Location','eastoutside');
% title(h,'N')
% set(gca,'FontSize',14,'FontName','Arial','YTick',[0 .5 1])
% axis([0 inf 0 1])
% 
% figure
% subplot(3,2,1)
% hold on
% plot(Nall,F50,'o')
% subplot(3,2,2)
% hold on
% plot(Nall,AvgPoly,'o')
% subplot(3,1,[2 3])
% FO_F50=fit(AvgPoly',F50','poly1');
% plot(1:(max(AvgPoly)+10),FO_F50(1:(max(AvgPoly)+10)),'--k')
% hold on
% plot(AvgPoly,F50,'o','LineWidth',1,'MarkerFaceColor','w')