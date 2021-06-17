clear;close all;clc;
N=2750;
% N=300;
n=0:N;
kdeg=25;
kon=[logspace(-1,5,21)];
kon(end)=[];
tmax=60;
Init=zeros(size(n));
Init(end)=1;

for i=1:length(kon)
    dt=.01/max(N*[kdeg kon(i)]);
    MM=eye(N+1);
    for j=2:N
        MM(j,j-1)=dt*kdeg*(j-1);
        MM(j,j+1)=dt*kon(i)*(N-j+1);
        MM(j,j)=1-MM(j,j+1)-MM(j,j-1);
    end
    MM(end,end-1)=dt*kdeg*N;
    MM(end,end)=1-MM(end,end-1);
    
    NumStep=round(tmax/dt);
    Frac=Init*MM^NumStep;
    
    FracCum=cumsum([0 Frac(2:end)]);
    [~,Ind]=unique(FracCum);
    FracCum=FracCum(Ind);
    nInterp=n(Ind);
    if length(Ind)>1
        FracCum=FracCum/sum(Frac(2:end));
        [~,Ind]=unique(FracCum);
        FracCum=FracCum(Ind);
        nInterp=n(Ind);
        LB(i)=interp1(FracCum,nInterp,[.05]);
        UB(i)=interp1(FracCum,nInterp,[.95]);
        MED(i)=interp1(FracCum,nInterp,[.5]);
    else
        LB(i)=0;
        UB(i)=0;
        MED(i)=0;
    end
    
    Bound(i)=1-Frac(1);
    AvgPoly(i)=sum(Frac(2:end).*n(2:end))/sum(Frac(2:end));
    
    i/length(kon)
    
end
figure
semilogx(kon,AvgPoly)
xlabel('k_{on} (s^{-1})')
ylabel('Steady-state polyvalency')
axis([-inf inf -inf inf])