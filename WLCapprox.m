function [x]=WLCapprox(F,N,k,P,kBT)
    
    ExponentTerm=exp((900*kBT./(F*P)).^.25);
    
    Term1=4/3*(1-1./sqrt(F*P/kBT + 1));
    Term2=10*ExponentTerm./(sqrt(F*P/kBT).*(ExponentTerm-1).^2);
    Term3=(F*P/kBT).^1.62./(3.55 + 3.8*(F*P/kBT).^2.2);
    
    x=N*k*(Term1-Term2+Term3);
    
end