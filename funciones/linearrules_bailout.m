function SS=linearrules_bailout(S)
    Nmid=floor(S.N/2+1);
    S.Tvec = zeros(S.N,1);S.T2vec = zeros(S.N,1);
    S.Tvec(Nmid:end,1) =  S.that*(S.Pvec(Nmid:end,1)-0.5) ;
    S.T2vec(1:Nmid,1)  = -S.that*(S.Pvec(1:Nmid,1)-0.5) ;
    
    S.Tvec(1:Nmid,1) =  -((1+S.T2vec(1:Nmid,1)).^(1-S.kappa)-1)./(1-S.kappa) ;
    S.T2vec(Nmid:end,1)  = -((1+S.Tvec(Nmid:end,1)).^(1-S.kappa)-1)./(1-S.kappa) ;   
    
    S.bailouted =0;
    if S.Tvec(1)>-S.bail*S.Cwp
        S.Tvec(1)   =-S.bail*S.Cwp;
        S.Tvec(end) = (S.rho./S.alpha).*((1+(1-S.kappa).*S.alpha.*S.Tvec(1)./S.rho).^(1/(1-S.kappa))-1);
        S.T2vec(end)=-S.bail*S.Cwp;
        S.T2vec(1)  = (S.rho./S.alpha).*((1+(1-S.kappa).*S.alpha.*S.Tvec(1)./S.rho).^(1/(1-S.kappa))-1);
        S.bailouted =1;
    end
    S.fracWP0=-S.Tvec(1)/S.Cwp;
    SS=S;
end