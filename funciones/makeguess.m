function SS=makeguess(S)
       % Consumption at the borders
        Call = S.rho*S.eta;         % Region 1's C(1) when they have all wealth.
        C2no= -S.T2vec(S.N);        % C'(1): Just eat transfers.
        C2all= S.rho*S.eta;         % Region 2's C'(0) when they have all wealth.
        Cno = - S.Tvec(1);        % C(0) : Just eat transfers.
    % Make linear guess for value functions.    
        S.Cmat(:,1)  = (1-S.Pvec).*Cno   + S.Pvec.*Call;    
        S.C2mat(:,1) = (1-S.Pvec).*C2all + S.Pvec.*C2no;
    % Get Pdot for our guess:
        aguess = Pdotfct( S );
    % Flow rates
        [fUpGs, fDnGs, fLvGs] = FlowRates(aguess,svec,dP);

        Fmat =   sparse(1:N-1,2:N  ,fUpGs(1:N-1),N,N) ...
               + sparse(2:N  ,1:N-1,fDnGs(2:N  ),N,N) ...
               - sparse(1:N  ,1:N  ,fLvGs       ,N,N)    ;
    % Flow utility of both agents:
        [U,U2] = flowutility(S);
    
    % Values
        FF = S.rho*speye(S.N) - Fmat;
        S.Vmat(:,1) = FF \  U(S);
        S.V2mat(:,1)= FF \ U2(S.Cmat(:,1),C2guess);
end