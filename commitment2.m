% Solves the commitment model by iterating backwards in time using the HJBs
% (seeing them as PDEs in P and t).
%clear all;clc;close all
function SS=commitment2(S)
N=S.N;
%% Options
    S.options.makePlot = false;   % Make 3D plots on the way or not
    S.options.PlotEvery =1;     % How often to make these plots
    S.options.eps = 10^(-7);      % Convergence criterion (on log-consumption difference,
                        % in yearly terms)
    S.options.T  = 10000;         % Maximal number of time iterations.                    
    T=S.options.T ;  
                  
%% Set up gridfor P
    S.dP = 1/(N-1);
    S.Pvec  = linspace(0,1,N)';   % Make grid for P: N-by-1 vector.
    
%     if S.constVol                 % To set up vector with volatility s(P):
%         S.svec = S.sigma*ones(N,1);         % Suppress Matlab's message.
%     else                          % Make it constant or put quadratic form.
        S.svec = sqrt(S.Pvec.*(1-S.Pvec))*S.sigma;  %#ok<*UNRCH>
        volPenaltyVec =  S.sigma^2/(2*S.rho).*ones(N,1);
%     end

% Function that takes the centered derivative:
    CentDiff = @(x) [ ( x(2    ) - x(1      ) )/    S.dP; ...
                      ( x(3:end) - x(1:end-2) )/(2*S.dP); ...
                      ( x(N    ) - x(N-1    ) )/    S.dP      ];
                  
%% Policy
S=feval(S.policyrules,S);
%     load T_smooth_symm.mat
%     S.Tvec=Tsymm;S.T2vec=-Tsymm;
%     load T_nocomm.mat
%     S.Tvec=T1;S.T2vec=T2;
%     S.fracWP0=-S.Tvec(1)/S.Cwp;
%     clear T1 T2
%% Start to construct guess for value function: get initial guess for 
    % Consumption at the borders
        Call = S.rho;         % Region 1's C(1) when they have all wealth.
        C2no= -S.T2vec(N);        % C'(1): Just eat transfers.
        C2all= S.rho;         % Region 2's C'(0) when they have all wealth.
        Cno = - S.Tvec(1);        % C(0) : Just eat transfers.
    % Make linear guess for value functions.    
        Cguess = (1-S.Pvec).*Cno   + S.Pvec.*Call;    
        C2guess= (1-S.Pvec).*C2all + S.Pvec.*C2no;
    % Get Pdot for our guess:
        Pdotfct = @(C,Cp,T1,T2) -(1-S.Pvec).*(C+T1) + S.Pvec.*(Cp+T2);
        %Pdot = Pdotfct( Cguess, C2guess, Pvec );
        aguess = Pdotfct( Cguess, C2guess, S.Tvec, S.T2vec );
    % Flow rates
        [fUpGs, fDnGs, fLvGs] = FlowRates(aguess,S.svec,S.dP);

        Fmat =   sparse(1:N-1,2:N  ,fUpGs(1:N-1),N,N) ...
               + sparse(2:N  ,1:N-1,fDnGs(2:N  ),N,N) ...
               - sparse(1:N  ,1:N  ,fLvGs       ,N,N)    ;
    % Flow utility of both agents:
        U = @(C,C2) log(C )  - (C+C2+S.Tvec+S.T2vec)/S.rho;
        U2= @(C,C2) log(C2)  - (C+C2+S.Tvec+S.T2vec)/S.rho;
    % Values
        FF = S.rho*speye(N) - Fmat;
        V = FF \(  U(Cguess,C2guess)- volPenaltyVec);
        V2= FF \ (U2(Cguess,C2guess)- volPenaltyVec);
   
%% Allocate and store guesses
    % Allocate
        Vmat = nan(N,T);  V2mat = nan(N,T);     % Set up storage for value and
        Cmat = nan(N,T);  C2mat = nan(N,T);     % consumption functions
        tvec = nan(1,T);                        % Time vector
        LnDCvec = nan(1,T);
    % Store guesses.
        Vmat(:,1) = V     ;   V2mat(:,1) = V2     ; 
        Cmat(:,1) = Cguess;   C2mat(:,1) = C2guess;
        tvec(  1) = 0;
    % Prepare everything for the while to start
        oldLnCC  = 1;
        LnDiffCC = 1;
        i = 2;
        abort = false;
    % Prepare command window
        if S.options.print
            disp('........................................................')
            fprintf('%1.0f grid points, that=%4.2f, sigma=%5.3f\n', [N, S.that, S.sigma])
        end
    %tic
    if S.options.makePlot, figure; end         

%% MAIN LOOP
    
    while ~abort && LnDiffCC>S.options.eps && i<=T
        % Differentiate values
            VP    = CentDiff(V );
            V2P   = CentDiff(V2);
        % Compute consumptions
            C     = 1./ ( 1/S.rho + (1-S.Pvec).*VP  );
            C(1)  = min(  C(1), - S.Tvec(1) );
            C2    = 1 ./ ( 1/S.rho -    S.Pvec .*V2P );
            C2(N) = min( C2(N), -S.T2vec(N) );
            if any( [C; C2] < 0 );abort = true; end
            if  ~isreal([C; C2]) ;abort = true; end
        % Get Pdot
            aa    = Pdotfct(C,C2, S.Tvec, S.T2vec); % -(1-Pvec).*(C+Tvec) + P.*(C2+T2vec);    % get drift
        % Flowrates
            [ffUp, ffDn, ffLv] = FlowRates( aa, S.svec, S.dP);
            dt  =1/max(ffLv);       % Get maximal time increment that still has all stable probabilities (but don't go beyond one year).
        % Continuation values
            Vtom  = dt*ffUp.*[ V(2:N);  0 ] + ...
                    dt*ffDn.*[ 0; V(1:N-1)] + ...
                   (1-dt*ffLv).* V;   % Get continuation value.
            V2tom = dt*ffUp.*[ V2(2:N);  0 ] + ...
                    dt*ffDn.*[ 0; V2(1:N-1)] + ...
                    (1-dt*ffLv).*V2;    
        % New values
            V =  U(C,C2)*dt + (1-S.rho*dt).* Vtom - volPenaltyVec;
            V2= U2(C,C2)*dt + (1-S.rho*dt).*V2tom - volPenaltyVec;
        % Store
            Vmat(:,i) = V;   V2mat(:,i) = V2;
            Cmat(:,i) = C;   C2mat(:,i) = C2;
            tvec(  i) = tvec(i-1)-dt;
        % ..
            newLnCC  = log( [C; C2] );
            LnDiffCC = max(abs( (newLnCC - oldLnCC)/dt ));
            LnDCvec(i) = LnDiffCC;

        % Plot stuff if nedeed
            if S.options.makePlot && mod(i,S.options.PlotEvery)==0
                subplot(1,3,1), surf(tvec,S.Pvec,Vmat);  title('Value')       
                xlabel('t'), ylabel('P'), zlabel('V')
                subplot(1,3,2), surf(tvec,S.Pvec,Cmat);  title('Consumption')
                xlabel('t'), ylabel('P'), zlabel('C')
                subplot(1,3,3), plot(tvec,LnDCvec); title('max.abs. cons. growth rate')
                xlabel('t'), ylabel('|Cdot|/C'), 
                lastdiff = max( LnDCvec(i-S.options.PlotEvery+1:i) );
                %ylim([0, lastdiff*3]);
                getframe;
            end
        % Update loop variables
            oldLnCC = newLnCC;              
            i = i+1;
    end
    %toc 

%% Display messages    
    if abort
        if S.options.print
            disp('NEGATIVE VALUE OF SAVINGS!')
            fprintf('Value-function iteration aborted after %1.0f iterations (%4.1f years)\n', ...
                    [i-1, abs(tvec(i-1))] );
        end
        %V = Vmat(:,i-2);  V2 = V2mat(:,i-2);  C = Cmat(:,i-2);  C2 = C2mat(:,i-2);
        S.V = NaN(N,1);  S.V2 = NaN(N,1);  S.C =NaN(N,1);  S.C2 = NaN(N,1);
    else
        if S.options.print
            fprintf('Value-function iteration converged after %1.0f iterations (%4.1f years)\n', ...
                    [i-1, abs(tvec(i-1))] );
            fprintf('Terminal log-consumption distance (yearly): %10.8f%%\n', LnDiffCC ); 
        end
        S.V = Vmat(:,i-1);  S.V2 = V2mat(:,i-1);  S.C = Cmat(:,i-1);  S.C2 = C2mat(:,i-1);
    end

%% Plot results:
        S.a = Pdotfct(C,C2,S.Tvec, S.T2vec );
        S.vbar= log(S.rho); % Constant Vbar from the paper.
        S.rhoVwp =  S.vbar - 1 ;          % rho*V^{wp}: value from wealth-pooling
        %S.C=C;S.C2=C2;
% if S.options.ploteach==1
%     figure;set(gcf,'units','normalized','position',[0.01,0.25,0.65,0.65])%set(gcf,'units','points','position',[10,50,800,650])
%     % 1. Value functions:
%         subplot(2,2,1)
%         
%         % Planner's solutions are hard to compare since they are non-stochastic...
%         %rhoVmuP = vbar + eta*log( Pvec );   % rho*V^{mu=P}: Value for Region 1 from 
%                                             % the planner's solution under mu=P.
%         %rhoV2muP= vbar + eta*log(1-Pvec);   % rho*V^{',mu=P}: Value for Region 2 
%                                             % that same planner's solution.
%         plot(S.Pvec,S.rho*S.V,'-b', S.Pvec,S.rho*S.V2,'--r',...     % Multiply value functions by
%              S.Pvec,S.rhoVwp*ones(N,1),'-.k'           )    % rho to get them to per-period 
%         %     Pvec,rhoVmuP,'-.k',Pvec,rhoV2muP,'-.k')   % level, as the others.
%         xlabel('P'), ylabel('V^1, V^2'), title('Value functions')
%         yticks([])
%         text(1.05,S.rhoVwp,'V^{wp}')
%     % 2. Drift (and volatility, if varying) of P:
%         subplot(2,2,2),
%         plot(S.Pvec,S.a,'--k',S.Pvec,0*S.Pvec,'-k')
%         if S.constVol
%             ylabel('drift a(P)')
%             text(0.5,0.7*max(S.a),sprintf('s(P)=%1.3f',S.sigma))
%         else
%             hold on, plot(S.Pvec,S.svec,'-.k') 
%             legend('drift','vol.')
%             ylabel('a(P),s(P)')
%         end
%         xlabel('P'), title('P: Law of motion')
%     % 3. Consumption functions:
%         subplot(2,2,3), plot(S.Pvec,S.C,'-b', S.Pvec,S.C2,'--r', S.Pvec,S.Cwp*ones(N,1),'-.k')
%         xlabel('P'), ylabel('C^1, C^2'), title('Consumption')
%         text(1.05,S.Cwp,'C^{wp}')
%     % 4. Government policy rules:
%         subplot(2,2,4), 
%         plot(S.Pvec,S.Gvec,'-.k',S.Pvec,S.Tvec,'-b',S.Pvec,S.T2vec,'--r')
%         title('Government policy')
%         text(1.05, S.Gvec(N),  'G(P)')
%         text(1.05, S.Tvec(N),'T^1(P)')
%         text(1.05,S.T2vec(N),'T^2(P)')
% 
%         ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%         
%         st1=['(g0=' char(string(S.g0)) ',g2=' char(string(S.g2))  ',that=' char(string(S.that)) ',fracWP0=' char(string(S.fracWP0)) ')'];
%         namefig0=['Parameters:' st1];
%         text(0.45, 0.98,namefig0)
% end

SS=S;
end
