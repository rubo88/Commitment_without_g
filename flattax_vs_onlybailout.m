% This code compares the commitment eq of a flat tax with the one of only
% bailouts(zero tax)

clear all
clc
close all
diary on
addpath('..','funciones');  
%%
%%
    S.options.ploteach=0;
    S.options.print=0;print=S.options.print;
%% Parameters
    S.N=101;N=S.N;
    S.constVol= false; % If true: s(P)=sigma is constant.% If false: s(P)= 4*P*(1-P)*sigma
    S.policyrules=@linearrules_bailout;
% Preference parameters:
    S.rho       = 0.04;         % Time discount rate.
    S.sigma     = 0.1;    
    S.sigma2    = S.sigma^2;
    S.Cwp       = S.rho;
    S.kappa     = 0.7;
    S.alpha     = 1;
    
%% Cases
    bail0=-(S.rho./S.alpha).*((1-(1-S.kappa).*S.alpha.*0.025/2./S.rho).^(1/(1-S.kappa))-1);
    that         =[0.025 0];K=length(that);
    bail         =[0 bail0/S.Cwp];L=length(bail);
    bail         =repmat(bail,[K,1]);
    that         =repmat(that,[L,1])';

%% Allocate variables to store alocations
    CONS1=nan(K,L,S.N);CONS2=nan(K,L,S.N);
    VAL1=nan(K,L,S.N);VAL2=nan(K,L,S.N);
    TR1=nan(K,L,S.N);TR2=nan(K,L,S.N);
    GOV=nan(K,L,S.N);
    
    VALWP=nan(K,L);CONSWP=nan(K,L);
    PLAW=nan(K,L,S.N);VOL=nan(K,L,S.N);
    fracWP0=nan(K,L);V05=nan(K,L);bailouted=nan(K,L);
%% Run commitment2.m for every case
      fprintf('Progress:\n');fprintf(['\n' repmat('.',1,K) '\n\n'])
   for k=1:K
        
       for l=1:L
            if print
                fprintf('**************************************************************************************************************\n')
                fprintf('Running with values that= %4.1f\n',[that(k,l) ])   
            end
            SS=S;
            SS.that=that(k,l);SS.bail=bail(k,l);
            SSS=commitment2(SS);
            CONS1(k,l,:)=SSS.C;CONS2(k,l,:)=SSS.C2;
            VAL1(k,l,:)=SSS.rho*SSS.V;VAL2(k,l,:)=SSS.rho*SSS.V2;
            TR1(k,l,:)=SSS.Tvec;TR2(k,l,:)=SSS.T2vec;
            PLAW(k,l,:)=SSS.a;VOL(k,l,:)=SSS.svec;
            VALWP(k,l)=SSS.rhoVwp;CONSWP(k,l)=SSS.Cwp;
            fracWP0(k,l)=SSS.fracWP0;
            V05(k,l)=SSS.V(ceil(SSS.N/2));
            bailouted(k,l)=SSS.bailouted;
            if print
                fprintf('fracWP0= %4.2f\n',fracWP0(k,l))
            end
       end
       fprintf('\b|\n');
    end

   diary([pwd '/figures/mesh/output.txt'])
    diary off
%% Plot results:
Pvec=SSS.Pvec;
figure;set(gcf,'units','normalized','position',[0.01,0.25,0.65,0.65])%set(gcf,'units','points','position',[10,50,800,650])
    % 1. Value functions:
        subplot(2,2,1)
        
        % Planner's solutions are hard to compare since they are non-stochastic...
        %rhoVmuP = vbar + eta*log( Pvec );   % rho*V^{mu=P}: Value for Region 1 from 
                                            % the planner's solution under mu=P.
        %rhoV2muP= vbar + eta*log(1-Pvec);   % rho*V^{',mu=P}: Value for Region 2 
                                            % that same planner's solution.
        plot(Pvec,squeeze(VAL1(1,1,:)),'-b', Pvec,squeeze(VAL1(2,2,:)),'--r',...     % Multiply value functions by
             Pvec,VALWP(1,1)*ones(N,1),'-.k'           )    % rho to get them to per-period 
        %     Pvec,rhoVmuP,'-.k',Pvec,rhoV2muP,'-.k')   % level, as the others.
        xlabel('P'), ylabel('V^1, V^2'), title('Value functions')
        yticks([])
        text(1.05,VALWP(1,1),'V^{wp}')
        legend('Flat tax','Only Bailout')
    % 2. Drift (and volatility, if varying) of P:
        subplot(2,2,2),
        plot(Pvec,squeeze(PLAW(1,1,:)),'-b',Pvec,squeeze(PLAW(2,2,:)),'--r',Pvec,0*Pvec,'-.k')
            %hold on, plot(S.Pvec,S.svec,'-.k') 
            %legend('drift','vol.')
            ylabel('a(P)')
        xlabel('P'), title('P: Law of motion')
    % 3. Consumption functions:
        subplot(2,2,3), plot(Pvec,squeeze(CONS1(1,1,:)),'-b', Pvec,squeeze(CONS1(2,2,:)),'--r', Pvec,SSS.Cwp*ones(N,1),'-.k')
        xlabel('P'), ylabel('C^1'), title('Consumption')
        text(1.05,SSS.Cwp,'C^{wp}')
    % 4. Government policy rules:
        subplot(2,2,4), 
        plot(Pvec,squeeze(TR1(1,1,:)),'-b', Pvec,squeeze(TR1(2,2,:)),'--r')
        title('Government policy')
%         text(1.05, TR1(1,1,N),'T^1(P)')
%         text(1.05, TR1(2,2,N),'T^2(P)')

        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        
        st1=['that=' char(string(that(1,1))) ',fracWP0=' char(string(fracWP0(1,1))) ')'];
        namefig0=['Parameters:' st1];
        text(0.45, 0.98,namefig0)

