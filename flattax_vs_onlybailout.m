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
% Preference parameters:
    S.rho       = 0.04;         % Time discount rate.
    S.sigma     = 0.1;    
    S.sigma2    = S.sigma^2;
    S.Cwp       = S.rho;
    S.kappa     = 0.7;
    S.alpha     = 1;
%% Tax schedule
    S.policyrules=@linearrules_bailout; % linear taxing with bailout
    
%% Cases
% Possible thats (second column zeros)
    that=[0.01:0.001:0.035];K=length(that);
    that=[that; zeros(1,K)]';
% Possible equivalent bailouts (first column zeros)    
    equivalent_bail=@(t) -(S.rho./S.alpha).*((1-(1-S.kappa).*S.alpha.*t/2./S.rho).^(1/(1-S.kappa))-1)./S.Cwp;    
    bail=equivalent_bail(that(:,1));
    bail=[zeros(K,1) bail];
    L=2;

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
   parfor k=1:K
        
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
            V075(k,l)=SSS.V(ceil(3.*SSS.N/4));
            bailouted(k,l)=SSS.bailouted;
            if print
                fprintf('fracWP0= %4.2f\n',fracWP0(k,l))
            end
       end
       fprintf('\b|\n');
    end

   diary([pwd '/figures/mesh/output.txt'])
    diary off
%% Plot C,V,LoM and tax schedules for flat tax vs bailout
pos=20;% Choose with position to plot

Pvec=linspace(0,1,N)';
figure;set(gcf,'units','normalized','position',[0.01,0.25,0.65,0.65])%set(gcf,'units','points','position',[10,50,800,650])
    % 1. Value functions:
        subplot(2,2,1)
        
        % Planner's solutions are hard to compare since they are non-stochastic...
        %rhoVmuP = vbar + eta*log( Pvec );   % rho*V^{mu=P}: Value for Region 1 from 
                                            % the planner's solution under mu=P.
        %rhoV2muP= vbar + eta*log(1-Pvec);   % rho*V^{',mu=P}: Value for Region 2 
                                            % that same planner's solution.
        plot(Pvec,squeeze(VAL1(pos,1,:)),'-b', Pvec,squeeze(VAL1(pos,2,:)),'--r',...     % Multiply value functions by
             Pvec,VALWP(pos,1)*ones(N,1),'-.k'           )    % rho to get them to per-period 
        %     Pvec,rhoVmuP,'-.k',Pvec,rhoV2muP,'-.k')   % level, as the others.
        xlabel('P'), ylabel('V^1, V^2'), title('Value functions')
        yticks([])
        text(1.05,VALWP(pos,1),'V^{wp}')
        legend('Flat tax','Only Bailout')
    % 2. Drift (and volatility, if varying) of P:
        subplot(2,2,2),
        plot(Pvec,squeeze(PLAW(pos,1,:)),'-b',Pvec,squeeze(PLAW(pos,2,:)),'--r',Pvec,0*Pvec,'-.k')
            %hold on, plot(S.Pvec,S.svec,'-.k') 
            %legend('drift','vol.')
            ylabel('a(P)')
        xlabel('P'), title('P: Law of motion')
    % 3. Consumption functions:
        subplot(2,2,3), plot(Pvec,squeeze(CONS1(pos,1,:)),'-b', Pvec,squeeze(CONS1(pos,2,:)),'--r', Pvec,S.Cwp*ones(N,1),'-.k')
        xlabel('P'), ylabel('C^1'), title('Consumption')
        text(1.05,S.Cwp,'C^{wp}')
    % 4. Government policy rules:
        subplot(2,2,4), 
        plot(Pvec,squeeze(TR1(pos,1,:)),'-b', Pvec,squeeze(TR1(pos,2,:)),'--r')
        title('Government policy')
%         text(1.05, TR1(1,1,N),'T^1(P)')
%         text(1.05, TR1(2,2,N),'T^2(P)')

        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        
        st1=['that=' char(string(that(1,1))) ',fracWP0=' char(string(fracWP0(1,1))) ')'];
        namefig0=['Parameters:' st1];
        text(0.45, 0.98,namefig0)
      saveas(gcf,[pwd '/figures/unidim/LoM_flattax_vs_onlybailout'])
      saveas(gcf,[pwd '/figures/unidim/LoM_flattax_vs_onlybailout.png'])
   
%% V(0.5) vs that
figure
hold on
plot(that(:,1),V05(:,1),'-','LineWidth',2)
plot(that(:,1),V05(:,2),'--','LineWidth',2)

xlabel('$\hat{t}$','Interpreter','Latex','FontSize',16), ylabel('V(0.5)','FontSize',16), title('Bailouts under commitment','FontSize',16)%xlabel('that'), ylabel('g2'), title('V05')
legend({'Flat tax','Only bailout'},'Interpreter','Latex','FontSize',14,'Location','southeast')
   saveas(gcf,[pwd '/figures/unidim/V05_flattax_vs_onlybailout'])
   saveas(gcf,[pwd '/figures/unidim/V05_flattax_vs_onlybailout.png'])
   
%% V(0.75) vs that
figure
hold on
plot(that(:,1),V075(:,1),'-','LineWidth',2)
plot(that(:,1),V075(:,2),'--','LineWidth',2)

xlabel('$\hat{t}$','Interpreter','Latex','FontSize',16), ylabel('V(0.75)','FontSize',16), title('Bailouts under commitment','FontSize',16)%xlabel('that'), ylabel('g2'), title('V05')
legend({'Flat tax','Only bailout'},'Interpreter','Latex','FontSize',14,'Location','southeast')
   saveas(gcf,[pwd '/figures/unidim/V075_flattax_vs_onlybailout'])
   saveas(gcf,[pwd '/figures/unidim/V075_flattax_vs_onlybailout.png'])
   
%% EV vs that
% Mean value
EV1=mean(VAL1,3);
figure
hold on
plot(that(:,1),EV1(:,1),'-','LineWidth',2)
plot(that(:,1),EV1(:,2),'--','LineWidth',2)
xlabel('$\hat{t}$','Interpreter','Latex','FontSize',16), ylabel('EV','FontSize',16), title('Bailouts under commitment','FontSize',16)%xlabel('that'), ylabel('g2'), title('V05')
legend({'Flat tax','Only bailout'},'Interpreter','Latex','FontSize',14,'Location','southeast')
      saveas(gcf,[pwd '/figures/unidim/EV_flattax_vs_onlybailout'])
   saveas(gcf,[pwd '/figures/unidim/EV_flattax_vs_onlybailout.png'])
   