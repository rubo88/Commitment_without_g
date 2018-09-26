%% This code plots V(0.5) over one paramter
%% Experiments

clear all
clc
close all
diary on
addpath('..','funciones');  
%% xaxis
%   =1  that
%   =2  g2
%   =3  bail
xaxis=1;
%%
    S.options.ploteach=1;
    S.options.print=0;print=S.options.print;
%% Parameters
    S.N=1001;N=S.N;
    S.constVol= false; % If true: s(P)=sigma is constant.% If false: s(P)= 4*P*(1-P)*sigma
    S.policyrules=@linearrules;
    S.X=@cara;
    % Preference parameters:
        S.rho = 0.04;         % Time discount rate.
        S.sigma   = S.rho/2;      % Strength of Brownian noise dP_t = ... + s(P_t) * dB_t
        S.sigma2  = S.sigma^2;
        S.Cwp = S.rho;
        S.kappa   = 0;%4/(S.rho);
%% Cases
L=1;
    if xaxis==1
        that=[0.01];K=length(that);
        bail=0;bail=repmat(bail,[K,1]);
    elseif xaxis==2
%         g2=[-0.5:0.01:0.5];K=length(g2);
%         that=2.7;that=repmat(that,[K,1]);
%         bail=0;bail=repmat(bail,[K,1]);
%         g0           =1;
    elseif xaxis==3
        bail=[0:0.001:0.5];K=length(bail);
        that=2.7;that=repmat(that,[K,1]);
        g2  =0;g2=repmat(g2,[K,1]);
        g0           =1;
    end

%% Allocate variables to store alocations
    CONS1=nan(K,L,S.N);CONS2=nan(K,L,S.N);
    VAL1=nan(K,L,S.N);VAL2=nan(K,L,S.N);
    TR1=nan(K,L,S.N);TR2=nan(K,L,S.N);
    
    VALWP=nan(K,L);CONSWP=nan(K,L);
    PLAW=nan(K,L,S.N);VOL=nan(K,L,S.N);
    fracWP0=nan(K,L);V05=nan(K,L);bailouted=nan(K,L);
%% Run commitment2.m for every case
      l=1;
      fprintf('Progress:\n');fprintf(['\n' repmat('.',1,K) '\n\n'])
   for k=1:K
            if print
                fprintf('**************************************************************************************************************\n')
                fprintf('Running with values that= %4.1f\n',that(k,l))   
            end
            SS=S;
            SS.that=that(k);SS.bail=bail(k);
            SSS=commitment2(SS);
            CONS1(k,l,:)=SSS.C;CONS2(k,l,:)=SSS.C2;
            VAL1(k,l,:)=SSS.rho*SSS.V;VAL2(k,l,:)=SSS.rho*SSS.V2;
            TR1(k,l,:)=SSS.Tvec;TR2(k,l,:)=SSS.T2vec;
            PLAW(k,l,:)=SSS.a;VOL(k,l,:)=SSS.svec;
            VALWP(k,l)=SSS.rhoVwp;CONSWP(k,l)=SSS.Cwp;
            fracWP0(k,l)=SSS.fracWP0;
            V05(k,l)=SSS.V(ceil(SSS.N/2));
            %bailouted(k,1)=SSS.bailouted;
            if print
                fprintf('fracWP0= %4.2f\n',fracWP0(k,l))
            end
       fprintf('\b|\n');
    end

   diary([pwd '/figures/mesh/output.txt'])
    diary off
%% Plot results:
figure;set(gcf,'units','normalized','position',[0.01,0.25,0.65,0.65])
if xaxis==1
   plot(that,V05);    
   xlabel('that'), ylabel('V05'), title('V05')
   saveas(gcf,[pwd '/figures/unidim/V05_vs_that'])
   saveas(gcf,[pwd '/figures/unidim/V05_vs_that.png'])
elseif xaxis==2
%    plot(g2,V05);    
%    xlabel('g2'), ylabel('V05'), title('V05')
%    saveas(gcf,[pwd '/figures/unidim/V05_vs_g2'])
%    saveas(gcf,[pwd '/figures/unidim/V05_vs_g2.png'])
elseif xaxis==3
   plot(bail,V05);    
   xlabel('bail'), ylabel('V05'), title('V05')
   saveas(gcf,[pwd '/figures/unidim/V05_vs_bail2'])
   saveas(gcf,[pwd '/figures/unidim/V05_vs_bail2.png'])
   save('data/unidim/bail1.mat')
end
