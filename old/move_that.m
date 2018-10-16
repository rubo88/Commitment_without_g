%% This code plots V(0.5) over that
 
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
    S.options.ploteach=0;
    S.options.print=0;print=S.options.print;
%% Parameters
    S.N=101;N=S.N;
    S.constVol= false; % If true: s(P)=sigma is constant.% If false: s(P)= 4*P*(1-P)*sigma
    %S.policyrules=@linearrules;
    S.policyrules=@linearrules_bailout;
    %S.X=@cara;
% Preference parameters:
    S.rho       = 0.04;         % Time discount rate.
    S.sigma     = 0.4;    
    S.sigma2    = S.sigma^2;
    S.Cwp       = S.rho;
    S.kappa     = 0.7;
    S.alpha     = 1;
%% Cases

    if xaxis==1
        that=[0.01:0.001:0.05];K=length(that);
        bail=0;bail=repmat(bail,[K,1]);
    elseif xaxis==2
%         g2=[-0.5:0.01:0.5];K=length(g2);
%         that=2.7;that=repmat(that,[K,1]);
%         bail=0;bail=repmat(bail,[K,1]);
%         g0           =1;
    elseif xaxis==3
        bail=[0:0.01:0.5];K=length(bail);
        that=0.025;that=repmat(that,[K,1]);
    end

%% Allocate variables to store alocations
    CONS1=nan(K,S.N);CONS2=nan(K,S.N);
    VAL1=nan(K,S.N);VAL2=nan(K,S.N);
    TR1=nan(K,S.N);TR2=nan(K,S.N);
    
    VALWP=nan(K);CONSWP=nan(K);
    PLAW=nan(K,S.N);VOL=nan(K,S.N);
    fracWP0=nan(K);V05=nan(K);bailouted=nan(K);
%% Run commitment2.m for every case
      l=1;
      fprintf('Progress:\n');fprintf(['\n' repmat('.',1,K) '\n\n'])
   parfor k=1:K
            if print
                fprintf('**************************************************************************************************************\n')
                fprintf('Running with values that= %4.1f\n',that(k))   
            end
            SS=S;
            SS.that=that(k);SS.bail=bail(k);
            SSS=commitment2(SS);
            CONS1(k,:)=SSS.C;CONS2(k,:)=SSS.C2;
            VAL1(k,:)=SSS.rho*SSS.V;VAL2(k,:)=SSS.rho*SSS.V2;
            TR1(k,:)=SSS.Tvec;TR2(k,:)=SSS.T2vec;
            PLAW(k,:)=SSS.a;VOL(k,:)=SSS.svec;
            VALWP(k)=SSS.rhoVwp;CONSWP(k)=SSS.Cwp;
            fracWP0(k)=SSS.fracWP0;
            V05(k)=SSS.V(ceil(SSS.N/2));
            %bailouted(k,1)=SSS.bailouted;
            if print
                fprintf('fracWP0= %4.2f\n',fracWP0(k))
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
