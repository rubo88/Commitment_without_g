%% This code plots V(0.5) over that and g2
%% Experiments
% 1 and 2 bail=0.2
% 3       bail=0.4
% 4       bail=0.6
% 5       fix g2 and move bail
% 6       fix g2 and move bail (finer)  
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
    S.sigma     = 0.4;    
    S.sigma2    = S.sigma^2;
    S.Cwp       = S.rho;
    S.kappa     = 0.7;
    S.alpha     = 1;
    
%% Cases
    that         =[0.01:0.001:0.05];K=length(that);
    bail         =[0:0.005:0.5];L=length(bail);
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
% Select cases to compare

V05real=real(V05);V05imag=imag(V05);V05imag_log=(V05imag~=0);%contour(that,g2,V05imag_log,80)
%mesh(that,g2,V05)
figure;set(gcf,'units','normalized','position',[0.01,0.25,0.65,0.65])
contour3(that,bail,V05real,80);%contour(that,g2,V05real,80)
xlabel('that'), ylabel('bail'), title('V05')%xlabel('that'), ylabel('g2'), title('V05')
saveas(gcf,[pwd '/figures/mesh/bailout7'])
saveas(gcf,[pwd '/figures/mesh/bailout7.png'])

% clear all 
% clc
% load('data/mesh/plotmesh_nobailout6') 
[fval_t,tmax_temp]=max(V05);
[fval_g,gmax_temp]=max(V05');
[fval_glo,gmax_ind]=max(fval_t');
[fval_glo2,tmax_ind]=max(fval_g);
that_max=that(tmax_ind,gmax_ind)
%g2_max  =g2(tmax_ind,gmax_ind)
save('data/mesh/plotmesh_bailout7.mat')


plot(bail(16,:),V05(16,:));hold all;plot(bail(46,:),V05(46,:));hold all;plot(bail(87,:),V05(87,:));
xlabel('bail'), ylabel('V(0.5)'), title('Bailouts')%xlabel('that'), ylabel('g2'), title('V05')
legend('that=2.3','that=2.7','that=3.1')