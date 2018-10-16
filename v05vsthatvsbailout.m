%% This code plots V(0.5) over that and bailout

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
    that         =[0.01:0.001:0.035];K=length(that);
    bail         =[0:0.1:0.5];L=length(bail);
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
%% Plot contour
% Select cases to compare

V05real=real(V05);V05imag=imag(V05);V05imag_log=(V05imag~=0);%contour(that,g2,V05imag_log,80)
%mesh(that,g2,V05)
figure;set(gcf,'units','normalized','position',[0.01,0.25,0.65,0.65])
contour3(that,bail,V05real,80);%contour(that,g2,V05real,80)
xlabel('that'), ylabel('bail'), title('V05')%xlabel('that'), ylabel('g2'), title('V05')
saveas(gcf,[pwd '/figures/mesh/bailout7'])
saveas(gcf,[pwd '/figures/mesh/bailout7.png'])

%% Plot V(0.5) and EV for different bailouts
% Choose which bailouts
pos1=3;pos2=4;
% Mean value
EV1=mean(VAL1,3);

figure
hold on
plot(that(:,1),V05(:,1),'-','LineWidth',2)
plot(that(:,pos1),V05(:,pos1),'--','LineWidth',2)
plot(that(:,pos2),V05(:,pos2),'-.','LineWidth',2)
xlabel('$\hat{t}$','Interpreter','Latex','FontSize',16), ylabel('V(0.5)','FontSize',16), title('Bailouts under commitment','FontSize',16)%xlabel('that'), ylabel('g2'), title('V05')
legend({'bailout=0','bailout=0.3$*C_{WP}$','bailout=0.4$*C_{WP}$'},'Interpreter','Latex','FontSize',14,'Location','southeast')
   saveas(gcf,[pwd '/figures/unidim/V05_vs_bail2'])
   saveas(gcf,[pwd '/figures/unidim/V05_vs_bail2.png'])

   figure
hold on
plot(that(:,1),V075(:,1),'-','LineWidth',2)
plot(that(:,pos1),V075(:,pos1),'--','LineWidth',2)
plot(that(:,pos2),V075(:,pos2),'-.','LineWidth',2)
xlabel('$\hat{t}$','Interpreter','Latex','FontSize',16), ylabel('V(0.75)','FontSize',16), title('Bailouts under commitment','FontSize',16)%xlabel('that'), ylabel('g2'), title('V05')
legend({'bailout=0','bailout=0.3$*C_{WP}$','bailout=0.4$*C_{WP}$'},'Interpreter','Latex','FontSize',14,'Location','southeast')
   saveas(gcf,[pwd '/figures/unidim/V075_vs_bail2'])
   saveas(gcf,[pwd '/figures/unidim/V075_vs_bail2.png'])
   
   
figure
hold on
plot(that(:,1),EV1(:,1),'-','LineWidth',2)
plot(that(:,pos1),EV1(:,pos1),'--','LineWidth',2)
plot(that(:,pos2),EV1(:,pos2),'-.','LineWidth',2)
xlabel('$\hat{t}$','Interpreter','Latex','FontSize',16), ylabel('EV','FontSize',16), title('Bailouts under commitment','FontSize',16)%xlabel('that'), ylabel('g2'), title('V05')
legend({'bailout=0','bailout=0.3$*C_{WP}$','bailout=0.4$*C_{WP}$'},'Interpreter','Latex','FontSize',14,'Location','southeast')
   saveas(gcf,[pwd '/figures/unidim/EV_vs_bail2'])
   saveas(gcf,[pwd '/figures/unidim/EV_vs_bail2.png'])
   
   
   
%    % clear all 
% % clc
% % load('data/mesh/plotmesh_nobailout6') 
% [fval_t,tmax_temp]=max(V05);
% [fval_g,gmax_temp]=max(V05');
% [fval_glo,gmax_ind]=max(fval_t');
% [fval_glo2,tmax_ind]=max(fval_g);
% that_max=that(tmax_ind,gmax_ind)
% %g2_max  =g2(tmax_ind,gmax_ind)
% save('data/mesh/plotmesh_bailout7.mat')
%    