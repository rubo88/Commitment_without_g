%% Simple code to call commitment2.m and experiment

clear all
clc
close all
diary on
addpath('..','funciones');  

%%
    S.options.ploteach=1;
    S.options.print=1;
%% Parameters
    S.N=101;N=S.N;
    S.constVol= false; % If true: s(P)=sigma is constant.% If false: s(P)= 4*P*(1-P)*sigma
    %S.policyrules=@linearrules;
    S.policyrules=@linearrules_bailout;
    %S.X=@cara;
% Preference parameters:
    S.rho       = 0.04;         % Time discount rate.
    S.sigma     = 0.02;    
    S.sigma2    = S.sigma^2;
    S.Cwp       = S.rho;
    S.kappa     = 0.7;
    S.alpha     = 1;
    S.that      =0.025;
    S.bail      =0;
    
    S=commitment2(S);