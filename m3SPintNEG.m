close all 
clear all
clc


%%
syms alpha_BS gamma_EB beta_BS mu_EB  alpha_PS gamma_BP gamma_EP beta_PS mu_EP  alpha_ES beta_ES
syms BS PS ES SL 

%% ODES
dBS=((alpha_BS*BS*SL)/(1+(ES*gamma_EB)))-(BS*(beta_BS+(mu_EB*ES))) ; %Tepozanes
dPS=((alpha_PS*PS*SL)/((1+(BS*gamma_BP))+(ES*gamma_EP)))   -(PS*(beta_PS+(mu_EP*ES))) ; %Palo locos
dES=(alpha_ES*ES*SL)-(beta_ES*ES); %Eucaliptos
dSL= -((alpha_BS*BS*SL)/(1+(ES*gamma_EB)))+(BS*(beta_BS+(mu_EB*ES))) -((alpha_PS*PS*SL)/((1+(BS*gamma_BP))+(ES*gamma_EP))) +(PS*(beta_PS+(mu_EP*ES))) -(alpha_ES*ES*SL)+(beta_ES*ES); %espacio disponible
%% checar que se cumplan ecuaciones de conservación
dBS+dPS+dES+dSL