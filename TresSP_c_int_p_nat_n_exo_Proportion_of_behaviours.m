%% MODELO 3 SP con int pos y neg QUE COMPITEN POR SL
close all 
clear all
clc

%% Function to sample parameters 
par=@(dummy)rand*(10^randi([-5, 2]));

%% Functions for the steady states
solBS =@(alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST)[                                                                                                                                                                                                                                                                                                                                                                                                                -(beta_BS - ST*alpha_BS)/alpha_BS;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0;
(alpha_ES*beta_PS*gamma_EB*mu_EB - alpha_PS*beta_ES*gamma_EB*mu_EB + alpha_BS*beta_ES*gamma_EP*mu_EP - alpha_ES*beta_BS*gamma_EP*mu_EP)/(alpha_PS*beta_ES*gamma_EB*mu_EB*teta_BP) - ((mu_EB + beta_BS*gamma_EB + ((alpha_ES*beta_BS^2*gamma_EB^2 - 2*alpha_ES*beta_BS*gamma_EB*mu_EB + 4*alpha_BS*beta_ES*gamma_EB*mu_EB + alpha_ES*mu_EB^2)/alpha_ES)^(1/2))*(alpha_ES*gamma_EB*mu_EB*mu_EP - alpha_ES*gamma_EP*mu_EB*mu_EP - alpha_ES*beta_BS*gamma_EB*gamma_EP*mu_EP + alpha_ES*beta_PS*gamma_EB*gamma_EP*mu_EB))/(2*alpha_PS*beta_ES*gamma_EB^2*mu_EB^2*teta_BP);
(alpha_ES*beta_PS*gamma_EB*mu_EB - alpha_PS*beta_ES*gamma_EB*mu_EB + alpha_BS*beta_ES*gamma_EP*mu_EP - alpha_ES*beta_BS*gamma_EP*mu_EP)/(alpha_PS*beta_ES*gamma_EB*mu_EB*teta_BP) - ((mu_EB + beta_BS*gamma_EB - ((alpha_ES*beta_BS^2*gamma_EB^2 - 2*alpha_ES*beta_BS*gamma_EB*mu_EB + 4*alpha_BS*beta_ES*gamma_EB*mu_EB + alpha_ES*mu_EB^2)/alpha_ES)^(1/2))*(alpha_ES*gamma_EB*mu_EB*mu_EP - alpha_ES*gamma_EP*mu_EB*mu_EP - alpha_ES*beta_BS*gamma_EB*gamma_EP*mu_EP + alpha_ES*beta_PS*gamma_EB*gamma_EP*mu_EB))/(2*alpha_PS*beta_ES*gamma_EB^2*mu_EB^2*teta_BP);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    (alpha_BS*beta_PS - alpha_PS*beta_BS)/(alpha_PS*beta_BS*teta_BP);
                                                                                                                                                                                                                                                                                                                           (ST*alpha_ES - beta_ES + (alpha_ES*(mu_EB + beta_BS*gamma_EB - ((alpha_ES*beta_BS^2*gamma_EB^2 - 2*alpha_ES*beta_BS*gamma_EB*mu_EB + 4*alpha_BS*beta_ES*gamma_EB*mu_EB + alpha_ES*mu_EB^2)/alpha_ES)^(1/2)))/(2*gamma_EB*mu_EB))/alpha_ES;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0;
                                                                                                                                                                                                                                                                                                                           (ST*alpha_ES - beta_ES + (alpha_ES*(mu_EB + beta_BS*gamma_EB + ((alpha_ES*beta_BS^2*gamma_EB^2 - 2*alpha_ES*beta_BS*gamma_EB*mu_EB + 4*alpha_BS*beta_ES*gamma_EB*mu_EB + alpha_ES*mu_EB^2)/alpha_ES)^(1/2)))/(2*gamma_EB*mu_EB))/alpha_ES;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0];
%%
solPS =@(alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST)[   0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  -(beta_PS - ST*alpha_PS)/alpha_PS;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  0;
(((mu_EB + beta_BS*gamma_EB + ((alpha_ES*beta_BS^2*gamma_EB^2 - 2*alpha_ES*beta_BS*gamma_EB*mu_EB + 4*alpha_BS*beta_ES*gamma_EB*mu_EB + alpha_ES*mu_EB^2)/alpha_ES)^(1/2))*(alpha_ES^2*gamma_EB*mu_EB*mu_EP - alpha_ES^2*gamma_EP*mu_EB*mu_EP - alpha_ES^2*beta_BS*gamma_EB*gamma_EP*mu_EP + alpha_ES^2*beta_PS*gamma_EB*gamma_EP*mu_EB + alpha_ES*alpha_PS*beta_ES*gamma_EB*mu_EB*teta_BP))/(2*gamma_EB*mu_EB) - alpha_ES^2*beta_PS*gamma_EB*mu_EB + alpha_ES^2*beta_BS*gamma_EP*mu_EP + alpha_ES*alpha_PS*beta_ES*gamma_EB*mu_EB - alpha_BS*alpha_ES*beta_ES*gamma_EP*mu_EP - alpha_PS*beta_ES^2*gamma_EB*mu_EB*teta_BP + ST*alpha_ES*alpha_PS*beta_ES*gamma_EB*mu_EB*teta_BP)/(alpha_ES*alpha_PS*beta_ES*gamma_EB*mu_EB*teta_BP);
(alpha_ES^2*beta_BS*gamma_EP*mu_EP - alpha_ES^2*beta_PS*gamma_EB*mu_EB + ((mu_EB + beta_BS*gamma_EB - ((alpha_ES*beta_BS^2*gamma_EB^2 - 2*alpha_ES*beta_BS*gamma_EB*mu_EB + 4*alpha_BS*beta_ES*gamma_EB*mu_EB + alpha_ES*mu_EB^2)/alpha_ES)^(1/2))*(alpha_ES^2*gamma_EB*mu_EB*mu_EP - alpha_ES^2*gamma_EP*mu_EB*mu_EP - alpha_ES^2*beta_BS*gamma_EB*gamma_EP*mu_EP + alpha_ES^2*beta_PS*gamma_EB*gamma_EP*mu_EB + alpha_ES*alpha_PS*beta_ES*gamma_EB*mu_EB*teta_BP))/(2*gamma_EB*mu_EB) + alpha_ES*alpha_PS*beta_ES*gamma_EB*mu_EB - alpha_BS*alpha_ES*beta_ES*gamma_EP*mu_EP - alpha_PS*beta_ES^2*gamma_EB*mu_EB*teta_BP + ST*alpha_ES*alpha_PS*beta_ES*gamma_EB*mu_EB*teta_BP)/(alpha_ES*alpha_PS*beta_ES*gamma_EB*mu_EB*teta_BP);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          -(alpha_BS^2*beta_PS - alpha_BS*alpha_PS*beta_BS + alpha_PS*beta_BS^2*teta_BP - ST*alpha_BS*alpha_PS*beta_BS*teta_BP)/(alpha_BS*alpha_PS*beta_BS*teta_BP);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          (ST*alpha_ES - beta_ES + (alpha_ES*(mu_EP + beta_PS*gamma_EP - ((alpha_ES*beta_PS^2*gamma_EP^2 - 2*alpha_ES*beta_PS*gamma_EP*mu_EP + 4*alpha_PS*beta_ES*gamma_EP*mu_EP + alpha_ES*mu_EP^2)/alpha_ES)^(1/2)))/(2*gamma_EP*mu_EP))/alpha_ES;
 0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          (ST*alpha_ES - beta_ES + (alpha_ES*(mu_EP + beta_PS*gamma_EP + ((alpha_ES*beta_PS^2*gamma_EP^2 - 2*alpha_ES*beta_PS*gamma_EP*mu_EP + 4*alpha_PS*beta_ES*gamma_EP*mu_EP + alpha_ES*mu_EP^2)/alpha_ES)^(1/2)))/(2*gamma_EP*mu_EP))/alpha_ES];
%%
solES =@(alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST)[  0;
                                                                                                                                                           -(beta_ES - ST*alpha_ES)/alpha_ES;
                                                                                                                                                                                           0;
                                                                                                                                                                                           0;
-(mu_EB + beta_BS*gamma_EB + ((alpha_ES*beta_BS^2*gamma_EB^2 - 2*alpha_ES*beta_BS*gamma_EB*mu_EB + 4*alpha_BS*beta_ES*gamma_EB*mu_EB + alpha_ES*mu_EB^2)/alpha_ES)^(1/2))/(2*gamma_EB*mu_EB);
-(mu_EB + beta_BS*gamma_EB - ((alpha_ES*beta_BS^2*gamma_EB^2 - 2*alpha_ES*beta_BS*gamma_EB*mu_EB + 4*alpha_BS*beta_ES*gamma_EB*mu_EB + alpha_ES*mu_EB^2)/alpha_ES)^(1/2))/(2*gamma_EB*mu_EB);
                                                                                                                                                                                           0;
-(mu_EB + beta_BS*gamma_EB - ((alpha_ES*beta_BS^2*gamma_EB^2 - 2*alpha_ES*beta_BS*gamma_EB*mu_EB + 4*alpha_BS*beta_ES*gamma_EB*mu_EB + alpha_ES*mu_EB^2)/alpha_ES)^(1/2))/(2*gamma_EB*mu_EB);
-(mu_EP + beta_PS*gamma_EP - ((alpha_ES*beta_PS^2*gamma_EP^2 - 2*alpha_ES*beta_PS*gamma_EP*mu_EP + 4*alpha_PS*beta_ES*gamma_EP*mu_EP + alpha_ES*mu_EP^2)/alpha_ES)^(1/2))/(2*gamma_EP*mu_EP);
-(mu_EB + beta_BS*gamma_EB + ((alpha_ES*beta_BS^2*gamma_EB^2 - 2*alpha_ES*beta_BS*gamma_EB*mu_EB + 4*alpha_BS*beta_ES*gamma_EB*mu_EB + alpha_ES*mu_EB^2)/alpha_ES)^(1/2))/(2*gamma_EB*mu_EB);
-(mu_EP + beta_PS*gamma_EP + ((alpha_ES*beta_PS^2*gamma_EP^2 - 2*alpha_ES*beta_PS*gamma_EP*mu_EP + 4*alpha_PS*beta_ES*gamma_EP*mu_EP + alpha_ES*mu_EP^2)/alpha_ES)^(1/2))/(2*gamma_EP*mu_EP)];
 

%% Functions for the Jacobian

Jacobian =@(BS, PS, ES, alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST)[[ - beta_BS - ES*mu_EB - (alpha_BS*(BS + ES + PS - ST))/(ES*gamma_EB + 1) - (BS*alpha_BS)/(ES*gamma_EB + 1),                                                                                                            -(BS*alpha_BS)/(ES*gamma_EB + 1),                                   (BS*alpha_BS*gamma_EB*(BS + ES + PS - ST))/(ES*gamma_EB + 1)^2 - (BS*alpha_BS)/(ES*gamma_EB + 1) - BS*mu_EB]
[ - (PS*alpha_PS*(BS*teta_BP + 1))/(ES*gamma_EP + 1) - (PS*alpha_PS*teta_BP*(BS + ES + PS - ST))/(ES*gamma_EP + 1), - beta_PS - ES*mu_EP - (PS*alpha_PS*(BS*teta_BP + 1))/(ES*gamma_EP + 1) - (alpha_PS*(BS*teta_BP + 1)*(BS + ES + PS - ST))/(ES*gamma_EP + 1), (PS*alpha_PS*gamma_EP*(BS*teta_BP + 1)*(BS + ES + PS - ST))/(ES*gamma_EP + 1)^2 - (PS*alpha_PS*(BS*teta_BP + 1))/(ES*gamma_EP + 1) - PS*mu_EP]
[                                                                                                     -ES*alpha_ES,                                                                                                                                -ES*alpha_ES,                                                                                        - beta_ES - ES*alpha_ES - alpha_ES*(BS + ES + PS - ST)]];



for jj=1:100000 % iterations 
    
    

alpha_BS=par(1);
gamma_EB=par(1);
beta_BS=par(1);
mu_EB=par(1);
alpha_PS=par(1);
teta_BP=par(1);
gamma_EP=par(1);
beta_PS=par(1);
mu_EP=par(1);
alpha_ES=par(1);
beta_ES=par(1);


ST=146.9;


%% % evaluate the steady state solutions
BSss=solBS(alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST);
PSss=solPS(alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST);
ESss=solES(alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST);

%% filter out those that are positive , real and stable
counter_stable=0;
counter_positive_real=1;

for ii=1:1:11
    
   Jaci= Jacobian(BSss(ii), PSss(ii), ESss(ii), alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST);
   
   
   
   if sum([BSss(ii), PSss(ii), ESss(ii)]>=0)==3 && isreal([BSss(ii), PSss(ii), ESss(ii)]) && sum(eig(Jaci)<=0)==3
       
             counter_stable=counter_stable+1;
  
       behaviour_matrix(jj,2+3*(counter_stable-1))=BSss(ii);
       behaviour_matrix(jj,3+3*(counter_stable-1))=PSss(ii);
       behaviour_matrix(jj,4+3*(counter_stable-1))=ESss(ii);
       
            
   end
        
    
end


behaviour_matrix(jj,1)=counter_stable; % number of steady states
end

%%
save('Behaviour_matrix_ModelIIB.mat', 'behaviour_matrix')


