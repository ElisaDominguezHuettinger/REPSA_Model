% %% Modelo_2_especies que cooperan
close all 
clear all
clc

%% En este codigo vamos a jugar con parametros, asi que lo primero es poner los valores
alpha_PS=.0011;
beta_PS=.29; 
teta_BP=.0995; 
alpha_BS=.00137; 
beta_BS=.07; 
ST=146.9; % hectareas de la repsa

%conjunto de parametros con buen ajuste
%{
alpha_PS=.0011;
beta_PS=.29; 
teta_BP=.0995; 
alpha_BS=.00137; 
beta_BS=.07; 
ST=146.9;

%}
%% cargar los datos
load('Indiv_tep_Palo_Euc.mat')
load('Cob_tep_Palo_Euc.mat')

%%
%figure;
subplot(1,2,1)
scatter(Cob_tep_Palo_Euc(:,1), Cob_tep_Palo_Euc(:,2)+Cob_tep_Palo_Euc(:,3), 'rs', 'filled');
hold on
ylabel('Tepozán (arboles + arbustos)')
xlabel('tiempo')

subplot(1,2,2)
scatter(Cob_tep_Palo_Euc(:,1), Cob_tep_Palo_Euc(:,4), 'db', 'filled');
hold on
ylabel('Palo loco')
xlabel('tiempo')



%% declarar funciones que usaremos mas adelante

% declara la ODE;
Modelo_2_especies=@(t,y,alpha_BS, ST, beta_BS, alpha_PS, beta_PS, teta_BP)[(alpha_BS*y(1)*(ST-y(1)-y(2)))-(beta_BS*y(1));(alpha_PS*y(2)*(ST-y(1)-y(2))*(1+teta_BP*y(1)))-(beta_PS*y(2))];



%% Integrar numéricamente
vec_tiempo=Cob_tep_Palo_Euc(:,1) 
vec_Var_estado_1=Cob_tep_Palo_Euc(:,2)+Cob_tep_Palo_Euc(:,3)  
vec_Var_estado_2=Cob_tep_Palo_Euc(:,4)

tspan=[vec_tiempo(2,1) vec_tiempo(4,1)]-vec_tiempo(2,1); % intervalo de integración
y0 = [vec_Var_estado_1(2) vec_Var_estado_2(2)]; % intervalo de integracion

%seleccionamos nuestro integrador favorito * (euler) (ode15, ode23,..)
[t,y] = ode45(@(t,y)Modelo_2_especies(t,y,alpha_BS, ST, beta_BS, alpha_PS, beta_PS, teta_BP),tspan,y0);
%%
t=t+vec_tiempo(2,1)
%% graficamos
subplot(1,2,1)
plot(t, y(:,1), 'r')
hold on

subplot(1,2,2)
plot(t, y(:,2), 'b')




%% Soluciones de equilibrio
solBS =@(alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)[ -(beta_BS - ST.*alpha_BS)./alpha_BS;
                                                               0;
                                                               0;
(alpha_BS.*beta_PS - alpha_PS.*beta_BS)./(alpha_PS.*beta_BS.*teta_BP)];

solPS =@(alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)[                                                                                            0;
                                                                                                                        -(beta_PS - ST.*alpha_PS)./alpha_PS;
                                                                                                                                                        0;
-(alpha_BS.^2.*beta_PS - alpha_BS.*alpha_PS.*beta_BS + alpha_PS.*beta_BS.^2.*teta_BP - ST.*alpha_BS.*alpha_PS.*beta_BS.*teta_BP)./(alpha_BS.*alpha_PS.*beta_BS.*teta_BP)];

%% Ceroclinas
% declarar la ceroclina como una funcion de PS
BSss2_fun=@(PS,alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP,ST)-(beta_BS + PS.*alpha_BS - ST.*alpha_BS)./alpha_BS;
% declarar la ceroclina como una funcion de BS
PSss2_fun=@(BS,alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)-(beta_PS + BS.*alpha_PS - ST.*alpha_PS + BS.^2.*alpha_PS*teta_BP - BS.*ST.*alpha_PS*teta_BP)./(alpha_PS + BS.*alpha_PS*teta_BP);

%% Jacobiano 
Jac=@(BS,PS,alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)[[                  - beta_BS - BS*alpha_BS - alpha_BS*(BS + PS - ST),                                                                        -BS*alpha_BS]
[- PS*alpha_PS*(BS*teta_BP + 1) - PS*alpha_PS*teta_BP*(BS + PS - ST), - beta_PS - PS*alpha_PS*(BS*teta_BP + 1) - alpha_PS*(BS*teta_BP + 1)*(BS + PS - ST)]];

 
%%
teta_BP_crit=-(beta_PS*alpha_BS^2 - alpha_PS*beta_BS*alpha_BS)/(alpha_PS*beta_BS^2 - ST*alpha_BS*alpha_PS*beta_BS)


%% Grafica el equilibrio que es estable

BSss=solBS(alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)
PSss=solPS(alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)

J1=Jac(BSss(1), PSss(1),alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST);
J2=Jac(BSss(2), PSss(2),alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST);
J3=Jac(BSss(3), PSss(3),alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST);
J4=Jac(BSss(4), PSss(4),alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST);



% valuamos estabilidad
if (sum(eig(J1)>0) >0)
   S(1)=0;
else
    S(1)=1;
end

if (sum(eig(J2)>0) >0)
   S(2)=0;
else
    S(2)=1;
end

if (sum(eig(J3)>0) >0)
   S(3)=0;
else
    S(3)=1;
end

if (sum(eig(J4)>0) >0)
   S(4)=0;
else
    S(4)=1;
end

S

subplot(1,2,1)
% grafica solo los equilbrios estables, es decir, donde la matriz S=1
line([t(1) t(end)], [BSss(S==1) BSss(S==1)])

subplot(1,2,2)
% grafica solo los equilbrios estables, es decir, donde la matriz S=1
%line([t(1) t(end)], [PSss(S==1) PSss(S==1)])


title(['Cobertura en 146.9 ha de REPSA, \alpha_{PS}=' num2str(alpha_PS)  '\beta_{PS}=' num2str(beta_PS)  ' \alpha_{BS}=' num2str(alpha_BS) '\beta_{BS}=' num2str(beta_BS) 'S_T=' num2str(ST)])

