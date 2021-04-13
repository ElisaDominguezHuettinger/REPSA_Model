%% MODELO 3 SP con int pos y neg QUE COMPITEN POR SL
%Modelo mas realista
close all 
clear all
clc

%%
%% cargar los datos
load('Indiv_tep_Palo_Euc.mat')
load('Cob_tep_Palo_Euc.mat')
%%

%Dos puntos de equilibrio... creo que mu_EB es uno de los determinantes....
%de ser así tendríamos que analizar más a detalle la interacción entre el
%tepozán y el eucalipto... lo cual es lógico ecológicamente... Si la
%Budleia es el dominante y es sombrilla del palo loco, entonces lo mnás
%relevante para conocer el futuro del sistema es cnocer su relación... Si
%varío el valor de muEB de .01 a .0001 pasa de biestabilidad a un estado
%dominado por los eucaliptos.

alpha_BS=.00137;
gamma_EB=.01;
beta_BS=.07;
mu_EB=.01;
alpha_PS=.0011;
teta_BP=.0995;
gamma_EP=.01;
beta_PS=.29;
mu_EP=.01;
alpha_ES=.001; 
beta_ES=.1;
ST=146.9;
%}
%% LO SIMBOLICO SE HACE UNA VEZ NOMAS POR ESO ESTA AHORA COMENTADO.
% YA LO GUARDAMOS EN FUNCIONES

%% declaramos las variables simbolicas:
%{
syms alpha_BS gamma_EB beta_BS mu_EB  alpha_PS teta_BP gamma_EP beta_PS mu_EP  alpha_ES beta_ES
syms BS PS ES SL  
%% ODES
dBS=((alpha_BS*BS*SL)/(1+(ES*gamma_EB)))-(BS*(beta_BS+(mu_EB*ES))) ; %Tepozanes
dPS=(((alpha_PS*PS*SL)*(1+teta_BP*BS))/(1+(ES*gamma_EP)))-(PS*(beta_PS+(mu_EP*ES))) ; %Palo locos
dES=(alpha_ES*ES*SL)-(beta_ES*ES); %Eucaliptos
dSL=-((alpha_BS*BS*SL)/(1+(ES*gamma_EB)))+(BS*(beta_BS+(mu_EB*ES)))-(((alpha_PS*PS*SL)*(1+teta_BP*BS))/(1+(ES*gamma_EP)))+(PS*(beta_PS+(mu_EP*ES))) -(alpha_ES*ES*SL)+(beta_ES*ES); %espacio disponible
%% checar que se cumplan ecuaciones de conservaciÃ³n
dBS+dPS+dES+dSL
%% reducir el sistema a dos dimensines
syms ST
SL=(ST-BS-PS-ES)
dBS=((alpha_BS*BS*SL)/(1+(ES*gamma_EB)))-(BS*(beta_BS+(mu_EB*ES))) ;
dPS=((alpha_PS*PS*SL)*(1+teta_BP*BS))/(1+(ES*gamma_EP))-(PS*(beta_PS+(mu_EP*ES))) ;
dES=(alpha_ES*ES*SL)-(beta_ES*ES) ;
%% sacamos los puntos de equilibrio

[solBS,solPS,solES] = solve(dBS== 0, dPS== 0, dES==0, [BS PS ES])

%}
%% y los acomodamos dentro de funciones
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
 

%% sacamos el jacobiano

%jacobian([dBS, dPS, dES], [BS, PS, ES])

%%

Jacobian =@(BS, PS, ES, alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST)[[ - beta_BS - ES*mu_EB - (alpha_BS*(BS + ES + PS - ST))/(ES*gamma_EB + 1) - (BS*alpha_BS)/(ES*gamma_EB + 1),                                                                                                            -(BS*alpha_BS)/(ES*gamma_EB + 1),                                   (BS*alpha_BS*gamma_EB*(BS + ES + PS - ST))/(ES*gamma_EB + 1)^2 - (BS*alpha_BS)/(ES*gamma_EB + 1) - BS*mu_EB]
[ - (PS*alpha_PS*(BS*teta_BP + 1))/(ES*gamma_EP + 1) - (PS*alpha_PS*teta_BP*(BS + ES + PS - ST))/(ES*gamma_EP + 1), - beta_PS - ES*mu_EP - (PS*alpha_PS*(BS*teta_BP + 1))/(ES*gamma_EP + 1) - (alpha_PS*(BS*teta_BP + 1)*(BS + ES + PS - ST))/(ES*gamma_EP + 1), (PS*alpha_PS*gamma_EP*(BS*teta_BP + 1)*(BS + ES + PS - ST))/(ES*gamma_EP + 1)^2 - (PS*alpha_PS*(BS*teta_BP + 1))/(ES*gamma_EP + 1) - PS*mu_EP]
[                                                                                                     -ES*alpha_ES,                                                                                                                                -ES*alpha_ES,                                                                                        - beta_ES - ES*alpha_ES - alpha_ES*(BS + ES + PS - ST)]];




%%
BSss=solBS(alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST);
PSss=solPS(alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST);
ESss=solES(alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST);

%%
counter_stable=0;
counter_positive_real=1;
%positive_real_solutions;
for ii=1:1:11
    
   Jaci= Jacobian(BSss(ii), PSss(ii), ESss(ii), alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST);
   
   
   
   if sum([BSss(ii), PSss(ii), ESss(ii)]>=0)==3 && isreal([BSss(ii), PSss(ii), ESss(ii)]) % positive 
   
       matrix_pos_real(counter_positive_real,1)=BSss(ii);
       matrix_pos_real(counter_positive_real,2)=PSss(ii);
       matrix_pos_real(counter_positive_real,3)=ESss(ii);
       
       
          
   if sum(eig(Jaci)<=0)==3 %
        stable_sol=ii;
        %disp(['stable' ii])
        counter_stable=counter_stable+1;
        
        matrix_pos_real(counter_positive_real,4)=1;
        
   else
       matrix_pos_real(counter_positive_real,4)=0;
    end;
        
    counter_positive_real=counter_positive_real+1;
   end;
end;
%%
[BSss(stable_sol), PSss(stable_sol), ESss(stable_sol)]


% bajo estos valores de parámetros, la solucion que domina es.

%% vamos a ver 

%[BSss, PSss, ESss]

 matrix_pos_real
 
 
 %% empezamos con la integración numérica
 
Modelo_3_especies=@(t,y,alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST)[((alpha_BS*y(1)*(ST-y(1)-y(2)-y(3)))/(1+(y(3)*gamma_EB)))-(y(1)*(beta_BS+(mu_EB*y(3)))) ;((alpha_PS*y(2)*(ST-y(1)-y(2)-y(3)))*(1+teta_BP*y(1)))/(1+(y(3)*gamma_EP))-(y(2)*(beta_PS+(mu_EP*y(3)))) ; (alpha_ES*y(3)*(ST-y(1)-y(2)-y(3)))-(beta_ES*y(3))] ;


%% Integrar numéricamente
vec_tiempo=Cob_tep_Palo_Euc(:,1) ;
vec_Var_estado_1=Cob_tep_Palo_Euc(:,2)+Cob_tep_Palo_Euc(:,3)  ;
vec_Var_estado_2=Cob_tep_Palo_Euc(:,4);
vec_Var_estado_3=Cob_tep_Palo_Euc(:,5);


%% dado que el primer (y unico) tiempo para el que tenemos datos para las tres especies es el segundo, partimos de ahi
% tanto para definir el intervalo de integración como las condiciones
% iniciales
tspan=[vec_tiempo(2,1) 20+vec_tiempo(end,1)]-vec_tiempo(2,1); % intervalo de integración
y0 = [vec_Var_estado_1(2) vec_Var_estado_2(2) .5*vec_Var_estado_3(2)]; % intervalo de integracion

%seleccionamos nuestro integrador favorito * (euler) (ode15, ode23,..)
[t,y] = ode45(@(t,y)Modelo_3_especies(t,y,alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST),tspan,y0);
%%
t=t+vec_tiempo(2,1)
%% graficamos
subplot(1,3,1)
plot(t, y(:,1), 'r')
hold on
scatter(Cob_tep_Palo_Euc(:,1), Cob_tep_Palo_Euc(:,2)+Cob_tep_Palo_Euc(:,3), 'rs', 'filled');
ylabel('Tepozán (cobertura, arbustos + arboles')

subplot(1,3,2)
plot(t, y(:,2), 'b')
hold on 
scatter(Cob_tep_Palo_Euc(:,1), Cob_tep_Palo_Euc(:,4), 'db', 'filled');
ylabel('Palo loco (cobertura)')

subplot(1,3,3)
plot(t, y(:,3), 'g')
hold on
ylabel('Eucalipto (cobertura')
scatter(Cob_tep_Palo_Euc(:,1), Cob_tep_Palo_Euc(:,5), 'pg', 'filled');
%%
index_stable_steady_state_sol=find(matrix_pos_real(:,4)==1); %encuentra los renglones de la matriz que tienen las soluciones de equilibrio estables

%Col = hsv(11);
for ii=1:size(index_stable_steady_state_sol) % por cada una de las sol estables, ...
    
  
if ii==1
    Col='m'
    Style='--'
elseif ii==2
    Col='k'
    Style=':'
elseif ii==3
    Col='c'
    Style='-'
    
end;
        
    subplot(1,3,1)
    hold on
    line([t(1) t(end)], [matrix_pos_real(index_stable_steady_state_sol(ii),1) matrix_pos_real(index_stable_steady_state_sol(ii),1) ], 'Color', Col);%s, 'LineStyle', 'Style')
    
    subplot(1,3,2)
    
    hold on
    line([t(1) t(end)], [matrix_pos_real(index_stable_steady_state_sol(ii),2) matrix_pos_real(index_stable_steady_state_sol(ii),2) ], 'Color', Col);%,'LineStyle', 'Style')
    
    subplot(1,3,3)
    line([t(1) t(end)], [matrix_pos_real(index_stable_steady_state_sol(ii),3) matrix_pos_real(index_stable_steady_state_sol(ii),3) ] ,'Color',  Col);%,'LineStyle', 'Style')
    
    
end;

%% Jugar con condiciones iniciales
y0 = [vec_Var_estado_1(2) vec_Var_estado_2(2)  vec_Var_estado_3(2)*.1] % intervalo de integracion
%

tspan2=[vec_tiempo(2,1) 1000+vec_tiempo(end,1)]-vec_tiempo(2,1)
[t2,y2] = ode45(@(t,y)Modelo_3_especies(t,y,alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST),tspan2,y0);
%%
t2=t2+vec_tiempo(2,1);
if ii>1 % si hay más de una solucion de equilibrio, ver a cual converge
 for jj=1:size(index_stable_steady_state_sol)
  
  Error=sum(abs( y2(end,:)-matrix_pos_real(index_stable_steady_state_sol(jj),1:3)))
  
  if Error <0.5
      sol_converge_a=jj
      
      
  end;
    
 end;   
end;

if sol_converge_a==1
    Col='m'
    Style='--'
elseif sol_converge_a==2
    Col='k'
    Style=':'
elseif sol_converge_a==3
    Col='c'
    Style='-'
    
end;

%% graficamos
subplot(1,3,1)
plot(t2, y2(:,1),'Color', Col, 'LineStyle', Style)
xlim([vec_tiempo(1)-10, t(end)])

subplot(1,3,2)
plot(t2, y2(:,2),'Color', Col, 'LineStyle', Style)
xlim([vec_tiempo(1)-10, t(end)])

subplot(1,3,3)
plot(t2, y2(:,3),'Color', Col, 'LineStyle', Style)
xlim([vec_tiempo(1)-10, t(end)])




%%












title(['\alpha_B_S, \gamma_E_B, \beta_B_S, \mu_E_B,  \alpha_P_S, \theta_B_P, \gamma_E_P, \beta_P_S, \mu_E_P,  \alpha_E_S, \beta_E_S, ST' num2str([alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, teta_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST])]);

