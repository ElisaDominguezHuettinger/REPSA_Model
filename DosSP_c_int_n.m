%% modelo 2 especies que compiten
close all 
clear all
clc

%%% TODA LA PARTE SIMBOLICA SOLO SE TIENE QUE HACER UNA VEZ - AHORA YA
%%% COMENTADO; LO MARCO COMO SYM
%% declaramos las variables simbolicas:
%%%%%Acabando de sacar lo que necito lo pongo como %% SYM SECCION
%%%% SACAR EXPRESIONES SIMBÓLICAS PARA LOS PUNTOS DE EQUIIBRIO
%
syms alpha_BS beta_BS alpha_PS beta_PS gamma_BP
syms BS PS SL  

% sistma de odes:
dBS=(alpha_BS*BS*SL)-(beta_BS*BS)  %tepozanes

dPS=(alpha_PS*PS*SL/(1+gamma_BP*BS))-(beta_PS*PS) %palolocos

dSL=-(alpha_BS*BS*SL+alpha_PS*PS*SL/(1+gamma_BP*BS))+(beta_BS*BS+beta_PS*PS) %espacio


%% checar que se cumplan ecuaciones de conservacion
dBS+dPS+dSL
%% reducir el sistema a dos dimensines
syms ST
SL=(ST-BS-PS)
dBS=(alpha_BS*BS*SL)-(beta_BS*BS)  ;
dPS=(alpha_PS*PS*SL/(1+gamma_BP*BS))-(beta_PS*PS);
%% sacamos los puntos de equilibrio

[solBSsym,solPSsym] = solve(dBS== 0, dPS== 0,[BS PS])
%}
%%
solBS =@(alpha_BS, beta_BS, alpha_PS, beta_PS, gamma_BP, ST)[  -(beta_BS - ST*alpha_BS)/alpha_BS;
                         0;
                                                                0;
    -(alpha_BS*beta_PS - alpha_PS*beta_BS)/(alpha_BS*beta_PS*gamma_BP)];

%%
solPS =@(alpha_BS, beta_BS, alpha_PS, beta_PS, gamma_BP, ST)[0;
    -(beta_PS - ST*alpha_PS)/alpha_PS;
    0;
(alpha_BS*beta_PS - alpha_PS*beta_BS - beta_BS*beta_PS*gamma_BP + ST*alpha_BS*beta_PS*gamma_BP)/(alpha_BS*beta_PS*gamma_BP)];


%% ya desde aca vemos que hay 4 posibes soluciones: una gana y mata a la otra, vice versa, ambas mueren o ¿?
 %%% aCA VA OTRA PARTE SIMBOLICA PARA SACAR LAS CEROCLINAS; ESTO ES
 %%% IMPORTNTE PARA TENER UNA NOCION GEOMETRICA DE LA ESTABILIDAD Y DEL
 %%% SISTEMA EN GRAL
 %
%% sacamos las ceroclinas (con la sig linea encontramos las soluciones de la eq.

% recordar que las ceroclinas del sistema son las curvas donde esa variable
% es==0; en la sig linea obtenemos la ceroclina de BS. Eso significa que en
% toda esa curva BS esta en equilibrio; pero nota que la curva  es una
% funciòn de PS!
sol1=solve(dBS==0,BS)
%% soluciones obtenidad de la línea aterior (copiar y pegar de linea de comando)
BSss1=0 % =sol1(1)
BSss2=-(beta_BS + PS*alpha_BS - ST*alpha_BS)/alpha_BS % = sol1(2)

%% lo mismo para la otra variable
sol2=solve(dPS==0,PS)
%% Soluciones obtenidas de la línea anterior
PSss1=0 % =sol2(1)
PSss2=-(beta_PS + BS*alpha_PS - ST*alpha_PS + BS*beta_PS*gamma_BP)/alpha_PS %sol2(2)
%}
 
 %% vamos a escribir la ceroclinas de PS como una funcion de BS y viceversa: (el  arroba se usa para declarar funciones "anonimas", es decir, en el mismo archivo en el que estamos trabajando
%% notar tambien que como queremos poder meter vectores a esta funcion (vectores de BS por ejemplo), necesitamos mdificar a funcion de la ceroclina ta que acepte
% vectores; en particular, sustitur * por .* y / por ./ y ^ por .^ (el
% punto indica multiplicacion/divicion / a la/  entrada por entrada de vector
%% declarar la ceroclina como una función de PS
BSss2_fun=@(PS,alpha_BS, beta_BS, alpha_PS, beta_PS,gamma_BP,ST)-(beta_BS + PS.*alpha_BS - ST.*alpha_BS)./alpha_BS
%% declarar la ceroclina como una función de BS
PSss2_fun=@(BS,alpha_BS, beta_BS, alpha_PS, beta_PS, gamma_BP, ST)-(beta_PS + BS.*alpha_PS - ST.*alpha_PS + BS.*beta_PS*gamma_BP)./alpha_PS





 
 
%%  sacamos el Jacobiano, lo vamos a necesitar más adelamnte para evaluar
% equilibrio
%%% SYM

jacobian([dBS, dPS], [BS, PS])
%}
%% recordar que la matriz jacobiana es una linearizacion del sistema; es decir, es el equivalente lineal del sistema no lineal .. esto es valido
% solo cerca de un punto ; entionces esta expresion se deberà evaluar en el
% punto de equilibrio que obtengamos mpas adelante. 
% una vez que evaluamos el jacobiano en el punto de equilbrio, hacemos un
% analisis deestabilidad local, es decir, sacamos los eigenvalores y
% eigenvectrs y pedimos que sean todos nevativos (convergencia asintòtica)
 
Jac=@(BS,PS,alpha_BS, beta_BS, alpha_PS, beta_PS, gamma_BP, ST)[[- beta_BS - BS*alpha_BS - alpha_BS*(BS + PS - ST),                                                                              -BS*alpha_BS]
[(PS*alpha_PS*gamma_BP*(BS + PS - ST))/(BS*gamma_BP + 1)^2 - (PS*alpha_PS)/(BS*gamma_BP + 1), - beta_PS - (alpha_PS*(BS + PS - ST))/(BS*gamma_BP + 1) - (PS*alpha_PS)/(BS*gamma_BP + 1)]]
 
%%%%%%% NICOL�S 


%% AHORA SI YA TENEMOS TODAS LAS FUNCIONES :D



%% ahora le metemos numeros a esto (Jugar con lo parametro)
alpha_PS=1; beta_PS=.2; gamma_BP=0.1; % tiene que ser mas grande de  1.3200
alpha_BS=1.1; beta_BS=.1; 
ST=1;
%%
% para graficar las ceroclinas necesitamos un vector de BS y otro para PS
BS_vec=0:0.01:ST;
PS_vec=0:0.01:ST;


%%
figure;
plot(BS_vec, PSss2_fun(BS_vec,alpha_BS, beta_BS, alpha_PS, beta_PS, gamma_BP, ST), 'r'); % función de PS
hold on
plot(BSss2_fun(PS_vec,alpha_BS, beta_BS, alpha_PS, beta_PS, gamma_BP, ST), PS_vec, 'b'); % función de BS
xlim([0,ST+0.1])
ylim([0, ST+0.1])
%%
axis square
xlabel('BS')
ylabel('PS')
legend('PSceroclina','BSceroclina')


%% TENEMOS 4 PUNTOS DE EQUILIBRIO;
BSeq=solBS(alpha_BS, beta_BS, alpha_PS, beta_PS, gamma_BP, ST)
PSeq=solPS(alpha_BS, beta_BS, alpha_PS, beta_PS, gamma_BP, ST)

%%
for jj=1:1:4


% grafica los puntos de equilibrio
BSeq(jj)
PSeq(jj)

% evaluar el jacobiano en el punto de equilibrio
Jac1=Jac(BSeq(jj), PSeq(jj),alpha_BS, beta_BS, alpha_PS, beta_PS, gamma_BP, ST)
eig(Jac1)

if (sum(eig(Jac1)>0) >0)
    disp('unstable')
    
   text(BSeq(jj), PSeq(jj), 'unstable ')
   scatter(BSeq(jj), PSeq(jj), 'k')
else
    disp('stable')
       text(BSeq(jj), PSeq(jj) , 'stable ')
       scatter(BSeq(jj), PSeq(jj), 'k', 'filled')
end


end
%%
line([ST, ST], [0, ST]);
line([0, ST], [ST, ST]);



%% ahora le voy a agregar unas trayectorias dinamicas

% declara la ODE;
Modelo_2_especies=@(t,y,alpha_BS, ST, beta_BS, alpha_PS, beta_PS, gamma_BP)[(alpha_BS*y(1)*(ST-y(1)-y(2)))-(beta_BS*y(1));(alpha_PS*y(2)*(ST-y(1)-y(2))/(1+gamma_BP*y(1)))-(beta_PS*y(2))];



y0=[.01 .01]; %condición inicial
tspan = [0 100]; % intervalo de integración
% 
%seleccionamos nuestro integrador favorito * (euler) (ode15, ode23,..)
[t,y] = ode45(@(t,y)Modelo_2_especies(t,y,alpha_BS, ST, beta_BS, alpha_PS, beta_PS, gamma_BP),tspan,y0);

plot(y(:,1), y(:,2), '--m')


% varias condiciones iniciales
%%
for ii=1:1:100
    
    % generate 2 random numbers
     %{
    if rand<0.25
    y0 = [0 rand];
    elseif rand<0.5
    y0 = [ rand 0];
    elseif rand<0.75
    y0 = [ST rand];
    elseif rand>.75
    %}
    
    % generate 3 random numbers
    r1=rand;
    r2=rand;
    r3=rand;
    
    if r1<0.5    
       if r2<0.5
    y0 = [0.05 r3];
       else
          y0 = [1 r3]; 
       end
    else
        if r2<0.5
            y0 = [r3 0.05];
        else
            
    y0 = [r3 1];
        end
    end
    
   
    

[t,y] = ode45(@(t,y)Modelo_2_especies(t,y,alpha_BS, ST, beta_BS, alpha_PS, beta_PS, gamma_BP),tspan,y0);



 plot(y(:,1),y(:,2),'c')
   
end
 
title(['2 Arboles vs suelo, \alpha_{PS}=' num2str(alpha_PS)  '\beta_{PS}=' num2str(beta_PS)  ' \alpha_{BS}=' num2str(alpha_BS) '\beta_{BS}=' num2str(beta_BS) 'S_T=' num2str(ST)])
legend('off')

%% Diagrama de bifurcaciones 

% vamos a evaluar cómo cambia la estabilidad de los 4 puntos de
% equilibrio cuando variamos un parámetro de bifurcación: beta_TS

%%
syms alpha_BS beta_BS alpha_PS beta_PS gamma_BP ST



solBSSYM =[ -(beta_BS - ST*alpha_BS)/alpha_BS;
                         0;
                                                                0;
    (alpha_BS*beta_PS - alpha_PS*beta_BS)/(alpha_PS*beta_BS*gamma_BP)];


solPSSYM =[0;
    -(beta_PS - ST*alpha_PS)/alpha_PS;
    0;
-(alpha_BS^2*beta_PS - alpha_BS*alpha_PS*beta_BS + alpha_PS*beta_BS^2*gamma_BP - ST*alpha_BS*alpha_PS*beta_BS*gamma_BP)/(alpha_BS*alpha_PS*beta_BS*gamma_BP)];
%%
jj=4
BS=solBSSYM(jj);
PS=solPSSYM(jj)
Jac=[[                    - beta_BS - BS*alpha_BS - alpha_BS*(BS + PS - ST),                                                                          -BS*alpha_BS]
[- PS*alpha_PS*(BS*gamma_BP + 1) - PS*alpha_PS*gamma_BP*(BS + PS - ST), - beta_PS - PS*alpha_PS*(BS*gamma_BP + 1) - alpha_PS*(BS*gamma_BP + 1)*(BS + PS - ST)]]


%%
eig(Jac)
% punto de equilibrio 3 es estable <-> ambos eigenvalores son < 0; vamos a
% ver cuándo beta_TS, nuestro parámetro de bifurcación, genera este cambio
%%

%%%%%%% vamos a ir guardando los egienvalores acà, son 4 pares, un par por
%%%%%%% punto de equilirbio. para que sea estable el punto deequlibrio
%%%%%%% correspondiente pedimos que amobos sean negativso; entnces
%%%%%%% dependiendo del parametro que te intereese puedes igualar los
%%%%%%% eigenvaloers a 0 y resolver / despejar para ese parametro;
eig1=[- beta_PS - alpha_PS*(ST + (beta_BS - ST*alpha_BS)/alpha_BS)*((gamma_BP*(beta_BS - ST*alpha_BS))/alpha_BS - 1);
                                                alpha_BS*(ST + (beta_BS - ST*alpha_BS)/alpha_BS) - ST*alpha_BS];
                                            
                                            % pr ejemplo si nos interesa
                                            % gamma_BP
                                            
                                            
                                            
eig2=[   alpha_BS*(ST + (beta_PS - ST*alpha_PS)/alpha_PS) - beta_BS;
alpha_PS*(ST + (beta_PS - ST*alpha_PS)/alpha_PS) - ST*alpha_PS]                                            
                                            

eig3=[ ST*alpha_BS - beta_BS;
        ST*alpha_PS - beta_PS];
    
    
    eig4=[(((ST^2*alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2*gamma_BP^2 + 2*ST*alpha_BS^4*alpha_PS*beta_BS^2*beta_PS^2*gamma_BP - 2*ST*alpha_BS^4*alpha_PS*beta_BS*beta_PS^3*gamma_BP - 2*ST*alpha_BS^3*alpha_PS^2*beta_BS^3*beta_PS*gamma_BP + 2*ST*alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2*gamma_BP - 4*ST*alpha_BS^2*alpha_PS^2*beta_BS^4*beta_PS*gamma_BP^2 - 2*ST*alpha_BS^2*alpha_PS^2*beta_BS^3*beta_PS^2*gamma_BP^2 + 4*ST*alpha_BS*alpha_PS^3*beta_BS^5*gamma_BP^2 + alpha_BS^5*beta_BS^2*beta_PS^2 - 2*alpha_BS^5*beta_BS*beta_PS^3 + alpha_BS^5*beta_PS^4 - 2*alpha_BS^4*alpha_PS*beta_BS^3*beta_PS + 4*alpha_BS^4*alpha_PS*beta_BS^2*beta_PS^2 - 2*alpha_BS^4*alpha_PS*beta_BS*beta_PS^3 + alpha_BS^3*alpha_PS^2*beta_BS^4 - 2*alpha_BS^3*alpha_PS^2*beta_BS^3*beta_PS + alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2 + 2*alpha_BS^3*alpha_PS*beta_BS^3*beta_PS^2*gamma_BP + 2*alpha_BS^3*alpha_PS*beta_BS^2*beta_PS^3*gamma_BP - 6*alpha_BS^2*alpha_PS^2*beta_BS^4*beta_PS*gamma_BP - 2*alpha_BS^2*alpha_PS^2*beta_BS^3*beta_PS^2*gamma_BP + 4*alpha_BS*alpha_PS^3*beta_BS^5*gamma_BP + 4*alpha_BS*alpha_PS^2*beta_BS^5*beta_PS*gamma_BP^2 + alpha_BS*alpha_PS^2*beta_BS^4*beta_PS^2*gamma_BP^2 - 4*alpha_PS^3*beta_BS^6*gamma_BP^2)/alpha_BS)^(1/2) + alpha_BS^2*beta_PS^2 + alpha_BS*alpha_PS*beta_BS^2 - alpha_BS^2*beta_BS*beta_PS + alpha_PS*beta_BS^2*beta_PS*gamma_BP - alpha_BS*alpha_PS*beta_BS*beta_PS - ST*alpha_BS*alpha_PS*beta_BS*beta_PS*gamma_BP)/(2*alpha_PS*beta_BS^2*gamma_BP);
-(((ST^2*alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2*gamma_BP^2 + 2*ST*alpha_BS^4*alpha_PS*beta_BS^2*beta_PS^2*gamma_BP - 2*ST*alpha_BS^4*alpha_PS*beta_BS*beta_PS^3*gamma_BP - 2*ST*alpha_BS^3*alpha_PS^2*beta_BS^3*beta_PS*gamma_BP + 2*ST*alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2*gamma_BP - 4*ST*alpha_BS^2*alpha_PS^2*beta_BS^4*beta_PS*gamma_BP^2 - 2*ST*alpha_BS^2*alpha_PS^2*beta_BS^3*beta_PS^2*gamma_BP^2 + 4*ST*alpha_BS*alpha_PS^3*beta_BS^5*gamma_BP^2 + alpha_BS^5*beta_BS^2*beta_PS^2 - 2*alpha_BS^5*beta_BS*beta_PS^3 + alpha_BS^5*beta_PS^4 - 2*alpha_BS^4*alpha_PS*beta_BS^3*beta_PS + 4*alpha_BS^4*alpha_PS*beta_BS^2*beta_PS^2 - 2*alpha_BS^4*alpha_PS*beta_BS*beta_PS^3 + alpha_BS^3*alpha_PS^2*beta_BS^4 - 2*alpha_BS^3*alpha_PS^2*beta_BS^3*beta_PS + alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2 + 2*alpha_BS^3*alpha_PS*beta_BS^3*beta_PS^2*gamma_BP + 2*alpha_BS^3*alpha_PS*beta_BS^2*beta_PS^3*gamma_BP - 6*alpha_BS^2*alpha_PS^2*beta_BS^4*beta_PS*gamma_BP - 2*alpha_BS^2*alpha_PS^2*beta_BS^3*beta_PS^2*gamma_BP + 4*alpha_BS*alpha_PS^3*beta_BS^5*gamma_BP + 4*alpha_BS*alpha_PS^2*beta_BS^5*beta_PS*gamma_BP^2 + alpha_BS*alpha_PS^2*beta_BS^4*beta_PS^2*gamma_BP^2 - 4*alpha_PS^3*beta_BS^6*gamma_BP^2)/alpha_BS)^(1/2) - alpha_BS^2*beta_PS^2 - alpha_BS*alpha_PS*beta_BS^2 + alpha_BS^2*beta_BS*beta_PS - alpha_PS*beta_BS^2*beta_PS*gamma_BP + alpha_BS*alpha_PS*beta_BS*beta_PS + ST*alpha_BS*alpha_PS*beta_BS*beta_PS*gamma_BP)/(2*alpha_PS*beta_BS^2*gamma_BP)];


 
    
%%
solve(eig4(1)==0, gamma_BP)
solve(eig4(2)==0, gamma_BP)

%
% el valor critic que obtenemos es
gamma_BP=-(beta_PS*alpha_BS^2 - alpha_PS*beta_BS*alpha_BS)/(alpha_PS*beta_BS^2 - ST*alpha_BS*alpha_PS*beta_BS)
% dado que solve(eig1(2)==0, gamma_BP) no da solucion, esto significa que
% gamma_BP no puede inducir un cambio deestambilidad en este punto de
% equilibrio


%% solve(Eig3_1==0, beta_TS)
% nuevamente declaramos parametros
alpha_PS=1; beta_PS=.2; gamma_BP=0.1; % tiene que ser mas grande de  1.3200
alpha_BS=1.1; beta_BS=.1; 
ST=1;
beta_BS_crit1=(alpha_BS*beta_PS)/alpha_PS;

beta_BS_crit2=ST*alpha_BS;

if beta_BS_crit1< beta_BS_crit2
    
    beta_BS_crit=beta_BS_crit1;
else
     beta_BS_crit=beta_BS_crit2;
end

    

%{
for beta_TS=0:0.01:.5
    
    %calcular los puntos de equilibrio
E1=Eq1(beta_BS);
E2=Eq2(beta_BS);
E3=Eq3(beta_BS);

J1=Jac(beta_BS, E1(1), E1(2));
J2=Jac(beta_BS, E2(1), E2(2));
J3=Jac(beta_BS, E3(1), E3(2));

% valuamos estabilidad
if (sum(eig(J1)>0) >0)
   S1=0;
else
    S1=1;
end

if (sum(eig(J2)>0) >0)
   S2=0;
else
    S2=1;
end

if (sum(eig(J3)>0) >0)
   S3=0;
else
    S3=1;
end


subplot(2,1,1)
xlabel('\beta_{BS}')
ylabel('BS_{ss}')
hold on
if S1==1
scatter(beta_BS, E1(1), 100, 'k', 'filled')
else
    scatter(beta_BS, E1(1), 100, 'k')
end


if S2==1
scatter(beta_BS, E2(1), 'r', 'filled')
else
    scatter(beta_BS, E2(1), 'r')
end


if S3==1
scatter(beta_BS, E3(1), 'b', 'filled')
else
    scatter(beta_BS, E3(1), 'b')
end


subplot(2,1,2)
xlabel('\beta_{BS}')
ylabel('PS_{ss}')
hold on
if S1==1
scatter(beta_BS, E1(2),100, 'k', 'filled')
else
    scatter(beta_BS, E1(2),100, 'k')
end


if S2==1
scatter(beta_BS, E2(2), 'r', 'filled')
else
    scatter(beta_BS, E2(2), 'r')
end


if S3==1
scatter(beta_BS, E3(2), 'b', 'filled')
else
    scatter(beta_BS, E3(2), 'b')
end


end

%%
subplot(2,1,1)
xlim([0, beta_BS])
ylim([0, ST])
axis square
line([ beta_BS_crit  beta_BS_crit], [0, ST])
subplot(2,1,2)
xlim([0, beta_BS])
ylim([0, ST])
axis square
line([ beta_BS_crit  beta_BS_crit], [0, ST])
%% daigrama de bifurcación de 2D, de manera analítica
 
% la estabilidad del punto de equilibio 3;PS=-(beta_TS + 0*alpha_TS - ST*alpha_TS)/alpha_TS
                                        % TS=0];
% depende de los eigenvalores del jacobiano evaluado en ese punto de
% equilibrio
% en particular, pedimos para que sea estable que los dos eigenvalores sean
% menor igual a 0, entonces podemos encontrar la curva de bifurcación, dada por los pares
% [alpha_TS*,beta_TS*] que dividen el plano  [alpha_TS,beta_TS] en una región en la que el punto de equlibio 3
%es estable, y una en la que es inestable.

%%
alpha_bS=0:0.01:1;

beta_BS_crit1=(alpha_BS.*beta_PS)./alpha_PS;
beta_BS_crit2=ST.*alpha_BS;

figure;
area(alpha_BS, beta_BS_crit1)%, 'r')
hold on
plot(alpha_BS, beta_BS_crit1, 'r')

%plot(alpha_TS, beta_TS_crit2, 'b')
% nos quedamos con la curva que sea más pequeña

xlabel('\alpha_{BS}')
ylabel('\beta_{BS}')
axis square
xlim([0, 1])
ylim([0,1])
line([0,1], [0,1], 'color', 'k')

%%



%}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





