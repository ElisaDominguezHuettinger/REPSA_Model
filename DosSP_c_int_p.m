%% Modelo_2_especies que cooperan
close all 
clear all
clc


%%% ESTE CODIGO JALA DE PRINICPIO A FIN %. EDH 24 NOV 2020 %%%

%%% TODA LA PARTE SIMBOLICA SLO SE TIENE QUE HACER UNA VEZ - AHORA YA
%%% COMENTADO; LO MARCO COMO SYM

%% declaramos las variables simbolicas:
%% SYM SECCION
%%%% SACAR EXPRESIONES SIMBÓLICAS PARA LOS PUNTOS DE EQUIIBRIO
%
syms alpha_BS beta_BS alpha_PS beta_PS teta_BP
syms BS PS SL  

% sistema de odes:
dBS=(alpha_BS*BS*SL)-(beta_BS*BS)  %tepozanes

dPS=(alpha_PS*PS*SL*(1+teta_BP*BS))-(beta_PS*PS) %palolocos

dSL=-(alpha_BS*BS*SL+alpha_PS*PS*SL*(1+teta_BP*BS))+(beta_BS*BS+beta_PS*PS) %espacio


%% checar que se cumplan ecuaciones de conservación
dBS+dPS+dSL
%% reducir el sistema a dos dimensines
syms ST
SL=(ST-BS-PS)
dBS=(alpha_BS*BS*SL)-(beta_BS*BS)  ;
dPS=(alpha_PS*PS*SL*(1+teta_BP*BS))-(beta_PS*PS);

%% sacamos los puntos de equilibrio

[solBS,solPS] = solve(dBS== 0, dPS== 0,[BS PS])
%}
%% y los acomodamos dentro de funciones
solBS =@(alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)[ -(beta_BS - ST.*alpha_BS)./alpha_BS;
                                                               0;
                                                               0;
(alpha_BS.*beta_PS - alpha_PS.*beta_BS)./(alpha_PS.*beta_BS.*teta_BP)];

%%
solPS =@(alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)[                                                                                            0;
                                                                                                                        -(beta_PS - ST.*alpha_PS)./alpha_PS;
                                                                                                                                                        0;
-(alpha_BS.^2.*beta_PS - alpha_BS.*alpha_PS.*beta_BS + alpha_PS.*beta_BS.^2.*teta_BP - ST.*alpha_BS.*alpha_PS.*beta_BS.*teta_BP)./(alpha_BS.*alpha_PS.*beta_BS.*teta_BP)];


%% ya desde aca vemos que hay 4 posibes soluciones: una gana y mata a la otra, vice versa, ambas mueren o coexistenica
 %%% aCA VA OTRA PARTE SIMBOLICA PARA SACAR LAS CEROCLINAS; ESTO ES
 %%% IMPORTNTE PARA TENER UNA NOCION GEOMETRICA DE LA ESTABILIDAD Y DEL
 %%% SISTEMA EN GRAL
% 
%% sacamos las ceroclinas (con la sig línea encontramos las soluciones de la eq.

% recordar que las ceroclinas del sistema son las curvas donde esa variable
% es==0; en la sig linea obtenemos la ceroclina de BS. Eso significa que en
% toda esa curva BS esta en equilibrio; pero nota que la curva  es una
% funcion de PS!
sol1=solve(dBS==0,BS)
%% soluciones obtenidad de la linea aterior (copiar y pegar de linea de comando)
BSss1=0 % =sol1(1)
BSss2=-(beta_BS + PS*alpha_BS - ST*alpha_BS)/alpha_BS % = sol1(2)

%% lo mismo para la otra variable
sol2=solve(dPS==0,PS)
%% Soluciones obtenidas de la linea anterior
PSss1=0 % =sol2(1)
PSss2=-(beta_PS + BS*alpha_PS - ST*alpha_PS + BS^2*alpha_PS*teta_BP - BS*ST*alpha_PS*teta_BP)/(alpha_PS + BS*alpha_PS*teta_BP) %sol2(2)
%}
 
 %% vamos a escribir la ceroclinas de PS como una funcion de BS y viceversa: (el  arroba se usa para declarar funciones "anonimas", es decir, en el mismo archivo en el que estamos trabajando
%% notar tambien que como queremos poder meter vectores a esta funcion (vectores de BS por ejemplo), necesitamos mdificar a funcion de la ceroclina ta que acepte
% vectores; en particular, sustitur * por .* y / por ./ y ^ por .^ (el
% punto indica multiplicacion/divicion / a la/  entrada por entrada de vector
%% declarar la ceroclina como una función de PS
BSss2_fun=@(PS,alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP,ST)-(beta_BS + PS.*alpha_BS - ST.*alpha_BS)./alpha_BS
%% declarar la ceroclina como una función de BS

PSss2_fun=@(BS,alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)-(beta_PS + BS.*alpha_PS - ST.*alpha_PS + BS.^2.*alpha_PS*teta_BP - BS.*ST.*alpha_PS*teta_BP)./(alpha_PS + BS.*alpha_PS*teta_BP)


 
 
%%  sacamos el Jacobiano, lo vamos a necesitar más adelamnte para evaluar
% equilibrio
%%% SYM
%
jacobian([dBS, dPS], [BS, PS])
%}
%% recordar que la matriz jacobiana es una linearizacion del sistema; es decir, es el equivalente lineal del sistema no lineal .. esto es valido
% solo cerca de un punto ; entionces esta expresion se deberà evaluar en el
% punto de equilibrio que obtengamos mpas adelante. 
% una vez que evaluamos el jacobiano en el punto de equilbrio, hacemos un
% analisis deestabilidad local, es decir, sacamos los eigenvalores y
% eigenvectrs y pedimos que sean todos nevativos (convergencia asintòtica)
 
Jac=@(BS,PS,alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)[[                  - beta_BS - BS*alpha_BS - alpha_BS*(BS + PS - ST),                                                                        -BS*alpha_BS]
[- PS*alpha_PS*(BS*teta_BP + 1) - PS*alpha_PS*teta_BP*(BS + PS - ST), - beta_PS - PS*alpha_PS*(BS*teta_BP + 1) - alpha_PS*(BS*teta_BP + 1)*(BS + PS - ST)]]

%% AHORA SI YA TENEMOS TODAS LAS FUNCIONES :D



%% ahora le metemos números a esto (Jugar con lo parámetro)
alpha_PS=1; beta_PS=.2; teta_BP=10; % tiene que ser mas grande de  1.3200
alpha_BS=1.1; beta_BS=.1; 
ST=1;
%%
% para graficar las ceroclinas necesitamos un vector de BS y otro para PS
BS_vec=0:0.01:ST;
PS_vec=0:0.01:ST;


%%
figure;
p1=plot(BS_vec, PSss2_fun(BS_vec,alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST), 'r'); % función de PS
hold on
p2=plot(BSss2_fun(PS_vec,alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST), PS_vec, 'b'); % función de BS
xlim([0,ST+0.1])
ylim([0, ST+0.1])
%%
axis square
xlabel('BS')
ylabel('PS')




%% TENEMOS 4 PUNTS DE EQUILIBRIO;
BSeq=solBS(alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)
PSeq=solPS(alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)

%%
for jj=1:1:4


% grafica los puntos de equilibrio
BSeq(jj)
PSeq(jj)

% evaluar el jacobiano en el punto de equilibrio
Jac1=Jac(BSeq(jj), PSeq(jj),alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)
eig(Jac1)

if (sum(eig(Jac1)>0) >0)
    disp('unstable')
    
   text(BSeq(jj), PSeq(jj), 'U')
   scatter(BSeq(jj), PSeq(jj), 'k')
else
    disp('stable')
       text(BSeq(jj), PSeq(jj) , 'S')
       scatter(BSeq(jj), PSeq(jj), 'k', 'filled')
end

end

%%
line([ST, ST], [0, ST]);
line([0, ST], [ST, ST]);



%% ahora le voy a agregar unas trayectorias dinámicas

% declara la ODE;
Modelo_2_especies=@(t,y,alpha_BS, ST, beta_BS, alpha_PS, beta_PS, teta_BP)[(alpha_BS*y(1)*(ST-y(1)-y(2)))-(beta_BS*y(1));(alpha_PS*y(2)*(ST-y(1)-y(2))*(1+teta_BP*y(1)))-(beta_PS*y(2))];



y0=[.01 .01]; %condición inicial
tspan = [0 100]; % intervalo de integración
% 
%seleccionamos nuestro integrador favorito * (euler) (ode15, ode23,..)
[t,y] = ode45(@(t,y)Modelo_2_especies(t,y,alpha_BS, ST, beta_BS, alpha_PS, beta_PS, teta_BP),tspan,y0);

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
    
 
    

[t,y] = ode45(@(t,y)Modelo_2_especies(t,y,alpha_BS, ST, beta_BS, alpha_PS, beta_PS, teta_BP),tspan,y0);



 plot(y(:,1),y(:,2),'k')  % si ponemos 'p' salen estrellitas en las trayectoria
   
end
 
title(['2 Arboles vs suelo, \alpha_{PS}=' num2str(alpha_PS)  '\beta_{PS}=' num2str(beta_PS)  ' \alpha_{BS}=' num2str(alpha_BS) '\beta_{BS}=' num2str(beta_BS) 'S_T=' num2str(ST)])
legend([p1 p2], 'PSceroclina','BSceroclina')
%% Diagrama de bifurcaciones 

% vamos a evaluar cómo cambia la estabilidad de los 4 puntos de
% equilibrio cuando variamos un parámetro de bifurcación: beta_TS

%%
syms alpha_BS beta_BS alpha_PS beta_PS teta_BP ST



solBSSYM =[ -(beta_BS - ST*alpha_BS)/alpha_BS;
                         0;
                                                                0;
    (alpha_BS*beta_PS - alpha_PS*beta_BS)/(alpha_PS*beta_BS*teta_BP)];


solPSSYM =[0;
    -(beta_PS - ST*alpha_PS)/alpha_PS;
    0;
-(alpha_BS^2*beta_PS - alpha_BS*alpha_PS*beta_BS + alpha_PS*beta_BS^2*teta_BP - ST*alpha_BS*alpha_PS*beta_BS*teta_BP)/(alpha_BS*alpha_PS*beta_BS*teta_BP)];

%%


jj=4;

BS=solBSSYM(jj);
PS=solPSSYM(jj)
Jac=[[                    - beta_BS - BS*alpha_BS - alpha_BS*(BS + PS - ST),                                                                          -BS*alpha_BS]
[- PS*alpha_PS*(BS*teta_BP + 1) - PS*alpha_PS*teta_BP*(BS + PS - ST), - beta_PS - PS*alpha_PS*(BS*teta_BP + 1) - alpha_PS*(BS*teta_BP + 1)*(BS + PS - ST)]]



eig(Jac)

% punto de equilibrio 3 es estable <-> ambos eigenvalores son < 0; vamos a
% ver cuando beta_TS, nuestro parámetro de bifurcación, genera este cambio
%%

%%%%%%% vamos a ir guardando los egienvalores acà, son 4 pares, un par por
%%%%%%% punto de equilirbio. para que sea estable el punto deequlibrio
%%%%%%% correspondiente pedimos que amobos sean negativso; entnces
%%%%%%% dependiendo del parametro que te intereese puedes igualar los
%%%%%%% eigenvaloers a 0 y resolver / despejar para ese parametro;
eig1=[     - beta_PS - alpha_PS*(ST + (beta_BS - ST*alpha_BS)/alpha_BS)*((teta_BP*(beta_BS - ST*alpha_BS))/alpha_BS - 1);
                                                alpha_BS*(ST + (beta_BS - ST*alpha_BS)/alpha_BS) - ST*alpha_BS];
                                         
                                            % pr ejemplo si nos interesa
                                            % gamma_BP
                                            
                                            
                                            
eig2=[  alpha_BS*(ST + (beta_PS - ST*alpha_PS)/alpha_PS) - beta_BS;
 alpha_PS*(ST + (beta_PS - ST*alpha_PS)/alpha_PS) - ST*alpha_PS]                                            
                                            

eig3=[  ST*alpha_BS - beta_BS;
 ST*alpha_PS - beta_PS];
    
    
eig4=[ (((ST^2*alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2*teta_BP^2 + 2*ST*alpha_BS^4*alpha_PS*beta_BS^2*beta_PS^2*teta_BP - 2*ST*alpha_BS^4*alpha_PS*beta_BS*beta_PS^3*teta_BP - 2*ST*alpha_BS^3*alpha_PS^2*beta_BS^3*beta_PS*teta_BP + 2*ST*alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2*teta_BP - 4*ST*alpha_BS^2*alpha_PS^2*beta_BS^4*beta_PS*teta_BP^2 - 2*ST*alpha_BS^2*alpha_PS^2*beta_BS^3*beta_PS^2*teta_BP^2 + 4*ST*alpha_BS*alpha_PS^3*beta_BS^5*teta_BP^2 + alpha_BS^5*beta_BS^2*beta_PS^2 - 2*alpha_BS^5*beta_BS*beta_PS^3 + alpha_BS^5*beta_PS^4 - 2*alpha_BS^4*alpha_PS*beta_BS^3*beta_PS + 4*alpha_BS^4*alpha_PS*beta_BS^2*beta_PS^2 - 2*alpha_BS^4*alpha_PS*beta_BS*beta_PS^3 + alpha_BS^3*alpha_PS^2*beta_BS^4 - 2*alpha_BS^3*alpha_PS^2*beta_BS^3*beta_PS + alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2 + 2*alpha_BS^3*alpha_PS*beta_BS^3*beta_PS^2*teta_BP + 2*alpha_BS^3*alpha_PS*beta_BS^2*beta_PS^3*teta_BP - 6*alpha_BS^2*alpha_PS^2*beta_BS^4*beta_PS*teta_BP - 2*alpha_BS^2*alpha_PS^2*beta_BS^3*beta_PS^2*teta_BP + 4*alpha_BS*alpha_PS^3*beta_BS^5*teta_BP + 4*alpha_BS*alpha_PS^2*beta_BS^5*beta_PS*teta_BP^2 + alpha_BS*alpha_PS^2*beta_BS^4*beta_PS^2*teta_BP^2 - 4*alpha_PS^3*beta_BS^6*teta_BP^2)/alpha_BS)^(1/2) + alpha_BS^2*beta_PS^2 + alpha_BS*alpha_PS*beta_BS^2 - alpha_BS^2*beta_BS*beta_PS + alpha_PS*beta_BS^2*beta_PS*teta_BP - alpha_BS*alpha_PS*beta_BS*beta_PS - ST*alpha_BS*alpha_PS*beta_BS*beta_PS*teta_BP)/(2*alpha_PS*beta_BS^2*teta_BP);
 -(((ST^2*alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2*teta_BP^2 + 2*ST*alpha_BS^4*alpha_PS*beta_BS^2*beta_PS^2*teta_BP - 2*ST*alpha_BS^4*alpha_PS*beta_BS*beta_PS^3*teta_BP - 2*ST*alpha_BS^3*alpha_PS^2*beta_BS^3*beta_PS*teta_BP + 2*ST*alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2*teta_BP - 4*ST*alpha_BS^2*alpha_PS^2*beta_BS^4*beta_PS*teta_BP^2 - 2*ST*alpha_BS^2*alpha_PS^2*beta_BS^3*beta_PS^2*teta_BP^2 + 4*ST*alpha_BS*alpha_PS^3*beta_BS^5*teta_BP^2 + alpha_BS^5*beta_BS^2*beta_PS^2 - 2*alpha_BS^5*beta_BS*beta_PS^3 + alpha_BS^5*beta_PS^4 - 2*alpha_BS^4*alpha_PS*beta_BS^3*beta_PS + 4*alpha_BS^4*alpha_PS*beta_BS^2*beta_PS^2 - 2*alpha_BS^4*alpha_PS*beta_BS*beta_PS^3 + alpha_BS^3*alpha_PS^2*beta_BS^4 - 2*alpha_BS^3*alpha_PS^2*beta_BS^3*beta_PS + alpha_BS^3*alpha_PS^2*beta_BS^2*beta_PS^2 + 2*alpha_BS^3*alpha_PS*beta_BS^3*beta_PS^2*teta_BP + 2*alpha_BS^3*alpha_PS*beta_BS^2*beta_PS^3*teta_BP - 6*alpha_BS^2*alpha_PS^2*beta_BS^4*beta_PS*teta_BP - 2*alpha_BS^2*alpha_PS^2*beta_BS^3*beta_PS^2*teta_BP + 4*alpha_BS*alpha_PS^3*beta_BS^5*teta_BP + 4*alpha_BS*alpha_PS^2*beta_BS^5*beta_PS*teta_BP^2 + alpha_BS*alpha_PS^2*beta_BS^4*beta_PS^2*teta_BP^2 - 4*alpha_PS^3*beta_BS^6*teta_BP^2)/alpha_BS)^(1/2) - alpha_BS^2*beta_PS^2 - alpha_BS*alpha_PS*beta_BS^2 + alpha_BS^2*beta_BS*beta_PS - alpha_PS*beta_BS^2*beta_PS*teta_BP + alpha_BS*alpha_PS*beta_BS*beta_PS + ST*alpha_BS*alpha_PS*beta_BS*beta_PS*teta_BP)/(2*alpha_PS*beta_BS^2*teta_BP)]
 
    
%% dado que el  punto de equilibrio 4 es el de coexistencia vemos en qu� momento este pierde estabilidad
solve(eig4(1)==0, teta_BP)
%%
solve(eig4(2)==0, teta_BP)

%%

teta_BP_crit=-(beta_PS*alpha_BS^2 - alpha_PS*beta_BS*alpha_BS)/(alpha_PS*beta_BS^2 - ST*alpha_BS*alpha_PS*beta_BS)

Jac=@(BS,PS,alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)[[                  - beta_BS - BS*alpha_BS - alpha_BS*(BS + PS - ST),                                                                        -BS*alpha_BS]
[- PS*alpha_PS*(BS*teta_BP + 1) - PS*alpha_PS*teta_BP*(BS + PS - ST), - beta_PS - PS*alpha_PS*(BS*teta_BP + 1) - alpha_PS*(BS*teta_BP + 1)*(BS + PS - ST)]]

%%
%% ahora le metemos numeros a esto (Jugar con lo parámetro)
alpha_PS=1; beta_PS=.2; teta_BP=10; % tiene que ser mas grande de  1.3200
alpha_BS=1.1; beta_BS=.1; 
ST=1;


teta_BP_crit=-(beta_PS*alpha_BS^2 - alpha_PS*beta_BS*alpha_BS)/(alpha_PS*beta_BS^2 - ST*alpha_BS*alpha_PS*beta_BS)
%% 
figure;
for teta_BP=0.5:0.1:teta_BP_crit+1
    
    %calcular los puntos de equilibrio
BSss=solBS(alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)%%
PSss=solPS(alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST)                                                                                            
     

J1=Jac(BSss(1), PSss(1),alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST);
J2=Jac(BSss(2), PSss(2),alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST);
J3=Jac(BSss(3), PSss(3),alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST);
J4=Jac(BSss(4), PSss(4),alpha_BS, beta_BS, alpha_PS, beta_PS, teta_BP, ST);
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

if (sum(eig(J4)>0) >0)
   S4=0;
else
    S4=1;
end


subplot(2,1,1)
xlabel('\theta_{BS}')
ylabel('BS_{ss}')
hold on
if S1==1
scatter(teta_BP, BSss(1), 100, 'ok', 'filled')
else
    scatter(teta_BP, BSss(1), 100, 'ok')
end

if S2==1
scatter(teta_BP, BSSss(2), 100, 'db', 'filled')
else
    scatter(teta_BP, BSss(2), 100, 'db')
end


if S3==1
scatter(teta_BP, BSss(3), 100, 'sc', 'filled')
else
    scatter(teta_BP, BSss(3), 100, 'sc')
end


if S4==1
scatter(teta_BP, BSss(4), 100, 'pg', 'filled')
else
    scatter(teta_BP, BSss(4), 100, 'pg')
end


subplot(2,1,2)
xlabel('\theta_{BS}')
ylabel('PS_{ss}')
hold on
if S1==1
scatter(teta_BP, PSss(1), 100, 'ok', 'filled')
else
    scatter(teta_BP, PSss(1), 100, 'ok')
end

if S2==1
scatter(teta_BP, PSss(2), 100, 'db', 'filled')
else
    scatter(teta_BP, PSss(2), 100, 'db')
end


if S3==1
scatter(teta_BP, PSss(3), 100, 'sr', 'filled')
else
    scatter(teta_BP, PSss(3), 100, 'sr')
end


if S4==1
scatter(teta_BP, PSss(4), 100, 'pg', 'filled')
else
    scatter(teta_BP, PSss(4), 100, 'pg')
end


end

%%
subplot(2,1,1)
ylim([0, ST])
axis square
line([ teta_BP_crit  teta_BP_crit], [0, ST])

subplot(2,1,2)
ylim([0, ST])
axis square
line([ teta_BP_crit  teta_BP_crit], [0, ST])


suptitle('black 1 BS only , blue 2 PS only, red 3 00, green 4 coex; filled stable')
hold on

