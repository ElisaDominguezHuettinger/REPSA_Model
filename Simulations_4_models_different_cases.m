close all 
clear all
clc

dir=fullfile('Simulations');
if ~isfolder(dir)
    mkdir(dir)
end

dir1=fullfile('Simulations','Model1');
if ~isfolder(dir1)
    mkdir(dir1)
end
dir2=fullfile('Simulations','Model2');
if ~isfolder(dir2)
    mkdir(dir2)
end
dir3=fullfile('Simulations','Model3');
if ~isfolder(dir3)
    mkdir(dir3)
end
dir4=fullfile('Simulations','Model4');
if ~isfolder(dir4)
    mkdir(dir4)
end

tspan=[0 10];
for i=1:2
    for j=1:4
        if i==1
            p=[4 2 4 1 1];
        elseif i==2
            if j==3
                p=[2 4 2 5 1];
            else
                p=[2 5 2 4 1];
            end
        end
        file_dir=fullfile('Simulations',['Model', num2str(i)]);
        file_name=['Case',num2str(j)];
        
        switch j                %i:Número modelo
            case 1              %j: Numero caso
                x0=[0.1;0.01;0.1];  %Puede modificar x0 para estudiar distintos comportamientos
            case 2
                x0=[0.01;1;0.5];
            case 3
                x0=[2;0.01;0.7];
            case 4
                x0=[1.5;2;0.5];
        end
        
        plot_case(i,j,p,tspan,x0,file_name,file_dir)
    end
end

for i=3:4
    for j=1:8
        if i==3
            if j==3
                p=[5; 1; 2; 1; 5; 1; 1; 1; 1; 5; 2];
            elseif j==5
                p=[5; 1; 3; 1; 2; 1; 1; 2; 2; 5; 5];
            else
                p=[5; 3; 7; 2; 10; 3; 3; 5; 2; 7; 2];
            end
        else
            if j==3
                p=[5; 1; 2; 1; 5; 1; 1; 1; 1; 5; 2];
            elseif j==5
                p=[5; 2; 3; 3; 2; 1; 2; 3; 1; 2; 3];
            else
                p=[5; 2; 3; 2; 2; 1; 3; 3; 2; 3; 3];
            end
        end
        file_dir=fullfile('Simulations',['Model', num2str(i)]);
        file_name=['Case',num2str(j)];
        
        switch j
            case 1
                x0=[0.1;0.1;0.1;0.1];       %Puede modificar x0 para estudiar distintos comportamientos
            case 2
                x0=[0.01;0.01;0.01;1];
            case 3
                x0=[0.1;2;0.1;0.4];
            case 4
                x0=[0.1; 2; 5; 10];
            case 5
                x0=[0.1;0.03;0.02;2];
            case 6
                x0=[3;0.02;5;8];
            case 7
                x0=[1;15;0.01;0.7];
            case 8
                x0=[10;15;17;14];
        end
        plot_case(i,j,p,tspan,x0,file_name,file_dir)
    end
end




function plot_case(model,model_case,p,tspan,x0,file_name,file_dir)
if model==1
    par=['alpha_1=',num2str(p(1)), '; beta_1=',num2str(p(2)), '; alpha_2=',num2str(p(3)), '; beta_2=',num2str(p(4)), '; gamma=',num2str(p(5))];
elseif model==2
    par=['alpha_1=',num2str(p(1)), '; beta_1=',num2str(p(2)), '; alpha_2=',num2str(p(3)), '; beta_2=',num2str(p(4)), '; theta=',num2str(p(5))];
elseif model==3
    par=['alpha_1=',num2str(p(1)), '; gamma_1=',num2str(p(2)), '; beta_1=',num2str(p(3)), '; mu_1=',num2str(p(4)), '; alpha_2=',num2str(p(5)), '; gamma_2=',num2str(p(6)), '; gamma_3=',num2str(p(7)), '; beta_2=',num2str(p(8)), '; mu_2=',num2str(p(9)), '; alpha_3=',num2str(p(10)), '; beta_3=',num2str(p(11))];
else
    par=['alpha_1=',num2str(p(1)), '; gamma_1=',num2str(p(2)), '; beta_1=',num2str(p(3)), '; mu_1=',num2str(p(4)), '; alpha_2=',num2str(p(5)), '; theta=',num2str(p(6)), '; gamma_2=',num2str(p(7)), '; beta_2=',num2str(p(8)), '; mu_2=',num2str(p(9)), '; alpha_3=',num2str(p(10)), '; beta_3=',num2str(p(11))];
end

[eq, x_free, stab_dir]=get_equilibrium(model,model_case,p);
S=set_simulation(model,p,tspan,x0);
t=S(1,:);
y=S(2:end,:);

figure('Name',['Model ',num2str(model),', case ',num2str(model_case)])

if model<3
    for i=1:3
        subplot(1,3,i)
        plot(t, y(i,:), 'k','LineWidth',2)
        hold on
        if x_free==i
            if length(stab_dir)<length(eq(1,:))
                if stab_dir==1
                    plot([t(end) t(end)], [eq(i,2) max(y(i,:))+eq(i,2)], '--b', 'LineWidth',2)
                    plot([t(end) t(end)], [0 eq(i,1)], '--r', 'LineWidth',2)
                elseif stab_dir==-1
                    plot([t(end) t(end)], [eq(i,2) max(y(i,:))+eq(i,2)], '--r', 'LineWidth',2)
                    plot([t(end) t(end)], [0 eq(i,1)], '--b', 'LineWidth',2)
                end
            else
                for j=1:length(eq(1,:))
                    if stab_dir(j)==1
                        plot([t(end) t(end)], [eq(i,j) max(y(i,:))+eq(i,j)], '--b', 'LineWidth',2)
                    elseif stab_dir(j)==-1
                        plot([t(end) t(end)], [0 eq(i,j)], '--b', 'LineWidth',2)
                    end
                end
            end
        else
            if stab_dir~=0
                plot([0 t(end)], [eq(i) eq(i)], '--g', 'LineWidth',2)
            end
        end
        xlim([0 t(end)]);
        ylabel(['x_',num2str(i)])
        xlabel('Time')
        hold on
        axis square
    end
else
    for i=1:4
        subplot(2,2,i)
        plot(t, y(i,:), 'k','LineWidth',2)
        hold on
        if x_free==i
            if length(stab_dir)<length(eq(1,:))
                if stab_dir==1
                    plot([t(end) t(end)], [eq(i,2) max(y(i,:))+eq(i,2)], '--b', 'LineWidth',2)
                    plot([t(end) t(end)], [0 eq(i,1)], '--r', 'LineWidth',2)
                elseif stab_dir==-1
                    plot([t(end) t(end)], [eq(i,2) max(y(i,:))+eq(i,2)], '--r', 'LineWidth',2)
                    plot([t(end) t(end)], [0 eq(i,1)], '--b', 'LineWidth',2)
                end
            else
                for j=1:length(eq(1,:))
                    if stab_dir(j)==1
                        plot([t(end) t(end)], [eq(i,j) max(y(i,:))], '--b', 'LineWidth',2)
                    elseif stab_dir(j)==-1
                        plot([t(end) t(end)], [0 eq(i,j)], '--b', 'LineWidth',2)
                    end
                end
            end
        else
            if stab_dir~=0
                plot([0 t(end)], [eq(i) eq(i)], '--g', 'LineWidth',2)
            end
        end
        xlim([0 t(end)]);
        ylabel(['x_',num2str(i)])
        xlabel('Time')
        hold on
        axis square
    end
end
h=suptitle({['Simulation model ',num2str(model), ', case ',num2str(model_case)], par});
set(h,'FontSize',8,'FontWeight','normal')

f=fullfile(file_dir,file_name);
print(gcf,f,'-dpdf','-r0')

close
end



function S=set_simulation(model,p,tspan,x0)
switch model
    case 1
        [t,y]=ode45(@(t,y)model1(t,y,p),tspan,x0);
    case 2
        [t,y]=ode45(@(t,y)model2(t,y,p),tspan,x0);
    case 3
        [t,y]=ode45(@(t,y)model3(t,y,p),tspan,x0);
    case 4
        [t,y]=ode45(@(t,y)model4(t,y,p),tspan,x0);
    otherwise
        warning('Variable "model" out of range.')
end
S=[t';y'];
end


function [eq,x_free,stab_dir]=get_equilibrium(model,model_case,p)
x=0;
switch model
    case 1
        alpha1=p(1); beta1=p(2); alpha2=p(3); beta2=p(4); gamma=p(5);
        switch model_case
            case 1
                x1=min([beta1/alpha1,beta2/alpha2]);
                x2=max([beta1/alpha1,beta2/alpha2]);
                eq=[0 0;0 0;x1 x2];
                x_free=3;
                stab_dir=-1;
            case 2
                eq=[0;0;beta2/alpha2];
                x_free=2;
                stab_dir=1;
            case 3
                x=max([0,(alpha2-alpha1)/(alpha1*gamma)]);
                eq=[x;0;beta1/alpha1];
                x_free=1;
                stab_dir=1;
            case 4
                eq=[(alpha2*beta1-alpha1*beta2)/(alpha1*beta2*gamma);0;beta1/alpha1];
                x_free=2;
                stab_dir=1;
            otherwise
                warning('Variable "model_case" out of range.')
        end
    case 2
        alpha1=p(1); beta1=p(2); alpha2=p(3); beta2=p(4); teta=p(5);
        switch model_case
            case 1
                x1=min([beta1/alpha1,beta2/alpha2]);
                x2=max([beta1/alpha1,beta2/alpha2]);
                eq=[0 0;0 0;x1 x2];
                x_free=3;
                stab_dir=-1;
            case 2
                x=max([beta1/alpha1,beta2/alpha2]);
                eq=[0;x;beta2/alpha2];
                x_free=2;
                stab_dir=1;
            case 3
                x=(alpha1*beta2-alpha2*beta1)/(alpha1*teta);
                eq=[x;0;beta1/alpha1];
                x_free=1;
                stab_dir=-1;
            case 4
                G=(alpha2*beta1-alpha1*beta2)*(2+((alpha1^2)*beta2)/(alpha2*(beta1^2)*teta));
                H=(alpha1^2)*(alpha2*beta1-alpha1*beta2)/((alpha2*beta1*teta)^2);
                x=(beta1^2)/((alpha1^2)*(beta2^2))*(G+sqrt((G^2)-H*(alpha1^2)*(beta2^2)/(beta1^2)));
                eq=[(alpha1*beta2-alpha2*beta1)/(alpha2*beta1*teta);x;beta1/alpha1];
                x_free=2;
                stab_dir=-1;
            otherwise
                warning('Variable "model_case" out of range.')
        end
    case 3
        alpha1=p(1); gamma1=p(2); beta1=p(3); mu1=p(4);  alpha2=p(5); 
        gamma2=p(6); gamma3=p(7); beta2=p(8); mu2=p(9);  alpha3=p(10); beta3=p(11);
        
        M2=(1/(2*mu2*gamma3))*(-(mu2+beta2*gamma3)+sqrt((mu2+beta2*gamma3)^2+4*alpha1*beta3*mu2*gamma3/alpha3));
        M1=(1/(2*mu1*gamma1))*(-(mu1+beta1*gamma1)+sqrt((mu1+beta1*gamma1)^2+4*alpha1*beta3*mu1*gamma1/alpha3));
        Q=(alpha2*beta3-alpha3*(beta2+mu2*M1)*(1+gamma3*M1))/(alpha3*gamma2*(beta2+mu2*M1));
        switch model_case
            case 1
                x1=min([beta1/alpha1,beta2/alpha2,beta3/alpha3]);
                x2=max([beta1/alpha1,beta2/alpha2,beta3/alpha3]);
                eq=[0 0;0 0;0 0;x1 x2];
                x_free=4;
                stab_dir=-1;
            case 2
                x=max([M1,M2,0]);
                eq=[0;0;x;beta3/alpha3];
                x_free=3;
                stab_dir=1;
            case 3
                eq=[0;0;0;beta2/alpha2];
                x_free=2;
                stab_dir=1;
            case 4
                eq=[0;x; M2;beta3/alpha3];
                x_free=2;
                stab_dir=0;
            case 5
                x=max([0,(alpha2*beta1-alpha1*beta2)/(alpha1*beta2*gamma2)]);
                eq=[x;0;0;beta1/alpha1];
                x_free=1;
                stab_dir=1;
            case 6
                eq=[x;0;M1;beta3/alpha3];
                x_free=1;
                stab_dir=0;
            case 7
                eq=[(alpha2*beta1-alpha1*beta2)/(alpha1*beta2*gamma2);0;0;beta1/alpha1];
                x_free=2;
                stab_dir=1;
            case 8
                eq=[Q;x;M1;beta3/alpha3];
                x_free=2;
                stab_dir=0;
            otherwise
                warning('Variable "model_case" out of range.')
        end
    case 4
        alpha1=p(1); gamma1=p(2); beta1=p(3); mu1=p(4);  alpha2=p(5); 
        teta=p(6); gamma2=p(7); beta2=p(8); mu2=p(9);  alpha3=p(10); beta3=p(11);
        
        M2=(1/(2*mu2*gamma2))*(-(mu2+beta2*gamma2)+sqrt((mu2+beta2*gamma2)^2+4*alpha2*beta3*mu2*gamma2/alpha3));
        M1=(1/(2*mu1*gamma1))*(-(mu1+beta1*gamma1)+sqrt((mu1+beta1*gamma1)^2+4*alpha1*beta3*mu1*gamma1/alpha3));
        K=(alpha3*beta2-alpha2*beta3+4*alpha1*beta3*gamma1*mu1)/alpha3+2*((beta1^2)*(gamma1^2)+mu1^2);
        M=beta2*gamma2+mu2-2*gamma2*mu2*(beta1*gamma1+mu1);
        N=(beta1*gamma1-mu1)^2+4*alpha1*beta3*gamma1*mu1/alpha3;
        Q=K+M*sqrt(N);
        
        switch model_case
            case 1
                x1=min([beta1/alpha1,beta2/alpha2,beta3/alpha3]);
                x2=max([beta1/alpha1,beta2/alpha2,beta3/alpha3]);
                eq=[0 0;0 0;0 0;x1 x2];
                x_free=4;
                stab_dir=-1;
            case 2
                x=max([0,M1,M2]);
                eq=[0;0;x;beta3/alpha3];
                x_free=3;
                stab_dir=1;
            case 3
                eq=[0;0;0;beta2/alpha2];
                x_free=2;
                stab_dir=1;
            case 4
                eq=[0;x;M2;beta3/alpha3];
                x_free=2;
                stab_dir=0;
            case 5
                x=(alpha1*beta2-alpha2*beta1)/(alpha2*beta1*teta);
                eq=[x;0;0;beta1/alpha1];
                x_free=1;
                stab_dir=-1;
            case 6
                eq=[x;0;M1;beta3/alpha3];
                x_free=1;
                stab_dir=0;
            case 7
                P=((beta1^2)/((alpha1^2)*(beta2^2)))*(-(2+(((alpha1^2)*beta2)/(alpha2*(beta1^2)*teta)))*(alpha2*beta1-alpha1*beta2)+2*sqrt((1+((alpha1^2)*beta2)/(alpha2*(beta1^2)*teta))*((alpha2*beta1-alpha1*beta2)^2)));
                x1=P;
                x2=max([P,(alpha1^2)*(alpha1*beta2-alpha2*beta1)/(alpha2*teta*((alpha1^2)*beta2+alpha2*(beta1^2)*teta))]);
                eq=[(alpha1*beta2-alpha2*beta1)/(alpha2*beta1*teta), (alpha1*beta2-alpha2*beta1)/(alpha2*beta1*teta);x1, x2;0, 0;beta1/alpha1, beta1/alpha1];
                x_free=2;
                stab_dir=[-1, 1];
            case 8
                eq=[Q;x;M1;beta3/alpha3];
                x_free=2;
                stab_dir=0;
            otherwise
                warning('Variable "model_case" out of range.')
        end
    otherwise
        warning('Variable "model" out of range.')
end
end


function dy=model1(t,y,p)
BS=y(1);
PS=y(2);
SL=y(3);

alpha_BS=p(1); beta_BS=p(2); alpha_PS=p(3); beta_PS=p(4); gamma_BP=p(5);

dBS=(alpha_BS*BS*SL)-(beta_BS*BS);  %tepozanes
dPS=(alpha_PS*PS*SL/(1+gamma_BP*BS))-(beta_PS*PS); %palolocos
dSL=-(alpha_BS*BS*SL+alpha_PS*PS*SL/(1+gamma_BP*BS))+(beta_BS*BS+beta_PS*PS); %espacio

dy=[dBS; dPS; dSL];
end


function dy=model2(t,y,p)
BS=y(1);
PS=y(2);
SL=y(3);

alpha_BS=p(1); beta_BS=p(2); alpha_PS=p(3); beta_PS=p(4); teta_BP=p(5);

dBS=(alpha_BS*BS*SL)-(beta_BS*BS);  %tepozanes
dPS=(alpha_PS*PS*SL*(1+teta_BP*BS))-(beta_PS*PS); %palolocos
dSL=-(alpha_BS*BS*SL+alpha_PS*PS*SL*(1+teta_BP*BS))+(beta_BS*BS+beta_PS*PS); %espacio

dy=[dBS; dPS; dSL];
end


function dy=model3(t,y,p)
BS=y(1);
PS=y(2);
ES=y(3);
SL=y(4);

alpha_BS=p(1); gamma_EB=p(2); beta_BS=p(3); mu_EB=p(4);  alpha_PS=p(5); 
gamma_BP=p(6); gamma_EP=p(7); beta_PS=p(8); mu_EP=p(9);  alpha_ES=p(10); beta_ES=p(11);

dBS=((alpha_BS*BS*SL)/(1+(ES*gamma_EB)))-(BS*(beta_BS+(mu_EB*ES))) ; %Tepozanes
dPS=((alpha_PS*PS*SL)/((1+(BS*gamma_BP))+(ES*gamma_EP)))   -(PS*(beta_PS+(mu_EP*ES))) ; %Palo locos
dES=(alpha_ES*ES*SL)-(beta_ES*ES); %Eucaliptos
dSL= -((alpha_BS*BS*SL)/(1+(ES*gamma_EB)))+(BS*(beta_BS+(mu_EB*ES))) -((alpha_PS*PS*SL)/((1+(BS*gamma_BP))+(ES*gamma_EP))) +(PS*(beta_PS+(mu_EP*ES))) -(alpha_ES*ES*SL)+(beta_ES*ES); %espacio disponible

dy=[dBS; dPS; dES; dSL];
end


function dy=model4(t,y,p)
BS=y(1);
PS=y(2);
ES=y(3);
SL=y(4);

alpha_BS=p(1); gamma_EB=p(2); beta_BS=p(3); mu_EB=p(4);  alpha_PS=p(5); 
teta_BP=p(6); gamma_EP=p(7); beta_PS=p(8); mu_EP=p(9);  alpha_ES=p(10); beta_ES=p(11);

dBS=((alpha_BS*BS*SL)/(1+(ES*gamma_EB)))-(BS*(beta_BS+(mu_EB*ES))) ; %Tepozanes
dPS=(((alpha_PS*PS*SL)*(1+teta_BP*BS))/(1+(ES*gamma_EP)))-(PS*(beta_PS+(mu_EP*ES))) ; %Palo locos
dES=(alpha_ES*ES*SL)-(beta_ES*ES); %Eucaliptos
dSL=-((alpha_BS*BS*SL)/(1+(ES*gamma_EB)))+(BS*(beta_BS+(mu_EB*ES)))-(((alpha_PS*PS*SL)*(1+teta_BP*BS))/(1+(ES*gamma_EP)))+(PS*(beta_PS+(mu_EP*ES))) -(alpha_ES*ES*SL)+(beta_ES*ES); %espacio disponible

dy=[dBS; dPS; dES; dSL];
end