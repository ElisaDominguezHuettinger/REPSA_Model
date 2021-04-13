close all 
clear all
clc

f=fullfile('Simulations','AtractingRegions');
if ~isfolder(f)
    mkdir(f)
end

%%Parameters
% alpha_BS=.00137;
% gamma_EB=.01;
% beta_BS=.07;
% mu_EB=.01;
% alpha_PS=.0011;
% teta_BP=.0995;
% gamma_EP=.01;
% beta_PS=.29;
% mu_EP=.01;
% alpha_ES=.001; 
% beta_ES=.1;
% ST=146.9;
% f=fullfile('Simulations','AtractingRegions','InitialCase');
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
%{
Optim=load('OptimModel4.mat');
p=Optim.p;
ST=Optim.ST;
alpha_BS=p(1); gamma_EB=p(2); beta_BS=p(3); mu_EB=p(4);  alpha_PS=p(5); 
teta_BP=p(6); gamma_EP=p(7); beta_PS=p(8); mu_EP=p(9);  alpha_ES=p(10); 
beta_ES=p(11);
%}
f=fullfile('Simulations','AtractingRegions','OptimCase');

n=50;   %Cantidad muestras por dimensión
%%

p=[alpha_BS;gamma_EB;beta_BS;mu_EB;alpha_PS;teta_BP;gamma_EP;beta_PS;mu_EP;alpha_ES;beta_ES];

[eq,ind]=get_equilibrium(p,ST);

c={'r','g','b','c','m','y','k'};
str='';
for i=1:length(ind)
    switch i
        case 1
            str=[str,'case ',num2str(ind(i)),': red'];
        case 2
            str=[str,'; case ',num2str(ind(i)),': green'];
        case 3
            str=[str,'; case ',num2str(ind(i)),': blue'];
        case 4
            str=[str,'; case ',num2str(ind(i)),': cyan'];
        case 5
            str=[str,'; case ',num2str(ind(i)),': magenta'];
        case 6
            str=[str,'; case ',num2str(ind(i)),': yellow'];
        case 7
            str=[str,'; case ',num2str(ind(i)),': black'];
    end
end
figure('Name','Atracting Regions')
for i=1:n
    j=1;
    while j<(n-i)
        m=sqrt(3)*ST;
        tspan=[0, 200];
        x0=[ST*i/n;ST*j/n;ST-ST*(i+j)/n];
        [~,y]=ode45(@(t,y)model4(t,y,p,ST),tspan,x0);
        y=y';
        for k=1:length(eq(1,:))
            if norm(y(:,end)-eq(:,k))<m
                m=norm(y(:,end)-eq(:,k));
                color=c(k);
                index=ind(k);
            end
        end
        ls=strcat('-',color);
        plot3(y(1,:),y(2,:),y(3,:),ls{1})
        hold on
        j=j+1;
    end
end
for k=1:length(eq(1,:))
    color=c(k);
    ls=strcat('o',color);
    plot3(eq(1,k),eq(2,k),eq(3,k),ls{1})
end
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')

par={['alpha_B_S=',num2str(p(1)), '; gamma_E_B=',num2str(p(2)), '; beta_B_S=',num2str(p(3)), '; mu_E_B=',num2str(p(4)), '; alpha_P_S=',num2str(p(5)), '; theta_B_P=',num2str(p(6)),';'], ['gamma_E_P=',num2str(p(7)), '; beta_P_S=',num2str(p(8)), '; mu_E_P=',num2str(p(9)), '; alpha_E_S=',num2str(p(10)), '; beta_E_S=',num2str(p(11)), '; ST=',num2str(ST)],[str]};
h=suptitle({['Simulation Model 4: Optim Parameters '], [par{1}],[par{2}],[par{3}]});
set(h,'FontSize',8,'FontWeight','normal')


print(gcf,f,'-dpdf','-r0')







function [eq,ind]=get_equilibrium(p,ST)
alpha1=p(1); gamma1=p(2); beta1=p(3); mu1=p(4);  alpha2=p(5); 
teta=p(6); gamma2=p(7); beta2=p(8); mu2=p(9);  alpha3=p(10); beta3=p(11);

M2=(1/(2*mu2*gamma2))*(-(mu2+beta2*gamma2)+sqrt((mu2+beta2*gamma2)^2+4*alpha2*beta3*mu2*gamma2/alpha3));
M1=(1/(2*mu1*gamma1))*(-(mu1+beta1*gamma1)+sqrt((mu1+beta1*gamma1)^2+4*alpha1*beta3*mu1*gamma1/alpha3));
K=(alpha3*beta2-alpha2*beta3+4*alpha1*beta3*gamma1*mu1)/alpha3+2*((beta1^2)*(gamma1^2)+mu1^2);
M=beta2*gamma2+mu2-2*gamma2*mu2*(beta1*gamma1+mu1);
N=(beta1*gamma1-mu1)^2+4*alpha1*beta3*gamma1*mu1/alpha3;
Q=K+M*sqrt(N);

eq=[];
ind=[];

% eq1=[0;0;0;ST];
eq1=[0;0;0];

x=ST-beta3/alpha3;
% eq2=[0;0;x;beta3/alpha3];
eq2=[0;0;x];

x=ST-beta2/alpha2;
% eq3=[0,x,0,beta2/alpha2];
eq3=[0;x;0];

x=ST-beta1/alpha1;
% eq5=[x,0,0,beta1/alpha1];
eq5=[x;0;0];

x=ST-M1-beta3/alpha3;
% eq6=[x;0;M1;beta3/alpha3];
eq6=[x;0;M1];

x=ST-(alpha1*beta2-alpha2*beta1)/(alpha2*beta1*teta)-beta1/alpha1;
% eq7=[(alpha1*beta2-alpha2*beta1)/(alpha2*beta1*teta);x;0;beta1/alpha1];
eq7=[(alpha1*beta2-alpha2*beta1)/(alpha2*beta1*teta);x;0];

x=ST-Q-M1-beta3/alpha3;
% eq8=[Q;x;M1;beta3/alpha3];
eq8=[Q;x;M1];

x=ST-beta3/alpha3;
if ST<min([beta1/alpha1,beta2/alpha2,beta3/alpha3]) && x>=0
    eq=[eq,eq1];
    ind=[ind,1];
end
x=ST-beta3/alpha3;
if ST-beta3/alpha3>max([M1,M2]) && x>0
    eq=[eq,eq2];
    ind=[ind,2];
end
x=ST-beta2/alpha2;
if alpha1*beta2<alpha2*beta1 && alpha3*beta2<alpha2*beta3 && x>0
    eq=[eq,eq3];
    ind=[ind,3];
end
x=ST-beta1/alpha1;
if alpha2*beta1<alpha1*beta2 && (ST-beta1/alpha1)<(alpha1*beta2-alpha2*beta1)/(alpha2*beta1*teta) && x>0
    eq=[eq,eq5];
    ind=[ind,5];
end
J22=alpha2*beta3*(1+teta*(ST-M1-beta3/alpha3))/(alpha3*(1+gamma2*M1))-(beta2+mu2*M1);
J44=-alpha3*M1-alpha1*(ST-M1-beta3/alpha3)/(1+gamma1*M1);
J42=beta2+mu2*M1-alpha2*beta3*(1+teta*(ST-M1-beta3/alpha3))/(alpha3*(1+gamma2*M1));
J43=(mu2+mu1+alpha1*beta3*gamma1/(alpha3*((1+gamma1*M1)^2)))*(ST-M1-beta3/alpha3);
lambda=roots([1, -(J22+J44), J22*J44-alpha3*M1*J43, alpha3*M1*(J22*J43+mu2*(ST-M1-beta3/alpha3)*J42)]);
x=ST-M1-beta3/alpha3;
if x>0 && M1>0 && real(lambda(1))<0 && real(lambda(2))<0 && real(lambda(3))<0
    eq=[eq,eq6];
    ind=[ind,6];
end
P=((beta1^2)/((alpha1^2)*(beta2^2)))*(-(2+(((alpha1^2)*beta2)/(alpha2*(beta1^2)*teta)))*(alpha2*beta1-alpha1*beta2)+2*sqrt((1+((alpha1^2)*beta2)/(alpha2*(beta1^2)*teta))*((alpha2*beta1-alpha1*beta2)^2)));
x=ST-(alpha1*beta2-alpha2*beta1)/(alpha2*beta1*teta)-beta1/alpha1;
if x>0 && alpha2*beta1<alpha1*beta2 && alpha3*beta1<alpha1*beta3 && ((x<P && x>0) || x>max([P, (alpha1^2)*(alpha1*beta2-alpha2*beta1)/(alpha2*teta*((alpha1^2)*beta2+alpha2*(beta1^2)*teta))]))
    eq=[eq,eq7];
    ind=[ind,7];
end
x=ST-Q-M1-beta3/alpha3;
J13=-(mu1+gamma1*(beta1+mu1*M1)/(1+gamma1*M1))*Q;
J14=Q*(beta1+mu1*M1)*alpha3/beta3;
J21=alpha2*beta3*teta*x/(alpha3*(1+gamma2*M1));
J43=Q*(mu2+mu1+gamma1*(beta1+mu1*M1)/(1+gamma1*M1))+gamma2*x*(beta2+mu2*M1)/(1+gamma2*M1);
J44=-alpha3*M1-(beta2+mu2*M1)*alpha3*x/beta2-alpha3*Q*(beta1+mu1*M1)/beta3;
lambda=roots([1, -J44, J14*J21-alpha3*M1*J43, alpha3*M1*J13*J21]);
if alpha3*beta1<alpha1*beta3 && Q>0 && M1>0 && x>0 && real(lambda(1))<0 && real(lambda(2))<0 && real(lambda(3))
    eq=[eq,eq8];
    ind=[ind,8];
end
end



function dy=model4(t,y,p,ST)
BS=y(1);
PS=y(2);
ES=y(3);

alpha_BS=p(1); gamma_EB=p(2); beta_BS=p(3); mu_EB=p(4);  alpha_PS=p(5); 
teta_BP=p(6); gamma_EP=p(7); beta_PS=p(8); mu_EP=p(9);  alpha_ES=p(10); 
beta_ES=p(11);

SL=ST-BS-PS-ES;     %Espacio disponible

dBS=((alpha_BS*BS*SL)/(1+(ES*gamma_EB)))-(BS*(beta_BS+(mu_EB*ES))) ; %Tepozanes
dPS=(((alpha_PS*PS*SL)*(1+teta_BP*BS))/(1+(ES*gamma_EP)))-(PS*(beta_PS+(mu_EP*ES))) ; %Palo locos
dES=(alpha_ES*ES*SL)-(beta_ES*ES); %Eucaliptos

dy=[dBS; dPS; dES];
end