close all 
clear all
clc
warning ('off','all');

%% Data Prosessing and Setting parameters
dataT=readmatrix('Series para las tres especies.xlsx','Sheet','Sheet1','Range','A12:A17');
dataX=readmatrix('Series para las tres especies.xlsx','Sheet','Sheet1','Range','B12:E17');
dataX=[dataX(:,1)+dataX(:,2), dataX(:,3:end)];
    
if isfile('OptimModel4.mat')
    Optim=load('OptimModel4.mat');
    p=Optim.p;
    ST=Optim.ST;
    Cost=Optim.Cost;
else
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

    p0=[alpha_BS;gamma_EB;beta_BS;mu_EB;alpha_PS;teta_BP;gamma_EP;beta_PS;mu_EP;alpha_ES;beta_ES];

    lb=zeros(size(p0));
    ub=[];
    min=@(p)CostFunction(p,ST,dataT,dataX);

    tic
    tStop=3600*2;
    
    gs=GlobalSearch('StartPointsToRun','bounds-ineqs','Display','iter','MaxTime',tStop);
    problem=createOptimProblem('fmincon','x0',p0,'objective',min,...
        'lb',lb,'ub',ub);
    [p,Cost]=run(gs,problem);
    toc
    

    save('OptimModel4.mat','p','ST','Cost');

end
file_dir='Simulations';
file_name='OptimSimulation';
SimulateOptim(p,ST,Cost,dataT,dataX,file_dir,file_name)




function SimulateOptim(p,ST,Cost,dataT,dataX,file_dir,file_name)
tspan2=[dataT(2),dataT(end)+20];
tspan1=[dataT(2),dataT(1)];
x0=dataX(2,:);

[t2,y2]=ode15s(@(t,y)model4(t,y,p,ST),tspan2,x0);
[t1,y1]=ode15s(@(t,y)model4(t,y,p,ST),tspan1,x0);
t=[t1(end:-1:2)', t2'];
y=[y1(end:-1:2,:)', y2'];

figure('Name','Model 4: Optim Simulation')
for i=1:3
    ind=~isnan(dataX(:,i));
    subplot(1,3,i)
        plot(t, y(i,:), 'k','LineWidth',2)
        hold on
        scatter(dataT(ind),dataX(ind,i), 'rs', 'filled');
        xlim([t(1) t(end)]);
        xlabel('year')
        ylabel(['x_',num2str(i),' (coverage)'])
        hold on
        axis square
end
par={['alpha_B_S=',num2str(p(1)), '; gamma_E_B=',num2str(p(2)), '; beta_B_S=',num2str(p(3)), '; mu_E_B=',num2str(p(4)), '; alpha_P_S=',num2str(p(5)), '; theta_B_P=',num2str(p(6)),';'], ['gamma_E_P=',num2str(p(7)), '; beta_P_S=',num2str(p(8)), '; mu_E_P=',num2str(p(9)), '; alpha_E_S=',num2str(p(10)), '; beta_E_S=',num2str(p(11)), '; ST=',num2str(ST),'; Cost=',num2str(Cost)]};
h=suptitle({['Simulation Model 4: Optim Parameters '], [par{1}],[par{2}]});
set(h,'FontSize',8,'FontWeight','normal')

f=fullfile(file_dir,file_name);
print(gcf,f,'-dpdf','-r0')
end



function Cost=CostFunction(p,ST,dataT,dataX)
dataX=dataX';
x0=dataX(:,2);
tspan2=[dataT(2),dataT(end)];
tspan1=[dataT(2),dataT(1)];

Sol2=ode15s(@(t,y)model4(t,y,p,ST),tspan2,x0);
Sol1=ode15s(@(t,y)model4(t,y,p,ST),tspan1,x0);
Y2=deval(Sol2,dataT(2:end));
Y1=Sol1.y(:,end);
Y=[Y1,Y2];
ind=~isnan(dataX);
Cost=sum((Y(ind)-dataX(ind)).^2,'all');
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