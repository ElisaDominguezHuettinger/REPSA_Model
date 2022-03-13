close all
clear all
clc

ST=146.9;
% Analyse the possible dynamic behaviours of the two 3D models:
 % matrices created in
 
 
%TresSP_c_int_p_nat_n_exo_Proportion_of_behaviours.m % MODEL II-A
% m3SPintNEG_Proportion_of_behaviours.m  % MODEL II-B

load('Behaviour_matrix_ModelIIA.mat')
 
behaviour_matrix_IIA=behaviour_matrix;

%%

load('Behaviour_matrix_ModelIIB.mat')
 
behaviour_matrix_IIB=behaviour_matrix;



%%
behaviour_matrix=behaviour_matrix_IIB;

jj=length(behaviour_matrix); %number of iterations

sum(behaviour_matrix(:,1)==1)/jj
%0.7253 are monostable
% let's have a look at those that are monostable
index_mono=find(behaviour_matrix(:,1)==1);


% should retrieve 1
multistab_IIB=[sum(behaviour_matrix(:,1)==2)/jj ; sum(behaviour_matrix(:,1)==1)/jj; sum(behaviour_matrix(:,1)==3)/jj];

sum(multistab_IIB)

%
number_mono=sum(behaviour_matrix(:,1)==1); 
 

% coexistence
Eu_only= sum(behaviour_matrix(index_mono,2)==0 & behaviour_matrix(index_mono,3)==0 & behaviour_matrix(index_mono,4)>0)/number_mono;

BS_only= sum(behaviour_matrix(index_mono,2)>0 & behaviour_matrix(index_mono,3)==0 & behaviour_matrix(index_mono,4)==0)/number_mono;

PS_only= sum(behaviour_matrix(index_mono,2)==0 & behaviour_matrix(index_mono,3)>0  & behaviour_matrix(index_mono,4)==0 )/number_mono;

coex= sum(behaviour_matrix(index_mono,2)>0 & behaviour_matrix(index_mono,3)>0)/number_mono;

all_extint= sum(behaviour_matrix(index_mono,2)==0 & behaviour_matrix(index_mono,3)==0 & behaviour_matrix(index_mono,4)==0)/number_mono;


Eu_only+BS_only+PS_only+coex+all_extint
mono_IIB=[all_extint; Eu_only; BS_only; PS_only; coex]

sum(mono_IIB)

%%
behaviour_matrix=behaviour_matrix_IIA;

jj=length(behaviour_matrix); %number of iterations

sum(behaviour_matrix(:,1)==1)/jj
%0.7253 are monostable
% let's have a look at those that are monostable
index_mono=find(behaviour_matrix(:,1)==1);


%% should retrieve 1 - Model IIA has up to 3 stable steady states!!!!!!
multistab_IIA=[sum(behaviour_matrix(:,1)==2)/jj ; sum(behaviour_matrix(:,1)==1)/jj  ; sum(behaviour_matrix(:,1)==3)/jj];


%%
number_mono=sum(behaviour_matrix(:,1)==1); 
 

%% coexistence
Eu_only= sum(behaviour_matrix(index_mono,2)==0 & behaviour_matrix(index_mono,3)==0 & behaviour_matrix(index_mono,4)>0)/number_mono;

BS_only= sum(behaviour_matrix(index_mono,2)>0 & behaviour_matrix(index_mono,3)==0 & behaviour_matrix(index_mono,4)==0)/number_mono;

PS_only= sum(behaviour_matrix(index_mono,2)==0 & behaviour_matrix(index_mono,3)>0  & behaviour_matrix(index_mono,4)==0 )/number_mono;

coex= sum(behaviour_matrix(index_mono,2)>0 & behaviour_matrix(index_mono,3)>0)/number_mono;

all_extint= sum(behaviour_matrix(index_mono,2)==0 & behaviour_matrix(index_mono,3)==0 & behaviour_matrix(index_mono,4)==0)/number_mono;


Eu_only+BS_only+PS_only+coex+all_extint
%

mono_IIA=[all_extint; Eu_only; BS_only; PS_only; coex];

 sum(mono_IIA)

%%
subplot(3,1,1)
x = [1,2,3];
y = [multistab_IIA,  multistab_IIB];
b=bar(x,y);
b(1).FaceColor = 'r';
b(2).FaceColor = 'b';
legend('II A (-)', 'II B (+)')
ylabel('Proportion')
xlabel('number of stable steady states')
axis square
 
 %%

subplot(3,1,2)
x = categorical({' all extint' ,'ES only', 'BS only', 'PS only', 'coex'});
y = [mono_IIA,  mono_IIB];
b=bar(x,y)
b(1).FaceColor = 'r';
b(2).FaceColor = 'b';

legend('II A (-)', 'II B (+)')
title('Mono-stable')
ylabel('proportion of behaviours')
xlabel('long term community configurations')
axis square



%% Now we will look at the bistable
for ii=[0, 3]
behaviour_matrix=behaviour_matrix_IIA;
% let's have a look at those that are bistable
index_bi=find(behaviour_matrix(:,1)==2);
number_bi=sum(behaviour_matrix(:,1)==2); 

Eu_only= sum(behaviour_matrix(index_bi,2+ii)==0 & behaviour_matrix(index_bi,3+ii)==0 & behaviour_matrix(index_bi,4+ii)>0)/number_bi;
BS_only= sum(behaviour_matrix(index_bi,2+ii)>0 & behaviour_matrix(index_bi,3+ii)==0 & behaviour_matrix(index_bi,4+ii)==0)/number_bi;
PS_only= sum(behaviour_matrix(index_bi,2+ii)==0 & behaviour_matrix(index_bi,3+ii)>0  & behaviour_matrix(index_bi,4+ii)==0 )/number_bi;
coex= sum(behaviour_matrix(index_bi,2+ii)>0 & behaviour_matrix(index_bi,3+ii)>0)/number_bi;
all_extint= sum(behaviour_matrix(index_bi,2+ii)==0 & behaviour_matrix(index_bi,3+ii)==0 & behaviour_matrix(index_bi,4+ii)==0)/number_bi;
Eu_only+BS_only+PS_only+coex+all_extint
bi_IIA=[all_extint; Eu_only; BS_only; PS_only; coex]
sum(bi_IIA)


%% 
behaviour_matrix=behaviour_matrix_IIB;
% let's have a look at those that are bistable
index_bi=find(behaviour_matrix(:,1)==2);
number_bi=sum(behaviour_matrix(:,1)==2); 

% only EU only or BS only!!!

%% First steady state
Eu_only= sum(behaviour_matrix(index_bi,2+ii)==0 & behaviour_matrix(index_bi,3+ii)==0 & behaviour_matrix(index_bi,4+ii)>0)/number_bi;
BS_only= sum(behaviour_matrix(index_bi,2+ii)>0 & behaviour_matrix(index_bi,3+ii)==0 & behaviour_matrix(index_bi,4+ii)==0)/number_bi;
PS_only= sum(behaviour_matrix(index_bi,2+ii)==0 & behaviour_matrix(index_bi,3+ii)>0  & behaviour_matrix(index_bi,4+ii)==0 )/number_bi;
coex= sum(behaviour_matrix(index_bi,2+ii)>0 & behaviour_matrix(index_bi,3+ii)>0)/number_bi;
all_extint= sum(behaviour_matrix(index_bi,2+ii)==0 & behaviour_matrix(index_bi,3+ii)==0 & behaviour_matrix(index_bi,4+ii)==0)/number_bi;
Eu_only+BS_only+PS_only+coex+all_extint
bi_IIB=[all_extint; Eu_only; BS_only; PS_only; coex]
sum(bi_IIB)



%%
%figure;
if ii==0
    uu=1;
else
    uu=2
end;
subplot(3,1,uu)
x = categorical({' all extint' ,'ES only', 'BS only', 'PS only', 'coex'});
y = [bi_IIA, bi_IIB];
b=bar(x, y)
b(1).FaceColor = 'r';
b(2).FaceColor = 'b';
title('Bi-stable')
ylabel('proportion of behaviours')
xlabel('long term community configurations')
if ii==0
    title('first steady state')
else
   title('second steady state')
end;
axis square
end



%% tristable
behaviour_matrix=behaviour_matrix_IIA;
% let's have a look at those that are bistable
index_tri=find(behaviour_matrix(:,1)==3)
number_tri=sum(behaviour_matrix(:,1)==3)

%%
figure
subplot(3,3,1)
histogram(behaviour_matrix(index_tri,2)./ST, 100, 'Normalization','probability', 'FaceColor', 'yellow')
hold on
axis square
ylabel('BS')
xlim([0, 1])

subplot(3,3,2)
histogram(behaviour_matrix(index_tri,2+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'cyan')
axis square
xlim([0, 1])

subplot(3,3,3)
histogram(behaviour_matrix(index_tri,2+6)./ST, 100, 'Normalization','probability', 'FaceColor', 'green')
xlim([0, 1])
axis square
%
subplot(3,3,1+3)
histogram(behaviour_matrix(index_tri,3)./ST, 100, 'Normalization','probability', 'FaceColor', 'yellow')
hold on
axis square
ylabel('PS')
xlim([0, 1])

subplot(3,3,2+3)
histogram(behaviour_matrix(index_tri,3+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'cyan')
axis square
xlim([0, 1])

subplot(3,3,3+3)
histogram(behaviour_matrix(index_tri,3+6)./ST, 100, 'Normalization','probability', 'FaceColor', 'green')
xlim([0, 1])
axis square


%
subplot(3,3,1+6)
histogram(behaviour_matrix(index_tri,4)./ST, 100, 'Normalization','probability', 'FaceColor', 'yellow')
hold on
axis square
ylabel('ES')
xlim([0, 1])

subplot(3,3,2+6)
histogram(behaviour_matrix(index_tri,4+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'cyan')
axis square
xlim([0, 1])

subplot(3,3,3+6)
histogram(behaviour_matrix(index_tri,4+6)./ST, 100, 'Normalization','probability', 'FaceColor', 'green')
xlim([0, 1])
axis square





%%




Eu_only= sum(behaviour_matrix(index_tri,2)==0 & behaviour_matrix(index_tri,3)==0 & behaviour_matrix(index_tri,4)>0)/number_tri

BS_only= sum(behaviour_matrix(index_bi,2)>0 & behaviour_matrix(index_bi,3)==0 & behaviour_matrix(index_bi,4)==0)/number_bi;
PS_only= sum(behaviour_matrix(index_bi,2)==0 & behaviour_matrix(index_bi,3)>0  & behaviour_matrix(index_bi,4)==0 )/number_bi;
coex= sum(behaviour_matrix(index_bi,2)>0 & behaviour_matrix(index_bi,3)>0)/number_bi;
all_extint= sum(behaviour_matrix(index_bi,2)==0 & behaviour_matrix(index_bi,3)==0 & behaviour_matrix(index_bi,4)==0)/number_bi;
Eu_only+BS_only+PS_only+coex+all_extint
bi_IIA=[all_extint; Eu_only; BS_only; PS_only; coex]
sum(bi_IIA)

















%% Now we look at the distrubutions of second steady states for the two cases
behaviour_matrix=behaviour_matrix_IIA
Eu_only_ind= find(behaviour_matrix(:,1)==2 & behaviour_matrix(:,2)==0 & behaviour_matrix(:,3)==0 & behaviour_matrix(:,4)>0);
BS_only_ind= find(behaviour_matrix(:,1)==2 & behaviour_matrix(:,2)>0 & behaviour_matrix(:,3)==0 & behaviour_matrix(:,4)==0);



ST=146.9;
figure;
subplot(3,2,1)
histogram(behaviour_matrix(Eu_only_ind,2+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red')
ylabel('BS')
xlim([0, 1]);
hold on
axis square
subplot(3,2,3)
histogram(behaviour_matrix(Eu_only_ind,3+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red')
ylabel('PS')
xlim([0, 1]);
axis square
subplot(3,2,5)
histogram(behaviour_matrix(Eu_only_ind,4+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red')
xlim([0, 1]);
ylabel('ES')
axis square

suptitle('Second steady state, bistable, Model II-A, First steady state is ES only')

%%

figure;
subplot(3,1,1)
histogram(behaviour_matrix(BS_only_ind,2+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red')
xlim([0, 1]);
hold on
subplot(3,1,2)
histogram(behaviour_matrix(BS_only_ind,3+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'yellow')
xlim([0, 1]);
subplot(3,1,3)
histogram(behaviour_matrix(BS_only_ind,4+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'blue')
xlim([0, 1]);

%% And the tristable
index_tri=find(behaviour_matrix(:,1)==3)
number_tri=sum(behaviour_matrix(:,1)==3) 

%% First steady state
Eu_only= sum(behaviour_matrix(index_tri,2)==0 & behaviour_matrix(index_tri,3)==0 & behaviour_matrix(index_tri,4)>0)/number_tri;
BS_only= sum(behaviour_matrix(index_tri,2)>0 & behaviour_matrix(index_tri,3)==0 & behaviour_matrix(index_tri,4)==0)/number_tri;

PS_only= sum(behaviour_matrix(index_tri,2)==0 & behaviour_matrix(index_tri,3)>0  & behaviour_matrix(index_tri,4)==0 )/number_tri;
coex= sum(behaviour_matrix(index_tri,2)>0 & behaviour_matrix(index_tri,3)>0)/number_tri;
all_extint= sum(behaviour_matrix(index_tri,2)==0 & behaviour_matrix(index_tri,3)==0 & behaviour_matrix(index_tri,4)==0)/number_tri;
Eu_only+BS_only+PS_only+coex+all_extint
tri_IIA=[all_extint; Eu_only; BS_only; PS_only; coex];
sum(tri_IIA)

% only  BS only!!!

%%
BS_only_ind= find(behaviour_matrix(:,1)==3 & behaviour_matrix(:,2)>0 & behaviour_matrix(:,3)==0 & behaviour_matrix(:,4)==0);



figure;
subplot(3,1,1)
histogram([behaviour_matrix(BS_only_ind,2+3) behaviour_matrix(BS_only_ind,2+6)]./ST, 100, 'Normalization','probability', 'FaceColor', 'red')
xlim([0, 1]);
hold on
subplot(3,1,2)
histogram([behaviour_matrix(BS_only_ind,3+3) behaviour_matrix(BS_only_ind,3+6)]./ST, 100, 'Normalization','probability', 'FaceColor', 'yellow')
xlim([0, 1]);
subplot(3,1,3)
histogram([behaviour_matrix(BS_only_ind,4+3) behaviour_matrix(BS_only_ind,4+6)]./ST, 100, 'Normalization','probability', 'FaceColor', 'blue')
xlim([0, 1]);

%%

% let's have a look at those that are bistable
index_bi_IIB=find(behaviour_matrix_IIB(:,1)==2);
% let's have a look at those that are bistable
index_bi_IIA=find(behaviour_matrix_IIA(:,1)==2);


figure
subplot(3,2,1)
histogram(behaviour_matrix_IIA(index_bi_IIA,2)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
ylabel('BS')
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])

subplot(3,2,2)
histogram(behaviour_matrix_IIA(index_bi_IIA,2+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])

subplot(3,2,3)
histogram(behaviour_matrix_IIA(index_bi_IIA,3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
ylabel('PS')
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])

subplot(3,2,4)
histogram(behaviour_matrix_IIA(index_bi_IIA,3+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])

subplot(3,2,5)
histogram(behaviour_matrix_IIA(index_bi_IIA,4)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
ylabel('ES')
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])
subplot(3,2,6)
histogram(behaviour_matrix_IIA(index_bi_IIA,4+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])

suptitle('Model II-A, bistable')


%%


figure
subplot(3,2,1)
histogram(behaviour_matrix_IIB(index_bi_IIB,2)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
ylabel('BS')
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])

subplot(3,2,2)
histogram(behaviour_matrix_IIB(index_bi_IIB,2+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])

subplot(3,2,3)
histogram(behaviour_matrix_IIB(index_bi_IIB,3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
ylabel('PS')
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])

subplot(3,2,4)
histogram(behaviour_matrix_IIB(index_bi_IIB,3+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])

subplot(3,2,5)
histogram(behaviour_matrix_IIB(index_bi_IIB,4)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
ylabel('ES')
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])
subplot(3,2,6)
histogram(behaviour_matrix_IIB(index_bi_IIB,4+3)./ST, 100, 'Normalization','probability', 'FaceColor', 'red', 'EdgeColor', 'red')
axis square
xlim([0, 1])
set(gca,'xtick',[])
ylim([0, 1]); 
set(gca,'ytick',[])

suptitle('Model II-B, bistable')

%% Now we will look at the bistable

% PLEASE NOTE: TO GET THE COMPUTATIONS FOR THESE LAST LINES I MANUALLY CHANGED SOME TERMS (NOT TO REPEAT..)THESE LAST LINES 
% in line 494: change between IIA and IIB
% in lines 511 and 512: change between behaviour_matrix(index_bi,4+ii)>0
% and behaviour_matrix(index_bi,4+ii)==0 and
% behaviour_matrix(:,2)==0  and behaviour_matrix(:,2)>0 

ii=0;

% first steady state 

behaviour_matrix=behaviour_matrix_IIB;
% let's have a look at those that are bistable
index_bi=find(behaviour_matrix(:,1)==2);
number_bi=sum(behaviour_matrix(:,1)==2); 

Eu_only= sum(behaviour_matrix(index_bi,2+ii)==0 & behaviour_matrix(index_bi,3+ii)==0 & behaviour_matrix(index_bi,4+ii)>0)/number_bi;
BS_only= sum(behaviour_matrix(index_bi,2+ii)>0 & behaviour_matrix(index_bi,3+ii)==0 & behaviour_matrix(index_bi,4+ii)==0)/number_bi;
PS_only= sum(behaviour_matrix(index_bi,2+ii)==0 & behaviour_matrix(index_bi,3+ii)>0  & behaviour_matrix(index_bi,4+ii)==0 )/number_bi;
coex= sum(behaviour_matrix(index_bi,2+ii)>0 & behaviour_matrix(index_bi,3+ii)>0)/number_bi;
all_extint= sum(behaviour_matrix(index_bi,2+ii)==0 & behaviour_matrix(index_bi,3+ii)==0 & behaviour_matrix(index_bi,4+ii)==0)/number_bi;
Eu_only+BS_only+PS_only+coex+all_extint
bi_IIA=[all_extint; Eu_only; BS_only; PS_only; coex]
sum(bi_IIA)
%%


%%
Eu_only_bi_index= find(behaviour_matrix(:,1)==2 & behaviour_matrix(:,2)==0 & behaviour_matrix(:,3)==0 & behaviour_matrix(:,4)>0)
Eu_only_bi_number= sum(behaviour_matrix(:,1)==2 & behaviour_matrix(:,2)==0 & behaviour_matrix(:,3)==0 & behaviour_matrix(:,4)>0)

ii=3;
Eu_only= sum(behaviour_matrix(Eu_only_bi_index,2+ii)==0 & behaviour_matrix(Eu_only_bi_index,3+ii)==0 & behaviour_matrix(Eu_only_bi_index,4+ii)>0)/Eu_only_bi_number;
BS_only= sum(behaviour_matrix(Eu_only_bi_index,2+ii)>0 & behaviour_matrix(Eu_only_bi_index,3+ii)==0 & behaviour_matrix(Eu_only_bi_index,4+ii)==0)/Eu_only_bi_number;
PS_only= sum(behaviour_matrix(Eu_only_bi_index,2+ii)==0 & behaviour_matrix(Eu_only_bi_index,3+ii)>0  & behaviour_matrix(Eu_only_bi_index,4+ii)==0 )/Eu_only_bi_number;
coex= sum(behaviour_matrix(Eu_only_bi_index,2+ii)>0 & behaviour_matrix(Eu_only_bi_index,3+ii)>0)/Eu_only_bi_number;
all_extint= sum(behaviour_matrix(Eu_only_bi_index,2+ii)==0 & behaviour_matrix(Eu_only_bi_index,3+ii)==0 & behaviour_matrix(Eu_only_bi_index,4+ii)==0)/Eu_only_bi_number;
%%
labels = {' all extint' ,'ES only', 'BS only', 'PS only', 'coex'};

pie([all_extint; Eu_only; BS_only; PS_only; coex],labels)
[all_extint; Eu_only; BS_only; PS_only; coex]
