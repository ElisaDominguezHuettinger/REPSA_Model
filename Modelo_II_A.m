function dydt = Modelo_II_A(t,y,alpha_BS, gamma_EB, beta_BS, mu_EB,  alpha_PS, gamma_BP, gamma_EP, beta_PS, mu_EP,  alpha_ES, beta_ES, ST)


dydt = zeros(3,1);
BS=y(1);
PS=y(2);
ES=y(3);

SL=ST-BS-PS-ES


dydt(1)=((alpha_BS*BS*SL)/(1+(ES*gamma_EB)))-(BS*(beta_BS+(mu_EB*ES))) ; %Tepozanes
dydt(2)=((alpha_PS*PS*SL)/((1+(BS*gamma_BP))+(ES*gamma_EP)))   -(PS*(beta_PS+(mu_EP*ES))) ; %Palo locos
dydt(3)=(alpha_ES*ES*SL)-(beta_ES*ES); %Eucaliptos
