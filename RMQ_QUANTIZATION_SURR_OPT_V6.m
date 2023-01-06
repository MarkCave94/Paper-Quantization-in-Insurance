clear;
clc;
toll=10e-10;
r=0.05;  
T=4; 
S0=100; 
sigma=0.15; 
q=0;
i_min=0.03;
i_tec=0.03;
beta=linspace(0.4,1,13);
C_0=100;
NS=120;
K=12;
dt=T/(K*T) ;


for g=1 :length(beta) 


[Gamma,Rmin,Rmax,Prob]= quantizzatore(NS,toll);

Gamma_X =S0+S0*r*dt+S0*sigma*Gamma*sqrt(dt); %Gamma_x uno %

Passi_successivi=T*K;

Quantizzatori_X=zeros(NS,Passi_successivi);

Probabilita_successive_X=zeros(Passi_successivi,NS);

Quantizzatori_X(:,1)=Gamma_X;

Probabilita_successive_X(1,:)=Prob;

for ripetizioni=2:Passi_successivi


  Gamma=Quantizzatori_X(:,ripetizioni-1);
  Gamma_confronto=zeros(NS,1);
  Prob=Probabilita_successive_X(ripetizioni-1,:);
  M_small= sigma*Gamma*sqrt(dt) ;   
  C_small= Gamma*(1+r*dt) ;  


while(norm(Gamma-Gamma_confronto)>toll)
  
  Gamma_confronto=Gamma;
  Rmin=zeros(NS,1);
  Rmin(1,1)=-inf;
  Rmin(2:end,1)=(Gamma(1:end-1,1)+Gamma(2:end,1))*0.5;
  Rmax=zeros(NS,1);
  Rmax(1:end-1,1)=(Gamma(1:end-1,1)+Gamma(2:end,1))*0.5;
  Rmax(NS,1)=+inf;
  R_norm_maximus_X= (Rmax(:,1)'-C_small(1:end,1))./M_small(1:end,1);
  R_norm_minus_X=    (Rmin(:,1)'-C_small(1:end,1))./M_small(1:end,1);
  DeltaGamma=Gamma(2:end,1)-Gamma(1:end-1,1) ;   
  Transition_matrix=sign(M_small).*(normcdf(R_norm_maximus_X)-normcdf(R_norm_minus_X));
  M=(normpdf(R_norm_minus_X)-normpdf(R_norm_maximus_X));
  f_small=(normpdf(R_norm_maximus_X(:,1:NS-1)));
  J=ones(1,NS);
  J_small=ones(1,NS-1);
  Gradient=2*Prob*((((Gamma*J)')-C_small*J).*Transition_matrix-M_small*J.*M); % è la versione trasposta! 
  Super_diagonal=-0.5*Prob*(((((M_small).^-1) *J_small)).*f_small.*((DeltaGamma*J)'));
  Sub_diagonal=-0.5*Prob*(((((M_small).^-1) *J_small)).*f_small.*((DeltaGamma*J)'));
  Diagonal=2*Prob*Transition_matrix+ [Super_diagonal 0]+[0 Super_diagonal];
  Hessian=diag(Sub_diagonal,-1)+diag(Diagonal,0)+diag(Super_diagonal,+1);
  Gamma=(Gamma-((inv(Hessian))*(Gradient')));

end

Quantizzatori_X(:,ripetizioni)=sort(Gamma);
Transition_matrix= sign(M_small).*(normcdf(R_norm_maximus_X)-normcdf(R_norm_minus_X));
Transition_matrix_tot{ripetizioni}=Transition_matrix;
Prob_new=Prob * Transition_matrix; % nuove
Probabilita_successive_X(ripetizioni,:)=Prob_new;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plot(1:48,Quantizzatori_X(:,2:end),Probabilita_successive_X');

surf(Quantizzatori_X(:,1:end)',(1:48),Probabilita_successive_X);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPTION PRICING Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strike=100;
% 
% euro_call_BS_quant = exp(-r*T)*Probabilita_successive_X(end,:)*max(Quantizzatori_X(:,end)-Strike,0);
% 
% [Call] = blsprice(S0, Strike, r, T, sigma);
% 
% delta_Vanilla_Option = (Call/euro_call_BS_quant-1) * 100 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SURRENDER OPTION PRICING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


Quantizzatori_X=[repmat(S0,NS,1) Quantizzatori_X] ;

year_end=(K:K:K*T);

Fund_Returns=Quantizzatori_X(:,year_end+1)./Quantizzatori_X(:,year_end+1-K) -1  ;

s_min=(i_min-i_tec)/(1+i_tec);

min_gar_rate=max(  ((beta(g).*Fund_Returns)-i_tec)/(1+i_tec),  s_min);


C_t=zeros(NS,T+1);

C_t(:,1)=C_0;

for i = 1 : NS

    for j = 1:T

     C_t(i,j+1)=C_t(i,j)*(1+min_gar_rate(i,j));

    end   

end

surf(C_t(:,end)',(1:48),Probabilita_successive_X);


Euro_Contract(g,1)= exp(-r*T)*Probabilita_successive_X(1,:)*C_t(:,end); 

H_small_matrix_call=zeros(NS,T+1);

H_small_matrix_call(:,:)=C_t;



for k=(T):-1:1
    for j=year_end(k):-1:year_end(k)
     Transition_matrix=zeros(NS,NS);
     Transition_matrix =Transition_matrix_tot{j}*Transition_matrix_tot{j-1};
    end    
    
    H_small_matrix_call(:,k)=max(H_small_matrix_call(:,k),exp(-r)*Transition_matrix*H_small_matrix_call(:,k+1));
   
end

AmerPrice(g,1)=exp(-r)*Probabilita_successive_X(1,:)*H_small_matrix_call(:,1);

surrender(g,1)=max(AmerPrice(g)-Euro_Contract(g),0);

end


Results=table(Euro_Contract,AmerPrice,surrender);
writetable(Results,'Results_Quantizzazione_Varying_Beta_4Y.xls')



%%
% AC_Input_data=importdata("AC_Input_data.xlsx");
% 
% AC_AMER=AC_Input_data.data(:,2);
% AC_EU=AC_Input_data.data(:,4);
% AC_SURR=AC_Input_data.data(:,6);
% 
% 
% Y=[AmerPrice,AC_AMER];
% 
% figure
% plot(beta,Y)
% 
% title('American Contract varying Beta')
% xlabel('Beta')
% ylabel('American Contract')
% legend('RMQ','A&C')
% 
% 
% Y=[Euro_Contract,AC_EU];
% 
% figure
% plot(beta,Y)
% 
% title('European Contract varying Beta')
% xlabel('Beta')
% ylabel('European Contract')
% legend('RMQ','A&C')
% 
% Y=[surrender,AC_SURR];
% 
% figure
% plot(beta,Y)
% 
% title('Surrender Option varying Beta')
% xlabel('Beta')
% ylabel('Surrender Option')
% legend('RMQ','A&C')
