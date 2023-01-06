clear
clc
r=0.05;  
T=10; 
M=1;
dt=1/M;
nsim=100000; 
S0=100; 
sigma=0.15; 
q=0;
i_min=0.03;
i_tec=0.03;
beta=linspace(0.4,1,13);
C_0=100;
surrender=zeros(length(beta),1);
surrender_SVR=zeros(length(beta),1);

for g=1 :length(beta) 

[S] = BS_Sim(sigma,S0,r,q,1,T,nsim);


I=S(:,2:end)./S(:,1:end-1)-1;

Stock_Index=S(:,2:end);

s_min=(i_min-i_tec)/(1+i_tec);

min_gar_rate=max(((beta(g).*I)-i_tec)/(1+i_tec),s_min);


C_t=zeros(nsim,(T*M)+1);

C_t(:,1)=C_0;

for i = 1 : nsim

    for j = 1:(T*M)

     C_t(i,j+1)=C_t(i,j)*(1+min_gar_rate(i,j));

    end   

end

Euro_Contract(g,1)=exp(-r*T)*mean(C_t(:,end));

mean(min_gar_rate(:,end))


%% 

% Longstaff-Schwartz method for surrender option


% Set the last cash flows to the intrinsic value.

CF=C_t(:,1:end);

for t = (T*M+1)-1:-1:2
 
 % Cash flows at time t+1, discounted one period
 Y = exp(-r*dt)*CF(:,t+1);
 
 X=CF(:,t);
 
Stock=Stock_Index(:,t);
 

V=[Stock,X,(1+ min_gar_rate(:,t))];

mdl2 = fitlm(V,Y,'poly020');
% Regression parameters and predicted cash flows

mdl2.Rsquared;

PredCF = mdl2.Fitted;
  
%  paths where continuation is optimal
J = find(CF(:,t) < PredCF); 
 
% Continued CF are discounted back one period.
CF(J,t) = exp(-r*dt)*CF(J,t+1);
end

% The American option price = cash flow at period 2 discounted one peroid
AmerPrice(g,1) = exp(-r*dt)*mean(CF(:,2));
surrender(g,1)=AmerPrice(g,1)-Euro_Contract(g,1);




end
%%
Results=[ AmerPrice Euro_Contract surrender];
writematrix(Results,'Results_LSMC_Varying_Beta_4Y.xls')