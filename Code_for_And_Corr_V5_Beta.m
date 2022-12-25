clear
clc
r=0.05;  
T=4; 
M=1;
dt=1/M;
nsim=400000; 
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
AmerPrice = exp(-r*dt)*mean(CF(:,2));
surrender(g,1)=AmerPrice-Euro_Contract(g,1);



%%
%% support vector regression 

% 
% CF2=C_t(:,1:end);
% 
% for t = (T*M+1)-1:-1:2
%  
%  % Cash flows at time t+1, discounted one period
%  Y = exp(-r*dt)*CF2(:,t+1);
%  
%  X=CF2(:,t);
%  
% Stock=Stock_Index(:,t);
%  
% 
% V=[X];
% 
% Mdl = fitrsvm(V,Y);
% 
% PredCF =predict(Mdl,Y);
%   
% %  paths where continuation is optimal
% J = find(CF2(:,t) < PredCF); 
%  
% % Continued CF are discounted back one period.
% CF2(J,t) = exp(-r*dt)*CF2(J,t+1);
% end
% 
% % The American option price = cash flow at period 2 discounted one peroid
% AmerPrice_SVR = exp(-r*dt)*mean(CF2(:,2));
% surrender_SVR(g,1)=AmerPrice_SVR-Euro_Contract;

end

