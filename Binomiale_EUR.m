% %% Binomial Tree approach 

clear
clc
beta_vec=0.4:0.05:1;

for i=1:size(beta_vec,2)
i_min=0.03;
i_tec=0.03;
beta=beta_vec(i);
r=0.05; 
C0=100;
T=4;
N=50;
delta=1/N;
sigma=0.15; 
u=exp(sigma*sqrt(delta));
d=1/u;
q = (((1+r)^delta)-d)/(u-d);
gamma=zeros(1,N+1);
reval_rate=zeros(1,N+1);
Q=zeros(1,N+1);

for j=0:N

gamma(j+1)= ((u^(N-j))*d^(j))-1;

Q(j+1)=binopdf(N-j,N,q);

reval_rate(j+1)=(   (beta*gamma(j+1))-i_tec    )/(1+i_tec);

end

n_proof=round(N/2+1-log(1+i_tec/beta)/(2*log(u))); % Dovrebbe essere uguale a riga 65
s_min=(i_min-i_tec)/(1+i_tec);
n=sum(reval_rate>s_min);
a=reval_rate(1:n);
b=Q(1:n);
c=a*b';

% lambda=(r-c)/(1+c)
% 
% prova=(C0/(1+c))*(1+lambda)^(-T)
% 
% prezzo=((C0)/(1+c))*((1+r)/(1+c))^(-T)

prezzo_binomiale(i)= C0*exp(-r*T)*(1+c)^(T)

end


