% American Contract


beta_vec=0.4:0.05:1;

for g=1:size(beta_vec,2)

i_min=0.03;
i_tec=0.03;
beta=beta_vec(g);
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

reval_rate(j+1)=max ((   (beta*gamma(j+1))-i_tec    )/(1+i_tec),s_min);

end

n_proof=round(N/2+1-log(1+i_tec/beta)/(2*log(u))); % Dovrebbe essere uguale a riga 65
s_min=(i_min-i_tec)/(1+i_tec);
n=sum(reval_rate>s_min);
a=reval_rate(1:n);
b=Q(1:n);
c=a*b';
prezzo_binomiale_eur= C0*exp(-r*T)*(1+c)^(T);
surrender_start=zeros(1,T);
surrender_start(T+1)=C0*(1+c)^(T);

for t=T:-1:1

c_t=C0*exp(-r*t)*(1+c)^(t);

surrender_start(t)=max(exp(-r)*surrender_start(t+1),c_t) ; 


end

prezzo_americano(g,1)=surrender_start(1);
surrender(g,1)=prezzo_americano(g,1)-prezzo_binomiale_eur;

surrender_start=zeros(1,T);

end