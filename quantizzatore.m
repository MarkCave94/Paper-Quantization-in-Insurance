function [Gamma,Rmin,Rmax,Prob]= quantizzatore(N,toll)
n=1:1:N;
Gamma= (5.5*n/(N+1)-2.75)';
Grad=1;
while any(Grad>toll)
Rmax = (Gamma(1:end-1,1)+Gamma(2:end,1))/2; Rmax=vertcat(Rmax,inf);
Rmin= (Gamma(1:end-1,1)+Gamma(2:end,1))/2;  Rmin=vertcat(-inf,Rmin);

f_x=normpdf(Rmax(1:end-1,1));

Delta_Gamma=Gamma(2:end,1)-Gamma(1:end-1,1);

M_x=normpdf(Rmin)-normpdf(Rmax);

Prob=(normcdf(Rmax)-normcdf(Rmin))';


Grad=2*Gamma.*(Prob')-2*M_x;
hoff=-0.5*((f_x.*Delta_Gamma));
super=[hoff; 0];
lower=[0 ;hoff];
hmain=diag(2*Prob')+diag(super)+diag(lower);
hessiano=hmain + diag(hoff,-1)+diag(hoff,1);
Gamma=Gamma-(inv(hessiano))*Grad;
end