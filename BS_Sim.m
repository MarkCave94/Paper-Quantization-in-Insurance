function [S] = BS_Sim(sigma,S0,r,q,M,T,N)

% Euler for BS
% INPUTS


%   params = Heston parameters
%   S0  = Spot Price
%   Mat = Maturity
%   r = riskf ree rate
%   q = dividend yield
%   T   = Number of time steps
%   N   = Number of stock price paths
% OUTPUTS
%   S = Vector of simulated stock prices
% Time increment
dt = 1/M;

% Initialize the variance and stock processes
S = zeros(T*M,N);

% Starting values for the variance and stock processes
S(1,:) = S0;       % Spot price 

% Generate the stock and volatility paths
for i=1:N;
	for t=2:T*M+1;
		% Generate two dependent N(0,1) variables with correlation rho
		Zv = randn(1);
		% Discretize the log stock price
		S(t,i) = S(t-1,i)*exp((r-1/2*sigma^2)*dt + sigma*sqrt(dt)*Zv);
	end
end

S=S';