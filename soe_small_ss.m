function [betta,c,cf,lambda,lambdaf,yy,w,wf,bl,b,bf,R,rho] = soe_small_ss


rho = 0.90; %% Add AR parameter
betta = 0.99;
R=1/betta;
b=5; %% steady-state bond holdings
yy = 1; %% steady-state endowment 
c=yy+b*(R-1);%% steady-state consumption
lambda=1/c;%% steady-state lagrangian multiplier

%take logs
c=log(c);
lambda=log(lambda);
b=log(b);
w=log(yy); %% Endowment at time t

%define t+1 and t-1 variables at the steady state
cf=c;
lambdaf=lambda;
bl=b;
bf=b;
wf=w; %% Endowment at time t+1
