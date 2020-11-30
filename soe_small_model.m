function [fx,fxp,fy,fyp,f] = soe_small_model
global approx nx betta R yy
%Define parameters
syms betta R yy rho
%Define variables 
syms c cf lambda lambdaf bl b bf w wf %% Add new parameter

%define utility function
utility=log(c);%%define the utility function here;
%define budget constraint of consumer
bc  = -c - b + w + R*bl;    %% define the budget constraint at time t here;
bcf = -cf - bf + wf + R*b;%% define the budget cosntraint at time t+1 here;

%define lagrangian
lagrangian_cons = utility +lambda *bc + lambdaf*betta*bcf;    %%s specify the lagrangian problem;
f1  = jacobian(lagrangian_cons,c);%% foc c
f2  = lambda - lambdaf*R*betta;%% foc b_{t}
f3  = -c - b + w + R*bl;%% foc lambda
f4  = (1-rho)*log(yy) + rho*log(w) - log(wf); %% foc endowment process

%create f
f = [f1;f2;f3;f4];
% Define the vector of controls, y, and states, x
x=[bl, w];  %% state variable at time t
y=[c lambda];
xp=[b, wf];  %% state variable at time t+1
yp=[cf lambdaf];
lv=[x,y,xp,yp];
for i=1:size(lv,2)
f = subs(f,lv(i),exp(lv(i)));
end
nx=size(x,2);
%Compute analytical derivatives of f
[fx,fxp,fy,fyp]=anal_deriv(f,x,y,xp,yp,approx); 

