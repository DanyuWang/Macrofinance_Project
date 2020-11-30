%MODEL_RUN.M
%Luisa Lambertini
%Date October 2018
clear all
close all
%addpath (Volumes/sfi-ll/SFI-LL/MFE/Fall 2019/Problem Sets/ps1/soe_small)
global s2 t2 approx nx %additional global variables
%degree of approximation in the linearization around the steady state
approx=1;
[fx,fxp,fy,fyp,f] = soe_small_model;

%Steady State and Parameter Values
[betta,c,cf,lambda,lambdaf,yy,w,wf,bl,b,bf,R,rho] = soe_small_ss;  %% update ss function

%Obtain numerical derivatives of f
num_eval
%evaluate f at the steady state
eval(f)

%First-order approximation
[gx,hx] = gx_hx_new(nfy,nfx,nfyp,nfxp);
%order eigenvalues
tt=sort(abs(diag(s2).\diag(t2)));
disp('stable eigenvalues')
tt(1:nx)
%shocks
x0 =[0, 1];  %% define the endowment shock

IR=ir(gx,hx,x0,100);

titlegraph=['c        ';...
            'lambda   ';...
            'b        ';...
            'y        '];

nrow = 2;  %number of rows in the plot;
ncol = 2;   %number of column in the plot;
hFig = figure(1);
for j=1:size(IR,2)
    subplot(nrow,ncol,j);
    plot(IR(:,j),'-b');
    title(titlegraph(j,:));
      xlabel('Quarters after 0');
      ylabel('% dev from stst');
end
set(gcf,'color','w')
set(hFig,'Units','normalized','Position',[0 0 .6 .6]);

% save figure
% print('IRF','-dpng')
