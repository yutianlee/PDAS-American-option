%% Primal-dual active-set method for American options with regime-switching
% default: two-regimes

clear;  
clc;
close all;
format long

%%%%%%%%%%%%%% Model parameters %%%%%%%%
T = 1.0;  % expiry time     

K1 = 10; 
K2 = 10; 
K = K1;  %strike price K

sigma1 = 0.5;  %volatility in regime 1
sigma2 = 0.4;  %volatility in regime 2
r1 = 0.2;  % risk-free interest rate in regime 1
r2 = 0.2;  %risk-free interest rate in regime 2
sigma = [sigma1 sigma2];
r = [r1 r2];

%generation matrix A
a11 = -0.05;
a12 = 0.05;
a21 = 0.15;
a22 = -0.15;
a = [a11 a12;a21 a22];  

%%%%%%%%%%% Domain truncation %%%%%%%%%%
% Truncate the spatial domain to [-L,L]
% find L0 and X as in Lemma 1, see Eq. (10)
epsilon = 1e-6;
L0 = max(-1.25*sigma.^2*T.*(r./sigma.^2-0.5)+0.5*sqrt(6.25*sigma.^4*(T^2).*(r./sigma.^2-0.5).^2-10*sigma.^2*T.*log(epsilon/sqrt(5*K))));
X = min(2*r./(2*r+sigma));
% calculate L using Eq. (11)
L = max([-log(K*X) L0+log(K)]); 

%%%%%%%%%%%%% FDM parameters %%%%%%%%%%%%
x0 = -L;  % left end-points of the spatial domain
t0 = 0;   % initial time
M = 2000;  % number of time steps
N = 1500;  % number of spatial nodes
dt = T/M;   % temporal step size 
dx = 2*L/N;  % spatial step size 

% defining the nodes
t = linspace(0,T,M+1);   % time levels
x = linspace(-L,L,N+1); % spatial nodes
S = exp(x);   % underlying asset price 

%solving the option values and optimal exercise boundries using FDM+PDAS
[V1,V2,b_1,b_2] = FDM_PDAS(M,N,dx,dt,x,t,L,sigma,r,a,x0,t0,K);

v1=V1(M+1,:);  %option values P(s,0) in regime 1
v2=V2(M+1,:);  %option values P(s,0) in regime 2



%%%%%%%%%% Plots %%%%%%%%%%%


hold on
plot_3D(V1,V2,t,S,K,b_1,b_2)  %three-dimensional graphs of P(S,t) values with two regimes 
figure
plot_v0(S,v1,v2,K)  %two-dimensional graphs of initial values P(S,0) values with two regimes 
figure
plot_boundary_2D(b_1,b_2,t)  %optimal exercise boundaries with two regimes



function plot_3D(V1,V2,t,S,K,b_1,b_2)
    P1 = flipud(V1);
    P2 = flipud(V2);
   
    
    subplot(1,2,1)
    mesh(t,S,P1')
    hold on
    plot3(t,b_1,K-b_1,'k-','LineWidth',3)
    set(gca,'YDir','reverse')
    axis([0 1 0 30 0 10])
    title('$Example \ 1: option \ value \ under \ regime \ 1$','Interpreter','latex','fontsize',18);
    xlabel('$t$','Interpreter','latex','fontsize',18)
    ylabel('$S$','Interpreter','latex','fontsize',18)
    zlabel('$P$','Interpreter','latex','fontsize',18,'rotation',1) 

    subplot(1,2,2)
    mesh(t,S,P2')
    hold on
    plot3(t,b_2,K-b_2,'k-','LineWidth',3)
    title('$Example \ 1: option \ value \ under \ regime \ 2$','Interpreter','latex','fontsize',18);
    set(gca,'YDir','reverse')
    axis([0 1 0 30 0 10])
    xlabel('$t$','Interpreter','latex','fontsize',18)
    ylabel('$S$','Interpreter','latex','fontsize',18)
    zlabel('$P$','Interpreter','latex','fontsize',18,'rotation',1) 
end
    
function plot_v0(s,v1,v2,K)
    hold on
    plot(s,max(K-s,0),'c--')
    plot(s,v1,'mo-','MarkerSize',3)
    plot(s,v2,'go-','MarkerSize',3)

    title('$Example \ 1: option \ value \ at \ t=0$','Interpreter','latex','fontsize',18);
    xlabel('$S$','Interpreter','latex','fontsize',18)
    ylabel('$P$','Interpreter','latex','fontsize',18,'rotation',0)
    legend('max(K-S,0)','PDAS-P1','PDAS-P2')
    axis([0 30 0 9]) 
end

function plot_boundary_2D(b_1,b_2,t)
    hold on
    plot(b_1,t,'go-');   
    plot(b_2,t,'mo-');
    
    title('$Example \ 1: optimal \ exercise \ boundary$','Interpreter','latex','fontsize',18);
    ylabel('$t$','Interpreter','latex','fontsize',18,'rotation',0)
    xlabel('$S$','Interpreter','latex','fontsize',18)
    axis([b_1(1)-1  10+1.2  0 1 ])
    legend('PDAS-\Gamma_1','PDAS-\Gamma_2') 
end




