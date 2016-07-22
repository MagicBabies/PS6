%% 
%Yedidya Moise
n=(0:1:100)/100;
A=3.9;
alpha=A*(1-n).*n;

plot(n,alpha,'b')

B=.1;
beta=n./(B+n);

hold on
plot(n,beta,'r')

%%
dn=alpha-beta;
plot(n,dn,'c')

r=zeros(size(n));
plot(n,r,'m')

%the equilibrium points are when n= 0, 0.2353, and 0.664

%%
clear all

dt=.01;
t=0:dt:10;
n=zeros(size(t));
n(1)=.23531;
for i=2:length(t)
    A=3.9;
    alpha=A*(1-n(i-1))*n(i-1);
    B=.1;
    beta=n(i-1)/(B+n(i-1));
    dn=alpha-beta;
    n(i)=n(i-1) + dt*dn;
end
figure(2)
hold on
plot(t,n)

%n goes to the two extreme equilibrium depending if it above or below the middle
% equilibrium point. The stable equilibrium points are 0 and 0.664. The
% unstable equilibrium point is 0.2353.