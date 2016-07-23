% 1. Phosphorylation **********
% Written by Yedidya Moise. Rechecked by Zehua Li.
clear all

dt = .01;
n = 0:dt:1;
A = 3.9;
alpha = A * (1 - n) .* n;

figure(1)
clf
plot(n, alpha, 'b')

% 2. Dephosphorylation **********
% Written by Yedidya Moise. Rechecked by Zehua Li.

B = .1;
beta = n ./ (B + n);

hold on
plot(n, beta, 'r')

%%
dn = alpha - beta;
plot(n, dn, 'c')

r=zeros(size(n));
plot(n, r, 'm')

%the equilibrium points are when n= 0, 0.2353, and 0.664

%%
clear all

dt = .01;
t = 0:dt:10;
n = zeros(size(t));
n(1) = .23531;
for i = 2:length(t)
    A = 3.9;
    alpha = A * (1 - n(i - 1)) * n(i - 1);
    B = .1;
    beta = n(i - 1) / (B + n(i - 1));
    dn = (alpha - beta) * dt;
    n(i) = n(i - 1) + dn;
end
figure(2)
hold on
plot(t, n)

% n goes to the two extreme equilibrium depending if it above or below the middle
% equilibrium point. The stable equilibrium points are 0 and 0.664. The
% unstable equilibrium point is 0.2353.

%% Section 2. Stochastic Models
% 1. Stochastic Model **********
% Written by Xinhang Chen. Rechecked by Zehua Li.
% initial n = 12
clear all
dt = .01;
t = 0:dt:10;
n = zeros(size(t));
p = zeros(size(t));
n(1) = 66.4;
for i = 2:length(t)
    A = 1;
    alpha = A * (100 - n(i - 1)) * n(i - 1);     
   
    B = 1;
    beta = (n(i - 1) / (B + n(i - 1))) * 100;
    p(i) = (alpha - beta) * dt;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This part seems weird. Each subunit should be randomly
    % phosphorylated or dephosphorylated at every time point.
    %                                             Simon
    if p(i) > rand(1)
        n(i) = n(i - 1) + 1;
    else
        n(i) = n(i - 1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
figure(3)
plot(t, n)
hold all
plot(t, p)
% initial n = 20
clear all
dt = .01;
t = 0:dt:10;
n = zeros(size(t));
p = zeros(size(t));
n(1) = 23.53;
for i = 2:length(t)
    A = 1;
    alpha = A * (100 - n(i - 1)) * n(i - 1);     
   
    B = 1;
    beta = (n(i - 1) / (B + n(i - 1))) * 100;
    p(i) = (alpha - beta) * dt;
    if p(i) >= rand(1)
        n(i) = n(i - 1) + 1;
    else
        n(i) = n(i - 1);
    end
end
figure(4)
plot(t, n, 'b')
hold all
plot(t, p, 'r')
%continue to increase until reach 100 subunits