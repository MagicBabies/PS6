%% Problem Set 6: Learning and Memory
% *Computer Modeling of the Brain 2016*
% 
% *Due July 26*
% Zehua Li, Xihang Chen, Yedidya Moise, Sherry Shi 

%% Section 1: Stable and Unstable Equilibrium 

%1. Phosphorylation **********
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

title('Rate of phosphorylation and dephophorylation');
legend('Alpha','Beta');
xlabel('n (Percentage of subunits phosphorylated)');
ylabel('Rate') ;

% 3. Creating Stable States  **********
dn = alpha - beta;
figure(2)
clf
plot(n, dn, 'c')
hold on 
r=zeros(size(n));
plot(n, r, 'm')

title('Alpha-Beta as a Function of n');
xlabel('n (Percentage of subunits phosphorylated)');
ylabel('Alpha-Beta')

%the equilibrium points are when n= 0, 0.2353, and 0.664


%% 4. Modeling Continuous Change
clear all
%Calculate the value of n overtime for 3 different initial values of n 
%Choose values of n between the three equilibrium point

%First initial value of n = .2
dt = .01;
t = 0:dt:10;
n = zeros(size(t));
n(1) = .2;
for i = 2:length(t)
    A = 3.9;
    alpha = A * (1 - n(i - 1)) * n(i - 1);
    B = .1;
    beta = n(i - 1) / (B + n(i - 1));
    dn = (alpha - beta) * dt;
    n(i) = n(i - 1) + dn;
end
figure(3)
clf
plot(t, n)
%Second initial value of n = .4

n(1) = .4;
for i = 2:length(t)
    A = 3.9;
    alpha = A * (1 - n(i - 1)) * n(i - 1);
    B = .1;
    beta = n(i - 1) / (B + n(i - 1));
    dn = (alpha - beta) * dt;
    n(i) = n(i - 1) + dn;
end
hold on
plot(t, n)

%Third initial value of n = .8
n(1) = .8;
for i = 2:length(t)
    A = 3.9;
    alpha = A * (1 - n(i - 1)) * n(i - 1);
    B = .1;
    beta = n(i - 1) / (B + n(i - 1));
    dn = (alpha - beta) * dt;
    n(i) = n(i - 1) + dn;
end

hold on
plot(t, n)
legend('0.2','0.4','0.8');
title('Value of n Over Time for 3 different n(1)');
xlabel('Time');
ylabel('n (Percentage of subunits phosphorylated)');

% What happens to the value of n over time?
%    Value of n(i) between the first two equilibrium points exponentially
%      decays to the first equilibrium point, 0 
%    Value of n(i) between the 2nd and 3rd equilibrium point and above the 
%      last both stabilize at the third equilibrium point, .664

% The stable equilibrium points are 0 and 0.664. 
% The unstable equilibrium point is 0.2353.


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


%%
%Yedidya Moise
dt=.01; maxt=10; t=0:dt:maxt; %setting up time

n=zeros(size(t));
subunits=12;
m=zeros(1,subunits);

for i=1:round(subunits/2)
    m(i)=1
end


n(1)=sum(m)/subunits; %numbers that are ones divided by numbers that are zeros
%so it is the initial proportions of subunits open
  B=.1;
  A=3.9;
for nindex=2:length(t)
    alpha=A*(1-n(nindex-1))*n(nindex-1);
    beta=n(nindex-1)/(B+n(nindex-1));
    
    pclose=beta*dt; %probability of subunit closing
    popen=alpha*dt; %probability of subunit opening (phosphorylation
   
    for mindex=1:subunits
        if m(mindex)==1
            if pclose>rand(1)
            m(mindex)=0; %I close m(i)
            end
        else % If the subunit is already closed
            if popen>rand(1);
                m(mindex)=1; % I open it
            end
        end
    end
    n(nindex)=sum(m)/subunits;
end
figure(1)
clf
plot(t,n)
