%% Problem Set 6: Learning and Memory
% *Computer Modeling of the Brain 2016*
% 
% *Due July 26*
% Zehua Li, Xihang Chen, Yedidya Moise, Sherry Shi 

%% Section 1: Stable and Unstable Equilibrium 

%1. Phosphorylation **********
% Written by Yedidya Moise. Rechecked by Zehua Li.
clear all
s=7
dt = .01;
n = 0:dt:1;
A = 3.9;
alpha = A * (1 - n) .*( n+s);

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

%% 3. Creating Stable States  **********
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


%%
%Yedidya Moise
%2.1)
dt=.01; maxt=10; t=0:dt:maxt; %setting up time

n=zeros(size(t)); %setting up the vector for the values of n
subunits=12; %setting up how many subunits you have
m=zeros(1,subunits); % setting up m(the subunits) as a vector (for each subunit)
%starting out as zero (starting out closed/dephosphorylated)

n(1)=sum(m)/subunits; %numbers that are ones divided by the number of subunits
%so it is the initial proportions of subunits open
 %setting up the variables for the equations
  B=.1; 
  A=3.9;

for i=1:round(subunits)/2 %starting out with the first half of subunits being
      %phosphorylated
   m(i)=1
end

  %setting up a for loop with the probabilities of subunits opening and
  %closing and with the equations for alpha and beta
  
  for nindex=2:length(t)
    alpha=A*(1-n(nindex-1))*n(nindex-1);
    beta=n(nindex-1)/(B+n(nindex-1));
    
    pclose=beta*dt; %probability of subunit closing as given in the assignment
    popen=alpha*dt; %probability of subunit opening (phosphorylation)
   %setting up a for loop with each value of the subunits 
    for mindex=1:subunits
        if m(mindex)==1 %if that subunit is open(phosphorylated) then if the 
            %probabilty for closing it is high enough it will close
            if pclose>rand(1)
            m(mindex)=0; 
            end
        else % If the subunit is already closed
            if popen>rand(1); %if the probability for opening is high enough
                m(mindex)=1; % It opens
            end
        end
    end
    n(nindex)=sum(m)/subunits;%sets up n(nindex) to be equal to numbers 
 %that are ones divided by the total number of subunits so it is the
 %proportions of subunits open
end
figure(1)
clf
plot(t,n)

%over time the number of subunits phosphorylated goes up and down according
%to the probabilities of individual subunits phosphorylating and
%dephosphorylating.

%%
%2.2)
dt=.01; maxt=500; t=0:dt:maxt; %setting up time

n=zeros(size(t));
subunits=12;
m=zeros(1,subunits);

for i=1:round(subunits) %starting out with all of the subunits being dephosphorylated
    m(i)=1
end


n(1)=sum(m)/subunits; 
  B=.05; %value of B that has 3 equilibrium points and makes all of the subunits 
  %eventually decay to zero
  A=3.9;
for nindex=2:length(t)
    alpha=A*(1-n(nindex-1))*n(nindex-1);
    beta=n(nindex-1)/(B+n(nindex-1));
    
    pclose=beta*dt; %probability of subunit closing
    popen=alpha*dt; %probability of subunit opening (phosphorylation)
   
    for mindex=1:subunits
        if m(mindex)==1
            if pclose>rand(1)
            m(mindex)=0; 
            end
        else % If the subunit is already closed
            if popen>rand(1);
                m(mindex)=1; 
            end
        end
    end
    n(nindex)=sum(m)/subunits;
end
figure(1)
clf
plot(t,n)

%%
%%PS6 s2 problem 3
dt=.01; maxt=1000; t=0:dt:maxt; %setting up time

n=zeros(size(t));
subunits=12;
m=zeros(1,subunits); %all subunits dephosphorilated
n(1)=sum(m)/subunits; %numbers that are ones divided by numbers that are zeros
%so it is the initial proportions of subunits open
  B=.05;
  A=3.9;
  s=.05; %setting up the value for the new variable 
for nindex=2:length(t)
    alpha=A*(1-n(nindex-1))*(n(nindex-1)+s); %altered equation as given in the assignment
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
%figure that shows that this s works is figurePS6section2part3
%% Extra Credit: Calcium Influx (Long term potentiation))
%created by Xinhang Chen
dt=.01; maxt=1000; t=0:dt:maxt; %setting up time

n=zeros(size(t));
subunits=1000;%use a stable number of enzymes 1000 is a very stable number
m=zeros(1,subunits); %all subunits dephosphorilated
n(1)=sum(m)/subunits; %numbers that are ones divided by numbers that are zeros
%so it is the initial proportions of subunits open
  B=.05;
  A=3.9;
  s=.05; %setting up the value for the new variable 
for nindex=2:round(length(t)/2)
    s=.01;%set s=.01 for the first half of t
    alpha=A*(1-n(nindex-1))*(n(nindex-1)+s); %altered equation as given in the assignment
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
for nindex=round((length(t)/2))+1:length(t)
    s=3;%set it equal to a value such that there is only 1 equilibrium point for the second half of t
    alpha=A*(1-n(nindex-1))*(n(nindex-1)+s); %altered equation as given in the assignment
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
%figure that shows that this works is Extra Credit 1
%Its has a relatively low percentage of phosphorylation for the first part
%of time and has a relatively high percentage of phosphrylation for the
%second part of time Via the frequent neuron stimulation, synapse between two neurons get stronger,which means the percentage of 
%phosphorylation for now(second part of time) is higher than previous time(first part of time)






