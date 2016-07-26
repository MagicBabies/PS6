%% Problem Set 6: Learning and Memory
% *Computer Modeling of the Brain 2016*
% 
% *Due July 26*
% Zehua Li, Xihang Chen, Yedidya Moise, Sherry Shi 

%% Section 1: Stable and Unstable Equilibrium 
% 1) Phosphorylation 
% Written by Yedidya Moise rechecked by Zehua Li 

n=(0:1:100)/100; %this sets up n to be a fraction of the subunits

A=3.9; %setting up variables according to the equations given in the assignment
alpha=A*(1-n).*n; %setting up the equation

figure(1)
plot(n,alpha,'b')

% 2) Desphosphorylation 
% Written by Yedidya Moise rechecked by Zehua Li 

B=.1; %setting up variables according to the equations given in the assignment
beta=n./(B+n); %setting up the equation

hold on
plot(n,beta,'r')

% 3) Creating Stable States 
% Written by Yedidya Moise 
dn=alpha-beta; %setting up the change in n to be according to the equation in the assignment
figure(2)
plot(n,dn,'c')

r=zeros(size(n)); %this creates a straight line, when dn intercepts with it
%that is the equilibrium point
hold on
plot(n,r,'m')

%the equilibrium points are when n= 0, 0.2353, and 0.664

% 4) Modeling Continuous Change 
clear all
close all
dt=.01; t=0:dt:100; %setting up time
n=zeros(size(t)); %setting up a vector for the value of n at each point
%n(1)=.23531; %starting n being equal to 
n(1)=.25;
%figurePS6 is the graph for when n is a variety of values
%setting up a for loop
A=3.9; %setting up variables
B=.1; %setting up variables
%setting up a for loop that calculates alpha beta and dn at every iteration
for i=2:length(t)
    alpha=A*(1-n(i-1))*n(i-1);
    beta=n(i-1)/(B+n(i-1));
    dn=alpha-beta;
    n(i)=n(i-1) + dt*dn; %predicts the value of n(i) using the previous value
    %rate of change and dt
end
figure(3)
hold on
plot(t,n)

%n goes to the two extreme stable equilibrium depending if it above or below the middle
% equilibrium point. The stable equilibrium points are 0 and 0.664. The
% unstable equilibrium point is 0.2353.


%% Section 2. Stochastic Models
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

%% 2.4 
dt=.01; maxt=2000; t=0:dt:maxt; %setting up time

n=zeros(size(t));
subunits=8;
m=zeros(1,subunits); %all subunits dephosphorilated
n(1)=sum(m)/subunits; %numbers that are ones divided by numbers that are zeros
%so it is the initial proportions of subunits open
 
B=.000000001;
A=6.3;
s=.0005; %setting up the value for the new variable 
 
  hi=zeros(size(n));
  yo=zeros(size(n));
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
    hi(nindex)=alpha;
    yo(nindex)=beta;
    
end
figure(1)
clf
plot(t,n)  


%%
%2.5) 
dt=.01; maxt=1500; t=0:dt:maxt; %setting up time 

n=zeros(size(t));
subunits=100;
m=zeros(1,subunits); %all subunits dephosphorilated
n(1)=sum(m)/subunits; %numbers that are ones divided by numbers that are zeros
%so it is the initial proportions of subunits open
 
 B=.0000001;
 A=4.9;
 s=.0005; %setting up the value for the new variable 
 
  hi=zeros(size(n));
  yo=zeros(size(n));
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
    hi(nindex)=alpha;
    yo(nindex)=beta;
    
end
figure(2)
clf
plot(t,n)  

%The range of fluctuations remains between .3 and .7, very stable
%There was no shift between UP and DOWN 

%% Extra Credit 6: Calcium Influx (Long term potentiation))
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
        elseif popen>rand(1); % If the subunit is already closed
            m(mindex)=1; % I open it
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
        elseif popen>rand(1); % If the subunit is already closed
            m(mindex)=1; % I open it
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

%% Extra Credit 7: PP1 (Long term depression)
dt=.01; maxt=1500; t=0:dt:maxt; %setting up time 

subunits=100;
m=ones(1,subunits); %all subunits dephosphorilated
n=ones(1, length(t)); %numbers that are ones divided by numbers that are zeros
% so it is the initial proportions of subunits open
 
B=0.1;
A=3;
s=.0005; %setting up the value for the new variable 
 
hi=zeros(size(n));
yo=zeros(size(n));
for nindex=2:length(n)
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
    hi(nindex)=alpha;
    yo(nindex)=beta;
    
    if nindex > length(t) / 2
        B=.001;
    end
end
figure(2)
clf
plot(t, n, 'r')