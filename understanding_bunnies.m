clear all
close all
clc

syms x(t) y(t)

%%
% increasing bunnies with initial condition and no "harvesting" rate
dx = diff(x, t, 1);
a = 0.1;
x0 = 5;
bunnies1DE = dx == a*x;
x1(t) = dsolve(bunnies1DE, x(0)==x0)
fplot(x1(t),[0 35])
grid on
title('Plotting bunnies')
xlabel('time');
ylabel('bunnies');
legend('x1(t) = 5*exp(t/10)');


%%
% increasing bunnies with initial condition and fixed "harvesting" rate
dx = diff(x, t, 1);
a = 0.1;
x0 = 5;
h1 = 0.6; % bunnies die off (harvesting rate is too strong)
h2 = 0.5; % bunnies stay at a population of 5 for all time
h3 = 0.4; % bunnies grow to infinity (harvesting rate is too weak)
bunnies2DE = dx == a*x - h3;
x2(t) = dsolve(bunnies2DE, x(0)==x0)
fplot(x2(t),[0 35])
grid on
title('Plotting bunnies');
xlabel('time');
ylabel('bunnies');
legend('x2(t) = exp(t/10) + 4'); % this will vary depending on which h you pick


%%
% logistic equation
dx = diff(x, t, 1);
a = 0.2;
x0 = 2;
K = 35 % carrying capacity
bunnies3DE = dx == a*x*(1 - x/K);
x3(t) = dsolve(bunnies3DE, x(0)==x0)
fplot(x3(t),[0 50])
grid on
title('Plotting bunnies')
xlabel('time');
ylabel('bunnies');
legend('x3(t) = 35/(exp(log(33/2) - t/5) + 1)')


%%
% adding foxes (no fox harvesting rate)
dy = diff(y, t, 1);
b = 0.5;
y0 = 2;
foxes1DE = dy == b*y;
foxes_initial = y(0) == y0;

% putting the effect of foxes on our bunnies
dx = diff(x, t, 1);
a = 0.1;
x0 = 5;
bunnies_initial = x(0) == x0;
h4 = 0.1*y; % A harvesting bunnies rate that is effected by how many foxes there are
bunnies4DE = dx == a*x - h4;

% solving the system of equations (no fox harvesting rate)
initial_conditions = [foxes_initial, bunnies_initial];
odes = [foxes1DE; bunnies4DE];
[x4(t), y1(t)] = dsolve(odes, initial_conditions) %x4(t) is bunnies, y1(t) is foxes

% plotting
fplot(@(t) x4(t), [0,5])
hold on
fplot(@(t) y1(t), [0,5])
hold off
grid on
title('Plotting Bunnies and Foxes')
xlabel('time');
ylabel('Bunnies and Foxes Populations');
legend('x4(t) = (11*exp(t/10))/2 - exp(t/2)/2', 'y1(t) = 2*exp(t/2)')


%%
% adding foxes (included fox harvesting rate)
dy = diff(y, t, 1);
b = 0.5;
y0 = 2;
k1 = 0.1*x; % harvesting foxes (if there are no bunnies, foxes starve and die off)
foxes2DE = dy == b*y - k1;
foxes_initial = y(0) == y0;

% putting the effect of foxes on our bunnies
dx = diff(x, t, 1);
a = 0.1;
x0 = 5;
bunnies_initial = x(0) == x0;
h4 = 0.1*y; % A harvesting bunnies rate that is effected by how many foxes there are
bunnies5DE = dx == a*x - h4;

% solving the system of equations (included fox harvesting rate)
initial_conditions = [foxes_initial, bunnies_initial];
odes = [foxes2DE; bunnies5DE];
[x5(t), y2(t)] = dsolve(odes, initial_conditions) %x5(t) is bunnies, y2(t) is foxes

% plotting
fplot(@(t) x5(t), [0,5])
hold on
fplot(@(t) y2(t), [0,5])
hold off
grid on
title('Plotting Bunnies and Foxes')
xlabel('time');
ylabel('Bunnies and Foxes Populations');
legend('x4(t)', 'y1(t)')


%%
% ideal case
syms p(t) m(t) l(t) T Y
dx = diff(x, t, 1)
dy = diff(y, t, 1)
Eqns = [dx == (450 - 6*x - 3*y)*x; 
        dy == (800 - 4*x - 12*y)*y];
[DEsys,Subs] = odeToVectorField(Eqns);
DEFcn = matlabFunction(DEsys, 'Vars',{T,Y});
tspan = [0,0.1];
initial_conditions = [49 49];
[T,Y] = ode45(DEFcn, tspan, initial_conditions)
figure(1)
plot(T,Y)
legend('x(t)','y(t)')
grid