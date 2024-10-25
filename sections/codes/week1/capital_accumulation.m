%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Script: capital_accumulation.m
%
% Author: Kiyea Jin
% Date: October 14, 2024
%
% Description:
% This script solves the continuous-time capital accumulation equation:
%   dK/dt = I - delta * K
% where {I} represents a constant, exogenously given investment, and 
% an initial condition for capital K0 is provided.
%
% The solutions are obtained using:
% 1. The integrating factor method for the analytical solution.
% 2. The finite difference method for the numerical solution.
%
% The script also:
% 3. Plots both the analytical and numerical solutions for comparison.
% 4. Explores how the numerical solution varies with grid fineness.
%
% Parameters:
% - Depreciation rate (delta): 0.05
% - Constant investment (I): 0 or 5
% - Initial boundary condition: K0 = 10
%
% Code Structure:
% 1. DEFINE PARAMETERS
% 2. INITIALIZE GRID POINTS
% 3. ANALYTICAL SOLUTION: INTEGRATING FACTOR METHOD
% 4. NUMERICAL SOLUTION: FINITE DIFFERENCE METHOD
% 5. PLOT
% (6. SOLVE USING ODE45)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
close all;
clc;

%% 1. DEFINE PARAMETERS

p = define_parameters();

%% 2. INITIALIZE GRID POINTS

t = linspace(p.tmin, p.tmax, p.I)';
dt = (p.tmax-p.tmin)/(p.I-1);

% Define investment series
I = p.Inv*ones(p.I,1);

%% 3. ANALYTICAL SOLUTION: INTEGRATING FACTOR METHOD

% 3-1. Pre-allocate arrays for solutions
K_analytical = zeros(p.I,1);

% 3-2. Analytical solution: K(t) = I(t)/delta*(1-exp^{-delta*t})+exp^{-delta*t}*K(0)
for i=1:p.I

    K_analytical(i) = (I(i)/p.delta) * (1 - exp(-p.delta*t(i))) + exp(-p.delta*t(i)) * p.K0;

end

%% 4. NUMERICAL SOLUTION: FINITE DIFFERENCE METHOD

% 4-1. Pre-allocate arrays for solutions
K_numerical  = zeros(p.I,1);
K_numerical(1) = p.K0;

% 4-2. Numerical solution: (K(t+dt)-K(t))/dt = I(t) - delta*K(t)

for i=1:p.I-1

    K_numerical(i+1) = K_numerical(i) + dt * (I(i) - p.delta*K_numerical(i));

end

%% 5. PLOT

% Compute the difference between two solutions
error = K_analytical - K_numerical;

% Plot both solutions for comparison
figure;
plot(t, K_analytical, 'r-', 'LineWidth', 2); hold on;
plot(t, K_numerical, 'bo--', 'LineWidth', 1.5);
xlabel('Time', 'FontSize', 11);
ylabel('Capital (K)', 'FontSize', 11);
legend('Analytical Solution', 'Numerical Solution (FDM)', 'FontSize', 12);
title('Comparison of Analytical and Numerical Solutions', 'FontSize', 15);
grid on;

%% (6. SOLVE USING ODE45)
% ode45  Solve non-stiff differential equations, medium order method.
% [TOUT,YOUT] = ode45(ODEFUN,TSPAN,Y0) integrates the system of
% differential equations y' = f(t,y) from time TSPAN(1) to TSPAN(end)
% with initial conditions Y0. Each row in the solution array YOUT
% corresponds to a time in the column vector TOUT. 

% Define the ODE as a function handle
ODE_capital = @(t, K) p.Inv - p.delta * K;

K_numerical2  = zeros(p.I,1);

% Call ode45 to solve the ODE
[t, K_numerical2] = ode45(ODE_capital, t, p.K0);

% Plot both solutions for comparison
figure;
plot(t, K_numerical, 'bo--', 'LineWidth', 2); hold on;
plot(t, K_numerical2, 'go--', 'LineWidth', 1.5);
xlabel('Time', 'FontSize', 11);
ylabel('Capital (K)', 'FontSize', 11);
legend('FD Scheme', 'ODE45', 'FontSize', 12);
title('Comparison of FD Scheme and ODE45', 'FontSize', 15);
grid on;
