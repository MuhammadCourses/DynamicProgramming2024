% MATLAB Code: HJB_Huggett_GE_Newton
%
% Author: Kiyea Jin
% Date: Nov 7, 2024
%
% Description:
% This MATLAB script solves the general equilibrium of the Huggett model,
% finding the equilibrium interest rate that clears the bond market using Newton’s method.
%
% Reference: Huggett_equilibrium_iterate.m by Benjamin Moll
%
% Notes:
% - CRRA utility function: U(c) = (c^(1 - sigma))/(1 - sigma)
% - Elasticity of intertemporal substitution (sigma): 2
% - Discount rate (rho): 0.05
% - Income: z = [z_u, z_e] = [0.1, 0.2];
% - Lambda: la = [la_u, la_e] = [1.2, 1.2];
% - Discrete grid of asset levels (a): -0.15 to 5
% - Borrowing constraint: a >= -0.15
% - Delta = 1000; (Can be arbitrarily large in implicit method)
%
% Code Structure:
% 1. DEFINE PARAMETERS
% 2. NEWTON’S METHOD
% 3. GRAPHS
clear all;
close all;
clc;

%% 1. DEFINE PARAMETERS
p = define_parameters_GE();

%% 2. NEWTON’S METHOD
% Guess an initial value of the interest rate
r0 = 0.03;

% Use fsolve to find the equilibrium interest rate that makes S(r) = 0
% Note: X = fsolve(FUN, X0, options) starts at the matrix X0 and tries to solve the equation
% Set OPTIONS = optimoptions('fsolve', 'Algorithm', 'trust-region'), and then pass OPTIONS
S = @(r) stationary(r, p);
options = optimoptions('fsolve', 'TolFun', p.tol, 'Display', 'iter');
r_equilibrium = fsolve(S, r0, options);

% Rerun stationary with the correct interest rate
[~, ss] = stationary(r_equilibrium, p);

%% 6. GRAPHS
% 6-1. Optimal consumption
set(gca, 'FontSize', 18)
plot(ss.a, ss.c(:,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(ss.a, ss.c(:,2), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
yy = get(gca, 'yLim');
plot([p.amin, p.amin], yy, '--k', 'LineWidth', 2)
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Consumption, c_j(a)', 'FontSize', 14)
xlim([p.amin - 0.1, p.amax])
legend(sprintf('Unemployed, r = %.4f', r_equilibrium), ...
    sprintf('Employed, r = %.4f', r_equilibrium), 'Location', 'best')
text(-0.08, 0.35, '$a = \underline{a}$', 'FontSize', 16, 'Interpreter', 'latex')

% 6-2. Optimal savings
set(gca, 'FontSize', 18)
plot(ss.a, ss.adot(:,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(ss.a, ss.adot(:,2), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
yy = get(gca, 'yLim');
plot([p.amin, p.amin], yy, '--k', 'LineWidth', 2)
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Saving, s_j(a)', 'FontSize', 14)
xlim([p.amin - 0.1, p.amax])
legend(sprintf('Unemployed, r = %.4f', r_equilibrium), ...
    sprintf('Employed, r = %.4f', r_equilibrium), 'Location', 'best')
text(-0.08, -0.1, '$a = \underline{a}$', 'FontSize', 16, 'Interpreter', 'latex')

% 6-3. Wealth distribution
set(gca, 'FontSize', 14)
plot(ss.a, ss.gg(:,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r')
hold on
plot(ss.a, ss.gg(:,2), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
hold on
yy = get(gca, 'yLim');
plot([p.amin, p.amin], yy, '--k', 'LineWidth', 2)
hold off
grid
xlabel('Wealth, a', 'FontSize', 14)
ylabel('Densities, g_j(a)', 'FontSize', 14)
xlim([p.amin - 0.05, 1])
legend(sprintf('Unemployed, r = %.4f', r_equilibrium), ...
    sprintf('Employed, r = %.4f', r_equilibrium), 'Location', 'best')
text(-0.12, 3, '$a = \underline{a}$', 'FontSize', 16, 'Interpreter', 'latex')
set(gca, 'FontSize', 14)