%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: main_GE
% 
% Author: Kiyea Jin
% Date: Nov 16, 2024
%
% Description:
% This MATLAB script solves the general equilibrium of the Huggett model.
% The implementation leverages the SparseEcon repository developed by 
% Andreas Schaab and Allen T. Zhang (available at https://github.com/schaab-lab/SparseEcon).
%
% Reference:
% - Huggett_equilibrium_iterate.m by Benjamin Moll
% - Codes developed during the 2024 Coding Workshop by Andreas Schaab
%
% Notes:
% - CRRA utility function: U(c) = (c^(1-sigma))/(1-sigma)
% - Elasticity of intertemporal substitution (sigma): 2
% - Discount rate (rho): 0.05
% - Income: z = [z_u, z_e] = [0.1, 0.2];
% - Lambda: la = [la_u, la_e] = [1.5, 1];
% - Discrete grid of asset levels (a): -0.15 to 5
% - Borrowing constraint: a>=-0.15
% - Delta = 1000; (Can be arbitrarily large in implicit method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;
clc;

%% 0. Add the path of SparseEcon folder

addpath(genpath('../SparseEcon/lib'));
% addpath(genpath('SparseEcon/lib'));

%% 1. DEFINE PARAMETERS

p = define_parameters();

%% 2. INITIALIZE GRID POINTS

% a = linspace(p.kmin, p.kmax, p.I)';
% da = (p.kmax-p.kmin)/(p.I-1);
% 
% aa = [a, a]; % I*2 matrix

%{
setup_grid.m: Initializes grid structure

function G = setup_grid(n, surplus, min, max, varargin)

INPUTS:
- n:       Level of sparse grid
- surplus: Dimension-specific level surplus e.g., [0, 1, -1]
- min:     Dimension-specific minimums in economic units
- max:     Dimension-specific maximums in economic units

VARIABLE INPUTS:
- NamedDims: Cell of vectors of named dimensions
- Names:     Names for named dimensions
- DxxDims:   Vector of dimensions to compute dxx operators
- DxyDims:   (n x 2) matrix of dimensions to compute dxy operators

%}

G = setup_grid(p.l, 0, p.kmin, p.kmax, 'NamedDims', {1}, 'Names', {'k'});

% Notes:  G.a is the grid points. G.J is the number of grid points.

%% NEWTON'S METHOD

% Guess an initial value of the interest rate
r0 = p.r;

%% SOLVE STATIONARY EQUILIBRIUM
K0 = 15; L0 = 1.2; %p.L;   % equally weight all observations. Here in first iteration I equally weight all observations i.,e every household supplies one unit inelastic labor.
; if exist('ss', 'var'), K0 = ss.K; L0 = ss.L; end; X0 = [K0, L0]; J0 = [];
% K0 = 10; if exist('ss', 'var'), K0 = ss.K; end; X0 = K0; J0 = [];       % because labor supply is inelastic

    % % Get better guess for value function:
    [diff,S, ss,G] = stationary(X0, G, p);
    
    % Solve for steady state prices: 
    % f = @(x, y) stationary(x, y, GDense, param); y0 = G;
    % [X, J0] = fsolve_newton(f, reshape(X0, [numel(X0), 1]), diff0, y0, J0, 5, 0);
    options = optimset('Display', 'off', 'UseParallel', true, 'TolX', 1e-12);
    X = fsolve(@(x) stationary(x,G, p), X0, options);

% % Rerun stationary using the equilibrium prices:
    [diff,S, ss,G] = stationary(X, G, p);

%% 6. GRAPHS 

% 6-1. Optimal consumption 
figure;
set(gca, 'FontSize', 18)
plot(G.k, ss.c, 'LineWidth', 2)
grid
xlabel('Capital, k','FontSize', 14)
ylabel('Consumption, c_j(k)','FontSize', 14)
xlim([p.kmin p.kmax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)
% save figure 
saveas(gcf, fullfile('..', 'Output', ['GE_optimal_consumption.png']));


% 6-2. Optimal Investment 
G.income = ss.w.*p.zz + ss.r.*G.k; 
kdot = G.income - ss.c;
figure;
set(gca, 'FontSize', 18)
plot(G.k, kdot, G.k, zeros(1,G.J), '--k', 'LineWidth', 2)
grid
xlabel('Capital, k', 'FontSize', 14)
ylabel('Investment, I_j(k)', 'FontSize', 14)
xlim([p.kmin p.kmax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)
% save figure
saveas(gcf, fullfile('..', 'Output', ['GE_optimal_investment.png']));

% 6-3. Value function
figure;
set(gca, 'FontSize', 18)
plot(G.k, ss.V, 'LineWidth', 2)
grid
xlabel('Capital, k', 'FontSize', 14)
ylabel('Value function, V_j(a)', 'FontSize', 14)
xlim([p.kmin p.kmax])
legend('Unemployed', 'Employed', 'Location', 'best', 'FontSize', 14)
% save figure
saveas(gcf, fullfile('..', 'Output', ['GE_value_function.png']));

% 6-4. Capital distribution
figure;
set(gca, 'FontSize', 14)
plot(G.k, ss.gg, 'LineWidth', 2)
grid
xlabel('Capital, k', 'FontSize', 14)
ylabel('Densities, g_j(k)', 'FontSize', 14)
yy = get(gca, 'yLim');
hold on
plot([0,0], yy, '--k', 'LineWidth', 2)
xlim([p.kmin p.kmax])
legend('Unemployed', 'Employed', 'No Shortsell Constraint', 'Location', 'best', 'FontSize', 14)
% save figure
saveas(gcf, fullfile('..', 'Output', ['GE_wealth_distribution.png']));