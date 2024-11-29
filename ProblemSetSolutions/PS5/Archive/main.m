%------------------------------------------------------------------------%
% 
% This code computes the stationary equilibrium of a model similar to
% Aiyagari (1994) using linear grids. Labor supply is inelastic 
% and earnings risk follows a two-state Markov chain.
% 
%------------------------------------------------------------------------%

clear all; close all; clc;

fprintf('Running algorithm:\n')
run_time = tic;
addpath(genpath('/Users/muhammadbashir/GitHub/MuhammadCourses/SparseEcon/lib/'))
%% 1. Load parameters
param = define_parameters();

%% 2. Setup grid
k = linspace(param.kmin, param.kmax, param.I)';
dk = (param.kmax-param.kmin)/(param.I-1);
kk = [k, k];

%% 3. Initial guess
r0 = 0.03;  % initial interest rate guess
w0 = ((1-param.alpha)*param.A*(param.alpha*param.A/(r0 + param.delta))^(param.alpha/(1-param.alpha)));

fprintf('Starting iterations for equilibrium interest rate...\n')

%% 4. Solve for equilibrium
for iter = 1:param.maxit
    % Get wage from interest rate
    w = ((1-param.alpha)*param.A*(param.alpha*param.A/(r0 + param.delta))^(param.alpha/(1-param.alpha)));
    
    % Solve household problem
    [V, policy] = VFI(kk, r0, w, param);
    
    % Find stationary distribution
    g = KF(policy.s, param);
    
    % Calculate aggregate capital supply and demand
    K_supply = sum(sum(k.*g))*dk;
    K_demand = ((r0 + param.delta)/(param.alpha*param.A))^(1/(param.alpha-1));
    
    % Check market clearing
    diff = K_supply - K_demand;
    
    fprintf('Iter %d: r = %.4f, K_supply = %.4f, K_demand = %.4f, diff = %.4f\n', ...
            iter, r0, K_supply, K_demand, diff);
    
    if abs(diff) < param.crit
        break;
    end
    
    % Update interest rate with dampening
    r0 = r0 - 0.1*diff;
end

%% 5. Store equilibrium results
eq.r = r0;
eq.w = w;
eq.K = K_supply;
eq.V = V;
eq.g = g;
eq.policy = policy;

%% 6. Output results
run_time = toc(run_time); 
fprintf('\nAlgorithm converged. Run-time: %.2f seconds.\n', run_time);
fprintf('Equilibrium: r = %.4f, K = %.4f\n', eq.r, eq.K);

%% 7. Plot results
fprintf('\nPlotting figures...\n');

% Value functions
figure('Name', 'Value Functions');
plot(k, V(:,1), 'b-', 'LineWidth', 2); hold on;
plot(k, V(:,2), 'r--', 'LineWidth', 2);
xlabel('Capital (k)');
ylabel('Value Function');
legend('Low Income State', 'High Income State', 'Location', 'southeast');
title('Value Functions');
saveas(gcf, 'output/value_functions.png');

% Policy functions
figure('Name', 'Policy Functions');
plot(k, policy.c(:,1), 'b-', 'LineWidth', 2); hold on;
plot(k, policy.c(:,2), 'r--', 'LineWidth', 2);
xlabel('Capital (k)');
ylabel('Consumption');
legend('Low Income State', 'High Income State', 'Location', 'southeast');
title('Consumption Policy Functions');
saveas(gcf, 'output/policy_functions.png');

% Stationary distribution
figure('Name', 'Stationary Distribution');
plot(k, g(:,1), 'b-', 'LineWidth', 2); hold on;
plot(k, g(:,2), 'r--', 'LineWidth', 2);
xlabel('Capital (k)');
ylabel('Density');
legend('Low Income State', 'High Income State', 'Location', 'northeast');
title('Stationary Distribution');
saveas(gcf, 'output/stationary_dist.png');

fprintf('Figures saved in output directory.\n');