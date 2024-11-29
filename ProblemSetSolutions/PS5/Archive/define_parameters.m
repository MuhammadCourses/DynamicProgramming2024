function [param] = define_parameters()
    % Economic Parameters
    param.A = 0.1;           % total factor productivity
    param.alpha = 1/3;       % capital share
    param.delta = 0.05;      % depreciation rate
    param.rho = 0.05;        % discount rate
    param.gamma = 2;         % risk aversion
    
    % Grid Parameters
    param.kmin = 0;        % minimum capital
    param.kmax = 20;         % maximum capital
    param.I = 1000;          % number of grid points
    param.k = linspace(param.kmin, param.kmax, param.I)';  % capital grid
    
    % Labor Market Parameters
    param.la1 = 1/3;         % transition rate 1
    param.la2 = 1/3;         % transition rate 2
    param.zz = [1, 2];   % productivity states
    param.S = 2;             % number of states
    
    % Utility Functions
    param.u = @(c) (c.^(1-param.gamma)-1)/(1-param.gamma);    % utility
    param.u1 = @(c) c.^(-param.gamma);                        % marginal utility
    param.u1inv = @(x) x.^(-1/param.gamma);                   % inverse marginal utility
    
    % Iteration Parameters
    param.maxit = 100;       % maximum iterations
    param.crit = 1e-6;       % convergence criterion
    param.Delta = 1000;      % time step for implicit method
end