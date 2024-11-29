

% Author:  Muhammad

function p = define_parameters_GE()
%% Economic Parameters
% Relative risk aversion
p.gamma = 2;
% Discount rate
p.rho = 0.05;
%% WE NO LONGER ASSUME EXOGENOUS INTEREST RATE
% Exogenous interest rate
% p . r = 0.035;
% Income process
p.z_u = 1;
p.z_e = 2;
p.zz = [ p . z_u , p . z_e ];
% P r o b a b i l i t y density
p.lambda_u = 1/3;
p.lambda_e = 1/3;
p.lambda = [ p.lambda_u , p.lambda_e ];
%% Economic Functions
% Utility funtion
p.u = @ ( c ) c.^(1 - p.gamma )/(1 - p.gamma );
% Marginal utility function
p.mu = @ ( c ) c.^(- p.gamma );
% FOC : mu ( c )= dV -> c = inv_mu ( dV )
p.inv_mu = @ ( dV ) dV.^(-1/p.gamma );
% Production function
p.alpha = 1/3;
p.delta = 0.05;
p.A = 0.1;
p.f = @ ( x ) p.A*x(1)^p.alpha*x(2)^(1 - p.alpha);
% Marginal product of capital
p.fk = @ ( x ) p.alpha*p.A*x(1)^(p.alpha - 1)*x(2)^(1 - p.alpha);
% Marginal product of labor
p.fl = @ ( x ) (1 - p.alpha)*p.A*x(1)^p.alpha*x(2)^(-p.alpha); 
%% Grid Parmaters
p.kmin = 0;
p.kmax = 20;
% The number of grid points
p.I = 1000;
%% Tuning parameters
% Step size : can be a r b i t r a r i l y large in implicit method
p.Delta = 1000;
% The maximum number of value function it er at io ns
p.maxit = 100;
% Tolerance for value function it er at io ns
p.tol = 10^( -6);
%% TUNING P AR AM ET ER S FOR INTEREST RATE ITERATION
% The maximum number of interest rate it er at io ns
p.Nr = 40;
% Tolerance for interest rate it er at io ns
p.tol_S = 10^( -6);

end 