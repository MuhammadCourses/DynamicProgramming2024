%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code: HJB_Huggett_assetsupply
% 
% Author: Kiyea Jin
% Date: Nov 7, 2024
%
% Description:
% This MATLAB function computes the asset supply S(r) for a given interest rate r.
% It solves the HJB equation and the KF equation of the Huggett model.
%
% Reference: Huggett_asset_supply.m by Benjamin Moll
%
% Notes:
% - CRRA utility function: U(c) = (c^(1-sigma))/(1-sigma)
% - Discount rate (rho): 0.05
% - Income: z = [z_u, z_e] = [0.1, 0.2];
% - Lambda: la = [la_u, la_e] = [1.2, 1.2];
% - Asset grid (a): -0.15 to 5
% - Borrowing constraint: a >= -0.15
% - Delta = 1000 (for implicit method)
%
% Code Structure:
% 1. VALUE FUNCTION ITERATION
% 2. KF EQUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = Huggett_assetsupply(r, p, aa, zz, Df, Db, A_switch, da)

    % Initial guess for the value function
    v0 = p.u(zz + max(r, 0.01) .* aa) ./ p.rho; 
    V = v0;

    % 1. VALUE FUNCTION ITERATION
    for n = 1:p.maxit

        % 1-1. Compute derivatives of the value function
        dVf = Df * V;
        dVb = Db * V;

        % 1-2. Boundary conditions
        dVb(1, :) = p.mu(zz(1, :) + r .* aa(1, :));      % Enforce borrowing constraint a >= a_min
        dVf(end, :) = p.mu(zz(end, :) + r .* aa(end, :)); % Enforce upper bound a <= a_max

        % 1-3. Compute optimal consumption
        cf = p.inv_mu(dVf);
        cb = p.inv_mu(dVb);

        % 1-4. Compute optimal savings
        sf = zz + r .* aa - cf;
        sb = zz + r .* aa - cb;

        % 1-5. Upwind scheme
        If = sf > 0;
        Ib = sb < 0;
        I0 = 1 - If - Ib;
        dV_Upwind = If .* dVf + Ib .* dVb + I0 .* p.mu(zz + r .* aa);
        c = p.inv_mu(dV_Upwind);

        % 1-6. Construct the matrix for the linear system
        c_stacked = c(:);
        V_stacked = V(:);

        % Drift terms
        A1 = spdiags(If(:, 1) .* sf(:, 1), 0, p.I, p.I) * Df + spdiags(Ib(:, 1) .* sb(:, 1), 0, p.I, p.I) * Db;
        A2 = spdiags(If(:, 2) .* sf(:, 2), 0, p.I, p.I) * Df + spdiags(Ib(:, 2) .* sb(:, 2), 0, p.I, p.I) * Db;
        A = [A1, sparse(p.I, p.I); sparse(p.I, p.I), A2];

        % Total generator
        P = A + A_switch;

        % Right-hand side and left-hand side for the linear system
        B = (p.rho + 1 / p.Delta) * speye(2 * p.I) - P;
        b = p.u(c_stacked) + V_stacked / p.Delta;

        % Solve for the updated value function
        V_new = B \ b;
        V_change = max(abs(V_new - V_stacked));

        % Update value function
        V = reshape(V_new, p.I, 2);

        % Check for convergence
        if V_change < p.tol
            break;
        end
    end

    % 2. KF EQUATION

    % Transpose of generator matrix
    P_transpose = P';

    % Solve P'*g = 0 with normalization condition
    g_stacked = zeros(2 * p.I, 1);
    i_fix = 1;
    P_transpose(i_fix, :) = zeros(1, 2 * p.I);
    P_transpose(i_fix, i_fix) = 1;
    g_stacked(i_fix) = 1;

    % Solve for the stationary distribution
    g = P_transpose \ g_stacked;

    % Normalize the distribution
    g = g / sum(g * da);

    % Reshape to original dimensions
    g = reshape(g, p.I, 2);

    % Compute asset supply S(r)
    S = sum(g(:, 1) .* aa(:, 1) * da) + sum(g(:, 2) .* aa(:, 2) * da);

end
