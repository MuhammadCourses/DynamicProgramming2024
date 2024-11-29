


function [S, ss] = stationary(r0, p)
% This function solves the HJB equation and KF equation of the Huggett model
% for given interest rates r0.

%% 2. INITIALIZE GRID POINTS
k = linspace(p.kmin, p.kmax, p.I)';
ka = (p.kmax - p.kmin) / (p.I - 1);
kk = [k, k]; % I * 2 matrix

%% 3. PRE-ITERATION INITIALIZATION
% 3-1. Construct the forward and backward differential operator
% Df such that Df * V = dVf and Db such that Db * V = dVb
Df = zeros(p.I, p.I);
for i = 1:p.I-1
    Df(i, i) = -1 / da;
    Df(i, i + 1) = 1 / da;
end
Df = sparse(Df);

Db = zeros(p.I, p.I);
for i = 2:p.I
    Db(i, i - 1) = -1 / da;
    Db(i, i) = 1 / da;
end
Db = sparse(Db);

% 3-2. Construct A_switch matrix
A_switch = [speye(p.I) * (-p.lambda(1)), speye(p.I) * p.lambda(2);
            speye(p.I) * p.lambda(1), speye(p.I) * (-p.lambda(2))];

% 3-3. Guess an initial value of the value function
zz = ones(p.I, 1) * p.zz; % I * 2 matrix
% The value function of "staying put"
r = r0;
v0 = p.u(zz + r .* aa) / p.rho;
v = v0;

%% 4. VALUE FUNCTION ITERATION
for n = 1:p.maxit
    V = v;
    V_u = V(:, 1);
    V_e = V(:, 2);

    % 4-1. Compute the derivative of the value function
    dVf = [Df * V_u, Df * V_e];
    dVb = [Db * V_u, Db * V_e];

    % 4-2. Boundary conditions
    dVb(1, :) = p.mu(zz(1, :) + r .* aa(1, :)); % a >= a_min is enforced (borrowing constraint)
    dVf(end, :) = p.mu(zz(end, :) + r .* aa(end, :)); % a <= a_max is enforced

    I_concave = dVb > dVf; % indicator whether value function is concave

    % 4-3. Compute the optimal consumption
    cf = p.inv_mu(dVf);
    cb = p.inv_mu(dVb);

    % 4-4. Compute the optimal savings
    sf = zz + r .* aa - cf;
    sb = zz + r .* aa - cb;

    % 4-5. Upwind scheme
    If = sf > 0;
    Ib = sb < 0;
    I0 = 1 - If - Ib;
    dV0 = p.mu(zz + r .* aa); % If sf <= 0 <= sb, set s = 0
    dV_upwind = If .* dVf + Ib .* dVb + I0 .* dV0;
    c = p.inv_mu(dV_upwind);
    s = zz + r .* aa - c;

    % 4-6. Update value function
    V_stacked = V(:); % 2I * 1 matrix
    c_stacked = c(:); % 2I * 1 matrix

    % A = SD
    Su = spdiags(If(:, 1) .* sf(:, 1), 0, p.I, p.I) * Df + spdiags(Ib(:, 1) .* sb(:, 1), 0, p.I, p.I) * Db;
    Se = spdiags(If(:, 2) .* sf(:, 2), 0, p.I, p.I) * Df + spdiags(Ib(:, 2) .* sb(:, 2), 0, p.I, p.I) * Db;
    A = [Su, sparse(p.I, p.I); sparse(p.I, p.I), Se]; % 2I * 2I matrix

    % P = A + A_switch
    P = A + A_switch;

    % B = [(rho + 1/Delta) * I - P]
    B = (p.rho + 1 / p.Delta) * speye(2 * p.I) - P;

    % b = u(c) + 1/Delta * V
    b = p.u(c_stacked) + (1 / p.Delta) * V_stacked;

    % Solve for the new value function
    V_update = B \ b; % 2I * 1 matrix
    V_change = V_update - V_stacked;
    v = reshape(V_update, p.I, 2); % I * 2 matrix

    % 3-6. Convergence criterion
    dist(n) = max(abs(V_change));
    if dist(n) < p.tol
        disp('Value function converged. Iteration =')
        disp(n)
        break
    end
end
toc;

%% 5. KF EQUATION
% 5-1. Solve for 0 = gdot = P' * g
PT = P';
gdot_stacked = zeros(2 * p.I, 1);
% need to fix one value, otherwise matrix is singular
i_fix = 1;
gdot_stacked(i_fix) = .1;
row_fix = [zeros(1, i_fix - 1), 1, zeros(1, 2 * p.I - i_fix)];
PT(i_fix, :) = row_fix;
g_stacked = PT \ gdot_stacked;

% 5-2. Normalization
g_sum = g_stacked' * ones(2 * p.I, 1) * da;
g_stacked = g_stacked / g_sum;

% 5-3. Reshape
gg = reshape(g_stacked, p.I, 2);

%% OUTPUT
S = gg(:, 1)' * a * da + gg(:, 2)' * a * da;
ss.a = a;
ss.gg = gg;
ss.adot = zz + r .* aa - c;
ss.V = V;
ss.dV = dV_upwind;
ss.c = c;
end