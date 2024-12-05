function [diff,S, ss,G] = stationary(x, G, p)

%% Aggregates
% K = x(1);
K = x;
L= G.L;
% L = x(2);
Y = p.A*K^p.alpha*L^(1-p.alpha);
rk = p.alpha*Y/K;
w = (1-p.alpha)*Y/L;
r= rk - p.delta;
%% HJB

G.r = r;
G.w = w;

[V, c, s, P] = HJB(G, p);

%% KF

[gg] = KF(P, G);

% %% OUTPUT

% S = gg(:,1)'*G.k*G.dk + gg(:,2)'*G.k*G.dk;

% ss.V = V;
% ss.gg = gg;
% ss.c = c;
% ss.s = s;

%% MARKET CLEARING
KH = sum(sum(G.k .* gg .* G.dk));           % first sum gives household level 
LH = sum(sum(p.zz .* 1 .* gg.*G.dk));
C  = sum(sum( (c .* gg .* G.dk)));
S  = sum(sum( (s .* gg .* G.dk)));

excess_supply  = Y - C - p.delta*KH;
excess_capital = K - KH;
excess_labor   = L - LH;
% size(LH)

%diff = %[excess_capital, excess_labor];
diff = excess_capital;

ss.V = V; ss.g = gg; ss.c = c; ss.s = s;
ss.K = K;  ss.C = C; ss.S = S; ss.r = r; ss.Y = Y; ss.w = w;  % ss.L = L;
ss.excess_supply = excess_supply; ss.excess_capital = excess_capital; ss.excess_labor = excess_labor; ss.gg = gg;

end