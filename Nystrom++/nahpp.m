function [trest] = nahpp(n,Afun,m)

S = randn(n,round(m/4)); R = randn(n,round(m/2)); G = randn(n,round(m/4));

Z = Afun(R); W = Afun(S);

M = S'*Z; 
[Q,U] = qr(M',0);
X = Z*Q; Y = (W/U)'; 

trest = trace(Y*X) + 4*(trace(G'*Afun(G)) - trace((G'*X)*(Y*G)))/m;

% S = randn(n,round(m/4)); R = randn(n,round(m/2)); G = randn(n,round(m/4));
% 
% Z = Afun(R); W = Afun(S);
% 
% M = pinv(S'*Z);
% 
% trest = trace(M*(W'*Z)) + 4*(trace(G'*Afun(G)) - trace((G'*Z)*(M*(W'*G))))/m;
%appr_err = trace(Afun(W)*(M'\Z'))+trace(Afun(Z)*(M\W'))-norm(Z*(M\W'),'fro')^2;

end