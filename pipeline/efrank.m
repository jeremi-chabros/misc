function erank = erank(A)

% This script is an adaptation of the effective rank as by definition in
% Roy & Vetterli (2007)

% Input:
% adjM - adjacency matrix of the size N x N, where N represents the number
% of nodes in a network. adjM needs to be symmetrical and real-valued

% Output:
% erank - value of effective rank

% A(isnan(A)) = 0;
% A(1:length(A)+1:end) = 1;

[U,D,V] = svd(A, 'econ');
sigma = diag(D);
lnorm = sum(abs(sigma));
p = sigma./lnorm;
s = p.*log(p);
s(isnan(s)) = 0;
H = -sum(s);
erank = exp(H);

end