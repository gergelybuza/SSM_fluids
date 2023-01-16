function [lambda,v] = inverse_iter(obj,A,M,mu,varargin)
% inverse iteration
% mu is the approximate eigenvalue
% varargin{1} is approximate eigenvector (if specified)

% initial eigenvector guess b_0
if ~isempty(varargin)
    b_0 = varargin{1};
else
    b_0 = ones(length(M),1);
end

err = 1;
tol = 1e-12;
C = (A-mu*M)\M; % predefine inverted matrix
while err > tol
    b = C*b_0; % iterate
    th = b_0'*b;
    magb_0 = b_0'*b_0;
    lam = magb_0/th+mu; % approx eigval
    err = norm(b-th*b_0)/abs(th); % error
    b_0 = b/norm(b); % next eigenvector
end
v = b_0;
lambda = lam;
end