function A = toreal(X)
% turns complex input vector or matrix into a real- version.
if size(X,2)>1
    A = [real(X), -imag(X); imag(X), real(X)];
else
    A = [real(X); imag(X)];
end
end