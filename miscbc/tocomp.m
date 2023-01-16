function A = tocomp(X)
% turns an input vector of the form (ReA,ImA) into A.
A = X(1:length(X)/2)+1i*X(length(X)/2+1:end);
end