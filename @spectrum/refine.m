function [V,Lambda,W] = refine(obj,A,M,tangentModes)
% refines eigenvalues/eigenvectors of the relevant modes
% via inverse iteration (will be more accurate than eig)

mus = obj.Lambda(tangentModes);
Vs = obj.V(:,tangentModes);
Ws = obj.W(:,tangentModes);

[V,Lambda,W] = obj.inv_iter_all(A,M,mus,Vs,Ws);
end