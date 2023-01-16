function spectrum = normalize_modes(obj,spectrum,M)
% Copied from SSMTool.
%NORMALIZE_MODES: This function normalizes the right and left eigenvectors  
% W, V withrespect to the matrix M. 

V = spectrum{1};
V = V*diag(1./vecnorm(V));
spectrum{1} = V;

if length(spectrum) > 2
    W = spectrum{3};
    mu = diag(W'*M*V);
    W = W*diag(1./(mu'));
    spectrum{3} = W;
end
end