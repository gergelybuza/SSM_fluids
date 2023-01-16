function X = stopinterf(obj,X)
% removes lines in preparation of BCs;
% also used for setting RHS zero in same lines.

bar = obj.bar;

if size(X,1) < obj.dim-1
    X(1:bar:end,:,:) = zeros;             % utop
else
    X(1:bar:obj.dim-2,:,:) = zeros;             % utop
end
X(obj.N:bar:end,:,:) = zeros;           % ubot
X(obj.N+1:bar:end,:,:) = zeros;         % vtop
X(2*obj.N:bar:end,:,:) = zeros;         % vbot

if isequal(obj.BC,'sym') && (obj.epsilon ~= 0) 
% these are symmetry conditions for T
X(3*obj.N+1:bar:end,:,:) = zeros; % txxtop
X(4*obj.N+1:bar:end,:,:) = zeros; % tyytop
X(5*obj.N+1:bar:end,:,:) = zeros; % txytop
end

end