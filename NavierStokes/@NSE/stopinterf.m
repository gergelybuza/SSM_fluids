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

end