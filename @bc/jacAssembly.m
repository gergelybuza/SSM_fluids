function jac = jacAssembly(obj,U,phi,phid,W,tangent)
% jacobian assembler

%% assemble gradient wrt parameter(s)
% find relevant directions
nonzero = find(obj.direction.values);
gradarr = zeros(W.dim_chopped,length(nonzero));

% compute finite differenced gradients
for j = 1:length(nonzero)
    gradarr(:,j) = obj.FDF(W,obj.direction.names{nonzero(j)},U,phi,phid);
end
gradp = gradarr*obj.direction.values(nonzero);

%% full jacobian
jac = [W.construct_linpart(U,phi,phid), gradp; tangent];

end