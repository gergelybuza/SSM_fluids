function varargout = coeffs_to_semip(obj,U)
% maps from coefficient space to semiphysical one
% not sure if worth in general.
if obj.TW
    U = tocomp(U(1:end-2));
else
    U = tocomp(U(1:end-1));
end
D0cell = repmat({sparse(obj.D0)}, 1, obj.nvar*obj.Q);
D0 = blkdiag(D0cell{:});
phi = D0*U;
varargout{1} = phi;
if nargout > 1
    D1cell = repmat({sparse(obj.D1)}, 1, obj.nvar*obj.Q);
    D1 = blkdiag(D1cell{:});
    phid = D1*U;
    varargout{2} = phid;
end
end