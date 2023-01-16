function [plot_struct,varargout] = Wztoplot_1D(z,W_0,W,U_b,varargin)
% computes a 'plot' array across the specified array 'z'
% for 1D (real or complex) reduced models only

if ~isempty(varargin) % specify whatever to reach 'polar' mode
    Wfun = @(x) reduced_to_full([x; conj(x)],W_0,[],0);
else
    Wfun = @(x) reduced_to_full(x,W_0,[],0);
end
plotfun = @(s) toplot(s,Wfun,W,U_b);
if nargout > 1
[plot_struct,varargout{1}] = arrayfun(plotfun,z);
else
plot_struct = arrayfun(plotfun,z);
end
end

function [plot,varargout] = toplot(z,Wfun,W,U_b)
U = U_b + readd(W,Wfun(z));
plot = W.energies(U);
if nargout > 1
    varargout{1} = norm(W.F(U));
end
end