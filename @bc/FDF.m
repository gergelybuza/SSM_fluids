function  DF = FDF(obj,W,vare,varargin)
% finite differences F or L, depending on what varargin you call it with
% vare contains the parameter (in 'pars') wrt which you want to diff.
% varargin should either contain arguments of L (om_r) for WNA
% or the arguments of F for BC

pars = W.pars;
grid_input = W.grid_input;
pars.(vare) = pars.(vare) + obj.FDgap;  
W_plus = feval(obj.Wname,grid_input,pars);
pars.(vare) = pars.(vare) - 2*obj.FDgap;
W_minus = feval(obj.Wname,grid_input,pars);
if length(varargin) == 1 % this is for WNA, FDs L
    DF = (W_plus.get_LnB(W_plus.LAM,W_plus.k,varargin{1})-...
        W_minus.get_LnB(W_minus.LAM,W_minus.k,varargin{1}))/2/obj.FDgap;
else % this is for BC
    DF = (W_plus.F(varargin{:})-W_minus.F(varargin{:}))/2/obj.FDgap;
end

end