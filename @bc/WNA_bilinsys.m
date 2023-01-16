function [phi1,phi20,phi22,c,d,varargout] = WNA_bilinsys(obj,W,varargin)
% very 'not-optimized', but it's rarely used anyway.
% varargin is 'Wi','Re', directions for which we want coeffs.

% compute spectrum
if isempty(W.small_spectrum.W)
    W.small_spectrum.do_adjoint = true;
    W.update_spectrum();
end

% coeff to semip matrices
d0cell = repmat({W.D0},1,W.nvar);
d1cell = repmat({W.D1},1,W.nvar);
D0 = blkdiag(d0cell{:});
D1 = blkdiag(d1cell{:});

om_r = 1i*(W.small_spectrum.Lambda(1)); % freq
phi1 = W.small_spectrum.V(:,1); % eigvect
psi = W.small_spectrum.W(:,1); % adjoint

% first eq at order 2
L0 = W.get_LnB(W.LAM,0,0); 
phi1c = conj(phi1); 
RHS = -W.bilin(D0*phi1,D1*phi1,W.k,D0*phi1c,D1*phi1c,-W.k)...
    -W.bilin(D0*phi1c,D1*phi1c,-W.k,D0*phi1,D1*phi1,W.k);
L = W.stopinterf(L0);
L = W.forceBCs(L,'normal');
B = W.stopinterf(RHS);

% eqs for v and p are removed
L(W.N+1:3*W.N,:) = [];
L(:,W.N+1:3*W.N) = [];
B(W.N+1:3*W.N) = [];

% the below searches for a body force that keeps volflux = const
bonuscol = zeros(length(B),1);
bonuscol(2:W.N-1) = ones;
B = [B; 0];
bonuseq = zeros(1,length(B));
WW = W.weights;
bonuseq(1:W.N) = sum((WW*ones(1,W.N)).*W.D0,1);
L = [L, bonuscol; bonuseq];
phi20 = L\B;
phi20 = phi20(1:end-1); % last element was the body force
phi20 = [phi20(1:W.N); zeros(2*W.N,1); phi20(W.N+1:end)];


% second eq at order 2 
Lk = W.get_LnB(W.LAM,2*W.k,2*om_r);
RHS = -W.bilin(D0*phi1,D1*phi1,W.k,D0*phi1,D1*phi1,W.k);
L = W.stopinterf(Lk);
L = W.forceBCs(L,'normal');
B = W.stopinterf(RHS);
phi22 = L\B;

% order 3
RHS = -W.bilin(D0*phi1,D1*phi1,W.k,D0*phi20,D1*phi20,0)...
      -W.bilin(D0*phi20,D1*phi20,0,D0*phi1,D1*phi1,W.k)...
      -W.bilin(D0*phi1c,D1*phi1c,-W.k,D0*phi22,D1*phi22,2*W.k)...
      -W.bilin(D0*phi22,D1*phi22,2*W.k,D0*phi1c,D1*phi1c,-W.k);
RHS = W.stopinterf(RHS);

% derivatives wrt parameters
DL = cell(1,length(varargin)); % varargin is of the form 'Wi','Re',
for j = 1:length(varargin)
    DL{j} = obj.FDF(W,varargin{j},om_r);
    DL{j} = W.stopinterf(DL{j});
end

Dcoeff = 1i*W.B*phi1;
Dcoeff = W.stopinterf(Dcoeff);

% varargout{1}*Wi1 = c|A|^2 + d*om2 (for sample par. Wi)
c = psi'*RHS;
d = psi'*Dcoeff;

% output coeffs of param derivatives
varargout = cell(1,length(varargin));
for j = 1:length(varargin)
    varargout{j} = psi'*DL{j}*phi1;
end
end
