function [U,W2] = initWNA(obj,initstep,W)
% gets an initial condition for the first point based on WNA

% Weakly nonlinear analysis
coefarr = cell(1,length(obj.direction.names));
coefarr2 = zeros(2,length(obj.direction.names));
[eigv,~,~,c,d,coefarr{:}] = obj.WNA_bilinsys(W,obj.direction.names{:});

if abs(real(W.small_spectrum.Lambda(1))) > 1e-5
    disp('Might not be close enough to the bifurcation.')
end

% initiate spectrum array (if requested)
if ~isempty(obj.options.specmodes)
    spec = spectrum();
    spec.Lambda = real(W.small_spectrum.Lambda(obj.options.specmodes));
%     spec.Lambda = 0;
    obj.spectrum = spec;
end

om_r = real(1i*(W.small_spectrum.Lambda(1)));

% putting init data in for bifurcation point.
Fx = W.Fx; 
obj.U = [W.LAM; zeros(W.dim-length(W.LAM)-2,1); Fx; om_r];

% solve WNA equations c|A|^2 + d \omega = coefarr{i}*par_1{i}
for j = 1:length(coefarr) 
coefarr2(:,j) =...
[real(c), real(d); imag(c), imag(d)]\[real(coefarr{j}); imag(coefarr{j})];
end

eps = 0.1;
initstep = initstep/eps^2; % (pars scale with eps^2)

% solution
amp = real(sqrt(initstep*coefarr2(1,:)*obj.direction.values)); 
om1 = initstep*coefarr2(2,:)*obj.direction.values;
om_r = om_r + eps^2*om1;

% generate initial condition
U = zeros(W.nvar*W.N*W.Q,1);
U(1:W.bar) = W.LAM;
z = W.D0(W.r,:)*eigv(1:W.N); % twist IC to match phase of BC
U(W.nvar*W.N+1:2*W.nvar*W.N) = eps*amp*eigv*conj(z/abs(z)); 
U = toreal(U);
U = [U; Fx; om_r];

% adjust pars to match those where IC is
pars = W.pars;
for j = 1:length(obj.direction.names)
    pars.(obj.direction.names{j}) = pars.(obj.direction.names{j}) + ...
        eps^2*initstep*obj.direction.values(j);
end

% save W
W2 = feval(obj.Wname,W.grid_input,pars);
obj.W = W2;
end