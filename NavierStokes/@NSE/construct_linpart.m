function A = construct_linpart(obj,U,varargin)
% computes DF(U) = Jacobian of F at U
% U is in (real) coefficient form
% varargin can be {phi, phid} (complex in semiphysical space)
% [phi,phid] = obj.coeffs_to_semip(U);

if length(varargin) > 1
    phi = varargin{1}; 
    phid = varargin{2};
else
    [phi,phid] = obj.coeffs_to_semip(U);
end
if obj.TW
    om_r = U(end);
else
    om_r = 0;
end
bar = obj.bar;

% contributions to Jacobian from the bilinear term
gradphi = constructjac(obj,phi,phid,obj.k); 

% D_{om_r} F (derivative wrt om_r)
gradom = zeros(length(phi),1);
acell = cell(1,obj.Q); % DF(0)
for j = 0:obj.Q-1 % loop over Fourier modes (0 has no om_r term)
    gradom(bar*j+1:bar*(j+1)) =...
        -1i*j*phi(bar*j+1:bar*(j+1));
    % remove the pressure part
    gradom(bar*j+1+2*obj.N:bar*j+3*obj.N) = zeros;
    acell{j+1} = obj.get_LnB(0,j*obj.k,j*om_r); % DF(0) terms
end
gradom = toreal(gradom);
gradphi = toreal(blkdiag(acell{:})) + gradphi; % linear + bilin parts

% assemble full gradient
[bonuseqc,bonuseqf,bonuscol] = obj.extralines();
if obj.TW
    A = [gradphi, bonuscol, gradom; bonuseqf; bonuseqc];
else
    A = [gradphi, bonuscol; bonuseqf];
end

% enforce BCs + reduce to chopped system
A = obj.stopinterf(A);
A = obj.forceBCs(A,'bc');
A = obj.remove_extra_lines(A); % removes trivial lines
end

function NLPART = constructjac(O,phi,phid,k)
% constructs bilinear part of the jacobian from small building blocks
% corresponding to each Fourier mode
M = length(phi)/O.nvar/O.N-1;
% if O.parlel
%     NLPART1 = cell(1,(M+1));
%     NLPART2 = cell(1,(M+1));
%     parfor q = 0:M 
%         [NLPART1{q+1}, NLPART2{q+1}] = kloop(O,phi,phid,k,q);
%         % if this is inefficient, try codistributing wrt rows
% %         NLPART1{q+1} = distributed(NLPART1{q+1}); 
% %         NLPART2{q+1} = distributed(NLPART2{q+1});
%     end
%     NLPART = [NLPART1, NLPART2];
% else
    NLPART = cell(1,2*(M+1));
    for q = 0:M 
        [NLPART{q+1}, NLPART{M+1+q+1}] = kloop(O,phi,phid,k,q);
    end
% end
NLPART = cat(1,NLPART{:});
end

function block = blockone(O,phi,phid,k,varargin)
% D Bilinpart wrt first argument
% you need the complementary phi's in here.
% k is n*k for n the complementary index
M = ones(1,O.N);
Z = zeros(O.N);
[u,v] = phi_parsern(O,phi);
[du,dv] = phi_parsern(O,phid);
block = ...
[ (1i*k*u*M).*O.D0 (du*M).*O.D0 Z;
(1i*k*v*M).*O.D0 (dv*M).*O.D0 Z;
Z Z Z];
end

function block = blocktwo(O,phi,k,varargin)
% D Bilin part wrt second argument
% you need the complementary phi's in here. (so first entries)
% k is n*k for n the primary index (wrt which I differentiated)
M = ones(1,O.N);
Z = zeros(O.N);
[u,v] = phi_parsern(O,phi);
block = ...
[ (1i*k*u*M).*O.D0+(v*M).*O.D1 Z Z;
Z (1i*k*u*M).*O.D0+(v*M).*O.D1 Z;
Z Z Z];
end

function [u,v] = phi_parsern(O,phi)
u = phi(1:O.N);
v = phi(O.N+1:2*O.N);
end

function [NLPARTr, NLPARTi] = kloop(O,phi,phid,k,q)

bar = O.N*O.nvar;
M = length(phi)/bar-1;

% these are horizontal bars correspondining to each 'k'
NLPARTr = zeros(bar,2*length(phi)); 
NLPARTi = zeros(bar,2*length(phi));
for index = -M+q:M
    cindex = q-index;
    if cindex < 0
        entry2 = conj(phi(-cindex*O.nvar*O.N+1:(-cindex+1)*O.nvar*O.N));
        entry2d = conj(phid(-cindex*O.nvar*O.N+1:(-cindex+1)*O.nvar*O.N));
    else
        entry2 = phi(cindex*O.nvar*O.N+1:(cindex+1)*O.nvar*O.N);
        entry2d = phid(cindex*O.nvar*O.N+1:(cindex+1)*O.nvar*O.N);
    end
    bone = blockone(O,entry2,entry2d,cindex*k);
    btwo = blocktwo(O,entry2,index*k);
    if index < 0
    NLPARTr(:,-index*bar+1:(-index+1)*bar) =...
        NLPARTr(:,-index*bar+1:(-index+1)*bar)+...
        real(bone)+real(btwo);
    NLPARTr(:,(M+1-index)*bar+1:(M+1-index+1)*bar) =...
        NLPARTr(:,(M+1-index)*bar+1:(M+1-index+1)*bar)...
        +imag(bone)+imag(btwo);
    NLPARTi(:,(-index)*bar+1:(-index+1)*bar) =...
        NLPARTi(:,(-index)*bar+1:(-index+1)*bar)...
        +imag(bone)+imag(btwo);
    NLPARTi(:,(M+1-index)*bar+1:(M+1-index+1)*bar) =...
        NLPARTi(:,(M+1-index)*bar+1:(M+1-index+1)*bar)...
        -real(bone)-real(btwo);
    else
    NLPARTr(:,index*bar+1:(index+1)*bar) =...
        NLPARTr(:,index*bar+1:(index+1)*bar)+...
        real(bone)+real(btwo);
    NLPARTr(:,(M+1+index)*bar+1:(M+1+index+1)*bar) =...
        NLPARTr(:,(M+1+index)*bar+1:(M+1+index+1)*bar)...
        -imag(bone)-imag(btwo);
    NLPARTi(:,(index)*bar+1:(index+1)*bar) =...
        NLPARTi(:,(index)*bar+1:(index+1)*bar)+...
        imag(bone)+imag(btwo);
    NLPARTi(:,(M+1+index)*bar+1:(M+1+index+1)*bar) =...
        NLPARTi(:,(M+1+index)*bar+1:(M+1+index+1)*bar)+...
        real(bone)+real(btwo);
    end
end
end

