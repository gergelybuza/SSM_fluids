function F = F(obj,U,varargin)
% computes F(U)
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
    Fx = U(end-1);
    om_r = U(end);
    U = tocomp(U(1:end-2));
else
    Fx = U(end);
    om_r = 0;
    U = tocomp(U(1:end-1));
end
k = obj.k;
bar = obj.bar;

F = zeros(size(phi)); 
for q = 0:obj.Q-1 % loop over all +ve Fourier modes
% Linear part
summation = obj.get_LnB(0,q*k,q*om_r)*U(q*bar+1:(q+1)*bar);
% Bilinear part
for index = 0:obj.Q-1
    cindex = q-index;
    if cindex < 0
        % conjugate, do bilin two ways! (sym)
        entry1 = phi(index*bar+1:(index+1)*bar);
        entry1d = phid(index*bar+1:(index+1)*bar);
        entry2 = conj(phi(-cindex*bar+1:(-cindex+1)*bar));
        entry2d = conj(phid(-cindex*bar+1:(-cindex+1)*bar));
        summation = summation ...
            + obj.bilin(entry1,entry1d,index*k,entry2,entry2d,cindex*k)+...
            obj.bilin(entry2,entry2d,cindex*k,entry1,entry1d,index*k);
    else
        % don't conj, do bilin but only one way.
        % (these will repeat naturally once the conjugate index is reached)
        entry1 = phi(index*bar+1:(index+1)*bar);
        entry1d = phid(index*bar+1:(index+1)*bar);
        entry2 = phi(cindex*bar+1:(cindex+1)*bar);
        entry2d = phid(cindex*bar+1:(cindex+1)*bar);
        summation = summation ...
            + obj.bilin(entry1,entry1d,index*k,entry2,entry2d,cindex*k);
    end      
end
F(q*bar+1:(q+1)*bar) = summation;
% BCs
switch obj.BC
    case 'sym'
        F(q*bar+1) = phid(q*bar+1);
    case 'full'
        F(q*bar+1) = phi(q*bar+1);
end
bcspots = [q*bar+obj.N q*bar+obj.N+1 q*bar+2*obj.N];
F(bcspots) = phi(bcspots);
end
F = toreal(F);

% extras
F(2:obj.N-1) = F(2:obj.N-1) - Fx;
FFx = dot(obj.weights,real(phi(1:obj.N))) - dot(obj.weights,obj.U);
Fomr = imag(phi(bar+obj.r));
F = [F; FFx; Fomr];
F = obj.remove_extra_lines(F);
end





