function V = bilin(obj,phi1,phi1d,k1,phi2,phi2d,k2)
% bilinear term in NSE

[u1,v1] = phi_parser2(obj,phi1);
[u2,v2] = phi_parser2(obj,phi2);
[u2d,v2d] = phi_parser2(obj,phi2d);

V = [ (1i*k2*u1.*u2+v1.*u2d);
      (1i*k2*u1.*v2+v1.*v2d);
      zeros(obj.N,1)];

end

function [u,v] = phi_parser2(O,phi)
u = phi(1:O.N);
v = phi(O.N+1:2*O.N);
end