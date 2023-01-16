function V = bilin(obj,phi1,phi1d,k1,phi2,phi2d,k2)
% bilinear term in OB

[u1,v1,txx1,tyy1,txy1] = phi_parser2(obj,phi1);
[u2,v2,txx2,tyy2,txy2] = phi_parser2(obj,phi2);
[u2d,v2d,txx2d,tyy2d,txy2d] = phi_parser2(obj,phi2d);

V = [ obj.Re*(1i*k2*u1.*u2+v1.*u2d);
      obj.Re*(1i*k2*u1.*v2+v1.*v2d);
      zeros(obj.N,1);
      (1i*k2*u1.*txx2+v1.*txx2d-2*1i*k2*txx1.*u2-2*txy1.*u2d);
      (1i*k2*u1.*tyy2+v1.*tyy2d-2*1i*k2*txy1.*v2-2*tyy1.*v2d);
      (1i*k2*u1.*txy2+v1.*txy2d-1i*k2*txx1.*v2-1i*k2*txy1.*u2-txy1.*v2d-tyy1.*u2d)];
end

function [u,v,txx,tyy,txy] = phi_parser2(obj,phi)
u = phi(1:obj.N);
v = phi(obj.N+1:2*obj.N);
txx = phi(3*obj.N+1:4*obj.N);
tyy = phi(4*obj.N+1:5*obj.N);
txy = phi(5*obj.N+1:end);
end