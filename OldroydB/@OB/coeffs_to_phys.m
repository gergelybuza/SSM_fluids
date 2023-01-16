function [u_1,u_2,t_11,t_22,t_12,du_1,du_2,Fx,c,varargout] = coeffs_to_phys(obj,U)
% maps the coefficient array 'U' into physical space

c = U(end)/obj.k;       % phase speed
Fx = U(end-1);          % forcing

% velocities, stresses and their derivatives
[U,DU] = obj.coeffs_to_semip(U);
U = reshape(U,obj.bar,[]);
DU = reshape(DU,obj.bar,[]);
u = U(:,1)*ones(size(obj.x));
du = DU(:,1)*ones(size(obj.x));
for j = 1:size(U,2)-1
    u = u + U(:,j+1)*exp(1i*j*obj.x*obj.k) + conj(U(:,j+1))*exp(-1i*j*obj.x*obj.k);
    du = du + U(:,j+1)*exp(1i*j*obj.x*obj.k) + conj(U(:,j+1))*exp(-1i*j*obj.x*obj.k);
end

switch obj.BC
    case 'sym'
        u_1 = [flipud(u(1:obj.N,:)); u(2:obj.N,:)];
        u_2 = [-flipud(u(obj.N+1:2*obj.N,:)); u(obj.N+2:2*obj.N,:)];
        t_11 = [flipud(u(3*obj.N+1:4*obj.N,:)); u(3*obj.N+2:4*obj.N,:)];
        t_22 = [flipud(u(4*obj.N+1:5*obj.N,:)); u(4*obj.N+2:5*obj.N,:)];
        t_12 = [-flipud(u(5*obj.N+1:6*obj.N,:)); u(5*obj.N+2:6*obj.N,:)];
        du_1 = [-flipud(du(1:obj.N,:)); du(2:obj.N,:)];
        du_2 = [flipud(du(obj.N+1:2*obj.N,:)); du(obj.N+2:2*obj.N,:)];
    case 'full'
        u_1 = u(1:obj.N,:);
        u_2 = u(obj.N+1:2*obj.N,:);
        du_1 = du(1:obj.N,:);
        du_2 = du(obj.N+1:2*obj.N,:);
        t_11 = u(3*obj.N+1:4*obj.N,:);
        t_22 = u(4*obj.N+1:5*obj.N,:);
        t_12 = u(5*obj.N+1:6*obj.N,:);
end

% compute stream function and vorticity, if requested.
if nargout > 9
    U_2 = U(obj.N+1:2*obj.N,:);
    DU_1 = DU(1:obj.N,:);
    if strcmp(obj.BC,'sym')
        DU_1 = [-flipud(DU_1); DU_1];
        U_2 = [-flipud(U_2); U_2];
    end
    vorticity = -DU_1(:,1)*ones(size(obj.x));
    stream_func = zeros(size(vorticity));
for j = 1:size(U,2)-1
    stream_func = stream_func ...
        - U_2(:,j+1)/(1i*j*obj.k)*exp(1i*j*obj.k*obj.x) ...
        - conj(U_2(:,j+1))/(-1i*j*obj.k)*exp(-1i*j*obj.k*obj.x);
    vorticity = vorticity ...
        + U_2(:,j+1)*(1i*j*obj.k)*exp(1i*j*obj.k*obj.x)...
        - DU_1(:,j+1)*exp(1i*j*obj.k*obj.x)...
        + conj(U_2(:,j+1))*(-1i*j*obj.k)*exp(-1i*j*obj.k*obj.x)...
        - conj(DU_1(:,j+1))*exp(-1i*j*obj.k*obj.x);
end
varargout{1} = real(stream_func);
varargout{2} = real(vorticity);
end

end