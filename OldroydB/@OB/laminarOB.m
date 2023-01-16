function [U,Txx,Txy,Fx] = laminarOB(obj)
% laminar solver (needed only due to the polymeric diffusion term)

switch obj.BC
    case 'sym'
        Qexp = 1;
    case 'full'
        Qexp = 2;
end
Fx = 3;

% initial data
Uin = 1.5*(1-obj.y.^2);
Up = -3*obj.y;
Txxin = 2*obj.Wi*(1-obj.beta)*Up.^2;
Txyin = (1-obj.beta)*Up;
x0 = [Uin; Txxin; Txyin; Fx];
err = 1;

% the below is to ensure bulk velocity = 1
bonuseq2 = zeros(1,length(x0));
WW = obj.weights(1:obj.N);
bonuseq2(1:obj.N) = sum((WW*ones(1,obj.N)).*obj.D0,1);
bonuscol = zeros(length(x0)-1,1);
bonuscol(2:obj.N-1) = -ones;

D0 = blkdiag(obj.D0,obj.D0,obj.D0);
while err > 1e-10 % Newton-Raphson
    U = x0(1:obj.N);
    U0 = obj.D0\U;
    Txx = x0(obj.N+1:2*obj.N);
    Txy = x0(2*obj.N+1:3*obj.N);
    Fx = x0(end);
    
    jac = DF(obj,U,Txx,Txy);
    jac2 = [jac, bonuscol; bonuseq2];
    F = eval(obj,U,Txx,Txy,Fx);
    F2 = [F; bonuseq2(1:obj.N)*U0-Qexp]; 
    
    mu = jac2\F2;
    xx = x0(1:end-1) - D0*mu(1:end-1);
    xx2 = [xx; x0(end)-mu(end)];
    err = norm(xx2-x0)/norm(x0);
%      disp(err)
    x0 = xx2;
end

U = x0(1:obj.N);
Txx = x0(obj.N+1:2*obj.N);
Txy = x0(2*obj.N+1:3*obj.N);
Fx = x0(end);

end



function jac = DF(obj,U,Txx,Txy)  
% computes jacobian (of laminar eqs. F = 0)

% preliminary stuff
Up = obj.dif*U;
M = ones(1,length(U));
Z = zeros(length(U));

% Polymer diffusion term
DDD = obj.Sctrickf(-obj.epsilon*obj.D2);

jac = [-obj.beta*obj.D2, Z, -obj.D1;
    -2*(Txy*M).*obj.D1, obj.D0/obj.Wi+DDD, (-2*(Up)*M).*obj.D0;
    -(1-obj.beta)/obj.Wi*obj.D1, Z, obj.D0/obj.Wi+DDD];

% boundary conditions
jac(1,:) = 0;
switch obj.BC % U_top
    case 'sym'
        jac(1,1:obj.N) = obj.D1(1,:);
    case 'full'
        jac(1,1:obj.N) = obj.D0(1,:);
end
jac(obj.N,:) = 0;
jac(obj.N,1:obj.N) = obj.D0(end,:); % U_bot

% symmetry conditions for T whenever necessary
if isequal(obj.BC,'sym') && (obj.epsilon ~= 0) 
    jac(obj.N+1,:) = 0;
    jac(obj.N+1,obj.N+1:2*obj.N) = obj.D1(1,:); % Txx
    jac(2*obj.N+1,:) = 0;
    jac(2*obj.N+1,2*obj.N+1:3*obj.N) = obj.D0(1,:); % Txy
end

end

function F = eval(obj,U,Txx,Txy,Fx)
% evaluates the laminar eqs F

% preliminary stuff
Txxp = obj.dif*Txx;
Txyp = obj.dif*Txy;
Up = obj.dif*U;
Upp = obj.D2*(obj.D0\U);

% Polymer diffusion term
DDD = obj.Sctrickf(-obj.epsilon*obj.D2);

F = [-Fx-obj.beta*Upp-Txyp;
    1/obj.Wi*Txx-2*Txy.*Up+DDD*(obj.D0\Txx);
    1/obj.Wi*Txy-(1-obj.beta)/obj.Wi*Up+DDD*(obj.D0\Txy)];

% boundary conditions
switch obj.BC 
    case 'sym'
        F(1) = Up(1); % U_top
    case 'full'
        F(1) = U(1);
end
F(obj.N) = U(end); % U_bottom

% symmetry conditions for T whenever necessary
if isequal(obj.BC,'sym') && (obj.epsilon ~= 0) 
    F(obj.N+1) = Txxp(1);
    F(2*obj.N+1) = Txy(1);
end

end
