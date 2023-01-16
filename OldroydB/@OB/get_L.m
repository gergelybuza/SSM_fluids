function L = get_L(obj,LAM,k)
% Linear operator.

% extract laminar state
if LAM == 0
    U = zeros(obj.N,1);
    Txx = zeros(obj.N,1);
    Txy = zeros(obj.N,1);
    Up = zeros(obj.N,1);
    Txxp = zeros(obj.N,1);
    Txyp = zeros(obj.N,1);
else
    U = obj.D0*LAM(1:obj.N);
    Txx = obj.D0*LAM(3*obj.N+1:4*obj.N);
    Txy = obj.D0*LAM(5*obj.N+1:end);
    Up = obj.D1*LAM(1:obj.N);
    Txxp = obj.D1*LAM(3*obj.N+1:4*obj.N);
    Txyp = obj.D1*LAM(5*obj.N+1:end);
end

M  = ones(1,length(U)); 
Z  = zeros(obj.N);

% polymeric dissipation term (Sctrick removes @ boundary)
DDD = obj.Sctrickf(-obj.epsilon*obj.D2+obj.epsilon*(k^2)*obj.D0);

% m'operator
L = ...
 ... %----u------------------------------v-----------------------------------------p----------------c11-----------------------c22------------------------c12---------------------%
[ -obj.beta*obj.D2+obj.Re*1i*k*(U*M).*obj.D0+obj.beta*k^2*obj.D0                    obj.Re*(Up*M).*obj.D0                                                   1i*k*obj.D0      -1i*k*obj.D0                               Z                                       -obj.D1                                 ;  % x momentum
   Z                                                                                -obj.beta*obj.D2+obj.Re*1i*k*(U*M).*obj.D0+obj.beta*k^2*obj.D0          obj.D1           Z                                          -obj.D1                                 -1i*k*obj.D0                            ;  % y momentum
   1i*k*obj.D0                                                                      obj.D1                                                                  Z                Z                                          Z                                       Z                                       ;  % continuity
   -2*(Txy*M).*obj.D1-2*1i*k*(Txx*M).*obj.D0-2*(1-obj.beta)/obj.Wi*1i*k*obj.D0      (Txxp*M).*obj.D0                                                        Z                obj.D0/obj.Wi+1i*k*(U*M).*obj.D0+DDD       Z                                       -2*(Up*M).*obj.D0                       ;  % c11 evolution
   Z                                                                                -2*(1-obj.beta)/obj.Wi*obj.D1-2*1i*k*(Txy*M).*obj.D0                    Z                Z                                          obj.D0/obj.Wi+1i*k*(U*M).*obj.D0+DDD    Z                                       ;  % c22 evolution
   -(1-obj.beta)/obj.Wi*obj.D1                                                      (Txyp*M).*obj.D0-1i*k*(Txx*M).*obj.D0-(1-obj.beta)/obj.Wi*1i*k*obj.D0   Z                Z                                          -(Up*M).*obj.D0                         obj.D0/obj.Wi+1i*k*(U*M).*obj.D0+DDD    ];  % c12 evolution    
end

