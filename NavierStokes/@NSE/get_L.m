function L = get_L(obj,LAM,k)
% Linear operator.

% extract laminar state
if LAM == 0
    U = zeros(obj.N,1);
    Up = zeros(obj.N,1);
else
    U = obj.D0*LAM(1:obj.N);
    Up = obj.D1*LAM(1:obj.N);
end
M  = ones(1,length(U)); 
Z  = zeros(obj.N);

% m'operator
L = ...
 ... %----u------------------------------------v--------------------------------------------p-----------
[  (k^2*obj.D0-obj.D2)/obj.Re+1i*k*(U*M).*obj.D0 (Up*M).*obj.D0                        1i*k*obj.D0  ;  % x momentum
   Z                                           (k^2*obj.D0-obj.D2)/obj.Re+1i*k*(U*M).*obj.D0  obj.D1      ;  % y momentum
   1i*k*obj.D0                                 obj.D1                                       Z           ];  % continuity]
    
end

