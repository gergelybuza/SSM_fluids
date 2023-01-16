function X = forceBCs(obj,X,type)
% enforces boundary conditions

bar = obj.bar;

switch type
    case 'bc'
        M = obj.Q-1;
        er = 1;
    case 'trick'
        M = -0.5;
        er = -200000*1i; % gets rid of spurious bdy modes for linstab
    case 'normal'
        M = -0.5;
        er = 1;
end

for j = 0:2*M+1 % goes through all Fmodes in real/imag

% setting u at top bdy:
switch obj.BC
    case 'sym'
        X(j*bar+1,j*bar+1:j*bar+obj.N) = er*obj.D1(1,:);
    case 'full'
        X(j*bar+1,j*bar+1:j*bar+obj.N) = er*obj.D0(1,:);
end

X(j*bar+obj.N,j*bar+1:j*bar+obj.N) = er*obj.D0(end,:); % u at bottom bdy
X(j*bar+obj.N+1,j*bar+obj.N+1:j*bar+2*obj.N) = er*obj.D0(1,:); % v top
X(j*bar+2*obj.N,j*bar+obj.N+1:j*bar+2*obj.N) = er*obj.D0(end,:); % v bot

end

end