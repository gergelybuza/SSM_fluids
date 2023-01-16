function [bonuseqc,bonuseqf,bonuscol] = extralines(obj)
% extra lines for Jacobian

bar = obj.bar;

% set bonus equation for phase degeneracy
bonuseqc = zeros(1,obj.dim);
bonuseqc(bar*(obj.Q+1)+1:bar*obj.Q+(obj.nvar+1)*obj.N) =...
    obj.D0(obj.r,:); 

% set bonus equation for body force (setting bulk vel = 0)
bonuseqf = zeros(1,obj.dim);
bonuseqf(1:obj.N) = sum((obj.weights*ones(1,obj.N)).*obj.D0,1);
bonuscol = zeros(obj.nvar*obj.N*obj.Q*2,1);
bonuscol(2:obj.N-1) = -ones; % DF of body force term
end