function fixtangents(obj,n)
% fixes the tangent n+1 points prior (so 0 => prev point)
% this is lazy coding but tangent(end) is one point behind here
ddist = 0;
vars = obj.direction.names;
for j = 1:length(vars) % full distance between points
    ddist = ddist + obj.direction.values(j)*...
        (obj.Wsave(end-n).(vars{j}) - obj.Wsave(end-n-1).(vars{j}));
end
obj.tangent(:,end-n) = ...
[obj.Wsave(end).remove_extra_lines(obj.U(:,end-n)-obj.U(:,end-n-1),'more'); ddist];
obj.tangent(:,end-n) = obj.tangent(:,end-n)/norm(obj.tangent(:,end-n));
end