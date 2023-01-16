function M = construct_mass(obj)
% 'mass' matrix, coefficient of du/dt
if obj.TW
    bcell = [repmat({(obj.B)},1,2*obj.Q) zeros(2)];
else
    bcell = [repmat({(obj.B)},1,2*obj.Q) zeros(1)];
end
M = blkdiag(bcell{:});
M = obj.stopinterf(M);
M = obj.remove_extra_lines(M);
end