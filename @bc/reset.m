function reset(obj)
% reset function for when NR fails

if size(obj.tangent,2) > 1
    obj.tangsave = [obj.tangsave obj.tangent(:,1:end-2)];
    obj.tangent = obj.tangent(:,end-1); % this makes sure that the tangent
    % is no longer interpolated (in case that was the issue)
else
    obj.tangent = [];
end
obj.stepsize = obj.stepsize/2; % halves stepsize
obj.count = []; % for consistent 'stepadapt's
end