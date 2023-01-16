function stepadapt(obj)
% adapts step
% you should change the exact values (14 and 24) empirically for 
% each problem

numpts = 3; % considers last 3 points
if length(obj.count)>numpts
    if sum(obj.count(end-numpts+1:end)) < 14 % if sum(iterationsteps)<14
        obj.stepsize = 2*obj.stepsize;
        obj.count = []; 
    elseif sum(obj.count(end-numpts+1:end)) > 24
        obj.stepsize = obj.stepsize/2; 
        obj.count = [];
    end
end
end