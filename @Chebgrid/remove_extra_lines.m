function X = remove_extra_lines(obj,X,varargin)
% removes trivial lines from the equations

if isempty(varargin)
    mode = 'single';
else
    mode = varargin{1};
end
relevant = [1:obj.N 3*obj.N+1:obj.Q*obj.bar (obj.Q+1)*obj.bar+1:obj.dim];
ndimsX = numel(find(size(X) == obj.dim)); 
relevant = repmat({relevant},1,ndimsX);
switch mode
    case 'single'
        X = X(relevant{:});
    case 'more'
        X = X(relevant{:},:);
end
end