function varargout = mainloop(obj,U,W,varargin)
% main branch continuation loop

% input parser. if no limit is specified it goes endlessly.
p = inputParser;
addParameter(p,'breakonfold',false,@islogical);
addParameter(p,'limdir','Re',@ischar);
addParameter(p,'limit',Inf);
addParameter(p,'spect',0);
parse(p,varargin{:});
breakonfold = p.Results.breakonfold;
limdir = p.Results.limdir;
limit = p.Results.limit;

% initial iteration
if ~isempty(obj.options.specmodes)
    U = obj.initITER(U,W,p.Results.spect);
else
    U = obj.initITER(U,W);
end

% main loop
limx2 = obj.Wsave(1).(limdir);
while (W.(limdir)-limit)*(limx2-limit) > 0 % break when swaps sign
    [U, W] = obj.branchcont(U,W);
    if isempty(U) % U = [] when NR fails
        obj.reset();
        U = obj.U(:,end);
        W = obj.Wsave(end);
    end
    if breakonfold && (obj.tangent(end,end) < 0)
        break
    end
end

if nargout > 0
    varargout{1} = U;
    varargout{2} = W;
end

end