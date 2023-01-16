classdef bc < handle

properties
    stepsize        % step size for continuation
    W               % wna object containing everything at equations level
    Wname           % class(W)
    direction       % Table with names (cell) and values (array)
    
    Wsave           % a save of W's across all continuation points
    U               % save of states across all continuation points
    tangent         % array of tangents (currently used)
    tangsave        % array of tangents that have been discontinued
    
    count           % #points since a stepsize change (used in adapt)
    options         % bc options struct
    FDgap = 1e-5    % finite differencing gap
    
    spectrum        % array of spectrum objects (for each point)
    
    to_plot         % array of variables to plot (for each point)
end

methods
    %% constructor
    function obj = bc(direction,stepsize,W,varargin)
        
        % input parser for technical parameters 
        p = inputParser;
        addParameter(p,'parlel',false,@islogical);
        addParameter(p,'tol',1e-10,@(x)validateattributes(x,{'numeric'},...
            {'nonempty','positive'}));
        addParameter(p,'maxiter',15,@(x)validateattributes(x,{'numeric'},...
            {'nonempty','integer','positive'}));
        addParameter(p,'tangentinterpol',true,@islogical);
        addParameter(p,'adapt',true,@islogical);
        addParameter(p,'specmodes',[]);
        parse(p,varargin{:});
        options.tangentinterpol = p.Results.tangentinterpol;
        options.adapt = p.Results.adapt;
        options.maxiter = p.Results.maxiter;
        options.tol = p.Results.tol;
        options.parlel = p.Results.parlel;
        options.numpts = 4;
        options.specmodes = p.Results.specmodes;
        obj.options = options;
        
        % specified things
        direction.values = direction.values/norm(direction.values);
        obj.direction = direction;
        obj.stepsize = stepsize;
        obj.W = W;
        obj.Wname = class(W);
    end
    
    %% set methods
    function set.W(obj,W)
        if ~isequal(obj.W, W)
            % saves previous values
            obj.W = W;
            obj.Wsave = [obj.Wsave, W];
            if size(obj.U,2) == length(obj.Wsave)
                st = W.energies(obj.U(:,end));
                obj.to_plot = [obj.to_plot st];
            end
        end
    end
    
    function set.U(obj,U)
        if ~isequal(obj.U,U)
            obj.U = U;
            if size(obj.U,2) == length(obj.Wsave)
                st = obj.Wsave(end).energies(U(:,end));
                obj.to_plot = [obj.to_plot st];
            end
        end
    end
    
    %% WNA methods
    [phi1,phi20,phi22,c,d,varargout] = WNA_bilinsys(obj,W,varargin)
    DF = FDF(obj,W,vare,varargin)
    
    %% core methods for branch continuation
    varargout = mainloop(obj,U,W,varargin)
    [U,W] = initWNA(obj,initstep,W)
    U = initITER(obj,U,W,varargin)
    [U, W] = branchcont(obj,initU,initW)
    updatetangent(obj,U,W)
    jac = jacAssembly(obj,U,phi,phid,W,tangent)
    stepadapt(obj)
    reset(obj)
    fixtangents(obj,n)
 
    %% plot / postproc methods
    varargout = compareplot(O,varargin)
    [xarr,plotarr] = postproc(obj,varargin)
    plotcompute(obj)
end
end