classdef Chebgrid < matlab.mixin.SetGetExactNames
    
properties
    %% grid size
    N           % #Chebyshevs
    Q           % #Fouriers
    y           % y grid
    
    %% coeff -> physical space operators
    D0          % 0th order Chebyshev matrix
    D1          % 1st ^^^^^^^^^^^^^^^^^^^^^^
    D2          % 2nd ^^^^^^^^^^^^^^^^^^^^^^
    dif         % Chebyshev differentiation matrix (D1/D0)
    weights     % Chebyshev integration weights
    
    %% other, for large system
    grid_input  % stores grid
    BC          % boundary condition
    nvar        % #vars in system, defined in subclass constructor
    bar         % N*nvar
    dim         % dimension of full system
    dim_chopped % dimension of chopped system
    r = 4       % point at which phase condition is set
    TW          % boolean of travelling wave formulation
    
    %% physical domain (for plotting only)
    x
    y_full
end

methods
    %% constructor
    function obj = Chebgrid(grid_input)
        set(obj,'grid_input',grid_input,'N',grid_input.N,...
            'Q',grid_input.Q,'BC',grid_input.BC)
        if isfield(grid_input,'TW')
            obj.TW = grid_input.TW;
        else
            obj.TW = true;
        end
        switch obj.BC
            case 'full'
                y_st = -1; y_end = 1;
            case 'sym'
                y_st = -1; y_end = 0;
            otherwise
                error(['unknown boundary condition: ' obj.BC])
        end
        
        % create grid
        obj.y = obj.gausslobatto(y_st,y_end);
        [obj.D0, obj.D1, obj.D2,~,~] = obj.cheb(y_st,y_end);
        obj.weights = obj.update_weights(y_st,y_end);
        obj.dif = obj.D1/obj.D0;
    end
    
    %% set methods
    function set.nvar(obj,nvar) 
        obj.nvar = nvar;
        if obj.TW
            obj.dim = obj.nvar*obj.N*obj.Q*2+2;
        else
            obj.dim = obj.nvar*obj.N*obj.Q*2+1;
        end
        obj.dim_chopped = obj.dim-2*obj.N-obj.nvar*obj.N;
    end
    
    %% get methods
    function x = get.x(obj)
        if isempty(obj.x)
            obj.x = 0:(2*pi/obj.k)/300:2*pi/obj.k;
        end
        x = obj.x;
    end
    
    function y_full = get.y_full(obj)
        if isempty(obj.y_full)
            switch obj.BC
                case 'sym'
                    obj.y_full = [-flipud(obj.y); obj.y(2:end)]; 
                case 'full'
                    obj.y_full = obj.y;
            end
        end
        y_full = obj.y_full;
    end
    
    %% grid methods
    y = gausslobatto(obj,y_st,y_end)
    [D0, D1, D2, D3, D4] = cheb(obj,y_st,y_end)
    wk = update_weights(obj,y_st,y_end)
    
    %% other methods
    varargout = coeffs_to_semip(obj,U)
    X = remove_extra_lines(obj,X,varargin)
    [bonuseqc,bonuseqf,bonuscol] = extralines(obj)
end
end






