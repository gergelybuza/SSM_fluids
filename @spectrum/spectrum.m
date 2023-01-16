classdef spectrum < handle
    
properties
    Lambda                  % eigenvalues
    V                       % eigenvectors
    W                       % adjoint eigenvectors
    do_adjoint = true       % whether or not to compute W
end

methods
    %% constructor
    function obj = spectrum(varargin)
        % varargin{1} = A -- linearized operator
        % varargin{2} = M -- 'mass' matrix
        % varargin{3} = spectrum -- only specify this if invers_iter
        if length(varargin) == 1 
            obj.compute_full_spectrum(varargin{1})
        elseif length(varargin) == 2 
            obj.compute_full_spectrum(varargin{1},varargin{2})
        elseif length(varargin) > 2 
            spec = varargin{3};
            mus = spec.Lambda;
            if ~isempty(spec.W)
            Vs = spec.V;
            Ws = spec.W;
            [obj.V,obj.Lambda,obj.W] =...
                obj.inv_iter_all(varargin{1},varargin{2},mus,Vs,Ws);
            else
            [obj.V,obj.Lambda,obj.W] =...
                obj.inv_iter_all(varargin{1},varargin{2},mus);
            end
        end
    end
    
    %% spectrum methods
    compute_full_spectrum(obj,A,varargin)
    [lambda,v] = inverse_iter(obj,A,M,mu,varargin)
    [V,Lambda,W] = refine(obj,A,M,tangentModes)
    spectrum = normalize_modes(obj,spectrum,M)
    [V,Lambda,W] = inv_iter_all(obj,A,M,mus,varargin)
end
end