classdef DynSysOB < OB
    
properties
    A                       % main linear operator
    M                       % 'mass' matrix in fromt of du/dt
    F                       % cell of polynomial orders of RHS: F = {A B}
    
    spectrum                % relevant part of spectrum(A,M)
    Options = DSOptions()   % Options for mfd computation from SSMTool
end
methods
    %% constructor
    function obj = DynSysOB(grid_input,pars,U,notation,varargin)
        % varargin is approximate spectrum, if known
        obj = obj@OB(grid_input,pars); % call OB constructor
        
        obj.A = -obj.construct_linpart(U);
        obj.M = obj.construct_mass();
        if isempty(varargin) % compute full spectrum
            obj.spectrum = spectrum(obj.A,obj.M);
        else % only refine the specified part in (varargin)
            obj.spectrum = spectrum(obj.A,obj.M,varargin{1});
        end
        
        obj.Options.notation = notation;
        % assemble F
        switch notation
            case 'tensor'
                obj.F = {(obj.A), -obj.nonlinearity()};
            case 'multiindex'
%                 obj.F = {tensor_to_multi_index(sptensor(obj.A)),...
%                     tensor_to_multi_index(-obj.nonlinearity())};
                obj.F = {tensor_to_multi_index(sptensor(obj.A)),...
                    obj.nonlinearity_multiindex()};
            otherwise
                error('notation should be tensor or multiindex')
        end
    end
    
    %% methods assembling the bilinear form 
    N = nonlinearity(obj)
    N = nonlinearity_multiindex(obj)
end    
end