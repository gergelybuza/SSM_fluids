classdef NSE < Chebgrid
    
properties  
    %% flow parameters
    pars            % struct of parameters (Re,k)
    Re              % Reynolds number
    k               % (streamwise) wavenumber
    
    %% laminar state
    U               % base (streamwise) velocity (in phys)
    LAM             % full laminar state (coeffs)
    Fx              % force sustaining the base state

    %% other 
    B               % diagonal coeff matrix (for RHS of linstab)
    small_spectrum  % spectrum object (for the 1-mode expansion)
end

methods
    %% constructor
    function obj = NSE(grid_input,pars)
        obj = obj@Chebgrid(grid_input); % call Chebgrid constructor
        set(obj,'pars',pars,'Re',pars.Re,'k',pars.k) % set pars

        % specify grid-dependent stuff
        obj.nvar = 3;
        obj.bar = obj.N*obj.nvar;
        obj.U = 1.5*(1 - obj.y.^2);
        obj.LAM = [obj.D0\obj.U; zeros(2*obj.N,1)];
        obj.Fx = 3/obj.Re;

        % coeff 'mass' matrix
        bcell = [repmat({obj.D0},1,2) zeros(obj.N)];
        obj.B = blkdiag(bcell{:});
        obj.small_spectrum = spectrum();
    end
    
    %% core methods for linear stability
    L = get_L(obj,LAM,k)
    LnB = get_LnB(obj,LAM,k,om_r)
    A = construct_linpart(obj,U,varargin)
    update_spectrum(obj)
    X = forceBCs(obj,X,type)
    X = stopinterf(obj,X)
    
    %% plots / postproc
    [u_1,u_2,du_1,du_2,Fx,c,varargout] = coeffs_to_phys(obj,U)
    to_plot = energies(obj,U)
    pointplot(obj,U)
    
    %% core methods for WNA/BC
    V = bilin(obj,phi1,phi1d,k1,phi2,phi2d,k2)
    F = F(obj,U,varargin)
    M = construct_mass(obj)
end
end