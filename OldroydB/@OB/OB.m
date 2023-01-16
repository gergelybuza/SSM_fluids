classdef OB < Chebgrid
    
properties  
    %% flow parameters
    pars            % struct of parameters (Re,k,Wi,beta,epsilon)
    Re              % Reynolds number
    k               % (streamwise) wavenumber
    Wi              % Weissenberg number
    beta            % viscosity ratio
    epsilon         % polymeric dissipation coefficient
    
    %% laminar state
    U               % base (streamwise) velocity (in phys)
    Txx             % (1,1) component of base stress tensor (in phys)
    Txy             % (1,2) component of base stress tensor (in phys)
    LAM             % full laminar state (coeffs)
    Fx              % force sustaining the base state

    %% other 
    B               % diagonal coeff matrix (for RHS of linstab)
    small_spectrum  % spectrum object (for the 1-mode expansion)
end

methods
    %% constructor
    function obj = OB(grid_input,pars)
        obj = obj@Chebgrid(grid_input); % call Chebgrid constructor
        set(obj,'pars',pars,'Re',pars.Re,'k',pars.k,'Wi',pars.Wi,...
            'beta',pars.beta,'epsilon',pars.epsilon) % set pars

        % specify grid-dependent stuff
        obj.nvar = 6;
        obj.bar = obj.N*obj.nvar;
        [obj.U,obj.Txx,obj.Txy,obj.Fx] = obj.laminarOB(); % compute laminar sol
        obj.LAM = [obj.D0\obj.U; zeros(2*obj.N,1);...
            obj.D0\obj.Txx; zeros(obj.N,1); obj.D0\obj.Txy];
        
        % coeff 'mass' matrix
        bcell = [repmat({obj.Re*obj.D0},1,2) zeros(obj.N) repmat({obj.D0},1,obj.nvar-3)];
        obj.B = blkdiag(bcell{:});
        obj.small_spectrum = spectrum();
    end
    
    %% core methods for linear stability
    [U,Txx,Txy,Fx] = laminarOB(obj)
    L = get_L(obj,LAM,k) 
    LnB = get_LnB(obj,LAM,k,om_r) 
    A = construct_linpart(obj,U,varargin) 
    update_spectrum(obj) 
    varargout = Sctrickf(obj,varargin) 
    X = forceBCs(obj,X,type) 
    X = stopinterf(obj,X)
    
    %% plots / postproc
    [u_1,u_2,t_11,t_22,t_12,du_1,du_2,Fx,c,varargout] = coeffs_to_phys(obj,U) %
    to_plot = energies(obj,U) 
    pointplot(obj,U)
    pointplot2(obj,varargin)
    
    %% core methods for WNA/BC
    V = bilin(obj,phi1,phi1d,k1,phi2,phi2d,k2) 
    F = F(obj,U,varargin)
    M = construct_mass(obj)
end
end