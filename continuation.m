clear
clc
close all
format long

%% define parameters
pars.Re = 3847.56;      % Reynolds number
pars.k = 1.02056;       % streamwise wavenumber
% pars.Re = 5427.2;       
% pars.k = 0.85;          

% domain size
grid_input.N = 30;      % #Chebyshev modes
grid_input.Q = 3;       % #Fourier modes
grid_input.BC = 'full'; % Boundary Condition: 'full' or 'sym'
grid_input.TW = true;   % enable travelling wave formulation

initstep = 0.01;        % step size for initial step along WNA
stepsize = 0.1;         % step size for all other steps
masterModes = [1];      % modes for which eigenvalues are continuated
% (eigenvalues are ordered according to their real part, 1 is least stable)

% specify direction (name-value array pairs)
names = {'Re'};         % anything in 'pars' works
values = [-1];          % means start downwards in 'Re'
direction = table(names,values);

%% init from bifurcation point
W = NSE(grid_input,pars); % NSE object
BC = bc(direction,stepsize,W,'specmodes',masterModes); % BC object
bcarr = []; % array of BCs
[U,W] = BC.initWNA(initstep,W); % initial WNA step
spect = BC.spectrum(end);       % spectrum init

%% init from saved data
% bcarr = load('plots/data/test.mat').bcarr; % load saved array of BCs
% BC = bcarr(end); % get last BC
% U = BC.U(:,end); % get last point's state on BC
% Win = BC.Wsave(end); % get last point's NSE object on BC
% parse = Win.pars; % get last point's parameters
% spect = BC.spectrum(end); % get last point's spectrum
% 
% W = NSE(grid_input,parse); % create new object with possibly new N,Q...
% BC = bc(direction,stepsize,W,'specmodes',masterModes); % init new BC
% U = adjustphi(U,Win,W); % adjust last point's state to match current discr.

%% continuation
% cont. will stop when 'limdir' reaches 'limit' or upon reaching a fold
% (if 'breakonfold' is true)
disp('Starting loop')
[U,W] = BC.mainloop(U,W,'limdir','Re','limit',3840,...
    'breakonfold',false,'spect',spect);
bcarr = [bcarr BC]; % pad BCs

% save('plots/data/test.mat','bcarr')

%% plots
BC.compareplot(20,bcarr) % first argument is WNA x-range in figure
Up = perturbation(U,W);
W.pointplot(Up) 
