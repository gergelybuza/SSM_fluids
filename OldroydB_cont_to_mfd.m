clear
clc
close all
format long

%% define parameters
pars.beta    = 0.8;         % viscosity ratio
pars.Re      = 0;           % Reynolds number
pars.Wi      = 6.53;        % Weissenberg number
pars.k       = 2.6;         % streamwise wavenumber
pars.epsilon = 1e-2;        % polymeric dissipation coefficient

grid_input.N = 35;          % #Chebyshev modes
grid_input.Q = 3;           % #Fourier modes
grid_input.BC = 'full';     % Boundary Condition: 'full' or 'sym'
grid_input.TW = true;       % enable travelling wave formulation

initstep = 0.01;             % step size for initial step along WNA
stepsize = 1;               % step size for all other steps
masterModes = [1];          % modes for which eigenvalues are continuated
% (eigenvalues are ordered according to their real part, 1 is least stable)
% here, mastermodes also specifies SSM to be computed

% specify direction (name-value array pairs)
names = {'Wi'};             % anything in 'pars' works
values = [1];               % means start upwards in 'Wi'
direction = table(names,values);

%% init from bifurcation point
W = OB(grid_input,pars); % OB object
BC = bc(direction,stepsize,W,'specmodes',masterModes); % BC object
bcarr = []; % array of BCs
[U,W] = BC.initWNA(initstep,W); % initial WNA step
spect = BC.spectrum(end);       % spectrum init

%% continuation
% cont. will stop when 'limdir' reaches 'limit' or upon reaching a fold
% (if 'breakonfold' is true)
disp('Starting loop')
[U,W] = BC.mainloop(U,W,'limdir','Wi','limit',6.5,...
    'breakonfold',false,'spect',spect);
bcarr = [bcarr BC]; % pad BCs

% BC.compareplot(0.1,bcarr);
% W.pointplot(U)

% compute mfd about point 'ind'
ind = 6;
W = BC.Wsave(ind);
U = BC.U(:,ind);

%% init DS, MFD
notation = 'multiindex'; % 'multiindex' or 'tensor' (slow+large memory req)
DS = DynSysOB(W.grid_input,W.pars,U,notation,BC.spectrum(ind)); % DSOB object
S = Manifold(DS); % manifold object for SSM computation
set(S.Options, 'reltol', 10, 'notation', notation) 
S.choose_E(masterModes); % specify SSM via its tangent space

%% mfds
order = 5; % SSM order
[W_0, R_0, errf] = S.compute_whisker(order); % compute SSM

% reduced dynamics (1D real)
r = [R_0{:}];
pol = [fliplr(full([r.coeffs])) 0]; % R polynomial
sol = roots(pol);
fps = sol((sol ~= 0) & (abs(sol)<10) & (abs(imag(sol))<1e-10));
fps = real(fps).'; % fixed points of R

% plot over z
% zarr = linspace(-max(abs(fps))*1.1,max(abs(fps))*1.1,100); 
tol = 3; % error tolerance of SSM
[zarr,err] = zarr_create(errf,1e-1,tol);

% compute plot variables
plot_full = Wztoplot_1D(zarr,W_0,W,U);
[plot_fps,fpchecks] = Wztoplot_1D(fps,W_0,W,U);
plot_U = W.energies(U);

%% plots

colors = {'#377eb8','#e41a1c','#4daf4a','#984ea3','#ff7f00','#f781bf'};
set(0,'DefaultFigureVisible','off');
[xarr,Aarr] = bcarr(end).postproc(bcarr,'A');
[~,Earr] = bcarr(end).postproc(bcarr,'E');
set(0,'DefaultFigureVisible','on');

% 2d, just mfd
figure
plot(1,0,'.b','markersize',30); hold on
plot(plot_U.A,plot_U.E,'.r','markersize',30);
plot([plot_full.A],[plot_full.E],'color',colors{3}); 
plot([plot_fps.A],[plot_fps.E],'.','color',colors{3},'markersize',20);
xlabel('$\mathcal{T}$')
ylabel('kinetic energy')

% 3d, mfd+branch
figure 
plot3(xarr,Aarr,Earr,'.-k','linewidth',1,'markersize',11); hold on
plot3(DS.Wi,1,0,'.b','markersize',30);
plot3(DS.Wi,plot_U.A,plot_U.E,'.r','markersize',30);
plot3(DS.Wi*ones(1,length(plot_full)),[plot_full.A],[plot_full.E],'color',colors{3}); 
plot3(DS.Wi*ones(1,length(plot_fps)),[plot_fps.A],[plot_fps.E],'.','color',colors{3},'markersize',20);
xlabel('Wi')
ylabel('$\mathcal{T}$')
zlabel('kinetic energy')

function [zarr,err] = zarr_create(errf,delt,tol)
% awful code that creates an array remaining in accepted tolerance
st = 10*tol;
err = 0;
while (err(1) < tol) || (err(end) < tol)
    zarr = -st:delt:st;
    err = zeros(length(errf{1}(:)),length(zarr));
    for j = 2:length(errf)+1
        err = err + (errf{j-1}(:))*(zarr.^j); 
    end
    err = vecnorm(err);
    st = 2*st;
end
zarr = zarr(err < tol);
err = err(err < tol);
end