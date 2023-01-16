clear
clc
close all
format long

%% define parameters
pars.beta    = 0.9;         % viscosity ratio
pars.Re      = 0;           % Reynolds number
pars.Wi      = 13.6;        % Weissenberg number
pars.k       = 2.3;         % streamwise wavenumber
pars.epsilon = 1e-2;        % polymeric dissipation coefficient

grid_input.N = 35;          % #Chebyshev modes
grid_input.Q = 3;           % #Fourier modes
grid_input.BC = 'full';     % Boundary Condition: 'full' or 'sym'
grid_input.TW = false;      % disable travelling wave formulation

W = OB(grid_input,pars);    % OB object
W.update_spectrum();        % compute (1-Fourier mode) spectrum @LAM

% base state
U_b = zeros(W.dim,1);
U_b(1:W.N*W.nvar) = W.LAM;
U_b(end) = W.Fx;

%% init DS, MFD
notation = 'multiindex'; % 'multiindex' or 'tensor' (slow+large memory req)
specte = spectrum(); % approximate spectrum 
specte.Lambda = [W.small_spectrum.Lambda(1) conj(W.small_spectrum.Lambda(1))];
DS = DynSysOB(grid_input,pars,U_b,notation,specte); % DSOB object
disp('first 2 eigenvalues')
disp(DS.spectrum.Lambda(1:2))
spectrum = DS.spectrum.Lambda;

S = Manifold(DS); % manifold object for SSM computation
set(S.Options, 'reltol', 1,'notation',notation,'paramStyle','graph')
masterModes = [1 2]; % should always be [1 2] for this example
S.choose_E(masterModes); % specify SSM via its tangent space

%% mfds
order = 6; % SSM order
[W_0, R_0] = S.compute_whisker(order); % compute SSM

% reduced dynamics (polar 2D)
[rpol,tpol] = twoDtoPOLAR(R_0); % rpol - radial polynomial, tpol for theta
sol = roots(rpol);
fps = sol((sol > 0) & (abs(sol)<10) & (abs(imag(sol))<1e-7));
fps = real(fps).'; % fixed points

% twist starting point along its periodic orbit to match BC phase condition
[U1,z1] = twistfp(min(fps),W,W_0,U_b,rpol,tpol);

%% initiate BC from obtained FP
grid_input.TW = true;       % re-enable TW formulation
W = OB(grid_input,pars);    % and create a new OB object with it

% specify direction (name-value array pairs)
names = {'Wi'};             % anything in 'pars' works
values = [1];               % means start upwards in 'Wi'
direction = table(names,values);
stepsize = 0.3;             % step size 

BC = bc(direction,stepsize,W); % BC object
bcarr = []; % array of BCs

%% continuation
% cont. will stop when 'limdir' reaches 'limit' or upon reaching a fold
% (if 'breakonfold' is true)
disp('Starting loop')
[U,W] = BC.mainloop(U1,W,'limdir','Wi','limit',13,'breakonfold',false);
bcarr = [bcarr BC]; 
% BC.compareplot(0.2,bcarr)
BC.postproc(bcarr,'A');
% W.pointplot(U)

function [U,z] = twistfp(fp,W,W_0,U_b,rpol,tpol)
% awful code that twists FP to match phase condition of BC
th = linspace(0,2*pi,10000);
Wfun = @(x) reduced_to_full([x; conj(x)],W_0,[],0);
pertr = zeros(1,length(th));
for i = 1:length(th)
    z = fp*exp(1i*th(i));
    eh = readd(W,Wfun(z));
    pert = W.coeffs_to_semip(eh);
    pertr(i) = imag(pert(W.bar+W.r));
    if i > 1
    if pertr(i)*pertr(i-1) < 0
%         U = U_b+[toreal(pert*conj(pertr(i)/abs(pertr(i)))); eh(end)];
        U = U_b + eh;
        break
    end
    end
end

[~, ome] = polareval(fp,rpol,tpol);
U = [U; ome];
end

function [yr, yth] = polareval(r,rpol,tpol)
yr = polyval(rpol,r);
yth = polyval(tpol,r);
end

function [rpol,tpol] = twoDtoPOLAR(R_0)
order = length(R_0);
% zarr = [z conj(z)];
rpol = zeros(1,length(R_0)+1);
tpol = zeros(1,length(R_0)+1);
rpol(end) = 0;
for j = 1:order
    if ~isempty(R_0{j}.coeffs)
    coeffs = full(R_0{j}.coeffs);
    ind = full(R_0{j}.ind);
    i = find(coeffs(1,:));
    if (length(i)>1) || (ind(i,1)-ind(i,2) ~= 1)
        error('no bueno')
    else
        rpol(end-j) = real(coeffs(1,i));
        tpol(end-j+1) = imag(coeffs(1,i));
    end
    end
end
end
