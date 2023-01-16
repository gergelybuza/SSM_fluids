clear
clc
close all
format long

%% define parameters
pars.Re = 3600;         % Reynolds number
pars.k = 1.02056;       % streamwise wavenumber

grid_input.N = 40;      % #Chebyshev modes
grid_input.Q = 3;       % #Fourier modes
grid_input.BC = 'full'; % Boundary Condition: 'full' or 'sym'
grid_input.TW = false;  % disable travelling wave formulation

W = NSE(grid_input,pars);   % NSE object
W.update_spectrum();        % compute (1-Fourier mode) spectrum @LAM

% base state
U_b = zeros(W.dim,1);
U_b(1:W.N*W.nvar) = W.LAM;
U_b(end) = W.Fx;

%% init DS, MFD
notation = 'multiindex'; % 'multiindex' or 'tensor' (slow+large memory req)
DS = DynSys(grid_input,pars,U_b,notation); % DS object
disp('first 10 eigenvalues')
disp(DS.spectrum.Lambda(1:2))
% spectrum = DS.spectrum.Lambda; 

S = Manifold(DS); % manifold object for SSM computation
set(S.Options, 'reltol', 1,'notation',notation,'paramStyle','graph') 
masterModes = [1 2]; % should always be [1 2] for this example
S.choose_E(masterModes); % specify SSM via its tangent space

%% mfds
order = 3; % SSM order
[W_0, R_0] = S.compute_whisker(order); % compute SSM

% reduced dynamics (polar 2D)
[rpol,tpol] = twoDtoPOLAR(R_0); % rpol - radial polynomial, tpol for theta
sol = roots(rpol);
fps = sol((sol > 0) & (abs(sol)<10) & (abs(imag(sol))<1e-7));
fps = real(fps).'; % fixed points
rmax = max(fps)*1.2;
r = linspace(0,rmax,15).';
th = linspace(0,2*pi,100).';

% grid on reduced model
x = linspace(-rmax,rmax,20);
y = linspace(-rmax,rmax,20);
[x,y] = meshgrid(x,y);
z = (x + 1i*y).';

% vector field R
Rfun = @(s) ...
subsref(reduced_to_full_complex([s; conj(s)],R_0,[],0),struct('type','()','subs',{{1}}));
R = arrayfun(Rfun,z); % evaluated at grid

% integrate 1 trajectory from zinit
zinit = z(7,7);
[t,zfun] = ode45(@(t,z) Rfun(z),[0 30],zinit);

%% plots

colors = {'#377eb8','#e41a1c','#4daf4a','#984ea3','#ff7f00','#f781bf'};
fg = figure;
% reduced model plot
subplot(1,2,1)
plot(0,0,'.','markersize',30,'color',colors{4}); hold on
plot(real(zfun(:,1)),imag(zfun(:,1)),'linewidth',0.1,'color',colors{5});
for i = 1:length(fps)
plot(fps(i)*cos(th),fps(i)*sin(th),...
    'color',colors{3},'linewidth',4);
end
quiver(real(z),imag(z),real(R),imag(R),'k')
axis tight
set(gca,'visible','off')

% compute plot variables
plot_full = Wztoplot_1D(z,W_0,W,U_b,'polar'); 
e = reshape([plot_full.E],size(plot_full));
u_1_mid_x0 = reshape([plot_full.u_1_mid_x0],size(plot_full));
u_2_avg_x0 = reshape([plot_full.u_2_avg_x0],size(plot_full));
plot_U = W.energies(U_b);
plot_fun = Wztoplot_1D(zfun(:,1),W_0,W,U_b,'polar');

% embedding plot
subplot(1,2,2)
surf(u_1_mid_x0,u_2_avg_x0,e,'FaceAlpha',0.7,...
    'edgecolor','none','facecolor','interp','facelighting','gouraud'); hold on
colormap('winter')
elim = min([e(1,:) e(:,1).']);
caxis([0 elim])
plot3([plot_fun.u_1_mid_x0],[plot_fun.u_2_avg_x0],[plot_fun.E],...
    'color',colors{5},'linewidth',0.1)
plot3(0,0,0,'.','markersize',30,'color',colors{4});
for j = 1:length(fps)
zfps = fps(j)*exp(1i*th);
plot_fps = Wztoplot_1D(zfps,W_0,W,U_b,'polar');
plot3([plot_fps.u_1_mid_x0],[plot_fps.u_2_avg_x0],[plot_fps.E],...
    'color',colors{3},'linewidth',4)
end
xlabel('$\mathrm{svf}$') % set to midline vel.; change in 'energies.m'
ylabel('$\mathrm{mwnv}$')
zlabel('$e$')
axis tight
zlim([0 elim])
ylim([-0.012 0.012])
xlim([-0.2 4.5]*1e-3)
fg.Position = [100 50 1120 450];
% save('/store/DAMTP/gb643/mfd/kcrit_40x30_Re3600r.mat','R_0','W_0','spectrum','W','U_b')


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

