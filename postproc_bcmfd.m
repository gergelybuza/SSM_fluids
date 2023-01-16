clear
clc
close all
format long

% this script loads a bcarr and computes/plots 1D mfds at point 'ind'

bcarr = load('plots/data/test.mat').bcarr;
% bcarr = load('/store/DAMTP/gb643/mfd/k085_40x30.mat').bcarr;

set(0,'DefaultFigureVisible','off');
[xarr,Earr] = bcarr(end).postproc(bcarr,'E');
[~,Darr] = bcarr(end).postproc(bcarr,'D');
set(0,'DefaultFigureVisible','on');

ind = 3; % index along all of bcarr 
disp(xarr(ind)) % Re, or whatever you've continuated along
disp(Earr(ind)) % kin energy
spect = [bcarr.spectrum];
disp(spect(ind))

% get W, U at ind
Warr = [bcarr.Wsave]; 
uarr = [bcarr.U];
W = Warr(ind);
U = uarr(:,ind);

%% init DS, MFD
notation = 'multiindex'; % 'multiindex' or 'tensor' (slow+large memory req)
% DS object:
% DS = DynSys(W.grid_input,W.pars,U,notation); % recomputes full spectrum
DS = DynSys(W.grid_input,W.pars,U,notation,spect(ind)); % uses spectrum from BC
% spectrum = DS.spectrum.Lambda;

S = Manifold(DS); % manifold object for SSM computation
set(S.Options, 'reltol', 5,'notation',notation,'paramStyle','graph') 
masterModes = [1]; % modes along which SSM is to be computed
S.choose_E(masterModes); % specify SSM via its tangent space

%% mfds
order = 3; % SSM order
[W_0, R_0, errf] = S.compute_whisker(order); % compute SSM

% reduced dynamics (1D real)
r = [R_0{:}];
pol = [fliplr(full([r.coeffs])) 0]; % R polynomial
sol = roots(pol);
fps = sol((sol ~= 0) & (abs(sol)<10) & (abs(imag(sol))<1e-7));
fps = real(fps).'; % fixed points of R
zarr = linspace(-max(abs(fps))*1.1,max(abs(fps))*1.1,100); 

% compute plot variables
plot_full = Wztoplot_1D(zarr,W_0,W,U);
[plot_fps,fpchecks] = Wztoplot_1D(fps,W_0,W,U);
plot_U = W.energies(U);

%% plots

colors = {'#377eb8','#e41a1c','#4daf4a','#984ea3','#ff7f00','#f781bf'};

% 2d, Re-E
figure 
plot(xarr,Earr,'-k','linewidth',1,'markersize',11); hold on
plot(DS.Re,0,'.b','markersize',30);
plot(DS.Re,plot_U.E,'.r','markersize',30);
plot(DS.Re*ones(1,length(plot_full)),[plot_full.E],'-','color',colors{3});
plot(DS.Re*ones(1,length(plot_fps)),[plot_fps.E],'.','color',colors{3},'markersize',20);
xlabel('Re')
ylabel('kinetic energy')
% ylim([0 0.06])

% 3d, mfd+branch
figure
plot3(xarr,Earr,Darr,'.-k','linewidth',1,'markersize',11); hold on
plot3(DS.Re,0,0,'.b','markersize',30);
plot3(DS.Re,plot_U.E,plot_U.D,'.r','markersize',30);
plot3(DS.Re*ones(1,length(plot_full)),[plot_full.E],[plot_full.D],'color',colors{3}); 
plot3(DS.Re*ones(1,length(plot_fps)),[plot_fps.E],[plot_fps.D],'.','color',colors{3},'markersize',20);
xlabel('Re')
ylabel('kinetic energy')
zlabel('dissipation')
% ylim([0 0.06])
% zlim([0 0.3])

% save('/store/DAMTP/gb643/mfd/k085_40x30_ind134_6.mat','bcarr','R_0','W_0','ind','errf','spectrum')