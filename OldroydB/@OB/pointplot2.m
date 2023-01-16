function pointplot2(obj,varargin)
% plots a state and some other stuff

x = obj.x;
y = obj.y_full;

trT = cell(1,length(varargin));
stream_func = cell(1,length(varargin));
%1236
colore = {'#e41a1c','#4daf4a','#f781bf','#377eb8'};
titel = {'$t = 0$' '$t = 4590$' '$t = 7475$' '$t \to \infty$'};
fg = figure;
for i = 1:length(varargin)
[~,~,t_11,t_22,~,~,~,~,~,stream_func{i},~] = obj.coeffs_to_phys(varargin{i});
trT{i} = t_11 + t_22;
subplot(2,2,i)
[sO sO] = contourf(x, y, trT{i}, 50);
set(sO,'linestyle','none')
hold on
[sV sV] = contour(x, y, stream_func{i}, linspace(0.001, 0.05, 10), 'color',[0.8 0.8 0.8]);
[sV sV] = contour(x, y, stream_func{i}, -linspace(0.001, 0.05, 10), '--','color',[0.8 0.8 0.8]);
% xlabel('$x$'); ylabel('$y$');
% title(['$\varepsilon =$' num2str(obj.epsilon) ', $Wi =$' num2str(obj.Wi)...
%     ', $Re=$' num2str(obj.Re) ', $\beta =$' num2str(obj.beta)])
title(titel{i})

% set(gcf, 'Position',  [100, 100, 560, 300])
set(gca,'xtick',[])
set(gca,'ytick',[])
plot(x(10),y(8),'.','markersize',30,'color',colore{i})

inf = inferno();
colormap(inf)

colorbar
% caxis([-0.5 0.5])
set(gcf,'renderer','OpenGL')
end
fg.Position = [100, 100, 1120, 600];

% % Chebyshev contributions
% chebcontr = zeros(1,obj.N);
% for j = 1:obj.N
%     chebcontr(j) = norm(U(j:obj.N:obj.dim-2));
% end
% figure
% semilogy((1:obj.N)-1,chebcontr,'k')
% title('Chebyshev contributions')
% xlabel('coefficient number')
% ylabel('contribution')
% 
% % Fourier contributions
% Fcontr = zeros(1,obj.Q);
% for j = 1:obj.Q
%     Fcontr(j) = norm(U(obj.bar*(j-1)+1:obj.bar*j));
% end
% figure
% semilogy((1:obj.Q)-1,Fcontr,'k')
% title('Fourier contributions')
% xlabel('coefficient number')
% ylabel('contribution')

end

