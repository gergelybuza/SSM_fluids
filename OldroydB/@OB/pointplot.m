function pointplot(obj,U)
% plots a state and some other stuff

[u_1,u_2,t_11,t_22,t_12,du_1,du_2,Fx,c,stream_func,vorticity] = obj.coeffs_to_phys(U);

x = obj.x;
y = obj.y_full;
colors = {'#377eb8','#e41a1c','#4daf4a','#984ea3','#ff7f00','#ffff33'};

trT = t_11 + t_22;
figure
[sO sO] = contourf(x, y, trT, 50);
set(sO,'linestyle','none')
hold on
[sV sV] = contour(x, y, stream_func, linspace(0.001, 0.05, 10), 'color',[0.8 0.8 0.8]);
[sV sV] = contour(x, y, stream_func, -linspace(0.001, 0.05, 10), '--','color',[0.8 0.8 0.8]);
% xlabel('$x$'); ylabel('$y$');
% title(['$\varepsilon =$' num2str(obj.epsilon) ', $Wi =$' num2str(obj.Wi)...
%     ', $Re=$' num2str(obj.Re) ', $\beta =$' num2str(obj.beta)])
plot(x(10),y(8),'.','markersize',30,'color',colors{4})
set(gcf, 'Position',  [100, 100, 560, 300])
set(gca,'xtick',[])
set(gca,'ytick',[])

inf = inferno();
colormap(inf)

colorbar
% caxis([-0.5 0.5])
set(gcf,'renderer','OpenGL')

% Chebyshev contributions
chebcontr = zeros(1,obj.N);
for j = 1:obj.N
    chebcontr(j) = norm(U(j:obj.N:obj.dim-2));
end
figure
semilogy((1:obj.N)-1,chebcontr,'k')
title('Chebyshev contributions')
xlabel('coefficient number')
ylabel('contribution')

% Fourier contributions
Fcontr = zeros(1,obj.Q);
for j = 1:obj.Q
    Fcontr(j) = norm(U(obj.bar*(j-1)+1:obj.bar*j));
end
figure
semilogy((1:obj.Q)-1,Fcontr,'k')
title('Fourier contributions')
xlabel('coefficient number')
ylabel('contribution')

end

