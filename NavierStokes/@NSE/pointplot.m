function pointplot(obj,U)
% plots a state and some other stuff

[u_1,~,~,~,~,c,stream_func,vorticity] = obj.coeffs_to_phys(U);

x = obj.x;
y = obj.y_full;

% figure
% plot(y,u_1(:,1)-obj.U)

fg = figure;
[sO sO] = contourf(x, y, vorticity, 50);
set(sO,'linestyle','none')
hold on
[sV sV] = contour(x, y, stream_func, linspace(0.004, 0.1, 10), 'color',[0.8 0.8 0.8]);
[sV sV] = contour(x, y, stream_func, -linspace(0.004, 0.1, 10), '--','color',[0.8 0.8 0.8]);
% xlabel('$x$'); ylabel('$y$');
% title(['$c =$' num2str(c) ', $Re =$' num2str(obj.Re)])

set(gcf, 'Position',  [100, 100, 560, 300])
set(gca,'xtick',[])
set(gca,'ytick',[])

vir = viridis();
colormap(vir)

colorbar
% caxis([-2 2])
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

