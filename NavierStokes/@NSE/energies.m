function to_plot = energies(obj,U)
% computes energy measures to plot

[u_1,u_2,du_1,du_2,to_plot.Fx,to_plot.c] = obj.coeffs_to_phys(U);
ulam = repmat(obj.U,1,length(obj.x));
dulam = repmat(obj.dif*obj.U,1,length(obj.x));
en = (u_1-ulam).^2 + u_2.^2;
dis = (du_1-dulam).^2 + du_2.^2;
E = -0.5*trapz(obj.y_full,trapz(obj.x,en,2))/(2*pi/obj.k*2);
D = sqrt(abs(-trapz(obj.y_full,trapz(obj.x,dis,2))/(2*pi/obj.k*2))); 
to_plot.E = E;
to_plot.D = D;

u_2_x0 = u_2(:,1);
u_1_mid_x0 = u_1(round(obj.N/2),1)-obj.U(round(obj.N/2));
u_2_avg_x0 = -trapz(obj.y_full,u_2_x0)/2; 
to_plot.u_1_mid_x0 = u_1_mid_x0;
to_plot.u_2_avg_x0 = u_2_avg_x0;

end