function wk = update_weights(obj,y_st,y_end)
% Chebyshev integration weights

vec = 0:1:(obj.N-1);
% computes \int_{-1}^1 T_j(y)dy 
jacobian = (y_end - y_st)/2;
intTj = ((-1).^vec+1)./(1-vec.^2)*jacobian;
intTj(2) = 0;

wk = 2/(obj.N-1)*cos(vec.'*vec*pi/(obj.N-1)).*intTj;
wk(:,1) = wk(:,1)/2;        
wk(:,end) = wk(:,end)/2;
wk = sum(wk,2);
wk(1) = wk(1)/2;
wk(end) = wk(end)/2;
end