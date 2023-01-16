function update_spectrum(obj)
% computes spectrum at the laminar state

A = obj.get_L(obj.LAM,obj.k);
A = obj.stopinterf(A);
A = -obj.forceBCs(A,'normal');
B = obj.stopinterf(obj.B);
obj.small_spectrum.compute_full_spectrum(A,B)
disp('Eigenvalues at the laminar state')
disp(obj.small_spectrum.Lambda(1:2))
end