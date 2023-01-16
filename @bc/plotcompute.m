function plotcompute(obj)
% computes plotted measures of amplitude across BC

if length(obj.Wsave) ~= size(obj.U,2)
    error('Something is wrong.')
end

for j = 1:size(obj.U,2)
    obj.to_plot(j) = obj.Wsave(j).energies(obj.U(:,j)); 
end

end