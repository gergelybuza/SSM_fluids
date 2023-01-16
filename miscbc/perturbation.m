function U = perturbation(U,W)
% extracts perturbation above the laminar state
U(1:W.bar) = U(1:W.bar) - W.LAM;
end