function U = initITER(obj,U,W,varargin)
% initial iteration: a simple Newton-Raphson

count = 0;
err = 1;
while err > obj.options.tol
    count = count + 1;

    if count>obj.options.maxiter || isnan(err)
        error('Newton iterations did not converge')
    end
    [phi,phid] = W.coeffs_to_semip(U); % U to semiphysical space
    dU = -W.construct_linpart(U,phi,phid)\W.F(U,phi,phid); % solve
    dU = readd(W,dU);  % readds 0s to trivial lines
    
    % error
    err = norm(dU)/norm(perturbation(U,W));
    disp(['Iteration ' num2str(count)...
        ': ||dx||/||x|| = '  num2str(err)])
    
    % correction
    U = U + dU;
end

% save
obj.U = [obj.U, U];

% compute spectrum, if requested
if ~isempty(obj.options.specmodes)
    [phi,phid] = W.coeffs_to_semip(U);
    spect = varargin{1};
    spect.V = [];
    spect.W = []; % in case a different discretization is used, recompute V,W
    spect = spectrum(-W.construct_linpart(U,phi,phid),W.construct_mass(),spect);
    obj.spectrum = [obj.spectrum spect];
end
end

