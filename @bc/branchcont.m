function [U, W] = branchcont(obj,initU,initW)
% main branch continuation routine

if obj.options.adapt
    obj.stepadapt();
end
disp(['stepsize: ' num2str(obj.stepsize)])

if obj.options.tangentinterpol && ~isempty(obj.tangent)
    obj.fixtangents(0); % fixes prev. tangent if we interpolate
end

obj.updatetangent(initU,initW); % computes tangent at current point

% update phi according to tangent
U = initU + obj.stepsize*readd(initW,obj.tangent(1:end-1,end));

% update parameters along tangent
pars = initW.pars;
initp = zeros(size(obj.direction.values));
for j = 1:length(obj.direction.names)
    initp(j) = pars.(obj.direction.names{j});
end
p = initp + obj.stepsize*obj.tangent(end,end)*obj.direction.values;
for j = 1:length(obj.direction.names)
    pars.(obj.direction.names{j}) = p(j);
    if j == 1
        disptext = [obj.direction.names{j} ' = ' num2str(p(j))];
    else
        disptext =...
            [disptext ', ' obj.direction.names{j} ' = ' num2str(p(j))];
    end
end
disp(disptext)
W = feval(obj.Wname,initW.grid_input,pars); % create wna object

count = 0;
err = 1;
while err > obj.options.tol % main NR loop
count = count + 1;
   
clear DF 
if count>obj.options.maxiter || isnan(err) || err>10
    warning('Newton iterations did not converge')
    U = [];
    W = [];
    return
end

[phi,phid] = W.coeffs_to_semip(U); % U to semiphysical space
DF = obj.jacAssembly(U,phi,phid,W,obj.tangent(:,end).'); % jacobian
% last element of RHS
lastelem = obj.stepsize-...
    obj.tangent(:,end).'*[W.remove_extra_lines(U-initU); (p-initp).'*obj.direction.values];
F = [W.F(U,phi,phid); -lastelem];
dUp = -DF\F;                     % solves system
dUp = readd(W,dUp);             % readds 0s to trivial lines

% error
err = norm(dUp)/norm([perturbation(U,W); p.'*obj.direction.values]); 
disp(['Iteration: ' num2str(count) ': ||dx||/||x|| = '  num2str(err)])

% correction
dU = dUp(1:end-1);
dp = dUp(end);
U = U + dU;
p = p + dp*obj.direction.values;

% update parameters
for j = 1:length(obj.direction.names) 
    pars.(obj.direction.names{j}) = p(j);
    if j == 1
        disptext = [obj.direction.names{j} ' = ' num2str(p(j))];
    else
        disptext =...
            [disptext ', ' obj.direction.names{j} ' = ' num2str(p(j))];
    end
end
disp(disptext)
W = feval(obj.Wname,W.grid_input,pars); % update wna object 
end

% save everything
obj.W = W;
obj.U = [obj.U, U];
obj.count = [obj.count count];

% compute spectrum, if requested
if ~isempty(obj.options.specmodes)
    [phi,phid] = W.coeffs_to_semip(U);
    obj.spectrum(end+1) =...
        spectrum(-W.construct_linpart(U,phi,phid),W.construct_mass(),obj.spectrum(end));
end
end