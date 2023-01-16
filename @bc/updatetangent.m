function updatetangent(obj,U,W)
% updates tangent

% check if next tangent vector is to be interpolated
not_interpol =...
    (~obj.options.tangentinterpol) || (size(obj.tangent,2)<obj.options.numpts);

% set prev. tangent
if isempty(obj.tangent) % if no tangent is available at prev. step
    tangent = [zeros(1,W.dim_chopped), 1];
elseif not_interpol
    tangent = obj.tangent(:,end).';
end

if not_interpol % standard
    [phi,phid] = W.coeffs_to_semip(U);
    Dp = obj.jacAssembly(U,phi,phid,W,tangent); % Jacobian
    vec = [zeros(W.dim_chopped,1); 1]; % RHS
    ker = Dp\vec; % solves system
    obj.tangent = [obj.tangent ker/norm(ker)]; % updates tangent array
else % interpolating, if enough points are available and enabled
    disp('interpolating tangent')
    % construct an arclength parameter:
    % (this is numpts+1 long, last one is where to interp. to)
    U = obj.U(:,end-obj.options.numpts:end);
    U = W.remove_extra_lines(U,'more');
    vars = obj.direction.names;
    vararr = zeros(length(vars),obj.options.numpts+1);
    for j = 1:length(vars)
        for i = 0:obj.options.numpts
            vararr(j,i+1) = obj.Wsave(end-obj.options.numpts+i).(vars{j});
        end
    end
    vec = [U; vararr];
    param = zeros(1,obj.options.numpts+1); % arclength parameter
    for j = 2:obj.options.numpts+1
        param(j) = norm(vec(:,j)-vec(:,j-1))+param(j-1);
    end

    % interpolate tangent:
    newtang = interp1(param(1:end-1),obj.tangent(:,end-obj.options.numpts+1:end).',...
        param(end),'linear','extrap');
    obj.tangent = [obj.tangent newtang.'/norm(newtang.')]; % update      
end
disp(['tangend: ' num2str(obj.tangent(end,end))])
end
