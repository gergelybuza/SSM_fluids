function varargout = Sctrickf(obj,varargin)
% Removes the polymeric diffusion term at the boundary.
% the input(s) is usually some form of the obj.D2 term

varargout = cell(1,length(varargin));
for i = 1:length(varargin)
    A = varargin{i};
    if isequal(obj.BC,'full')
        A(1,:) = 0;
    end
    A(end,:) = 0;
    varargout{i} = A;
end
end
