function [V,Lambda,W] = inv_iter_all(obj,A,M,mus,varargin)
% does inverse iteration over all requested eigenvalues in mu
% opt args: varargin{1} - eigenvectors, 
%           varargin{2} - adjoint eigenvectors

Lambda = zeros(size(mus));
V = zeros(length(M),length(mus));
W = zeros(size(V));
for j = 1:length(mus)
    switch nargin
        case 4
        [Lambda(j),V(:,j)] = obj.inverse_iter(A,M,mus(j)); % normal
        [~,W(:,j)] = obj.inverse_iter(A',M',conj(mus(j))); % adjoint
        case 5
        [Lambda(j),V(:,j)] = obj.inverse_iter(A,M,mus(j),varargin{1}(:,j));
        [~,W(:,j)] = obj.inverse_iter(A',M',conj(mus(j)));
        case 6
        [Lambda(j),V(:,j)] = obj.inverse_iter(A,M,mus(j),varargin{1}(:,j));
        [~,W(:,j)] = obj.inverse_iter(A',M',conj(mus(j)),varargin{2}(:,j));
    end
end
spectr = obj.normalize_modes({V Lambda W},M);
V = spectr{1};
W = spectr{3};
end