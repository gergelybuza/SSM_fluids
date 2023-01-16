function compute_full_spectrum(obj,A,varargin)
% computes the full spectrum (will be slow for large systems)
% varargin{1} = M -- mass matrix, if specified

if obj.do_adjoint
    spectrum = cell(1,3); % V D W
else
    spectrum = cell(1,2); % V D 
end

% get mass matrix and compute eigs
if ~isempty(varargin)
    M = varargin{1};
    [spectrum{:}] = eig(A,M);
else
    [spectrum{:}] = eig(A);
    M = eye(size(A));
end

% organize
spectrum = cleanse(spectrum); % remove spurious modes
spectrum = sort_modes(spectrum);            
spectrum = obj.normalize_modes(spectrum,M);

% output
obj.V = spectrum{1};
obj.Lambda = spectrum{2};
if length(spectrum) > 2
    obj.W = spectrum{3};
end

end

function spectrum = cleanse(spectrum)
% removes spurious modes

LAMBDA = diag(spectrum{2});
sp = find((abs(LAMBDA)>20)|(abs(LAMBDA)<1e-10)|(real(LAMBDA)>0.5)...
    | isnan(LAMBDA));
LAMBDA(sp) = [];
spectrum{1}(:,sp) = [];
if length(spectrum) > 2
    spectrum{3}(:,sp) = [];
end
spectrum{2} = diag(LAMBDA);
end

function spectrum = sort_modes(spectrum)
% Copied from SSMTool.
%SORT_MODES: This function sorts the eigenvectors (V) in descending order of
%the real parts of the corresponding eigenvalues (D). The resulting
%eigenvectors are also normalized to unit magnitude. 

% obtain the eigenvalues as a vector instead of a diagonal matrix
Lambda = diag(spectrum{2}); 
if ~iscolumn(Lambda)
    Lambda = transpose(Lambda);
end
% sort eigenvalues in the descending order of real parts, incase of tie by
% ascending order of magnitude of imaginary parts
[Lambda_sorted,I] = sortrows([real(Lambda), abs(imag(Lambda)) sign(imag(Lambda))],[1 2],{'descend' 'ascend'});
Lambda = Lambda_sorted(:,1) + 1i * Lambda_sorted(:,2).*Lambda_sorted(:,3);
% arrange eigenvectors accordingly
V = spectrum{1}(:,I);
if length(spectrum) > 2
    W = spectrum{3}(:,I);
end

% ensure positive imaginary part first in every complex pair
skip = false;
for j = 1:length(Lambda)-1
    if skip 
        skip = false;
        continue;
    end    
    if ~isreal(Lambda(j))&& abs(Lambda(j)-conj(Lambda(j+1)))<1e-8*abs(Lambda(j))
        % extract complex eigenpair
        Lambda0 = Lambda(j:j+1);
        % sort eigenvalues in the descending order of imaginary parts
        [~,I] = sort(imag(Lambda0),'descend','ComparisonMethod','real');        
        % rearrange the ordre of the pair
        Lambda([j,j+1]) = Lambda0(I);
        V0 = V(:,j:j+1);
        V(:,[j,j+1]) = V0(:,I);
        if length(spectrum) > 2
            W0 = W(:,j:j+1);
            W(:,[j,j+1]) = W0(:,I);
        end
%         Lambda(j) = Lambda0(I(1));
%         V(:,j) = V0(:,I(1));
%         W(:,j) = W0(:,I(1));
        % ensure complex conjugate eigenvalues and eigenvectors - not true
        % if A and B are not real
%         Lambda(j+1) = conj(Lambda(j));                
%         V(:,j+1) = conj(V(:,j));
%         W(:,j+1) = conj(W(:,j));
        % move to the next pair of eigenvalues
        skip = true; 
    end
end
spectrum = {V Lambda W};
% D = diag(Lambda);
end

