function varargout = compareplot(obj,dista,varargin)
% plots the WNA-BC comparison
% very old (rushed) code that I don't care to fix

if isempty(varargin)
    bcarr = obj;
else
    bcarr = varargin{1};
end

% WNA
W = bcarr(1).Wsave(1);
coefarr = cell(1,length(obj.direction.names));
coefarr2 = zeros(2,length(obj.direction.names));
[~,~,~,c,d,coefarr{:}] = obj.WNA_bilinsys(W,obj.direction.names{:});
for j = 1:length(coefarr) 
coefarr2(:,j) =...
[real(c), real(d); imag(c), imag(d)]\[real(coefarr{j}); imag(coefarr{j})];
end

Warr = [bcarr.Wsave]; 
xarr = zeros(length(obj.direction.names),length(Warr));
initpoint = zeros(length(obj.direction.names),1);
for j = 1:length(obj.direction.names)
    xarr(j,:) = [Warr.(obj.direction.names{j})];
    initpoint(j) = W.(obj.direction.names{j});
end

eps = 0.1;
dist2 = linspace(0,dista/eps^2,1000);

arrwna = initpoint + obj.direction.values*eps^2*dist2;
wnaamp = eps*real(sqrt(dist2*(coefarr2(1,:)*obj.direction.values))); 

inds = (max(xarr,[],2)-min(xarr,[],2)) ~= 0;
if isempty(inds)
    disp('There is nothing to plot here, boss.')
    return
else
xarr = xarr(inds,:);
arrwna = arrwna(inds,:);
names = plotnames(obj.direction.names);
names = names(inds);
colors = {'#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'};

uarr = [bcarr.U];
unorm = zeros(1,size(uarr,2));
for j = 1:length(unorm)
unorm(j) = vecnorm(tocomp(uarr(1:end-2,j)-uarr(1:end-2,1)));
end
if length(names) < 2
    figure
    plot(arrwna,wnaamp,'k','linewidth',1.5); hold on
    plot(xarr,unorm,'+','color',colors{5},'linewidth',1.2);
    xlabel(names{1})
    ylabel('amplitude')
else
    figure
    plot3(arrwna(1,:),arrwna(2,:),wnaamp,'k','linewidth',1.5); hold on
    plot3(xarr(1,:),xarr(2,:),unorm,'color',colors{5},'linewidth',1.2);
    xlabel(names{1})
    ylabel(names{2})
    zlabel('amplitude')
    grid on
end
end

if nargout > 1
    varargout{1} = {arrwna(:) wnaamp(:)};
    varargout{2} = {xarr(:) uarr(:)};
end       

end

function names = plotnames(dirnames)
parnames = {'l' 'beta' 'Re' 'Wi' 'Sc' 'epsilon' 'lambda' 'k'};
displaynames = {'$L_{max}$' '$\beta$' '$Re$' '$Wi$' '$Sc$'...
    '$\varepsilon$' '$\lambda$' '$k$'};
names = cell(size(dirnames));
for j = 1:length(dirnames)
ind = find(strcmp(parnames,dirnames{j}));
names{j} = displaynames{ind};
end
end
